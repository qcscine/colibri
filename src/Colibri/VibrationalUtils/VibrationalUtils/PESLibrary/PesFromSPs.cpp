/**
 * @file PesFromSPs.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#include "PesFromSPs.h"
#include "Core/Interfaces/Calculator.h"
#include "Utils/Constants.h"
#include <Core/Log.h>
#include <Core/ModuleManager.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Geometry/GeometryUtilities.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/DoubleDescriptor.h>
#include <Utils/UniversalSettings/ParametrizedOptionValue.h>
#include <Utils/UniversalSettings/SettingPopulator.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#ifdef MPI_PARALLEL
  #include <mpi.h>
#endif

namespace Scine::Colibri {

PesFromSPs::PesFromSPs(VibrationalParameters& vibParms) {
  int id = 0;
#ifdef MPI_PARALLEL
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  double startOpt = MPI_Wtime();
#endif
  if (id == 0) {
    std::cout << "Now setting up PES from precalculated SinglePoints" << std::endl;
  }

  auto wavenumbers = readHarmFreqs(vibParms.harmFreqFname_, vibParms.numModes_);
  if (id == 0) {
    std::cout << "Finished reading in harmonic frequencies" << std::endl;
  }
  double harmZPVE = 0;
  if (id == 0) {
    std::cout << std::endl
              << "  Normal Frequencies:\n  [Rot. and trans. freq. removed, "
                 "imaginary freq. shown as negative.]\n"
              << std::endl;
  }
  if (id == 0) {
    printf("  %4s %8s\n", "#", "cm^-1");
  }
  for (unsigned int i = 0; i < wavenumbers.size(); i++) {
    vibParms.modeFreqs_[i] = (wavenumbers[i] / (Utils::Constants::invCentimeter_per_hartree));
    if (id == 0) {
      printf("  %4d %+13.6f\n", i + 1, wavenumbers[i]);
    }
    harmZPVE += 0.5 * wavenumbers[i];
  }
  if (id == 0) {
    std::cout << std::endl
              << "Harmonic ZPVE is: " << harmZPVE / Utils::Constants::invCentimeter_per_hartree << " Hartree" << std::endl;
  }
  if (id == 0) {
    std::cout << "Harmonic ZPVE is: " << harmZPVE << " cm^-1" << std::endl;
  }

  int numSPs = readSinglePoints(vibParms.singlePointsFname_);
  if (id == 0) {
    std::cout << "Finished reading in SinglePoints" << std::endl;
  }
  if (id == 0) {
    std::cout << "Num terms in 1B: " << calc1BSPs_.size() << ", in 2B: " << calc2BSPs_.size()
              << ", in 3B: " << calc3BSPs_.size() << std::endl;
  }

  refEnergy_ = vibParms.refEnergy_;
  if (refEnergy_ == 0) {
    throw std::runtime_error("You need to provide the reference energy in the input file for a PES "
                             "constructed from precalculated SinglePoints.");
  }

  if (id == 0) {
    std::cout << std::endl
              << "Finished preparing for PES from precalculated SinglePoints "
                 "calculation"
              << std::endl;
  }
#ifdef MPI_PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
  double setUp = MPI_Wtime();
  if (id == 0)
    std::cout << "Setup took " << setUp - startOpt << " seconds." << std::endl;
#endif
  // All set and ready to go!
}

std::vector<double> PesFromSPs::readHarmFreqs(const std::string& filename, int numModes) {
  // Open the integral file and do a loop over this file
  std::ifstream freqFile;
  std::string someString;
  freqFile.open(filename);
  if (!freqFile) {
    throw std::runtime_error("Error: input file " + filename + " could not be opened!");
  }

  freqFile >> someString >> someString; // Header, first line is skipped

  int mode = 0;
  double freq = 0.0;
  std::vector<double> harmFreqs;
  for (unsigned int i = 1; i <= numModes; i++) {
    freqFile >> mode >> freq;
    if (mode != i) {
      throw std::runtime_error("Error: something went wrong with reading in "
                               "the harmonic frequencies!");
    }
    harmFreqs.push_back(freq);
  }

  freqFile.close();
  return harmFreqs;
}

int PesFromSPs::readSinglePoints(const std::string& filename) {
  // Open the single point results file and do a loop over this file
  std::ifstream singlePointStream;
  std::string lineString;
  singlePointStream.open(filename);
  if (!singlePointStream) {
    throw std::runtime_error("Error: input file " + filename + " could not be opened!");
  }
  while (getline(singlePointStream, lineString)) {
    // Extract the data
    parseLineOfSinglePoint(lineString);
  }
  singlePointStream.close();
  return calc1BSPs_.size() + calc2BSPs_.size() + calc3BSPs_.size();
}

void PesFromSPs::parseLineOfSinglePoint(std::string lineString) {
  // Trim leading and final spaces in the string.
  lineString.erase(lineString.begin(),
                   std::find_if(lineString.begin(), lineString.end(), [&](int ch) { return std::isspace(ch) == 0; }));
  lineString.erase(std::find_if(lineString.rbegin(), lineString.rend(), [&](int ch) { return std::isspace(ch) == 0; }).base(),
                   lineString.end());

  // Split the string into modes, displacments and energy
  std::string delimiter = " --- ";
  std::string modeTerm = lineString.substr(0, lineString.find(delimiter));
  lineString.erase(0, lineString.find(delimiter) + delimiter.length());

  std::string dispTerm = lineString.substr(0, lineString.find(delimiter));
  lineString.erase(0, lineString.find(delimiter) + delimiter.length());

  std::string energyTerm = lineString;

  // Get order
  int modeCouplingOrder = 1;
  std::string::size_type pos = 0;
  std::string target = "and";
  while ((pos = modeTerm.find(target, pos)) != std::string::npos) {
    ++modeCouplingOrder;
    pos += target.length();
  }

  std::vector<std::string> energyTermVec;
  boost::split(energyTermVec, energyTerm, boost::is_any_of(":"), boost::token_compress_on);
  double energyOfSP = std::atof(energyTermVec[energyTermVec.size() - 1].c_str());

  std::vector<std::string> modeTermVec;
  boost::split(modeTermVec, modeTerm, boost::is_any_of(":"), boost::token_compress_on);
  modeTerm = modeTermVec[1];

  std::vector<std::string> dispTermVec;
  boost::split(dispTermVec, dispTerm, boost::is_any_of(":"), boost::token_compress_on);
  dispTerm = dispTermVec[1];

  if (modeCouplingOrder == 1) {
    int mode = std::atoi(modeTerm.c_str());
    int disp = convertToIntKey(std::atof(dispTerm.c_str()));
    std::array<std::pair<int, int>, 1> key = {std::make_pair(mode, disp)};
    auto find = calc1BSPs_.find(key);
    if (find == calc1BSPs_.end()) {
      calc1BSPs_[key] = energyOfSP;
    } else {
      throw std::runtime_error("Tried to add 1B integral term which was "
                               "already in map. Are there dublicates?!");
    }
  } else {
    std::string seperator = " and ";
    std::vector<int> modes;
    std::vector<int> disps;
    std::string mode;
    std::string disp;
    for (int i = 1; i < modeCouplingOrder; i++) {
      mode = modeTerm.substr(0, modeTerm.find(seperator));
      modes.push_back(std::atoi(mode.c_str()));
      modeTerm.erase(0, modeTerm.find(seperator) + seperator.length());
      disp = dispTerm.substr(0, dispTerm.find(seperator));
      disps.push_back(convertToIntKey(std::atof(disp.c_str())));
      dispTerm.erase(0, dispTerm.find(seperator) + seperator.length());
    }
    modes.push_back(std::atoi(modeTerm.c_str()));                  // removed all previous mode indices and "and"s
    disps.push_back(convertToIntKey(std::atof(dispTerm.c_str()))); // dito

    if (modeCouplingOrder == 2) {
      std::array<std::pair<int, int>, 2> key = {std::make_pair(modes[0], disps[0]), std::make_pair(modes[1], disps[1])};
      std::sort(key.begin(), key.end(), [](std::pair<int, int> i, std::pair<int, int> j) {
        return ((i.first < j.first) || (i.first == j.first && i.second < j.second));
      });
      auto find = calc2BSPs_.find(key);
      if (find == calc2BSPs_.end()) {
        calc2BSPs_[key] = energyOfSP;
      } else {
        throw std::runtime_error("Tried to add 2B integral term which was "
                                 "already in map. Are there dublicates?!");
      }
    } else if (modeCouplingOrder == 3) {
      std::array<std::pair<int, int>, 3> key = {std::make_pair(modes[0], disps[0]), std::make_pair(modes[1], disps[1]),
                                                std::make_pair(modes[2], disps[2])};
      std::sort(key.begin(), key.end(), [](std::pair<int, int> i, std::pair<int, int> j) {
        return ((i.first < j.first) || (i.first == j.first && i.second < j.second));
      });
      auto find = calc3BSPs_.find(key);
      if (find == calc3BSPs_.end()) {
        calc3BSPs_[key] = energyOfSP;
      } else {
        throw std::runtime_error("Tried to add 3B integral term which was "
                                 "already in map. Are there dublicates?!");
      }
    } else {
      throw std::runtime_error("Modecouplingorder of term in singlepoints file "
                               "not recognized. Only 1,2,3 are valid. Abort");
    }
  }
}

PesFromSPs::PesFromSPs(const PesFromSPs& rhs) : CloneInterface(rhs) {
  refEnergy_ = rhs.refEnergy_;
}

PesFromSPs& PesFromSPs::operator=(const PesFromSPs& rhs) {
  if (this == &rhs) {
    return *this;
  }
  refEnergy_ = rhs.refEnergy_;
  return *this;
}

// One-mode variation
double PesFromSPs::getPES(double q, int mode) const {
  if (std::abs(q) < 0.001) {
    return 0;
  }
  double Vmi = 0;
  std::array<std::pair<int, int>, 1> key = {std::make_pair(mode, convertToIntKey(q))};
  auto find = calc1BSPs_.find(key);
  if (find != calc1BSPs_.end()) {
    Vmi = find->second;
  } else {
    throw std::runtime_error("Could not find 1B entry for key " + std::to_string(key[0].first) + " " +
                             std::to_string(key[0].second));
  }
  Vmi -= refEnergy_;
  return Vmi;
}

// Two-mode variation
double PesFromSPs::getPES(double qi, double qj, int modei, int modej) const {
  if (std::abs(qi) < 0.001 || std::abs(qj) < 0.001) {
    return 0;
  }
  double Vmimj = 0;
  std::array<std::pair<int, int>, 2> key = {std::make_pair(modei, convertToIntKey(qi)),
                                            std::make_pair(modej, convertToIntKey(qj))};
  std::sort(key.begin(), key.end(), [](std::pair<int, int> i, std::pair<int, int> j) {
    return ((i.first < j.first) || (i.first == j.first && i.second < j.second));
  });
  auto find = calc2BSPs_.find(key);
  if (find != calc2BSPs_.end()) {
    Vmimj = find->second;
  } else {
    throw std::runtime_error("Could not find 2B entry for key " + std::to_string(key[0].first) + " " +
                             std::to_string(key[0].second) + " " + std::to_string(key[1].first) + " " +
                             std::to_string(key[1].second));
  }
  Vmimj -= refEnergy_;
  Vmimj = Vmimj - getPES(qi, modei) - getPES(qj, modej);
  return Vmimj;
}

// Three-mode variation
double PesFromSPs::getPES(double qi, double qj, double qk, int modei, int modej, int modek) const {
  if (std::abs(qi) < 0.001 || std::abs(qj) < 0.001 || std::abs(qk) < 0.001) {
    return 0;
  }
  double Vmimjmk = 0;
  std::array<std::pair<int, int>, 3> key = {std::make_pair(modei, convertToIntKey(qi)),
                                            std::make_pair(modej, convertToIntKey(qj)),
                                            std::make_pair(modek, convertToIntKey(qk))};
  std::sort(key.begin(), key.end(), [](std::pair<int, int> i, std::pair<int, int> j) {
    return ((i.first < j.first) || (i.first == j.first && i.second < j.second));
  });
  auto find = calc3BSPs_.find(key);
  if (find != calc3BSPs_.end()) {
    Vmimjmk = find->second;
  } else {
    throw std::runtime_error("Could not find 3B entry for key " + std::to_string(key[0].first) + " " +
                             std::to_string(key[0].second) + " " + std::to_string(key[1].first) + " " +
                             std::to_string(key[1].second) + " " + std::to_string(key[2].first) + " " +
                             std::to_string(key[2].second));
  }
  Vmimjmk -= refEnergy_;
  Vmimjmk = Vmimjmk - getPES(qi, qj, modei, modej) - getPES(qi, qk, modei, modek) - getPES(qj, qk, modej, modek);
  Vmimjmk = Vmimjmk - getPES(qi, modei) - getPES(qj, modej) - getPES(qk, modek);
  return (Vmimjmk);
}

int PesFromSPs::convertToIntKey(double value) {
  value = round(value * 1000.0) / 1000.0;
  int converted = (int)(value * 100);
  return converted;
}

double PesFromSPs::round(double val) {
  if (val < 0) {
    return ceil(val - 0.5);
  }
  return floor(val + 0.5);
}

} // namespace Scine::Colibri
