/**
 * @file HybridOTFfromSPs.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "HybridOTFfromSPs.h"
#include "Core/Interfaces/Calculator.h"
#include "Utils/Constants.h"
#include <Core/Log.h>
#include <Core/ModuleManager.h>
#include <Utils/ExternalQC/Turbomole/TurbomoleCalculator.h>
#include <Utils/GeometricDerivatives/NormalModeAnalysis.h>
#include <Utils/GeometricDerivatives/NumericalHessianCalculator.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Geometry/GeometryUtilities.h>
#include <Utils/GeometryOptimization/GeometryOptimization.h>
#include <Utils/GeometryOptimization/GeometryOptimizer.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Optimizer/GradientBased/Bfgs.h>
#include <Utils/Optimizer/GradientBased/GradientBasedCheck.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/DoubleDescriptor.h>
#include <Utils/UniversalSettings/ParametrizedOptionValue.h>
#include <Utils/UniversalSettings/SettingPopulator.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#ifdef MPI_PARALLEL
  #include <mpi.h>
#endif

namespace Scine::Colibri {

HybridPes::HybridPes(VibrationalParameters& vibParms) {
  int id = 0;
#ifdef MPI_PARALLEL
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  double startOpt = MPI_Wtime();
#endif
  if (id == 0) {
    std::cout << "Now setting up hybrid PES combining precalculated "
                 "SinglePoints and on-the-fly calculations"
              << std::endl;
  }

  // First, read in molecular geometry
  setGeometry(vibParms.xyzFname_);
  if (id == 0) {
    std::cout << "Structure is set to: ";
    Utils::XyzStreamHandler::write(std::cout, referenceGeometry_);
    std::cout << std::endl;
  }

  // Read in harmonic frequencies for later check
  std::vector<double> wavenumbers = readHarmFreqs(vibParms.harmFreqFname_, vibParms.numModes_);
  if (id == 0) {
    std::cout << "Finished reading in harmonic frequencies" << std::endl;
  }

  // Get Calculator
  auto& manager = Core::ModuleManager::getInstance();
  try {
    pesCalculator_ =
        manager.get<Core::Calculator>(Core::Calculator::supports(vibParms.otfMethodFamily_), vibParms.otfProgram_);
  }
  catch (...) {
    if (id == 0) {
      std::cout << "No module named '" << vibParms.otfProgram_ << "' providing '" << vibParms.otfMethodFamily_
                << "' is currently loaded." << std::endl;
    }
    if (id == 0) {
      std::cout << "Please add the module to the MODULE_PATH in order for it "
                   "to be accessible.\n"
                << std::endl;
    }
    exit(1);
  }

  // Set structure
  pesCalculator_->setStructure(referenceGeometry_);

  // Set other calculator settings
  pesCalculator_->setLog(Core::Log::silent());
  if (vibParms.otfProgram_ == "Turbomole") {
    pesCalculator_->settings().modifyBool(Utils::SettingsNames::scfDamping, true);
    pesCalculator_->settings().modifyInt(Utils::SettingsNames::maxScfIterations, 200);
    pesCalculator_->settings().modifyString(Utils::SettingsNames::basisSet, "aug-cc-pVDZ");
    pesCalculator_->settings().modifyString(Utils::SettingsNames::method, "b3lyp-d3bj");
    pesCalculator_->settings().modifyDouble(Utils::SettingsNames::selfConsistenceCriterion, 1e-8);
  }

  Utils::HessianMatrix hessian =
      Eigen::MatrixXd::Zero(referenceGeometry_.getPositions().size(), referenceGeometry_.getPositions().size());

  // ¨Only do this part once on the master thread
  if (id == 0) {
    // Calculate energy and get hessian at minimum
    pesCalculator_->setRequiredProperties(Utils::Property::Hessian | Utils::Property::Thermochemistry | Utils::Property::Energy);
    pesCalculator_->calculate("hess");
    hessian = pesCalculator_->results().get<Utils::Property::Hessian>();
    refEnergy_ = pesCalculator_->results().get<Utils::Property::Energy>();
    refEnergy_ = vibParms.refEnergy_;
    if (refEnergy_ == 0) {
      throw std::runtime_error("You need to provide the reference energy in "
                               "the input file for a PES "
                               "constructed from precalculated SinglePoints.");
    }
    if (std::abs(refEnergy_ - vibParms.refEnergy_) > 1e-10) {
      throw std::runtime_error("Provided reference energy does not agree with "
                               "calculated one. Abort.");
    }
    std::cout << std::endl << "Reference energy: " << refEnergy_ << std::endl << std::endl;
  }

#ifdef MPI_PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&refEnergy_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(hessian.data(), hessian.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // Get the normal frequencies and unnormalized MW normal modes
  normalModes_ = Utils::NormalModeAnalysis::calculateNormalModes(hessian, referenceGeometry_.getElements(),
                                                                 referenceGeometry_.getPositions(), false);
  std::vector<double> wavenumbersCalc = normalModes_.getWaveNumbers();

  // Compare wavenumbers
  for (unsigned int i = 0; i < wavenumbers.size(); i++) {
    if (std::abs(wavenumbers[i] - wavenumbersCalc[i]) > 1e-4) {
      if (id == 0) {
        std::cout << "Wavenumbers - provided vs. calculated:" << std::endl;
        for (unsigned int i = 0; i < wavenumbersCalc.size(); i++) {
          printf("  %4d %+13.6f %+13.6f\n", i + 1, wavenumbers[i], wavenumbersCalc[i]);
        }
      }
      throw std::runtime_error("Calculated wavenumbers deviate from provided ones. Abort.");
    }
  }

  double harmZPVE = 0;
  if (id == 0) {
    std::cout << "  Normal Frequencies:\n  [Rot. and trans. freq. removed, "
                 "imaginary freq. shown as negative.]\n"
              << std::endl;
  }
  if (id == 0) {
    printf("  %4s %8s\n", "#", "cm^-1");
  }
  for (unsigned int i = 0; i < wavenumbersCalc.size(); i++) {
    vibParms.modeFreqs_[i] = wavenumbersCalc[i] / (Utils::Constants::invCentimeter_per_hartree);
    if (id == 0) {
      printf("  %4d %+13.6f\n", i + 1, wavenumbersCalc[i]);
    }
    harmZPVE += 0.5 * wavenumbersCalc[i];
  }
  if (id == 0) {
    std::cout << std::endl
              << "Harmonic ZPVE is: " << harmZPVE / Utils::Constants::invCentimeter_per_hartree << " Hartree" << std::endl;
  }
  if (id == 0) {
    std::cout << "Harmonic ZPVE is: " << harmZPVE << " cm^-1" << std::endl;
  }

  pesCalculator_->setRequiredProperties(Utils::Property::Energy);
  numSPsForVSCF_ = 0;

  int numSPs = readSinglePoints(vibParms.singlePointsFname_);
  if (id == 0) {
    std::cout << "Finished reading in provided precalculated SinglePoints" << std::endl;
  }
  if (id == 0) {
    std::cout << "Num terms in 1B: " << calc1BSPs_.size() << ", in 2B: " << calc2BSPs_.size()
              << ", in 3B: " << calc3BSPs_.size() << std::endl;
  }

  if (id == 0) {
    std::cout << std::endl
              << "Finished preparing for hybrid PES combining precalculated "
                 "SinglePoints with on-the-fly calculations."
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

HybridPes::HybridPes(const HybridPes& rhs) : CloneInterface(rhs) {
  normalModes_ = rhs.normalModes_;
  referenceGeometry_ = rhs.referenceGeometry_;
  refEnergy_ = rhs.refEnergy_;
  pesCalculator_ = rhs.pesCalculator_->clone();
  numSPsForVSCF_ = rhs.numSPsForVSCF_;
}

HybridPes& HybridPes::operator=(const HybridPes& rhs) {
  if (this == &rhs) {
    return *this;
  }
  normalModes_ = rhs.normalModes_;
  referenceGeometry_ = rhs.referenceGeometry_;
  refEnergy_ = rhs.refEnergy_;
  pesCalculator_ = rhs.pesCalculator_->clone();
  numSPsForVSCF_ = rhs.numSPsForVSCF_;
  return *this;
}

double HybridPes::getEnergy(const Utils::AtomCollection& displacedGeom) const {
  assert(pesCalculator_);
  auto localPesCalculator = pesCalculator_->clone();
  localPesCalculator->setStructure(displacedGeom);
  auto nOldThreads = omp_get_num_threads();
  omp_set_num_threads(1);
  Utils::Results r = localPesCalculator->calculate("");
  r = localPesCalculator->results();
  omp_set_num_threads(nOldThreads);
  auto en = r.template get<Utils::Property::Energy>();
  {
    auto* calcBase = localPesCalculator.get();
    if (dynamic_cast<Utils::ExternalQC::TurbomoleCalculator*>(calcBase)) {
      std::string dirName = (dynamic_cast<Utils::ExternalQC::TurbomoleCalculator*>(calcBase))->getCalculationDirectory();
      boost::filesystem::path dirPath(dirName);
      if (boost::filesystem::exists(dirPath)) {
        boost::filesystem::remove_all(dirPath);
      }
    }
  }
  return en;
}

// One-mode variation
double HybridPes::getPES(double q, int mode) const {
  if (std::abs(q) < 0.001) {
    return 0;
  }
  double Vmi = 0;
  std::array<std::pair<int, int>, 1> key = {std::make_pair(mode, convertToIntKey(q))};
  auto find = calc1BSPs_.find(key);
  if (find != calc1BSPs_.end()) {
    Vmi = find->second;
  } else {
    Utils::AtomCollection displacedGeom = this->getGeometryFromDisplacement({mode}, {q});
    Vmi = getEnergy(displacedGeom);
#pragma omp critical
    {
      find = calc1BSPs_.find(key);
      if (find == calc1BSPs_.end()) {
        calc1BSPs_[key] = Vmi;
        ++numSPsForVSCF_;
        std::ofstream outfile;
        outfile.open("pes1ModeCuts.out", std::ios_base::app);
        outfile << "Mode: " << std::setw(2) << mode << " --- Displacement: " << std::fixed << std::setprecision(5)
                << std::setw(9) << q;
        outfile << " --- Energy:" << std::setprecision(10) << std::setw(16) << Vmi << std::endl;
        outfile.close();
      }
    }
  }
  Vmi -= refEnergy_;
  return Vmi;
}

// Two-mode variation
double HybridPes::getPES(double qi, double qj, int modei, int modej) const {
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
    Utils::AtomCollection displacedGeom = this->getGeometryFromDisplacement({modei, modej}, {qi, qj});
    Vmimj = getEnergy(displacedGeom);
#pragma omp critical
    {
      find = calc2BSPs_.find(key);
      if (find == calc2BSPs_.end()) {
        calc2BSPs_[key] = Vmimj;
        ++numSPsForVSCF_;
        std::ofstream outfile2;
        outfile2.open("pes2ModeCuts.out", std::ios_base::app);
        outfile2 << "Modes: " << std::setw(2) << modei << " and " << std::setw(2) << modej;
        outfile2 << " --- Displacements: " << std::fixed << std::setprecision(5) << std::setw(9) << qi << " and "
                 << std::setprecision(5) << std::setw(9) << qj;
        outfile2 << " --- Energy: " << std::setprecision(10) << std::setw(16) << Vmimj << std::endl;
        outfile2.close();
      }
    }
  }
  Vmimj -= refEnergy_;
  Vmimj = Vmimj - getPES(qi, modei) - getPES(qj, modej);
  return Vmimj;
}

// Three-mode variation
double HybridPes::getPES(double qi, double qj, double qk, int modei, int modej, int modek) const {
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
    Utils::AtomCollection displacedGeom = this->getGeometryFromDisplacement({modei, modej, modek}, {qi, qj, qk});
    Vmimjmk = getEnergy(displacedGeom);
#pragma omp critical
    {
      find = calc3BSPs_.find(key);
      if (find == calc3BSPs_.end()) {
        calc3BSPs_[key] = Vmimjmk;
        ++numSPsForVSCF_;
        std::ofstream outfile3;
        outfile3.open("pes3ModeCuts.out", std::ios_base::app);
        outfile3 << "Modes: " << std::setw(2) << modei << " and " << std::setw(2) << modej << " and " << std::setw(2) << modek;
        outfile3 << " --- Displacements: " << std::fixed << std::setprecision(5) << std::setw(9) << qi << " and "
                 << std::setw(9) << qj << " and " << std::setw(9) << qk;
        outfile3 << " --- Energy:" << std::setprecision(10) << std::setw(16) << Vmimjmk << std::endl;
        outfile3.close();
      }
    }
  }
  Vmimjmk -= refEnergy_;
  Vmimjmk = Vmimjmk - getPES(qi, qj, modei, modej) - getPES(qi, qk, modei, modek) - getPES(qj, qk, modej, modek);
  Vmimjmk = Vmimjmk - getPES(qi, modei) - getPES(qj, modej) - getPES(qk, modek);
  return (Vmimjmk);
}

void HybridPes::setGeometry(std::string filename) {
  boost::filesystem::path filepath(filename);
  if (!boost::filesystem::exists(filepath)) {
    std::cout << "XYZ-File for reference geometry not found!" << std::endl;
    exit(1);
  }
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cout << "XYZ-File for reference geometry could not be opened!" << std::endl;
    exit(1);
  }
  referenceGeometry_ = Utils::XyzStreamHandler::read(file);
  file.close();
}

Utils::AtomCollection HybridPes::getGeometryFromDisplacement(std::vector<int>&& idxModes, std::vector<double>&& steps) const {
  Utils::AtomCollection currentGeometry = referenceGeometry_;
  assert(idxModes.size() == steps.size());
  Eigen::Matrix<double, Eigen::Dynamic, 3> finalPositions = referenceGeometry_.getPositions();
  for (std::size_t iModes = 0; iModes < idxModes.size(); iModes++) {
    auto modeNoWeight = normalModes_.getMode(idxModes[iModes]);
    int nAtoms = pesCalculator_->getStructure()->size();
    for (int iAtoms = 0; iAtoms < nAtoms; iAtoms++) {
      modeNoWeight(iAtoms, 0) *= std::sqrt(Utils::Constants::u_per_electronRestMass);
      modeNoWeight(iAtoms, 1) *= std::sqrt(Utils::Constants::u_per_electronRestMass);
      modeNoWeight(iAtoms, 2) *= std::sqrt(Utils::Constants::u_per_electronRestMass);
    }
    modeNoWeight *= steps[iModes];
    finalPositions += modeNoWeight;
  }
  currentGeometry.setPositions(finalPositions);
  return currentGeometry;
}

int HybridPes::getNumSPsForVSCF() const {
  return numSPsForVSCF_;
}

void HybridPes::setNumSPsForVSCF(int num) {
  numSPsForVSCF_ = num;
}

std::vector<double> HybridPes::readHarmFreqs(const std::string& filename, int numModes) {
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

int HybridPes::readSinglePoints(const std::string& filename) {
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

void HybridPes::parseLineOfSinglePoint(std::string lineString) {
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

int HybridPes::convertToIntKey(double value) {
  value = round(value * 1000.0) / 1000.0;
  int converted = (int)(value * 100);
  return converted;
}

double HybridPes::round(double val) {
  if (val < 0) {
    return ceil(val - 0.5);
  }
  return floor(val + 0.5);
}

} // namespace Scine::Colibri