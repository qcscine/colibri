/**
 * @file OnTheFlyPes.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "OnTheFlyPes.h"
#include "Core/Interfaces/Calculator.h"
#include "Utils/Constants.h"
#include "Utils/UniversalSettings/OptimizationSettingsNames.h"
#include <Core/Log.h>
#include <Core/ModuleManager.h>
#include <Utils/ExternalQC/Orca/OrcaCalculator.h>
#include <Utils/ExternalQC/Orca/OrcaCalculatorSettings.h>
#include <Utils/ExternalQC/Orca/OrcaHessianOutputParser.h>
#include <Utils/ExternalQC/Turbomole/TurbomoleCalculator.h>
#include <Utils/ExternalQC/Turbomole/TurbomoleCalculatorSettings.h>
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
#include <boost/filesystem.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#ifdef MPI_PARALLEL
  #include <mpi.h>
#endif

namespace Scine::Colibri {

OTFPes::OTFPes(VibrationalParameters& vibParms) {
  int id = 0;
#ifdef MPI_PARALLEL
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  double startOpt = MPI_Wtime();
#endif
  if (id == 0) {
    std::cout << "Now setting up on-the-fly PES" << std::endl;
  }

  // First, read in starting geometry and set up calculator
  setStartGeometry(vibParms.xyzFname_);
  Utils::PositionCollection refPositions = referenceGeometry_.getPositions();
  auto& manager = Core::ModuleManager::getInstance();

  // Get Calculator
  try {
    pesCalculator_ =
        manager.get<Core::Calculator>(Core::Calculator::supports(vibParms.otfMethodFamily_), vibParms.otfProgram_);
  }
  catch (...) {
    if (id == 0) {
      std::cout << "No module named '" << vibParms.otfProgram_ << "' providing '" << vibParms.otfMethodFamily_
                << "' is currently loaded." << std::endl;
      std::cout << "Please add the module to the MODULE_PATH in order for it to be accessible." << std::endl;
    }
    exit(1);
  }

  pesCalculator_->setLog(Core::Log::silent());
  if (vibParms.otfProgram_ == "Turbomole") {
    pesCalculator_->settings().modifyBool(Utils::SettingsNames::scfDamping, true);
    pesCalculator_->settings().modifyInt(Utils::SettingsNames::maxScfIterations, 200);
    pesCalculator_->settings().modifyString(Utils::SettingsNames::basisSet, "aug-cc-pVDZ");
    if (vibParms.otfMethodFamily_ == "DFT") {
      pesCalculator_->settings().modifyString(Utils::SettingsNames::method, "b3lyp-d3bj");
    } else if (vibParms.otfMethodFamily_ == "HF") {
      pesCalculator_->settings().modifyString(Utils::SettingsNames::method, "HF");
      auto* calcBase = pesCalculator_.get();
      if (dynamic_cast<Utils::ExternalQC::TurbomoleCalculator*>(calcBase)) {
        pesCalculator_->settings().modifyBool(Utils::ExternalQC::SettingsNames::enableRi, false);
      }
    } else {
      throw std::runtime_error("The specified method family is not (yet) acciessible via the "
                               "turbomole calculator. Abort.");
    }

    pesCalculator_->settings().modifyDouble(Utils::SettingsNames::selfConsistenceCriterion, 1e-8);
  } else if (vibParms.otfProgram_ == "Orca") {
    if (vibParms.otfMethodFamily_ == "CC") {
      pesCalculator_->settings().modifyString(Utils::SettingsNames::method, vibParms.otfMethod_);
      pesCalculator_->settings().modifyString(Utils::SettingsNames::basisSet, "cc-pvtz-f12");
      pesCalculator_->settings().modifyString(Scine::Utils::ExternalQC::SettingsNames::orcaAuxCBasisSet, "aug-cc-pvtz");
      pesCalculator_->settings().modifyString(Scine::Utils::ExternalQC::SettingsNames::orcaCabsBasisSet,
                                              "cc-pvtz-f12-cabs");
      pesCalculator_->settings().modifyInt(Utils::SettingsNames::externalProgramNProcs, 1);
      pesCalculator_->settings().modifyInt(Utils::SettingsNames::externalProgramMemory, 8000);
      pesCalculator_->settings().modifyString(Scine::Utils::ExternalQC::SettingsNames::specialOption,
                                              "VeryTightSCF TightPNO");
      pesCalculator_->settings().modifyDouble(Utils::SettingsNames::selfConsistenceCriterion, 1E-9);
    } else {
      throw std::runtime_error("The specified method family is not (yet) acciessible via the "
                               "colibri module with the orca calculator. Abort.");
    }
  }

  Utils::HessianMatrix hessian = Eigen::MatrixXd::Zero(refPositions.size(), refPositions.size());

  // For the geometry optimization, only do this part once on the master thread
  if (id == 0) {
    // Set initial structure
    pesCalculator_->setStructure(referenceGeometry_);
    std::cout << "Starting structure: ";
    Utils::XyzStreamHandler::write(std::cout, *(pesCalculator_->getStructure()));
    std::cout << std::endl;

    // Second, optimize the structure with the given calculator
    if (vibParms.doTransitionStateOptimization_) {
      Utils::GeometryOptimizer<Utils::Bofill> optimizer(*pesCalculator_);
      auto settings = optimizer.getSettings();
      optimizer.setSettings(applyVeryThightProfile(settings));
      optimizer.optimize(*(pesCalculator_->getStructure()), pesCalculator_->getLog());
    } else if (vibParms.otfMethodFamily_ == "CC") {
      std::cout << "For CC methods, the geometry is not optimized, the "
                   "reference energy is read from the input, and the Hessian "
                   "is read in from file."
                << std::endl;
      refEnergy_ = vibParms.refEnergy_;
    } else {
      Utils::GeometryOptimizer<Utils::Bfgs> optimizer(*pesCalculator_);
      auto settings = optimizer.getSettings();
      optimizer.setSettings(applyVeryThightProfile(settings));
      optimizer.optimize(*(pesCalculator_->getStructure()), pesCalculator_->getLog());
    }

    referenceGeometry_ = *(pesCalculator_->getStructure());
    std::cout << "Optimized reference structure: ";
    Utils::XyzStreamHandler::write(std::cout, referenceGeometry_);
    auto centerOfMass = Utils::Geometry::Properties::getCenterOfMass(referenceGeometry_);
    std::cout << std::endl << "COM: " << centerOfMass << std::endl;
    refPositions = referenceGeometry_.getPositions();
    if (vibParms.otfMethodFamily_ != "CC") {
      refEnergy_ = pesCalculator_->results().get<Utils::Property::Energy>();
    }
    std::cout << std::endl << "Reference energy: " << refEnergy_ << std::endl << std::endl;
    if (vibParms.otfProgram_ == "Orca" && vibParms.otfMethodFamily_ == "CC") {
      auto* calcBase = pesCalculator_.get();
      if (dynamic_cast<Utils::ExternalQC::OrcaCalculator*>(calcBase)) {
        hessian = Utils::ExternalQC::OrcaHessianOutputParser::getHessian(vibParms.hessianFname_);
      }
    } else {
      // Third, calculate energy and get hessian at minimum
      pesCalculator_->setRequiredProperties(Utils::Property::Hessian | Utils::Property::Energy);
      pesCalculator_->calculate("hess");
      hessian = pesCalculator_->results().get<Utils::Property::Hessian>();
    }
  }

#ifdef MPI_PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
  double endOpt = MPI_Wtime();
  if (id == 0 && vibParms.otfMethodFamily_ != "CC")
    std::cout << "Optimisation took " << endOpt - startOpt << " seconds." << std::endl;
  MPI_Bcast(&refEnergy_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(refPositions.data(), refPositions.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(hessian.data(), hessian.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  referenceGeometry_.setPositions(refPositions);
  pesCalculator_->setStructure(referenceGeometry_);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // Fourth, get the normal frequencies and unnormalized MW normal modes
  normalModes_ = Utils::NormalModeAnalysis::calculateNormalModes(hessian, referenceGeometry_.getElements(),
                                                                 referenceGeometry_.getPositions(), false);
  auto wavenumbers = normalModes_.getWaveNumbers();
  double harmZPVE = 0;
  if (id == 0) {
    std::cout << "  Normal Frequencies:\n  [Rot. and trans. freq. removed, "
                 "imaginary freq. shown as negative.]\n"
              << std::endl;
  }
  if (id == 0) {
    printf("  %4s %8s\n", "#", "cm^-1");
  }
  for (unsigned int i = 0; i < wavenumbers.size(); i++) {
    vibParms.modeFreqs_.push_back(std::abs(wavenumbers[i]) / Utils::Constants::invCentimeter_per_hartree);
    if (id == 0) {
      printf("  %4d %+13.6f\n", i + 1, wavenumbers[i]);
    }
    if (wavenumbers[i] > 0) {
      harmZPVE += 0.5 * wavenumbers[i];
    }
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
  if (id == 0) {
    std::cout << std::endl << "Finished preparing for on-the-fly PES calculation" << std::endl;
  }

#ifdef MPI_PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
  double setUp = MPI_Wtime();
  if (id == 0)
    std::cout << "Setup took " << setUp - startOpt << " seconds." << std::endl;
#endif
  // All set and ready to go!
}

OTFPes::OTFPes(const OTFPes& rhs) : CloneInterface(rhs) {
  normalModes_ = rhs.normalModes_;
  referenceGeometry_ = rhs.referenceGeometry_;
  refEnergy_ = rhs.refEnergy_;
  pesCalculator_ = rhs.pesCalculator_->clone();
  numSPsForVSCF_ = rhs.numSPsForVSCF_;
}

OTFPes& OTFPes::operator=(const OTFPes& rhs) {
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

double OTFPes::getEnergy(const Utils::AtomCollection& displacedGeom) const {
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
    std::string dirName;
    if (dynamic_cast<Utils::ExternalQC::TurbomoleCalculator*>(calcBase)) {
      dirName = (dynamic_cast<Utils::ExternalQC::TurbomoleCalculator*>(calcBase))->getCalculationDirectory();
    } else if (dynamic_cast<Utils::ExternalQC::OrcaCalculator*>(calcBase)) {
      dirName = (dynamic_cast<Utils::ExternalQC::OrcaCalculator*>(calcBase))->getCalculationDirectory();
    }
    if (!dirName.empty()) {
      boost::filesystem::path dirPath(dirName);
      if (boost::filesystem::exists(dirPath)) {
        boost::filesystem::remove_all(dirPath);
      }
    }
  }
  return en;
}

// One-mode variation
double OTFPes::getPES(double q, int mode) const {
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
double OTFPes::getPES(double qi, double qj, int modei, int modej) const {
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
double OTFPes::getPES(double qi, double qj, double qk, int modei, int modej, int modek) const {
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

void OTFPes::setStartGeometry(std::string filename) {
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

Utils::AtomCollection OTFPes::getGeometryFromDisplacement(std::vector<int>&& idxModes, std::vector<double>&& steps) const {
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

Utils::Settings OTFPes::applyVeryThightProfile(Utils::Settings& settings) {
  settings.modifyDouble(Utils::SettingsNames::Optimizations::Convergence::stepMaxCoeff, 2.0e-5);
  settings.modifyDouble(Utils::SettingsNames::Optimizations::Convergence::stepRMS, 1.0e-5);
  settings.modifyDouble(Utils::SettingsNames::Optimizations::Convergence::gradMaxCoeff, 2.0e-5);
  settings.modifyDouble(Utils::SettingsNames::Optimizations::Convergence::gradRMS, 1.0e-5);
  settings.modifyDouble(Utils::SettingsNames::Optimizations::Convergence::deltaValue, 1.0e-7);
  settings.modifyInt(Utils::SettingsNames::Optimizations::Convergence::maxIter, 1000);
  settings.modifyInt(Utils::SettingsNames::Optimizations::Convergence::requirement, 4);
  return settings;
}

int OTFPes::getNumSPsForVSCF() const {
  return numSPsForVSCF_;
}

void OTFPes::setNumSPsForVSCF(int num) {
  numSPsForVSCF_ = num;
}

int OTFPes::convertToIntKey(double value) {
  value = round(value * 1000.0) / 1000.0;
  int converted = (int)(value * 100);
  return converted;
}

double OTFPes::round(double val) {
  if (val < 0) {
    return ceil(val - 0.5);
  }
  return floor(val + 0.5);
}

} // namespace Scine::Colibri
