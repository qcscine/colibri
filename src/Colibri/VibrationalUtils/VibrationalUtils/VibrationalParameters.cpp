/**
 * @file VibrationalParameters.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#include "VibrationalParameters.h"
#include <boost/filesystem.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <iostream>
#include <set>
#include <sstream>
#include <string>

namespace Scine::Colibri {

template<int rank>

using tensor = std::unordered_map<std::array<int, rank>, double, boost::hash<std::array<int, rank>>>;

// This function is probably faulty as the force constants are not added but
// overwritten... Maybe delete the second part
void VibrationalParameters::setParameters(const std::string& path) {
  boost::property_tree::iptree pt;
  boost::property_tree::ini_parser::read_ini(path, pt);
  numModes_ = pt.get<int>("numModes");
  vscfIter_ = pt.get_optional<int>("vscfIter").get_value_or(1);
  nModePotentialOrder_ = pt.get_optional<int>("nModePotentialOrder").get_value_or(1);

  dumpOnlyCoupledModes_ = pt.get_optional<bool>("dumpOnlyCoupledModes").get_value_or(false);

  std::string nMaxStr = pt.get_optional<std::string>("nmax").get_value_or("");
  int nm = 0;
  std::stringstream issNmax(nMaxStr);
  while (issNmax >> nm) {
    nMax_.push_back(nm);
  }
  if (nMax_.empty()) {
    for (int mode = 0; mode < numModes_; mode++) {
      nMax_.push_back(6);
    }
  } else if (nMax_.size() == 1) {
    for (int mode = 1; mode < numModes_; mode++) {
      nMax_.push_back(nMax_[0]);
    }
  }

  auto occVector = pt.get<std::string>("ONVector");
  fcidumpFname_ = pt.get<std::string>("fciDumpName");
  int oc = 0;
  std::stringstream iss(occVector);
  while (iss >> oc) {
    ONVector_.push_back(oc);
  }

  std::string coupledModesIds = pt.get_optional<std::string>("coupledModes").get_value_or("");
  int cM = 0;
  std::stringstream issCoup(coupledModesIds);
  while (issCoup >> cM) {
    coupledModes_.push_back(cM);
  }
  if (coupledModes_.empty()) {
    for (int mode = 0; mode < numModes_; mode++) {
      coupledModes_.push_back(mode);
    }
  }

  coeffsFname_ = pt.get_optional<std::string>("coeffsFname").get_value_or("potAndCoeffs.out");

  vscfEnTol_ = pt.get_optional<double>("vscfEnTol").get_value_or(1.0E-8);
  vscfCoeffTol_ = pt.get_optional<double>("vscfCoeffTol").get_value_or(1.0E-10);
  twoBodyTol_ = pt.get_optional<double>("twoBodyTol").get_value_or(1.0E-8);
  numqp_ = pt.get_optional<int>("numqp").get_value_or(8);
  numGridPoints_.resize(numModes_);
  modeFreqs_.resize(numModes_);
  quantumNums_.resize(numModes_);
  for (int mode = 0; mode < numModes_; mode++) {
    numGridPoints_[mode] = pt.get_optional<int>("Mode" + std::to_string(mode) + ".numGrid").get_value_or(10);
    quantumNums_[mode] = pt.get_optional<double>("Mode" + std::to_string(mode) + ".nquantum").get_value_or(5.0);
    modeFreqs_[mode] = pt.get_optional<double>("Mode" + std::to_string(mode) + ".frequency").get_value_or(1.0);
  }
  std::string forceConst1Tmp = pt.get_optional<std::string>("ForceConstants.firstOrder").get_value_or(" ");
  linForceConst_ = getForceConstants<1>(forceConst1Tmp);
  std::string forceConst2Tmp = pt.get_optional<std::string>("ForceConstants.secondOrder").get_value_or(" ");
  quadraticForceConst_ = getForceConstants<2>(forceConst2Tmp);
  std::string forceConst3Tmp = pt.get_optional<std::string>("ForceConstants.thirdOrder").get_value_or(" ");
  cubicForceConst_ = getForceConstants<3>(forceConst3Tmp);
  std::string forceConst4Tmp = pt.get_optional<std::string>("ForceConstants.fourthOrder").get_value_or(" ");
  quarticForceConst_ = getForceConstants<4>(forceConst4Tmp);
}

void VibrationalParameters::setParametersN(const std::string& path) {
  boost::property_tree::iptree pt;
  boost::property_tree::ini_parser::read_ini(path, pt);
  numModes_ = pt.get<int>("numModes");
  vscfIter_ = pt.get_optional<int>("vscfIter").get_value_or(1);
  nModePotentialOrder_ = pt.get_optional<int>("nModePotentialOrder").get_value_or(2);

  dumpOnlyCoupledModes_ = pt.get_optional<bool>("dumpOnlyCoupledModes").get_value_or(false);

  std::string nMaxStr = pt.get_optional<std::string>("nmax").get_value_or("");
  int nm = 0;
  std::stringstream issNmax(nMaxStr);
  while (issNmax >> nm) {
    nMax_.push_back(nm);
  }
  if (nMax_.empty()) {
    for (int mode = 0; mode < numModes_; mode++) {
      nMax_.push_back(6);
    }
  } else if (nMax_.size() == 1) {
    for (int mode = 1; mode < numModes_; mode++) {
      nMax_.push_back(nMax_[0]);
    }
  }

  taylorExpOrder_ = pt.get_optional<int>("taylorExpOrder").get_value_or(6);
  auto occVector = pt.get<std::string>("ONVector");
  fcidumpFname_ = pt.get_optional<std::string>("fciDumpName").get_value_or("");
  inputPESFname_ = pt.get<std::string>("inputPESFile");
  inputPESType_ = pt.get<std::string>("inputPESType");
  primitiveBasisType_ = pt.get_optional<std::string>("primBasisType").get_value_or("DG");
  int oc = 0;
  std::stringstream iss(occVector);
  while (iss >> oc) {
    ONVector_.push_back(oc);
  }

  std::string coupledModesIds = pt.get_optional<std::string>("coupledModes").get_value_or("");
  int cM = 0;
  std::stringstream issCoup(coupledModesIds);
  while (issCoup >> cM) {
    coupledModes_.push_back(cM);
  }
  if (coupledModes_.empty()) {
    for (int mode = 0; mode < numModes_; mode++) {
      coupledModes_.push_back(mode);
    }
  }

  coeffsFname_ = pt.get_optional<std::string>("coeffsFname").get_value_or("potAndCoeffs.out");

  vscfEnTol_ = pt.get_optional<double>("vscfEnTol").get_value_or(1.0E-8);
  vscfCoeffTol_ = pt.get_optional<double>("vscfCoeffTol").get_value_or(1.0E-10);
  twoBodyTol_ = pt.get_optional<double>("twoBodyTol").get_value_or(1.0E-8);
  numqp_ = pt.get_optional<int>("numqp").get_value_or(8);
  numGridPoints_.resize(numModes_);
  modeFreqs_.resize(numModes_);
  quantumNums_.resize(numModes_);
  for (int mode = 0; mode < numModes_; mode++) {
    numGridPoints_[mode] = pt.get_optional<int>("Mode" + std::to_string(mode) + ".numGrid").get_value_or(10);
    quantumNums_[mode] = pt.get_optional<double>("Mode" + std::to_string(mode) + ".nquantum").get_value_or(5.0);
    modeFreqs_[mode] = pt.get_optional<double>("Mode" + std::to_string(mode) + ".frequency").get_value_or(1.0);
  }
}

void VibrationalParameters::setParametersAlsoVCI(const std::string& path) {
  boost::property_tree::iptree pt;
  boost::property_tree::ini_parser::read_ini(path, pt);
  // General settings
  numModes_ = pt.get<int>("numModes");
  nModePotentialOrder_ = pt.get_optional<int>("nModePotentialOrder").get_value_or(2);

  dumpOnlyCoupledModes_ = pt.get_optional<bool>("dumpOnlyCoupledModes").get_value_or(false);

  std::string nMaxStr = pt.get_optional<std::string>("nmax").get_value_or("");
  int nm = 0;
  std::stringstream issNmax(nMaxStr);
  while (issNmax >> nm) {
    nMax_.push_back(nm);
  }
  if (nMax_.empty()) {
    for (int mode = 0; mode < numModes_; mode++) {
      nMax_.push_back(6);
    }
  } else if (nMax_.size() == 1) {
    for (int mode = 1; mode < numModes_; mode++) {
      nMax_.push_back(nMax_[0]);
    }
  }

  // if calculation is started from FCIDUMP
  taylorExpOrder_ = pt.get_optional<int>("taylorExpOrder").get_value_or(6);
  inputPESFname_ = pt.get_optional<std::string>("inputPESFile").get_value_or("");
  inputPESType_ = pt.get_optional<std::string>("inputPESType").get_value_or("");

  // if calculation is started from precalculated singlepoints
  harmFreqFname_ = pt.get_optional<std::string>("harmFreqs").get_value_or("harmfreqs.txt");
  singlePointsFname_ = pt.get_optional<std::string>("singlePoints").get_value_or("");
  refEnergy_ = pt.get_optional<double>("refEnergy").get_value_or(0);

  // if calculation uses OTF calculations
  xyzFname_ = pt.get_optional<std::string>("XYZFile").get_value_or("");
  hessianFname_ = pt.get_optional<std::string>("HessianFile").get_value_or("");
  otfProgram_ = pt.get_optional<std::string>("program").get_value_or("");
  otfMethodFamily_ = pt.get_optional<std::string>("method_family").get_value_or("");
  if (!otfProgram_.empty()) {
    for (auto& x : otfProgram_) {
      x = std::tolower(x);
    }
    otfProgram_[0] = std::toupper(otfProgram_[0]);

    for (auto& x : otfMethodFamily_) {
      x = std::toupper(x);
    }
  }
  otfMethod_ = pt.get_optional<std::string>("method").get_value_or("");

  // VSCF settings
  primitiveBasisType_ = pt.get_optional<std::string>("primBasisType").get_value_or("DG");
  fcidumpFname_ = pt.get_optional<std::string>("fciDumpName").get_value_or("");
  vscfIter_ = pt.get_optional<int>("vscfIter").get_value_or(1);
  vscfEnTol_ = pt.get_optional<double>("vscfEnTol").get_value_or(1.0E-8);
  vscfCoeffTol_ = pt.get_optional<double>("vscfCoeffTol").get_value_or(1.0E-10);
  twoBodyTol_ = pt.get_optional<double>("twoBodyTol").get_value_or(1.0E-8);
  numqp_ = pt.get_optional<int>("numqp").get_value_or(8);
  auto occVector = pt.get<std::string>("ONVector");
  int oc = 0;
  std::stringstream iss(occVector);
  while (iss >> oc) {
    ONVector_.push_back(oc);
  }

  if (ONVector_.size() != numModes_) {
    std::cout << "WARNING!!! The entered ONV is not equal to the number of modes!" << std::endl;
  }
  coeffsFname_ = pt.get_optional<std::string>("coeffsFname").get_value_or("potAndCoeffs.out");

  // Parameters of all the modes
  numGridPoints_.resize(numModes_);
  modeFreqs_.resize(numModes_);
  quantumNums_.resize(numModes_);
  for (int mode = 0; mode < numModes_; mode++) {
    numGridPoints_[mode] = pt.get_optional<int>("Mode" + std::to_string(mode) + ".numGrid").get_value_or(10);
    quantumNums_[mode] = pt.get_optional<double>("Mode" + std::to_string(mode) + ".nquantum").get_value_or(5.0);
    modeFreqs_[mode] = pt.get_optional<double>("Mode" + std::to_string(mode) + ".frequency").get_value_or(1.0);
  }

  std::string coupledModesIds = pt.get_optional<std::string>("coupledModes").get_value_or("");
  int cM = 0;
  std::stringstream issCoup(coupledModesIds);
  while (issCoup >> cM) {
    coupledModes_.push_back(cM);
  }
  if (coupledModes_.empty()) {
    for (int mode = 0; mode < numModes_; mode++) {
      coupledModes_.push_back(mode);
    }
  }

  // General VCI settings
  vciExModes_ = pt.get_optional<int>("vciExModes").get_value_or(2);
  vciTotEx_ = pt.get_optional<int>("vciTotEx").get_value_or(6);
  vciEVSolver_ = pt.get_optional<std::string>("vciEVSolver").get_value_or("SelfAdjoint");
  vciNumStates_ = pt.get_optional<int>("vciNumStates").get_value_or(numModes_ + 1);

  // State-specific VCI settings
  vciDoStateSpecific_ = pt.get_optional<bool>("vciStateSpecific").get_value_or(false);
  if (vciDoStateSpecific_) {
    std::string occupRef = pt.get_optional<std::string>("vciRefDet").get_value_or("");
    vciRefDet_.clear();
    int ocRef = 0;
    std::stringstream issRef(occupRef);
    while (issRef >> ocRef) {
      vciRefDet_.push_back(ocRef);
    }

    std::string occupMin = pt.get_optional<std::string>("vciOccupMin").get_value_or("");
    vciOccupMin_.clear();
    int ocMin = 0;
    std::stringstream issMin(occupMin);
    while (issMin >> ocMin) {
      vciOccupMin_.push_back(ocMin);
    }

    std::string occupMax = pt.get_optional<std::string>("vciOccupMax").get_value_or("");
    vciOccupMax_.clear();
    int ocMax = 0;
    std::stringstream issMax(occupMax);
    while (issMax >> ocMax) {
      vciOccupMax_.push_back(ocMax);
    }
  }
}

void VibrationalParameters::setParametersOnlyVCI(const std::string& path) {
  boost::property_tree::iptree pt;
  boost::property_tree::ini_parser::read_ini(path, pt);
  // General settings
  numModes_ = pt.get<int>("numModes");
  nModePotentialOrder_ = pt.get_optional<int>("nModePotentialOrder").get_value_or(2);

  std::string nMaxStr = pt.get_optional<std::string>("nmax").get_value_or("");
  int nm = 0;
  std::stringstream issNmax(nMaxStr);
  while (issNmax >> nm) {
    nMax_.push_back(nm);
  }
  if (nMax_.empty()) {
    for (int mode = 0; mode < numModes_; mode++) {
      nMax_.push_back(6);
    }
  } else if (nMax_.size() == 1) {
    for (int mode = 1; mode < numModes_; mode++) {
      nMax_.push_back(nMax_[0]);
    }
  }

  fcidumpFname_ = pt.get_optional<std::string>("fciDumpName").get_value_or("");

  // General VCI settings
  vciExModes_ = pt.get_optional<int>("vciExModes").get_value_or(2);
  vciTotEx_ = pt.get_optional<int>("vciTotEx").get_value_or(6);
  vciEVSolver_ = pt.get_optional<std::string>("vciEVSolver").get_value_or("SelfAdjoint");
  vciNumStates_ = pt.get_optional<int>("vciNumStates").get_value_or(numModes_ + 1);

  // State-specific VCI settings
  vciDoStateSpecific_ = pt.get_optional<bool>("vciStateSpecific").get_value_or(false);
  if (vciDoStateSpecific_) {
    std::string occupRef = pt.get_optional<std::string>("vciRefDet").get_value_or("");
    vciRefDet_.clear();
    int ocRef = 0;
    std::stringstream issRef(occupRef);
    while (issRef >> ocRef) {
      vciRefDet_.push_back(ocRef);
    }

    std::string occupMin = pt.get_optional<std::string>("vciOccupMin").get_value_or("");
    vciOccupMin_.clear();
    int ocMin = 0;
    std::stringstream issMin(occupMin);
    while (issMin >> ocMin) {
      vciOccupMin_.push_back(ocMin);
    }

    std::string occupMax = pt.get_optional<std::string>("vciOccupMax").get_value_or("");
    vciOccupMax_.clear();
    int ocMax = 0;
    std::stringstream issMax(occupMax);
    while (issMax >> ocMax) {
      vciOccupMax_.push_back(ocMax);
    }
  }
}

void VibrationalParameters::setParametersOTF(const std::string& path) {
  boost::property_tree::iptree pt;
  boost::property_tree::ini_parser::read_ini(path, pt);
  // General information
  numModes_ = pt.get<int>("numModes");
  nModePotentialOrder_ = pt.get_optional<int>("nModePotentialOrder").get_value_or(2);

  dumpOnlyCoupledModes_ = pt.get_optional<bool>("dumpOnlyCoupledModes").get_value_or(false);

  doTransitionStateOptimization_ = pt.get_optional<bool>("transitionState").get_value_or(false);

  // OTF settings
  xyzFname_ = pt.get<std::string>("XYZFile");
  refEnergy_ = pt.get_optional<double>("refEnergy").get_value_or(0);
  hessianFname_ = pt.get_optional<std::string>("HessianFile").get_value_or("");
  otfProgram_ = pt.get<std::string>("program");
  for (auto& x : otfProgram_) {
    x = std::tolower(x);
  }
  otfProgram_[0] = std::toupper(otfProgram_[0]);

  otfMethodFamily_ = pt.get<std::string>("method_family");
  for (auto& x : otfMethodFamily_) {
    x = std::toupper(x);
  }
  otfMethod_ = pt.get_optional<std::string>("method").get_value_or("");

  otfBasisSet_ = pt.get_optional<std::string>("basis_set").get_value_or("cc-pvdz");

  stepForNumericalDifferences_ = pt.get_optional<double>("stepForNumericalDifferences").get_value_or(0.002);

  // VSCF settings
  primitiveBasisType_ = pt.get_optional<std::string>("primBasisType").get_value_or("DVR");
  fcidumpFname_ = pt.get_optional<std::string>("fciDumpName").get_value_or("");
  vscfIter_ = pt.get_optional<int>("vscfIter").get_value_or(1);
  vscfEnTol_ = pt.get_optional<double>("vscfEnTol").get_value_or(1.0E-8);
  vscfCoeffTol_ = pt.get_optional<double>("vscfCoeffTol").get_value_or(1.0E-10);
  twoBodyTol_ = pt.get_optional<double>("twoBodyTol").get_value_or(1.0E-8);
  numqp_ = pt.get_optional<int>("numqp").get_value_or(8);

  coeffsFname_ = pt.get_optional<std::string>("coeffsFname").get_value_or("potAndCoeffs.out");

  std::string nMaxStr = pt.get_optional<std::string>("nmax").get_value_or("");
  int nm = 0;
  std::stringstream issNmax(nMaxStr);
  while (issNmax >> nm) {
    nMax_.push_back(nm);
  }
  if (nMax_.empty()) {
    for (int mode = 0; mode < numModes_; mode++) {
      nMax_.push_back(6);
    }
  } else if (nMax_.size() == 1) {
    for (int mode = 1; mode < numModes_; mode++) {
      nMax_.push_back(nMax_[0]);
    }
  }

  auto occVector = pt.get<std::string>("ONVector");
  int oc = 0;
  std::stringstream iss(occVector);
  while (iss >> oc) {
    ONVector_.push_back(oc);
  }

  if (ONVector_.size() != numModes_) {
    std::cout << "WARNING!!! The entered ONV is not equal to the number of modes!" << std::endl;
  }

  std::string coupledModesIds = pt.get_optional<std::string>("coupledModes").get_value_or("");
  int cM = 0;
  std::stringstream issCoup(coupledModesIds);
  while (issCoup >> cM) {
    coupledModes_.push_back(cM);
  }
  if (coupledModes_.empty()) {
    for (int mode = 0; mode < numModes_; mode++) {
      coupledModes_.push_back(mode);
    }
  }

  // Parameters of the individual modes
  numGridPoints_.resize(numModes_);
  quantumNums_.resize(numModes_);
  for (int mode = 0; mode < numModes_; mode++) {
    numGridPoints_[mode] = pt.get_optional<int>("Mode" + std::to_string(mode) + ".numGrid").get_value_or(10);
    quantumNums_[mode] = pt.get_optional<double>("Mode" + std::to_string(mode) + ".nquantum").get_value_or(5.0);
  }

  // General VCI settings
  vciExModes_ = pt.get_optional<int>("vciExModes").get_value_or(2);
  vciTotEx_ = pt.get_optional<int>("vciTotEx").get_value_or(6);
  vciEVSolver_ = pt.get_optional<std::string>("vciEVSolver").get_value_or("SelfAdjoint");
  vciNumStates_ = pt.get_optional<int>("vciNumStates").get_value_or(numModes_ + 1);

  // State-specific VCI settings
  vciDoStateSpecific_ = pt.get_optional<bool>("vciStateSpecific").get_value_or(false);
  if (vciDoStateSpecific_) {
    std::string occupRef = pt.get_optional<std::string>("vciRefDet").get_value_or("");
    vciRefDet_.clear();
    int ocRef = 0;
    std::stringstream issRef(occupRef);
    while (issRef >> ocRef) {
      vciRefDet_.push_back(ocRef);
    }

    std::string occupMin = pt.get_optional<std::string>("vciOccupMin").get_value_or("");
    vciOccupMin_.clear();
    int ocMin = 0;
    std::stringstream issMin(occupMin);
    while (issMin >> ocMin) {
      vciOccupMin_.push_back(ocMin);
    }

    std::string occupMax = pt.get_optional<std::string>("vciOccupMax").get_value_or("");
    vciOccupMax_.clear();
    int ocMax = 0;
    std::stringstream issMax(occupMax);
    while (issMax >> ocMax) {
      vciOccupMax_.push_back(ocMax);
    }
  }
}

void VibrationalParameters::setPES(std::string integral_file) {
  if (!boost::filesystem::exists(integral_file)) {
    throw std::runtime_error("integral_file " + integral_file + " does not exist\n");
  }
  std::ifstream orb_file;
  std::string line_string;
  orb_file.open(integral_file.c_str());
  while (getline(orb_file, line_string)) {
    std::istringstream ss(line_string);

    std::string val;
    std::string ind;
    std::string forceConstTmp;

    int num_ind = 0;

    ss >> val;

    while (ss >> ind && ind != "0") {
      forceConstTmp += std::to_string(std::stoi(ind) - 1);
      forceConstTmp += " ";
      num_ind++;
    }
    forceConstTmp += val;
    switch (num_ind) {
      case 1:
        addForceConstants<1>(linForceConst_, forceConstTmp);
        break;
      case 2:
        addForceConstants<2>(quadraticForceConst_, forceConstTmp);
        break;
      case 3:
        addForceConstants<3>(cubicForceConst_, forceConstTmp);
        break;
      case 4:
        addForceConstants<4>(quarticForceConst_, forceConstTmp);
        break;
      default:
        std::cout << "Integral term could not be read in. Abort." << std::endl;
        exit(1);
    }
  }
}

void VibrationalParameters::setPesMc(std::string integral_file) {
  if (!boost::filesystem::exists(integral_file)) {
    throw std::runtime_error("integral_file " + integral_file + " does not exist\n");
  }
  std::ifstream orb_file;
  std::string line_string;
  bool warn = true;
  orb_file.open(integral_file.c_str());
  while (getline(orb_file, line_string)) {
    std::istringstream ss(line_string);

    std::string val;
    std::string ind;
    std::string forceConstTmp;

    ss >> val;
    std::set<int> modes;

    while (ss >> ind && ind != "0") {
      forceConstTmp += std::to_string(std::stoi(ind) - 1);
      forceConstTmp += " ";
      modes.insert(std::stoi(ind) - 1);
    }
    forceConstTmp += val;
    int mcorder = modes.size();

    if (mcorder <= nModePotentialOrder_) {
      switch (mcorder) {
        case 1:
          addForceConstants<1>(oneModeCoeffs_, forceConstTmp);
          break;
        case 2:
          addForceConstants<2>(twoModeCoeffs_, forceConstTmp);
          break;
        case 3:
          addForceConstants<3>(threeModeCoeffs_, forceConstTmp);
          break;
        case 4:
          addForceConstants<4>(fourModeCoeffs_, forceConstTmp);
          break;
        case 5:
          addForceConstants<5>(fiveModeCoeffs_, forceConstTmp);
          break;
        case 6:
          addForceConstants<6>(sixModeCoeffs_, forceConstTmp);
          break;
        default:
          std::cout << "Integral term could not be read in. Abort." << std::endl;
          exit(1);
      }
    } else if (mcorder > nModePotentialOrder_ && warn) {
      std::cout << "Integral file seems to contain higher-order couplings than "
                   "what the nModePotentialOder is set to."
                << std::endl;
      std::cout << "WARNING! These higher-order coupling terms will be disregarded!" << std::endl;
      warn = false;
    }
  }
}

void VibrationalParameters::setPesMCSP(std::string integral_file) {
  if (!boost::filesystem::exists(integral_file)) {
    throw std::runtime_error("integral_file " + integral_file + " does not exist\n");
  }
  std::ifstream orb_file;
  std::string line_string;
  bool warn = true;
  orb_file.open(integral_file.c_str());
  while (getline(orb_file, line_string)) {
    std::istringstream ss(line_string);

    std::string val;
    std::string ind;
    std::string forceConstTmp;

    ss >> val;
    std::set<int> modes;

    while (ss >> ind && ind != "0") {
      forceConstTmp += std::to_string(std::stoi(ind) - 1);
      forceConstTmp += " ";
      modes.insert(std::stoi(ind) - 1);
    }
    forceConstTmp += val;
    int mcorder = modes.size();

    if (mcorder <= nModePotentialOrder_) {
      switch (mcorder) {
        case 1:
          addForceConstants<1>(oneModeTerms_, forceConstTmp);
          break;
        case 2:
          addForceConstants<2>(twoModeTerms_, forceConstTmp);
          break;
        case 3:
          addForceConstants<3>(threeModeTerms_, forceConstTmp);
          addSet3B(forceConstTmp);
          break;
        case 4:
          addForceConstants<4>(fourModeTerms_, forceConstTmp);
          break;
        case 5:
          addForceConstants<5>(fiveModeTerms_, forceConstTmp);
          break;
        case 6:
          addForceConstants<6>(sixModeTerms_, forceConstTmp);
          break;
        default:
          std::cout << "Integral term could not be read in. Abort." << std::endl;
          exit(1);
      }
    } else if (mcorder > nModePotentialOrder_ && warn) {
      std::cout << "Integral file seems to contain higher-order couplings than "
                   "what the nModePotentialOder is set to."
                << std::endl;
      std::cout << "WARNING! These higher-order coupling terms will be disregarded!" << std::endl;
      warn = false;
    }
  }
}

template<std::size_t rank>
void VibrationalParameters::addForceConstants(tensor<rank>& K, const std::string& inp) {
  std::istringstream iss(inp);
  std::vector<std::string> parsed{std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>()};
  for (std::size_t k = 0; k < parsed.size(); k += (rank + 1)) {
    std::array<int, rank> inds;
    for (std::size_t i = 0; i < rank; i++) {
      inds[i] = std::stoi(parsed[k + i]);
    }
    K[inds] = std::stod(parsed[k + rank]);
  }
}

template<std::size_t mcorder>
void VibrationalParameters::addForceConstants(tensorMC<mcorder>& K, const std::string& inp) {
  std::istringstream iss(inp);
  std::array<std::pair<int, int>, mcorder> inds;
  std::set<int> modes;
  std::multiset<int> modesOcc;
  std::vector<std::string> parsed{std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>()};
  for (std::size_t k = 0; k < parsed.size() - 1; k++) {
    modes.insert(std::stoi(parsed[k]));
    modesOcc.insert(std::stoi(parsed[k]));
  }
  int order = 0;
  auto it = modes.begin();
  while (it != modes.end()) {
    inds[order] = std::make_pair(*it, modesOcc.count(*it));
    order++;
    it++;
  }
  K[inds] = std::stod(parsed.back());
}

// ALB: 3B
template<std::size_t mcorder>
void VibrationalParameters::addForceConstants(tensorMCSP<mcorder>& K, const std::string& inp) {
  std::istringstream iss(inp);
  std::set<int> modes;
  std::multiset<int> modesOcc;
  std::vector<std::string> parsed{std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>()};

  for (std::size_t k = 0; k < parsed.size() - 1; k++) {
    modes.insert(std::stoi(parsed[k]));
    modesOcc.insert(std::stoi(parsed[k]));
  }

  std::array<std::pair<int, int>, mcorder> inds_exps;
  auto it = modes.begin();
  int order = 0;
  while (it != modes.end()) {
    inds_exps[order] = std::make_pair(*it, modesOcc.count(*it));
    order++;
    it++;
  }

  std::sort(inds_exps.begin(), inds_exps.end(),
            [](std::pair<int, int> i, std::pair<int, int> j) { return i.first < j.first; });

  std::array<int, mcorder> inds;
  std::array<int, mcorder> exps;
  for (int i = 0; i < inds_exps.size(); i++) {
    inds[i] = inds_exps[i].first;
    exps[i] = inds_exps[i].second;
  }

  bool foundAlready = false;
  auto find = K.find(inds);
  if (find != K.end()) {
    for (auto it = find->second.begin(); it != find->second.end(); it++) {
      if (it->first == exps) {
        std::cout << "WARNING! Integral file provided multiple lines for the "
                     "same mode coupling term."
                  << std::endl;
        std::cout << "Did not add the following line to the PES: " << inp << std::endl;
        foundAlready = true;
        break;
      }
    }
    if (!foundAlready) {
      find->second.push_back(std::make_pair(exps, std::stod(parsed.back())));
    }
  } else {
    std::vector<std::pair<std::array<int, mcorder>, double>> new_entry{std::make_pair(exps, std::stod(parsed.back()))};
    K[inds] = new_entry;
  }
}

void VibrationalParameters::addSet3B(const std::string& inp) {
  std::istringstream iss(inp);
  std::set<int> modes;
  std::vector<std::string> parsed{std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>()};

  for (std::size_t k = 0; k < parsed.size() - 1; k++) {
    modes.insert(std::stoi(parsed[k]));
  }

  std::array<int, 3> inds;
  auto it = modes.begin();
  int order = 0;
  while (it != modes.end()) {
    inds[order] = *it;
    order++;
    it++;
  }

  std::sort(inds.begin(), inds.end());
  set3B_.insert(inds);
}

template<std::size_t rank>
tensor<rank> VibrationalParameters::getForceConstants(const std::string& inp) {
  tensor<rank> K;
  std::istringstream iss(inp);
  std::vector<std::string> parsed{std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>()};
  for (std::size_t k = 0; k < parsed.size(); k += (rank + 1)) {
    std::array<int, rank> inds;
    for (std::size_t i = 0; i < rank; i++) {
      inds[i] = std::stoi(parsed[k + i]);
    }
    K[inds] = std::stod(parsed[k + rank]);
  }
  return K;
}

void VibrationalParameters::printParameters() {
  std::cout << "PARAMETERS OF CALCULATION" << std::endl;
  std::cout << "Number of modes: " << numModes_ << std::endl;
  std::cout << "Number of SCF iterations: " << vscfIter_ << std::endl;
  std::cout << "Order of potential: " << nModePotentialOrder_ << std::endl;

  std::cout << "Occupation number vector: " << std::endl;
  for (int i = 0; i < numModes_; i++) {
    std::cout << ONVector_[i];
  }
  std::cout << std::endl << std::endl;

  std::cout << "Coupled modes: " << std::endl;
  for (int coupledMode : coupledModes_) {
    std::cout << coupledMode << " ";
  }
  std::cout << std::endl << std::endl;

  if (dumpOnlyCoupledModes_) {
    std::cout << "Only dumpling coupled modes!" << std::endl;
  }

  std::cout << "Number of modals to dump: " << std::endl;
  for (int i = 0; i < numModes_; i++) {
    std::cout << nMax_[i] << " ";
  }
  std::cout << std::endl << std::endl;

  std::cout << "FCIDUMP file: " << fcidumpFname_ << std::endl;
  std::cout << "PotAndCoeffs file: " << coeffsFname_ << std::endl;
  std::cout << "SCF tolerance energy: " << vscfEnTol_ << std::endl;
  std::cout << "SCF tolerance eigenvectors: " << vscfCoeffTol_ << std::endl;
  std::cout << "Two-body integral tolerance: " << twoBodyTol_ << std::endl;
  std::cout << "Number of quadrature points: " << numqp_ << std::endl;
}

void VibrationalParameters::printParametersN() {
  std::cout << "PARAMETERS OF CALCULATION" << std::endl << std::endl;

  std::cout << "Number of modes: " << numModes_ << std::endl;
  std::cout << "Order of potential to consider: " << nModePotentialOrder_ << std::endl;
  std::cout << "Expansion order of Taylor Polynomials to consider: " << taylorExpOrder_ << std::endl;
  std::cout << "Input integrals taken from: " << inputPESFname_ << std::endl;
  std::cout << "Two- and more-body integral tolerance: " << twoBodyTol_ << std::endl << std::endl;

  std::cout << "Occupation number vector: " << std::endl;
  for (int i = 0; i < numModes_; i++) {
    std::cout << ONVector_[i] << " ";
  }
  std::cout << std::endl << std::endl;

  std::cout << "Coupled modes: " << std::endl;
  for (int coupledMode : coupledModes_) {
    std::cout << coupledMode << " ";
  }
  std::cout << std::endl << std::endl;

  if (dumpOnlyCoupledModes_) {
    std::cout << "Only dumpling coupled modes!" << std::endl;
  }

  std::cout << "Primitive basis type: " << primitiveBasisType_ << std::endl;
  std::cout << "Number of SCF iterations: " << vscfIter_ << std::endl;
  std::cout << "SCF tolerance energy: " << vscfEnTol_ << std::endl;
  std::cout << "SCF tolerance eigenvectors: " << vscfCoeffTol_ << std::endl;
  std::cout << "Number of quadrature points: " << numqp_ << std::endl << std::endl;

  std::cout << "Print out integrals to FCIDUMP file: " << fcidumpFname_ << std::endl;
  std::cout << "PotAndCoeffs file: " << coeffsFname_ << std::endl;

  std::cout << "Number of modals to dump: " << std::endl;
  for (int i = 0; i < numModes_; i++) {
    std::cout << nMax_[i] << " ";
  }
  std::cout << std::endl << std::endl;
}

void VibrationalParameters::printParametersVSCF() {
  std::cout << "PARAMETERS OF VSCF CALCULATION" << std::endl << std::endl;

  std::cout << "Number of modes: " << numModes_ << std::endl;
  std::cout << "Order of potential to consider: " << nModePotentialOrder_ << std::endl;
  std::cout << "Expansion order of Taylor Polynomials to consider: " << taylorExpOrder_ << std::endl;
  std::cout << "Input integrals taken from: " << inputPESFname_ << std::endl;
  std::cout << "Two- and more-body integral tolerance: " << twoBodyTol_ << std::endl << std::endl;

  std::cout << "Occupation number vector: " << std::endl;
  for (int i = 0; i < numModes_; i++) {
    std::cout << ONVector_[i] << " ";
  }
  std::cout << std::endl << std::endl;

  std::cout << "Coupled modes: " << std::endl;
  for (int coupledMode : coupledModes_) {
    std::cout << coupledMode << " ";
  }
  std::cout << std::endl << std::endl;

  if (dumpOnlyCoupledModes_) {
    std::cout << "Only dumpling coupled modes!" << std::endl;
  }

  std::cout << "Primitive basis type: " << primitiveBasisType_ << std::endl;
  std::cout << "Number of SCF iterations: " << vscfIter_ << std::endl;
  std::cout << "SCF tolerance energy: " << vscfEnTol_ << std::endl;
  std::cout << "SCF tolerance eigenvectors: " << vscfCoeffTol_ << std::endl;
  std::cout << "Number of quadrature points: " << numqp_ << std::endl << std::endl;

  std::cout << "Print out integrals to FCIDUMP file: " << fcidumpFname_ << std::endl;
  std::cout << "PotAndCoeffs file: " << coeffsFname_ << std::endl;

  std::cout << "Number of modals to dump: " << std::endl;
  for (int i = 0; i < numModes_; i++) {
    std::cout << nMax_[i] << " ";
  }
  std::cout << std::endl << std::endl;
}

void VibrationalParameters::printParametersOTF() {
  std::cout << "PARAMETERS OF ON-THE-FLY VSCF CALCULATION" << std::endl << std::endl;

  std::cout << "Number of modes: " << numModes_ << std::endl;
  std::cout << "Order of potential to consider: " << nModePotentialOrder_ << std::endl;
  std::cout << "Starting structure taken from: " << xyzFname_ << std::endl;
  std::cout << "Transition state optimization: " << doTransitionStateOptimization_ << std::endl;
  std::cout << "ES program: " << otfProgram_ << std::endl;
  std::cout << "ES method family: " << otfMethodFamily_ << std::endl;
  std::cout << "ES method: " << otfMethod_ << std::endl;
  std::cout << "ES basis set: " << otfBasisSet_ << std::endl;

  if (otfProgram_ == "Pipnn") {
    std::cout << "Step size for numerical differences: " << stepForNumericalDifferences_ << std::endl;
  }

  std::cout << "Occupation number vector: " << std::endl;
  for (int i = 0; i < numModes_; i++) {
    std::cout << ONVector_[i] << " ";
  }
  std::cout << std::endl << std::endl;

  std::cout << "Coupled modes: " << std::endl;
  for (int coupledMode : coupledModes_) {
    std::cout << coupledMode << " ";
  }
  std::cout << std::endl << std::endl;

  if (dumpOnlyCoupledModes_) {
    std::cout << "Only dumpling coupled modes!" << std::endl;
  }

  std::cout << "Primitive basis type: " << primitiveBasisType_ << std::endl;
  std::cout << "Number of SCF iterations: " << vscfIter_ << std::endl;
  std::cout << "SCF tolerance energy: " << vscfEnTol_ << std::endl;
  std::cout << "SCF tolerance eigenvectors: " << vscfCoeffTol_ << std::endl;
  std::cout << "Number of quadrature points: " << numqp_ << std::endl << std::endl;

  std::cout << "Print out integrals to FCIDUMP file: " << fcidumpFname_ << std::endl;
  std::cout << "PotAndCoeffs file: " << coeffsFname_ << std::endl;

  std::cout << "Number of modals to dump: " << std::endl;
  for (int i = 0; i < numModes_; i++) {
    std::cout << nMax_[i] << " ";
  }
  std::cout << std::endl << std::endl;
}

void VibrationalParameters::printParametersFromSPs() {
  std::cout << "PARAMETERS OF VSCF CALCULATION BASED ON PRECALCULATED SINGLEPOINTS" << std::endl << std::endl;

  std::cout << "Number of modes: " << numModes_ << std::endl;
  std::cout << "Harmonic frequencies taken from: " << harmFreqFname_ << std::endl;
  std::cout << "Singlepoints read from: " << singlePointsFname_ << std::endl;
  std::cout << "Order of potential to consider: " << nModePotentialOrder_ << std::endl;

  std::cout << "Occupation number vector: " << std::endl;
  for (int i = 0; i < numModes_; i++) {
    std::cout << ONVector_[i] << " ";
  }
  std::cout << std::endl << std::endl;

  std::cout << "Coupled modes: " << std::endl;
  for (int coupledMode : coupledModes_) {
    std::cout << coupledMode << " ";
  }
  std::cout << std::endl << std::endl;

  if (dumpOnlyCoupledModes_) {
    std::cout << "Only dumpling coupled modes!" << std::endl;
  }

  std::cout << "Primitive basis type: " << primitiveBasisType_ << std::endl;
  std::cout << "Number of SCF iterations: " << vscfIter_ << std::endl;
  std::cout << "SCF tolerance energy: " << vscfEnTol_ << std::endl;
  std::cout << "SCF tolerance eigenvectors: " << vscfCoeffTol_ << std::endl;
  std::cout << "Number of quadrature points: " << numqp_ << std::endl << std::endl;

  std::cout << "Print out integrals to FCIDUMP file: " << fcidumpFname_ << std::endl;
  std::cout << "PotAndCoeffs file: " << coeffsFname_ << std::endl;

  std::cout << "Number of modals to dump: " << std::endl;
  for (int i = 0; i < numModes_; i++) {
    std::cout << nMax_[i] << " ";
  }
  std::cout << std::endl << std::endl;
}

void VibrationalParameters::printParametersHybrid() {
  std::cout << "PARAMETERS OF VSCF CALCULATION BASED ON PRECALCULATED SINGLEPOINTS" << std::endl << std::endl;

  std::cout << "Number of modes: " << numModes_ << std::endl;
  std::cout << "Reference structure taken from: " << xyzFname_ << std::endl;
  std::cout << "Harmonic frequencies taken from: " << harmFreqFname_ << std::endl;
  std::cout << "Singlepoints read from: " << singlePointsFname_ << std::endl;
  std::cout << "Order of potential to consider: " << nModePotentialOrder_ << std::endl;

  std::cout << "ES program: " << otfProgram_ << std::endl;
  std::cout << "ES method family: " << otfMethodFamily_ << std::endl;
  std::cout << "ES method: " << otfMethod_ << std::endl;
  std::cout << "ES basis set: " << otfBasisSet_ << std::endl;

  std::cout << "Occupation number vector: " << std::endl;
  for (int i = 0; i < numModes_; i++) {
    std::cout << ONVector_[i] << " ";
  }
  std::cout << std::endl << std::endl;

  std::cout << "Coupled modes: " << std::endl;
  for (int coupledMode : coupledModes_) {
    std::cout << coupledMode << " ";
  }
  std::cout << std::endl << std::endl;

  if (dumpOnlyCoupledModes_) {
    std::cout << "Only dumpling coupled modes!" << std::endl;
  }

  std::cout << "Primitive basis type: " << primitiveBasisType_ << std::endl;
  std::cout << "Number of SCF iterations: " << vscfIter_ << std::endl;
  std::cout << "SCF tolerance energy: " << vscfEnTol_ << std::endl;
  std::cout << "SCF tolerance eigenvectors: " << vscfCoeffTol_ << std::endl;
  std::cout << "Number of quadrature points: " << numqp_ << std::endl << std::endl;

  std::cout << "Print out integrals to FCIDUMP file: " << fcidumpFname_ << std::endl;
  std::cout << "PotAndCoeffs file: " << coeffsFname_ << std::endl;

  std::cout << "Number of modals to dump: " << std::endl;
  for (int i = 0; i < numModes_; i++) {
    std::cout << nMax_[i] << " ";
  }
  std::cout << std::endl << std::endl;
}

void VibrationalParameters::printParametersVCI() {
  std::cout << "PARAMETERS OF VCI CALCULATION" << std::endl << std::endl;

  std::cout << "Number of modes that can simultaneously be excited: " << vciExModes_ << std::endl;
  if (!vciDoStateSpecific_) {
    std::cout << "Maximum excitation degree of each mode: ";
    for (int i = 0; i < numModes_; i++) {
      std::cout << nMax_[i] << " ";
    }
    std::cout << std::endl << std::endl;
  }
  std::cout << "Maximum total excitation degree: " << vciTotEx_ << std::endl;
  std::cout << "Solver used for eigenvalue problem: " << vciEVSolver_ << std::endl;
  std::cout << "Number of states of interest: " << vciNumStates_ << std::endl;

  if (vciDoStateSpecific_) {
    std::cout << std::endl << "Settings for state-specific VCI calculation:" << std::endl;
    if (vciRefDet_.size() != numModes_) {
      std::cout << std::endl
                << "WARNING!!! The entered reference determinant length is not "
                   "equal to the number of modes!"
                << std::endl;
      std::cout << "Warning: The reference will be set to zero for all modes" << std::endl;
      vciRefDet_.clear();
      std::vector<int> vect(numModes_, 0);
      vciRefDet_ = vect;
    }

    if (vciOccupMin_.size() != numModes_) {
      std::cout << "Warning: The minimum occupation will be set to zero for all modes" << std::endl;
      vciOccupMin_.clear();
      std::vector<int> vect(numModes_, 0);
      vciOccupMin_ = vect;
    }

    if (vciOccupMax_.size() != numModes_) {
      std::cout << "Warning: The maximum occupation will be set to nMax-1 for "
                   "all modes"
                << std::endl;
      vciOccupMax_.clear();
      vciOccupMax_.resize(numModes_);
      for (int i = 0; i < numModes_; i++) {
        vciOccupMax_[i] = nMax_[i] - 1;
      }
    }

    bool changed = false;
    for (int i = 0; i < numModes_; i++) {
      if (vciOccupMax_[i] < vciOccupMin_[i]) {
        changed = true;
        vciOccupMax_[i] = nMax_[i] - 1;
        vciOccupMin_[i] = 0;
      }
    }
    if (changed) {
      std::cout << std::endl
                << "WARNING!!! The minimum and maximum occupations are "
                   "incompatible!!!"
                << std::endl;
      std::cout << "Warning! They have been modified accordingly." << std::endl;
    }

    std::cout << "Occupation number vector of reference det: ";
    for (int i = 0; i < numModes_; i++) {
      std::cout << vciRefDet_[i];
    }
    std::cout << std::endl;
    std::cout << "Minimum occupation number vector: ";
    for (int i = 0; i < numModes_; i++) {
      std::cout << vciOccupMin_[i];
    }
    std::cout << std::endl;
    std::cout << "Maximum occupation number vector ";
    for (int i = 0; i < numModes_; i++) {
      std::cout << vciOccupMax_[i];
    }
    std::cout << std::endl;
  }
}

void VibrationalParameters::printParametersOnlyVCI() {
  std::cout << "PARAMETERS OF VCI CALCULATION" << std::endl << std::endl;

  std::cout << "Integrals read from FCIDUMP file: " << fcidumpFname_ << std::endl;

  std::cout << "Number of modes that can simultaneously be excited: " << vciExModes_ << std::endl;
  if (!vciDoStateSpecific_) {
    std::cout << "Maximum excitation degree of each mode: ";
    for (int i = 0; i < numModes_; i++) {
      std::cout << nMax_[i] << " ";
    }
    std::cout << std::endl << std::endl;
  }
  std::cout << "Maximum total excitation degree: " << vciTotEx_ << std::endl;
  std::cout << "Solver used for eigenvalue problem: " << vciEVSolver_ << std::endl;
  std::cout << "Number of states of interest: " << vciNumStates_ << std::endl;

  if (vciDoStateSpecific_) {
    std::cout << std::endl << "Settings for state-specific VCI calculation:" << std::endl;
    if (vciRefDet_.size() != numModes_) {
      std::cout << std::endl
                << "WARNING!!! The entered reference determinant length is not "
                   "equal to the number of modes!"
                << std::endl;
      std::cout << "Warning: The reference will be set to zero for all modes" << std::endl;
      vciRefDet_.clear();
      std::vector<int> vect(numModes_, 0);
      vciRefDet_ = vect;
    }

    if (vciOccupMin_.size() != numModes_) {
      std::cout << "Warning: The minimum occupation will be set to zero for all modes" << std::endl;
      vciOccupMin_.clear();
      std::vector<int> vect(numModes_, 0);
      vciOccupMin_ = vect;
    }

    if (vciOccupMax_.size() != numModes_) {
      std::cout << "Warning: The maximum occupation will be set to nMax-1 for "
                   "all modes"
                << std::endl;
      vciOccupMax_.clear();
      vciOccupMax_.resize(numModes_);
      for (int i = 0; i < numModes_; i++) {
        vciOccupMax_[i] = nMax_[i] - 1;
      }
    }

    bool changed = false;
    for (int i = 0; i < numModes_; i++) {
      if (vciOccupMax_[i] < vciOccupMin_[i]) {
        changed = true;
        vciOccupMax_[i] = nMax_[i] - 1;
        vciOccupMin_[i] = 0;
      }
    }
    if (changed) {
      std::cout << std::endl
                << "WARNING!!! The minimum and maximum occupations are "
                   "incompatible!!!"
                << std::endl;
      std::cout << "Warning! They have been modified accordingly." << std::endl;
    }

    std::cout << "Occupation number vector of reference det: ";
    for (int i = 0; i < numModes_; i++) {
      std::cout << vciRefDet_[i];
    }
    std::cout << std::endl;
    std::cout << "Minimum occupation number vector: ";
    for (int i = 0; i < numModes_; i++) {
      std::cout << vciOccupMin_[i];
    }
    std::cout << std::endl;
    std::cout << "Maximum occupation number vector ";
    for (int i = 0; i < numModes_; i++) {
      std::cout << vciOccupMax_[i];
    }
    std::cout << std::endl;
  }
}

} // namespace Scine::Colibri
