/**
 * @file VibrationalCI.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "VibrationalCI.h"
#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace Scine::Colibri {

VibrationalCI::VibrationalCI(const VibrationalParameters& vibParms) {
  nModes_ = vibParms.numModes_;
  potentialOrder_ = vibParms.nModePotentialOrder_;
  exModes_ = vibParms.vciExModes_;
  totEx_ = vibParms.vciTotEx_;
  nMax_ = vibParms.nMax_;
  doStateSpecific_ = vibParms.vciDoStateSpecific_;
  refDet_ = vibParms.vciRefDet_;
  occupMin_ = vibParms.vciOccupMin_;
  occupMax_ = vibParms.vciOccupMax_;
  evSolver_ = vibParms.vciEVSolver_;
  numStatesToCalculate_ = vibParms.vciNumStates_;
  fcidumpFile_ = vibParms.fcidumpFname_;
}

double VibrationalCI::CI() {
  std::cout << std::endl << "----- Started VCI -----" << std::endl << std::endl;
  std::cout << "VCI: Generating states";
  if (doStateSpecific_) {
    numStates_ = this->generateStatesFromDet();
  } else {
    numStates_ = this->generateStatesFromGS();
  }
  std::cout << "  -->  generated " << numStates_ << std::endl;

  std::cout << "VCI: Setting up the Hamiltonian matrix";
  try {
    hamiltonianMatrix_.resize(numStates_, numStates_);
  }
  catch (const std::bad_alloc&) {
    std::cout << "Running out of memory for Hamiltonian matrix." << std::endl;
    exit(1);
  }
  double h_00 = this->fillHamiltonianMat();
  if (evSolver_ == "Davidson") {
    sortStatesWithIncreasingEnergy();
    h_00 = this->fillHamiltonianMat(); // Matrix needs to be recomputed with new
                                       // order
  }
  std::cout << "  -->  H_00 is " << h_00 << std::endl;

  // printStates();

  std::cout << "VCI: Computing EVs with " << evSolver_ << std::endl;
  if (evSolver_ == "SelfAdjoint") {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(hamiltonianMatrix_);
    eVals_ = es.eigenvalues();
    eVecs_ = es.eigenvectors();
  } else if (evSolver_ == "Davidson") {
    this->doDavidson();
  } else {
    std::cout << std::endl << "ERROR!!! This eigensolver is not known!" << std::endl;
    exit(1);
  }

  // Print relvant results
  printEValsMinusGs();
  printEVecsLargestCoeffAndDet();

  return eVals_[0];
}

std::vector<std::pair<double, Eigen::VectorXd>> VibrationalCI::sortEigenPairs(const Eigen::MatrixXd& eVecs,
                                                                              const Eigen::VectorXd& eVals) {
  int nBasis = eVals.size();
  std::vector<std::pair<double, Eigen::VectorXd>> evalPairs(nBasis);
  for (int i = 0; i < nBasis; i++) {
    evalPairs[i].first = eVals(i);
    evalPairs[i].second = eVecs.col(i);
  }
  std::sort(evalPairs.begin(), evalPairs.end(), [](auto a, auto b) { return a.first < b.first; });
  return evalPairs;
}

void VibrationalCI::printEValsMinusGs() {
  std::cout << std::endl << "The VCI energies are:" << std::endl;
  std::cout << "GS energy: " << std::fixed << eVals_[0] * 219474.63 << " cm^-1" << std::endl;
  std::cout << "Ex energies relative to GS:" << std::endl;
  Eigen::VectorXd exEnergies = eVals_.array() - eVals_[0];
  for (int i = 1; (i < numStatesToCalculate_ && i < eVals_.size()); i++) {
    std::cout << i << " " << exEnergies[i] * 219474.63 << " cm^-1" << std::endl;
  }
  std::cout << std::endl;
}

void VibrationalCI::printEVecs() {
  std::cout << "GS coeffs: " << std::endl << eVecs_.col(0) << std::endl;
  for (int i = 1; (i < numStatesToCalculate_ && i < eVecs_.size()); i++) {
    std::cout << "First 20 CI coeffs for ex state " << i << std::endl << eVecs_.col(i).head(20) << std::endl;
  }
}

void VibrationalCI::printEVecsLargestCoeffAndDet() {
  std::cout << "The dominating VCI coefficients are:" << std::endl;
  // to get location of maximum coefficient
  int maxRow = 0;
  double max = eVecs_.col(0).array().abs().maxCoeff(&maxRow);
  std::cout << "GS max coeff is " << max << " for det ";
  for (std::vector<int>::iterator it = states_[maxRow].begin(); it != states_[maxRow].end(); it++) {
    std::cout << *it;
  }
  std::cout << std::endl;
  for (int i = 1; (i < numStatesToCalculate_ && i < eVecs_.size()); i++) {
    max = eVecs_.col(i).array().abs().maxCoeff(&maxRow);
    std::cout << "Ex state " << i << " max coeff is " << max << " for det ";
    for (std::vector<int>::iterator it = states_[maxRow].begin(); it != states_[maxRow].end(); it++) {
      std::cout << *it;
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void VibrationalCI::printStates() {
  std::cout << "States: " << std::endl;
  for (const auto& state : states_) {
    for (auto it2 : state) {
      std::cout << it2;
    }
    std::cout << std::endl;
  }
  std::cout << "End states." << std::endl;
}

void VibrationalCI::sortStatesWithIncreasingEnergy() {
  std::vector<std::pair<double, std::vector<int>>> energyOfEachState(numStates_);
  for (int i = 0; i < numStates_; i++) {
    energyOfEachState[i].first = hamiltonianMatrix_(i, i);
    energyOfEachState[i].second = states_[i];
  }
  std::sort(energyOfEachState.begin(), energyOfEachState.end(), [](auto a, auto b) { return a.first < b.first; });

  for (int i = 0; i < numStates_; i++) {
    states_[i] = energyOfEachState[i].second;
  }
}

} // namespace Scine::Colibri
