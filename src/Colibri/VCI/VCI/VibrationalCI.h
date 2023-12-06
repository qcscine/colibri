/**
 * @file VibrationalCI.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef VIB_CI_H
#define VIB_CI_H

#include <VibrationalUtils/VibrationalParameters.h>
#include <Eigen/Dense>
#include <iostream>

namespace Scine {
namespace Colibri {
class VibrationalCI {
  /**
   * @brief Base class for vibrational CI-implementations.
   * @param nModes Number of vibrational modes.
   * @param potentialOrder 1: only one-body terms, 2: up to two-body terms, etc.
   * @param refDet ONV to start excitations from.
   * @param occupMin minimum occupation, no determinants generated below this
   * one.
   * @param occupMax maximum occupation, no determinants generated above this
   * one.
   * @param exModes number of modes that can be excited simulatenously
   * @param totEx total number of excitations that can simultaneously occur
   * @param numStatesToCalculate number of states that the user is interested in
   * (only those solved in Davidson)
   * @param nMax maximum excitations of each mode
   * @param numStates number of states generated --> size of CI-space
   * @param doStateSpecific if false, start from ground state, if true start
   * from refDet
   * @param states vector with all states of the CI-space
   * @param evSolver can be either SelfAdjoint for exact diagonalisation or
   * Davidson
   * @param fcidumpFile fcidump file name from which integrals are read if only
   * VCI is done
   */

 public:
  VibrationalCI(const VibrationalParameters& vibParms);
  virtual ~VibrationalCI() = default;

  double CI();

 protected:
  int nModes_, potentialOrder_, exModes_, totEx_, numStates_, numStatesToCalculate_;
  bool doStateSpecific_;
  std::string fcidumpFile_;
  std::vector<int> refDet_, nMax_, occupMin_, occupMax_;
  std::vector<std::vector<int>> states_;
  Eigen::MatrixXd hamiltonianMatrix_;
  Eigen::VectorXd eVals_;
  Eigen::MatrixXd eVecs_;

 private:
  std::string evSolver_;

  virtual int generateStatesFromGS() = 0;
  virtual int generateStatesFromDet() {
    return 0;
  };
  virtual double fillHamiltonianMat() = 0;
  virtual double energyOfState(std::vector<int> /*state*/) {
    return 0;
  };
  virtual void doDavidson(){};

  void printEValsMinusGs();
  void printEVecs();
  void printEVecsLargestCoeffAndDet();
  void printStates();

  void sortStatesWithIncreasingEnergy();

  static std::vector<std::pair<double, Eigen::VectorXd>> sortEigenPairs(const Eigen::MatrixXd& eVecs,
                                                                        const Eigen::VectorXd& eVals);
};
} // namespace Colibri
} // namespace Scine

#endif
