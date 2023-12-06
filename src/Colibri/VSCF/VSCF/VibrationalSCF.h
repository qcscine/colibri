/**
 * @file VibrationalSCF.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef VIB_SCF_H
#define VIB_SCF_H

#include <VibrationalUtils/VibrationalParameters.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

namespace Scine {
namespace Colibri {
class VibrationalSCF {
  /**
   * @brief Base class for vibrational SCF-implementations
   * for calculating the modal-coefficients.
   * @param nModes Number of vibrational modes.
   * @param maxIter Maximum number of SCF iterations.
   * @param potentialOrder 0: harmonic, 1: only one-body terms, 2: up to
   * two-body terms.
   * @param occupation ONV of targeted state in SCF.
   * @param integralTol tolerance below which two-body integrals are neglected.
   * @param nMax number of modals of each mode written in the integral dump.
   * @param isConverged bool for checking whether SCF converged.
   */

 public:
  VibrationalSCF(const VibrationalParameters& vibParms);
  virtual ~VibrationalSCF() = default;

  double SCF();
  virtual Eigen::MatrixXd getMeanFieldOperator(int mode, const std::vector<Eigen::VectorXd>& coeffs) = 0;
  virtual Eigen::MatrixXd getMeanFieldOperator3(int mode, const std::vector<Eigen::VectorXd>& coeffs){};
  virtual Eigen::MatrixXd getKineticOperator(int mode){};
  virtual Eigen::MatrixXd getMeanFieldOperator(int mode) = 0;
  virtual Eigen::MatrixXd getOverlap(int mode) = 0;
  virtual std::vector<Eigen::VectorXd> getStartingGuess(bool randomize = false) = 0;
  virtual void dumpIntegrals(const std::vector<std::vector<std::pair<double, Eigen::VectorXd>>>& evalPairs) = 0;
  virtual double getEnergy(const std::vector<Eigen::VectorXd>& coeffs, const std::vector<double>& evals) = 0;
  virtual double getEnergy3(const std::vector<Eigen::VectorXd>& coeffs, const std::vector<double>& evals){};
  virtual double getEnergy(const std::vector<double>& evals) = 0;
  virtual void printAllEnergies(const std::vector<Eigen::VectorXd>& coeffs, std::vector<Eigen::VectorXd>& allEvals){};
  virtual int getBasisSize(int mode) = 0;
  virtual tensor<3> getOneBodyModals(){};
  virtual tensor<6> getTwoBodyModals(){};
  virtual tensor<9> getThreeBodyModals(){};

 protected:
  int nModes, maxIter, potentialOrder;
  bool dumpOnlyCoupledModes_;
  std::string integralDumpFile, coeffsDumpFile_;
  std::vector<int> occupation, nMax_;
  std::vector<int> selectedModes;
  std::vector<int> coupledModes_;
  double integralTol_;

 private:
  bool isConverged_;
  double scfEnergyTolerance_, scfCoeffTolerance_;
  std::vector<Eigen::VectorXd> coeffs_;
  std::vector<double> evals_;
  // the sorted eigenvalues - coefficients of the modals of each mode are stored
  std::vector<std::vector<std::pair<double, Eigen::VectorXd>>> allModalPairs_;
  std::vector<std::vector<Eigen::VectorXd>> evecsIterations_;
  std::vector<double> evalsIterations_;

  static std::vector<std::pair<double, Eigen::VectorXd>> sortEigenPairs(const Eigen::MatrixXd& eVecs,
                                                                        const Eigen::VectorXd& eVals);
  bool evecConverged(std::vector<Eigen::VectorXd>& thisIter, std::vector<Eigen::VectorXd>& prevIter, double crit);
};
} // namespace Colibri
} // namespace Scine

#endif
