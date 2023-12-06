/**
 * @file VSCF.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef VSCF_H
#define VSCF_H
#include "VibrationalSCF.h"
#include <VibrationalUtils/IntegralDumper.h>
#include <VibrationalUtils/PESLibrary/PES.h>
#include <VibrationalUtils/Tensor.h>
#include <VibrationalUtils/VibrationalBases/ModalBasis.h>
#include <VibrationalUtils/VibrationalParameters.h>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <memory>
#include <vector>

namespace Scine {
namespace Colibri {
class VSCF : public VibrationalSCF, public IntegralDumper {
  /**
   * @brief Implementation of the VSCF algorithm.
   * Stores the integrals which are used in the base-class SCF.
   */
 public:
  using Base = VibrationalSCF;
  VSCF(const VibrationalParameters& parms, std::shared_ptr<ModalBasis> basis);
  ~VSCF() override = default;

  /*Matrix representation of the operator in basis of mode mode*/
  Eigen::MatrixXd getMeanFieldOperator(int mode, const std::vector<Eigen::VectorXd>& coefficients) override;
  /*Matrix representation of the operator in basis of mode mode*/
  Eigen::MatrixXd getMeanFieldOperator3(int mode, const std::vector<Eigen::VectorXd>& coefficients) override;
  /*Matrix representation of the kinetic operator in basis of mode mode*/
  Eigen::MatrixXd getKineticOperator(int mode) override;
  /*Matrix representation of the harmonic operator in basis of mode mode*/
  Eigen::MatrixXd getMeanFieldOperator(int mode) override;
  /*Overlap of the basis object of mode mode*/
  Eigen::MatrixXd getOverlap(int mode) override;
  /*Vector with size of the number of modes, each containing the initial
   * coefficients of the modals*/
  std::vector<Eigen::VectorXd> getStartingGuess(bool randomize) override;
  /*Contracts coefficients and writes dump-file*/
  void dumpIntegrals(const std::vector<std::vector<std::pair<double, Eigen::VectorXd>>>& evalPairs) override;
  /*Returns VSCF-energy*/
  double getEnergy(const std::vector<Eigen::VectorXd>& coeffs, const std::vector<double>& evals) override;
  /*Returns VSCF-energy*/
  double getEnergy3(const std::vector<Eigen::VectorXd>& coeffs, const std::vector<double>& evals) override;
  /*One-body VSCF energy*/
  double getEnergy(const std::vector<double>& evals) override;
  /*Print all modal energies*/
  void printAllEnergies(const std::vector<Eigen::VectorXd>& coeffs, std::vector<Eigen::VectorXd>& allEvals) override;
  /*Returns the number of basis functions for mode mode*/
  int getBasisSize(int mode) override;

  /*Return the Modal parameters needed for a subsequent VCI calc*/
  tensor<3> getOneBodyModals() override;
  tensor<6> getTwoBodyModals() override;
  tensor<9> getThreeBodyModals() override;

 private:
  bool storeIntegrals();
  void storeHarmonicIntegrals();
  void storeOneBodyIntegrals();
  void storeTwoBodyIntegrals();
  void storeThreeBodyIntegrals();
  Eigen::MatrixXd contractTwoBodyIntegral(int mode, const std::vector<Eigen::VectorXd>& C);
  Eigen::MatrixXd contractThreeBodyIntegral(int mode, const std::vector<Eigen::VectorXd>& C);
  std::pair<int, int> getTwoBodyIndex(int i, int j, int k, int l, int modei, int modej) const;
  std::pair<int, int> getThreeBodyIndex(int i, int j, int k, int l, int m, int n, int mode1, int mode2, int mode3) const;
  int getDensityIndex(int k, int l, int modej) const;
  Eigen::VectorXd ConstructDensity(int mode, const Eigen::VectorXd& coeffLeft, const Eigen::VectorXd& coeffRight);
  Eigen::VectorXd ConstructDensity3(int mode2, const Eigen::VectorXd& coeffLeft2, const Eigen::VectorXd& coeffRight2,
                                    int mode3, const Eigen::VectorXd& coeffLeft3, const Eigen::VectorXd& coeffRight3);
  std::shared_ptr<ModalBasis> basis_;
  std::vector<Eigen::MatrixXd> oneBodyIntegrals_, overlaps_, oneBodyHarmonicIntegrals_;
  std::unordered_map<std::pair<int, int>, Eigen::SparseMatrix<double>, boost::hash<std::pair<int, int>>> twoBodyIntegrals_;

  std::unordered_map<std::tuple<int, int, int>, Eigen::SparseMatrix<double>, boost::hash<std::tuple<int, int, int>>> threeBodyIntegrals_;

  tensor<3> oneBodyModals_;
  tensor<6> twoBodyModals_;
  tensor<9> threeBodyModals_;

  setNbody<3> set3B_;
};

} // namespace Colibri
} // namespace Scine

#endif
