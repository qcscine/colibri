/**
 * @file DistributedGaussians.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#ifndef DISTR_GAUSS_H
#define DISTR_GAUSS_H
#include "ModalBasis.h"
#include "ModeGrid.h"
#include <Utils/Technical/CloneInterface.h>
#include <VibrationalUtils/GaussHermite.h>
#include <VibrationalUtils/PESLibrary/PES.h>
#include <VibrationalUtils/VibrationalParameters.h>
#include <cmath>
#include <memory>
#include <vector>

namespace Scine {
namespace Colibri {

/**
 * @brief Distributed Gaussian basis for one-mode basis functions.
 * Inherits from virtual base class ModalBasis.
 *
 * @param qp_ number of quadrature points used for integration
 * @param gauss_hermite object representing gauss-hermite quadrature
 */
class DistributedGaussians : public Utils::CloneInterface<DistributedGaussians, ModalBasis> {
 public:
  // Types definition
  using Base = ModalBasis;
  using Base::getBasisSize;

  /* Class constructors */
  DistributedGaussians(const VibrationalParameters& parms, std::shared_ptr<PES> pes);
  DistributedGaussians(const std::vector<ModeGrid>& grids, std::shared_ptr<PES> pes, int nQuad);
  DistributedGaussians(const DistributedGaussians& rhs);

  /**
   * @brief Default destructor.
   */
  ~DistributedGaussians() override = default;

  /**
   * @brief Overlap calculator
   */
  double getOverlap(int i, int j, int mode) const final;

  /**
   * @brief Returns the tensor element (V_{modei, modej})ijkl in DG basis
   *
   * @param i index of the gaussian
   * @param j index of the gaussian
   * @param k index of the gaussian
   * @param l index of the gaussian
   * @param modei integer ID for identfying the mode
   * @param modej integer ID for identfying the mode
   * @return double
   */
  // double getTwoBodyIntegrals(int i, int j, int k, int l, int modei,
  //                           int modej) const;

  // double getHarmonicIntegrals(int i, int j, int mode) const;

 private:
  /**
   * @brief Get the Kinetic integral.
   */
  double getKinetic(int i, int j, int mode) const override;

  /**
   * @brief One-body potential integral calculator.
   */
  double integrateOneBodyPES(int i, int j, int mode) const final;

  /**
   * @brief Harmonic integral calculator
   */
  double integrateHarmonicPES(int i, int j, int mode) const final;

  /**
   * @brief Two-body potential integral calculator.
   */
  double integrateTwoBodyPES(int i, int j, int k, int l, int modei, int modej) const final;

  double integrateThreeBodyPES(int i, int j, int k, int l, int m, int n, int mode1, int mode2, int mode3) const final;

  /**
   * @brief Simple helper function to get the midpoint between two Distributed
   * Gaussians.
   * @param i index of the first DG.
   * @param j index of the second DG.
   * @param mode index of the mode.
   * @return double midpoint between the i-th and the j-th basis functions.
   */
  double getMidpoint(int i, int j, int mode) const;

  /* Object storing the Gauss-Hermite quadrature points */
  GaussHermiteQuad gauss_hermite;
  int qp_;
};

} // namespace Colibri
} // namespace Scine

#endif
