/**
 * @file DVR.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#ifndef DVR_H
#define DVR_H

#include "ModalBasis.h"
#include "ModeGrid.h"
#include "VibrationalUtils/PESLibrary/PES.h"
#include "VibrationalUtils/VibrationalParameters.h"
#include <Utils/Technical/CloneInterface.h>
#include <memory>

namespace Scine {
namespace Colibri {

class DVR : public Utils::CloneInterface<DVR, ModalBasis> {
 public:
  // Types definition
  using Base = ModalBasis;
  using Base::getBasisSize;

  /**
   * @brief Class constructor
   * @param parms object storing the parameters of the calculation.
   * @param pes pointer to the PES of the Hamiltonian.
   */
  DVR(const VibrationalParameters& parms, std::shared_ptr<PES> pes);
  DVR(const DVR& rhs);
  DVR(const std::vector<ModeGrid>& grids, std::shared_ptr<PES> pes);

  /**
   * @brief Default destructor (overrides the virtual one of the base class)
   */
  ~DVR() override = default;

  /**
   * @brief Overlap calculator.
   */
  double getOverlap(int i, int j, int mode) const final;

  double getHarmonicIntegrals(int i, int j, int mode) const;

 private:
  /**
   * @brief Kinetic energy calculator.
   */
  double getKinetic(int i, int j, int mode) const final;

  /**
   * @brief One-body potential integral calculator.
   */
  double integrateOneBodyPES(int i, int j, int mode) const final;

  /**
   * @brief Harmonic potential integrator.
   */
  double integrateHarmonicPES(int i, int j, int mode) const final;

  /**
   * @brief Two-body potential integral calculator.
   */
  double integrateTwoBodyPES(int i, int j, int k, int l, int modei, int modej) const final;

  /**
   * @brief Thre-body potential integral calculator.
   */
  double integrateThreeBodyPES(int i, int j, int k, int l, int m, int n, int modei, int modej, int modek) const final;
};
} // namespace Colibri
} // namespace Scine

#endif
