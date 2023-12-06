/**
 * @file ModalBasis.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#ifndef MODAL_BASIS_H
#define MODAL_BASIS_H

#include "ModeGrid.h"
#include "VibrationalUtils/PESLibrary/PES.h"
#include "VibrationalUtils/VibrationalParameters.h"

namespace Scine {
namespace Colibri {

/**
 * @brief Virtual base class for 1D basis functions.
 *
 * This class is the interface for any one-mode basis that can be used
 * in conjunction with a VSCF/MCTDH simulation. It must implement functions
 * for the calculation of the overlap, of the kinetic energy and for the
 * potential energy.
 *
 * @param grids_ vector which contains parameters of the numerical grid
 * corresponding to a given mode.
 * @param pes_ pointer to general PES object.
 * @param nModes_ total number of modes
 */

class ModalBasis {
 public:
  /**
   * @brief Constructor from a vector of ModalGrid objects
   * @param grid vector of ModalGrid objects, one per mode.
   * @param pes pointer to the underlying PES.
   */
  ModalBasis() = default;
  void initBasis(const std::vector<ModeGrid>& grid, std::shared_ptr<PES> pes);
  void initBasis(const VibrationalParameters& parms, std::shared_ptr<PES> pes);

  std::shared_ptr<ModalBasis> clone() const {
    return std::shared_ptr<ModalBasis>(this->cloneImpl());
  };

  /**
   * @brief Class constructor.
   * @param parms object storing the input parameters for a calculation.
   * @param pes underlying PES object.
   */

  /**
   * @brief Pure virtual destructor.
   */
  virtual ~ModalBasis() = default;

  /**
   * @brief Overlap between two modals
   * @param i index of the bra basis function.
   * @param j index of the ket basis function.
   * @param mode mode for which the overlap must be calculated.
   * @return double < i | j >
   */
  virtual double getOverlap(int i, int j, int mode) const = 0;

  /**
   * @brief Calculats the one-body integral for a given mode and pair of modals.
   *
   * Note that this method will return the sum of the kinetic and the potential
   * parts.
   *
   * @param i index of the bra basis function.
   * @param j index of the ket basis function.
   * @param mode mode for which the integral must be calculated.
   * @return double < i | h^{mode}(Q_{mode}) | j >
   */
  double getOneBodyIntegrals(int i, int j, int mode) const;

  /**
   * @brief Calculates the two-body integral of the hamiltonian.
   * @param i index of the bra basis function for the first mode.
   * @param j index of the bra basis function for the second mode.
   * @param k index of the ket basis function for the first mode.
   * @param l index of the ket basis function for the second mode.
   * @param modei first mode for which the integral must be calculated.
   * @param modej second mode for which the integral must be calculated.
   * @return double < i | h^{mode}(Q_{mode}) | j >
   */
  double getTwoBodyIntegrals(int i, int j, int k, int l, int modei, int modej) const;

  double getThreeBodyIntegrals(int i, int j, int k, int l, int m, int n, int mode1, int mode2, int mode3) const;

  /**
   * @brief Get the Harmonic Integral between two modals
   * @param i index of the bra basis function.
   * @param j index of the ket basis function.
   * @param mode mode for which the integral must be calculated.
   * @return double < i | H_{harm}^{mode}(Q_{mode}) | j >
   */
  double getHarmonicIntegrals(int i, int j, int mode) const;

  /**
   * @brief Gets the size of the modal basis for a given mode.
   * @param modei index of the mode
   * @return int size of the modal basis
   */
  int getBasisSize(int modei) const;

 protected:
  /**
   * @brief Calculate the integral of the kinetic energy operator.
   *
   * For a given mode modei and pair of modal bases (i,j), this function
   * should return the integral - < i | d^2 / dQ^2 | j > / 2.
   *
   * @param i index of the bra modal.
   * @param j index of the ket modal.
   * @param modei index of the mode.
   * @return double value of the kinetic energy integral.
   */
  virtual double getKinetic(int i, int j, int modei) const = 0;

  /**
   * @brief Calculate the integral of the one-body potential.
   *
   * For a given mode modei and pair of modal bases (i,j), this function
   * should return the integral < i | V^{modei}(Q_{modei}) | j >
   *
   * @param i index of the bra modal.
   * @param j index of the ket modal.
   * @param modei index of the mode.
   * @return double value of the one-body potential integral.
   */
  virtual double integrateOneBodyPES(int i, int j, int mode) const = 0;

  /**
   * @brief Calculate the integral of the one-body potential.
   *
   * For a given mode modei and pair of modal bases (i,j), this function
   * should return the integral < i | V^{modei}(Q_{modei}) | j > where V
   * represents the harmonic part of the PES.
   *
   * @param i index of the bra modal.
   * @param j index of the ket modal.
   * @param modei index of the mode.
   * @return double value of the harmonic potential integral.
   */
  virtual double integrateHarmonicPES(int i, int j, int mode) const = 0;

  /**
   * @brief Calculate the integral of the two-body potential.
   *
   * For a given mode modei and pair of modal bases (i,j), this function
   * should return the integral < i | V^{modei, mode_j} (Q_{modei}, Q_{modej}) |
   * j >
   *
   * @param i index of the bra modal for the first mode.
   * @param j index of the bra modal for the second mode.
   * @param k index of the ket modal for the first mode.
   * @param l index of the ket modal for the second mode.
   * @param modei index of the first mode.
   * @param modej index of the second mode.
   * @return double value of the two-body potential integral.
   */
  virtual double integrateTwoBodyPES(int i, int j, int k, int l, int modei, int modej) const = 0;

  virtual double integrateThreeBodyPES(int i, int j, int k, int l, int m, int n, int mode1, int mode2, int mode3) const = 0;

  /**
   * @brief Initializes the grid for each mode.
   */
  void setGrids(const VibrationalParameters& parms);

  /* Grid for each mode */
  std::vector<ModeGrid> grids_;
  /* PES underlying the basis set definition */
  std::shared_ptr<PES> pes_;
  int nModes_;

 private:
  virtual std::shared_ptr<ModalBasis> cloneImpl() const = 0;
  friend class VSCF;
};

} // namespace Colibri
} // namespace Scine
#endif
