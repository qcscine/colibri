/**
 * @file ModeGrid.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#ifndef MODE_GRID_H
#define MODE_GRID_H

#include <VibrationalUtils/VibrationalMode.h>
#include <Eigen/Dense>
#include <vector>

namespace Scine {
namespace Colibri {

/**
 * @brief Small class representing the integration grid for a given mode.
 *
 * In practice, when performing any VSCF/VCI calculation, an interval [-Qmax,
 * Qmax] is set to perform any integration. This class contains the information
 * regarding this integral.
 * The grid is constructed as in DOI:10.1002/cphc.201402251.
 */
class ModeGrid : public VibrationalMode {
 public:
  // -- Types declaration --
  using ListOfDisplacementsType = Eigen::MatrixXd;

  /* Default class constructor */
  ModeGrid() = default;

  /**
   * @brief Constructor for a given frequency
   *
   * @param freq Harmonic Frequency (in Hartree)
   * @param nquantum maximum number of quanta of excitation
   * @param gridsize number of points in the grid.
   */
  ModeGrid(double freq, double nquantum, int gridsize);

  /* Default destructors */
  ~ModeGrid() = default;

  /* Getter of the maximum displacement */
  double getQMax() const {
    return qMax_;
  }

  /* Getter for the basis grid */
  std::vector<double> getBasisGrid() const;
  /* Getter for the size */
  int getSize() const;

  /* The basis objects can access the private members */
  friend class ModalBasis;
  friend class DVR;
  friend class DistributedGaussians;

 private:
  void setParams();
  double setQstep();
  double qMax_, qStep_, paramB_, paramA_, ampl_;
  Eigen::MatrixXd paramC_;
  int gridSize_;
  std::vector<double> gridVals_;
  static const double paramc_;
};

} // namespace Colibri
} // namespace Scine

#endif
