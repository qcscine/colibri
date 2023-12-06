/**
 * @file ModeGrid.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#include "ModeGrid.h"

namespace Scine::Colibri {

const double ModeGrid::paramc_ = 0.7;

ModeGrid::ModeGrid(double freq, double nquantum, int gridsize)
  : VibrationalMode{freq, static_cast<int>(nquantum)}, gridSize_(gridsize) {
  qMax_ = sqrt(2. * (nquantum_ + 0.5) / freq_);
  qStep_ = 2. * qMax_ / (gridSize_ - 1);
  gridVals_ = getBasisGrid();
  // sets parameters A, B, C for Distributed Gaussians
  // TODO make it DG-specific
  setParams();
}

std::vector<double> ModeGrid::getBasisGrid() const {
  std::vector<double> gridVals(gridSize_);
  for (int k = 0; k < gridSize_; k++) {
    gridVals[k] = -qMax_ + k * qStep_;
  }
  return gridVals;
}

int ModeGrid::getSize() const {
  return gridSize_;
}

void ModeGrid::setParams() {
  ampl_ = pow(paramc_, 2) / pow(qStep_, 2);
  paramA_ = pow((4 * ampl_ * ampl_ / pow(M_PI, 2)), 1. / 4);
  paramB_ = pow((2 * ampl_), 1. / 2);
  paramC_.resize(gridSize_, gridSize_);
  for (int i = 0; i < gridSize_; i++) {
    for (int j = 0; j < gridSize_; j++) {
      double r = gridVals_[i] - gridVals_[j];
      paramC_(i, j) = (pow(ampl_, 2) / pow(paramB_, 2)) * pow(r, 2);
    }
  }
}

} // namespace Scine::Colibri
