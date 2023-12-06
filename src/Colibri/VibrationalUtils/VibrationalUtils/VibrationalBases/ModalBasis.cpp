/**
 * @file ModalBasis.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#include "ModalBasis.h"
#include <utility>

namespace Scine::Colibri {

void ModalBasis::initBasis(const std::vector<ModeGrid>& grid, std::shared_ptr<PES> pes) {
  pes_ = std::move(pes);
  nModes_ = grid.size();
  grids_.reserve(nModes_);
  for (const auto& iGrid : grid) {
    grids_.push_back(iGrid);
  }
}

void ModalBasis::initBasis(const VibrationalParameters& parms, std::shared_ptr<PES> pes) {
  nModes_ = parms.numModes_;
  pes_ = std::move(pes);
  setGrids(parms);
}

void ModalBasis::setGrids(const VibrationalParameters& parms) {
  grids_.reserve(nModes_);
  for (int mode = 0; mode < nModes_; mode++) {
    ModeGrid g(parms.modeFreqs_[mode], parms.quantumNums_[mode], parms.numGridPoints_[mode]);
    grids_.emplace_back(std::move(g));
  }
}

int ModalBasis::getBasisSize(int modei) const {
  return grids_[modei].gridSize_;
}

double ModalBasis::getOneBodyIntegrals(int i, int j, int mode) const {
  return getKinetic(i, j, mode) + integrateOneBodyPES(i, j, mode);
}

double ModalBasis::getTwoBodyIntegrals(int i, int j, int k, int l, int modei, int modej) const {
  return integrateTwoBodyPES(i, j, k, l, modei, modej);
}

double ModalBasis::getThreeBodyIntegrals(int i, int j, int k, int l, int m, int n, int mode1, int mode2, int mode3) const {
  return integrateThreeBodyPES(i, j, k, l, m, n, mode1, mode2, mode3);
}

double ModalBasis::getHarmonicIntegrals(int i, int j, int mode) const {
  return integrateHarmonicPES(i, j, mode) + getKinetic(i, j, mode);
}

} // namespace Scine::Colibri
