/**
 * @file DVR.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#include "DVR.h"
#include <cmath>
#include <utility>

namespace Scine::Colibri {

DVR::DVR(const VibrationalParameters& parms, std::shared_ptr<PES> pes) {
  Base::initBasis(parms, std::move(pes));
}

DVR::DVR(const std::vector<ModeGrid>& grids, std::shared_ptr<PES> pes) {
  Base::initBasis(grids, std::move(pes));
}

DVR::DVR(const DVR& rhs) : CloneInterface(rhs) {
  // pes needs to be cloned
  Base::grids_ = rhs.grids_;
  Base::pes_ = rhs.pes_->clone();
}

double DVR::getOverlap(int i, int j, int /*mode*/) const {
  return (i == j) ? 1.0 : 0.0;
}

double DVR::getKinetic(int i, int j, int mode) const {
  double ba = std::abs(2 * (grids_[mode].qMax_ + grids_[mode].qStep_));
  int nGrid = getBasisSize(mode) + 1;
  double tmpt = NAN;
  i = i + 1;
  j = j + 1;
  if (i != j) {
    tmpt = (pow(M_PI, 2) * pow((-1.), (i - j))) / (4 * pow(ba, 2));
    tmpt = tmpt * (1. / (pow((sin(M_PI * (i - j) / (2 * nGrid))), 2)) - 1. / (pow((sin(M_PI * (i + j) / (2 * nGrid))), 2)));
  } else {
    tmpt = pow(M_PI, 2) / (4 * pow((ba), 2));
    tmpt = tmpt * ((2 * pow(nGrid, 2) + 1.) / 3. - 1. / pow(sin(M_PI * i / nGrid), 2));
  }
  return tmpt;
}

double DVR::integrateOneBodyPES(int i, int j, int mode) const {
  return (i == j) ? pes_->getOneBodyPES(grids_[mode].gridVals_[i], mode) : 0.;
}

double DVR::integrateTwoBodyPES(int i, int j, int k, int l, int modei, int modej) const {
  double qi = NAN;
  double qj = NAN;
  double res = 0.;
  if (i == j && k == l) {
    qi = grids_[modei].gridVals_[i];
    qj = grids_[modej].gridVals_[k];
    res = pes_->getTwoBodyPES(qi, qj, modei, modej);
  }
  return res;
}

double DVR::integrateThreeBodyPES(int i, int j, int k, int l, int m, int n, int modei, int modej, int modek) const {
  double qi = NAN;
  double qj = NAN;
  double qk = NAN;
  double res = 0.;
  if (i == j && k == l && m == n) {
    qi = grids_[modei].gridVals_[i];
    qj = grids_[modej].gridVals_[k];
    qk = grids_[modek].gridVals_[m];
    res = pes_->getThreeBodyPES(qi, qj, qk, modei, modej, modek);
  }
  return res;
}

double DVR::integrateHarmonicPES(int /*i*/, int /*j*/, int /*mode*/) const {
  return 0.0;
}

} // namespace Scine::Colibri
