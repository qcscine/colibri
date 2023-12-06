/**
 * @file DistributedGaussians.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#include "DistributedGaussians.h"
#include <utility>

namespace Scine::Colibri {

DistributedGaussians::DistributedGaussians(const VibrationalParameters& parms, std::shared_ptr<PES> pes)
  : qp_(parms.numqp_) {
  Base::initBasis(parms, std::move(pes));
}

DistributedGaussians::DistributedGaussians(const std::vector<ModeGrid>& grids, std::shared_ptr<PES> pes, int nQuad)
  : qp_(nQuad) {
  Base::initBasis(grids, std::move(pes));
}

DistributedGaussians::DistributedGaussians(const DistributedGaussians& rhs) : CloneInterface(rhs) {
  qp_ = rhs.qp_;
  gauss_hermite = rhs.gauss_hermite;
  Base::grids_ = rhs.grids_;
  Base::pes_ = rhs.pes_->clone();
}

double DistributedGaussians::getOverlap(int i, int j, int mode) const {
  return sqrt(M_PI) * grids_[mode].paramA_ / grids_[mode].paramB_ * exp(-grids_[mode].paramC_(i, j));
}

double DistributedGaussians::getKinetic(int i, int j, int mode) const {
  return getOverlap(i, j, mode) * pow(grids_[mode].ampl_, 2) / pow(grids_[mode].paramB_, 2) *
         (1. - 2. * grids_[mode].paramC_(i, j));
}

double DistributedGaussians::integrateOneBodyPES(int i, int j, int mode) const {
  const double prefact = grids_[mode].paramA_ * exp(-grids_[mode].paramC_(i, j));
  const double midpoint = getMidpoint(i, j, mode);
  // original paper takes b^2, actually b
  const double param = grids_[mode].paramB_;
  auto f = [this, param, midpoint, mode](const double& qi) { return pes_->getOneBodyPES(qi / param + midpoint, mode); };
  return prefact / param * gauss_hermite.oneDimIntegrate(f, qp_);
}

double DistributedGaussians::integrateTwoBodyPES(int i, int j, int k, int l, int modei, int modej) const {
  // calculates two body integrals (Vijkl)
  // i,j: mode1; k,l:mode2
  const double midpointi = getMidpoint(i, j, modei);
  const double parami = grids_[modei].paramB_;
  const double prefacti = grids_[modei].paramA_ * exp(-grids_[modei].paramC_(i, j));

  const double midpointj = getMidpoint(k, l, modej);
  const double paramj = grids_[modej].paramB_;
  const double prefactj = grids_[modej].paramA_ * exp(-grids_[modej].paramC_(k, l));

  auto g = [this, parami, paramj, midpointi, midpointj, modei, modej](const double& qi, const double& qj) {
    return pes_->getTwoBodyPES(qi / parami + midpointi, qj / paramj + midpointj, modei, modej);
  };
  double twobody = prefactj / paramj * prefacti / parami * gauss_hermite.twoDimIntegrate(g, qp_);
  return twobody;
}

double DistributedGaussians::integrateThreeBodyPES(int i, int j, int k, int l, int m, int n, int mode1, int mode2,
                                                   int mode3) const {
  // calculates three body integrals (Vijklmn)
  // i,j: mode1; k,l: mode2; m,n: mode3
  const double midpoint1 = getMidpoint(i, j, mode1);
  const double param1 = grids_[mode1].paramB_;
  const double prefact1 = grids_[mode1].paramA_ * exp(-grids_[mode1].paramC_(i, j));

  const double midpoint2 = getMidpoint(k, l, mode2);
  const double param2 = grids_[mode2].paramB_;
  const double prefact2 = grids_[mode2].paramA_ * exp(-grids_[mode2].paramC_(k, l));

  const double midpoint3 = getMidpoint(m, n, mode3);
  const double param3 = grids_[mode3].paramB_;
  const double prefact3 = grids_[mode3].paramA_ * exp(-grids_[mode3].paramC_(m, n));

  auto p = [this, param1, param2, param3, midpoint1, midpoint2, midpoint3, mode1, mode2,
            mode3](const double& q1, const double& q2, const double& q3) {
    return pes_->getThreeBodyPES(q1 / param1 + midpoint1, q2 / param2 + midpoint2, q3 / param3 + midpoint3, mode1, mode2, mode3);
  };
  double threebody = prefact1 / param1 * prefact2 / param2 * prefact3 / param3 * gauss_hermite.threeDimIntegrate(p, qp_);
  return threebody;
}

// NG: ALB: why is this weighted by the amplitude? Doesn't it divide out?
double DistributedGaussians::getMidpoint(int i, int j, int mode) const {
  return (grids_[mode].ampl_ * grids_[mode].gridVals_[i] + grids_[mode].ampl_ * grids_[mode].gridVals_[j]) /
         (2 * grids_[mode].ampl_);
}

double DistributedGaussians::integrateHarmonicPES(int i, int j, int mode) const {
  double prefact = grids_[mode].paramA_ * exp(-grids_[mode].paramC_(i, j));
  double midpoint = getMidpoint(i, j, mode);
  double param = grids_[mode].paramB_;
  auto f = [this, param, midpoint, mode](const double& qi) { return pes_->getHarmonicPES(qi / param + midpoint, mode); };
  return prefact / param * gauss_hermite.oneDimIntegrate(f, qp_);
}

} // namespace Scine::Colibri
