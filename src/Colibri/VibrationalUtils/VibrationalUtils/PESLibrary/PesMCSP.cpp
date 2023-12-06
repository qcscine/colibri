/**
 * @file PesMCSP.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#include "PesMCSP.h"

namespace Scine::Colibri {

PesMCSP::PesMCSP(const VibrationalParameters& vibParms) {
  oneModeForceConstants_ = vibParms.oneModeTerms_;
  twoModeForceConstants_ = Base::sortForceConstants(vibParms.twoModeTerms_);
  threeModeForceConstants_ = Base::sortForceConstants(vibParms.threeModeTerms_);
  fourModeForceConstants_ = Base::sortForceConstants(vibParms.fourModeTerms_);
  fiveModeForceConstants_ = Base::sortForceConstants(vibParms.fiveModeTerms_);
  sixModeForceConstants_ = Base::sortForceConstants(vibParms.sixModeTerms_);
}

PesMCSP::PesMCSP(const PesMCSP& rhs) : CloneInterface(rhs) {
  oneModeForceConstants_ = rhs.oneModeForceConstants_;
  twoModeForceConstants_ = rhs.twoModeForceConstants_;
  threeModeForceConstants_ = rhs.threeModeForceConstants_;
  fourModeForceConstants_ = rhs.fourModeForceConstants_;
  fiveModeForceConstants_ = rhs.fiveModeForceConstants_;
  sixModeForceConstants_ = rhs.sixModeForceConstants_;
}

void PesMCSP::setReferenceGeometry(tensorMCSP<1>& k1, tensorMCSP<2>& k2, tensorMCSP<3>& k3, tensorMCSP<4>& k4,
                                   tensorMCSP<5>& k5, tensorMCSP<6>& k6) {
  oneModeForceConstants_ = k1;
  twoModeForceConstants_ = Base::sortForceConstants(k2);
  threeModeForceConstants_ = Base::sortForceConstants(k3);
  fourModeForceConstants_ = Base::sortForceConstants(k4);
  fiveModeForceConstants_ = Base::sortForceConstants(k5);
  sixModeForceConstants_ = Base::sortForceConstants(k6);
}

// One-mode variation
double PesMCSP::getPES(double q, int mode) const {
  double Vmi = 0.;
  std::array<int, 1> key = {mode};
  auto find = oneModeForceConstants_.find(key);
  if (find != oneModeForceConstants_.end()) {
    for (const auto& it : find->second) {
      Vmi += it.second * pow(q, it.first[0]);
    }
  }
  return Vmi;
}

// Two-mode variation
double PesMCSP::getPES(double qi, double qj, int modei, int modej) const {
  double Vmimj = 0.;
  std::array<int, 2> key;
  std::array<double, 2> qs;

  std::array<std::pair<int, double>, 2> inds_qs = {std::make_pair(modei, qi), std::make_pair(modej, qj)};
  std::sort(inds_qs.begin(), inds_qs.end(),
            [](std::pair<int, double> i, std::pair<int, double> j) { return i.first < j.first; });
  for (int i = 0; i < inds_qs.size(); i++) {
    key[i] = inds_qs[i].first;
    qs[i] = inds_qs[i].second;
  }

  auto find = twoModeForceConstants_.find(key);
  if (find != twoModeForceConstants_.end()) {
    for (const auto& it : find->second) {
      Vmimj += it.second * pow(qs[0], it.first[0]) * pow(qs[1], it.first[1]);
    }
  }
  return Vmimj;
}

// this should be truly harmonic, not first- and second-order as in PesQFF
double PesMCSP::getHarmonicPES(double q, int mode) const {
  double Vmi = 0.;
  std::array<int, 1> key = {mode};
  auto find = oneModeForceConstants_.find(key);
  if (find != oneModeForceConstants_.end()) {
    for (const auto& it : find->second) {
      if (it.first[0] == 2) {
        Vmi += it.second * pow(q, 2);
        break;
      }
    }
  }
  return Vmi;
}

// Three-mode variation
double PesMCSP::getPES(double qi, double qj, double qk, int modei, int modej, int modek) const {
  double Vmimjmk = 0.;

  std::array<int, 3> key;
  std::array<double, 3> qs;

  std::array<std::pair<int, double>, 3> inds_qs = {std::make_pair(modei, qi), std::make_pair(modej, qj),
                                                   std::make_pair(modek, qk)};
  std::sort(inds_qs.begin(), inds_qs.end(),
            [](std::pair<int, double> i, std::pair<int, double> j) { return i.first < j.first; });
  for (int i = 0; i < inds_qs.size(); i++) {
    key[i] = inds_qs[i].first;
    qs[i] = inds_qs[i].second;
  }

  auto find = threeModeForceConstants_.find(key);
  if (find != threeModeForceConstants_.end()) {
    for (const auto& it : find->second) {
      Vmimjmk += it.second * pow(qs[0], it.first[0]) * pow(qs[1], it.first[1]) * pow(qs[2], it.first[2]);
    }
  }
  return Vmimjmk;
}

} // namespace Scine::Colibri
