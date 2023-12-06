/**
 * @file PESTaylorParam.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#include "PesTaylorParam.h"

namespace Scine::Colibri {

PesTayParam::PesTayParam(const VibrationalParameters& vibParms) {
  oneModeForceConstants_ = vibParms.oneModeCoeffs_;
  twoModeForceConstants_ = vibParms.twoModeCoeffs_;
  threeModeForceConstants_ = Base::sortForceConstants(vibParms.threeModeCoeffs_);
  fourModeForceConstants_ = Base::sortForceConstants(vibParms.fourModeCoeffs_);
  fiveModeForceConstants_ = Base::sortForceConstants(vibParms.fiveModeCoeffs_);
  sixModeForceConstants_ = Base::sortForceConstants(vibParms.sixModeCoeffs_);
  taylorExpansionOrder_ = vibParms.taylorExpOrder_;
}

PesTayParam::PesTayParam(const PesTayParam& rhs) : CloneInterface(rhs) {
  oneModeForceConstants_ = rhs.oneModeForceConstants_;
  twoModeForceConstants_ = rhs.twoModeForceConstants_;
  threeModeForceConstants_ = rhs.threeModeForceConstants_;
  fourModeForceConstants_ = rhs.fourModeForceConstants_;
  fiveModeForceConstants_ = rhs.fiveModeForceConstants_;
  sixModeForceConstants_ = rhs.sixModeForceConstants_;
  taylorExpansionOrder_ = rhs.taylorExpansionOrder_;
}

void PesTayParam::setReferenceGeometry(tensorMC<1>& k1, tensorMC<2>& k2, tensorMC<3>& k3, tensorMC<4>& k4,
                                       tensorMC<5>& k5, tensorMC<6>& k6, int expOrder) {
  oneModeForceConstants_ = k1;
  twoModeForceConstants_ = k2;
  threeModeForceConstants_ = Base::sortForceConstants(k3);
  fourModeForceConstants_ = Base::sortForceConstants(k4);
  fiveModeForceConstants_ = Base::sortForceConstants(k5);
  sixModeForceConstants_ = Base::sortForceConstants(k6);
  taylorExpansionOrder_ = expOrder;
}

void PesTayParam::setReferenceGeometry(tensor<1>& /*k1*/, tensor<2>& /*k2*/, tensor<3>& /*k3*/, tensor<4>& /*k4*/) {
  exit(2);
}

// One-mode variation
double PesTayParam::getPES(double q, int mode) const {
  double Vmi = 0.;
  std::array<std::pair<int, int>, 1> key;
  for (int i = 1; i <= taylorExpansionOrder_; i++) {
    key = {std::make_pair(mode, i)};
    auto find = oneModeForceConstants_.find(key);
    if (find != oneModeForceConstants_.end()) {
      Vmi += find->second * pow(q, i);
    }
  }
  return Vmi;
}

// Two-mode variation
double PesTayParam::getPES(double qi, double qj, int modei, int modej) const {
  double Vmi = getPES(qi, modei);
  double Vmj = getPES(qj, modej);

  double Vmimj = 0.;
  std::array<std::pair<int, int>, 2> key;
  std::array<std::pair<int, int>, 2> keySwitched;

  for (int i = 1; i <= taylorExpansionOrder_ - 1; i++) {
    for (int j = 1; j <= taylorExpansionOrder_ - i; j++) {
      key = {std::make_pair(modei, i), std::make_pair(modej, j)};
      auto find = twoModeForceConstants_.find(key);
      if (find != twoModeForceConstants_.end()) {
        Vmimj += find->second * pow(qi, i) * pow(qj, j);
      } else {
        keySwitched = {std::make_pair(modej, j), std::make_pair(modei, i)};
        auto findS = twoModeForceConstants_.find(keySwitched);
        if (findS != twoModeForceConstants_.end()) {
          Vmimj += findS->second * pow(qi, i) * pow(qj, j);
        }
      }
    }
  }
  return (Vmimj + Vmi + Vmj);
}

// this should be truly harmonic, not first- and second-order as in PesQFF
double PesTayParam::getHarmonicPES(double q, int mode) const {
  double Vmi = 0.;
  std::array<std::pair<int, int>, 1> key = {std::make_pair(mode, 2)};
  auto find = oneModeForceConstants_.find(key);
  if (find != oneModeForceConstants_.end()) {
    Vmi += find->second * pow(q, 2);
  }
  return Vmi;
}

// Three-mode variation
double PesTayParam::getPES(double qi, double qj, double qk, int modei, int modej, int modek) const {
  double Vmimjmk = 0.;

  std::array<std::pair<int, int>, 3> key;

  for (int i = 1; i <= taylorExpansionOrder_ - 2; i++) {
    for (int j = 1; j <= taylorExpansionOrder_ - i - 1; j++) {
      for (int k = 1; k <= taylorExpansionOrder_ - i - j; k++) {
        key = {std::make_pair(modei, i), std::make_pair(modej, j), std::make_pair(modek, k)};

        std::sort(key.begin(), key.end(), sortInds);
        do {
          auto find = threeModeForceConstants_.find(key);
          if (find != threeModeForceConstants_.end()) {
            Vmimjmk += find->second * pow(qi, i) * pow(qj, j) * pow(qk, k);
            break;
          }
        } while (std::next_permutation(key.begin(), key.end(), sortInds));
      }
    }
  }
  return (Vmimjmk);
}

bool PesTayParam::sortInds(std::pair<int, int> ind1, std::pair<int, int> ind2) {
  return (ind1.first < ind2.first);
}

} // namespace Scine::Colibri
