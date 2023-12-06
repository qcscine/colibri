/**
 * @file PesChristoffel.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#include "PesChristoffel.h"

namespace Scine::Colibri {

PesChristoffel::PesChristoffel(const VibrationalParameters& vibParms) {
  secondOrderForceConstants_ = vibParms.quadraticForceConst_;
  thirdOrderForceConstants_ = Base::sortForceConstants(vibParms.cubicForceConst_);
}

PesChristoffel::PesChristoffel(const PesChristoffel& rhs) : CloneInterface(rhs) {
  secondOrderForceConstants_ = rhs.secondOrderForceConstants_;
  thirdOrderForceConstants_ = rhs.thirdOrderForceConstants_;
}

void PesChristoffel::setReferenceGeometry(tensor<1>& /*k1*/, tensor<2>& k2, tensor<3>& k3, tensor<4>& /*k4*/) {
  secondOrderForceConstants_ = k2;
  thirdOrderForceConstants_ = Base::sortForceConstants(k3);
}

double PesChristoffel::getPES(double q, int mode) const {
  double Vi = 0;
  std::array<int, 2> key1 = {mode, mode};
  auto find1 = secondOrderForceConstants_.find(key1);
  if (find1 != secondOrderForceConstants_.end()) {
    Vi += 1. / 2 * find1->second * pow(q, 2);
  }
  std::array<int, 3> key2 = {mode, mode, mode};
  auto find2 = thirdOrderForceConstants_.find(key2);
  if (find2 != thirdOrderForceConstants_.end()) {
    Vi += find2->second * pow(q, 3);
  }
  return Vi;
}

double PesChristoffel::getPES(double qi, double qj, int modei, int modej) const {
  double Vij = 0;
  std::array<int, 3> key1 = {modei, modej, modej};
  std::sort(key1.begin(), key1.end());
  auto find1 = thirdOrderForceConstants_.find(key1);
  std::array<int, 3> key2 = {modej, modei, modei};
  std::sort(key2.begin(), key2.end());
  auto find2 = thirdOrderForceConstants_.find(key2);
  if (find1 != thirdOrderForceConstants_.end()) {
    Vij += find1->second * qi * pow(qj, 2);
  }
  if (find2 != thirdOrderForceConstants_.end()) {
    Vij += find2->second * qj * pow(qi, 2);
  }
  return Vij;
}

double PesChristoffel::getHarmonicPES(double q, int mode) const {
  double Vi = 0;
  std::array<int, 2> key1 = {mode, mode};
  auto find1 = secondOrderForceConstants_.find(key1);
  if (find1 != secondOrderForceConstants_.end()) {
    Vi += 1. / 2 * find1->second * pow(q, 2);
  }
  return Vi;
}

} // namespace Scine::Colibri
