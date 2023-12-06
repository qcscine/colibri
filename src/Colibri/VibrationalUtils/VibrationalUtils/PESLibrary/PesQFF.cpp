/**
 * @file PesQFF.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#include "PesQFF.h"

namespace Scine::Colibri {

PesQFF::PesQFF(const VibrationalParameters& vibParms) {
  firstOrderForceConstants_ = vibParms.linForceConst_;
  secondOrderForceConstants_ = vibParms.quadraticForceConst_;
  thirdOrderForceConstants_ = Base::sortForceConstants(vibParms.cubicForceConst_);
  fourthOrderForceConstants_ = Base::sortForceConstants(vibParms.quarticForceConst_);
}

PesQFF::PesQFF(const PesQFF& rhs) : CloneInterface(rhs) {
  firstOrderForceConstants_ = rhs.firstOrderForceConstants_;
  secondOrderForceConstants_ = rhs.secondOrderForceConstants_;
  thirdOrderForceConstants_ = rhs.thirdOrderForceConstants_;
  fourthOrderForceConstants_ = rhs.fourthOrderForceConstants_;
}

void PesQFF::setReferenceGeometry(tensor<1>& k1, tensor<2>& k2, tensor<3>& k3, tensor<4>& k4) {
  firstOrderForceConstants_ = k1;
  secondOrderForceConstants_ = k2;
  thirdOrderForceConstants_ = Base::sortForceConstants(k3);
  fourthOrderForceConstants_ = Base::sortForceConstants(k4);
}

double PesQFF::getPES(double q, int mode) const {
  double Vmi = 0.;
  std::array<int, 1> key3 = {mode};
  auto find3 = firstOrderForceConstants_.find(key3);
  std::array<int, 2> key4 = {mode, mode};
  auto find4 = secondOrderForceConstants_.find(key4);
  if (find3 != firstOrderForceConstants_.end()) {
    Vmi += find3->second * q;
  }
  if (find4 != secondOrderForceConstants_.end()) {
    Vmi += 1. / 2 * find4->second * pow(q, 2);
  }

  std::array<int, 3> key1 = {mode, mode, mode};
  auto find1 = thirdOrderForceConstants_.find(key1);
  std::array<int, 4> key2 = {mode, mode, mode, mode};
  auto find2 = fourthOrderForceConstants_.find(key2);
  // check whether key present and add term
  if (find1 != thirdOrderForceConstants_.end()) {
    Vmi += 1. / 6 * find1->second * pow(q, 3);
  }
  if (find2 != fourthOrderForceConstants_.end()) {
    Vmi += 1. / 24 * find2->second * pow(q, 4);
  }
  return Vmi;
}

double PesQFF::getPES(double qi, double qj, int modei, int modej) const {
  // if not normal coords, add: secondOrderForceConst(modei, modej)*qi*qj;
  double Vmimj = 0.;
  std::array<int, 3> key1 = {modei, modej, modej};
  std::sort(key1.begin(), key1.end());
  auto find1 = thirdOrderForceConstants_.find(key1);
  std::array<int, 3> key2 = {modej, modei, modei};
  std::sort(key2.begin(), key2.end());
  auto find2 = thirdOrderForceConstants_.find(key2);
  std::array<int, 4> key3 = {modej, modei, modei, modei};
  std::sort(key3.begin(), key3.end());
  auto find3 = fourthOrderForceConstants_.find(key3);
  std::array<int, 4> key4 = {modei, modej, modej, modej};
  std::sort(key4.begin(), key4.end());
  auto find4 = fourthOrderForceConstants_.find(key4);
  std::array<int, 4> key5 = {modei, modei, modej, modej};
  std::sort(key5.begin(), key5.end());
  auto find5 = fourthOrderForceConstants_.find(key5);

  if (find1 != thirdOrderForceConstants_.end()) {
    Vmimj += 1. / 2 * find1->second * qi * pow(qj, 2);
  }
  if (find2 != thirdOrderForceConstants_.end()) {
    Vmimj += 1. / 2 * find2->second * qj * pow(qi, 2);
  }
  if (find3 != fourthOrderForceConstants_.end()) {
    Vmimj += 1. / 6 * find3->second * qj * pow(qi, 3);
  }
  if (find4 != fourthOrderForceConstants_.end()) {
    Vmimj += 1. / 6 * find4->second * qi * pow(qj, 3);
  }
  if (find5 != fourthOrderForceConstants_.end()) {
    Vmimj += 1. / 4. * find5->second * pow(qi, 2) * pow(qj, 2);
  }
  return Vmimj;
}

double PesQFF::getHarmonicPES(double q, int mode) const {
  double Vmi = 0.;
  std::array<int, 1> key3 = {mode};
  auto find3 = firstOrderForceConstants_.find(key3);
  std::array<int, 2> key4 = {mode, mode};
  auto find4 = secondOrderForceConstants_.find(key4);
  if (find3 != firstOrderForceConstants_.end()) {
    Vmi += find3->second * q;
  }
  if (find4 != secondOrderForceConstants_.end()) {
    Vmi += 1. / 2 * find4->second * pow(q, 2);
  }
  return Vmi;
}

} // namespace Scine::Colibri
