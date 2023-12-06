/**
 * @file PesCoupledOscillator.cpp
 *
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 *
 * This model is taken from Leclerc and Carrington, JCP 2014
 */
#include "PesCoupledOscillator.h"
#include <iostream>

namespace Scine::Colibri {

PesCoupledOscillator::PesCoupledOscillator(const VibrationalParameters& vibParms, ScalarType alpha) {
  coupling_ = alpha;
  int numberOfModes = vibParms.numModes_;
  for (int iMode = 0; iMode < numberOfModes; iMode++) {
    // Harmonic part
#pragma omp critical
    {
      std::array<int, 2> tmp = {iMode, iMode};
      secondOrderForceConstants_[tmp] = static_cast<ScalarType>((iMode + 1.) / 2.);
      for (int jMode = 0; jMode < iMode; jMode++) {
        std::array<int, 2> tmp = {jMode, iMode};
        secondOrderForceConstants_[tmp] =
            static_cast<ScalarType>(std::sqrt(std::sqrt((iMode + 1.) / 2.) * std::sqrt((jMode + 1.) / 2.)) * coupling_);
      }
    }
  }
}

PesCoupledOscillator::PesCoupledOscillator(const PesCoupledOscillator& rhs) : CloneInterface(rhs) {
  secondOrderForceConstants_ = rhs.secondOrderForceConstants_;
}

double PesCoupledOscillator::getPES(double q, int mode) const {
  double Vi = 0;
  std::array<int, 2> key1 = {mode, mode};
#pragma omp critical
  {
    auto find1 = secondOrderForceConstants_.find(key1);
    if (find1 != secondOrderForceConstants_.end()) {
      Vi += 1. / 2. * find1->second * pow(q, 2);
    }
  }
  return Vi;
}

double PesCoupledOscillator::getPES(double qi, double qj, int modei, int modej) const {
  double Vij = 0.;
  std::array<int, 2> key1;
  if (modei < modej) {
    key1 = {modei, modej};
  } else {
    key1 = {modej, modei};
  };
#pragma omp critical
  {
    auto find1 = secondOrderForceConstants_.find(key1);
    if (find1 != secondOrderForceConstants_.end()) {
      Vij += find1->second * qi * qj;
    }
  }
  return Vij;
}

double PesCoupledOscillator::getHarmonicPES(double q, int mode) const {
  double Vi = 0;
  std::array<int, 2> key1 = {mode, mode};
#pragma omp critical
  {
    auto find1 = secondOrderForceConstants_.find(key1);
    if (find1 != secondOrderForceConstants_.end()) {
      Vi += 1. / 2 * find1->second * pow(q, 2);
    }
  }
  return Vi;
}

} // namespace Scine::Colibri
