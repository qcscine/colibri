/**
 * @file PesQFF.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#ifndef PES_QFF_H
#define PES_QFF_H

#include "PES.h"
#include "VibrationalUtils/VibrationalParameters.h"
#include <Utils/Technical/CloneInterface.h>
#include <cmath>
#include <vector>

namespace Scine {
namespace Colibri {

class PesQFF : public Utils::CloneInterface<PesQFF, PES> {
 public:
  using Base = PES;
  PesQFF() = default;
  ~PesQFF() override = default;
  PesQFF(const VibrationalParameters& vibParms);
  PesQFF(const PesQFF& rhs);
  void setReferenceGeometry(tensor<1>& k1, tensor<2>& k2, tensor<3>& k3, tensor<4>& k4) override;
  double getPES() const override {
    return 0.0;
  };

 private:
  double getPES(double q, int modei) const override;
  double getPES(double qi, double qj, int modei, int modej) const override;
  double getHarmonicPES(double q, int mode) const override;
  tensor<1> firstOrderForceConstants_;
  tensor<2> secondOrderForceConstants_;
  tensor<3> thirdOrderForceConstants_;
  tensor<4> fourthOrderForceConstants_;
};

} // namespace Colibri
} // namespace Scine

#endif
