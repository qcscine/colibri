/**
 * @file PesTaylorParam.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#ifndef PES_TayParam_H
#define PES_TayParam_H

#include "PES.h"
#include "Utils/Technical/CloneInterface.h"
#include "VibrationalUtils/VibrationalParameters.h"
#include <cmath>
#include <vector>

namespace Scine {
namespace Colibri {

class PesTayParam : public Utils::CloneInterface<PesTayParam, PES> {
 public:
  using Base = PES;
  PesTayParam() = default;
  ~PesTayParam() override = default;
  PesTayParam(const VibrationalParameters& vibParms);
  PesTayParam(const PesTayParam& rhs);
  void setReferenceGeometry(tensorMC<1>& k1, tensorMC<2>& k2, tensorMC<3>& k3, tensorMC<4>& k4, tensorMC<5>& k5,
                            tensorMC<6>& k6, int expOrder) override;
  void setReferenceGeometry(tensor<1>& k1, tensor<2>& k2, tensor<3>& k3, tensor<4>& k4) override;
  double getPES() const override {
    return 0.0;
  };

 private:
  int taylorExpansionOrder_;
  double getPES(double q, int modei) const override;
  double getPES(double qi, double qj, int modei, int modej) const override;
  double getPES(double qi, double qj, double qk, int modei, int modej, int modek) const override;
  double getHarmonicPES(double q, int mode) const override;
  tensorMC<1> oneModeForceConstants_;
  tensorMC<2> twoModeForceConstants_;
  tensorMC<3> threeModeForceConstants_;
  tensorMC<4> fourModeForceConstants_;
  tensorMC<5> fiveModeForceConstants_;
  tensorMC<6> sixModeForceConstants_;
  static bool sortInds(std::pair<int, int> ind1, std::pair<int, int> ind2);
};

} // namespace Colibri
} // namespace Scine

#endif
