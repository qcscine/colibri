/**
 * @file PesMCSP.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#ifndef PES_MCSP_H
#define PES_MCSP_H

#include "PES.h"
#include <Utils/Technical/CloneInterface.h>
#include <VibrationalUtils/VibrationalParameters.h>
#include <cmath>
#include <vector>

namespace Scine {
namespace Colibri {

class PesMCSP : public Utils::CloneInterface<PesMCSP, PES> {
 public:
  using Base = PES;
  PesMCSP() = default;
  ~PesMCSP() override = default;
  PesMCSP(const VibrationalParameters& vibParms);
  PesMCSP(const PesMCSP& rhs);
  void setReferenceGeometry(tensorMCSP<1>& k1, tensorMCSP<2>& k2, tensorMCSP<3>& k3, tensorMCSP<4>& k4,
                            tensorMCSP<5>& k5, tensorMCSP<6>& k6) override;
  double getPES() const override {
    return 0.0;
  };

 private:
  double getPES(double q, int modei) const override;
  double getPES(double qi, double qj, int modei, int modej) const override;
  double getPES(double qi, double qj, double qk, int modei, int modej, int modek) const override;
  double getHarmonicPES(double q, int mode) const override;
  tensorMCSP<1> oneModeForceConstants_;
  tensorMCSP<2> twoModeForceConstants_;
  tensorMCSP<3> threeModeForceConstants_;
  tensorMCSP<4> fourModeForceConstants_;
  tensorMCSP<5> fiveModeForceConstants_;
  tensorMCSP<6> sixModeForceConstants_;
};

} // namespace Colibri
} // namespace Scine

#endif
