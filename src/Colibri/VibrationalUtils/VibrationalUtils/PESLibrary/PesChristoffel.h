/**
 * @file PesChristoffel.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#ifndef PES_CH_H
#define PES_CH_H

#include "PES.h"
#include <Utils/Technical/CloneInterface.h>
#include <VibrationalUtils/Tensor.h>
#include <VibrationalUtils/VibrationalParameters.h>
#include <cmath>

namespace Scine {
namespace Colibri {

class PesChristoffel : public Utils::CloneInterface<PesChristoffel, PES> {
 public:
  /**
   * @brief Implementation of the model potential of Christoffel
   *
   * @param secondOrder Eigen matrix storing the quadratic force constants
   * @param thirdOrder  unordered map storing the cubic force constants. The key
   * is an array of indices.
   */
  using Base = PES;
  PesChristoffel() = default;
  ~PesChristoffel() override = default;
  PesChristoffel(const VibrationalParameters& vibParms);
  PesChristoffel(const PesChristoffel& rhs);

  double getPES() const override {
    return 0.0;
  };

  void setReferenceGeometry(tensor<1>& k1, tensor<2>& k2, tensor<3>& k3, tensor<4>& k4) override;

 private:
  double getPES(double q, int modei) const override;
  double getPES(double qi, double qj, int modei, int modej) const override;
  double getHarmonicPES(double q, int mode) const override;
  tensor<2> secondOrderForceConstants_;
  tensor<3> thirdOrderForceConstants_;
};

} // namespace Colibri
} // namespace Scine
#endif
