/**
 * @file PesCoupledOscillator.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#ifndef PES_COUPLED_OSC_H
#define PES_COUPLED_OSC_H

#include "PES.h"
#include "VibrationalUtils/Tensor.h"
#include "VibrationalUtils/VibrationalParameters.h"
#include <Utils/Technical/CloneInterface.h>
#include <cmath>

namespace Scine {
namespace Colibri {

class PesCoupledOscillator : public Utils::CloneInterface<PesCoupledOscillator, PES> {
 public:
  // Types declaration
  using Base = PES;
  using ScalarType = double;
  /* Default constructor */
  PesCoupledOscillator() = default;
  /* Default destructor */
  ~PesCoupledOscillator() override = default;
  /* Constructor from input parameters */
  PesCoupledOscillator(const VibrationalParameters& vibParms, ScalarType alpha);
  /* Copy constructor */
  PesCoupledOscillator(const PesCoupledOscillator& rhs);

  /**
   * @brief Getter for the energy of the equilibrium structure
   * @return 0. in all cases
   */
  double getPES() const override {
    return 0.0;
  };

  void setReferenceGeometry(tensor<1>& /*k1*/, tensor<2>& /*k2*/, tensor<3>& /*k3*/, tensor<4>& /*k4*/) override {
    throw std::runtime_error("setReferenceGeometry not available for the PesCoupledOscillator");
  };

 private:
  /* Private getters for the energy along given modes */
  double getPES(double q, int modei) const override;
  double getPES(double qi, double qj, int modei, int modej) const override;
  double getHarmonicPES(double q, int mode) const override;
  /* Container for the second-order force constants */
  tensor<2> secondOrderForceConstants_;
  /* Scaling factor for the coupling term */
  ScalarType coupling_;
};

} // namespace Colibri
} // namespace Scine
#endif
