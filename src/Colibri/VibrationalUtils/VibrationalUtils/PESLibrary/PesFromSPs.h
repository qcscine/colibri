/**
 * @file PesFromSPs.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#ifndef PES_FROM_SPS_H
#define PES_FROM_SPS_H

#include "PES.h"
#include "VibrationalUtils/VibrationalParameters.h"
#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/CalculatorBasics/StatesHandler.h>
#include <Utils/GeometricDerivatives/NormalModesContainer.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Settings.h>
#include <Utils/Technical/CloneInterface.h>
#include <cmath>
#include <vector>

#ifdef MPI_PARALLEL
  #include "VSCF/VSCF.h"
#endif

namespace Scine {
namespace Colibri {

class PesFromSPs : public Utils::CloneInterface<PesFromSPs, PES> {
 public:
  // Types definition
  using Base = PES;
  PesFromSPs() = default;
  ~PesFromSPs() override = default;
  PesFromSPs(VibrationalParameters& vibParms);
  PesFromSPs(const PesFromSPs& rhs);
  PesFromSPs& operator=(const PesFromSPs& rhs);
  double getPES() const override {
    return 0.0;
  };

 private:
  double getPES(double q, int mode) const override;
  double getPES(double qi, double qj, int modei, int modej) const override;
  double getPES(double qi, double qj, double qk, int modei, int modej, int modek) const override;
  double getHarmonicPES(double /*qi*/, int /*modei*/) const override {
    return 0.0;
  };

  static std::vector<double> readHarmFreqs(const std::string& filename, int numModes);
  int readSinglePoints(const std::string& filename);
  void parseLineOfSinglePoint(std::string lineString);
  static int convertToIntKey(double value);
  static double round(double value);

  double refEnergy_;

  mutable calculatedSPs<1> calc1BSPs_;
  mutable calculatedSPs<2> calc2BSPs_;
  mutable calculatedSPs<3> calc3BSPs_;

  friend class VSCF;
};

} // namespace Colibri
} // namespace Scine

#endif
