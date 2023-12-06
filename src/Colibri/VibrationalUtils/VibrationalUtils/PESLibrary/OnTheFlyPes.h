/**
 * @file OnTheFlyPes.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef OTF_PES_H
#define OTF_PES_H

#include "PES.h"
#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/CalculatorBasics/StatesHandler.h>
#include <Utils/GeometricDerivatives/NormalModesContainer.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Settings.h>
#include <Utils/Technical/CloneInterface.h>
#include <VibrationalUtils/VibrationalParameters.h>
#include <cmath>
#include <vector>

#ifdef MPI_PARALLEL
  #include "VSCF/VSCF.h"
#endif

namespace Scine {
namespace Colibri {

class OTFPes : public Utils::CloneInterface<OTFPes, PES> {
 public:
  // Types definition
  using GeometryType = Eigen::Matrix<double, Eigen::Dynamic, 3>;
  using HessianType = Eigen::MatrixXd;
  using NormalModeContainer = Utils::NormalModesContainer;
  using Base = PES;
  using Base::referenceGeometry_;
  OTFPes() = default;
  ~OTFPes() override = default;
  OTFPes(VibrationalParameters& vibParms);
  OTFPes(const OTFPes& rhs);
  OTFPes& operator=(const OTFPes& rhs);
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
  double getEnergy(const Utils::AtomCollection& geom) const;

  int getNumSPsForVSCF() const override;
  void setNumSPsForVSCF(int num) override;

  void setStartGeometry(std::string filename);

  static Utils::Settings applyVeryThightProfile(Utils::Settings& settings);

  Utils::AtomCollection getGeometryFromDisplacement(std::vector<int>&& idxModes, std::vector<double>&& steps) const;

  double refEnergy_;
  mutable int numSPsForVSCF_;

  static int convertToIntKey(double value);
  static double round(double value);

  Utils::NormalModesContainer normalModes_;

  /* Pointer to a Core calculator */
  std::shared_ptr<Core::Calculator> pesCalculator_;

  mutable calculatedSPs<1> calc1BSPs_;
  mutable calculatedSPs<2> calc2BSPs_;
  mutable calculatedSPs<3> calc3BSPs_;

  friend class VSCF;
};

} // namespace Colibri
} // namespace Scine

#endif
