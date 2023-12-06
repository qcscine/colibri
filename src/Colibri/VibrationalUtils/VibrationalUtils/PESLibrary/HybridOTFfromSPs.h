/**
 * @file HybridOTFfromSPs.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef HYBRID_PES_H
#define HYBRID_PES_H

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

class HybridPes : public Utils::CloneInterface<HybridPes, PES> {
 public:
  // Types definition
  using GeometryType = Eigen::Matrix<double, Eigen::Dynamic, 3>;
  using HessianType = Eigen::MatrixXd;
  using NormalModeContainer = Utils::NormalModesContainer;
  using Base = PES;
  using Base::referenceGeometry_;
  HybridPes() = default;
  ~HybridPes() override = default;
  HybridPes(VibrationalParameters& vibParms);
  HybridPes(const HybridPes& rhs);
  HybridPes& operator=(const HybridPes& rhs);
  double getPES() const override {
    return refEnergy_;
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

  // Stuff needed for OTF calcs
  void setGeometry(std::string filename);

  Utils::AtomCollection getGeometryFromDisplacement(std::vector<int>&& idxModes, std::vector<double>&& steps) const;

  Utils::NormalModesContainer normalModes_;

  /* Pointer to a Core calculator */
  std::shared_ptr<Core::Calculator> pesCalculator_;

  // Stuff needed for precalculatedSPs
  static std::vector<double> readHarmFreqs(const std::string& filename, int numModes);
  int readSinglePoints(const std::string& filename);
  void parseLineOfSinglePoint(std::string lineString);
  static int convertToIntKey(double value);
  static double round(double value);

  double refEnergy_;
  mutable int numSPsForVSCF_;

  mutable calculatedSPs<1> calc1BSPs_;
  mutable calculatedSPs<2> calc2BSPs_;
  mutable calculatedSPs<3> calc3BSPs_;

  friend class VSCF;
};

} // namespace Colibri
} // namespace Scine

#endif
