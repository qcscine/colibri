/**
 * @file PesCartesian.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#ifndef PES_CARTESIAN_H
#define PES_CARTESIAN_H

#include "PES.h"
#include "VibrationalUtils/Tensor.h"
#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/CalculatorBasics/StatesHandler.h>
#include <Utils/GeometricDerivatives/NormalModesContainer.h>
#include <Utils/Settings.h>
#include <Utils/Technical/CloneInterface.h>
#include <Eigen/Dense>

namespace Scine {
namespace Colibri {

/**
 * @brief Class representing a generic PES in Cartesian normal coordinates.
 *
 * This PES calculates the energy as a generic function of the Cartesian
 * displacements (i.e., any function V(x_1, y_1, ..., z_N)).
 *
 * @tparam EnergyComputer functor class representing the calculator used to
 * evaluate the function V(x_1, y_1, ..., z_N). It must implement a
 * setStructure method, a calculate method and have a nAtoms member.
 */
class PesCartesian : public Utils::CloneInterface<PesCartesian, PES> {
 public:
  // Types definition
  using GeometryType = Eigen::Matrix<double, Eigen::Dynamic, 3>;
  using HessianType = Eigen::MatrixXd;
  using NormalModeContainer = Utils::NormalModesContainer;
  using Base = PES;
  using Base::referenceGeometry_;
  /* Default constructor/descructor */
  PesCartesian() = delete;
  PesCartesian(std::shared_ptr<Core::Calculator> pesCalculator);
  PesCartesian(const PesCartesian& rhs);
  ~PesCartesian() override = default;
  /* Getter for the energy at the reference geometry */
  double getPES() const override;
  /* Getter for the energy for a one-mode displacement */
  double getPES(double qi, int modei) const override;
  /* Getter for the energy for a two-modes dispalcement */
  double getPES(double qi, double qj, int modei, int modej) const override;
  /* Print the normal modes */
  void printNormalModeInformation(std::ostream& stream) const;
  double printHarmonicZPVE() const;
  // TODO ALB This is not something of a generic PES.
  void setReferenceGeometry(tensor<1>& k1, tensor<2>& k2, tensor<3>& k3, tensor<4>& k4) override;
  double getHarmonicPES(double qi, int modei) const override;

 private:
  double getEnergy(const Utils::AtomCollection& geom) const;
  void calculateNormalModesAndHarmonicFrequencies();
  Utils::AtomCollection updateReferenceGeometryFromDisplacement(std::vector<int>&& idxModes, std::vector<double>&& steps) const;
  /* Pointer to a Core calculator */
  std::shared_ptr<Core::Calculator> pesCalculator_;
  HessianType hessian_;
  NormalModeContainer normalModesContainer_;
  std::vector<double> masses_;
  mutable Utils::AtomCollection currentGeometry_;
  /* Step-size for the numeric evaluation of the Hessian/gradient */
  static const double stepForNumericalDifferences;
};

} // namespace Colibri
} // namespace Scine

#endif
