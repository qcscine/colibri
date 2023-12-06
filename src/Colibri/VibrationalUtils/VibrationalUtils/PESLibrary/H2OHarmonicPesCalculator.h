/**
 * @file H20HarmonicPesCalculator.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef H2O_HARMONIC_PES_CALCULATOR_H
#define H2O_HARMONIC_PES_CALCULATOR_H

#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/CalculatorBasics/StatesHandler.h>
#include <Utils/Settings.h>
#include <Utils/Technical/CloneInterface.h>
#include <Eigen/Dense>

namespace Scine {
namespace Colibri {

class H2OHarmonicPesCalculator : public Utils::CloneInterface<Colibri::H2OHarmonicPesCalculator, Core::Calculator> {
 public:
  static constexpr int nAtoms = 3;
  using GeometryType = Eigen::Matrix<double, nAtoms, 3>;
  using DistanceMatrixType = Eigen::Matrix<double, nAtoms*(nAtoms - 1) / 2, 1>;
  using HessianMatrixType = Eigen::Matrix<double, 9, 9>;

  /* Default constructor and destructor */
  H2OHarmonicPesCalculator();
  ~H2OHarmonicPesCalculator() final = default;

  /* Copy constructor */
  H2OHarmonicPesCalculator(const H2OHarmonicPesCalculator& rhs);

  std::string name() const final;

  Utils::Settings& settings() final;

  const Utils::Settings& settings() const final;

  Utils::Results& results() final;

  const Utils::Results& results() const final;

  const Utils::Results& calculate(std::string description) final;

  void setRequiredProperties(const Utils::PropertyList& requiredProperties) final;

  Utils::PropertyList getRequiredProperties() const final;

  Utils::PropertyList possibleProperties() const override;

  void setStructure(const Utils::AtomCollection& structure) final;

  std::unique_ptr<Utils::AtomCollection> getStructure() const final;

  void modifyPositions(Utils::PositionCollection newPositions) final;

  const Utils::PositionCollection& getPositions() const final;

  bool supportsMethodFamily(const std::string& methodFamily) const final;

  bool allowsPythonGILRelease() const override {
    return true;
  };

  // Public constants
  static constexpr const char* model = "N2OPes";

 protected:
  void loadState(std::shared_ptr<Core::State> state) final{};
  std::shared_ptr<Core::State> getState() const final {
    return std::make_shared<Core::State>();
  };
  Utils::Results results_;
  Utils::PropertyList requiredProperties_;
  std::unique_ptr<Utils::Settings> settings_;
  Utils::ElementTypeCollection elements_;
  Utils::PositionCollection positions_, referencePositions_;

 private:
  // Private members
  HessianMatrixType Hessian_;
};

} // namespace Colibri
} // namespace Scine

#endif
