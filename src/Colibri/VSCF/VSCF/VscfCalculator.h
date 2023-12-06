/**
 * @file VscfCalculator.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef VSCF_CALCULATOR_H
#define VSCF_CALCULATOR_H

#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/CalculatorBasics/StatesHandler.h>
#include <Utils/Settings.h>
#include <Utils/Technical/CloneInterface.h>
#include <VSCF/VSCF.h>
#include <VibrationalUtils/PESLibrary/PES.h>
#include <VibrationalUtils/PESLibrary/PesChristoffel.h>
#include <VibrationalUtils/PESLibrary/PesQFF.h>
#include <VibrationalUtils/VibrationalBases/DVR.h>
#include <VibrationalUtils/VibrationalBases/DistributedGaussians.h>
#include <VibrationalUtils/VibrationalBases/ModalBasis.h>
#include <VibrationalUtils/VibrationalParameters.h>

namespace Scine {
namespace Colibri {
// TODO include NMA
// Model names
class PesQFF;
class DistributedGaussians;
class DVR;
class PesChristoffel;

template<class PesType, class BasisType>
struct NameTrait {
  static constexpr const char* name = "";
};
template<>
struct NameTrait<PesQFF, DistributedGaussians> {
  static constexpr const char* name = "PesQffDistGauss";
};
template<>
struct NameTrait<PesQFF, DVR> {
  static constexpr const char* name = "PesQffDVR";
};
template<>
struct NameTrait<PesChristoffel, DistributedGaussians> {
  static constexpr const char* name = "PesChristoffelDistGauss";
};
template<>
struct NameTrait<PesChristoffel, DVR> {
  static constexpr const char* name = "PesChristoffelDVR";
};

template<typename pesType, typename basisType>
class VSCFCalculator : public Utils::CloneInterface<Colibri::VSCFCalculator<pesType, basisType>, Core::Calculator> {
 public:
  VSCFCalculator();
  VSCFCalculator(const VSCFCalculator<pesType, basisType>& rhs);
  ~VSCFCalculator() override;
  std::string name() const final;
  static constexpr const char* model = NameTrait<pesType, basisType>::name;
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

 protected:
  void loadState(std::shared_ptr<Core::State> state) final{};
  std::shared_ptr<Core::State> getState() const final {
    return std::make_shared<Core::State>();
  };

  Utils::Results results_;
  Utils::PropertyList requiredProperties_;
  std::unique_ptr<Utils::Settings> settings_;
  Utils::ElementTypeCollection elementTypes_;
  Utils::PositionCollection positions_;

 private:
  VibrationalParameters vibParms_;
  std::shared_ptr<PES> pes_;
  std::shared_ptr<ModalBasis> basis_;
  std::unique_ptr<VSCF> vscf_;
};

} // namespace Colibri
} // namespace Scine
#endif
