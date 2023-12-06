/**
 * @file VscfCalculator.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "VscfCalculator.h"
#include "VSCF.h"
#include "VscfSettings.h"
#include <cmath>

#ifdef MPI_PARALLEL
  #include <mpi.h>
#endif

namespace Scine::Colibri {
// TODO redo calculator

template class VSCFCalculator<PesQFF, DistributedGaussians>;
template class VSCFCalculator<PesQFF, DVR>;
template class VSCFCalculator<PesChristoffel, DistributedGaussians>;
template class VSCFCalculator<PesChristoffel, DVR>;

template<typename pesType, typename basisType>
VSCFCalculator<pesType, basisType>::VSCFCalculator() {
  settings_ = std::make_unique<VscfSettings>("vscf");
}

template<typename pesType, typename basisType>
VSCFCalculator<pesType, basisType>::VSCFCalculator(const VSCFCalculator<pesType, basisType>& rhs)
  : VSCFCalculator<pesType, basisType>() {
  *(this->settings_) = rhs.settings();
  setStructure(*rhs.getStructure());
}

template<typename pesType, typename basisType>
std::unique_ptr<Utils::AtomCollection> VSCFCalculator<pesType, basisType>::getStructure() const {
  return std::make_unique<Utils::AtomCollection>(elementTypes_, positions_);
}

template<typename pesType, typename basisType>
Utils::Settings& VSCFCalculator<pesType, basisType>::settings() {
  return *settings_;
}

template<typename pesType, typename basisType>
const Utils::Settings& VSCFCalculator<pesType, basisType>::settings() const {
  return *settings_;
}

template<typename pesType, typename basisType>
Utils::Results& VSCFCalculator<pesType, basisType>::results() {
  return results_;
}

template<typename pesType, typename basisType>
const Utils::Results& VSCFCalculator<pesType, basisType>::results() const {
  return results_;
}

template<typename pesType, typename basisType>
const Utils::Results& VSCFCalculator<pesType, basisType>::calculate(std::string description) {
  vibParms_.setParameters(description);
  pes_ = std::make_shared<pesType>(vibParms_);
  basis_ = std::make_shared<basisType>(vibParms_, pes_);
  vscf_ = std::make_unique<VSCF>(vibParms_, basis_);
  results_ = Utils::Results();
  int id = 0;
#ifdef MPI_PARALLEL
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
#endif
  double energy = NAN;
  if (id == 0) {
    energy = vscf_->SCF();
  }
  results_.set<Utils::Property::Energy>(energy);
  results_.set<Utils::Property::SuccessfulCalculation>(true);
  return results_;
}

template<typename pesType, typename basisType>
void VSCFCalculator<pesType, basisType>::setRequiredProperties(const Utils::PropertyList& requiredProperties) {
  requiredProperties_ = requiredProperties;
}

template<typename pesType, typename basisType>
Utils::PropertyList VSCFCalculator<pesType, basisType>::getRequiredProperties() const {
  return requiredProperties_;
}

template<typename pesType, typename basisType>
Utils::PropertyList VSCFCalculator<pesType, basisType>::possibleProperties() const {
  return Utils::Property::Energy;
}

template<typename pesType, typename basisType>
void VSCFCalculator<pesType, basisType>::setStructure(const Utils::AtomCollection& structure) {
}

template<typename pesType, typename basisType>
VSCFCalculator<pesType, basisType>::~VSCFCalculator() = default;

template<typename pesType, typename basisType>
bool VSCFCalculator<pesType, basisType>::supportsMethodFamily(const std::string& methodFamily) const {
  return methodFamily == name();
}

template<typename pesType, typename basisType>
std::string VSCFCalculator<pesType, basisType>::name() const {
  return (*this).model;
}

template<typename pesType, typename basisType>
void VSCFCalculator<pesType, basisType>::modifyPositions(Utils::PositionCollection newPositions) {
}

template<typename pesType, typename basisType>
const Utils::PositionCollection& VSCFCalculator<pesType, basisType>::getPositions() const {
  return positions_;
}

} // namespace Scine::Colibri
