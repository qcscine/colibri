/**
 * @file ColibriModule.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#include "ColibriModule.h"
#include <Core/DerivedModule.h>
#include <VSCF/VscfCalculator.h>
#include <VibrationalUtils/PESLibrary/PesChristoffel.h>
#include <VibrationalUtils/PESLibrary/PesQFF.h>
#include <VibrationalUtils/VibrationalBases/DVR.h>
#include <VibrationalUtils/VibrationalBases/DistributedGaussians.h>
#include <boost/mpl/list.hpp>
#include <boost/mpl/map.hpp>

namespace Scine::Colibri {

std::string ColibriModule::name() const noexcept {
  return "Colibri";
}

using InterfaceModelMap = boost::mpl::map<boost::mpl::pair<
    Core::Calculator, boost::mpl::list<Colibri::VSCFCalculator<Colibri::PesQFF, Colibri::DistributedGaussians>,
                                       Colibri::VSCFCalculator<Colibri::PesQFF, Colibri::DVR>, Colibri::VSCFCalculator<Colibri::PesChristoffel, Colibri::DistributedGaussians>,
                                       Colibri::VSCFCalculator<Colibri::PesChristoffel, Colibri::DVR>>>>;

boost::any ColibriModule::get(const std::string& interface, const std::string& model) const {
  boost::any resolved = Scine::Core::DerivedModule::resolve<InterfaceModelMap>(interface, model);
  // Throw an exception if we could not match an interface or model
  if (resolved.empty()) {
    throw Scine::Core::ClassNotImplementedError();
  }
  return resolved;
}

bool ColibriModule::has(const std::string& interface, const std::string& model) const noexcept {
  return Core::DerivedModule::has<InterfaceModelMap>(interface, model);
}

std::vector<std::string> ColibriModule::announceInterfaces() const noexcept {
  return Core::DerivedModule::announceInterfaces<InterfaceModelMap>();
}

std::vector<std::string> ColibriModule::announceModels(const std::string& interface) const noexcept {
  return Core::DerivedModule::announceModels<InterfaceModelMap>(interface);
}

std::shared_ptr<Core::Module> ColibriModule::make() {
  return std::make_shared<ColibriModule>();
}

std::vector<std::shared_ptr<Scine::Core::Module>> moduleFactory() {
  return {ColibriModule::make()};
}

} // namespace Scine::Colibri
