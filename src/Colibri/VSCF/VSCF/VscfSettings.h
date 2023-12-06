/**
 * @file VscfSettings.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef VSCF_SETTINGS_H
#define VSCF_SETTINGS_H
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingPopulator.h>

namespace Scine {
namespace Colibri {

class VscfSettings : public Utils::Settings {
 public:
  explicit VscfSettings(const std::string& name) : Utils::Settings(name) {
    Utils::UniversalSettings::FileDescriptor parameterFilePath("The path to the parameter file");
    parameterFilePath.setDefaultValue("");
    _fields.push_back("parameter_file", std::move(parameterFilePath));
    resetToDefaults();
  }
};
} // namespace Colibri
} // namespace Scine
#endif
