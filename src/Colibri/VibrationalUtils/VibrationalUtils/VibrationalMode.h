/**
 * @file VibrationalMode.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#ifndef VIB_MODE_H
#define VIB_MODE_H

namespace Scine {
namespace Colibri {

/**
 * @brief Representation of a vibrational mode
 * @param freq frequency of the mode
 * @param nquantum quantum number of the mode
 */
class VibrationalMode {
 public:
  VibrationalMode() = default;
  VibrationalMode(double freq, int nquantum) : freq_(freq), nquantum_(nquantum){};
  ~VibrationalMode() = default;

 protected:
  double freq_;
  int nquantum_;
};

} // namespace Colibri
} // namespace Scine

#endif