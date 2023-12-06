/**
 * @file IntegralDumper.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#ifndef INTEGRAL_DUMPER_H
#define INTEGRAL_DUMPER_H

#include <boost/functional/hash.hpp>
#include <unordered_map>

namespace Scine {
namespace Colibri {

class IntegralDumper {
 public:
  using modalOneBodyMap = typename std::unordered_map<std::array<int, 3>, double, boost::hash<std::array<int, 3>>>;
  using modalTwoBodyMap = typename std::unordered_map<std::array<int, 6>, double, boost::hash<std::array<int, 6>>>;
  using modalThreeBodyMap = typename std::unordered_map<std::array<int, 9>, double, boost::hash<std::array<int, 9>>>;

  IntegralDumper() = default;
  ~IntegralDumper() = default;

 protected:
  static void writeIntegrals(modalOneBodyMap& oneBody, modalTwoBodyMap& twoBody, const std::vector<int>& nMax, int nModes,
                             const std::string& filename, bool onlyDumpCoupledModes, const std::vector<int>& coupledModes);
  static void writeIntegrals(modalOneBodyMap& oneBody, modalTwoBodyMap& twoBody, modalThreeBodyMap& threeBody,
                             const std::vector<int>& nMax, int nModes, const std::string& filename,
                             bool onlyDumpCoupledModes, const std::vector<int>& coupledModes);
  static void writeIntegrals(modalOneBodyMap& oneBody, std::vector<int> nMax, int nModes, const std::string& filename,
                             bool onlyDumpCoupledModes, std::vector<int> coupledModes);

 private:
  static bool fileExists(const std::string& filename);
  static void writeOneBody(std::ofstream& fciDump, modalOneBodyMap& oneBody, int nModes, std::vector<int> nMax,
                           bool onlyDumpCoupledModes, std::vector<int> coupledModes);
  static void writeTwoBody(std::ofstream& fciDump, modalTwoBodyMap& twoBody, int nModes, std::vector<int> nMax,
                           bool onlyDumpCoupledModes, std::vector<int> coupledModes);
  static void writeThreeBody(std::ofstream& fciDump, modalThreeBodyMap& threeBody, int nModes, std::vector<int> nMax,
                             bool onlyDumpCoupledModes, std::vector<int> coupledModes);
};

} // namespace Colibri
} // namespace Scine

#endif
