/**
 * @file VCI.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef VCI_H
#define VCI_H
#include "VibrationalCI.h"
#include <VSCF/VibrationalSCF.h>
#include <VibrationalUtils/VibrationalParameters.h>
#include <memory>

namespace Scine {
namespace Colibri {

class VCI : public VibrationalCI {
  /**
   * @brief Implementation of the VCI algorithm.
   * Needs a prior VSCF calculation to get modal basis or a FCIDUMP file to read
   * from
   */
 public:
  using Base = VibrationalCI;
  VCI(const VibrationalParameters& parms, const std::shared_ptr<VibrationalSCF>& vscf);
  VCI(const VibrationalParameters& parms);
  ~VCI() override;

  int generateStatesFromGS() override;
  int generateStatesFromDet() override;
  double fillHamiltonianMat() override;
  double energyOfState(std::vector<int> state) override;
  void doDavidson() override;

 private:
  void generateAllStatesWithEx(std::vector<std::vector<int>>& exCombs);
  void generateAllStatesWithChange(std::vector<std::vector<int>>& exCombs);
  void generateAllExCombinations(std::vector<std::vector<int>>& exCombs, int possEx[], int chosen[], int numEle,
                                 int index, int start);
  void generateAllChangeCombinations(std::vector<std::vector<int>>& exCombs, int possEx[], int chosen[], int numEle,
                                     int index, int start);
  void parseLineFromFcidump(std::string line_string);

  tensor<3> oneBody_;
  tensor<6> twoBody_;
  tensor<9> threeBody_;
};

} // namespace Colibri
} // namespace Scine

#endif
