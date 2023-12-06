/**
 * @file VibrationalParameters.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#ifndef VIB_PARMS_H
#define VIB_PARMS_H

#include "Tensor.h"
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <vector>

namespace Scine {
namespace Colibri {

class VibrationalParameters {
  /**
   * @brief Class storing the parameters for a frequency calculation
   *
   * @param numModes_ number of vibrational modes
   * @param coupledModes_ vector holding the indices of modes which should be
   * considered coupled
   * @param dumpOnlyCoupledModes_ bool to dump only coupled modes
   * @param nMax_ number of modals to be considered per mode
   * @param vscfIter_ maximum number of VSCF iterations
   * @param numqp_ number of quadrature points used for the integration
   * @param nModePotentialOrder_ order of potential in n-mode expansion
   * @param quantumNums_ vector holding the quantum number for each mode
   * @param ONVector_ occupation number vector of considered wave-function
   * @param modeFreqs_ vector holding the frequency of each mode
   * @param vscfEnTol_ tolerance for energy convergence in VSCF
   * @param vscfCoeffTol_ tolerance for coefficients in VSCF
   * @param twoBodyTol_ tolerance for two-body integrals
   * @param xyzFname_ path/filename of xyz file to read reference geometry from
   * @param hessianFname_ path/filename of Hessian file to read Hessian from
   * @param doTransitionStateOptimization if set to true, TS will be optimized,
   * otherwise regular optimisation
   * @param harmFreqFname_ path/filename of text file to read harmonic
   * frequencies from
   * @param refEnergy_ reference energy in hartree to use for PES from
   * precalculated singlepoints
   * @param singlePointsFname_ path/filename of text file to read precalculated
   * singlepoints from
   * @param fcidumpFname_ path/filename of FCIDUMP file to dump to/read from
   * @param coeffsFname_ path/filename of file to dump VSCF coefficients to
   * @param inputPESFname_ path/filename of integral file to read in PES from
   * @param inputPESType_ type of PES provided in the input file
   * @param primitiveBasisType_ type of primitve basis to use, either DG or DVR
   * @param taylorExpOrder_ order of the polynomial of the input PES
   * @param otfProgram_ program to perform on-the-fly calculation, same keyword
   * as for readuct
   * @param otfMethodFamily_ method family to perform on-the-fly calculation,
   * same keyword as for readuct
   * @param otfMethod_ specific method to perform on-the-fly calculation,
   * same keyword as for readuct
   * @param otfBasisSet_ basis set specifier to perform on-the-fly calculation
   * @param stepForNumericalDifferences_ step size for the numerical calculation
   * of differences for gradients or hessians
   * @param vciExModes_ number of Modes to be simultaneously excited in VCI
   * @param vciTotEx_ number of total excitations to be simultaneously excited
   * in VCI
   * @param vciDoStateSpecific_ whether to do standard (GS-based) or
   * state-specific VCI
   * @param vciRefDet_ determinant from which to generate states for
   * state-specific VCI
   * @param vciOccupMin_ minmal occupation of generated states for
   * state-specific VCI
   * @param vciOccupMax_ maximal occupation of generated states for
   * state-specific VCI
   * @param vciEVSolver_ eigensolver to be used (Eigen::SelfAdjointEigenSolver
   * or Davidson)
   * @param vciNumStates_ number of states that the user is interested in (only
   * those solved in Davidson)
   *
   */
 public:
  /**
   * @brief Construct a new Vibrational Parameters object
   *
   * @param infile path to input file in .ini format
   */
  VibrationalParameters() = default;
  ~VibrationalParameters() = default;

  void setParameters(const std::string& path);
  void setParametersN(const std::string& path);
  void setParametersAlsoVCI(const std::string& path);
  void setParametersOnlyVCI(const std::string& path);
  void setParametersOTF(const std::string& path);
  void printParameters();
  void printParametersN();
  void printParametersVSCF();
  void printParametersVCI();
  void printParametersOnlyVCI();
  void printParametersOTF();
  void printParametersFromSPs();
  void printParametersHybrid();
  void setPES(std::string integral_file);

  tensor<1> linForceConst_;
  tensor<2> quadraticForceConst_;
  tensor<3> cubicForceConst_;
  tensor<4> quarticForceConst_;

  template<std::size_t rank>
  tensor<rank> getForceConstants(const std::string& inp);
  template<std::size_t rank>
  void addForceConstants(tensor<rank>& K, const std::string& inp);

  // Alternatively for a more efficient storage:
  void setPesMc(std::string integral_file);
  tensorMC<1> oneModeCoeffs_;
  tensorMC<2> twoModeCoeffs_;
  tensorMC<3> threeModeCoeffs_;
  tensorMC<4> fourModeCoeffs_;
  tensorMC<5> fiveModeCoeffs_;
  tensorMC<6> sixModeCoeffs_;
  template<std::size_t mcorder>
  void addForceConstants(tensorMC<mcorder>& K, const std::string& inp);

  // Alternatively for an even more efficient storage:
  void setPesMCSP(std::string integral_file);
  tensorMCSP<1> oneModeTerms_;
  tensorMCSP<2> twoModeTerms_;
  tensorMCSP<3> threeModeTerms_;
  tensorMCSP<4> fourModeTerms_;
  tensorMCSP<5> fiveModeTerms_;
  tensorMCSP<6> sixModeTerms_;
  template<std::size_t mcorder>
  void addForceConstants(tensorMCSP<mcorder>& K, const std::string& inp);

  // Colibri parameters
  bool vciDoStateSpecific_, dumpOnlyCoupledModes_, doTransitionStateOptimization_;
  int numModes_, vscfIter_, numqp_, nModePotentialOrder_, taylorExpOrder_;
  int vciExModes_, vciTotEx_, vciNumStates_;
  double vscfEnTol_, vscfCoeffTol_, refEnergy_;
  double twoBodyTol_, stepForNumericalDifferences_;
  std::vector<int> numGridPoints_, ONVector_, vciRefDet_, vciOccupMin_, vciOccupMax_, coupledModes_, nMax_;
  std::vector<double> modeFreqs_, quantumNums_;
  std::string fcidumpFname_, inputPESFname_, xyzFname_, hessianFname_, harmFreqFname_, singlePointsFname_, coeffsFname_;
  std::string inputPESType_, primitiveBasisType_;
  std::string vciEVSolver_;
  std::string otfProgram_, otfMethodFamily_, otfMethod_, otfBasisSet_;

  setNbody<3> set3B_;
  void addSet3B(const std::string& inp);
};

} // namespace Colibri
} // namespace Scine
#endif