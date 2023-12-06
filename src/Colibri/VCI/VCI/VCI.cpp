/**
 * @file VCI.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "VCI.h"
#include "Davidson.h"
#include "Utils/Constants.h"
#include <omp.h>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

namespace Scine::Colibri {

// constructor for VCI afte VSCF
VCI::VCI(const VibrationalParameters& parms, const std::shared_ptr<VibrationalSCF>& vscf) : Base(parms) {
  // Get integrals from previous vscf calculation
  oneBody_ = vscf->getOneBodyModals();
  twoBody_ = vscf->getTwoBodyModals();
  threeBody_ = vscf->getThreeBodyModals();
}

// Read in integrals from FCIDUMP
VCI::VCI(const VibrationalParameters& parms) : Base(parms) {
  // Open the integral file and do a loop over this file
  if (!boost::filesystem::exists(fcidumpFile_)) {
    throw std::runtime_error("integral file " + fcidumpFile_ + " does not exist\n");
  }
  std::ifstream int_file;
  std::string line_string;
  int_file.open(fcidumpFile_.c_str());
  while (getline(int_file, line_string)) {
    // Extract the data
    parseLineFromFcidump(line_string);
  }
  int_file.close();
  std::cout << "Finished reading in FCIDUMP" << std::endl;
}

VCI::~VCI() = default;

void VCI::generateAllExCombinations(std::vector<std::vector<int>>& exCombs, int possEx[], int chosen[], int numEle,
                                    int index, int start) {
  if (index == numEle) {
    std::vector<int> exComb;
    exComb.reserve(numEle);
    for (int i = 0; i < numEle; i++) {
      exComb.push_back(possEx[chosen[i]]);
    }
    int exQ = std::accumulate(exComb.begin(), exComb.end(), 0);
    if (exQ <= totEx_) {
      bool valid = true;
      for (int i = 0; i < numEle; i++) {
        if (exComb[i] >= nMax_[i]) {
          valid = false;
        }
      }
      if (valid) {
        exCombs.push_back(exComb);
      }
    }
    return;
  }

  for (int i = start; i < *std::max_element(nMax_.begin(), nMax_.end()) - 1; i++) {
    chosen[index] = i;
    generateAllExCombinations(exCombs, possEx, chosen, numEle, index + 1, i);
  }
}

void VCI::generateAllChangeCombinations(std::vector<std::vector<int>>& exCombs, int possEx[], int chosen[], int numEle,
                                        int index, int start) {
  if (index == numEle) {
    std::vector<int> exComb;
    std::vector<int> exCombAbs;
    for (int i = 0; i < numEle; i++) {
      exComb.push_back(possEx[chosen[i]]);
      exCombAbs.push_back(abs(possEx[chosen[i]]));
    }
    int exQ = std::accumulate(exCombAbs.begin(), exCombAbs.end(), 0);
    if (exQ <= totEx_) {
      bool valid = true;
      for (int i = 0; i < numEle; i++) {
        if (exComb[i] >= nMax_[i]) {
          valid = false;
        }
      }
      if (valid) {
        exCombs.push_back(exComb);
      }
    }
    return;
  }

  for (int i = start; i < 2 * (*std::max_element(nMax_.begin(), nMax_.end()) - 1); i++) {
    chosen[index] = i;
    generateAllChangeCombinations(exCombs, possEx, chosen, numEle, index + 1, i);
  }
}

void VCI::generateAllStatesWithEx(std::vector<std::vector<int>>& exCombs) {
  for (auto state : exCombs) {
    for (int i = state.size(); i < nModes_; i++) {
      state.push_back(0);
    }
    std::sort(state.begin(), state.end());
    do {
      std::vector<int> stateSorted = state;
      std::reverse(stateSorted.begin(), stateSorted.end());
      states_.push_back(stateSorted);
    } while (std::next_permutation(state.begin(), state.end()));
  }
}

void VCI::generateAllStatesWithChange(std::vector<std::vector<int>>& exCombs) {
  std::sort(exCombs.begin(), exCombs.end());
  // Loop over all possible (de)excitation combinations
  for (auto changes : exCombs) {
    std::vector<bool> selectedModes(nModes_);
    std::fill(selectedModes.begin(), selectedModes.begin() + changes.size(), true);
    std::sort(changes.begin(), changes.end());
    // Loop over all possible change permutatoins
    do {
      std::vector<int> modes;
      for (int i = 0; i < nModes_; ++i) {
        if (selectedModes[i]) {
          modes.push_back(i);
        }
      }
      // Loop over all permutations of that mode combination
      do {
        // Generate new state from refDet with those changes
        std::vector<int> state = refDet_;
        int i = 0;
        for (int& mode : modes) {
          state[mode] = state[mode] + changes[i];
          i++;
        }
        // Check if the generated state fullfills criteria
        bool isAccepted = true;
        for (int n = 0; n < nModes_; n++) {
          if (state[n] > occupMax_[n]) {
            isAccepted = false;
          }
          if (state[n] < occupMin_[n]) {
            isAccepted = false;
          }
        }
        if (isAccepted) {
          states_.push_back(state);
        }
      } while (std::next_permutation(changes.begin(), changes.end()));
    } while (std::prev_permutation(selectedModes.begin(), selectedModes.end()));
  }
}

int VCI::generateStatesFromGS() {
  states_.clear();
  // First: ground state
  std::vector<int> gs(nModes_, 0);
  states_.push_back(gs);
  // Second: all excited states for nmex simultaneously excited modes
  for (int nmex = 1; nmex <= exModes_; nmex++) {
    std::vector<std::vector<int>> exCombs;
    int possEx[*std::max_element(nMax_.begin(), nMax_.end()) - 1]; // array with all possible excitation quanta
    for (int i = 1; i < *std::max_element(nMax_.begin(), nMax_.end()); i++) {
      possEx[i - 1] = i;
    }
    int chosen[nmex]; // helper array for finding all combinations
    generateAllExCombinations(exCombs, possEx, chosen, nmex, 0,
                              0);     // all possible excitations combinations are generated
    generateAllStatesWithEx(exCombs); // for all excitation combinations, all possible determinant
                                      // permutations are generated
  }
  return states_.size();
}

int VCI::generateStatesFromDet() {
  states_.clear();
  // First: reference state
  states_.push_back(refDet_);
  // Create helper arrays for (de)excitation generator
  int maxNmax = *std::max_element(nMax_.begin(), nMax_.end());
  int possEx[2 * (maxNmax - 1)]; // array with all possible excitation quanta
  for (int i = 1; i < maxNmax; i++) {
    possEx[i - 1] = -(maxNmax - i);
    possEx[i - 2 + maxNmax] = i;
  }
  // Second: all (de)excited states for nmex simultaneously (de)excited modes
  for (int nmex = 1; nmex <= exModes_; nmex++) {
    std::vector<std::vector<int>> exCombs;
    int chosen[nmex]; // helper array for finding all combinations
    generateAllChangeCombinations(exCombs, possEx, chosen, nmex, 0,
                                  0);     // all possible combinations of (de)excitations are generated
    generateAllStatesWithChange(exCombs); // for all (de)excitation combinations, all possible
                                          // determinant permutations are generated
  }
  return states_.size();
}

double VCI::fillHamiltonianMat() {
  hamiltonianMatrix_.setZero();
#pragma omp parallel for collapse(2) schedule(dynamic)
  for (std::size_t a = 0; a < states_.size(); a++) {
    for (std::size_t b = 0; b < states_.size(); b++) {
      std::vector<int> state_a = states_[a];
      std::vector<int> state_b = states_[b];
      std::vector<int> delta(nModes_, 0);
      for (std::vector<int>::size_type i = 0; i != nModes_; i++) {
        delta[i] = std::abs(state_a[i] - state_b[i]);
      }
      int ss = std::accumulate(delta.begin(), delta.end(), 0);
      for (int n = 0; n < nModes_; n++) {
        if (ss - delta[n] == 0) {
          std::array<int, 3> inds = {n, state_a[n], state_b[n]};
#pragma omp critical
          { hamiltonianMatrix_(a, b) = hamiltonianMatrix_(a, b) + oneBody_[inds]; }
        }
        if (potentialOrder_ > 1) {
          for (int m = n + 1; m < nModes_; m++) {
            if (ss - delta[n] - delta[m] == 0) {
              std::array<int, 6> inds = {n, state_a[n], state_b[n], m, state_a[m], state_b[m]};
#pragma omp critical
              { hamiltonianMatrix_(a, b) = hamiltonianMatrix_(a, b) + twoBody_[inds]; }
            }
            if (potentialOrder_ > 2) {
              for (int l = m + 1; l < nModes_; l++) {
                if (ss - delta[n] - delta[m] - delta[l] == 0) {
                  std::array<int, 9> inds = {n,          state_a[n], state_b[n], m,         state_a[m],
                                             state_b[m], l,          state_a[l], state_b[l]};
#pragma omp critical
                  { hamiltonianMatrix_(a, b) = hamiltonianMatrix_(a, b) + threeBody_[inds]; }
                }
              }
            }
          }
        }
      }
    }
  }
  return hamiltonianMatrix_(0, 0);
}

double VCI::energyOfState(std::vector<int> state) {
  double energy = 0;
  for (int n = 0; n < nModes_; n++) {
    std::array<int, 3> inds = {n, state[n], state[n]};
    energy += oneBody_[inds];
    if (potentialOrder_ > 1) {
      for (int m = n + 1; m < nModes_; m++) {
        std::array<int, 6> inds = {n, state[n], state[n], m, state[m], state[m]};
        energy += twoBody_[inds];
        if (potentialOrder_ > 2) {
          for (int l = m + 1; l < nModes_; l++) {
            std::array<int, 9> inds = {n, state[n], state[n], m, state[m], state[m], l, state[l], state[l]};
            energy += threeBody_[inds];
          }
        }
      }
    }
  }
  return energy;
}

void VCI::doDavidson() {
  // Straightforward evaluator for applying the hamiltonian
  class Evaluator {
   private:
    Eigen::MatrixXd A_;

   public:
    Evaluator(Eigen::MatrixXd A) : A_(std::move(A)) {
    }
    Eigen::MatrixXd evaluate(const Eigen::MatrixXd& P) {
      return A_ * P;
    }
  };

  int numberOfEigenpairs = numStatesToCalculate_;
  int maxSubspaceDim = 5 * numberOfEigenpairs;

  if (maxSubspaceDim > numStates_) {
    std::cout << "ERROR! The maximum subspace dimension is larger than the "
                 "number of states!"
              << std::endl;
    std::cout << "In this case it makes no sense to do a Davidson calculation "
                 "--> Aborting!"
              << std::endl;
    exit(1);
  }

  Davidson davidson(numberOfEigenpairs, maxSubspaceDim, numStates_);
  davidson.setMaxIterations(1000);
  davidson.setConvThresh(1e-10);
  davidson.setDiagonalPreconditioner(hamiltonianMatrix_.diagonal().asDiagonal());
  davidson.setGuess(Eigen::MatrixXd::Random(numStates_, numberOfEigenpairs));

  Evaluator evaluator(hamiltonianMatrix_);
  davidson.compute<Evaluator>(evaluator, true);

  eVals_ = davidson.getEigvals();
  eVecs_ = davidson.getEigvecs();
}

//
// +---------------------+
//   PARSE_INTEGRAL_NMODE
// +---------------------+
//
// Reads the FF from an input file and creates the Hamiltonian for the n-mode
// representation. The potential is given in the following form: i-n i-m float_1
// --> 1-body term (the first index is the mode, the second the basis set) i-n
// i-m j-p j-q   float_2   --> 2-body term (again, the indexes of the mode are
// the same.
//

void VCI::parseLineFromFcidump(std::string line_string) {
  std::vector<std::string> line_splitted;
  // Trim leading and final spaces in the string.
  line_string.erase(line_string.begin(),
                    std::find_if(line_string.begin(), line_string.end(), [&](int ch) { return std::isspace(ch) == 0; }));
  line_string.erase(
      std::find_if(line_string.rbegin(), line_string.rend(), [&](int ch) { return std::isspace(ch) == 0; }).base(),
      line_string.end());

  // Split the string into i-n, i-m, float_1 or i-n, i-m, j-p, j-q, float_2 or
  // ............
  boost::split(line_splitted, line_string, boost::is_any_of(" "), boost::token_compress_on);
  double integral_term =
      std::atof(line_splitted[line_splitted.size() - 1].c_str()) / Utils::Constants::invCentimeter_per_hartree;

  // Check the degree of coupling
  if (!(line_splitted.size() == 3 || line_splitted.size() == 5 || line_splitted.size() == 7)) {
    if (line_splitted.size() >= 9) {
      throw std::runtime_error("VCI with more than three-body terms NYI");
    }
    throw std::runtime_error("FCIDUMP seems to contain lines not adhering to correct format");
  }
  if (potentialOrder_ < 1) {
    throw std::runtime_error("No 0th order VCI possible");
  }

  std::vector<std::string> term1, term2;
  boost::split(term1, line_splitted[0], boost::is_any_of("-"));
  boost::split(term2, line_splitted[1], boost::is_any_of("-"));
  // Checks data consistency
  assert(term1.size() == 2 && term2.size() == 2 && term1[0] == term2[0]);

  if (line_splitted.size() == 3) {
    // One-body term
    std::array<int, 3> ind = {std::stoi(term1[0]) - 1, std::stoi(term1[1]), std::stoi(term2[1])};
    oneBody_[ind] = integral_term;
    return;
  }

  // Now we get to higher-order terms
  if (potentialOrder_ < 2) {
    throw std::runtime_error("FCIDUMP integral file contains higher order "
                             "integral terms than set mode-coupling order");
  }

  std::vector<std::string> term3, term4;
  boost::split(term3, line_splitted[2], boost::is_any_of("-"));
  boost::split(term4, line_splitted[3], boost::is_any_of("-"));
  // Checks data consistency
  assert(term3.size() == 2 && term4.size() == 2 && term3[0] == term4[0]);

  if (line_splitted.size() == 5) {
    // Two-body term
    std::array<int, 6> ind = {std::stoi(term1[0]) - 1, std::stoi(term1[1]), std::stoi(term2[1]),
                              std::stoi(term3[0]) - 1, std::stoi(term3[1]), std::stoi(term4[1])};
    twoBody_[ind] = integral_term;
    return;
  }

  // If we reach here, we have a 3-body term
  if (potentialOrder_ != 3) {
    throw std::runtime_error("FCIDUMP integral file contains higher order "
                             "integral terms than set mode-coupling order");
  }

  std::vector<std::string> term5, term6;
  boost::split(term5, line_splitted[4], boost::is_any_of("-"));
  boost::split(term6, line_splitted[5], boost::is_any_of("-"));
  // Checks data consistency
  assert(term5.size() == 2 && term6.size() == 2 && term5[0] == term6[0]);

  // Three-body term
  std::array<int, 9> ind = {std::stoi(term1[0]) - 1, std::stoi(term1[1]), std::stoi(term2[1]),
                            std::stoi(term3[0]) - 1, std::stoi(term3[1]), std::stoi(term4[1]),
                            std::stoi(term5[0]) - 1, std::stoi(term5[1]), std::stoi(term6[1])};
  threeBody_[ind] = integral_term;
}

} // namespace Scine::Colibri
