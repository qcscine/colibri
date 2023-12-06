/**
 * @file VibrationalSCF.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "VibrationalSCF.h"
#include <sys/stat.h>
#include <Eigen/SparseCore>
#include <cmath>
#include <iomanip>

namespace Scine::Colibri {

VibrationalSCF::VibrationalSCF(const VibrationalParameters& vibParms) {
  nModes = vibParms.numModes_;
  maxIter = vibParms.vscfIter_;
  occupation = vibParms.ONVector_;
  integralDumpFile = vibParms.fcidumpFname_;
  coeffsDumpFile_ = vibParms.coeffsFname_;
  scfEnergyTolerance_ = vibParms.vscfEnTol_;
  scfCoeffTolerance_ = vibParms.vscfCoeffTol_;
  potentialOrder = vibParms.nModePotentialOrder_;
  integralTol_ = vibParms.twoBodyTol_;
  nMax_ = vibParms.nMax_;
  coupledModes_ = vibParms.coupledModes_;
  dumpOnlyCoupledModes_ = vibParms.dumpOnlyCoupledModes_;
  isConverged_ = false;
}

std::vector<std::pair<double, Eigen::VectorXd>> VibrationalSCF::sortEigenPairs(const Eigen::MatrixXd& eVecs,
                                                                               const Eigen::VectorXd& eVals) {
  int nBasis = eVals.size();
  std::vector<std::pair<double, Eigen::VectorXd>> evalPairs(nBasis);
  for (int i = 0; i < nBasis; i++) {
    evalPairs[i].first = eVals(i);
    evalPairs[i].second = eVecs.col(i);
  }
  std::sort(evalPairs.begin(), evalPairs.end(), [](auto a, auto b) { return a.first < b.first; });
  return evalPairs;
}

bool VibrationalSCF::evecConverged(std::vector<Eigen::VectorXd>& thisIter, std::vector<Eigen::VectorXd>& prevIter, double crit) {
  std::vector<double> vecDiffs(nModes);
  for (int i = 0; i < nModes; i++) {
    Eigen::VectorXd diff = thisIter[i] - prevIter[i];
    double norm = thisIter[i].transpose() * this->getOverlap(i) * thisIter[i];
    double vecDiff = diff.transpose() * this->getOverlap(i) * diff;
    vecDiffs[i] = vecDiff / norm;
  }
  auto res = std::max_element(std::begin(vecDiffs), std::end(vecDiffs));
  return (*res < crit);
}

double VibrationalSCF::SCF() {
  coeffs_.resize(nModes);
  evals_.resize(nModes);
  evecsIterations_.reserve(maxIter);
  evalsIterations_.reserve(maxIter);
  allModalPairs_.resize(nModes);
  coeffs_ = getStartingGuess();
  std::ofstream outfile;
  outfile.open(coeffsDumpFile_, std::ios_base::app);
  for (int mode = 0; mode < nModes; mode++) {
    outfile << "Starting guess coefficient vector for mode " << mode << " is: " << coeffs_[mode](0);
    for (int i = 1; i < coeffs_[mode].size(); i++) {
      outfile << ", " << coeffs_[mode](i);
    }
    outfile << std::endl;
  }
  std::cout << "STARTING SCF" << std::endl;
  double vscfEnergy = NAN;
  for (int vscfIt = 0; vscfIt < maxIter; vscfIt++) {
    std::vector<Eigen::VectorXd> allEvals;
    allEvals.resize(nModes);
    for (int mode = 0; mode < nModes; mode++) {
      int Nbasis = this->getBasisSize(mode);
      Eigen::MatrixXd F = Eigen::MatrixXd::Zero(Nbasis, Nbasis);
      if (potentialOrder == 2 || potentialOrder == 6 || potentialOrder == 7) {
        F = this->getMeanFieldOperator(mode, coeffs_);
      } else if (potentialOrder == 5 || potentialOrder == 1 || potentialOrder == 0) {
        F = this->getMeanFieldOperator(mode);
      } else if (potentialOrder == 3) {
        F = this->getMeanFieldOperator3(mode, coeffs_);
      } else {
        std::cout << "Potential order not known." << std::endl;
      }
      Eigen::MatrixXd pot = Eigen::MatrixXd::Zero(Nbasis, Nbasis);
      pot = F - this->getKineticOperator(mode);
      outfile << "Potential matrix of mode " << mode << " at iteration " << vscfIt << " is: " << pot(0, 0);
      for (int i = 1; i < Nbasis; i++) {
        outfile << ", " << pot(i, i);
      }
      outfile << std::endl;

      Eigen::MatrixXd S = Eigen::MatrixXd::Zero(Nbasis, Nbasis);
      S = getOverlap(mode);
      Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> scf;
      scf.compute(F, S);
      Eigen::VectorXd evals = scf.eigenvalues();
      Eigen::MatrixXd evecs = scf.eigenvectors();
      allEvals[mode] = evals;
      std::vector<std::pair<double, Eigen::VectorXd>> evalPairs(evals.size());
      evalPairs = sortEigenPairs(evecs, evals);
      int conf = occupation[mode];
      coeffs_[mode] = evalPairs[conf].second;
      evals_[mode] = evalPairs[conf].first;
      allModalPairs_[mode] = evalPairs;
      outfile << "The coefficient vector for mode " << mode << " at iteration " << vscfIt << " is: " << coeffs_[mode](0);
      for (int i = 1; i < coeffs_[mode].size(); i++) {
        outfile << ", " << coeffs_[mode](i);
      }
      outfile << std::endl;
    }
    evecsIterations_.push_back(coeffs_);
    if (potentialOrder == 2 || potentialOrder == 6 || potentialOrder == 7) {
      vscfEnergy = getEnergy(coeffs_, evals_);
    } else if (potentialOrder == 3) {
      vscfEnergy = getEnergy3(coeffs_, evals_);
    } else if (potentialOrder == 0 || potentialOrder == 1 || potentialOrder == 5) {
      vscfEnergy = getEnergy(evals_);
    }
    evalsIterations_.push_back(vscfEnergy);
    std::cout << "ITERATION: " << vscfIt << "   CONVERGED ENERGY: " << std::setprecision(10) << std::fixed << vscfEnergy
              << std::endl;
    if (vscfIt > 1 && std::abs(evalsIterations_.end()[-1] - evalsIterations_.end()[-2]) < scfEnergyTolerance_ &&
        this->evecConverged(evecsIterations_.end()[-1], evecsIterations_.end()[-2], scfCoeffTolerance_)) {
      std::cout << "CONVERGED." << std::endl;
      struct stat buf;
      if (!integralDumpFile.empty() && !(stat(integralDumpFile.c_str(), &buf) != -1)) {
        std::cout << "DUMPING INTEGRALS." << std::endl;
        this->dumpIntegrals(allModalPairs_);
      }
      if (dumpOnlyCoupledModes_) {
        printAllEnergies(coeffs_, allEvals);
      }
      isConverged_ = true;
      break;
    }
    if (vscfIt == (maxIter - 1)) {
      std::cout << "SCF NOT CONVERGED WITHIN MAX. ITERATIONS." << std::endl;
      break;
    }
  }
  outfile.close();
  return vscfEnergy;
}

} // namespace Scine::Colibri
