/**
 * @file PesCartesian.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#include "PesCartesian.h"
#include <Utils/GeometricDerivatives/NormalModeAnalysis.h>
#include <Utils/GeometricDerivatives/NumericalHessianCalculator.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Geometry/GeometryUtilities.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/IO/FormattedIOUtils.h>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <utility>

namespace Scine::Colibri {

const double PesCartesian::stepForNumericalDifferences = 0.001;

PesCartesian::PesCartesian(const PesCartesian& rhs) : CloneInterface(rhs) {
  hessian_ = rhs.hessian_;
  normalModesContainer_ = rhs.normalModesContainer_;
  masses_ = rhs.masses_;
  referenceGeometry_ = rhs.referenceGeometry_;
  currentGeometry_ = rhs.currentGeometry_;
  pesCalculator_ = rhs.pesCalculator_->clone();
}

PesCartesian::PesCartesian(std::shared_ptr<Core::Calculator> pesCalculator)
  : pesCalculator_(std::move(std::move(pesCalculator))) {
  referenceGeometry_ = *(pesCalculator_->getStructure());
  masses_ = Utils::Geometry::Properties::getMasses(referenceGeometry_.getElements());
  calculateNormalModesAndHarmonicFrequencies();
  // printNormalModeInformation(std::cout);
}

double PesCartesian::getPES() const {
  double energy = 0;
#pragma omp critical
  { energy = getEnergy(referenceGeometry_); }
  return energy;
}

double PesCartesian::getPES(double qi, int modei) const {
  double energy = 0;
#pragma omp critical
  {
    assert(pesCalculator_);
    Utils::AtomCollection updatedGeom = this->updateReferenceGeometryFromDisplacement({modei}, {qi});
    energy = getEnergy(updatedGeom);
  }
  return energy - getPES();
}

double PesCartesian::getPES(double qi, double qj, int modei, int modej) const {
  double energy = 0;
#pragma omp critical
  {
    assert(pesCalculator_);
    Utils::AtomCollection updatedGeom = this->updateReferenceGeometryFromDisplacement({modei, modej}, {qi, qj});
    energy = getEnergy(updatedGeom);
  }
  return energy - getPES() - getPES(qi, modei) - getPES(qj, modej);
}

double PesCartesian::getHarmonicPES(double /*qi*/, int /*modei*/) const {
  throw std::runtime_error("Harmonic PES for Cartesian PES NYI.");
}

double PesCartesian::getEnergy(const Utils::AtomCollection& updatedGeom) const {
  pesCalculator_->setStructure(updatedGeom);
  Utils::Results r = pesCalculator_->calculate("");
  auto en = r.template get<Utils::Property::Energy>();
  return en;
}

void PesCartesian::calculateNormalModesAndHarmonicFrequencies() {
  pesCalculator_->setStructure(referenceGeometry_);
  Scine::Utils::NumericalHessianCalculator hessianCalc(*pesCalculator_);
  hessian_ = hessianCalc.calculateFromEnergyDifferences(stepForNumericalDifferences);
  normalModesContainer_ = Scine::Utils::NormalModeAnalysis::calculateNormalModes(
      hessian_, referenceGeometry_.getElements(), referenceGeometry_.getPositions(), false);
}

void PesCartesian::printNormalModeInformation(std::ostream& stream) const {
  Utils::matrixPrettyPrint(stream, normalModesContainer_, referenceGeometry_.getElements());
}

double PesCartesian::printHarmonicZPVE() const {
  auto wavenumbers = normalModesContainer_.getWaveNumbers();
  double harmZPVE = 0;
  std::cout << "  Normal Frequencies:\n  [Rot. and trans. freq. removed, "
               "imaginary freq. shown as negative.]\n"
            << std::endl;
  printf("  %4s %8s\n", "#", "cm^-1");
  for (unsigned int i = 0; i < wavenumbers.size(); i++) {
    printf("  %4d %+13.6f\n", i + 1, wavenumbers[i]);
    harmZPVE += 0.5 * wavenumbers[i];
  }
  std::cout << std::endl
            << "Harmonic ZPVE is: " << harmZPVE / Utils::Constants::invCentimeter_per_hartree << " Hartree" << std::endl;
  std::cout << "Harmonic ZPVE is: " << harmZPVE << " cm^-1" << std::endl;
  return harmZPVE;
}

void PesCartesian::setReferenceGeometry(tensor<1>& /*k1*/, tensor<2>& /*k2*/, tensor<3>& /*k3*/, tensor<4>& /*k4*/) {
  throw std::runtime_error("Force constants not available for the PesCartesian class");
}

Utils::AtomCollection PesCartesian::updateReferenceGeometryFromDisplacement(std::vector<int>&& idxModes,
                                                                            std::vector<double>&& steps) const {
  Utils::AtomCollection currentGeometry = referenceGeometry_;
  assert(idxModes.size() == steps.size());
  Eigen::Matrix<double, Eigen::Dynamic, 3> finalPositions = referenceGeometry_.getPositions();
  for (std::size_t iModes = 0; iModes < idxModes.size(); iModes++) {
    auto modeNoWeight = normalModesContainer_.getMode(idxModes[iModes]);
    int nAtoms = pesCalculator_->getStructure()->size();
    for (int iAtoms = 0; iAtoms < nAtoms; iAtoms++) {
      modeNoWeight(iAtoms, 0) /= std::sqrt(1822.888486209);
      modeNoWeight(iAtoms, 1) /= std::sqrt(1822.888486209);
      modeNoWeight(iAtoms, 2) /= std::sqrt(1822.888486209);
    }
    modeNoWeight *= steps[iModes];
    finalPositions += modeNoWeight;
  }
  currentGeometry.setPositions(finalPositions);
  return currentGeometry;
}

} // namespace Scine::Colibri
