/**
 * @file H20HarmonicPesCalculator.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "H2OHarmonicPesCalculator.h"
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <iomanip>
#include <iostream>

namespace Scine::Colibri {

H2OHarmonicPesCalculator::H2OHarmonicPesCalculator() {
  // Sets the reference geometry
  auto h2oOpt = std::stringstream("3\n\n"
                                  "O    0.000000    0.000000    0.128335\n"
                                  "H    0.000000    0.756661   -0.476433\n"
                                  "H    0.000000   -0.756661   -0.476433\n");
  referencePositions_ = Utils::XyzStreamHandler::read(h2oOpt).getPositions();
  // Load the Hessian matrix
  Hessian_ = HessianMatrixType::Zero();
  Hessian_(0, 0) = 0.134940E-04;
  Hessian_(0, 3) = -0.674698E-05;
  Hessian_(3, 0) = -0.674698E-05;
  Hessian_(0, 6) = -0.674698E-05;
  Hessian_(6, 0) = -0.674698E-05;
  Hessian_(1, 1) = 0.641571E+00;
  Hessian_(1, 4) = -0.320786E+00;
  Hessian_(4, 1) = -0.320786E+00;
  Hessian_(1, 5) = 0.256380E+00;
  Hessian_(5, 1) = 0.256380E+00;
  Hessian_(1, 7) = -0.320786E+00;
  Hessian_(7, 1) = -0.320786E+00;
  Hessian_(1, 8) = -0.256380E+00;
  Hessian_(8, 1) = -0.256380E+00;
  Hessian_(2, 2) = 0.451973E+00;
  Hessian_(2, 4) = 0.193587E+00;
  Hessian_(4, 2) = 0.193587E+00;
  Hessian_(2, 7) = -0.193587E+00;
  Hessian_(7, 2) = -0.193587E+00;
  Hessian_(2, 5) = -0.225986E+00;
  Hessian_(5, 2) = -0.225986E+00;
  Hessian_(2, 8) = -0.225986E+00;
  Hessian_(8, 2) = -0.225986E+00;
  Hessian_(3, 3) = 0.705273E-05;
  Hessian_(4, 4) = 0.353999E+00;
  Hessian_(4, 5) = -0.224983E+00;
  Hessian_(5, 4) = -0.224983E+00;
  Hessian_(4, 7) = -0.332134E-01;
  Hessian_(7, 4) = -0.332134E-01;
  Hessian_(4, 8) = 0.313966E-01;
  Hessian_(8, 4) = 0.313966E-01;
  Hessian_(5, 5) = 0.215450E+00;
  Hessian_(5, 7) = -0.313966E-01;
  Hessian_(7, 5) = -0.313966E-01;
  Hessian_(5, 8) = 0.105359E-01;
  Hessian_(8, 5) = 0.105359E-01;
  Hessian_(6, 6) = 0.705273E-05;
  Hessian_(7, 7) = 0.353999E+00;
  Hessian_(7, 8) = 0.224983E+00;
  Hessian_(8, 7) = 0.224983E+00;
  Hessian_(8, 8) = 0.215450E+00;
}

H2OHarmonicPesCalculator::H2OHarmonicPesCalculator(const H2OHarmonicPesCalculator& rhs) : CloneInterface(rhs) {
  // Standard operations
  setStructure(*rhs.getStructure());
  results_ = rhs.results();
  setRequiredProperties(rhs.getRequiredProperties());
  *(this->settings_) = rhs.settings();
  Hessian_ = rhs.Hessian_;
  referencePositions_ = rhs.referencePositions_;
}

std::unique_ptr<Utils::AtomCollection> H2OHarmonicPesCalculator::getStructure() const {
  return std::make_unique<Utils::AtomCollection>(elements_, positions_);
}

Utils::Settings& H2OHarmonicPesCalculator::settings() {
  return *settings_;
}

const Utils::Settings& H2OHarmonicPesCalculator::settings() const {
  return *settings_;
}

Utils::Results& H2OHarmonicPesCalculator::results() {
  return results_;
}

const Utils::Results& H2OHarmonicPesCalculator::results() const {
  return results_;
}

const Utils::Results& H2OHarmonicPesCalculator::calculate(std::string /*description*/) {
  results_ = Utils::Results();
  double E = 0.;
  for (int i = 0; i < nAtoms; i++) {
    auto dri = (positions_.row(i) - referencePositions_.row(i));
    for (int j = 0; j < nAtoms; j++) {
      auto drj = (positions_.row(j) - referencePositions_.row(j));
      for (int iX = 0; iX < 3; iX++) {
        for (int jX = 0; jX < 3; jX++) {
          E += dri[iX] * drj[jX] * Hessian_(3 * i + iX, 3 * j + jX) / 2.;
        }
      }
    }
  }
  results_.set<Utils::Property::Energy>(E);
  results_.set<Utils::Property::SuccessfulCalculation>(true);
  return results_;
}

void H2OHarmonicPesCalculator::setRequiredProperties(const Utils::PropertyList& requiredProperties) {
  requiredProperties_ = requiredProperties;
}

Utils::PropertyList H2OHarmonicPesCalculator::getRequiredProperties() const {
  return requiredProperties_;
}

Utils::PropertyList H2OHarmonicPesCalculator::possibleProperties() const {
  return Utils::Property::Energy;
}

void H2OHarmonicPesCalculator::setStructure(const Utils::AtomCollection& structure) {
  elements_ = structure.getElements();
  positions_ = structure.getPositions();
}

bool H2OHarmonicPesCalculator::supportsMethodFamily(const std::string& methodFamily) const {
  return methodFamily == name();
}

std::string H2OHarmonicPesCalculator::name() const {
  return Scine::Colibri::H2OHarmonicPesCalculator::model;
}

void H2OHarmonicPesCalculator::modifyPositions(Utils::PositionCollection newPositions) {
  positions_ = newPositions;
}

const Utils::PositionCollection& H2OHarmonicPesCalculator::getPositions() const {
  return positions_;
}

} // namespace Scine::Colibri
