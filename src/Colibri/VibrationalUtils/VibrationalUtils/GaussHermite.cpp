/**
 * @file GaussHermite.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#include "GaussHermite.h"
#include <omp.h>
#include <iostream>
#include <numeric>

namespace Scine::Colibri {

GaussHermiteQuad::GaussHermiteQuad() {
  map_qp.insert(std::pair<int, std::vector<quadpoint>>(6, gh_6));
  map_qp.insert(std::pair<int, std::vector<quadpoint>>(8, gh_8));
  map_qp.insert(std::pair<int, std::vector<quadpoint>>(12, gh_12));
  map_qp.insert(std::pair<int, std::vector<quadpoint>>(16, gh_16));
  map_qp.insert(std::pair<int, std::vector<quadpoint>>(32, gh_32));
}

GaussHermiteQuad::GaussHermiteQuad(const GaussHermiteQuad& rhs) {
  map_qp = rhs.map_qp;
}

double GaussHermiteQuad::oneDimIntegrate(std::function<double(const double)> const& f, int key) const {
  std::vector<quadpoint> contour = getContour(key);
  double res = 0.;
  for (auto& it : contour) {
    res += it.weight * f(it.node);
  }
  return res;
}

double GaussHermiteQuad::twoDimIntegrate(std::function<double(const double, const double)> const& f, int key) const {
  std::vector<quadpoint> contour = getContour(key);
  double res = 0.;
  for (auto it = contour.begin(); it != contour.end(); it++) {
    auto g = [f, contour, it](double qj) { return f(it->node, qj); };
    for (auto& it2 : contour) {
      res += it->weight * it2.weight * g(it2.node);
    }
  }
  return res;
}

double GaussHermiteQuad::threeDimIntegrate(std::function<double(const double, const double, const double)> const& f,
                                           int key) const {
  std::vector<quadpoint> contour = getContour(key);
  double res = 0.;
  for (auto it = contour.begin(); it != contour.end(); ++it) {
    auto g = [f, contour, it](const double& qj, const double& qk) { return f(it->node, qj, qk); };
    for (auto it2 = contour.begin(); it2 != contour.end(); ++it2) {
      auto h = [g, contour, it2](const double& qk) { return g(it2->node, qk); };
      for (auto& it3 : contour) {
        res += it->weight * it2->weight * it3.weight * h(it3.node);
      }
    }
  }
  return res;
}

std::vector<GaussHermiteQuad::quadpoint> GaussHermiteQuad::getContour(int key) const {
  if (key == 6 || key == 8 || key == 12 || key == 16 || key == 32) {
    return map_qp.at(key);
  }
  throw std::invalid_argument("Number of quadrature points (numqp) not implemented.");
}

} // namespace Scine::Colibri
