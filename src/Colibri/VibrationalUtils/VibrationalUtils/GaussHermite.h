/**
 * @file GaussHermite.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#ifndef GAUSS_HERM_H
#define GAUSS_HERM_H

#include <omp.h>
#include <functional>
#include <map>
#include <utility>
#include <vector>

namespace Scine {
namespace Colibri {

class GaussHermiteQuad {
  /**
   * @brief Implementation of Gauss-Hermite quadrature rule.
   * @param map_qp Maps the key value of number of nodes (8,12,16 or 32) to the
   * corresponding quadrature node-weight pairs.
   * @param quadpoint is a struct for holding a node-weight pair
   * @param gh_8 vector of node-weight pairs for 8-point Gauss-Hermite
   * quadrature rule
   * @param gh_12 vector of node-weight pairs for 12-point Gauss-Hermite
   * quadrature rule
   * @param gh_16 vector of node-weight pairs for 16-point Gauss-Hermite
   * quadrature rule
   * @param gh_32 vector of node-weight pairs for 32-point Gauss-Hermite
   * quadrature rule
   *
   */
 public:
  /**
   * @brief Constructs all quadrature rules and the corresponding maps
   *
   */
  GaussHermiteQuad();
  GaussHermiteQuad(const GaussHermiteQuad& rhs);
  ~GaussHermiteQuad() = default;

  /**
   * @brief Evaluates an integral over (-infinity, +infinity) of e^(-x^2)f(x)
   *
   * @param lambda_f lambda function of form f(x)
   * @param key number of nodes
   * @return double
   */
  double oneDimIntegrate(std::function<double(const double)> const& lambda_f, int key) const;

  /**
   * @brief Evaluates integral over (-infinity, +infinity) of
   * e^(-x^2)e^(-y^2)f(x,y)
   *
   * @param lambda_g lambda function of form f(x, y)
   * @param key number of nodes
   * @return double
   */
  double twoDimIntegrate(std::function<double(const double, const double)> const& f, int key) const;
  double threeDimIntegrate(std::function<double(const double, const double, const double)> const& f, int key) const;

 protected:
  struct quadpoint {
    long double node;
    long double weight;
    quadpoint(const long double quad_node, const long double quad_weight) : node(quad_node), weight(quad_weight){};
  };

 private:
  std::map<int, std::vector<quadpoint>> map_qp;

  std::vector<quadpoint> getContour(int key) const;

  std::vector<quadpoint> gh_6 = {quadpoint(-2.35060497367449222281e+00, 4.53000990550884564102e-03),
                                 quadpoint(-4.36077411927616508688e-01, 7.24629595224392524086e-01),
                                 quadpoint(-1.33584907401369694976e+00, 1.57067320322856643914e-01),
                                 quadpoint(2.35060497367449222281e+00, 4.53000990550884564102e-03),
                                 quadpoint(4.36077411927616508688e-01, 7.24629595224392524086e-01),
                                 quadpoint(1.33584907401369694976e+00, 1.57067320322856643914e-01)};

  std::vector<quadpoint> gh_8 = {quadpoint(3.81186990207322116844e-01, 6.61147012558241291042e-01),
                                 quadpoint(1.15719371244678019474e+00, 2.07802325814891879546e-01),
                                 quadpoint(1.98165675669584292584e+00, 1.70779830074134754563e-02),
                                 quadpoint(2.93063742025724401920e+00, 1.99604072211367619211e-04),
                                 quadpoint(-3.81186990207322116844e-01, 6.61147012558241291042e-01),
                                 quadpoint(-1.15719371244678019474e+00, 2.07802325814891879546e-01),
                                 quadpoint(-1.98165675669584292584e+00, 1.70779830074134754563e-02),
                                 quadpoint(-2.93063742025724401920e+00, 1.99604072211367619211e-04)};

  std::vector<quadpoint> gh_12 = {quadpoint(3.14240376254359111269e-01, 5.70135236262479578326e-01),
                                  quadpoint(9.47788391240163743685e-01, 2.60492310264161129222e-01),
                                  quadpoint(1.59768263515260479672e+00, 5.16079856158839299912e-02),
                                  quadpoint(2.27950708050105990015e+00, 3.90539058462906185994e-03),
                                  quadpoint(3.02063702512088977175e+00, 8.57368704358785865472e-05),
                                  quadpoint(3.88972489786978191926e+00, 2.65855168435630160614e-07),
                                  quadpoint(-3.14240376254359111269e-01, 5.70135236262479578326e-01),
                                  quadpoint(-9.47788391240163743685e-01, 2.60492310264161129222e-01),
                                  quadpoint(-1.59768263515260479672e+00, 5.16079856158839299912e-02),
                                  quadpoint(-2.27950708050105990015e+00, 3.90539058462906185994e-03),
                                  quadpoint(-3.02063702512088977175e+00, 8.57368704358785865472e-05),
                                  quadpoint(-3.88972489786978191926e+00, 2.65855168435630160614e-07)};

  std::vector<quadpoint> gh_16 = {quadpoint(2.73481046138152452172e-01, 5.07929479016613741923e-01),
                                  quadpoint(8.22951449144655892596e-01, 2.80647458528533675357e-01),
                                  quadpoint(1.38025853919888079639e+00, 8.38100413989858294132e-02),
                                  quadpoint(1.95178799091625397740e+00, 1.28803115355099736832e-02),
                                  quadpoint(2.54620215784748136221e+00, 9.32284008624180529895e-04),
                                  quadpoint(3.17699916197995602682e+00, 2.71186009253788151199e-05),
                                  quadpoint(3.86944790486012269869e+00, 2.32098084486521065344e-07),
                                  quadpoint(4.68873893930581836465e+00, 2.65480747401118224476e-10),
                                  quadpoint(-2.73481046138152452172e-01, 5.07929479016613741923e-01),
                                  quadpoint(-8.22951449144655892596e-01, 2.80647458528533675357e-01),
                                  quadpoint(-1.38025853919888079639e+00, 8.38100413989858294132e-02),
                                  quadpoint(-1.95178799091625397740e+00, 1.28803115355099736832e-02),
                                  quadpoint(-2.54620215784748136221e+00, 9.32284008624180529895e-04),
                                  quadpoint(-3.17699916197995602682e+00, 2.71186009253788151199e-05),
                                  quadpoint(-3.86944790486012269869e+00, 2.32098084486521065344e-07),
                                  quadpoint(-4.68873893930581836465e+00, 2.65480747401118224476e-10)};

  std::vector<quadpoint> gh_32 = {quadpoint(1.94840741569399326713e-01, 3.75238352592802392864e-01),
                                  quadpoint(5.84978765435932448449e-01, 2.77458142302529898131e-01),
                                  quadpoint(9.76500463589682838499e-01, 1.51269734076642482578e-01),
                                  quadpoint(1.37037641095287183817e+00, 6.04581309559126141860e-02),
                                  quadpoint(1.76765410946320160465e+00, 1.75534288315734303030e-02),
                                  quadpoint(2.16949918360611217335e+00, 3.65489032665442807915e-03),
                                  quadpoint(2.57724953773231745414e+00, 5.36268365527972045989e-04),
                                  quadpoint(2.99249082500237420621e+00, 5.41658406181998255789e-05),
                                  quadpoint(3.41716749281857073593e+00, 3.65058512956237605727e-06),
                                  quadpoint(3.85375548547144464390e+00, 1.57416779254559402923e-07),
                                  quadpoint(4.30554795335119844506e+00, 4.09883216477089661816e-09),
                                  quadpoint(4.77716450350259639289e+00, 5.93329146339663861478e-11),
                                  quadpoint(5.27555098651588012760e+00, 4.21501021132644757306e-13),
                                  quadpoint(5.81222594951591383294e+00, 1.19734401709284866582e-15),
                                  quadpoint(6.40949814926966041214e+00, 9.23173653651829223381e-19),
                                  quadpoint(7.12581390983072757292e+00, 7.31067642738416239302e-23),
                                  quadpoint(-1.94840741569399326713e-01, 3.75238352592802392864e-01),
                                  quadpoint(-5.84978765435932448449e-01, 2.77458142302529898131e-01),
                                  quadpoint(-9.76500463589682838499e-01, 1.51269734076642482578e-01),
                                  quadpoint(-1.37037641095287183817e+00, 6.04581309559126141860e-02),
                                  quadpoint(-1.76765410946320160465e+00, 1.75534288315734303030e-02),
                                  quadpoint(-2.16949918360611217335e+00, 3.65489032665442807915e-03),
                                  quadpoint(-2.57724953773231745414e+00, 5.36268365527972045989e-04),
                                  quadpoint(-2.99249082500237420621e+00, 5.41658406181998255789e-05),
                                  quadpoint(-3.41716749281857073593e+00, 3.65058512956237605727e-06),
                                  quadpoint(-3.85375548547144464390e+00, 1.57416779254559402923e-07),
                                  quadpoint(-4.30554795335119844506e+00, 4.09883216477089661816e-09),
                                  quadpoint(-4.77716450350259639289e+00, 5.93329146339663861478e-11),
                                  quadpoint(-5.27555098651588012760e+00, 4.21501021132644757306e-13),
                                  quadpoint(-5.81222594951591383294e+00, 1.19734401709284866582e-15),
                                  quadpoint(-6.40949814926966041214e+00, 9.23173653651829223381e-19),
                                  quadpoint(-7.12581390983072757292e+00, 7.31067642738416239302e-23)};
};

} // namespace Colibri
} // namespace Scine
#endif
