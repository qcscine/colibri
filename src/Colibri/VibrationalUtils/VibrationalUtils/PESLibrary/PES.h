/**
 * @file PES.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#ifndef PES_H
#define PES_H

#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <VibrationalUtils/Tensor.h>
#include <VibrationalUtils/VibrationalParameters.h>
#include <Eigen/Dense>
#include <iomanip>
#include <iostream>

namespace Scine {
namespace Colibri {

/**
 * @brief PES Base class all PES objects should derive from.
 *
 * The PES objects must implement the getPES method (for various override) that
 * calculate the energy for a given displacement.
 */

class PES {
 public:
  /**
   * @brief Default constructor.
   * The referenceGeometry_ is initialized with 0 so that we can check
   * afterwards if it has been initialized explicitly or not.
   */
  PES() : referenceGeometry_(0){};

  std::shared_ptr<PES> clone() const {
    return std::shared_ptr<PES>(this->cloneImpl());
  };

  /**
   * @brief Default destructor
   * This is required to avoid memory leaks when creating pointer-to-base class
   * from the constructor of the derived classes.
   */
  virtual ~PES() = default;

  /**
   * @brief This virtual methods should implement the getter for the energy at
   * the referenceGeometry_.
   * @return double
   */
  virtual double getPES() const = 0;

  /**
   * @brief returns PES at displacement of one coordinate
   * @param qi displacement of i-th normal coordinate
   * @param modei i-th mode
   * @return double
   */
  virtual double getPES(double qi, int modei) const = 0;

  /**
   * @brief returns PES at displacement of two coordinates
   * @param qi displacement of i-th normal coordinate
   * @param qj displacement of j-th normal coordinate
   * @param modei i-th mode
   * @param modej j-th mode
   * @return double
   */
  virtual double getPES(double qi, double qj, int modei, int modej) const = 0;

  /**
   * @brief returns PES at displacement of three coordinates
   * @param qi displacement of i-th normal coordinate
   * @param qj displacement of j-th normal coordinate
   * @param qk displacement of k-th normal coordinate
   * @param modei i-th mode
   * @param modej j-th mode
   * @param modek k-th mode
   * @return double
   */
  virtual double getPES(double qi, double qj, double qk, int modei, int modej, int modek) const {};

  /**
   * @brief Setter for the force constants.
   * @param k1 First-order constants.
   * @param k2 Second-order constants.
   * @param k3 Third-order constants.
   * @param k4 Fourth-order constants.
   */
  virtual void setReferenceGeometry(tensor<1>& k1, tensor<2>& k2, tensor<3>& k3, tensor<4>& k4){};

  /**
   * @brief Setter for the force constants in the MC-storage.
   * @param k1 one-mode coeffs.
   * @param k2 two-mode coeffs.
   * @param k3 three-mode coeffs.
   * @param k4 four-mode coeffs.
   * @param k5 five-mode coeffs.
   * @param k6 six-mode coeffs.
   */
  virtual void setReferenceGeometry(tensorMC<1>& k1, tensorMC<2>& k2, tensorMC<3>& k3, tensorMC<4>& k4, tensorMC<5>& k5,
                                    tensorMC<6>& k6, int taylorExpOrder){};

  virtual void setReferenceGeometry(tensorMCSP<1>& k1, tensorMCSP<2>& k2, tensorMCSP<3>& k3, tensorMCSP<4>& k4,
                                    tensorMCSP<5>& k5, tensorMCSP<6>& k6){};

  /**
   * @brief returns solely the one-body contribution to the potential
   * @param q displacement of the normal coordinate
   * @param mode the normal mode
   * @return double one-body PES contribution
   */
  double getOneBodyPES(const double q, int mode) const {
    try {
      return getPES(q, mode);
    }
    catch (std::exception& e) {
      std::cout << "Exception thrown in getOneBodyPes!" << std::endl;
      std::cout << e.what() << std::endl;
      return 0.0;
    }
  };

  /**
   * @brief returns solely the two-body contribution to the potential
   *
   * @param qi displacement of i-th normal coordinate
   * @param qj displacement of j-th normal coordinate
   * @param modei i-th mode
   * @param modej j-th mode
   * @return double two-body PES contribution
   */
  double getTwoBodyPES(double qi, double qj, int modei, int modej) const {
    try {
      return getPES(qi, qj, modei, modej);
    }
    catch (std::exception& e) {
      std::cout << "Exception thrown in getTwoBodyPes!" << std::endl;
      std::cout << e.what() << std::endl;
      return 0.0;
    }
  };

  /**
   * @brief returns solely the three-body contribution to the potential
   *
   * @param qi displacement of i-th normal coordinate
   * @param qj displacement of j-th normal coordinate
   * @param qk displacement of k-th normal coordinate
   * @param modei i-th mode
   * @param modej j-th mode
   * @param modek k-th mode
   * @return double three-body PES contribution
   */
  double getThreeBodyPES(double qi, double qj, double qk, int modei, int modej, int modek) const {
    try {
      return getPES(qi, qj, qk, modei, modej, modek);
    }
    catch (std::exception& e) {
      std::cout << "Exception thrown in getThreeBodyPes!" << std::endl;
      std::cout << e.what() << std::endl;
      return 0.0;
    }
  };

  virtual double getHarmonicPES(double qi, int modei) const = 0;

  virtual int getNumSPsForVSCF() const {
    return 0;
  };

  virtual void setNumSPsForVSCF(int /*num*/){};

  /**
   * @brief Sets the reference geometry for the PES in Cartesian coordinates.
   * @param referenceGeometry input stream of the xyz file defining the
   * reference geometry.
   */
  void setReferenceGeometry(std::stringstream& referenceGeometry) {
    referenceGeometry_ = Utils::XyzStreamHandler::read(referenceGeometry);
  }

  Utils::AtomCollection getReferenceGeometry() {
    return referenceGeometry_;
  }

 protected:
  template<std::size_t rank>
  tensor<rank> sortForceConstants(const tensor<rank>& forceConst) {
    tensor<rank> tensorSorted;
    for (const auto& ki : forceConst) {
      std::array<int, rank> inds = ki.first;
      std::sort(inds.begin(), inds.end());
      tensorSorted[inds] = ki.second;
    };
    return tensorSorted;
  }

  template<std::size_t mcorder>
  tensorMC<mcorder> sortForceConstants(const tensorMC<mcorder>& forceConst) {
    tensorMC<mcorder> tensorSorted;
    for (const auto& ki : forceConst) {
      std::array<std::pair<int, int>, mcorder> inds = ki.first;
      std::sort(inds.begin(), inds.end(), [](std::pair<int, int> i, std::pair<int, int> j) { return i.first < j.first; });
      tensorSorted[inds] = ki.second;
    };
    return tensorSorted;
  }

  // for MCSP not only the mode indices, but also the exponents need to be
  // sorted!
  template<std::size_t mcorder>
  tensorMCSP<mcorder> sortForceConstants(const tensorMCSP<mcorder>& forceConst) {
    tensorMCSP<mcorder> tensorSorted;
    for (const auto& ki : forceConst) {
      std::array<int, mcorder> inds = ki.first;
      std::array<std::pair<int, int>, mcorder> inds_exps;
      std::vector<std::pair<std::array<int, mcorder>, double>> new_entry;
      for (auto it = ki.second.begin(); it != ki.second.end(); it++) {
        for (int i = 0; i < inds.size(); i++) {
          inds_exps[i] = std::make_pair(inds[i], it->first[i]);
        }
        std::sort(inds_exps.begin(), inds_exps.end(),
                  [](std::pair<int, int> i, std::pair<int, int> j) { return i.first < j.first; });
        std::array<int, mcorder> exps;
        for (int i = 0; i < inds.size(); i++) {
          exps[i] = inds_exps[i].second;
        }
        new_entry.push_back(std::make_pair(exps, it->second));
      }
      std::sort(inds.begin(),
                inds.end()); // This line needs to go after the other sorting
      tensorSorted[inds] = new_entry;
    };
    return tensorSorted;
  }

  static std::vector<std::pair<double, Eigen::VectorXd>> sortEigenPairs(const Eigen::MatrixXd& eVecs,
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

  /* Members */
  Utils::AtomCollection referenceGeometry_;

 private:
  virtual std::shared_ptr<PES> cloneImpl() const = 0;
};

} // namespace Colibri
} // namespace Scine

#endif
