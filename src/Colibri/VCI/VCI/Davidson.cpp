/**
 * @file Davidson.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef DAVIDSON_VCI_CPP
#define DAVIDSON_VCI_CPP

#include "Davidson.h"
#include <algorithm>

namespace Scine::Colibri {

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

Davidson::Davidson(int numberOfEigenPairs, int maxSubspaceDimension, int fullDimension)
  : numEigpairs_(numberOfEigenPairs), maxSubspaceDim_(maxSubspaceDimension), fullDim_(fullDimension) {
  eigvalsFinal_.resize(numEigpairs_);
  eigvecsFinal_.resize(fullDim_, numEigpairs_);
}

int Davidson::getSubspaceDim() const {
  return subspaceDim_;
}

int Davidson::getNumUnconvergedEVs() const {
  return numUnconvergedEVs_;
}

int Davidson::getIterations() const {
  return iterations_;
}

const Vector& Davidson::getEigvals() const {
  return eigvalsFinal_;
}

const Matrix& Davidson::getEigvecs() const {
  return eigvecsFinal_;
}

const Matrix& Davidson::getProjector() const {
  return projector_;
}

double Davidson::getError() const {
  return error_;
}

double Davidson::getFinalError() const {
  return errorFinal_;
}

const Matrix& Davidson::getSubspaceEigvecs() const {
  return subspaceEigvecs_;
}

void Davidson::setMaxIterations(int maxIt) {
  maxIterations_ = maxIt;
}

void Davidson::setConvThresh(double threshold) {
  thresh_ = threshold;
}

void Davidson::setDiagonalPreconditioner(Eigen::DiagonalMatrix<double, -1, -1>&& preCond) {
  diagonal_ = preCond;
  useDiagonalPreconditioner_ = true;
}

void Davidson::setGuess(Eigen::MatrixXd&& guess) {
  projector_ = guess;
  guessProvided_ = true;
  if (projector_.cols() < numEigpairs_) {
    throw std::runtime_error("Provide at least as many trial vectors, as eigenpairs are requested!");
  }
}

void Davidson::initiateModifiedGramSchmidt() {
  if (subspaceDim_ + numEigpairs_ > maxSubspaceDim_) {
    // collapse subspace
    subspaceDim_ = 2 * numEigpairs_;
    projector_.resize(fullDim_, subspaceDim_);

    for (auto i = 0; i < numEigpairs_; ++i) {
      projector_.col(i) = ritzVectors_.col(i);
      projector_.col(i + numEigpairs_) = newDirection_.col(i);
    }
  } else {
    // enlarge subspace
    int oldSubspaceDim = subspaceDim_;
    subspaceDim_ = subspaceDim_ + numEigpairs_;
    projector_.conservativeResize(fullDim_, subspaceDim_);

    for (int i = 0; i < numEigpairs_; ++i) {
      projector_.col(oldSubspaceDim + i) = newDirection_.col(i);
    }
  }
  projector_ = modifiedGramSchmidt();
}

void Davidson::initiateModifiedGramSchmidtAgainstEVecs() {
  if ((foundNumEVsInIt_ > 0) || (subspaceDim_ + numUnconvergedEVs_ > maxSubspaceDim_ - (numEigpairs_ - numUnconvergedEVs_))) {
    // collapse subspace
    subspaceDim_ = 2 * numUnconvergedEVs_;
    projector_.resize(fullDim_, subspaceDim_);

    for (auto i = 0; i < numUnconvergedEVs_; ++i) {
      projector_.col(i) = ritzVectors_.col(i + foundNumEVsInIt_);
      projector_.col(i + numUnconvergedEVs_) = newDirection_.col(i);
    }

    foundNumEVsInIt_ = 0;
  } else {
    // enlarge subspace
    int oldSubspaceDim = subspaceDim_;
    subspaceDim_ = subspaceDim_ + numUnconvergedEVs_;
    projector_.conservativeResize(fullDim_, subspaceDim_);

    for (int i = 0; i < numUnconvergedEVs_; ++i) {
      projector_.col(oldSubspaceDim + i) = newDirection_.col(i);
    }
  }
  projector_ = modifiedGramSchmidtAgainstEVecs();
}

Eigen::MatrixXd Davidson::modifiedGramSchmidt() {
  // Calculate modified Gram-Schmidt QR decomposition
  Matrix Q(projector_.rows(), projector_.cols());

  for (int j = 0; j < Q.cols(); ++j) {
    Q.col(j) = projector_.col(j);
    for (int i = 0; i < j; ++i) {
      Q.col(j) -= Q.col(i) * (Q.col(i).adjoint() * Q.col(j));
    }
    Q.col(j).normalize();
  }

  return Q;
}

// calculate Gram-Schmidt also against previously found EVecs
Eigen::MatrixXd Davidson::modifiedGramSchmidtAgainstEVecs() {
  // Concatenate EVS and Projector
  Matrix Q(projector_.rows(), numEigpairs_ - numUnconvergedEVs_ + projector_.cols());
  Q << eigvecsFinal_.block(0, 0, subspaceDim_, numEigpairs_ - numUnconvergedEVs_), projector_;

  // Calculate modified Gram-Schmidt QR decomposition

  for (int j = 0; j < Q.cols(); ++j) {
    for (int i = 0; i < j; ++i) {
      Q.col(j) -= Q.col(i) * (Q.col(i).adjoint() * Q.col(j));
    }
    Q.col(j).normalize();
  }

  // return only the projector
  return Q.block(0, numEigpairs_ - numUnconvergedEVs_, projector_.rows(), projector_.cols());
}

} // namespace Scine::Colibri

#endif