/**
 * @file Davidson.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

/**
 * https://epubs.siam.org/doi/abs/10.1137/0915004
 * SIAM J. Sci. Comput., 15(1), 62–76. (15 pages) 1994
 * Code core originally by Robin, heavily modified by Nina 2022
 * See also:
 * https://dl.acm.org/doi/10.1145/2543696
 * ACM Transactions on Mathematical Software, Vol. 40, No. 2, Article 13, 2014
 */

#ifndef DAVIDSON_VCI_H
#define DAVIDSON_VCI_H
#include <Eigen/Eigenvalues>
#include <chrono>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <vector>

namespace Scine {
namespace Colibri {

class Davidson {
  /**
   * @brief class for iterative diagonalisation with generalized hermitian
   * Davidson
   * @param maxIterations maximum number of Davidson iterations
   * @param numEigpairs number of Eigenpairs which are calculated
   * @param maxSubspaceDim maximum dimension of diagonalized subspace, if it
   * grows above this the subspace is collapsed
   * @param fullDim dimension of the full problem
   * @param thresh convergence threshold, applied to each individual eigenvector
   * residual
   */

  using Matrix = Eigen::MatrixXd;
  using Vector = Eigen::VectorXd;

 private:
  int maxIterations_ = 100;
  int iterations_;
  int numEigpairs_;
  int maxSubspaceDim_;
  int fullDim_;
  int subspaceDim_;
  int numUnconvergedEVs_;
  int foundNumEVsInIt_;
  double thresh_ = 1e-10;
  double error_;
  double errorFinal_;
  bool useDiagonalPreconditioner_ = false;
  bool guessProvided_ = false;

  Vector eigvals_;
  Vector eigvalsFinal_;
  Matrix projector_;
  Matrix sigmaVectors_;
  Matrix interactionMatrix_;
  Matrix subspaceEigvecs_;
  Matrix ritzVectors_;
  Matrix eigvecsFinal_;
  Matrix residuals_;
  Matrix newDirection_;
  Eigen::DiagonalMatrix<double, -1, -1> preConditioner_;
  Eigen::DiagonalMatrix<double, -1, -1> diagonal_;

  void initiateModifiedGramSchmidt();
  void initiateModifiedGramSchmidtAgainstEVecs();
  Eigen::MatrixXd modifiedGramSchmidt();
  Eigen::MatrixXd modifiedGramSchmidtAgainstEVecs();

 public:
  Davidson(int numberOfEigenPairs, int maxSubspaceDimension, int fullDimension);

  int getSubspaceDim() const;
  int getNumUnconvergedEVs() const;
  int getIterations() const;

  const Vector& getEigvals() const;
  const Matrix& getEigvecs() const;
  const Matrix& getProjector() const;
  const Matrix& getSubspaceEigvecs() const;

  double getError() const;
  double getFinalError() const;

  void setMaxIterations(int maxIt);
  void setConvThresh(double threshold);
  void setDiagonalPreconditioner(Eigen::DiagonalMatrix<double, -1, -1>&& preCond);
  void setGuess(Eigen::MatrixXd&& guess);

// The main routine had to be put in this file! Otherwise linker issues arise
// due to templating
#include "Davidson.inl"
  // template<class SigmaVectorEvaluator> void compute(SigmaVectorEvaluator&
  // evaluator, bool verbose = true); template<class SigmaVectorEvaluator> void
  // computeSequential(SigmaVectorEvaluator& evaluator, bool verbose = true);
};

} // namespace Colibri
} // namespace Scine

#endif
