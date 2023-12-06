/**
 * @file Davidson.inl
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef DAVIDSON_VCI_INL
#define DAVIDSON_VCI_INL

/**
 * Main routine. This is needed here to prevent linking errors from the template
 */
template<class SigmaVectorEvaluator>
void compute(SigmaVectorEvaluator& evaluator, bool verbose) {

  if (verbose) {
    std::cout << "----------------------------\n";
    std::cout << " Davidson\n";
    std::cout << "----------------------------\n";
    std::cout << " Iteration       Error      \n";
    std::cout << "----------------------------";
  }

  iterations_ = 0;
  errorFinal_ = 0;

  if (!guessProvided_) {
    projector_ = Eigen::MatrixXd::Identity(fullDim_, numEigpairs_);
    subspaceDim_ = numEigpairs_;
  } else {
    subspaceDim_ = projector_.cols();
    projector_ = modifiedGramSchmidt();
  }

  Eigen::SelfAdjointEigenSolver<Matrix> subspaceDiagonalizer;

  Eigen::DiagonalMatrix<double, -1, -1> Id;
  Id.setIdentity(fullDim_);

  while (true) {
    residuals_.resize(fullDim_, numEigpairs_);
    ritzVectors_.resize(fullDim_, numEigpairs_);

    sigmaVectors_ = evaluator.evaluate(projector_);
    interactionMatrix_ = projector_.transpose() * sigmaVectors_;

    // Diagonalization of subspace
    subspaceDiagonalizer.compute(interactionMatrix_);
    eigvals_ = subspaceDiagonalizer.eigenvalues().block(0, 0, numEigpairs_, 1);
    subspaceEigvecs_ = subspaceDiagonalizer.eigenvectors().block(0, 0, subspaceDim_, numEigpairs_);

    for (int i = 0; i < numEigpairs_; ++i) {
      ritzVectors_.col(i) = projector_ * subspaceEigvecs_.col(i);
      residuals_.col(i) = eigvals_(i) * ritzVectors_.col(i) - sigmaVectors_ * subspaceEigvecs_.col(i);
    }

    iterations_++;
    error_ = residuals_.norm();
    
    if (verbose) {
      std::cout << std::endl << std::right << std::setw(10) << iterations_ << std::string(3, ' ') << std::setw(15)
                << std::scientific << std::setprecision(5) << error_;
    }

    // Check if convergence has been reached, either all at once or each EVec individually
    if (residuals_.norm()<thresh_) { // Converge all Evecs at once
      for (int i = 0; i < numEigpairs_; ++i) {
        eigvalsFinal_(i) = eigvals_(i);
        eigvecsFinal_.col(i) = ritzVectors_.col(i);
      }
      errorFinal_ = residuals_.norm();
      if (verbose) { std::cout << std::endl << "----------------------------\n";}
      std::cout << "All eigenvectors have converged." << std::endl;
      std::cout << "Total iterations:  " << iterations_ << std::endl;
      break;
    }

    if (iterations_ == maxIterations_) {
      if (verbose) { std::cout << std::endl << "----------------------------\n";}
      std::cout << "WARNING! Davidson diagonalization did not converge in " << maxIterations_ << " iterations." << std::endl;
      errorFinal_ = error_;
      break;
    }

    /* Insert preconditioner here */
    //
    newDirection_.resize(fullDim_, numEigpairs_);
    for (int i = 0; i < numEigpairs_; ++i) {
      if (useDiagonalPreconditioner_ && (eigvals_(i)-diagonal_.diagonal()(i)!=0.0)) {
        preConditioner_ = (eigvals_(i) * Id.diagonal() - diagonal_.diagonal()).cwiseInverse().asDiagonal();
        newDirection_.col(i) = preConditioner_ * residuals_.col(i);
      } else {
        newDirection_.col(i) = residuals_.col(i);
      }
    }

    initiateModifiedGramSchmidt();
  }

  std::cout << "Final total error: " << errorFinal_ << std::endl;
}

/**
 * Main routine, sequential version
 */
template<class SigmaVectorEvaluator>
void computeSequential(SigmaVectorEvaluator& evaluator, bool verbose) {

  if (verbose) {
    std::cout << "----------------------------\n";
    std::cout << " Davidson (Sequential Conv.)\n";
    std::cout << "----------------------------\n";
    std::cout << " Iteration       Error      \n";
    std::cout << "----------------------------";
  }

  iterations_ = 0;
  errorFinal_ = 0;
  numUnconvergedEVs_ = numEigpairs_;

  if (!guessProvided_) {
    projector_ = Eigen::MatrixXd::Identity(fullDim_, numEigpairs_);
    subspaceDim_ = numEigpairs_;
  } else {
    subspaceDim_ = projector_.cols();
  }

  Eigen::SelfAdjointEigenSolver<Matrix> subspaceDiagonalizer;

  Eigen::DiagonalMatrix<double, -1, -1> Id;
  Id.setIdentity(fullDim_);

  while (true) {
    residuals_.resize(fullDim_, numUnconvergedEVs_);
    ritzVectors_.resize(fullDim_, numUnconvergedEVs_);

    sigmaVectors_ = evaluator.evaluate(projector_);
    interactionMatrix_ = projector_.transpose() * sigmaVectors_;

    // Diagonalization of subspace
    subspaceDiagonalizer.compute(interactionMatrix_);
    eigvals_ = subspaceDiagonalizer.eigenvalues().block(0, 0, numUnconvergedEVs_, 1);
    subspaceEigvecs_ = subspaceDiagonalizer.eigenvectors().block(0, 0, subspaceDim_, numUnconvergedEVs_);

    for (int i = 0; i < numUnconvergedEVs_; ++i) {
      ritzVectors_.col(i) = projector_ * subspaceEigvecs_.col(i);
      residuals_.col(i) = eigvals_(i) * ritzVectors_.col(i) - sigmaVectors_ * subspaceEigvecs_.col(i);
    }

    iterations_++;

    error_ = residuals_.norm() + errorFinal_;
    
    if (verbose) {
      std::cout << std::endl << std::right << std::setw(10) << iterations_ << std::string(3, ' ') << std::setw(15)
                << std::scientific << std::setprecision(5) << error_;
    }

    // Check if convergence has been reached
    foundNumEVsInIt_ = 0;
    for (int i = 0; i < residuals_.cols(); ++i) {
      if (residuals_.col(i).norm()<thresh_) { // Converge individual EVecs
        if (verbose) { std::cout << "     -->  Ev #" << numEigpairs_-numUnconvergedEVs_ << " is " << eigvals_(i);}
        eigvalsFinal_(numEigpairs_-numUnconvergedEVs_) = eigvals_(i);
        eigvecsFinal_.col(numEigpairs_-numUnconvergedEVs_) = ritzVectors_.col(i);
        errorFinal_+=residuals_.col(i).norm();
        ++foundNumEVsInIt_;
        --numUnconvergedEVs_;
      } else {
        break;
      }
    }
    
    
    if (numUnconvergedEVs_==0) {
      if (verbose) { std::cout << std::endl << "----------------------------\n";}
      std::cout << "All eigenvectors have converged." << std::endl;
      std::cout << "Total iterations:  " << iterations_ << std::endl;
      break;
    }

    if (iterations_ == maxIterations_) {
      if (verbose) { std::cout << std::endl << "----------------------------\n";}
      std::cout << "WARNING! Davidson diagonalization did not converge in " << maxIterations_ << " iterations." << std::endl;
      errorFinal_ = error_;
      break;
    }

    /* Insert preconditioner here */
    //
    newDirection_.resize(fullDim_, numUnconvergedEVs_);
    for (int i = 0; i < numUnconvergedEVs_; ++i) {
      if (useDiagonalPreconditioner_ && (eigvals_(i+foundNumEVsInIt_)-diagonal_.diagonal()(i+foundNumEVsInIt_)!=0.0)) {
        preConditioner_ = (eigvals_(i+foundNumEVsInIt_) * Id.diagonal() - diagonal_.diagonal()).cwiseInverse().asDiagonal();
        newDirection_.col(i) = preConditioner_ * residuals_.col(i+foundNumEVsInIt_);
      } else {
        newDirection_.col(i) = residuals_.col(i+foundNumEVsInIt_);
      }
    }

    initiateModifiedGramSchmidtAgainstEVecs();
  }

  std::cout << "Final total error: " << errorFinal_ << std::endl;
}

#endif