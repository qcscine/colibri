/**
 * @file VCITest.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#include <VCI/Davidson.h>
#include <VCI/VCI.h>
#include <VSCF/VSCF.h>
#include <VibrationalUtils/PESLibrary/PesCoupledOscillator.h>
#include <VibrationalUtils/PESLibrary/PesMCSP.h>
#include <VibrationalUtils/VibrationalBases/DVR.h>
#include <VibrationalUtils/VibrationalBases/DistributedGaussians.h>
#include <VibrationalUtils/VibrationalParameters.h>
#include <gmock/gmock.h>
#include <cmath>
#include <cstdlib>
#include <string>

#ifdef MPI_PARALLEL
  #include <mpi.h>
#endif

#include "writeInpFiles.h"

namespace Scine::Colibri {

using namespace testing;
using namespace VCITestInputs;

class VCITest : public Test {
 protected:
  std::string inputFileEthene = "tmp_ethene.inp";
  std::string inputFcidump = "tmp_FCIDUMPh_a2_p6";
  std::string inputFileCoupledOscillator = "tmp_coupledOsc.inp";

  class Evaluator {
    Eigen::MatrixXd A_;

   public:
    Evaluator(Eigen::MatrixXd A) : A_(std::move(A)) {
    }
    Eigen::MatrixXd evaluate(const Eigen::MatrixXd& P) {
      return A_ * P;
    }
  };

  int size = 1, id = 0;

  void finalizeMPI() {
#ifdef MPI_PARALLEL
    MPI_Finalize();
#endif
  }

  void initializeMPI() {
#ifdef MPI_PARALLEL
    MPI_Init(NULL, NULL);
#endif
  }

  void writeInputFiles(const std::string& molecule) {
#ifdef MPI_PARALLEL
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (id == 0) {
      if (molecule == "ethene") {
        writeInputEthene(inputFileEthene);
        writeFcidumpEthene(inputFcidump);
      } else if (molecule == "coupledOscillator") {
        writeInputCoupledOsc(inputFileCoupledOscillator);
      }
    }
#ifdef MPI_PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  void removeInputFiles(const std::string& molecule) {
#ifdef MPI_PARALLEL
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (id == 0) {
      if (molecule == "ethene") {
        auto successEtheneInp = std::remove(inputFileEthene.c_str());
        auto successFcidumpInp = std::remove(inputFcidump.c_str());
        auto successFcidumpOut = std::remove("tmp_FCIDUMPv_a2_p6");
        if ((successEtheneInp != 0) || (successFcidumpInp != 0) || (successFcidumpOut != 0)) {
          throw std::runtime_error("Error while removing files");
        }
      } else if (molecule == "coupledOscillator") {
        auto successInp = std::remove(inputFileCoupledOscillator.c_str());
        auto successFcidumpOut = std::remove("tmp_FCIDUMPv_coupOsc");
        if ((successInp != 0) || (successFcidumpOut != 0)) {
          throw std::runtime_error("Error while removing files");
        }
      }
    }
#ifdef MPI_PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
};

TEST_F(VCITest, Davidson) {
  initializeMPI();
  int dim = 50;
  int maxIt = 200;
  double sparsity = 0.001;
  Eigen::MatrixXd A;
  A = Eigen::MatrixXd::Random(dim, dim);

  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < i; ++j) {
      A(i, j) *= sparsity;
      A(j, i) = A(i, j);
    }
  }

  int numberOfEigenpairs = 2;
  int maxSubspaceDim = 10;

  Davidson davidson(numberOfEigenpairs, maxSubspaceDim, dim);

  davidson.setMaxIterations(maxIt);
  davidson.setDiagonalPreconditioner(A.diagonal().asDiagonal());

  Evaluator evaluator(A);
  davidson.compute<Evaluator>(evaluator, true);

  Eigen::MatrixXd eigvecs = davidson.getEigvecs();
  Eigen::VectorXd eigvals = davidson.getEigvals();

  if (id == 0) {
    std::cout << "Davidson eigenvalues:\n" << eigvals.block(0, 0, numberOfEigenpairs, 1) << std::endl;
  }

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig;
  eig.compute(A);

  Eigen::MatrixXd refVec = eig.eigenvectors();
  Eigen::VectorXd refVal = eig.eigenvalues();

  if (id == 0) {
    std::cout << "Reference eigenvalues:\n" << refVal.block(0, 0, numberOfEigenpairs, 1) << std::endl;
  }

  if (id == 0) {
    EXPECT_THAT((eigvals.block(0, 0, numberOfEigenpairs, 1) - refVal.block(0, 0, numberOfEigenpairs, 1)).norm(),
                DoubleNear(0.0, 1e-10));
  }
}

TEST_F(VCITest, DavidsonLarge) {
  int dim = 100;
  double sparsity = 0.001;
  Eigen::MatrixXd A;

  A = Eigen::MatrixXd::Random(dim, dim);

  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < i; ++j) {
      A(i, j) *= sparsity;
      A(j, i) = A(i, j);
    }
  }

  int numberOfEigenpairs = 5;
  int maxSubspaceDim = 10;

  Davidson davidson(numberOfEigenpairs, maxSubspaceDim, dim);

  davidson.setMaxIterations(1000);
  davidson.setDiagonalPreconditioner(A.diagonal().asDiagonal());

  Evaluator evaluator(A);
  davidson.compute<Evaluator>(evaluator, false);

  Eigen::MatrixXd eigvecs = davidson.getEigvecs();
  Eigen::VectorXd eigvals = davidson.getEigvals();

  if (id == 0) {
    std::cout << "Davidson eigenvalues:\n" << eigvals.block(0, 0, numberOfEigenpairs, 1) << std::endl;
  }

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig;
  eig.compute(A);

  Eigen::MatrixXd refVec = eig.eigenvectors();
  Eigen::VectorXd refVal = eig.eigenvalues();

  if (id == 0) {
    std::cout << "Reference eigenvalues:\n" << refVal.block(0, 0, numberOfEigenpairs, 1) << std::endl;
  }

  if (id == 0) {
    EXPECT_THAT((eigvals.block(0, 0, numberOfEigenpairs, 1) - refVal.block(0, 0, numberOfEigenpairs, 1)).norm(),
                DoubleNear(0.0, 1e-10));
  }
}

TEST_F(VCITest, VCI_Ethene_SelfAdjoint) {
  writeInputFiles("ethene");
  VibrationalParameters parms;
  parms.setParametersAlsoVCI(inputFileEthene);
  if (id == 0) {
    parms.printParametersVSCF();
  }
  parms.setPesMCSP(inputFcidump);
  std::shared_ptr<PES> pes = std::make_shared<PesMCSP>();

  pes->setReferenceGeometry(parms.oneModeTerms_, parms.twoModeTerms_, parms.threeModeTerms_, parms.fourModeTerms_,
                            parms.fiveModeTerms_, parms.sixModeTerms_);

  std::shared_ptr<ModalBasis> basis;
  basis = std::make_shared<DVR>(parms, pes);

  std::shared_ptr<VibrationalSCF> vscf = std::make_shared<VSCF>(parms, basis);

  if (id == 0) {
    double vscfEnergy = vscf->SCF();
    parms.printParametersVCI();

    std::shared_ptr<VibrationalCI> vci = std::make_shared<VCI>(parms, vscf);
    double vciEnergy = vci->CI() * 219474.63;

    EXPECT_THAT(vciEnergy, DoubleNear(10985.61864, 1e-05));
  }
  removeInputFiles("ethene");
}

TEST_F(VCITest, VCI_Ethene_Davidson) {
  writeInputFiles("ethene");
  VibrationalParameters parms;
  parms.setParametersAlsoVCI(inputFileEthene);
  parms.vciEVSolver_ = "Davidson";
  if (id == 0) {
    parms.printParametersVSCF();
  }
  parms.setPesMCSP(inputFcidump);
  std::shared_ptr<PES> pes = std::make_shared<PesMCSP>();

  pes->setReferenceGeometry(parms.oneModeTerms_, parms.twoModeTerms_, parms.threeModeTerms_, parms.fourModeTerms_,
                            parms.fiveModeTerms_, parms.sixModeTerms_);

  std::shared_ptr<ModalBasis> basis;
  basis = std::make_shared<DVR>(parms, pes);

  std::shared_ptr<VibrationalSCF> vscf = std::make_shared<VSCF>(parms, basis);

  if (id == 0) {
    double vscfEnergy = vscf->SCF();

    parms.printParametersVCI();

    std::shared_ptr<VibrationalCI> vci = std::make_shared<VCI>(parms, vscf);
    double vciEnergy = vci->CI() * 219474.63;

    EXPECT_THAT(vciEnergy, DoubleNear(10985.61864, 1e-05));
  }
  removeInputFiles("ethene");
}

/**
 * @brief Test of VSCF-VCI on the linearly coupled Harmonic oscillator model
 * Hamiltonian with coupling parameter = 0.1
 * The VCI result should be close to the analytical paper reference
 * Leclerc and Carrington, JCP 2014
 */

TEST_F(VCITest, GaussiansOscillatorCoupled) {
  writeInputFiles("coupledOscillator");
  VibrationalParameters parms;
  parms.setParametersAlsoVCI(inputFileCoupledOscillator);
  parms.vciEVSolver_ = "Davidson";
  if (id == 0) {
    parms.printParametersVSCF();
  }

  // Coupled version with 0.1
  std::shared_ptr<PES> pes = std::make_shared<PesCoupledOscillator>(parms, 0.1);
  std::shared_ptr<ModalBasis> basis = std::make_shared<DistributedGaussians>(parms, pes);
  std::shared_ptr<VibrationalSCF> vscf = std::make_shared<VSCF>(parms, basis);
  double vscfEnergy = NAN;
  if (id == 0) {
    double vscfEnergy = vscf->SCF();
    EXPECT_THAT(vscfEnergy, DoubleNear(3.8296274263021495, 1.0E-05));
    parms.printParametersVCI();

    std::shared_ptr<VibrationalCI> vci = std::make_shared<VCI>(parms, vscf);
    double vciEnergy = vci->CI();

    EXPECT_THAT(vciEnergy, DoubleNear(3.8164041, 1e-05));
  }
  removeInputFiles("coupledOscillator");
  finalizeMPI();
}

} // namespace Scine::Colibri
