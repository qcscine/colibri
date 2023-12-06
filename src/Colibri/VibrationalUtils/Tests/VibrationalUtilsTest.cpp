/**
 * @file VibrationalUtilsTest.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <VibrationalUtils/GaussHermite.h>
#include <VibrationalUtils/PESLibrary/H2OHarmonicPesCalculator.h>
#include <VibrationalUtils/PESLibrary/NMAPesCalculator.h>
#include <VibrationalUtils/PESLibrary/PesCartesian.h>
#include <VibrationalUtils/VibrationalBases/DVR.h>
#include <VibrationalUtils/VibrationalBases/DistributedGaussians.h>
#include <VibrationalUtils/VibrationalBases/ModeGrid.h>
#include <VibrationalUtils/VibrationalParameters.h>
#include <gmock/gmock.h>
#include <omp.h>
#include <Eigen/Eigenvalues>

#ifdef MPI_PARALLEL
  #include <mpi.h>
#endif

namespace Scine::Colibri {

using namespace testing;

class VibrationalUtilsTest : public Test {
 protected:
  // Class members
  ModeGrid modeGrid100CmM1With10NMax, modeGrid100CmM1With20NMax, modeGridLastModeNma;
  std::vector<ModeGrid> modeGridNMA, modeGridH2O;
  std::shared_ptr<NMAPesCalculator> nmaPesCalculator_;
  std::shared_ptr<H2OHarmonicPesCalculator> h2oPesCalculator_;
  std::shared_ptr<PES> pesCartesianNma_, pesCartesianH2O_;
  std::stringstream nmaOpt, h2oOpt;

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

  void getRank() {
#ifdef MPI_PARALLEL
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
#else
    id = 0;
#endif
  }

  int id = 0;

 private:
  // Setup function
  void SetUp() final {
    nmaOpt = std::stringstream("12\n\n"
                               "H    2.599018    0.636239   -0.691697\n"
                               "H    1.886641   -0.948031   -0.249315\n"
                               "H    2.495868    0.149889    1.029947\n"
                               "C    1.984051    0.102848    0.053387\n"
                               "N    0.643357    0.650678    0.127818\n"
                               "H    0.525182    1.619283    0.394361\n"
                               "O   -0.422052   -1.272056   -0.479170\n"
                               "C   -0.472330   -0.094752   -0.149604\n"
                               "C   -1.794371    0.647203   -0.019215\n"
                               "H   -2.320364    0.600899   -0.985014\n"
                               "H   -2.422247    0.120503    0.715590\n"
                               "H   -1.691736    1.700268    0.284410\n");
    h2oOpt = std::stringstream("3\n\n"
                               "O    0.000000    0.000000    0.120954\n"
                               "H    0.000000    0.756661   -0.483815\n"
                               "H    0.000000   -0.756661   -0.483815\n");
    modeGrid100CmM1With10NMax = ModeGrid(100., 10, 100);
    modeGrid100CmM1With20NMax = ModeGrid(100., 20, 100);
    modeGridLastModeNma = ModeGrid(0.016539497070800, 20, 100);
    int maxQuantum = 5;
    int nGridPoints = 30;
    modeGridNMA = {
        ModeGrid(0.000236929434623, maxQuantum, nGridPoints), ModeGrid(0.000373619493059, maxQuantum, nGridPoints),
        ModeGrid(0.000706231968588, maxQuantum, nGridPoints), ModeGrid(0.001307668225708, maxQuantum, nGridPoints),
        ModeGrid(0.001895435476984, maxQuantum, nGridPoints), ModeGrid(0.001963780506202, maxQuantum, nGridPoints),
        ModeGrid(0.002843153215476, maxQuantum, nGridPoints), ModeGrid(0.002875047562445, maxQuantum, nGridPoints),
        ModeGrid(0.003959455359373, maxQuantum, nGridPoints), ModeGrid(0.004487990251994, maxQuantum, nGridPoints),
        ModeGrid(0.004720363351336, maxQuantum, nGridPoints), ModeGrid(0.005057532162146, maxQuantum, nGridPoints),
        ModeGrid(0.005175996879457, maxQuantum, nGridPoints), ModeGrid(0.005294461596769, maxQuantum, nGridPoints),
        ModeGrid(0.005745538789609, maxQuantum, nGridPoints), ModeGrid(0.006310524364479, maxQuantum, nGridPoints),
        ModeGrid(0.006438101752353, maxQuantum, nGridPoints), ModeGrid(0.006593017151914, maxQuantum, nGridPoints),
        ModeGrid(0.006665918516414, maxQuantum, nGridPoints), ModeGrid(0.006693256528101, maxQuantum, nGridPoints),
        ModeGrid(0.006720594539788, maxQuantum, nGridPoints), ModeGrid(0.007080545027004, maxQuantum, nGridPoints),
        ModeGrid(0.008105720465277, maxQuantum, nGridPoints), ModeGrid(0.013751019878699, maxQuantum, nGridPoints),
        ModeGrid(0.013860371925448, maxQuantum, nGridPoints), ModeGrid(0.013992505648603, maxQuantum, nGridPoints),
        ModeGrid(0.014247660424351, maxQuantum, nGridPoints), ModeGrid(0.014252216759632, maxQuantum, nGridPoints),
        ModeGrid(0.014352456135819, maxQuantum, nGridPoints), ModeGrid(0.016539497070800, maxQuantum, nGridPoints)};
    nmaPesCalculator_ = std::make_shared<NMAPesCalculator>();
    auto refGeom = Utils::XyzStreamHandler::read(nmaOpt);
    nmaPesCalculator_->setStructure(refGeom);
    pesCartesianNma_ = std::make_shared<PesCartesian>(nmaPesCalculator_);
    //
    int nGridPointsDVR = 30;
    int maxQuantumDVR = 10;
    modeGridH2O = {ModeGrid(0.00755852692, maxQuantumDVR, nGridPointsDVR),
                   ModeGrid(0.01709230311, maxQuantumDVR, nGridPointsDVR),
                   ModeGrid(0.01755421162, maxQuantumDVR, nGridPointsDVR)};
    refGeom = Utils::XyzStreamHandler::read(h2oOpt);
    h2oPesCalculator_ = std::make_shared<H2OHarmonicPesCalculator>();
    h2oPesCalculator_->setStructure(refGeom);
    pesCartesianH2O_ = std::make_shared<PesCartesian>(h2oPesCalculator_);
  }
};

/**
 * @brief Checks that the maximum displacement for a given grid is proportional
 * to the square root of the maximum quantum of excitation (+1/2).
 */
TEST_F(VibrationalUtilsTest, CheckGridQMax) {
  initializeMPI();
  auto qMax10 = modeGrid100CmM1With10NMax.getQMax();
  auto qMax20 = modeGrid100CmM1With20NMax.getQMax();
  auto ratio = qMax20 / qMax10;
  EXPECT_THAT(ratio, DoubleNear(std::sqrt(20.5 / 10.5), 1.0E-15));
}

/**
 * @brief Checks that the basis size is correctly gotten through the dvr basis
 */
TEST_F(VibrationalUtilsTest, CheckBasisSetDimension) {
  std::vector<ModeGrid> tmpGridVec = {modeGrid100CmM1With10NMax};
  auto dvrBasis = DVR(tmpGridVec, pesCartesianNma_);
  ASSERT_EQ(dvrBasis.getBasisSize(0), 100);
  tmpGridVec = {modeGrid100CmM1With10NMax, modeGrid100CmM1With20NMax};
  dvrBasis = DVR(tmpGridVec, pesCartesianNma_);
  ASSERT_EQ(dvrBasis.getBasisSize(1), 100);
}

/** @brief Check for one-dim numerical integration */
TEST_F(VibrationalUtilsTest, oneDimIntegration) {
  auto f = [](double /*x*/) { return 1.; };
  GaussHermiteQuad quad;
  auto res = quad.oneDimIntegrate(f, 8);
  EXPECT_THAT(res, DoubleNear(std::sqrt(M_PI), 1.0E-8));
}

/** @brief Check for two-dim numerical integration */
TEST_F(VibrationalUtilsTest, twoDimIntegration) {
  auto f = [](double x, double y) { return pow(x / sqrt(3.) + 1., 2) * pow(y - 2., 2); };
  GaussHermiteQuad quad;
  double res = 1. / sqrt(3.) * quad.twoDimIntegrate(f, 8);
  EXPECT_THAT(res, DoubleNear(9.522446662228603, 1.0E-5));
}

/**
 * @brief Checks that the one-body DVR matrix is symmetric.
 * Note that this test does not ensure that the calculation of the Hamiltonian
 * matrix (per se) is correct.
 */
TEST_F(VibrationalUtilsTest, CheckConsistencyOneBodyDVR) {
  std::vector<ModeGrid> tmpGridVec = {modeGrid100CmM1With10NMax};
  auto dvrBasis = DVR(tmpGridVec, pesCartesianNma_);
  int size = dvrBasis.getBasisSize(0);
  Eigen::MatrixXd hamiltonianMatrix = Eigen::MatrixXd::Zero(size, size);
  for (int iRow = 0; iRow < size; iRow++) {
    for (int iCol = 0; iCol < size; iCol++) {
      hamiltonianMatrix(iRow, iCol) = dvrBasis.getOneBodyIntegrals(iRow, iCol, 0);
    }
  }
  for (int iRow = 0; iRow < size; iRow++) {
    for (int iCol = iRow + 1; iCol < size; iCol++) {
      EXPECT_THAT(hamiltonianMatrix(iRow, iCol), DoubleNear(hamiltonianMatrix(iCol, iRow), 1.0E-16));
    }
  }
}

/**
 * @brief Checks that the one-body DG matrix is symmetric.
 * Note that this test does not ensure that the calculation of the Hamiltonian
 * matrix (per se) is correct.
 */
TEST_F(VibrationalUtilsTest, CheckConsistencyOneBodyDG) {
  std::vector<ModeGrid> tmpGridVec = {modeGrid100CmM1With10NMax};
  auto dgBasis = DistributedGaussians(tmpGridVec, pesCartesianNma_, 16);
  int size = dgBasis.getBasisSize(0);
  Eigen::MatrixXd hamiltonianMatrix = Eigen::MatrixXd::Zero(size, size);
  for (int iRow = 0; iRow < size; iRow++) {
    for (int iCol = 0; iCol < size; iCol++) {
      hamiltonianMatrix(iRow, iCol) = dgBasis.getOneBodyIntegrals(iRow, iCol, 0);
    }
  }
  for (int iRow = 0; iRow < size; iRow++) {
    for (int iCol = iRow + 1; iCol < size; iCol++) {
      EXPECT_THAT(hamiltonianMatrix(iRow, iCol), DoubleNear(hamiltonianMatrix(iCol, iRow), 1.0E-16));
    }
  }
}

/**
 * @brief Checks that the DVR basis is orthogonal
 */
TEST_F(VibrationalUtilsTest, CheckOverlapForDVRBasis) {
  std::vector<ModeGrid> tmpGridVec = {modeGrid100CmM1With10NMax};
  auto dvrBasis = DVR(tmpGridVec, pesCartesianNma_);
  for (int iRow = 0; iRow < 100; iRow++) {
    for (int iCol = iRow + 1; iCol < 100; iCol++) {
      double overlap = dvrBasis.getOverlap(iRow, iCol, 0);
      EXPECT_THAT(overlap, DoubleNear(0., 1.0E-16));
    }
  }
}

/** @brief Checks that the reference energy of H2O is 0. (it's a parabola) */
TEST_F(VibrationalUtilsTest, CheckReferenceEnergyOfH2O) {
  double referenceEnergy = pesCartesianH2O_->getPES();
  EXPECT_THAT(referenceEnergy, DoubleNear(0., 1.0E-10));
}

/**
 * @brief Last excitation energy of NMA with a DVR basis.
 */
TEST_F(VibrationalUtilsTest, CheckLastExcitationOfNMAWithOneBody) {
  getRank();
  std::vector<ModeGrid> tmpGridVec = std::vector<ModeGrid>(30, modeGridLastModeNma);
  int idxMode = 29;
  auto dvrBasis = DVR(tmpGridVec, pesCartesianNma_);
  int size = dvrBasis.getBasisSize(idxMode);
  Eigen::MatrixXd hamiltonianMatrix = Eigen::MatrixXd::Zero(size, size);
  double referenceEnergy = pesCartesianNma_->getPES();
  for (int iRow = 0; iRow < size; iRow++) {
    for (int iCol = 0; iCol < size; iCol++) {
      hamiltonianMatrix(iRow, iCol) = dvrBasis.getOneBodyIntegrals(iRow, iCol, idxMode);
    }
  }
  // Diagonalization of the Hamiltonian
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  es.compute(hamiltonianMatrix);
  Eigen::VectorXd eigenvalues = es.eigenvalues();
  std::vector<double> diffs;
  if (id == 0) {
    std::cout << "Modal excitation energies of one-body NMA Hamiltonian of mode 29" << std::endl;
  }
  for (int idx = 1; idx < size; idx++) {
    if (id == 0) {
      std::cout << (eigenvalues(idx - 1) - referenceEnergy) * 219474.63;
    }
    double tmpval = (eigenvalues(idx) - eigenvalues(idx - 1)) * 219474.63;
    if (id == 0) {
      std::cout << " " << tmpval << std::endl;
    }
    diffs.push_back(tmpval);
  }
  // Although one-body potential is not strictly harmonic, we'd expect close to
  // equidistant modal energies
  EXPECT_THAT(diffs[1], DoubleNear(diffs[0], 300));
}

/**
 * @brief Excitation frequencies of H2O
 */
TEST_F(VibrationalUtilsTest, CheckExcitationsWithHarmonicH2OPES) {
  auto dvrBasis = DVR(modeGridH2O, pesCartesianH2O_);
  std::vector<double> vibrationalExcitations(3);
  for (int iMode = 0; iMode < 3; iMode++) {
    int size = dvrBasis.getBasisSize(iMode);
    Eigen::MatrixXd hamiltonianMatrix = Eigen::MatrixXd::Zero(size, size);
    for (int iRow = 0; iRow < size; iRow++) {
      for (int iCol = 0; iCol < size; iCol++) {
        hamiltonianMatrix(iRow, iCol) = dvrBasis.getOneBodyIntegrals(iRow, iCol, iMode);
      }
    }
    // Diagonalization of the Hamiltonian
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute(hamiltonianMatrix);
    Eigen::VectorXd eigenvalues = es.eigenvalues();
    vibrationalExcitations[iMode] = (eigenvalues(1) - eigenvalues(0)) * 219474.63;
  }
  // Although one-body potential is not strictly harmonic, we'd expect it to be
  // close
  EXPECT_THAT(vibrationalExcitations[0], DoubleNear(1658.83, 1.0E+2));
  EXPECT_THAT(vibrationalExcitations[1], DoubleNear(3751.17, 1.0E+2));
  EXPECT_THAT(vibrationalExcitations[2], DoubleNear(3852.53, 2.0E+2));
}

TEST_F(VibrationalUtilsTest, CheckConsistencyDistributedGaussianAndDvr) {
  getRank();
  // Sets up the DVR Hamiltonian
  auto dvrBasis = DVR(modeGridNMA, pesCartesianNma_);
  auto dgBasis = DistributedGaussians(modeGridNMA, pesCartesianNma_, 16);
  double au2cm1 = 219474.63;
  std::vector<int> lstModesToCheck = {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14,
                                      15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};
  for (const auto& iMode : lstModesToCheck) {
    if (id == 0) {
      std::cout << "Mode " << iMode << std::endl;
    }
    auto size = modeGridNMA[iMode].getSize();
    Eigen::MatrixXd hamiltonianMatrixDVR = Eigen::MatrixXd::Zero(size, size);
    Eigen::MatrixXd hamiltonianMatrixDG = Eigen::MatrixXd::Zero(size, size);
    Eigen::MatrixXd overlapMatrixDG = Eigen::MatrixXd::Zero(size, size);
    for (int iRow = 0; iRow < size; iRow++) {
      for (int iCol = 0; iCol < size; iCol++) {
        hamiltonianMatrixDVR(iRow, iCol) = dvrBasis.getOneBodyIntegrals(iRow, iCol, iMode);
        hamiltonianMatrixDG(iRow, iCol) = dgBasis.getOneBodyIntegrals(iRow, iCol, iMode);
        overlapMatrixDG(iRow, iCol) = dgBasis.getOverlap(iRow, iCol, iMode);
      }
    }
    // Diagonalization of the Hamiltonian
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esDVR;
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> esDG;
    esDVR.compute(hamiltonianMatrixDVR);
    esDG.compute(hamiltonianMatrixDG, overlapMatrixDG);
    auto energyDiffDVR = (esDVR.eigenvalues()(1) - esDVR.eigenvalues()(0)) * au2cm1;
    auto energyDiffDG = (esDG.eigenvalues()(1) - esDG.eigenvalues()(0)) * au2cm1;
    if (id == 0) {
      std::cout.precision(4);
      std::cout << std::fixed;
      std::cout << " Distributed Gaussian " << energyDiffDG << std::endl;
      std::cout << " DVR " << energyDiffDVR << std::endl;
    }
    double energyRatioDGtoDVR = energyDiffDG / energyDiffDVR;
    EXPECT_THAT(energyRatioDGtoDVR, DoubleNear(1, 1.0E-2)); // Allow a proportional error of up to 1 percent
  }
}

TEST_F(VibrationalUtilsTest, CheckConsistencyDistributedGaussianAndDvrForH2O) {
  getRank();
  // Sets up the DVR Hamiltonian
  auto dvrBasis = DVR(modeGridH2O, pesCartesianH2O_);
  auto dgBasis = DistributedGaussians(modeGridH2O, pesCartesianH2O_, 16);
  double au2cm1 = 219474.63;
  std::vector<int> lstModesToCheck = {0, 1, 2};
  for (const auto& iMode : lstModesToCheck) {
    if (id == 0) {
      std::cout << "Mode " << iMode << std::endl;
    }
    auto size = modeGridH2O[iMode].getSize();
    Eigen::MatrixXd hamiltonianMatrixDVR = Eigen::MatrixXd::Zero(size, size);
    Eigen::MatrixXd hamiltonianMatrixDG = Eigen::MatrixXd::Zero(size, size);
    Eigen::MatrixXd overlapMatrixDG = Eigen::MatrixXd::Zero(size, size);
    for (int iRow = 0; iRow < size; iRow++) {
      for (int iCol = 0; iCol < size; iCol++) {
        hamiltonianMatrixDVR(iRow, iCol) = dvrBasis.getOneBodyIntegrals(iRow, iCol, iMode);
        hamiltonianMatrixDG(iRow, iCol) = dgBasis.getOneBodyIntegrals(iRow, iCol, iMode);
        overlapMatrixDG(iRow, iCol) = dgBasis.getOverlap(iRow, iCol, iMode);
      }
    }
    // Diagonalization of the Hamiltonian
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esDVR;
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> esDG;
    esDVR.compute(hamiltonianMatrixDVR);
    esDG.compute(hamiltonianMatrixDG, overlapMatrixDG);
    auto energyDiffDVR = (esDVR.eigenvalues()(1) - esDVR.eigenvalues()(0)) * au2cm1;
    auto energyDiffDG = (esDG.eigenvalues()(1) - esDG.eigenvalues()(0)) * au2cm1;
    EXPECT_THAT(energyDiffDVR, DoubleNear(energyDiffDG, 1.0E-3));
  }
  finalizeMPI();
}

} // namespace Scine::Colibri
