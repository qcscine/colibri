/**
 * @file ModuleTest.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#include <Core/Interfaces/Calculator.h>
#include <Core/ModuleManager.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/Settings.h>
#include <gmock/gmock.h>
#include <boost/filesystem.hpp>

#ifdef MPI_PARALLEL
  #include <mpi.h>
#endif

namespace Scine::Colibri {

using namespace testing;

class ModuleTest : public Test {
 protected:
  // Members
  std::ofstream inputFileCO2;
  std::string inputFileCO2Name;
  int size = 1, id = 0; // Needed for potentially parallel MPI run

  void mpiSetup() {
#ifdef MPI_PARALLEL
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  }

  void mpiTeardown() {
#ifdef MPI_PARALLEL
    MPI_Finalize();
#endif
  }

  void writeFiles() {
#ifdef MPI_PARALLEL
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    inputFileCO2Name = "co2.ini";
    if (id == 0) {
      // Input file for the CO2 calculation
      inputFileCO2.open(inputFileCO2Name);
      inputFileCO2 << "numModes = 4 " << std::endl;
      inputFileCO2 << "vscfIter = 10 " << std::endl;
      inputFileCO2 << "ONVector = 0 0 0 0 " << std::endl;
      inputFileCO2 << "vscfEnTol = 10E-8 " << std::endl;
      inputFileCO2 << "vscfCoeffTol = 10E-12 " << std::endl;
      inputFileCO2 << "nModePotentialOrder = 2 " << std::endl;
      inputFileCO2 << "nMax = 8 " << std::endl;
      inputFileCO2 << "twoBodyTol = 1E-8 " << std::endl;
      inputFileCO2 << "fciDumpName = co2_test " << std::endl;
      inputFileCO2 << "numqp = 8 " << std::endl;
      inputFileCO2 << "[Mode0] " << std::endl;
      inputFileCO2 << "numGrid = 10 " << std::endl;
      inputFileCO2 << "frequency = 0.00310191311 " << std::endl;
      inputFileCO2 << "nquantum = 5 " << std::endl;
      inputFileCO2 << "[Mode1] " << std::endl;
      inputFileCO2 << "numGrid = 10 " << std::endl;
      inputFileCO2 << "frequency = 0.00310191311 " << std::endl;
      inputFileCO2 << "nquantum = 5 " << std::endl;
      inputFileCO2 << "[Mode2] " << std::endl;
      inputFileCO2 << "numGrid = 10 " << std::endl;
      inputFileCO2 << "frequency = 0.00611968448 " << std::endl;
      inputFileCO2 << "nquantum = 5 " << std::endl;
      inputFileCO2 << "[Mode3] " << std::endl;
      inputFileCO2 << "numGrid = 10 " << std::endl;
      inputFileCO2 << "frequency = 0.0108343253 " << std::endl;
      inputFileCO2 << "nquantum = 5 " << std::endl;
      inputFileCO2 << "[ForceConstants] " << std::endl;
      inputFileCO2 << "firstOrder = "
                   << "0  0.00000000000000E+00 "
                   << "1  0.00000000000000E+00 "
                   << "2 -3.79376328358405E-05 "
                   << "3  0.00000000000000E+00 " << std::endl;
      inputFileCO2 << "secondOrder = "
                   << "0 0 9.62186495432437E-06 "
                   << "1 1 9.62186495432437E-06 "
                   << "2 2 3.74505381530253E-05 "
                   << "3 3 0.00011738260475E+00 " << std::endl;
      inputFileCO2 << "thirdOrder = "
                   << "2 2 2  5.83894731065161E-07 "
                   << "2 0 0 -1.60960891738734E-07 "
                   << "2 1 1 -1.60960891738734E-07 "
                   << "2 3 3  1.95710501210382E-06 " << std::endl;
      inputFileCO2 << "fourthOrder = "
                   << "0 0 0 0  2.8412590870909E-09 "
                   << "1 1 1 1  2.8412590870909E-09 "
                   << "2 2 2 2  7.5342922652451E-09 "
                   << "3 3 3 3  8.8291418176427E-08 "
                   << "2 2 0 0 -3.8347558276766E-10 "
                   << "2 2 1 1 -3.8347558276766E-10 "
                   << "2 2 3 3  2.5417941550711E-08 "
                   << "0 0 1 1  9.8448211457855E-10 "
                   << "0 0 3 3 -1.6636940667766E-08 "
                   << "1 1 3 3 -1.6636940667766E-08 " << std::endl;
      inputFileCO2.close();
    }
#ifdef MPI_PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
  void removeFiles() {
#ifdef MPI_PARALLEL
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (id == 0) {
      auto successCO2 = std::remove(inputFileCO2Name.c_str());
      auto successFCIDUMP = std::remove("co2_test");
      if (successCO2 != 0 || successFCIDUMP != 0) {
        throw std::runtime_error("Error while removing input files");
      }
    }
#ifdef MPI_PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
};

TEST_F(ModuleTest, CalculateCO2ThroughInterfaceDG) {
  mpiSetup();
  writeFiles();
  auto& managerRef = Core::ModuleManager::getInstance();
  if (!managerRef.moduleLoaded("Colibri")) {
    managerRef.load("../libcolibri.so");
  }
  auto vscfCalc = managerRef.get<Core::Calculator>("PesQffDistGauss");
  auto results = vscfCalc->calculate(inputFileCO2Name);
  double vscfEnergy = results.get<Utils::Property::Energy>();
  if (id == 0) {
    EXPECT_THAT(vscfEnergy, DoubleNear(0.0115909, 1.0E-4));
  }
  removeFiles();
}

TEST_F(ModuleTest, CalculateCO2ThroughInterfaceDVR) {
  writeFiles();
  auto& managerRef = Core::ModuleManager::getInstance();
  if (!managerRef.moduleLoaded("Colibri")) {
    managerRef.load("../libcolibri.so");
  }
  auto vscfCalc = managerRef.get<Core::Calculator>("PesQffDVR");
  auto results = vscfCalc->calculate(inputFileCO2Name);
  double vscfEnergy = results.get<Utils::Property::Energy>();
  if (id == 0) {
    EXPECT_THAT(vscfEnergy, DoubleNear(0.0115909, 1.0E-4));
  }
  removeFiles();
  mpiTeardown();
}

} // namespace Scine::Colibri