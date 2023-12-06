/**
 * @file OTFVSCFTest.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Core/ModuleManager.h>
#include <VSCF/VSCF.h>
#include <VibrationalUtils/PESLibrary/HybridOTFfromSPs.h>
#include <VibrationalUtils/PESLibrary/OnTheFlyPes.h>
#include <VibrationalUtils/PESLibrary/PesFromSPs.h>
#include <VibrationalUtils/VibrationalBases/DVR.h>
#include <VibrationalUtils/VibrationalParameters.h>
#include <gmock/gmock.h>
#include <cmath>
#include <cstdlib>
#include <string>

#ifdef MPI_PARALLEL
  #include <mpi.h>
#endif

#include "writeOTFInpFiles.h"

namespace Scine::Colibri {

using namespace testing;
using namespace OTFVSCFTestInputs;

class OTFVSCFTest : public Test {
 protected:
  std::string inputFileEthene = "tmp_ethene_otf.inp";
  std::string inputXYZEthene = "tmp_ethene_struct.xyz";
  std::string inputFileH2 = "tmp_h2_otf.inp";
  std::string inputXYZH2 = "tmp_h2_struct.xyz";
  std::string inputFileWater = "tmp_water_otf.inp";
  std::string inputXYZWater = "tmp_water_struct.xyz";
  std::string harmfreqs = "harmfreqs.txt";

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

  void getMPIrank() {
    id = 0;
#ifdef MPI_PARALLEL
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  void createInputFiles(const std::string& molecule) {
#ifdef MPI_PARALLEL
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (id == 0) {
      if (molecule == "ethene") {
        writeInputEthene(inputFileEthene);
        writeXYZEthene(inputXYZEthene);
        writeHarmFreqsEthene(harmfreqs);
      } else if (molecule == "h2") {
        writeInputH2(inputFileH2);
        writeXYZforH2(inputXYZH2);
        writeHarmFreqsH2(harmfreqs);
      } else if (molecule == "water") {
        writeInputH2O(inputFileWater);
        writeXYZforH20(inputXYZWater);
        writeHarmFreqsH2O(harmfreqs);
      }
    }
#ifdef MPI_PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  void removeFiles(const std::string& molecule) {
#ifdef MPI_PARALLEL
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (id == 0) {
      if (molecule == "ethene") {
        auto successEtheneInp = std::remove(inputFileEthene.c_str());
        auto successFcidumpOut = std::remove("tmp_FCIDUMPv_otf2");
        auto successEtheneXYZ = std::remove(inputXYZEthene.c_str());
        auto successPes1ModeCuts = std::remove("pes1ModeCuts.out");
        auto successPes2ModeCuts = std::remove("pes2ModeCuts.out");
        auto successPesCuts = std::remove("pesCuts.out");
        auto successHarmFreqs = std::remove(harmfreqs.c_str());
        if ((successEtheneInp != 0) || (successFcidumpOut != 0) || (successEtheneXYZ != 0) || (successPes1ModeCuts != 0) ||
            (successPes2ModeCuts != 0) || (successPesCuts != 0) || (successHarmFreqs != 0)) {
          throw std::runtime_error("Error while removing input files");
        }
      } else if (molecule == "h2") {
        auto successH2Inp = std::remove(inputFileH2.c_str());
        auto successFcidumpH2 = std::remove("FCIDUMPv_otf1");
        auto successH2XYZ = std::remove(inputXYZH2.c_str());
        auto successPes1ModeCuts = std::remove("pes1ModeCuts.out");
        auto successHarmFreqs = std::remove(harmfreqs.c_str());
        if ((successH2Inp != 0) || (successFcidumpH2 != 0) || (successH2XYZ != 0) || (successPes1ModeCuts != 0) ||
            (successHarmFreqs != 0)) {
          throw std::runtime_error("Error while removing input files");
        }
      } else if (molecule == "water") {
        auto successInp = std::remove(inputFileWater.c_str());
        auto successFcidump = std::remove("tmp_FCIDUMPv_h2o");
        auto successXYZ = std::remove(inputXYZWater.c_str());
        auto successPes1ModeCuts = std::remove("pes1ModeCuts.out");
        auto successPes2ModeCuts = std::remove("pes2ModeCuts.out");
        auto successPes12ModeCuts = std::remove("pes12ModeCuts.out");
        auto successPes3ModeCuts = std::remove("pes3ModeCuts.out");
        auto successPes123ModeCuts = std::remove("pes123ModeCuts.out");
        auto successRefXYZ = std::remove("h2o_ref.xyz");
        auto successHarmFreqs = std::remove(harmfreqs.c_str());
        if ((successInp != 0) || (successFcidump != 0) || (successXYZ != 0) || (successRefXYZ != 0) ||
            (successPes1ModeCuts != 0) || (successPes2ModeCuts != 0) || (successPes12ModeCuts != 0) ||
            (successPes3ModeCuts != 0) || (successPes123ModeCuts != 0) || (successHarmFreqs != 0)) {
          throw std::runtime_error("Error while removing input files");
        }
      }
    }
#ifdef MPI_PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
};

TEST_F(OTFVSCFTest, OTFVSCF_Ethene_Sparrow) {
  initializeMPI();
  auto& manager = Core::ModuleManager::getInstance();
  if (manager.moduleLoaded("Sparrow")) {
    createInputFiles("ethene");
    VibrationalParameters parms;
    parms.setParametersOTF(inputFileEthene);

    std::shared_ptr<PES> pes = std::make_shared<OTFPes>(parms);
    std::shared_ptr<ModalBasis> basis = std::make_shared<DVR>(parms, pes);
    std::shared_ptr<VibrationalSCF> vscf = std::make_shared<VSCF>(parms, basis);

    if (id == 0) {
      EXPECT_EQ(pes->getNumSPsForVSCF(), 6720);
    }

    // concatenate pesCuts as input for next test
    std::ifstream pes1mode("pes1ModeCuts.out", std::ios_base::binary);
    std::ifstream pes2mode("pes2ModeCuts.out", std::ios_base::binary);
    std::ofstream pesCuts("pesCuts.out", std::ios_base::binary);

    pesCuts << pes1mode.rdbuf() << pes2mode.rdbuf();

    pes1mode.close();
    pes2mode.close();
    pesCuts.close();

    double vscfEnergy = NAN;
    if (id == 0) {
      vscfEnergy = vscf->SCF() * 219474.63;
    }
    if (id == 0) {
      EXPECT_THAT(vscfEnergy, DoubleNear(10868.5932158505, 1e-05));
    }
  } else {
    if (id == 0) {
      std::cout << "Sparrow OTFVSCF test was not run as the Sparrow module is "
                   "not loaded."
                << std::endl;
    }
  }
}

TEST_F(OTFVSCFTest, VSCF_Ethene_fromPreviousSparrowSPs) {
  getMPIrank();
  auto& manager = Core::ModuleManager::getInstance();
  if (manager.moduleLoaded("Sparrow")) {
    VibrationalParameters parms;
    parms.setParametersAlsoVCI(inputFileEthene);
    parms.singlePointsFname_ = "pesCuts.out";
    parms.refEnergy_ = -4.9007076061;
    if (id == 0) {
      parms.printParametersFromSPs();
    }

    std::shared_ptr<PES> pes = std::make_shared<PesFromSPs>(parms);
    std::shared_ptr<ModalBasis> basis = std::make_shared<DVR>(parms, pes);
    std::shared_ptr<VibrationalSCF> vscf = std::make_shared<VSCF>(parms, basis);

    if (id == 0) {
      EXPECT_EQ(pes->getNumSPsForVSCF(), 0);
    }

    double vscfEnergy = NAN;
    if (id == 0) {
      vscfEnergy = vscf->SCF() * 219474.63;
    }
    if (id == 0) {
      EXPECT_THAT(vscfEnergy, DoubleNear(10868.5932158505, 1e-04));
    }
    removeFiles("ethene");
  } else {
    if (id == 0) {
      std::cout << "Sparrow OTFVSCF test was not run as the Sparrow module is "
                   "not loaded."
                << std::endl;
    }
  }
}

TEST_F(OTFVSCFTest, OTFVSCF_H2_Turbomole) {
  getMPIrank();
  auto& manager = Core::ModuleManager::getInstance();
  if (manager.moduleLoaded("Turbomole")) {
    createInputFiles("h2");
    VibrationalParameters parms;
    parms.setParametersOTF(inputFileH2);

    std::shared_ptr<PES> pes = std::make_shared<OTFPes>(parms);
    std::shared_ptr<ModalBasis> basis = std::make_shared<DVR>(parms, pes);
    std::shared_ptr<VibrationalSCF> vscf = std::make_shared<VSCF>(parms, basis);

    if (id == 0) {
      EXPECT_EQ(pes->getNumSPsForVSCF(), 10);
    }

    double vscfEnergy = NAN;
    if (id == 0) {
      vscfEnergy = vscf->SCF() * 219474.63;
    }
    if (id == 0) {
      EXPECT_THAT(vscfEnergy, DoubleNear(2152.0750247592, 1e-02));
    }
  } else {
    if (id == 0) {
      std::cout << "Turbomole OTFVSCF test was not run as the Turbomole module "
                   "is not loaded."
                << std::endl;
    }
  }
}

TEST_F(OTFVSCFTest, VSCF_H2_fromPreviousTurbomoleSPs) {
  getMPIrank();
  auto& manager = Core::ModuleManager::getInstance();
  if (manager.moduleLoaded("Sparrow")) {
    VibrationalParameters parms;
    parms.setParametersAlsoVCI(inputFileH2);
    parms.singlePointsFname_ = "pes1ModeCuts.out";
    parms.refEnergy_ = -1.1675754883;
    if (id == 0) {
      parms.printParametersFromSPs();
    }

    std::shared_ptr<PES> pes = std::make_shared<PesFromSPs>(parms);
    std::shared_ptr<ModalBasis> basis = std::make_shared<DVR>(parms, pes);
    std::shared_ptr<VibrationalSCF> vscf = std::make_shared<VSCF>(parms, basis);

    if (id == 0) {
      EXPECT_EQ(pes->getNumSPsForVSCF(), 0);
    }

    double vscfEnergy = NAN;
    if (id == 0) {
      vscfEnergy = vscf->SCF() * 219474.63;
    }
    if (id == 0) {
      EXPECT_THAT(vscfEnergy, DoubleNear(2152.0750247592, 1e-02));
    }
    removeFiles("h2");
  } else {
    if (id == 0) {
      std::cout << "Sparrow OTFVSCF test was not run as the Sparrow module is "
                   "not loaded."
                << std::endl;
    }
  }
}

TEST_F(OTFVSCFTest, OTFVSCF_Water_Sparrow) {
  getMPIrank();
  auto& manager = Core::ModuleManager::getInstance();
  if (manager.moduleLoaded("Sparrow")) {
    createInputFiles("water");
    VibrationalParameters parms;
    parms.setParametersOTF(inputFileWater);
    parms.nModePotentialOrder_ = 2;

    std::shared_ptr<PES> pes = std::make_shared<OTFPes>(parms);
    std::shared_ptr<ModalBasis> basis = std::make_shared<DVR>(parms, pes);
    std::shared_ptr<VibrationalSCF> vscf = std::make_shared<VSCF>(parms, basis);

    // concatenate pesCuts as input for next test
    std::ifstream pes1mode("pes1ModeCuts.out", std::ios_base::binary);
    std::ifstream pes2mode("pes2ModeCuts.out", std::ios_base::binary);
    std::ofstream pesCuts("pes12ModeCuts.out", std::ios_base::binary);

    pesCuts << pes1mode.rdbuf() << pes2mode.rdbuf();

    pes1mode.close();
    pes2mode.close();
    pesCuts.close();

    std::ofstream refXYZ("h2o_ref.xyz");
    Utils::XyzStreamHandler::write(refXYZ, pes->getReferenceGeometry());
    refXYZ.close();

    if (id == 0) {
      EXPECT_EQ(pes->getNumSPsForVSCF(), 330);
    }

    double vscfEnergy = NAN;
    if (id == 0) {
      vscfEnergy = vscf->SCF() * 219474.63;
    }
    if (id == 0) {
      EXPECT_THAT(vscfEnergy, DoubleNear(4432.14, 1e-02));
    }
  } else {
    if (id == 0) {
      std::cout << "Sparrow OTFVSCF test was not run as the Sparrow module is "
                   "not loaded."
                << std::endl;
    }
  }
}

TEST_F(OTFVSCFTest, VSCF_Water_HybridFromPreviousAndOTFSparrowSPs) {
  getMPIrank();
  auto& manager = Core::ModuleManager::getInstance();
  if (manager.moduleLoaded("Sparrow")) {
    VibrationalParameters parms;
    parms.setParametersAlsoVCI(inputFileWater);
    parms.singlePointsFname_ = "pes12ModeCuts.out";
    parms.xyzFname_ = "h2o_ref.xyz";
    parms.refEnergy_ = -4.0715756440;
    parms.nModePotentialOrder_ = 3;

    if (id == 0) {
      parms.printParametersHybrid();
    }

    std::shared_ptr<PES> pes = std::make_shared<HybridPes>(parms);
    std::shared_ptr<ModalBasis> basis = std::make_shared<DVR>(parms, pes);
    std::shared_ptr<VibrationalSCF> vscf = std::make_shared<VSCF>(parms, basis);

    // concatenate pesCuts as input for next test
    std::ifstream pes12mode("pes12ModeCuts.out", std::ios_base::binary);
    std::ifstream pes3mode("pes3ModeCuts.out", std::ios_base::binary);
    std::ofstream pesCuts("pes123ModeCuts.out", std::ios_base::binary);

    pesCuts << pes12mode.rdbuf() << pes3mode.rdbuf();

    pes12mode.close();
    pes3mode.close();
    pesCuts.close();

    if (id == 0) {
      EXPECT_EQ(pes->getNumSPsForVSCF(), 1000);
    }

    double vscfEnergy = NAN;
    if (id == 0) {
      vscfEnergy = vscf->SCF() * 219474.63;
    }

    if (id == 0) {
      EXPECT_THAT(vscfEnergy, DoubleNear(4431.32, 1e-02)); // ref from OTF calc
    }
  } else {
    if (id == 0) {
      std::cout << "Sparrow OTFVSCF test was not run as the Sparrow module is "
                   "not loaded."
                << std::endl;
    }
  }
}

TEST_F(OTFVSCFTest, VSCF_Water_fromPreviousSparrowSPs) {
  getMPIrank();
  auto& manager = Core::ModuleManager::getInstance();
  if (manager.moduleLoaded("Sparrow")) {
    VibrationalParameters parms;
    parms.setParametersAlsoVCI(inputFileWater);
    parms.refEnergy_ = -4.0715756440;
    parms.nModePotentialOrder_ = 3;
    parms.singlePointsFname_ = "pes123ModeCuts.out";
    if (id == 0) {
      parms.printParametersFromSPs();
    }

    std::shared_ptr<PES> pes = std::make_shared<PesFromSPs>(parms);
    std::shared_ptr<ModalBasis> basis = std::make_shared<DVR>(parms, pes);
    std::shared_ptr<VibrationalSCF> vscf = std::make_shared<VSCF>(parms, basis);

    if (id == 0) {
      EXPECT_EQ(pes->getNumSPsForVSCF(), 0);
    }

    double vscfEnergy = NAN;
    if (id == 0) {
      vscfEnergy = vscf->SCF() * 219474.63;
    }

    removeFiles("water");

    if (id == 0) {
      EXPECT_THAT(vscfEnergy, DoubleNear(4431.32, 1e-02)); // ref from OTF calc
    }
  } else {
    if (id == 0) {
      std::cout << "Sparrow OTFVSCF test was not run as the Sparrow module is "
                   "not loaded."
                << std::endl;
    }
  }
  finalizeMPI();
}

} // namespace Scine::Colibri
