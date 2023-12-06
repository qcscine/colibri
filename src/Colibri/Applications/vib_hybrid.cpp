/**
 * @file vib_hybrid.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#include "Utils/Constants.h"
#include <VCI/VCI.h>
#include <VCI/VibrationalCI.h>
#include <VSCF/VSCF.h>
#include <VSCF/VibrationalSCF.h>
#include <VibrationalUtils/PESLibrary/HybridOTFfromSPs.h>
#include <VibrationalUtils/PESLibrary/PES.h>
#include <VibrationalUtils/VibrationalBases/DVR.h>
#include <VibrationalUtils/VibrationalParameters.h>
#include <sys/time.h>
#include <iostream>

#ifdef MPI_PARALLEL
  #include <mpi.h>
#endif

using namespace Scine;
using namespace Colibri;

/**************************************
 * Main routine for VCI calculation with prior VSCF
 **************************************/

int main(int argc, char** argv) {
#ifdef MPI_PARALLEL
  MPI_Init(&argc, &argv);
  int size, id;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
#else
  int id = 0;
#endif

  if (id == 0) {
    std::cout << "------ VCI with prior VSCF routine ------" << std::endl << std::endl;
  }

  timeval start, first, then, before_vci, after_vci;
  gettimeofday(&start, nullptr);
  std::string inputFileName = std::string(argv[1]);
  if (id == 0) {
    std::cout << "Input file: " << inputFileName << std::endl << std::endl;
  }

  std::cout << std::setprecision(8);

  try {
    if (id == 0) {
      std::cout << "------ First step: VSCF routine ------" << std::endl << std::endl;
    }
    VibrationalParameters parms;
    parms.setParametersAlsoVCI(inputFileName);
    if (id == 0) {
      parms.printParametersHybrid();
    }

    std::shared_ptr<PES> pes = std::make_shared<HybridPes>(parms);
    std::shared_ptr<ModalBasis> basis = std::make_shared<DVR>(parms, pes);
    std::shared_ptr<VibrationalSCF> vscf = std::make_shared<VSCF>(parms, basis);

    if (id == 0) {
      std::cout << "Number of SP Calculations needed to construct hybrid PES: " << pes->getNumSPsForVSCF() << std::endl;
      gettimeofday(&first, nullptr);
      double initializationTime = first.tv_sec - start.tv_sec + 1e-6 * (first.tv_usec - start.tv_usec);
      std::cout << "Initialization took " << initializationTime << " seconds." << std::endl << std::endl;

      double vscfEnergy = vscf->SCF();
      std::cout << std::endl << "Final VSCF energy is: " << vscfEnergy << " Hartree" << std::endl;
      std::cout << "Final VSCF energy is: " << vscfEnergy * Utils::Constants::invCentimeter_per_hartree << " cm^-1"
                << std::endl;
      // Final operations
      gettimeofday(&then, nullptr);
      double elapsed = then.tv_sec - start.tv_sec + 1e-6 * (then.tv_usec - start.tv_usec);
      std::cout << std::endl << "VSCF took " << elapsed << " seconds in total." << std::endl << std::endl;

      std::cout << "------ Second step: VCI routine ------" << std::endl << std::endl;

      parms.printParametersVCI();
      std::shared_ptr<VibrationalCI> vci = std::make_shared<VCI>(parms, vscf);
      double vciEnergy = vci->CI();

      std::cout << "Minimal VCI energy is: " << vciEnergy << " Hartree" << std::endl;
      std::cout << "Minimal VCI energy is: " << vciEnergy * Utils::Constants::invCentimeter_per_hartree << " cm^-1"
                << std::endl
                << std::endl;

      // Final operations
      gettimeofday(&after_vci, nullptr);
      elapsed = after_vci.tv_sec - then.tv_sec + 1e-6 * (after_vci.tv_usec - then.tv_usec);
      std::cout << "VCI took " << elapsed << " seconds in total." << std::endl;
    }
  }
  catch (std::exception& e) {
    std::cout << "Exception thrown!" << std::endl;
    std::cout << e.what() << std::endl;
  }
#ifdef MPI_PARALLEL
  MPI_Finalize();
#endif
}
