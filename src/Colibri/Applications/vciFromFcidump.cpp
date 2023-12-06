/**
 * @file vciFromFcidump.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#include "Utils/Constants.h"
#include <VCI/VCI.h>
#include <VCI/VibrationalCI.h>
#include <VSCF/VSCF.h>
#include <VSCF/VibrationalSCF.h>
#include <VibrationalUtils/PESLibrary/OnTheFlyPes.h>
#include <VibrationalUtils/PESLibrary/PES.h>
#include <VibrationalUtils/Tensor.h>
#include <VibrationalUtils/VibrationalBases/DVR.h>
#include <VibrationalUtils/VibrationalBases/DistributedGaussians.h>
#include <VibrationalUtils/VibrationalParameters.h>
#include <sys/time.h>
#include <iostream>

#ifdef MPI_PARALLEL
  #include <mpi.h>
#endif

using namespace Scine;
using namespace Colibri;
/**************************************
 * Main routine for VCI calculation based on FCIDUMP from previous calculation
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
    std::cout << "------ VCI based on FCIDUMP ------" << std::endl << std::endl;
  }

  timeval start, after_vci;
  gettimeofday(&start, nullptr);
  std::string inputFileName = std::string(argv[1]);
  if (id == 0) {
    std::cout << "Input file: " << inputFileName << std::endl << std::endl;
  }

  std::cout << std::setprecision(8);

  try {
    VibrationalParameters parms;
    parms.setParametersOnlyVCI(inputFileName);
    if (id == 0) {
      parms.printParametersOnlyVCI();

      std::cout << std::endl << "------ Setting up VCI ------" << std::endl << std::endl;
      std::shared_ptr<VibrationalCI> vci = std::make_shared<VCI>(parms);
      double vciEnergy = vci->CI();

      std::cout << "Minimal VCI energy is: " << vciEnergy << " Hartree" << std::endl;
      std::cout << "Minimal VCI energy is: " << vciEnergy * Utils::Constants::invCentimeter_per_hartree << " cm^-1"
                << std::endl
                << std::endl;

      // Final operations
      gettimeofday(&after_vci, nullptr);
      double elapsed = after_vci.tv_sec - start.tv_sec + 1e-6 * (after_vci.tv_usec - start.tv_usec);
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
