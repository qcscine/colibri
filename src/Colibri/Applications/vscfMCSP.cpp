/**
 * @file vscfMCSP.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#include <VSCF/VSCF.h>
#include <VSCF/VibrationalSCF.h>
#include <VibrationalUtils/PESLibrary/PES.h>
#include <VibrationalUtils/PESLibrary/PesMCSP.h>
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
 * Main routine for VSCF calculations
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
    std::cout << "------ VSCF routine ------" << std::endl << std::endl;
  }

  timeval now, first, then;
  gettimeofday(&now, nullptr);
  std::string inputFileName = std::string(argv[1]);
  if (id == 0) {
    std::cout << "Input file: " << inputFileName << std::endl << std::endl;
  }

  std::cout << std::setprecision(8);

  try {
    VibrationalParameters parms;
    parms.setParametersN(inputFileName);
    if (id == 0) {
      parms.printParametersN();
    }
    parms.setPesMCSP(parms.inputPESFname_);
    std::shared_ptr<PES> pes = std::make_shared<PesMCSP>();

    pes->setReferenceGeometry(parms.oneModeTerms_, parms.twoModeTerms_, parms.threeModeTerms_, parms.fourModeTerms_,
                              parms.fiveModeTerms_, parms.sixModeTerms_);
    std::shared_ptr<ModalBasis> basis;
    if (parms.primitiveBasisType_ == "DG") {
      basis = std::make_shared<DistributedGaussians>(parms, pes);
    } else if (parms.primitiveBasisType_ == "DVR") {
      basis = std::make_shared<DVR>(parms, pes);
    } else {
      if (id == 0) {
        std::cout << "The basis specifier is not know. Please set BasisType "
                     "either to DG or DVR. Abort";
      }
      exit(2);
    }
    std::shared_ptr<VibrationalSCF> vscf = std::make_shared<VSCF>(parms, basis);
    if (id == 0) {
      gettimeofday(&first, nullptr);
      double initializationTime = first.tv_sec - now.tv_sec + 1e-6 * (first.tv_usec - now.tv_usec);
      std::cout << "Initialization took " << initializationTime << " seconds." << std::endl;

      double vscfEnergy = vscf->SCF();
      std::cout << "Final VSCF Energy is: " << vscfEnergy << " Hartree" << std::endl;
      std::cout << "Final VSCF Energy is: " << vscfEnergy * 219474.63 << " cm^-1" << std::endl;
      // Final operations
      gettimeofday(&then, nullptr);
      double elapsed = then.tv_sec - now.tv_sec + 1e-6 * (then.tv_usec - now.tv_usec);
      std::cout << "VSCF took " << elapsed << " seconds in total." << std::endl;
    }
  }
  catch (std::exception& e) {
    std::cout << "Exception thrown!" << std::endl;
    std::cout << e.what() << std::endl;
    exit(1);
  }
#ifdef MPI_PARALLEL
  MPI_Finalize();
#endif
}
