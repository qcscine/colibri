/**
 * @file vscfMC.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#include <VSCF/VSCF.h>
#include <VSCF/VibrationalSCF.h>
#include <VibrationalUtils/PESLibrary/PES.h>
#include <VibrationalUtils/PESLibrary/PesTaylorParam.h>
#include <VibrationalUtils/Tensor.h>
#include <VibrationalUtils/VibrationalBases/DVR.h>
#include <VibrationalUtils/VibrationalBases/DistributedGaussians.h>
#include <VibrationalUtils/VibrationalParameters.h>
#include <sys/time.h>
#include <iostream>

using namespace Scine;
using namespace Colibri;

/**************************************
 * Main routine for VSCF calculations
 **************************************/

int main(int /*argc*/, char** argv) {
  std::cout << "------ VSCF routine ------" << std::endl << std::endl;

  timeval now;
  timeval then;
  gettimeofday(&now, nullptr);
  std::string inputFileName = std::string(argv[1]);
  std::cout << "Input file: " << inputFileName << std::endl << std::endl;

  try {
    VibrationalParameters parms;
    parms.setParametersN(inputFileName);
    parms.printParametersN();
    parms.setPesMc(parms.inputPESFname_);
    std::shared_ptr<PES> pes = std::make_shared<PesTayParam>();

    pes->setReferenceGeometry(parms.oneModeCoeffs_, parms.twoModeCoeffs_, parms.threeModeCoeffs_, parms.fourModeCoeffs_,
                              parms.fiveModeCoeffs_, parms.sixModeCoeffs_, parms.taylorExpOrder_);
    std::shared_ptr<ModalBasis> basis;
    if (parms.primitiveBasisType_ == "DG") {
      basis = std::make_shared<DistributedGaussians>(parms, pes);
    } else if (parms.primitiveBasisType_ == "DVR") {
      basis = std::make_shared<DVR>(parms, pes);
    } else {
      std::cout << "The basis specifier is not know. Please set BasisType "
                   "either to DG or DVR. Abort";
      exit(2);
    }
    std::shared_ptr<VibrationalSCF> vscf = std::make_shared<VSCF>(parms, basis);
    double vscfEnergy = vscf->SCF();
    std::cout << "Final VSCF Energy is: " << vscfEnergy << std::endl;
    // Final operations
    gettimeofday(&then, nullptr);
    double elapsed = then.tv_sec - now.tv_sec + 1e-6 * (then.tv_usec - now.tv_usec);
    std::cout << "Task took " << elapsed << " seconds." << std::endl;
  }
  catch (std::exception& e) {
    std::cout << "Exception thrown!" << std::endl;
    std::cout << e.what() << std::endl;
    exit(1);
  }
}
