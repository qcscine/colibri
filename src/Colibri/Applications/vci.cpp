/**
 * @file vci.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <VCI/VCI.h>
#include <VCI/VibrationalCI.h>
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

using namespace Scine;
using namespace Colibri;
/**************************************
 * Main routine for VCI calculation with prior VSCF
 **************************************/

int main(int /*argc*/, char** argv) {
  std::cout << "------ VCI with prior VSCF routine ------" << std::endl << std::endl;

  timeval now, first, then, before_vci, after_vci;
  gettimeofday(&now, nullptr);
  std::string inputFileName = std::string(argv[1]);
  std::cout << "Input file: " << inputFileName << std::endl << std::endl;

  std::cout << std::setprecision(8);

  try {
    std::cout << "------ First step: VSCF routine ------" << std::endl << std::endl;
    VibrationalParameters parms;
    parms.setParametersAlsoVCI(inputFileName);
    parms.printParametersVSCF();
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
      std::cout << "The basis specifier is not know. Please set BasisType "
                   "either to DG or DVR. Abort";
      exit(2);
    }
    std::shared_ptr<VibrationalSCF> vscf = std::make_shared<VSCF>(parms, basis);

    gettimeofday(&first, nullptr);
    double initializationTime = first.tv_sec - now.tv_sec + 1e-6 * (first.tv_usec - now.tv_usec);
    std::cout << "Initialization took " << initializationTime << " seconds." << std::endl;

    double vscfEnergy = vscf->SCF();
    std::cout << "Final VSCF energy is: " << vscfEnergy << " Hartree" << std::endl;
    std::cout << "Final VSCF energy is: " << vscfEnergy * 219474.63 << " cm^-1" << std::endl;
    // Final operations
    gettimeofday(&then, nullptr);
    double elapsed = then.tv_sec - now.tv_sec + 1e-6 * (then.tv_usec - now.tv_usec);
    std::cout << "VSCF took " << elapsed << " seconds in total." << std::endl << std::endl;

    std::cout << "------ Second step: VCI routine ------" << std::endl << std::endl;

    parms.printParametersVCI();
    std::shared_ptr<VibrationalCI> vci = std::make_shared<VCI>(parms, vscf);
    double vciEnergy = vci->CI();

    std::cout << "Minimal VCI energy is: " << vciEnergy << " Hartree" << std::endl;
    std::cout << "Minimal VCI energy is: " << vciEnergy * 219474.63 << " cm^-1" << std::endl << std::endl;

    // Final operations
    gettimeofday(&after_vci, nullptr);
    elapsed = after_vci.tv_sec - then.tv_sec + 1e-6 * (after_vci.tv_usec - then.tv_usec);
    std::cout << "VCI took " << elapsed << " seconds in total." << std::endl;
  }
  catch (std::exception& e) {
    std::cout << "Exception thrown!" << std::endl;
    std::cout << e.what() << std::endl;
    exit(1);
  }
}
