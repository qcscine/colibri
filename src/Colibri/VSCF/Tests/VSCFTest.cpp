/**
 * @file VSCFTest.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <VSCF/VSCF.h>
#include <VSCF/VibrationalSCF.h>
#include <VibrationalUtils/PESLibrary/NMAPesCalculator.h>
#include <VibrationalUtils/PESLibrary/PesCartesian.h>
#include <VibrationalUtils/PESLibrary/PesChristoffel.h>
#include <VibrationalUtils/PESLibrary/PesCoupledOscillator.h>
#include <VibrationalUtils/PESLibrary/PesQFF.h>
#include <VibrationalUtils/Tensor.h>
#include <VibrationalUtils/VibrationalBases/DVR.h>
#include <VibrationalUtils/VibrationalBases/DistributedGaussians.h>
#include <VibrationalUtils/VibrationalBases/ModeGrid.h>
#include <VibrationalUtils/VibrationalParameters.h>
#include <gmock/gmock.h>
#include <cmath>
#include <cstdlib>
#include <memory>

#ifdef MPI_PARALLEL
  #include <mpi.h>
#endif

namespace Scine::Colibri {

using namespace testing;

class VSCFTest : public Test {
 public:
  // Class members
  tensor<1> K1, k1Co2, k1Ch, k1c;
  tensor<2> K2, k2Co2, k2Ch, k2c;
  tensor<3> K3, k3Co2, k3Ch, k3ChD, k3c;
  tensor<4> K4, k4Co2, k4Ch, k4c;

  std::shared_ptr<NMAPesCalculator> nmaPesCalculator_;
  std::shared_ptr<PES> pesCartesianNma_;

  std::ofstream inputFileNMA, inputFileCO2, inputFileEthylene, inputFileChristoffel, inputFileCoupledOscillator;
  std::string inputFileNMAName, inputFileCO2Name, inputFileEthyleneName, inputFileChristoffelName,
      inputFileCoupledOscillatorName;

 protected:
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

  void SetUp() final {
    inputFileNMAName = "nma.ini";
    inputFileCO2Name = "co2.ini";
    inputFileEthyleneName = "ethylene.ini";
    inputFileChristoffelName = "christoffel.ini";
    inputFileCoupledOscillatorName = "coupled.ini";

    // Creates the PES calculators
    nmaPesCalculator_ = std::make_shared<NMAPesCalculator>();

    // Sets up the PES objects.
    // Force constants in Eh/a0^2u
    // CO2 Hirata
    k1Co2[{2}] = -3.79376328358405E-05;
    k2Co2[{0, 0}] = 9.62186495432437E-06;
    k2Co2[{1, 1}] = 9.62186495432437E-06;
    k2Co2[{2, 2}] = 3.74505381530253E-05;
    k2Co2[{3, 3}] = 0.000117382604756;
    k3Co2[{2, 2, 2}] = 5.83894731065161E-07;
    k4Co2[{0, 0, 0, 0}] = 2.8412590870909E-09;
    k4Co2[{1, 1, 1, 1}] = 2.8412590870909E-09;
    k4Co2[{2, 2, 2, 2}] = 7.53429226524512E-09;
    k4Co2[{3, 3, 3, 3}] = 8.82914181764278E-08;
    k3Co2[{2, 0, 0}] = -1.60960891738734E-07;
    k3Co2[{2, 1, 1}] = -1.60960891738734E-07;
    k3Co2[{2, 3, 3}] = 1.95710501210382E-06;
    k4Co2[{2, 2, 0, 0}] = -3.83475582767667E-10;
    k4Co2[{2, 2, 1, 1}] = -3.83475582767667E-10;
    k4Co2[{2, 2, 3, 3}] = 2.5417941550711E-08;
    k4Co2[{0, 0, 1, 1}] = 9.84482114578557E-10;
    k4Co2[{0, 0, 3, 3}] = -1.66369406677665E-08;
    k4Co2[{1, 1, 3, 3}] = -1.66369406677665E-08;
    // Christoffel
    double lambda = -0.10;
    double mu = -0.10;
    double eta = 0.10;
    double zeta = 0.10;
    k2Ch[{0, 0}] = 0.49;
    k2Ch[{1, 1}] = 1.00;
    k2Ch[{2, 2}] = 1.69;
    // orig
    k3Ch[{0, 0, 0}] = lambda * eta;
    k3Ch[{2, 2, 2}] = mu * zeta;
    k3Ch[{0, 2, 2}] = lambda;
    k3Ch[{2, 1, 1}] = mu;
    // mod
    double damping = 0.5;
    k3ChD[{0, 0, 0}] = lambda * eta;
    k3ChD[{2, 2, 2}] = mu * zeta;
    k3ChD[{0, 2, 2}] = lambda * damping;
    k3ChD[{2, 1, 1}] = mu * damping;

    // ETHYLENE
    // Quadratic
    K2[{0, 0}] = 0.00021566870545;
    K2[{1, 1}] = 0.000211878018278;
    K2[{2, 2}] = 0.000204598362885;
    K2[{3, 3}] = 0.000202804506582;
    K2[{4, 4}] = 5.9488005339e-05;
    K2[{5, 5}] = 4.54004732742e-05;
    K2[{6, 6}] = 3.96568416263e-05;
    K2[{7, 7}] = 3.2256498653e-05;
    K2[{8, 8}] = 2.36437940807e-05;
    K2[{9, 9}] = 2.008899627e-05;
    K2[{10, 10}] = 1.99024791009e-05;
    K2[{11, 11}] = 1.45099385948e-05;
    // CUBIC
    K3[{1, 1, 1}] = 9.12259082591e-09;
    K3[{2, 2, 2}] = 7.24179526943e-06;
    K3[{4, 4, 4}] = 1.033936456e-06;
    K3[{6, 6, 6}] = -1.49302120278e-07;
    K3[{1, 0, 0}] = 6.55284693129e-09;
    K3[{2, 0, 0}] = 7.89733693695e-06;
    K3[{2, 1, 1}] = 7.94577660937e-06;
    K3[{2, 2, 1}] = -2.69823108935e-09;
    K3[{3, 1, 0}] = -7.88114755042e-06;
    K3[{3, 2, 0}] = 2.56974389462e-09;
    K3[{3, 3, 2}] = 7.26826363155e-06;
    K3[{4, 0, 0}] = -3.94455687825e-08;
    K3[{4, 1, 1}] = -7.65783680597e-08;
    K3[{4, 2, 2}] = 4.7771539001e-07;
    K3[{4, 3, 3}] = 3.89059225646e-07;
    K3[{4, 4, 2}] = -5.97850917084e-07;
    K3[{5, 1, 0}] = -5.04569213709e-07;
    K3[{5, 2, 0}] = 6.42435973656e-10;
    K3[{5, 3, 2}] = -7.3109213802e-08;
    K3[{5, 4, 3}] = 4.81313031463e-07;
    K3[{5, 5, 2}] = -5.63416348896e-07;
    K3[{5, 5, 4}] = 1.86948868334e-07;
    K3[{6, 0, 0}] = -4.7527413331e-07;
    K3[{6, 1, 1}] = -4.29532691986e-07;
    K3[{6, 2, 2}] = -1.02404294201e-07;
    K3[{6, 3, 3}] = -6.84836747917e-08;
    K3[{6, 4, 2}] = -1.60865967803e-07;
    K3[{6, 4, 4}] = -4.61012054695e-07;
    K3[{6, 5, 3}] = 3.71456479968e-07;
    K3[{6, 5, 5}] = 1.31056938626e-07;
    K3[{6, 6, 2}] = -2.96805419829e-07;
    K3[{6, 6, 4}] = 3.3432368069e-07;
    K3[{7, 1, 1}] = -5.13948778924e-10;
    K3[{7, 2, 1}] = -2.0532253718e-07;
    K3[{7, 3, 0}] = 1.46732376383e-07;
    K3[{7, 4, 1}] = 6.13012406062e-07;
    K3[{7, 4, 2}] = -5.13948778924e-10;
    K3[{7, 5, 0}] = 7.15416700263e-07;
    K3[{7, 5, 3}] = 5.13948778924e-10;
    K3[{7, 6, 1}] = 4.23236819444e-07;
    K3[{7, 7, 2}] = -5.35920089223e-07;
    K3[{7, 7, 4}] = -1.11269910637e-07;
    K3[{7, 7, 6}] = -1.58810172688e-07;
    K3[{8, 1, 0}] = -7.70923168387e-10;
    K3[{8, 3, 2}] = 7.70923168387e-10;
    K3[{8, 8, 1}] = -1.92730792097e-09;
    K3[{8, 8, 2}] = -8.34010381e-07;
    K3[{8, 8, 4}] = 3.46915425774e-08;
    K3[{8, 8, 6}] = -6.28302382235e-08;
    K3[{9, 2, 1}] = 7.70923168387e-10;
    K3[{9, 3, 0}] = -7.70923168387e-10;
    K3[{9, 8, 0}] = 9.57615062331e-07;
    K3[{9, 8, 3}] = -6.42435973656e-10;
    K3[{9, 9, 2}] = -8.98639439949e-07;
    K3[{9, 9, 4}] = 8.45445741331e-08;
    K3[{9, 9, 6}] = 3.84176712246e-08;
    K3[{10, 2, 0}] = -8.99410363118e-10;
    K3[{10, 3, 1}] = 8.99410363118e-10;
    K3[{10, 8, 1}] = -9.68407986688e-07;
    K3[{10, 8, 7}] = 4.33001846244e-08;
    K3[{10, 9, 0}] = -5.13948778924e-10;
    K3[{10, 9, 3}] = -9.74189910451e-07;
    K3[{10, 9, 5}] = -1.53156736119e-07;
    K3[{10, 9, 8}] = 5.13948778924e-10;
    K3[{10, 10, 1}] = 1.15638475258e-09;
    K3[{10, 10, 2}] = -1.01055178656e-06;
    K3[{10, 10, 4}] = 1.03689166148e-07;
    K3[{10, 10, 6}] = 1.20906450242e-07;
    K3[{11, 1, 0}] = -5.13948778924e-10;
    K3[{11, 2, 0}] = -1.09214115521e-07;
    K3[{11, 3, 1}] = 1.09599577106e-07;
    K3[{11, 4, 0}] = 6.998697497e-07;
    K3[{11, 5, 1}] = 8.97997003976e-07;
    K3[{11, 5, 2}] = -7.70923168387e-10;
    K3[{11, 6, 0}] = 6.15839124346e-07;
    K3[{11, 7, 3}] = 7.99318838422e-07;
    K3[{11, 7, 5}] = -1.11398397832e-07;
    K3[{11, 9, 8}] = -2.44125669989e-08;
    K3[{11, 11, 2}] = -1.0169761463e-06;
    K3[{11, 11, 4}] = -8.12039070701e-08;
    K3[{11, 11, 6}] = -1.25018040473e-07;
    K4[{0, 0, 0, 0}] = 2.83148348759e-07;
    K4[{1, 1, 1, 1}] = 2.86263076896e-07;
    K4[{2, 2, 2, 2}] = 2.37008240741e-07;
    K4[{3, 3, 3, 3}] = 2.39108801359e-07;
    K4[{4, 4, 4, 4}] = 2.1788049964e-08;
    K4[{6, 6, 6, 6}] = 4.61942772027e-09;
    K4[{7, 7, 7, 7}] = 6.27459726173e-09;
    K4[{8, 8, 8, 8}] = 1.07977242087e-08;
    K4[{9, 9, 9, 9}] = 1.66269303938e-08;
    K4[{10, 10, 10, 10}] = 2.00395890484e-08;
    K4[{11, 11, 11, 11}] = 2.1986670309e-08;
    K4[{1, 0, 0, 0}] = 2.10657941641e-11;
    K4[{1, 1, 0, 0}] = 2.84303958038e-07;
    K4[{1, 1, 1, 0}] = -2.70845924967e-11;
    K4[{2, 1, 0, 0}] = -4.84513265774e-10;
    K4[{2, 1, 1, 1}] = -4.03259488284e-10;
    K4[{2, 2, 0, 0}] = 2.62594152453e-07;
    K4[{2, 2, 1, 0}] = -4.21315883282e-11;
    K4[{2, 2, 1, 1}] = 2.6281985739e-07;
    K4[{2, 2, 2, 1}] = -7.52349791574e-10;
    K4[{3, 0, 0, 0}] = 4.81503866607e-10;
    K4[{3, 1, 1, 0}] = 4.03259488284e-10;
    K4[{3, 2, 0, 0}] = -3.61127899956e-11;
    K4[{3, 2, 1, 1}] = 4.21315883282e-11;
    K4[{3, 2, 2, 0}] = 7.52349791574e-10;
    K4[{3, 2, 2, 2}] = 6.62067816585e-11;
    K4[{3, 3, 0, 0}] = 2.6328932366e-07;
    K4[{3, 3, 1, 0}] = 3.31033908293e-11;
    K4[{3, 3, 1, 1}] = 2.64020607658e-07;
    K4[{3, 3, 2, 1}] = -6.71096014084e-10;
    K4[{3, 3, 2, 2}] = 2.37586045381e-07;
    K4[{3, 3, 3, 0}] = 6.71096014084e-10;
    K4[{3, 3, 3, 2}] = -4.81503866607e-11;
    K4[{4, 1, 0, 0}] = -7.52349791574e-11;
    K4[{4, 1, 1, 1}] = -7.52349791574e-11;
    K4[{4, 2, 0, 0}] = 5.74494300846e-09;
    K4[{4, 2, 1, 1}] = 6.44011421587e-09;
    K4[{4, 2, 2, 1}] = -9.32913741552e-11;
    K4[{4, 2, 2, 2}] = 1.66058645996e-08;
    K4[{4, 3, 2, 2}] = 1.50469958315e-11;
    K4[{4, 3, 3, 1}] = -9.32913741552e-11;
    K4[{4, 3, 3, 2}] = 1.47972157007e-08;
    K4[{4, 4, 0, 0}] = -3.38467124233e-08;
    K4[{4, 4, 1, 0}] = -2.10657941641e-11;
    K4[{4, 4, 1, 1}] = -3.44696580508e-08;
    K4[{4, 4, 2, 1}] = -2.97930517463e-10;
    K4[{4, 4, 2, 2}] = -2.34823416946e-08;
    K4[{4, 4, 3, 0}] = 3.03949315796e-10;
    K4[{4, 4, 3, 2}] = 3.31033908293e-11;
    K4[{4, 4, 3, 3}] = -2.57634662627e-08;
    K4[{4, 4, 4, 1}] = -1.08338369987e-10;
    K4[{4, 4, 4, 2}] = -5.71785841596e-09;
    K4[{5, 0, 0, 0}] = -7.22255799911e-11;
    K4[{5, 1, 1, 0}] = -7.52349791574e-11;
    K4[{5, 2, 0, 0}] = 2.40751933304e-11;
    K4[{5, 2, 1, 1}] = -3.0093991663e-11;
    K4[{5, 2, 2, 0}] = -8.125377749e-11;
    K4[{5, 2, 2, 2}] = -4.21315883282e-11;
    K4[{5, 3, 0, 0}] = 1.1146814512e-08;
    K4[{5, 3, 1, 1}] = 1.13905758444e-08;
    K4[{5, 3, 2, 2}] = -6.41002022421e-09;
    K4[{5, 3, 3, 0}] = -8.42631766563e-11;
    K4[{5, 3, 3, 2}] = 2.70845924967e-11;
    K4[{5, 3, 3, 3}] = -3.64438239039e-09;
    K4[{5, 4, 0, 0}] = 1.50469958315e-11;
    K4[{5, 4, 1, 1}] = -1.50469958315e-11;
    K4[{5, 4, 2, 2}] = -2.10657941641e-11;
    K4[{5, 4, 3, 3}] = 1.50469958315e-11;
    K4[{5, 4, 4, 0}] = -1.23385365818e-10;
    K4[{5, 4, 4, 2}] = -2.10657941641e-11;
    K4[{5, 4, 4, 3}] = -4.33353479947e-10;
    K4[{5, 4, 4, 4}] = -1.20375966652e-11;
    K4[{5, 5, 0, 0}] = -6.09583895125e-08;
    K4[{5, 5, 1, 0}] = 2.70845924967e-11;
    K4[{5, 5, 1, 1}] = -6.13435926058e-08;
    K4[{5, 5, 2, 1}] = -2.61817727468e-10;
    K4[{5, 5, 2, 2}] = -4.77681929666e-08;
    K4[{5, 5, 3, 0}] = 2.678365258e-10;
    K4[{5, 5, 3, 2}] = -6.01879833259e-11;
    K4[{5, 5, 3, 3}] = -4.73258112892e-08;
    K4[{5, 5, 4, 1}] = -9.32913741552e-11;
    K4[{5, 5, 4, 2}] = 3.40363045708e-09;
    K4[{5, 5, 4, 3}] = -1.80563949978e-11;
    K4[{5, 5, 4, 4}] = 3.82494634036e-09;
    K4[{5, 5, 5, 0}] = -1.02319571654e-10;
    K4[{5, 5, 5, 2}] = 4.81503866607e-11;
    K4[{5, 5, 5, 3}] = -6.31973824922e-09;
    K4[{5, 5, 5, 4}] = 3.31033908293e-11;
    K4[{6, 1, 0, 0}] = -4.21315883282e-11;
    K4[{6, 1, 1, 1}] = -4.51409874944e-11;
    K4[{6, 2, 0, 0}] = -1.15591021977e-08;
    K4[{6, 2, 1, 1}] = -1.12822374744e-08;
    K4[{6, 2, 2, 1}] = -4.21315883282e-11;
    K4[{6, 2, 2, 2}] = 1.64915074313e-09;
    K4[{6, 3, 0, 0}] = -1.50469958315e-11;
    K4[{6, 3, 1, 1}] = 1.80563949978e-11;
    K4[{6, 3, 2, 2}] = 2.70845924967e-11;
    K4[{6, 3, 3, 1}] = -4.51409874944e-11;
    K4[{6, 3, 3, 2}] = -5.38682450767e-10;
    K4[{6, 3, 3, 3}] = -1.80563949978e-11;
    K4[{6, 4, 0, 0}] = -3.11051497828e-08;
    K4[{6, 4, 1, 1}] = -3.04340537688e-08;
    K4[{6, 4, 2, 2}] = -2.50923702486e-08;
    K4[{6, 4, 3, 3}] = -2.40601463345e-08;
    K4[{6, 4, 4, 1}] = -7.82443783237e-11;
    K4[{6, 4, 4, 2}] = 4.47798595945e-09;
    K4[{6, 4, 4, 3}] = 1.50469958315e-11;
    K4[{6, 4, 4, 4}] = -9.69929351297e-09;
    K4[{6, 5, 0, 0}] = 3.0093991663e-11;
    K4[{6, 5, 1, 1}] = -3.61127899956e-11;
    K4[{6, 5, 2, 2}] = -4.81503866607e-11;
    K4[{6, 5, 3, 3}] = 3.0093991663e-11;
    K4[{6, 5, 4, 4}] = -2.40751933304e-11;
    K4[{6, 5, 5, 1}] = -6.31973824922e-11;
    K4[{6, 5, 5, 2}] = 4.07773587033e-09;
    K4[{6, 5, 5, 3}] = -3.0093991663e-11;
    K4[{6, 5, 5, 4}] = 3.13278453211e-09;
    K4[{6, 5, 5, 5}] = 6.01879833259e-11;
    K4[{6, 6, 0, 0}] = -2.93265948756e-08;
    K4[{6, 6, 1, 1}] = -2.95583186114e-08;
    K4[{6, 6, 2, 1}] = -3.37052706625e-10;
    K4[{6, 6, 2, 2}] = -2.27089261089e-08;
    K4[{6, 6, 3, 0}] = 3.40062105792e-10;
    K4[{6, 6, 3, 3}] = -2.37050372329e-08;
    K4[{6, 6, 4, 1}] = -1.14357168319e-10;
    K4[{6, 6, 4, 2}] = 4.90532064106e-10;
    K4[{6, 6, 4, 4}] = 7.13227602412e-09;
    K4[{6, 6, 5, 0}] = -1.32413563317e-10;
    K4[{6, 6, 5, 3}] = -3.44275264624e-09;
    K4[{6, 6, 5, 5}] = 1.70933872646e-09;
    K4[{6, 6, 6, 1}] = -8.42631766563e-11;
    K4[{6, 6, 6, 2}] = 3.30732968376e-09;
    K4[{6, 6, 6, 4}] = -1.27598524651e-09;
    K4[{7, 0, 0, 0}] = 1.50469958315e-11;
    K4[{7, 1, 0, 0}] = 6.8132797125e-09;
    K4[{7, 1, 1, 0}] = -1.80563949978e-11;
    K4[{7, 1, 1, 1}] = 6.99685306164e-09;
    K4[{7, 2, 0, 0}] = 7.22255799911e-11;
    K4[{7, 2, 1, 1}] = 7.22255799911e-11;
    K4[{7, 2, 2, 0}] = -2.70845924967e-11;
    K4[{7, 2, 2, 1}] = -9.84374467296e-09;
    K4[{7, 2, 2, 2}] = 1.0532897082e-10;
    K4[{7, 3, 3, 0}] = 1.80563949978e-11;
    K4[{7, 3, 3, 1}] = -7.12625722579e-09;
    K4[{7, 3, 3, 2}] = 1.02319571654e-10;
    K4[{7, 4, 2, 2}] = -2.40751933304e-11;
    K4[{7, 4, 3, 3}] = -1.50469958315e-11;
    K4[{7, 4, 4, 0}] = -1.50469958315e-11;
    K4[{7, 4, 4, 1}] = 2.09755121891e-09;
    K4[{7, 4, 4, 2}] = 1.0532897082e-10;
    K4[{7, 4, 4, 4}] = 3.61127899956e-11;
    K4[{7, 5, 5, 0}] = 3.0093991663e-11;
    K4[{7, 5, 5, 1}] = -2.34131255138e-09;
    K4[{7, 5, 5, 2}] = 8.72725758226e-11;
    K4[{7, 5, 5, 4}] = 5.1159785827e-11;
    K4[{7, 6, 2, 2}] = -2.70845924967e-11;
    K4[{7, 6, 3, 3}] = -1.80563949978e-11;
    K4[{7, 6, 4, 4}] = 3.31033908293e-11;
    K4[{7, 6, 5, 5}] = 5.1159785827e-11;
    K4[{7, 6, 6, 1}] = -8.125377749e-10;
    K4[{7, 6, 6, 2}] = 1.14357168319e-10;
    K4[{7, 6, 6, 4}] = 4.81503866607e-11;
    K4[{7, 6, 6, 6}] = 4.21315883282e-11;
    K4[{7, 7, 0, 0}] = -4.29772294939e-08;
    K4[{7, 7, 1, 1}] = -4.25438760139e-08;
    K4[{7, 7, 2, 1}] = -2.76864723299e-10;
    K4[{7, 7, 2, 2}] = -3.9414100881e-08;
    K4[{7, 7, 3, 0}] = 2.79874122466e-10;
    K4[{7, 7, 3, 3}] = -3.93840068893e-08;
    K4[{7, 7, 4, 1}] = -9.63007733215e-11;
    K4[{7, 7, 4, 2}] = -2.63623366968e-09;
    K4[{7, 7, 4, 4}] = 1.65516954146e-09;
    K4[{7, 7, 5, 0}] = -1.11347769153e-10;
    K4[{7, 7, 5, 3}] = 4.25529042114e-09;
    K4[{7, 7, 5, 5}] = 7.36099036076e-09;
    K4[{7, 7, 6, 1}] = -7.22255799911e-11;
    K4[{7, 7, 6, 2}] = -4.49303295528e-09;
    K4[{7, 7, 6, 4}] = 4.22519642948e-09;
    K4[{7, 7, 6, 6}] = 2.97629577547e-09;
    K4[{7, 7, 7, 1}] = 5.85930017678e-09;
    K4[{7, 7, 7, 2}] = 9.02819749889e-11;
    K4[{7, 7, 7, 4}] = 4.81503866607e-11;
    K4[{7, 7, 7, 6}] = 4.51409874944e-11;
    K4[{8, 0, 0, 0}] = 1.56488756647e-10;
    K4[{8, 1, 1, 0}] = 1.56488756647e-10;
    K4[{8, 2, 2, 0}] = 1.95610945809e-10;
    K4[{8, 3, 0, 0}] = -2.40751933304e-10;
    K4[{8, 3, 1, 1}] = -2.37742534137e-10;
    K4[{8, 3, 2, 2}] = -2.34733134971e-10;
    K4[{8, 3, 3, 0}] = 1.89592147477e-10;
    K4[{8, 3, 3, 3}] = -2.37742534137e-10;
    K4[{8, 4, 4, 0}] = 4.81503866607e-11;
    K4[{8, 4, 4, 3}] = 3.31033908293e-11;
    K4[{8, 5, 0, 0}] = 7.82443783237e-11;
    K4[{8, 5, 1, 1}] = 7.82443783237e-11;
    K4[{8, 5, 2, 2}] = 6.31973824922e-11;
    K4[{8, 5, 3, 3}] = 6.31973824922e-11;
    K4[{8, 5, 5, 0}] = 5.71785841596e-11;
    K4[{8, 5, 5, 3}] = 2.40751933304e-11;
    K4[{8, 6, 6, 0}] = 6.31973824922e-11;
    K4[{8, 6, 6, 3}] = 4.21315883282e-11;
    K4[{8, 6, 6, 5}] = -1.20375966652e-11;
    K4[{8, 7, 7, 0}] = 4.81503866607e-11;
    K4[{8, 7, 7, 3}] = 2.70845924967e-11;
    K4[{8, 8, 0, 0}] = -5.43136361533e-08;
    K4[{8, 8, 1, 0}] = 2.70845924967e-11;
    K4[{8, 8, 1, 1}] = -5.37298127151e-08;
    K4[{8, 8, 2, 1}] = 2.40751933304e-11;
    K4[{8, 8, 2, 2}] = -4.98446783914e-08;
    K4[{8, 8, 3, 0}] = -2.10657941641e-11;
    K4[{8, 8, 3, 2}] = -4.51409874944e-11;
    K4[{8, 8, 3, 3}] = -5.21709439469e-08;
    K4[{8, 8, 4, 1}] = -1.20375966652e-11;
    K4[{8, 8, 4, 2}] = -6.56049018253e-10;
    K4[{8, 8, 4, 3}] = -1.20375966652e-11;
    K4[{8, 8, 4, 4}] = 1.38733301566e-09;
    K4[{8, 8, 5, 0}] = -2.40751933304e-11;
    K4[{8, 8, 5, 2}] = 3.0093991663e-11;
    K4[{8, 8, 5, 3}] = 1.03824271237e-09;
    K4[{8, 8, 5, 4}] = 1.50469958315e-11;
    K4[{8, 8, 5, 5}] = 1.34520142733e-09;
    K4[{8, 8, 6, 1}] = -2.10657941641e-11;
    K4[{8, 8, 6, 2}] = -7.61377989073e-10;
    K4[{8, 8, 6, 3}] = -1.80563949978e-11;
    K4[{8, 8, 6, 4}] = 8.06518976567e-10;
    K4[{8, 8, 6, 5}] = 3.61127899956e-11;
    K4[{8, 8, 6, 6}] = 1.35422962483e-10;
    K4[{8, 8, 7, 0}] = 2.10657941641e-11;
    K4[{8, 8, 7, 1}] = 1.79962070145e-09;
    K4[{8, 8, 7, 2}] = 2.40751933304e-11;
    K4[{8, 8, 7, 4}] = 2.10657941641e-11;
    K4[{8, 8, 7, 6}] = 2.10657941641e-11;
    K4[{8, 8, 7, 7}] = 2.56400808968e-09;
    K4[{8, 8, 8, 0}] = 1.26394764984e-10;
    K4[{8, 8, 8, 3}] = -1.53479357481e-10;
    K4[{8, 8, 8, 5}] = 4.81503866607e-11;
    K4[{9, 1, 0, 0}] = -2.49780130803e-10;
    K4[{9, 1, 1, 1}] = -2.49780130803e-10;
    K4[{9, 2, 0, 0}] = -1.80563949978e-10;
    K4[{9, 2, 1, 1}] = -1.80563949978e-10;
    K4[{9, 2, 2, 1}] = -2.4376133247e-10;
    K4[{9, 2, 2, 2}] = -2.16676739973e-10;
    K4[{9, 3, 3, 1}] = -2.49780130803e-10;
    K4[{9, 3, 3, 2}] = -2.10657941641e-10;
    K4[{9, 4, 0, 0}] = 7.82443783237e-11;
    K4[{9, 4, 1, 1}] = 7.82443783237e-11;
    K4[{9, 4, 2, 2}] = 6.62067816585e-11;
    K4[{9, 4, 3, 3}] = 6.62067816585e-11;
    K4[{9, 4, 4, 1}] = 3.61127899956e-11;
    K4[{9, 4, 4, 2}] = -4.81503866607e-11;
    K4[{9, 5, 5, 1}] = 3.0093991663e-11;
    K4[{9, 5, 5, 2}] = -6.01879833259e-11;
    K4[{9, 6, 0, 0}] = 9.32913741552e-11;
    K4[{9, 6, 1, 1}] = 9.32913741552e-11;
    K4[{9, 6, 2, 2}] = 8.42631766563e-11;
    K4[{9, 6, 3, 3}] = 8.125377749e-11;
    K4[{9, 6, 6, 1}] = 4.81503866607e-11;
    K4[{9, 6, 6, 2}] = -6.01879833259e-11;
    K4[{9, 7, 0, 0}] = 7.82443783237e-11;
    K4[{9, 7, 1, 1}] = 7.82443783237e-11;
    K4[{9, 7, 2, 2}] = 6.01879833259e-11;
    K4[{9, 7, 3, 3}] = 6.31973824922e-11;
    K4[{9, 7, 6, 6}] = -1.20375966652e-11;
    K4[{9, 7, 7, 1}] = 3.0093991663e-11;
    K4[{9, 7, 7, 2}] = -5.1159785827e-11;
    K4[{9, 8, 8, 1}] = -1.68526353313e-10;
    K4[{9, 8, 8, 2}] = -1.3843236165e-10;
    K4[{9, 8, 8, 6}] = 2.10657941641e-11;
    K4[{9, 8, 8, 7}] = 4.51409874944e-11;
    K4[{9, 9, 0, 0}] = -5.98960716068e-08;
    K4[{9, 9, 1, 0}] = 1.50469958315e-11;
    K4[{9, 9, 1, 1}] = -6.14338745808e-08;
    K4[{9, 9, 2, 1}] = 9.63007733215e-11;
    K4[{9, 9, 2, 2}] = -5.27968989735e-08;
    K4[{9, 9, 3, 0}] = -9.32913741552e-11;
    K4[{9, 9, 3, 2}] = -2.10657941641e-11;
    K4[{9, 9, 3, 3}] = -5.2878152751e-08;
    K4[{9, 9, 4, 1}] = -1.80563949978e-11;
    K4[{9, 9, 4, 2}] = 6.19936228257e-10;
    K4[{9, 9, 4, 4}] = -5.1159785827e-11;
    K4[{9, 9, 5, 0}] = -3.31033908293e-11;
    K4[{9, 9, 5, 3}] = -2.34131255138e-09;
    K4[{9, 9, 5, 5}] = 2.85291040965e-09;
    K4[{9, 9, 6, 1}] = -2.70845924967e-11;
    K4[{9, 9, 6, 2}] = 1.23385365818e-09;
    K4[{9, 9, 6, 4}] = 2.2058895889e-09;
    K4[{9, 9, 6, 5}] = 1.20375966652e-11;
    K4[{9, 9, 6, 6}] = 4.15297084949e-10;
    K4[{9, 9, 7, 1}] = 1.9019402731e-09;
    K4[{9, 9, 7, 2}] = 2.10657941641e-11;
    K4[{9, 9, 7, 4}] = 2.70845924967e-11;
    K4[{9, 9, 7, 6}] = 2.40751933304e-11;
    K4[{9, 9, 7, 7}] = 2.6873934555e-09;
    K4[{9, 9, 8, 0}] = 1.26394764984e-10;
    K4[{9, 9, 8, 3}] = -1.56488756647e-10;
    K4[{9, 9, 8, 5}] = 5.1159785827e-11;
    K4[{9, 9, 8, 8}] = 1.26394764984e-08;
    K4[{9, 9, 9, 1}] = -1.68526353313e-10;
    K4[{9, 9, 9, 2}] = -1.3843236165e-10;
    K4[{9, 9, 9, 4}] = 1.50469958315e-11;
    K4[{9, 9, 9, 6}] = 2.40751933304e-11;
    K4[{9, 9, 9, 7}] = 4.51409874944e-11;
    K4[{10, 0, 0, 0}] = 2.61817727468e-10;
    K4[{10, 1, 1, 0}] = 2.61817727468e-10;
    K4[{10, 2, 2, 0}] = 2.49780130803e-10;
    K4[{10, 3, 0, 0}] = -1.83573349144e-10;
    K4[{10, 3, 1, 1}] = -1.83573349144e-10;
    K4[{10, 3, 2, 2}] = -2.1968613914e-10;
    K4[{10, 3, 3, 0}] = 2.55798929135e-10;
    K4[{10, 3, 3, 3}] = -2.13667340807e-10;
    K4[{10, 4, 4, 0}] = -3.61127899956e-11;
    K4[{10, 4, 4, 3}] = -5.1159785827e-11;
    K4[{10, 5, 0, 0}] = -1.29404164151e-10;
    K4[{10, 5, 1, 1}] = -1.32413563317e-10;
    K4[{10, 5, 2, 2}] = -1.14357168319e-10;
    K4[{10, 5, 3, 3}] = -1.17366567486e-10;
    K4[{10, 5, 5, 0}] = -3.0093991663e-11;
    K4[{10, 5, 5, 3}] = -6.31973824922e-11;
    K4[{10, 5, 5, 5}] = -1.50469958315e-11;
    K4[{10, 6, 6, 0}] = -5.1159785827e-11;
    K4[{10, 6, 6, 3}] = -6.31973824922e-11;
    K4[{10, 7, 7, 0}] = -3.0093991663e-11;
    K4[{10, 7, 7, 3}] = -5.1159785827e-11;
    K4[{10, 8, 0, 0}] = -1.98620344976e-10;
    K4[{10, 8, 1, 1}] = -2.13667340807e-10;
    K4[{10, 8, 2, 2}] = -1.08338369987e-10;
    K4[{10, 8, 3, 3}] = -1.23385365818e-10;
    K4[{10, 8, 4, 4}] = -4.21315883282e-11;
    K4[{10, 8, 5, 5}] = -6.01879833259e-11;
    K4[{10, 8, 6, 6}] = -3.61127899956e-11;
    K4[{10, 8, 7, 7}] = -4.81503866607e-11;
    K4[{10, 8, 8, 0}] = 1.74545151645e-10;
    K4[{10, 8, 8, 3}] = -1.3843236165e-10;
    K4[{10, 8, 8, 5}] = -3.0093991663e-11;
    K4[{10, 8, 8, 8}] = -3.31033908293e-11;
    K4[{10, 9, 0, 0}] = -2.40751933304e-11;
    K4[{10, 9, 1, 1}] = 3.0093991663e-11;
    K4[{10, 9, 2, 2}] = 3.91221891619e-11;
    K4[{10, 9, 3, 3}] = -2.70845924967e-11;
    K4[{10, 9, 4, 4}] = 2.10657941641e-11;
    K4[{10, 9, 5, 5}] = -5.41691849933e-11;
    K4[{10, 9, 8, 8}] = -3.0093991663e-11;
    K4[{10, 9, 9, 0}] = 1.74545151645e-10;
    K4[{10, 9, 9, 3}] = -1.3843236165e-10;
    K4[{10, 9, 9, 5}] = -3.61127899956e-11;
    K4[{10, 9, 9, 8}] = -5.71785841596e-11;
    K4[{10, 10, 0, 0}] = -6.6233866251e-08;
    K4[{10, 10, 1, 1}] = -6.46930538779e-08;
    K4[{10, 10, 2, 1}] = 1.68526353313e-10;
    K4[{10, 10, 2, 2}] = -5.70160766047e-08;
    K4[{10, 10, 3, 0}] = -1.68526353313e-10;
    K4[{10, 10, 3, 3}] = -5.71966405546e-08;
    K4[{10, 10, 4, 1}] = -2.40751933304e-11;
    K4[{10, 10, 4, 2}] = 6.47020820754e-10;
    K4[{10, 10, 4, 4}] = 1.02018631737e-09;
    K4[{10, 10, 5, 0}] = -3.91221891619e-11;
    K4[{10, 10, 5, 3}] = -3.96939750035e-09;
    K4[{10, 10, 5, 5}] = 2.26306817306e-09;
    K4[{10, 10, 6, 1}] = -3.0093991663e-11;
    K4[{10, 10, 6, 2}] = 3.56914741123e-09;
    K4[{10, 10, 6, 4}] = 1.20375966652e-09;
    K4[{10, 10, 6, 6}] = 1.23686305735e-09;
    K4[{10, 10, 7, 1}] = 7.46330993242e-10;
    K4[{10, 10, 7, 2}] = 2.40751933304e-11;
    K4[{10, 10, 7, 4}] = 3.31033908293e-11;
    K4[{10, 10, 7, 6}] = 3.0093991663e-11;
    K4[{10, 10, 7, 7}] = 2.83184461549e-09;
    K4[{10, 10, 8, 0}] = 1.29404164151e-10;
    K4[{10, 10, 8, 3}] = -1.59498155814e-10;
    K4[{10, 10, 8, 5}] = 5.1159785827e-11;
    K4[{10, 10, 8, 8}] = 1.28651814359e-08;
    K4[{10, 10, 9, 1}] = -1.71535752479e-10;
    K4[{10, 10, 9, 2}] = -1.41441760816e-10;
    K4[{10, 10, 9, 4}] = 1.20375966652e-11;
    K4[{10, 10, 9, 6}] = 2.10657941641e-11;
    K4[{10, 10, 9, 7}] = 4.51409874944e-11;
    K4[{10, 10, 9, 9}] = 1.67924473479e-08;
    K4[{10, 10, 10, 0}] = 1.77554550812e-10;
    K4[{10, 10, 10, 3}] = -1.44451159982e-10;
    K4[{10, 10, 10, 5}] = -3.31033908293e-11;
    K4[{10, 10, 10, 8}] = -8.72725758226e-11;
    K4[{11, 0, 0, 0}] = 1.14116416386e-08;
    K4[{11, 1, 1, 0}] = 1.05750286704e-08;
    K4[{11, 2, 2, 0}] = -8.72123878393e-09;
    K4[{11, 3, 0, 0}] = -8.42631766563e-11;
    K4[{11, 3, 1, 1}] = -8.125377749e-11;
    K4[{11, 3, 2, 2}] = -1.14357168319e-10;
    K4[{11, 3, 3, 0}] = -6.44914241337e-09;
    K4[{11, 3, 3, 3}] = -1.11347769153e-10;
    K4[{11, 4, 4, 0}] = 5.95861034927e-10;
    K4[{11, 4, 4, 3}] = -1.20375966652e-10;
    K4[{11, 5, 2, 2}] = -4.51409874944e-11;
    K4[{11, 5, 3, 3}] = -3.0093991663e-11;
    K4[{11, 5, 4, 4}] = 5.41691849933e-11;
    K4[{11, 5, 5, 0}] = -4.94745222939e-09;
    K4[{11, 5, 5, 3}] = -1.02319571654e-10;
    K4[{11, 5, 5, 5}] = 7.82443783237e-11;
    K4[{11, 6, 6, 0}] = -3.16889732211e-09;
    K4[{11, 6, 6, 3}] = -1.29404164151e-10;
    K4[{11, 6, 6, 5}] = 6.92161808248e-11;
    K4[{11, 7, 7, 0}] = 4.79999167024e-09;
    K4[{11, 7, 7, 3}] = -1.0532897082e-10;
    K4[{11, 7, 7, 5}] = 7.22255799911e-11;
    K4[{11, 8, 0, 0}] = 1.08338369987e-10;
    K4[{11, 8, 1, 1}] = 1.08338369987e-10;
    K4[{11, 8, 2, 2}] = 9.63007733215e-11;
    K4[{11, 8, 3, 3}] = 9.63007733215e-11;
    K4[{11, 8, 6, 6}] = -1.50469958315e-11;
    K4[{11, 8, 8, 0}] = 1.32714503234e-09;
    K4[{11, 8, 8, 3}] = -2.70845924967e-11;
    K4[{11, 8, 8, 5}] = 3.0093991663e-11;
    K4[{11, 8, 8, 8}] = 1.50469958315e-11;
    K4[{11, 9, 9, 0}] = 2.82883521632e-10;
    K4[{11, 9, 9, 3}] = -2.70845924967e-11;
    K4[{11, 9, 9, 5}] = 4.21315883282e-11;
    K4[{11, 9, 9, 8}] = 2.40751933304e-11;
    K4[{11, 10, 0, 0}] = -8.72725758226e-11;
    K4[{11, 10, 1, 1}] = -8.72725758226e-11;
    K4[{11, 10, 2, 2}] = -6.92161808248e-11;
    K4[{11, 10, 3, 3}] = -7.22255799911e-11;
    K4[{11, 10, 6, 6}] = 1.80563949978e-11;
    K4[{11, 10, 7, 7}] = 1.20375966652e-11;
    K4[{11, 10, 8, 8}] = -5.41691849933e-11;
    K4[{11, 10, 9, 9}] = -5.71785841596e-11;
    K4[{11, 10, 10, 0}] = -9.96111124044e-10;
    K4[{11, 10, 10, 3}] = -3.0093991663e-11;
    K4[{11, 10, 10, 5}] = 5.41691849933e-11;
    K4[{11, 10, 10, 8}] = 1.50469958315e-11;
    K4[{11, 10, 10, 10}] = -6.01879833259e-11;
    K4[{11, 11, 0, 0}] = -5.60229748798e-08;
    K4[{11, 11, 1, 0}] = -1.50469958315e-11;
    K4[{11, 11, 1, 1}] = -5.70702457896e-08;
    K4[{11, 11, 2, 1}] = -2.58808328302e-10;
    K4[{11, 11, 2, 2}] = -5.33235438276e-08;
    K4[{11, 11, 3, 0}] = 2.61817727468e-10;
    K4[{11, 11, 3, 2}] = 2.10657941641e-11;
    K4[{11, 11, 3, 3}] = -5.35582769626e-08;
    K4[{11, 11, 4, 1}] = -8.72725758226e-11;
    K4[{11, 11, 4, 2}] = -4.32450660197e-09;
    K4[{11, 11, 4, 4}] = 7.09315383496e-09;
    K4[{11, 11, 5, 0}] = -9.93101724878e-11;
    K4[{11, 11, 5, 2}] = -1.50469958315e-11;
    K4[{11, 11, 5, 3}] = 2.77165663216e-09;
    K4[{11, 11, 5, 5}] = 1.26906362843e-08;
    K4[{11, 11, 6, 1}] = -6.31973824922e-11;
    K4[{11, 11, 6, 2}] = -1.9109684706e-09;
    K4[{11, 11, 6, 4}] = 7.18343580995e-09;
    K4[{11, 11, 6, 5}] = -1.80563949978e-11;
    K4[{11, 11, 6, 6}] = 9.23885544053e-09;
    K4[{11, 11, 7, 1}] = 3.75573015954e-09;
    K4[{11, 11, 7, 2}] = 8.125377749e-11;
    K4[{11, 11, 7, 4}] = 5.41691849933e-11;
    K4[{11, 11, 7, 6}] = 5.1159785827e-11;
    K4[{11, 11, 7, 7}] = 9.70531231131e-09;
    K4[{11, 11, 8, 0}] = 4.81503866607e-11;
    K4[{11, 11, 8, 3}] = 1.80563949978e-11;
    K4[{11, 11, 8, 8}] = 4.22519642948e-09;
    K4[{11, 11, 9, 1}] = 2.40751933304e-11;
    K4[{11, 11, 9, 2}] = -4.81503866607e-11;
    K4[{11, 11, 9, 9}] = 4.55021153944e-09;
    K4[{11, 11, 10, 0}] = -2.40751933304e-11;
    K4[{11, 11, 10, 3}] = -4.81503866607e-11;
    K4[{11, 11, 10, 8}] = -5.71785841596e-11;
    K4[{11, 11, 10, 9}] = 1.50469958315e-11;
    K4[{11, 11, 10, 10}] = 5.6606798318e-09;
    K4[{11, 11, 11, 0}] = 1.45053039815e-09;
    K4[{11, 11, 11, 3}] = -9.63007733215e-11;
    K4[{11, 11, 11, 5}] = 7.82443783237e-11;
  }

  void writeInputFiles() {
#ifdef MPI_PARALLEL
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (id == 0) {
      // Input file for the NMA calculation
      inputFileNMA.open(inputFileNMAName);
      inputFileNMA << "numModes = 30" << std::endl;
      inputFileNMA << "vscfIter = 50" << std::endl;
      inputFileNMA << "ONVector = 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
                      "0 0 0 0 0 0 0 "
                   << std::endl;
      inputFileNMA << "modeIds =1 2 3 4 5 6 7 8 91 2 3 4 5 6 7 8 9 10 11 12 13 "
                      "14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29"
                   << std::endl;
      inputFileNMA << "vscfEnTol = 10E-8" << std::endl;
      inputFileNMA << "vscfCoeffTol = 10E-12" << std::endl;
      inputFileNMA << "nMax = 6" << std::endl;
      inputFileNMA << "twoBodyTol = 1E-8" << std::endl;
      inputFileNMA << "numqp = 6" << std::endl;
      inputFileNMA << "fciDumpName = " << std::endl;
      inputFileNMA << "[Mode0]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency = 0.000173140740686" << std::endl;
      inputFileNMA << "[Mode1]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency =0.00037653737" << std::endl;
      inputFileNMA << "[Mode2]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency =0.0007083825588" << std::endl;
      inputFileNMA << "[Mode3]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency =0.0013093677" << std::endl;
      inputFileNMA << "[Mode4]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency =0.0018967887" << std::endl;
      inputFileNMA << "[Mode5]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency =0.00196409945" << std::endl;
      inputFileNMA << "[Mode6]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency =0.0028437865" << std::endl;
      inputFileNMA << "[Mode7]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency =0.002876278" << std::endl;
      inputFileNMA << "[Mode8]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency =0.0039633055" << std::endl;
      inputFileNMA << "[Mode9]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency =0.00449005883" << std::endl;
      inputFileNMA << "[Mode10]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency =0.0047207734" << std::endl;
      inputFileNMA << "[Mode11]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency = 0.0050602204" << std::endl;
      inputFileNMA << "[Mode12]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency = 0.0051769993" << std::endl;
      inputFileNMA << "[Mode13]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency = 0.0052947805" << std::endl;
      inputFileNMA << "[Mode14]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency =0.00574604" << std::endl;
      inputFileNMA << "[Mode15]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency =0.0063127114" << std::endl;
      inputFileNMA << "[Mode16]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency =0.0064383296" << std::endl;
      inputFileNMA << "[Mode17]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency =0.0065956598" << std::endl;
      inputFileNMA << "[Mode18]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency =0.006669564" << std::endl;
      inputFileNMA << "[Mode19]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency =0.0066978" << std::endl;
      inputFileNMA << "[Mode20]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency =0.0067252" << std::endl;
      inputFileNMA << "[Mode21]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency =0.0070805" << std::endl;
      inputFileNMA << "[Mode22]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency =0.0081057" << std::endl;
      inputFileNMA << "[Mode23]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency =0.013737" << std::endl;
      inputFileNMA << "[Mode24]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency = 0.013851" << std::endl;
      inputFileNMA << "[Mode25]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency =0.013979" << std::endl;
      inputFileNMA << "[Mode26]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency =0.014234" << std::endl;
      inputFileNMA << "[Mode27]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency =0.014243" << std::endl;
      inputFileNMA << "[Mode28]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency =0.014348" << std::endl;
      inputFileNMA << "[Mode29]" << std::endl;
      inputFileNMA << "numGrid = 10" << std::endl;
      inputFileNMA << "nquantum = 5" << std::endl;
      inputFileNMA << "frequency =0.01653" << std::endl;
      inputFileNMA.close();
      // Input file for CO2 calculation
      inputFileCO2.open(inputFileCO2Name);
      inputFileCO2 << "numModes = 4" << std::endl;
      inputFileCO2 << "vscfIter = 10" << std::endl;
      inputFileCO2 << "ONVector = 0 0 0 0" << std::endl;
      inputFileCO2 << "vscfEnTol = 10E-8" << std::endl;
      inputFileCO2 << "vscfCoeffTol = 10E-12" << std::endl;
      inputFileCO2 << "nModePotentialOrder = 2" << std::endl;
      inputFileCO2 << "nMax = 6" << std::endl;
      inputFileCO2 << "twoBodyTol = 1E-8" << std::endl;
      inputFileCO2 << "fciDumpName = " << std::endl;
      inputFileCO2 << "numqp = 8" << std::endl;
      inputFileCO2 << "[Mode0]" << std::endl;
      inputFileCO2 << "numGrid = 10" << std::endl;
      inputFileCO2 << "frequency = 0.00310191311" << std::endl;
      inputFileCO2 << "nquantum = 5" << std::endl;
      inputFileCO2 << "[Mode1]" << std::endl;
      inputFileCO2 << "numGrid = 10" << std::endl;
      inputFileCO2 << "frequency = 0.00310191311" << std::endl;
      inputFileCO2 << "nquantum = 5" << std::endl;
      inputFileCO2 << "[Mode2]" << std::endl;
      inputFileCO2 << "numGrid = 10" << std::endl;
      inputFileCO2 << "frequency = 0.00611968448" << std::endl;
      inputFileCO2 << "nquantum = 5" << std::endl;
      inputFileCO2 << "[Mode3]" << std::endl;
      inputFileCO2 << "numGrid = 10" << std::endl;
      inputFileCO2 << "frequency = 0.0108343253" << std::endl;
      inputFileCO2 << "nquantum = 5" << std::endl;
      inputFileCO2.close();
      // Input file for the Ethylene calculation
      inputFileEthylene.open(inputFileEthyleneName);
      inputFileEthylene << "numModes = 12" << std::endl;
      inputFileEthylene << "nModePotentialOrder = 2" << std::endl;
      inputFileEthylene << "vscfIter = 10" << std::endl;
      inputFileEthylene << "ONVector = 0 0 0 0 0 0 0 0 0 0 0 0" << std::endl;
      inputFileEthylene << "vscfEnTol = 10E-8" << std::endl;
      inputFileEthylene << "vscfCoeffTol = 10E-12" << std::endl;
      inputFileEthylene << "nMax = 6" << std::endl;
      inputFileEthylene << "twoBodyTol = 1E-8" << std::endl;
      inputFileEthylene << "numqp = 8" << std::endl;
      inputFileEthylene << "fciDumpName = " << std::endl;
      inputFileEthylene << "[Mode0]" << std::endl;
      inputFileEthylene << "numGrid = 10" << std::endl;
      inputFileEthylene << "nquantum = 5" << std::endl;
      inputFileEthylene << "frequency =0.014685663262168" << std::endl;
      inputFileEthylene << "[Mode1]" << std::endl;
      inputFileEthylene << "numGrid = 10" << std::endl;
      inputFileEthylene << "nquantum = 5" << std::endl;
      inputFileEthylene << "frequency =0.014556030306307" << std::endl;
      inputFileEthylene << "[Mode2]" << std::endl;
      inputFileEthylene << "numGrid = 10" << std::endl;
      inputFileEthylene << "nquantum = 5" << std::endl;
      inputFileEthylene << "frequency =0.014303788410238" << std::endl;
      inputFileEthylene << "[Mode3]" << std::endl;
      inputFileEthylene << "numGrid = 10" << std::endl;
      inputFileEthylene << "nquantum = 5" << std::endl;
      inputFileEthylene << "frequency =0.014240944722251" << std::endl;
      inputFileEthylene << "[Mode4]" << std::endl;
      inputFileEthylene << "numGrid = 10" << std::endl;
      inputFileEthylene << "nquantum = 5" << std::endl;
      inputFileEthylene << "frequency =0.00771284677269" << std::endl;
      inputFileEthylene << "[Mode5]" << std::endl;
      inputFileEthylene << "numGrid = 10" << std::endl;
      inputFileEthylene << "nquantum = 5" << std::endl;
      inputFileEthylene << "frequency =0.006737987331111" << std::endl;
      inputFileEthylene << "[Mode6]" << std::endl;
      inputFileEthylene << "numGrid = 10" << std::endl;
      inputFileEthylene << "nquantum = 5" << std::endl;
      inputFileEthylene << "frequency =0.006297367833175" << std::endl;
      inputFileEthylene << "[Mode7]" << std::endl;
      inputFileEthylene << "numGrid = 10" << std::endl;
      inputFileEthylene << "nquantum = 5" << std::endl;
      inputFileEthylene << "frequency =0.005679480491472" << std::endl;
      inputFileEthylene << "[Mode8]" << std::endl;
      inputFileEthylene << "numGrid = 10" << std::endl;
      inputFileEthylene << "nquantum = 5" << std::endl;
      inputFileEthylene << "frequency =0.004862488465868" << std::endl;
      inputFileEthylene << "[Mode9]" << std::endl;
      inputFileEthylene << "numGrid = 10" << std::endl;
      inputFileEthylene << "nquantum = 5" << std::endl;
      inputFileEthylene << "frequency =0.004482074996021" << std::endl;
      inputFileEthylene << "[Mode10]" << std::endl;
      inputFileEthylene << "numGrid = 10" << std::endl;
      inputFileEthylene << "nquantum = 5" << std::endl;
      inputFileEthylene << "frequency =0.004461219463429" << std::endl;
      inputFileEthylene << "[Mode11]" << std::endl;
      inputFileEthylene << "numGrid = 10" << std::endl;
      inputFileEthylene << "nquantum = 5" << std::endl;
      inputFileEthylene << "frequency =0.003809191330816" << std::endl;
      inputFileEthylene.close();
      // Input file for Christoffel calculation
      inputFileChristoffel.open(inputFileChristoffelName);
      inputFileChristoffel << "numModes = 3" << std::endl;
      inputFileChristoffel << "nMax =  5" << std::endl;
      inputFileChristoffel << "vscfIter = 100" << std::endl;
      inputFileChristoffel << "ONVector = 0 0 0" << std::endl;
      inputFileChristoffel << "vscfEnTol = 10E-8" << std::endl;
      inputFileChristoffel << "vscfCoeffTol = 10E-12" << std::endl;
      inputFileChristoffel << "twoBodyTol = 10E-9" << std::endl;
      inputFileChristoffel << "nModePotentialOrder = 2" << std::endl;
      inputFileChristoffel << "fciDumpName  = " << std::endl;
      inputFileChristoffel << "numqp = 8" << std::endl;
      inputFileChristoffel << "[Mode0]" << std::endl;
      inputFileChristoffel << "numGrid = 10" << std::endl;
      inputFileChristoffel << "frequency = 0.7" << std::endl;
      inputFileChristoffel << "nquantum = 5" << std::endl;
      inputFileChristoffel << "[Mode1]" << std::endl;
      inputFileChristoffel << "numGrid = 10" << std::endl;
      inputFileChristoffel << "frequency = 1.0" << std::endl;
      inputFileChristoffel << "nquantum = 5" << std::endl;
      inputFileChristoffel << "[Mode2]" << std::endl;
      inputFileChristoffel << "numGrid = 10" << std::endl;
      inputFileChristoffel << "frequency = 1.3" << std::endl;
      inputFileChristoffel << "nquantum = 5" << std::endl;
      inputFileChristoffel.close();
      // Input for the coupled oscillator calculation
      inputFileCoupledOscillator.open(inputFileCoupledOscillatorName);
      inputFileCoupledOscillator << "numModes = 6" << std::endl;
      inputFileCoupledOscillator << "nMax =  6" << std::endl;
      inputFileCoupledOscillator << "vscfIter = 100" << std::endl;
      inputFileCoupledOscillator << "ONVector = 0 0 0 0 0 0" << std::endl;
      inputFileCoupledOscillator << "vscfEnTol = 10E-8" << std::endl;
      inputFileCoupledOscillator << "vscfCoeffTol = 10E-12" << std::endl;
      inputFileCoupledOscillator << "twoBodyTol = 10E-9" << std::endl;
      inputFileCoupledOscillator << "nModePotentialOrder = 2" << std::endl;
      inputFileCoupledOscillator << "fciDumpName  = " << std::endl;
      inputFileCoupledOscillator << "numqp = 8" << std::endl;
      inputFileCoupledOscillator << "[Mode0]" << std::endl;
      inputFileCoupledOscillator << "numGrid = 20" << std::endl;
      inputFileCoupledOscillator << "frequency = 0.7071067811865476" << std::endl;
      inputFileCoupledOscillator << "nquantum = 5" << std::endl;
      inputFileCoupledOscillator << "[Mode1]" << std::endl;
      inputFileCoupledOscillator << "numGrid = 20" << std::endl;
      inputFileCoupledOscillator << "frequency = 1." << std::endl;
      inputFileCoupledOscillator << "nquantum = 5" << std::endl;
      inputFileCoupledOscillator << "[Mode2]" << std::endl;
      inputFileCoupledOscillator << "numGrid = 20" << std::endl;
      inputFileCoupledOscillator << "frequency = 1.224744871391589" << std::endl;
      inputFileCoupledOscillator << "nquantum = 5" << std::endl;
      inputFileCoupledOscillator << "[Mode3]" << std::endl;
      inputFileCoupledOscillator << "numGrid = 20" << std::endl;
      inputFileCoupledOscillator << "frequency = 1.4142135623730951" << std::endl;
      inputFileCoupledOscillator << "nquantum = 5" << std::endl;
      inputFileCoupledOscillator << "[Mode4]" << std::endl;
      inputFileCoupledOscillator << "numGrid = 20" << std::endl;
      inputFileCoupledOscillator << "frequency = 1.5811388300841898" << std::endl;
      inputFileCoupledOscillator << "nquantum = 5" << std::endl;
      inputFileCoupledOscillator << "[Mode5]" << std::endl;
      inputFileCoupledOscillator << "numGrid = 20" << std::endl;
      inputFileCoupledOscillator << "frequency = 1.7320508075688772" << std::endl;
      inputFileCoupledOscillator << "nquantum = 5" << std::endl;
    }
#ifdef MPI_PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  void removeInputFiles() {
#ifdef MPI_PARALLEL
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (id == 0) {
      auto successChristoffel = std::remove(inputFileChristoffelName.c_str());
      auto successCO2 = std::remove(inputFileCO2Name.c_str());
      if (successChristoffel != 0 || successCO2 != 0) {
        throw std::runtime_error("Error while removing input files");
      }
    }
#ifdef MPI_PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
};

TEST_F(VSCFTest, NMAComparisonWithReferencePaperEnergy) {
  initializeMPI();
  writeInputFiles();
  std::stringstream nmaPaper("12\n\n"
                             "H    2.723437   0.892771   0.013188 \n"
                             "H    2.130298  -0.539827  -0.886672 \n"
                             "H    2.125245  -0.541351   0.888778 \n"
                             "C    1.976395   0.084411   0.000797 \n"
                             "N    0.647510   0.681271  -0.004950 \n"
                             "H    0.549997   1.692012   0.009007 \n"
                             "O   -0.402564  -1.321377  -0.001241 \n"
                             "C   -0.479342  -0.095309   0.004434 \n"
                             "C   -1.806384   0.648329  -0.000575 \n"
                             "H   -2.370866   0.340581  -0.880167 \n"
                             "H   -2.374076   0.337436   0.891188 \n"
                             "H   -1.708633   1.734024  -0.002290 \n");
  nmaPesCalculator_ = std::make_shared<NMAPesCalculator>();
  auto refGeom = Utils::XyzStreamHandler::read(nmaPaper);
  nmaPesCalculator_->setStructure(refGeom);
  pesCartesianNma_ = std::make_shared<PesCartesian>(nmaPesCalculator_);
  std::shared_ptr<PES> pes = std::make_shared<PesCartesian>(nmaPesCalculator_);
  auto tmp = pes->getPES();
  if (id == 0) {
    EXPECT_THAT(tmp, DoubleNear(-248.38757324, 1.0E-4));
  }
  removeInputFiles();
}

TEST_F(VSCFTest, NMAComparisonWithReferenceEnergy) {
  writeInputFiles();
  VibrationalParameters parms;
  parms.setParameters(inputFileNMAName);
  std::stringstream nmaOpt("12\n\n"
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
  nmaPesCalculator_ = std::make_shared<NMAPesCalculator>();
  auto refGeom = Utils::XyzStreamHandler::read(nmaOpt);
  nmaPesCalculator_->setStructure(refGeom);
  pesCartesianNma_ = std::make_shared<PesCartesian>(nmaPesCalculator_);
  std::shared_ptr<PES> pes = std::make_shared<PesCartesian>(nmaPesCalculator_);
  auto tmp = pes->getPES();
  if (id == 0) {
    EXPECT_THAT(tmp, DoubleNear(-248.38880037, 1.0E-5));
  }
  removeInputFiles();
}

TEST_F(VSCFTest, NMAvscfGaussians) {
  writeInputFiles();
  VibrationalParameters parms;
  parms.setParameters(inputFileNMAName);
  parms.nModePotentialOrder_ = 1;
  if (id == 0) {
    parms.printParameters();
  }
  std::stringstream nmaOpt("12\n\n"
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
  nmaPesCalculator_ = std::make_shared<NMAPesCalculator>();
  auto refGeom = Utils::XyzStreamHandler::read(nmaOpt);
  nmaPesCalculator_->setStructure(refGeom);
  pesCartesianNma_ = std::make_shared<PesCartesian>(nmaPesCalculator_);
  std::shared_ptr<PES> pes = std::make_shared<PesCartesian>(nmaPesCalculator_);
  std::shared_ptr<ModalBasis> basis = std::make_shared<DistributedGaussians>(parms, pes);
  std::shared_ptr<VibrationalSCF> vscf = std::make_shared<VSCF>(parms, basis);
  auto* pesBase = pes.get();
  if (id == 0) {
    (dynamic_cast<PesCartesian*>(pesBase))->printHarmonicZPVE();
  }
  double vscfEnergy = 0;
  if (id == 0) {
    vscfEnergy = vscf->SCF();
  }
  // expect approx. half of sum of harmonic frequencies up to 1 percent
  // deviation
  if (id == 0) {
    EXPECT_THAT(vscfEnergy, DoubleNear(0.20182468333 / 2, 1E-2));
  }
  removeInputFiles();
}

TEST_F(VSCFTest, NMAvscfDVR) {
  writeInputFiles();
  VibrationalParameters parms;
  parms.setParameters(inputFileNMAName);
  parms.nModePotentialOrder_ = 2;
  if (id == 0) {
    parms.printParameters();
  }
  std::stringstream nmaOpt("12\n\n"
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
  nmaPesCalculator_ = std::make_shared<NMAPesCalculator>();
  auto refGeom = Utils::XyzStreamHandler::read(nmaOpt);
  nmaPesCalculator_->setStructure(refGeom);
  pesCartesianNma_ = std::make_shared<PesCartesian>(nmaPesCalculator_);
  std::shared_ptr<PES> pes = std::make_shared<PesCartesian>(nmaPesCalculator_);
  std::shared_ptr<ModalBasis> basis = std::make_shared<DVR>(parms, pes);
  std::shared_ptr<VibrationalSCF> vscf = std::make_shared<VSCF>(parms, basis);
  double vscfEnergy = NAN;
  if (id == 0) {
    vscfEnergy = vscf->SCF();
  }
  // expect approx. half of sum of harmonic frequencies up to 1 percent
  // deviation
  if (id == 0) {
    EXPECT_THAT(vscfEnergy, DoubleNear(0.20182468333 / 2, 1E-2));
  }
  removeInputFiles();
}

/** @brief VSCF energy with Christoffel potential + DVR */
TEST_F(VSCFTest, DVRChristoffel) {
  writeInputFiles();
  VibrationalParameters parms;
  parms.setParameters(inputFileChristoffelName);
  if (id == 0) {
    parms.printParameters();
  }
  std::shared_ptr<PES> pes = std::make_shared<PesChristoffel>();
  pes->setReferenceGeometry(k1Ch, k2Ch, k3Ch, k4Ch);
  std::shared_ptr<ModalBasis> basis = std::make_shared<DVR>(parms, pes);
  std::shared_ptr<VibrationalSCF> vscf = std::make_shared<VSCF>(parms, basis);
  double vscfEnergy = NAN;
  if (id == 0) {
    vscfEnergy = vscf->SCF();
  }
  if (id == 0) {
    EXPECT_THAT(vscfEnergy, DoubleNear(1.49503, 1.0E-5));
  }
  removeInputFiles();
}

/** @brief VSCF energy with Christoffel potential + DG */
TEST_F(VSCFTest, GaussiansChristoffel) {
  writeInputFiles();
  VibrationalParameters parms;
  parms.setParameters(inputFileChristoffelName);
  if (id == 0) {
    parms.printParameters();
  }
  std::shared_ptr<PES> pes = std::make_shared<PesChristoffel>();
  pes->setReferenceGeometry(k1Ch, k2Ch, k3Ch, k4Ch);
  std::shared_ptr<ModalBasis> basis = std::make_shared<DistributedGaussians>(parms, pes);
  std::shared_ptr<VibrationalSCF> vscf = std::make_shared<VSCF>(parms, basis);
  double vscfEnergy = NAN;
  if (id == 0) {
    vscfEnergy = vscf->SCF();
  }
  if (id == 0) {
    EXPECT_THAT(vscfEnergy, DoubleNear(1.49503, 1.0E-5));
  }
  removeInputFiles();
}

/** @brief VSCF energy with damped Christoffel potential + DG */
TEST_F(VSCFTest, GaussiansChristoffelDamped) {
  writeInputFiles();
  VibrationalParameters parms;
  parms.setParameters(inputFileChristoffelName);
  if (id == 0) {
    parms.printParameters();
  }
  std::shared_ptr<PES> pes = std::make_shared<PesChristoffel>();
  pes->setReferenceGeometry(k1Ch, k2Ch, k3ChD, k4Ch);
  std::shared_ptr<ModalBasis> basis = std::make_shared<DistributedGaussians>(parms, pes);
  std::shared_ptr<VibrationalSCF> vscf = std::make_shared<VSCF>(parms, basis);
  double vscfEnergy = NAN;
  if (id == 0) {
    vscfEnergy = vscf->SCF();
  }
  if (id == 0) {
    EXPECT_THAT(vscfEnergy, DoubleNear(1.4977876435, 1.0E-5));
  }
  removeInputFiles();
}

/** @brief VSCF energy of CO2 with Distributed Gaussians */
TEST_F(VSCFTest, GaussiansCO2) {
  writeInputFiles();
  VibrationalParameters parms;
  parms.setParameters(inputFileCO2Name);
  if (id == 0) {
    parms.printParameters();
  }
  std::shared_ptr<PES> pes = std::make_shared<PesQFF>();
  pes->setReferenceGeometry(k1Co2, k2Co2, k3Co2, k4Co2);
  std::shared_ptr<ModalBasis> basis = std::make_shared<DistributedGaussians>(parms, pes);
  std::shared_ptr<VibrationalSCF> vscf = std::make_shared<VSCF>(parms, basis);
  double vscfEnergy = NAN;
  if (id == 0) {
    vscfEnergy = vscf->SCF();
  }
  if (id == 0) {
    EXPECT_THAT(vscfEnergy, DoubleNear(0.0115909, 1.0E-7));
  }
  removeInputFiles();
}

/** @brief VSCF energy of ethylene + DG */
TEST_F(VSCFTest, GaussiansEthylene) {
  writeInputFiles();
  VibrationalParameters parms;
  parms.setParameters(inputFileEthyleneName);
  if (id == 0) {
    parms.printParameters();
  }
  std::shared_ptr<PES> pes = std::make_shared<PesQFF>();
  pes->setReferenceGeometry(K1, K2, K3, K4);
  std::shared_ptr<ModalBasis> basis = std::make_shared<DistributedGaussians>(parms, pes);
  std::shared_ptr<VibrationalSCF> vscf = std::make_shared<VSCF>(parms, basis);
  double vscfEnergy = NAN;
  if (id == 0) {
    vscfEnergy = vscf->SCF();
  }
  if (id == 0) {
    EXPECT_THAT(vscfEnergy, DoubleNear(0.0505168, 1.0E-7));
  }
  removeInputFiles();
}

/** @brief VSCF energy of ethylene + DVR */
TEST_F(VSCFTest, DVREthylene) {
  writeInputFiles();
  VibrationalParameters parms;
  parms.setParameters(inputFileEthyleneName);
  if (id == 0) {
    parms.printParameters();
  }
  std::shared_ptr<PES> pes = std::make_shared<PesQFF>();
  pes->setReferenceGeometry(K1, K2, K3, K4);
  std::shared_ptr<ModalBasis> basis = std::make_shared<DVR>(parms, pes);
  std::shared_ptr<VibrationalSCF> vscf = std::make_shared<VSCF>(parms, basis);
  double vscfEnergy = NAN;
  if (id == 0) {
    vscfEnergy = vscf->SCF();
  }
  if (id == 0) {
    EXPECT_THAT(vscfEnergy, DoubleNear(0.0505168, 1.0E-7));
  }
  removeInputFiles();
}

/** @brief VSCF energy of ethylene + DVR to check mode couplig settings */
TEST_F(VSCFTest, DVREthyleneCheckCoupledModes) {
  writeInputFiles();
  VibrationalParameters parms;
  parms.setParameters(inputFileEthyleneName);
  // First, perform a MC=1 (so uncoupled) calculation
  parms.nModePotentialOrder_ = 1;
  if (id == 0) {
    parms.printParameters();
  }
  std::shared_ptr<PES> pes = std::make_shared<PesQFF>();
  pes->setReferenceGeometry(K1, K2, K3, K4);
  std::shared_ptr<ModalBasis> basis = std::make_shared<DVR>(parms, pes);
  std::shared_ptr<VibrationalSCF> vscf = std::make_shared<VSCF>(parms, basis);
  double vscfEnergy = NAN;
  if (id == 0) {
    vscfEnergy = vscf->SCF();
  }

  // Now, perform a MC=2 (two-mode coupled) calculation, but only add one single
  // mode to the coupled ones
  parms.nModePotentialOrder_ = 2;
  parms.coupledModes_.resize(0);
  parms.coupledModes_.push_back(3);

  if (id == 0) {
    parms.printParameters();
  }
  std::shared_ptr<PES> pes2 = std::make_shared<PesQFF>();
  pes2->setReferenceGeometry(K1, K2, K3, K4);
  std::shared_ptr<ModalBasis> basis2 = std::make_shared<DVR>(parms, pes2);
  std::shared_ptr<VibrationalSCF> vscf2 = std::make_shared<VSCF>(parms, basis2);
  double vscfEnergy2 = NAN;
  if (id == 0) {
    vscfEnergy2 = vscf2->SCF();
  }

  // Now check that this results in the same energy
  if (id == 0) {
    EXPECT_THAT(vscfEnergy2, DoubleNear(vscfEnergy, 1.0E-7));
  }

  // And finally, perform a fully coupled calc with MC=2
  parms.coupledModes_.resize(0);
  for (int mode = 0; mode < parms.numModes_; mode++) {
    parms.coupledModes_.push_back(mode);
  }

  if (id == 0) {
    parms.printParameters();
  }
  std::shared_ptr<PES> pes3 = std::make_shared<PesQFF>();
  pes3->setReferenceGeometry(K1, K2, K3, K4);
  std::shared_ptr<ModalBasis> basis3 = std::make_shared<DVR>(parms, pes3);
  std::shared_ptr<VibrationalSCF> vscf3 = std::make_shared<VSCF>(parms, basis3);
  double vscfEnergy3 = NAN;
  if (id == 0) {
    vscfEnergy3 = vscf3->SCF();
  }
  if (id == 0) {
    EXPECT_THAT(vscfEnergy3, DoubleNear(0.0505168, 1.0E-7));
  }
  removeInputFiles();
}

/**
 * @brief Test of VSCF on the linearly coupled Harmonic oscillator model
 * Hamiltonian with no coupling parameter.
 */
TEST_F(VSCFTest, GaussiansOscillatorUncoupled) {
  writeInputFiles();
  VibrationalParameters parms;
  parms.setParameters(inputFileCoupledOscillatorName);
  if (id == 0) {
    parms.printParameters();
  }
  // Uncoupled version
  std::shared_ptr<PES> pes = std::make_shared<PesCoupledOscillator>(parms, 0.);
  std::shared_ptr<ModalBasis> basis = std::make_shared<DistributedGaussians>(parms, pes);
  std::shared_ptr<VibrationalSCF> vscf = std::make_shared<VSCF>(parms, basis);
  double vscfEnergy = NAN;
  if (id == 0) {
    vscfEnergy = vscf->SCF();
  }
  if (id == 0) {
    EXPECT_THAT(vscfEnergy, DoubleNear(3.8296274263021495, 1.0E-05));
  }
  removeInputFiles();
}

/**
 * @brief Test of VSCF on the linearly coupled Harmonic oscillator model
 * Hamiltonian with coupling parameter = 0.1
 */

TEST_F(VSCFTest, DVROscillatorCoupled) {
  writeInputFiles();
  VibrationalParameters parms;
  parms.setParameters(inputFileCoupledOscillatorName);
  if (id == 0) {
    parms.printParameters();
  }
  // Coupled version
  std::shared_ptr<PES> pes = std::make_shared<PesCoupledOscillator>(parms, 0.1);
  std::shared_ptr<ModalBasis> basis = std::make_shared<DVR>(parms, pes);
  std::shared_ptr<VibrationalSCF> vscf = std::make_shared<VSCF>(parms, basis);
  double vscfEnergy = NAN;
  if (id == 0) {
    vscfEnergy = vscf->SCF();
  }
  if (id == 0) {
    EXPECT_THAT(vscfEnergy, DoubleNear(3.8296274263021495, 1.0E-04));
  }
  removeInputFiles();
}

/**
 * @brief Test of VSCF on the linearly coupled Harmonic oscillator model
 * Hamiltonian with coupling parameter = 0.1
 */

TEST_F(VSCFTest, GaussiansOscillatorCoupled) {
  writeInputFiles();
  VibrationalParameters parms;
  parms.setParameters(inputFileCoupledOscillatorName);
  if (id == 0) {
    parms.printParameters();
  }
  // Coupled version
  std::shared_ptr<PES> pes = std::make_shared<PesCoupledOscillator>(parms, 0.1);
  std::shared_ptr<ModalBasis> basis = std::make_shared<DistributedGaussians>(parms, pes);
  std::shared_ptr<VibrationalSCF> vscf = std::make_shared<VSCF>(parms, basis);
  double vscfEnergy = NAN;
  if (id == 0) {
    vscfEnergy = vscf->SCF();
  }
  if (id == 0) {
    EXPECT_THAT(vscfEnergy, DoubleNear(3.8296274263021495, 1.0E-05));
  }
  removeInputFiles();
  finalizeMPI();
}

} // namespace Scine::Colibri
