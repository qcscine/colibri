/**
 * @file writeOTFInpFiles.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

namespace Scine {
namespace Colibri {
namespace OTFVSCFTestInputs {

// Create input file for ethylene
inline void writeInputEthene(const std::string& fileName) {
  std::ofstream inputFile;
  inputFile.open(fileName);
  inputFile << "numModes = 12" << std::endl;
  inputFile << "nModePotentialOrder = 2" << std::endl;
  inputFile << "XYZFile = tmp_ethene_struct.xyz" << std::endl;
  inputFile << "program = Sparrow" << std::endl;
  inputFile << "method_family = DFTB3" << std::endl;
  inputFile << "vscfIter = 10" << std::endl;
  inputFile << "ONVector = 0 0 0 0 0 0 0 0 0 0 0 0" << std::endl;
  inputFile << "vscfEnTol = 10E-8" << std::endl;
  inputFile << "vscfCoeffTol = 10E-12" << std::endl;
  inputFile << "nMax = 6" << std::endl;
  inputFile << "twoBodyTol = 1E-9" << std::endl;
  inputFile << "numqp = 8" << std::endl;
  inputFile << "inputPESType = QFF" << std::endl;
  inputFile << "primBasisType = DVR" << std::endl;
  inputFile << "fciDumpName = tmp_FCIDUMPv_otf2" << std::endl;
  inputFile << "[Mode0]" << std::endl;
  inputFile << "numGrid = 10" << std::endl;
  inputFile << "nquantum = 5" << std::endl;
  inputFile << "[Mode1]" << std::endl;
  inputFile << "numGrid = 10" << std::endl;
  inputFile << "nquantum = 5" << std::endl;
  inputFile << "[Mode2]" << std::endl;
  inputFile << "numGrid = 10" << std::endl;
  inputFile << "nquantum = 5" << std::endl;
  inputFile << "[Mode3]" << std::endl;
  inputFile << "numGrid = 10" << std::endl;
  inputFile << "nquantum = 5" << std::endl;
  inputFile << "[Mode4]" << std::endl;
  inputFile << "numGrid = 10" << std::endl;
  inputFile << "nquantum = 5" << std::endl;
  inputFile << "[Mode5]" << std::endl;
  inputFile << "numGrid = 11" << std::endl;
  inputFile << "nquantum = 5" << std::endl;
  inputFile << "[Mode6]" << std::endl;
  inputFile << "numGrid = 10" << std::endl;
  inputFile << "nquantum = 5" << std::endl;
  inputFile << "[Mode7]" << std::endl;
  inputFile << "numGrid = 10" << std::endl;
  inputFile << "nquantum = 5" << std::endl;
  inputFile << "[Mode8]" << std::endl;
  inputFile << "numGrid = 10" << std::endl;
  inputFile << "nquantum = 5" << std::endl;
  inputFile << "[Mode9]" << std::endl;
  inputFile << "numGrid = 10" << std::endl;
  inputFile << "nquantum = 5" << std::endl;
  inputFile << "[Mode10]" << std::endl;
  inputFile << "numGrid = 10" << std::endl;
  inputFile << "nquantum = 5" << std::endl;
  inputFile << "[Mode11]" << std::endl;
  inputFile << "numGrid = 10" << std::endl;
  inputFile << "nquantum = 5" << std::endl;
  inputFile.close();
}

// Create XYZ file for ethylene
inline void writeXYZEthene(const std::string& fileName) {
  std::ofstream inputFile;
  inputFile.open(fileName);
  inputFile << "6" << std::endl << std::endl;
  inputFile << "C    -0.6488329   -0.1535148   -0.0114205 " << std::endl;
  inputFile << "C     0.6488347    0.1535232    0.0114132 " << std::endl;
  inputFile << "H    -1.3348602    0.3030364   -0.7318331 " << std::endl;
  inputFile << "H    -1.0812490   -0.8747297    0.6892348 " << std::endl;
  inputFile << "H     1.0812711    0.8747384   -0.6892399 " << std::endl;
  inputFile << "H     1.3348363   -0.3030536    0.7318453" << std::endl;
  inputFile.close();
}

// Create harmfreqs file for ethylene
inline void writeHarmFreqsEthene(const std::string& fileName) {
  std::ofstream inputFile;
  inputFile.open(fileName);
  inputFile << "     #    cm^-1" << std::endl;
  inputFile << "     1   +837.767722" << std::endl;
  inputFile << "     2   +838.260541" << std::endl;
  inputFile << "     3   +917.898192" << std::endl;
  inputFile << "     4  +1042.493807" << std::endl;
  inputFile << "     5  +1220.403878" << std::endl;
  inputFile << "     6  +1284.435391" << std::endl;
  inputFile << "     7  +1360.888196" << std::endl;
  inputFile << "     8  +1774.984947" << std::endl;
  inputFile << "     9  +2953.621095" << std::endl;
  inputFile << "    10  +2955.428413" << std::endl;
  inputFile << "    11  +3050.908270" << std::endl;
  inputFile << "    12  +3081.329876" << std::endl;
  inputFile.close();
}

// Create input file for hydrogen
inline void writeInputH2(const std::string& fileName) {
  std::ofstream inputFile;
  inputFile.open(fileName);
  inputFile << "numModes = 1" << std::endl;
  inputFile << "nModePotentialOrder = 1" << std::endl;
  inputFile << "XYZFile = tmp_h2_struct.xyz" << std::endl;
  inputFile << "program = Turbomole" << std::endl;
  inputFile << "method_family = DFT" << std::endl;
  inputFile << "vscfIter = 10" << std::endl;
  inputFile << "ONVector = 0" << std::endl;
  inputFile << "vscfEnTol = 10E-8" << std::endl;
  inputFile << "vscfCoeffTol = 10E-12" << std::endl;
  inputFile << "nMax = 6" << std::endl;
  inputFile << "numqp = 8" << std::endl;
  inputFile << "primBasisType = DVR" << std::endl;
  inputFile << "fciDumpName = FCIDUMPv_otf1" << std::endl;
  inputFile << "vciExModes = 1" << std::endl;
  inputFile << "vciTotEx = 6" << std::endl;
  inputFile << "[Mode0]" << std::endl;
  inputFile << "numGrid = 10" << std::endl;
  inputFile << "nquantum = 5" << std::endl;
  inputFile.close();
}

// Create XYZ file for hydrogen
inline void writeXYZforH2(const std::string& fileName) {
  std::ofstream inputFile;
  inputFile.open(fileName);
  inputFile << "2" << std::endl << std::endl;
  inputFile << "H     0.3700000    0.0000000    0.0000000 " << std::endl;
  inputFile << "H    -0.3700000    0.0000000    0.0000000" << std::endl;
  inputFile.close();
}

// Create harmfreqs file for h2
inline void writeHarmFreqsH2(const std::string& fileName) {
  std::ofstream inputFile;
  inputFile.open(fileName);
  inputFile << "     #    cm^-1" << std::endl;
  inputFile << "   1  +4343.184861" << std::endl;
  inputFile.close();
}

// Create input file for H2O
inline void writeInputH2O(const std::string& fileName) {
  std::ofstream inputFile;
  inputFile.open(fileName);
  inputFile << "numModes = 3" << std::endl;
  inputFile << "nModePotentialOrder = 3" << std::endl;
  inputFile << "XYZFile = tmp_water_struct.xyz" << std::endl;
  inputFile << "program = Sparrow" << std::endl;
  inputFile << "method_family = DFTB3" << std::endl;
  inputFile << "vscfIter = 10" << std::endl;
  inputFile << "ONVector = 0 0 0" << std::endl;
  inputFile << "vscfEnTol = 10E-8" << std::endl;
  inputFile << "vscfCoeffTol = 10E-12" << std::endl;
  inputFile << "nMax = 6" << std::endl;
  inputFile << "numqp = 8" << std::endl;
  inputFile << "primBasisType = DVR" << std::endl;
  inputFile << "fciDumpName = tmp_FCIDUMPv_h2o" << std::endl;
  inputFile << "vciExModes = 1" << std::endl;
  inputFile << "vciTotEx = 6" << std::endl;
  inputFile << "[Mode0]" << std::endl;
  inputFile << "numGrid = 10" << std::endl;
  inputFile << "nquantum = 5" << std::endl;
  inputFile << "[Mode1]" << std::endl;
  inputFile << "numGrid = 10" << std::endl;
  inputFile << "nquantum = 5" << std::endl;
  inputFile << "[Mode2]" << std::endl;
  inputFile << "numGrid = 10" << std::endl;
  inputFile << "nquantum = 5" << std::endl;
  inputFile.close();
}

// Create XYZ file for water
inline void writeXYZforH20(const std::string& fileName) {
  std::ofstream inputFile;
  inputFile.open(fileName);
  inputFile << "3" << std::endl << std::endl;
  inputFile << "O     0.0000000    0.0000000    0.0000000 " << std::endl;
  inputFile << "H     0.5884000    0.7585000    0.0000000 " << std::endl;
  inputFile << "H     0.5884000   -0.7585000    0.0000000" << std::endl;
  inputFile.close();
}

// Create harmfreqs file for water
inline void writeHarmFreqsH2O(const std::string& fileName) {
  std::ofstream inputFile;
  inputFile.open(fileName);
  inputFile << "     #    cm^-1" << std::endl;
  inputFile << "     1  +1385.052523" << std::endl;
  inputFile << "     2  +3630.680661" << std::endl;
  inputFile << "     3  +3903.446308" << std::endl;
  inputFile.close();
}

} // namespace OTFVSCFTestInputs
} // namespace Colibri
} // namespace Scine
