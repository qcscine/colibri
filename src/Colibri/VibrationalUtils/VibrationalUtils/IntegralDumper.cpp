/**
 * @file IntegralDumper.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#include "IntegralDumper.h"
#include "Utils/Constants.h"
#include <sys/stat.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <utility>

namespace Scine::Colibri {

void IntegralDumper::writeIntegrals(modalOneBodyMap& oneBody, modalTwoBodyMap& twoBody, const std::vector<int>& nMax,
                                    int nModes, const std::string& filename, bool onlyDumpCoupledModes,
                                    const std::vector<int>& coupledModes) {
  // if file does not exist, open
  std::ofstream fciDump;
  if (!fileExists(filename)) {
    try {
      fciDump.open(filename, std::ofstream::out | std::ofstream::app);
      fciDump << std::fixed << std::setprecision(8);
      writeOneBody(fciDump, oneBody, nModes, nMax, onlyDumpCoupledModes, coupledModes);
      std::cout << "OneBody written " << std::endl;
      writeTwoBody(fciDump, twoBody, nModes, nMax, onlyDumpCoupledModes, coupledModes);
      std::cout << "TwoBody written " << std::endl;
      std::cout << "Integral dump file written." << std::endl;
      fciDump.close();
    }
    catch (const std::ofstream::failure& e) {
      std::cout << "Error creating FCIDUMP file." << std::endl;
    }
  } else {
    std::cout << "File " << filename << " already exists." << std::endl;
  }
}

void IntegralDumper::writeIntegrals(modalOneBodyMap& oneBody, modalTwoBodyMap& twoBody, modalThreeBodyMap& threeBody,
                                    const std::vector<int>& nMax, int nModes, const std::string& filename,
                                    bool onlyDumpCoupledModes, const std::vector<int>& coupledModes) {
  // if file does not exist, open
  std::ofstream fciDump;
  if (!fileExists(filename)) {
    try {
      fciDump.open(filename, std::ofstream::out | std::ofstream::app);
      fciDump << std::fixed << std::setprecision(8);
      writeOneBody(fciDump, oneBody, nModes, nMax, onlyDumpCoupledModes, coupledModes);
      std::cout << "OneBody written " << std::endl;
      writeTwoBody(fciDump, twoBody, nModes, nMax, onlyDumpCoupledModes, coupledModes);
      std::cout << "TwoBody written " << std::endl;
      writeThreeBody(fciDump, threeBody, nModes, nMax, onlyDumpCoupledModes, coupledModes);
      std::cout << "ThreeBody written " << std::endl;
      std::cout << "Integral dump file written." << std::endl;
      fciDump.close();
    }
    catch (const std::ofstream::failure& e) {
      std::cout << "Error creating FCIDUMP file." << std::endl;
    }
  } else {
    std::cout << "File " << filename << " already exists." << std::endl;
  }
}

void IntegralDumper::writeIntegrals(modalOneBodyMap& oneBody, std::vector<int> nMax, int nModes,
                                    const std::string& filename, bool onlyDumpCoupledModes, std::vector<int> coupledModes) {
  // if file does not exist, open
  std::ofstream fciDump;
  if (!fileExists(filename)) {
    try {
      fciDump.open(filename, std::ofstream::out | std::ofstream::app);
      fciDump << std::fixed << std::setprecision(8);
      writeOneBody(fciDump, oneBody, nModes, std::move(nMax), onlyDumpCoupledModes, std::move(coupledModes));
      fciDump.close();
    }
    catch (const std::ofstream::failure& e) {
      std::cout << "Error creating FCIDUMP file." << std::endl;
    }
  } else {
    std::cout << "File " << filename << " already exists." << std::endl;
  }
}

bool IntegralDumper::fileExists(const std::string& filename) {
  struct stat buf;
  return stat(filename.c_str(), &buf) != -1;
}

void IntegralDumper::writeOneBody(std::ofstream& fciDump, modalOneBodyMap& oneBody, int nModes, std::vector<int> nMax,
                                  bool onlyDumpCoupledModes, std::vector<int> coupledModes) {
  int modeNumber = 0;
  for (int mode = 0; mode < nModes; mode++) {
    if (!onlyDumpCoupledModes || std::find(coupledModes.begin(), coupledModes.end(), mode) != coupledModes.end()) {
      modeNumber++;
      for (int modal1 = 0; modal1 < nMax[mode]; modal1++) {
        for (int modal2 = 0; modal2 < nMax[mode]; modal2++) {
          std::array<int, 3> inds = {mode, modal1, modal2};
          if (std::abs(oneBody[inds]) > 1.0E-12) {
            fciDump << std::to_string(modeNumber) + "-" + std::to_string(inds[1]) << "    ";
            fciDump << std::to_string(modeNumber) + "-" + std::to_string(inds[2]) << "    ";
            fciDump << oneBody[inds] * Utils::Constants::invCentimeter_per_hartree << std::endl;
          }
        }
      }
    }
  }
}

void IntegralDumper::writeTwoBody(std::ofstream& fciDump, modalTwoBodyMap& twoBody, int nModes, std::vector<int> nMax,
                                  bool onlyDumpCoupledModes, std::vector<int> coupledModes) {
  int modeiNumber = 0;
  for (int modei = 0; modei < nModes; modei++) {
    if (!onlyDumpCoupledModes || std::find(coupledModes.begin(), coupledModes.end(), modei) != coupledModes.end()) {
      modeiNumber++;
      int modejNumber = 0;
      for (int modej = 0; modej < nModes; modej++) {
        if (!onlyDumpCoupledModes || std::find(coupledModes.begin(), coupledModes.end(), modej) != coupledModes.end()) {
          modejNumber++;
          if (modej != modei) {
            for (int modali1 = 0; modali1 < nMax[modei]; modali1++) {
              for (int modali2 = 0; modali2 < nMax[modei]; modali2++) {
                for (int modalj1 = 0; modalj1 < nMax[modej]; modalj1++) {
                  for (int modalj2 = 0; modalj2 < nMax[modej]; modalj2++) {
                    std::array<int, 6> ind = {modei, modali1, modali2, modej, modalj1, modalj2};
                    if (std::abs(twoBody[ind]) > 1.0E-12) {
                      fciDump << std::fixed << std::setprecision(8);
                      fciDump << std::to_string(modeiNumber) + "-" + std::to_string(modali1) + "  " +
                                     std::to_string(modeiNumber) + "-" + std::to_string(modali2);
                      fciDump << "  " << std::to_string(modejNumber) + "-" + std::to_string(modalj1) + "  "
                              << std::to_string(modejNumber) + "-" + std::to_string(modalj2);
                      fciDump << "  ";
                      fciDump << twoBody[ind] * Utils::Constants::invCentimeter_per_hartree << std::endl;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

void IntegralDumper::writeThreeBody(std::ofstream& fciDump, modalThreeBodyMap& threeBody, int nModes,
                                    std::vector<int> nMax, bool onlyDumpCoupledModes, std::vector<int> coupledModes) {
  int modeiNumber = 0;
  for (int modei = 0; modei < nModes; modei++) {
    if (!onlyDumpCoupledModes || std::find(coupledModes.begin(), coupledModes.end(), modei) != coupledModes.end()) {
      modeiNumber++;
      int modejNumber = 0;
      for (int modej = 0; modej < nModes; modej++) {
        if (!onlyDumpCoupledModes || std::find(coupledModes.begin(), coupledModes.end(), modej) != coupledModes.end()) {
          modejNumber++;
          int modekNumber = 0;
          for (int modek = 0; modek < nModes; modek++) {
            if (!onlyDumpCoupledModes || std::find(coupledModes.begin(), coupledModes.end(), modek) != coupledModes.end()) {
              modekNumber++;
              if (modej != modei && modej != modek && modek != modei) {
                for (int modali1 = 0; modali1 < nMax[modei]; modali1++) {
                  for (int modali2 = 0; modali2 < nMax[modei]; modali2++) {
                    for (int modalj1 = 0; modalj1 < nMax[modej]; modalj1++) {
                      for (int modalj2 = 0; modalj2 < nMax[modej]; modalj2++) {
                        for (int modalk1 = 0; modalk1 < nMax[modek]; modalk1++) {
                          for (int modalk2 = 0; modalk2 < nMax[modek]; modalk2++) {
                            std::array<int, 9> ind = {modei,   modali1, modali2, modej,  modalj1,
                                                      modalj2, modek,   modalk1, modalk2};
                            if (std::abs(threeBody[ind]) > 1.0E-12) {
                              fciDump << std::fixed << std::setprecision(8);
                              fciDump << std::to_string(modeiNumber) + "-" + std::to_string(modali1) + "  " +
                                             std::to_string(modeiNumber) + "-" + std::to_string(modali2);
                              fciDump << "  " << std::to_string(modejNumber) + "-" + std::to_string(modalj1) + "  "
                                      << std::to_string(modejNumber) + "-" + std::to_string(modalj2);
                              fciDump << "  " << std::to_string(modekNumber) + "-" + std::to_string(modalk1) + "  "
                                      << std::to_string(modekNumber) + "-" + std::to_string(modalk2);
                              fciDump << "  ";
                              fciDump << threeBody[ind] * Utils::Constants::invCentimeter_per_hartree << std::endl;
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

} // namespace Scine::Colibri
