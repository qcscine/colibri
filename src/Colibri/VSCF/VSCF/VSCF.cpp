/**
 * @file VSCF.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "VSCF.h"

#ifdef OTF
  #include "VibrationalUtils/PESLibrary/HybridOTFfromSPs.h"
  #include "VibrationalUtils/PESLibrary/OnTheFlyPes.h"
#endif

#include "Utils/Constants.h"
#include "VibrationalUtils/PESLibrary/PesFromSPs.h"
#include <omp.h>
#include <Eigen/Eigenvalues>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <random>
#include <utility>

#ifdef MPI_PARALLEL
  #include <mpi.h>
  #include <cstddef>
  #include <queue>
#endif

namespace Scine::Colibri {

VSCF::VSCF(const VibrationalParameters& parms, std::shared_ptr<ModalBasis> basis)
  : Base(parms), basis_(std::move(std::move(basis))) {
  set3B_ = parms.set3B_;
  bool areTheyStored = storeIntegrals();
}

bool VSCF::storeIntegrals() {
  oneBodyIntegrals_.reserve(nModes);
  oneBodyHarmonicIntegrals_.reserve(nModes);
  overlaps_.reserve(nModes);
  bool success = 0;
  if (potentialOrder == 1) {
    try {
      this->storeOneBodyIntegrals();
      success = true;
    }
    catch (const std::bad_alloc&) {
      std::cout << "Running out of memory." << std::endl;
      success = false;
    }
  } else if (potentialOrder == 0) {
    try {
      this->storeHarmonicIntegrals();
      return true;
    }
    catch (const std::bad_alloc&) {
      std::cout << "Running out of memory." << std::endl;
      success = false;
    }
  } else if (potentialOrder == 2) {
    try {
      this->storeOneBodyIntegrals();
      this->storeTwoBodyIntegrals();
      success = true;
    }
    catch (const std::bad_alloc&) {
      std::cout << "Running out of memory." << std::endl;
      success = false;
    }
  } else if (potentialOrder == 3) {
    try {
      this->storeOneBodyIntegrals();
      this->storeTwoBodyIntegrals();
      this->storeThreeBodyIntegrals();
      success = true;
    }
    catch (const std::bad_alloc&) {
      std::cout << "Running out of memory." << std::endl;
      success = false;
    }
  } else if (potentialOrder == 5) {
    try {
      this->storeOneBodyIntegrals();
      this->storeHarmonicIntegrals();
      this->storeTwoBodyIntegrals();
      success = true;
    }
    catch (const std::bad_alloc&) {
      std::cout << "Running out of memory." << std::endl;
      success = false;
    }
  } else {
    std::cout << "Order of coupling in the potential not implemented." << std::endl;
    success = false;
  }
#ifdef MPI_PARALLEL
  #ifdef OTF
  // Only do this part if we are doing an MPI-parallel on-the-fly PES
  // calculation
  if (std::dynamic_pointer_cast<OTFPes>(basis_->pes_) || std::dynamic_pointer_cast<HybridPes>(basis_->pes_)) {
    int numSpCalcs;
    int numSpCalcsLocal = basis_->pes_->getNumSPsForVSCF();
    MPI_Allreduce(&numSpCalcsLocal, &numSpCalcs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    basis_->pes_->setNumSPsForVSCF(numSpCalcs);
  }
  #endif
#endif
  return success;
}

void VSCF::storeHarmonicIntegrals() {
  for (int mode = 0; mode < nModes; mode++) {
    int nBasis = basis_->getBasisSize(mode);
    Eigen::MatrixXd harmonic = Eigen::MatrixXd::Zero(nBasis, nBasis);
    for (int i = 0; i < nBasis; i++) {
      for (int j = 0; j < nBasis; j++) {
        if (i >= j) {
          harmonic(i, j) = basis_->getHarmonicIntegrals(i, j, mode);
        }
      }
    }
    Eigen::MatrixXd transpH = harmonic.triangularView<Eigen::StrictlyLower>();
    transpH.transposeInPlace();
    oneBodyHarmonicIntegrals_.emplace_back(harmonic + transpH);
  }
}

#ifndef MPI_PARALLEL
void VSCF::storeOneBodyIntegrals() {
  overlaps_.resize(nModes);
  oneBodyIntegrals_.resize(nModes);
  #pragma omp parallel for
  for (int mode = 0; mode < nModes; mode++) {
    int nBasis = basis_->getBasisSize(mode);
    Eigen::MatrixXd onebody = Eigen::MatrixXd::Zero(nBasis, nBasis);
    Eigen::MatrixXd overlap = Eigen::MatrixXd::Zero(nBasis, nBasis);
    for (int i = 0; i < nBasis; i++) {
      for (int j = 0; j < nBasis; j++) {
        if (i >= j) {
          onebody(i, j) = basis_->getOneBodyIntegrals(i, j, mode);
          overlap(i, j) = basis_->getOverlap(i, j, mode);
        }
      }
    }
    Eigen::MatrixXd transp = onebody.triangularView<Eigen::StrictlyLower>();
    transp.transposeInPlace();
  #pragma omp critical
    { oneBodyIntegrals_[mode] = (onebody + transp); }

    Eigen::MatrixXd transpS = overlap.triangularView<Eigen::StrictlyLower>();
    transpS.transposeInPlace();
  #pragma omp critical
    { overlaps_[mode] = (overlap + transpS); }
  }
}
#else
void VSCF::storeOneBodyIntegrals() {
  double start1B = MPI_Wtime();
  int size, id;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size == 1) { // No need to do complicated task splitting for one-thread MPI job
    std::cout << "One-thread MPI job" << std::endl;
    overlaps_.resize(nModes);
    oneBodyIntegrals_.resize(nModes);
    for (int mode = 0; mode < nModes; mode++) {
      int nBasis = basis_->getBasisSize(mode);
      Eigen::MatrixXd onebody = Eigen::MatrixXd::Zero(nBasis, nBasis);
      Eigen::MatrixXd overlap = Eigen::MatrixXd::Zero(nBasis, nBasis);
      for (int i = 0; i < nBasis; i++) {
        for (int j = 0; j < nBasis; j++) {
          if (i >= j) {
            onebody(i, j) = basis_->getOneBodyIntegrals(i, j, mode);
            overlap(i, j) = basis_->getOverlap(i, j, mode);
          }
        }
      }
      Eigen::MatrixXd transp = onebody.triangularView<Eigen::StrictlyLower>();
      transp.transposeInPlace();
      oneBodyIntegrals_[mode] = (onebody + transp);

      Eigen::MatrixXd transpS = overlap.triangularView<Eigen::StrictlyLower>();
      transpS.transposeInPlace();
      overlaps_[mode] = (overlap + transpS);
    }
  } else {
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    // Prepare so all threads have everything they need
    std::vector<Eigen::MatrixXd> oneBodyIntegralsL, overlapsL;
    overlaps_.resize(nModes);
    oneBodyIntegrals_.resize(nModes);
    overlapsL.resize(nModes);
    oneBodyIntegralsL.resize(nModes);
    for (int mode = 0; mode < nModes; mode++) {
      int nBasis = basis_->getBasisSize(mode);
      Eigen::MatrixXd dummy = Eigen::MatrixXd::Zero(nBasis, nBasis);
      oneBodyIntegralsL[mode] = dummy;
      oneBodyIntegrals_[mode] = dummy;
      overlapsL[mode] = dummy;
      overlaps_[mode] = dummy;
    }
    double beforeCalcs;
    if (id == 0) { // Master thread is distributing work
      // Generate queue of tasks. Each calculation is specified by {mode,
      // modali, modalj}
      std::queue<std::tuple<int, int, int>> calcsToDo;
      for (int mode = 0; mode < nModes; mode++) {
        int nBasis = basis_->getBasisSize(mode);
        for (int i = 0; i < nBasis; i++) {
          for (int j = i; j < nBasis; j++) {
            std::tuple<int, int, int> calc = std::make_tuple(mode, i, j);
            calcsToDo.push(calc);
          }
        }
      }
      // Preliminary operator
      std::vector<bool> finished(size, false);
      finished[0] = true;
      int threadToSendCalc = -1;
      // Actual work distribution - loop until all threads have receive
      // termination signal
      while (!std::all_of(finished.begin(), finished.end(), [](bool i) { return i; })) {
        MPI_Status status;
        MPI_Recv(&threadToSendCalc, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        if (threadToSendCalc < 0) {
          finished[status.MPI_SOURCE] = true;
        } else {
          assert(threadToSendCalc == status.MPI_SOURCE);
          std::tuple<int, int, int> calc;
          if (!calcsToDo.empty()) {
            calc = calcsToDo.front();
            calcsToDo.pop();
          } else {
            calc = std::make_tuple(-1, -1, -1);
          }
          MPI_Send(&calc, 3, MPI_INT, threadToSendCalc, 1, MPI_COMM_WORLD);
        }
      }
    } else { // all other threads are doing calculations until they receive the
             // termination code
      bool doCalcs = true;
      while (doCalcs) {
        MPI_Send(&id, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        std::tuple<int, int, int> calc;
        MPI_Recv(&calc, 3, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        int mode, i, j;
        std::tie(mode, i, j) = calc;
        if (mode < 0) { // Termination code! Work of that thread is done
          int done = -1;
          MPI_Send(&done, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
          doCalcs = false;
          break;
        }
        double onebody = basis_->getOneBodyIntegrals(i, j, mode);
        double overlap = basis_->getOverlap(i, j, mode);

        oneBodyIntegralsL[mode](i, j) = onebody;
        oneBodyIntegralsL[mode](j, i) = onebody;

        overlapsL[mode](i, j) = overlap;
        overlapsL[mode](j, i) = overlap;
      }
    } // Now all calcs should be done
    for (int mode = 0; mode < nModes; mode++) {
      int nBasis = basis_->getBasisSize(mode);
      MPI_Allreduce(oneBodyIntegralsL[mode].data(), oneBodyIntegrals_[mode].data(), nBasis * nBasis, MPI_DOUBLE,
                    MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(overlapsL[mode].data(), overlaps_[mode].data(), nBasis * nBasis, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  #ifdef OTF
    // Only do this part if we are doing an on-the-fly PES calculation for
    // higher-order MC
    if ((std::dynamic_pointer_cast<OTFPes>(basis_->pes_) || std::dynamic_pointer_cast<HybridPes>(basis_->pes_)) &&
        potentialOrder > 1) {
      struct entryCalc1BSP {
        int mode;
        int disp;
        double res;
      };
      MPI_Aint mpidisps[3] = {offsetof(entryCalc1BSP, mode), offsetof(entryCalc1BSP, disp), offsetof(entryCalc1BSP, res)};
      MPI_Datatype entryCalc1BSPtype;
      int mpilengths[3] = {1, 1, 1};
      MPI_Datatype mpitypes[3] = {MPI_INT, MPI_INT, MPI_DOUBLE};
      MPI_Type_create_struct(3, mpilengths, mpidisps, mpitypes, &entryCalc1BSPtype);
      MPI_Type_commit(&entryCalc1BSPtype);
      calculatedSPs<1> calc1BSPsOfSendThread;
      if (std::dynamic_pointer_cast<OTFPes>(basis_->pes_)) {
        calc1BSPsOfSendThread = std::dynamic_pointer_cast<OTFPes>(basis_->pes_)->calc1BSPs_;
      } else if (std::dynamic_pointer_cast<HybridPes>(basis_->pes_)) {
        calc1BSPsOfSendThread = std::dynamic_pointer_cast<HybridPes>(basis_->pes_)->calc1BSPs_;
      }
      int numSPsOnThread = calc1BSPsOfSendThread.size();
      std::vector<std::pair<int, int>> keyVec;
      keyVec.reserve(numSPsOnThread);
      std::vector<double> valVec;
      valVec.reserve(numSPsOnThread);
      for (calculatedSPs<1>::iterator iter = calc1BSPsOfSendThread.begin(); iter != calc1BSPsOfSendThread.end(); ++iter) {
        keyVec.push_back(iter->first[0]);
        valVec.push_back(iter->second);
      }
      for (int send = 0; send < size; send++) {
        int numSPstoBCast = 0;
        if (id == send) {
          numSPstoBCast = numSPsOnThread;
        }
        MPI_Bcast(&numSPstoBCast, 1, MPI_INT, send, MPI_COMM_WORLD);
        if (numSPstoBCast > 0) {
          for (int i = 0; i < numSPstoBCast; i++) {
            entryCalc1BSP entry;
            if (id == send) {
              entry.mode = keyVec[i].first;
              entry.disp = keyVec[i].second;
              entry.res = valVec[i];
            }
            MPI_Bcast(&entry, 1, entryCalc1BSPtype, send, MPI_COMM_WORLD);
            if (id != send) {
              std::array<std::pair<int, int>, 1> key = {std::make_pair(entry.mode, entry.disp)};
              if (std::dynamic_pointer_cast<OTFPes>(basis_->pes_)) {
                std::dynamic_pointer_cast<OTFPes>(basis_->pes_)->calc1BSPs_[key] = entry.res;
              } else if (std::dynamic_pointer_cast<HybridPes>(basis_->pes_)) {
                std::dynamic_pointer_cast<HybridPes>(basis_->pes_)->calc1BSPs_[key] = entry.res;
              }
            }
            MPI_Barrier(MPI_COMM_WORLD);
          }
        }
      }
    }
  #endif
  }
  double after1B = MPI_Wtime();
  if (id == 0)
    std::cout << "Total 1B part took " << after1B - start1B << " seconds." << std::endl;
}
#endif

#ifndef MPI_PARALLEL
void VSCF::storeTwoBodyIntegrals() {
  #pragma omp parallel for collapse(2)
  for (int modei = 0; modei < nModes; modei++) {
    for (int modej = 0; modej < nModes; modej++) {
      if (std::find(coupledModes_.begin(), coupledModes_.end(), modei) != coupledModes_.end()) {
        if (std::find(coupledModes_.begin(), coupledModes_.end(), modej) != coupledModes_.end()) {
          if (modej != modei) {
            int nBasisi = basis_->getBasisSize(modei);
            int nBasisj = basis_->getBasisSize(modej);
            Eigen::SparseMatrix<double> tmpMat(nBasisi * nBasisi, nBasisj * nBasisj);
            tmpMat.reserve(nBasisi * nBasisi * nBasisj * nBasisj);
            for (int i = 0; i < nBasisi; i++) {
              for (int j = i; j < nBasisi; j++) {
                for (int k = 0; k < nBasisj; k++) {
                  for (int l = k; l < nBasisj; l++) {
                    double tijkl = basis_->getTwoBodyIntegrals(i, j, k, l, modei, modej);
                    if (std::abs(tijkl) > integralTol_) {
                      std::pair<int, int> inds = getTwoBodyIndex(i, j, k, l, modei, modej);
                      tmpMat.coeffRef(inds.first, inds.second) = tijkl;
                      std::pair<int, int> inds2 = getTwoBodyIndex(j, i, k, l, modei, modej);
                      tmpMat.coeffRef(inds2.first, inds2.second) = tijkl;
                      std::pair<int, int> inds3 = getTwoBodyIndex(j, i, l, k, modei, modej);
                      tmpMat.coeffRef(inds3.first, inds3.second) = tijkl;
                      std::pair<int, int> inds4 = getTwoBodyIndex(i, j, l, k, modei, modej);
                      tmpMat.coeffRef(inds4.first, inds4.second) = tijkl;
                    }
                  }
                }
              }
            }
  #pragma omp critical
            {
              tmpMat.makeCompressed();
              std::pair<int, int> ij = std::make_pair(modei, modej);
              twoBodyIntegrals_[ij] = tmpMat;
            }
          }
        }
      }
    }
  }
}
#else
void VSCF::storeTwoBodyIntegrals() {
  double start2B = MPI_Wtime();
  int size, id;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  // No need to do complicated task splitting for one-thread MPI job
  if (size == 1) {
    std::cout << "One-thread MPI job also for 2body" << std::endl;
    for (int modei = 0; modei < nModes; modei++) {
      if (std::find(coupledModes_.begin(), coupledModes_.end(), modei) != coupledModes_.end()) {
        for (int modej = 0; modej < nModes; modej++) {
          if (std::find(coupledModes_.begin(), coupledModes_.end(), modej) != coupledModes_.end() && modej != modei) {
            int nBasisi = basis_->getBasisSize(modei);
            int nBasisj = basis_->getBasisSize(modej);
            Eigen::SparseMatrix<double> tmpMat(nBasisi * nBasisi, nBasisj * nBasisj);
            tmpMat.reserve(nBasisi * nBasisi * nBasisj * nBasisj);
            for (int i = 0; i < nBasisi; i++) {
              for (int j = i; j < nBasisi; j++) {
                for (int k = 0; k < nBasisj; k++) {
                  for (int l = k; l < nBasisj; l++) {
                    double tijkl = basis_->getTwoBodyIntegrals(i, j, k, l, modei, modej);
                    if (std::abs(tijkl) > integralTol_) {
                      std::pair<int, int> inds = getTwoBodyIndex(i, j, k, l, modei, modej);
                      tmpMat.coeffRef(inds.first, inds.second) = tijkl;
                      std::pair<int, int> inds2 = getTwoBodyIndex(j, i, k, l, modei, modej);
                      tmpMat.coeffRef(inds2.first, inds2.second) = tijkl;
                      std::pair<int, int> inds3 = getTwoBodyIndex(j, i, l, k, modei, modej);
                      tmpMat.coeffRef(inds3.first, inds3.second) = tijkl;
                      std::pair<int, int> inds4 = getTwoBodyIndex(i, j, l, k, modei, modej);
                      tmpMat.coeffRef(inds4.first, inds4.second) = tijkl;
                    }
                  }
                }
              }
            }
            tmpMat.makeCompressed();
            std::pair<int, int> ij = std::make_pair(modei, modej);
            twoBodyIntegrals_[ij] = tmpMat;
          }
        }
      }
    }
  } else {
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    // Prepare so all threads have everything they need
    int numCalcs = 0;
    if (id == 0) { // Master thread is distributing work
      // Generate queue of tasks. Each calculation is specified by {calcIdx,
      // modei, modej, modali, modalj, modalk, modall} modali & modalj are of
      // modei, modalk & modall of modej
      std::queue<std::tuple<int, int, int, int, int, int, int>> calcsToDo;
      int calcIdx = 0;
      for (int modei = 0; modei < nModes; modei++) {
        if (std::find(coupledModes_.begin(), coupledModes_.end(), modei) != coupledModes_.end()) {
          for (int modej = modei + 1; modej < nModes; modej++) {
            if (std::find(coupledModes_.begin(), coupledModes_.end(), modej) != coupledModes_.end()) {
              int nBasisi = basis_->getBasisSize(modei);
              int nBasisj = basis_->getBasisSize(modej);
              for (int i = 0; i < nBasisi; i++) {
                for (int j = i; j < nBasisi; j++) {
                  for (int k = 0; k < nBasisj; k++) {
                    for (int l = k; l < nBasisj; l++) {
                      std::tuple<int, int, int, int, int, int, int> calc = std::make_tuple(calcIdx, modei, modej, i, j, k, l);
                      calcsToDo.push(calc);
                      ++calcIdx;
                    }
                  }
                }
              }
            }
          }
        }
      }
      numCalcs = calcIdx;
      // Preliminary operator
      std::vector<bool> finished(size, false);
      std::vector<int> numCalcsPerThread(size, 0);
      finished[0] = true;
      int threadToSendCalc = -1;
      // Actual work distribution - loop until all threads have receive
      // termination signal
      while (!std::all_of(finished.begin(), finished.end(), [](bool i) { return i; })) {
        MPI_Status status;
        MPI_Recv(&threadToSendCalc, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        if (threadToSendCalc < 0) {
          finished[status.MPI_SOURCE] = true;
        } else {
          assert(threadToSendCalc == status.MPI_SOURCE);
          std::tuple<int, int, int, int, int, int, int> calc;
          if (!calcsToDo.empty()) {
            calc = calcsToDo.front();
            calcsToDo.pop();
            numCalcsPerThread[threadToSendCalc] += 1;
          } else {
            calc = std::make_tuple(-1, -1, -1, -1, -1, -1, -1);
          }
          MPI_Send(&calc, 7, MPI_INT, threadToSendCalc, 1, MPI_COMM_WORLD);
        }
      }
      // Now all calcs should be done
      MPI_Barrier(MPI_COMM_WORLD);
      // Now retreive all the data -- loop over threads
      std::vector<int> allCalcIdxs;
      std::vector<double> allCalcResults;

      for (int worker = 1; worker < size; worker++) {
        int numCalcsOfWorker = numCalcsPerThread[worker];
        std::vector<int> calcIdxs;
        std::vector<double> calcResults;
        calcIdxs.resize(numCalcsOfWorker);
        calcResults.resize(numCalcsOfWorker);
        MPI_Recv(&calcIdxs[0], numCalcsOfWorker, MPI_INT, worker, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&calcResults[0], numCalcsOfWorker, MPI_DOUBLE, worker, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        allCalcIdxs.insert(allCalcIdxs.end(), calcIdxs.begin(), calcIdxs.end());
        allCalcResults.insert(allCalcResults.end(), calcResults.begin(), calcResults.end());
      }

      MPI_Barrier(MPI_COMM_WORLD);

      std::vector<int> indexOrdering(allCalcIdxs.size());
      std::iota(indexOrdering.begin(), indexOrdering.end(), 0);
      std::sort(indexOrdering.begin(), indexOrdering.end(),
                [&](int A, int B) -> bool { return allCalcIdxs[A] < allCalcIdxs[B]; });

      calcIdx = 0;
      for (int modei = 0; modei < nModes; modei++) {
        if (std::find(coupledModes_.begin(), coupledModes_.end(), modei) != coupledModes_.end()) {
          for (int modej = modei + 1; modej < nModes; modej++) {
            if (std::find(coupledModes_.begin(), coupledModes_.end(), modej) != coupledModes_.end()) {
              int nBasisi = basis_->getBasisSize(modei);
              int nBasisj = basis_->getBasisSize(modej);
              Eigen::SparseMatrix<double> tmpMat(nBasisi * nBasisi, nBasisj * nBasisj);
              tmpMat.reserve(nBasisi * nBasisi * nBasisj * nBasisj);
              for (int i = 0; i < nBasisi; i++) {
                for (int j = i; j < nBasisi; j++) {
                  for (int k = 0; k < nBasisj; k++) {
                    for (int l = k; l < nBasisj; l++) {
                      double tijkl = allCalcResults[indexOrdering[calcIdx]];
                      if (std::abs(tijkl) > integralTol_) {
                        std::pair<int, int> inds = getTwoBodyIndex(i, j, k, l, modei, modej);
                        tmpMat.coeffRef(inds.first, inds.second) = tijkl;
                        std::pair<int, int> inds2 = getTwoBodyIndex(j, i, k, l, modei, modej);
                        tmpMat.coeffRef(inds2.first, inds2.second) = tijkl;
                        std::pair<int, int> inds3 = getTwoBodyIndex(j, i, l, k, modei, modej);
                        tmpMat.coeffRef(inds3.first, inds3.second) = tijkl;
                        std::pair<int, int> inds4 = getTwoBodyIndex(i, j, l, k, modei, modej);
                        tmpMat.coeffRef(inds4.first, inds4.second) = tijkl;
                      }
                      ++calcIdx;
                    }
                  }
                }
              }
              Eigen::SparseMatrix<double> tmpMatTrans = tmpMat.transpose();
              tmpMat.makeCompressed();
              tmpMatTrans.makeCompressed();
              std::pair<int, int> ij = std::make_pair(modei, modej);
              twoBodyIntegrals_[ij] = tmpMat;
              std::pair<int, int> ji = std::make_pair(modej, modei);
              twoBodyIntegrals_[ji] = tmpMatTrans;
            }
          }
        }
      }
    } else { // all other threads are doing calculations until they receive the
             // termination code
      std::vector<int> calcIdxs;
      std::vector<double> calcResults;
      bool doCalcs = true;
      while (doCalcs) {
        MPI_Send(&id, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        std::tuple<int, int, int, int, int, int, int> calc;
        MPI_Recv(&calc, 7, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        int calcIdx, modei, modej, i, j, k, l;
        std::tie(calcIdx, modei, modej, i, j, k, l) = calc;
        if (calcIdx < 0) { // Termination code! Work of that thread is done
          int done = -1;
          MPI_Send(&done, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
          doCalcs = false;
          break;
        }
        double twobody = basis_->getTwoBodyIntegrals(i, j, k, l, modei, modej);
        calcIdxs.push_back(calcIdx);
        calcResults.push_back(twobody);
      }
      MPI_Barrier(MPI_COMM_WORLD);

      MPI_Send(&calcIdxs[0], calcIdxs.size(), MPI_INT, 0, 2, MPI_COMM_WORLD);
      MPI_Send(&calcResults[0], calcResults.size(), MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);

      MPI_Barrier(MPI_COMM_WORLD);

    } // Now all calcs should be done
    MPI_Barrier(MPI_COMM_WORLD);
  #ifdef OTF
    // Only do this part if we are doing an on-the-fly PES calculation for
    // higher-order MC
    if ((std::dynamic_pointer_cast<OTFPes>(basis_->pes_) || std::dynamic_pointer_cast<HybridPes>(basis_->pes_)) &&
        potentialOrder > 2) {
      struct entryCalc2BSP {
        int modes[2];
        int disps[2];
        double res;
      };
      MPI_Aint mpidisps[3] = {offsetof(entryCalc2BSP, modes), offsetof(entryCalc2BSP, disps), offsetof(entryCalc2BSP, res)};
      MPI_Datatype entryCalc2BSPtype;
      int mpilengths[3] = {2, 2, 1};
      MPI_Datatype mpitypes[3] = {MPI_INT, MPI_INT, MPI_DOUBLE};
      MPI_Type_create_struct(3, mpilengths, mpidisps, mpitypes, &entryCalc2BSPtype);
      MPI_Type_commit(&entryCalc2BSPtype);
      calculatedSPs<2> calc2BSPsOfSendThread;
      if (std::dynamic_pointer_cast<OTFPes>(basis_->pes_)) {
        calc2BSPsOfSendThread = std::dynamic_pointer_cast<OTFPes>(basis_->pes_)->calc2BSPs_;
      } else if (std::dynamic_pointer_cast<HybridPes>(basis_->pes_)) {
        calc2BSPsOfSendThread = std::dynamic_pointer_cast<HybridPes>(basis_->pes_)->calc2BSPs_;
      }
      int numSPsOnThread = calc2BSPsOfSendThread.size();
      std::vector<std::array<std::pair<int, int>, 2>> keyVec;
      keyVec.reserve(numSPsOnThread);
      std::vector<double> valVec;
      valVec.reserve(numSPsOnThread);
      for (calculatedSPs<2>::iterator iter = calc2BSPsOfSendThread.begin(); iter != calc2BSPsOfSendThread.end(); ++iter) {
        keyVec.push_back(iter->first);
        valVec.push_back(iter->second);
      }
      for (int send = 0; send < size; send++) {
        int numSPstoBCast = 0;
        if (id == send) {
          numSPstoBCast = numSPsOnThread;
        }
        MPI_Bcast(&numSPstoBCast, 1, MPI_INT, send, MPI_COMM_WORLD);
        if (numSPstoBCast > 0) {
          for (int i = 0; i < numSPstoBCast; i++) {
            entryCalc2BSP entry;
            if (id == send) {
              entry.modes[0] = keyVec[i][0].first;
              entry.modes[1] = keyVec[i][1].first;
              entry.disps[0] = keyVec[i][0].second;
              entry.disps[1] = keyVec[i][1].second;
              entry.res = valVec[i];
            }
            MPI_Bcast(&entry, 1, entryCalc2BSPtype, send, MPI_COMM_WORLD);
            if (id != send) {
              std::array<std::pair<int, int>, 2> key = {std::make_pair(entry.modes[0], entry.disps[0]),
                                                        std::make_pair(entry.modes[1], entry.disps[1])};
              if (std::dynamic_pointer_cast<OTFPes>(basis_->pes_)) {
                std::dynamic_pointer_cast<OTFPes>(basis_->pes_)->calc2BSPs_[key] = entry.res;
              } else if (std::dynamic_pointer_cast<HybridPes>(basis_->pes_)) {
                std::dynamic_pointer_cast<HybridPes>(basis_->pes_)->calc2BSPs_[key] = entry.res;
              }
            }
            MPI_Barrier(MPI_COMM_WORLD);
          }
        }
      }
    }
  #endif
  }
  double after2B = MPI_Wtime();
  if (id == 0)
    std::cout << "Total 2B part took " << after2B - start2B << " seconds." << std::endl;
}
#endif

// NG: think about this!
// this is currently probably the bottleneck, make it faster
// don't we store the same int 3! times here? Couldn't that loop be made more
// compact?
#ifndef MPI_PARALLEL
void VSCF::storeThreeBodyIntegrals() {
  #pragma omp parallel for collapse(3)
  for (int mode1 = 0; mode1 < nModes; mode1++) {
    for (int mode2 = 0; mode2 < nModes; mode2++) {
      for (int mode3 = 0; mode3 < nModes; mode3++) {
        if (std::find(coupledModes_.begin(), coupledModes_.end(), mode1) != coupledModes_.end()) {
          if (std::find(coupledModes_.begin(), coupledModes_.end(), mode2) != coupledModes_.end()) {
            if (std::find(coupledModes_.begin(), coupledModes_.end(), mode3) != coupledModes_.end()) {
              if ((mode2 != mode1) && (mode3 != mode1) && (mode3 != mode2)) {
                std::array<int, 3> key = {mode1, mode2, mode3};
                std::sort(key.begin(), key.end());
                auto find = set3B_.find(key);
  #ifdef OTF
                if (std::dynamic_pointer_cast<OTFPes>(basis_->pes_) || std::dynamic_pointer_cast<PesFromSPs>(basis_->pes_) ||
                    std::dynamic_pointer_cast<HybridPes>(basis_->pes_) || find != set3B_.end()) {
  #else
                if (std::dynamic_pointer_cast<PesFromSPs>(basis_->pes_) || find != set3B_.end()) {
  #endif
                  int nBasis1 = basis_->getBasisSize(mode1);
                  int nBasis2 = basis_->getBasisSize(mode2);
                  int nBasis3 = basis_->getBasisSize(mode3);
                  Eigen::SparseMatrix<double> tmpMat(nBasis1 * nBasis1, nBasis2 * nBasis2 * nBasis3 * nBasis3);
                  tmpMat.reserve(nBasis1 * nBasis1 * nBasis2 * nBasis2 * nBasis3 * nBasis3);
                  for (int i = 0; i < nBasis1; i++) {
                    for (int j = i; j < nBasis1; j++) {
                      for (int k = 0; k < nBasis2; k++) {
                        for (int l = k; l < nBasis2; l++) {
                          for (int m = 0; m < nBasis3; m++) {
                            for (int n = m; n < nBasis3; n++) {
                              double tijklmn = basis_->getThreeBodyIntegrals(i, j, k, l, m, n, mode1, mode2, mode3);
                              if (std::abs(tijklmn) > integralTol_) {
                                std::pair<int, int> inds = getThreeBodyIndex(i, j, k, l, m, n, mode1, mode2, mode3);
                                tmpMat.coeffRef(inds.first, inds.second) = tijklmn;
                                std::pair<int, int> inds2 = getThreeBodyIndex(i, j, k, l, n, m, mode1, mode2, mode3);
                                tmpMat.coeffRef(inds2.first, inds2.second) = tijklmn;
                                std::pair<int, int> inds3 = getThreeBodyIndex(i, j, l, k, m, n, mode1, mode2, mode3);
                                tmpMat.coeffRef(inds3.first, inds3.second) = tijklmn;
                                std::pair<int, int> inds4 = getThreeBodyIndex(i, j, l, k, n, m, mode1, mode2, mode3);
                                tmpMat.coeffRef(inds4.first, inds4.second) = tijklmn;
                                std::pair<int, int> inds5 = getThreeBodyIndex(j, i, k, l, m, n, mode1, mode2, mode3);
                                tmpMat.coeffRef(inds5.first, inds5.second) = tijklmn;
                                std::pair<int, int> inds6 = getThreeBodyIndex(j, i, k, l, n, m, mode1, mode2, mode3);
                                tmpMat.coeffRef(inds6.first, inds6.second) = tijklmn;
                                std::pair<int, int> inds7 = getThreeBodyIndex(j, i, l, k, m, n, mode1, mode2, mode3);
                                tmpMat.coeffRef(inds7.first, inds7.second) = tijklmn;
                                std::pair<int, int> inds8 = getThreeBodyIndex(j, i, l, k, n, m, mode1, mode2, mode3);
                                tmpMat.coeffRef(inds8.first, inds8.second) = tijklmn;
                                // probably there is a smarter way than storing
                                // this 2^3 * 3! times...
                              }
                            }
                          }
                        }
                      }
                    }
                  }
  #pragma omp critical
                  {
                    tmpMat.makeCompressed();
                    std::tuple<int, int, int> tup = std::make_tuple(mode1, mode2, mode3);
                    threeBodyIntegrals_[tup] = tmpMat;
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
#else
void VSCF::storeThreeBodyIntegrals() {
  double start3B = MPI_Wtime();
  int size, id;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  // No need to do complicated task splitting for one-thread MPI job
  if (size == 1) {
    std::cout << "One-thread MPI job also for 3body" << std::endl;
    for (int mode1 = 0; mode1 < nModes; mode1++) {
      if (std::find(coupledModes_.begin(), coupledModes_.end(), mode1) != coupledModes_.end()) {
        for (int mode2 = 0; mode2 < nModes; mode2++) {
          if (std::find(coupledModes_.begin(), coupledModes_.end(), mode2) != coupledModes_.end()) {
            for (int mode3 = 0; mode3 < nModes; mode3++) {
              if (std::find(coupledModes_.begin(), coupledModes_.end(), mode3) != coupledModes_.end()) {
                if ((mode2 != mode1) && (mode3 != mode1) && (mode3 != mode2)) {
                  std::array<int, 3> key = {mode1, mode2, mode3};
                  std::sort(key.begin(), key.end());
                  auto find = set3B_.find(key);
  #ifdef OTF
                  if (std::dynamic_pointer_cast<OTFPes>(basis_->pes_) || std::dynamic_pointer_cast<PesFromSPs>(basis_->pes_) ||
                      std::dynamic_pointer_cast<HybridPes>(basis_->pes_) || find != set3B_.end()) {
  #else
                  if (std::dynamic_pointer_cast<PesFromSPs>(basis_->pes_) || find != set3B_.end()) {
  #endif
                    int nBasis1 = basis_->getBasisSize(mode1);
                    int nBasis2 = basis_->getBasisSize(mode2);
                    int nBasis3 = basis_->getBasisSize(mode3);
                    Eigen::SparseMatrix<double> tmpMat(nBasis1 * nBasis1, nBasis2 * nBasis2 * nBasis3 * nBasis3);
                    tmpMat.reserve(nBasis1 * nBasis1 * nBasis2 * nBasis2 * nBasis3 * nBasis3);
                    for (int i = 0; i < nBasis1; i++) {
                      for (int j = i; j < nBasis1; j++) {
                        for (int k = 0; k < nBasis2; k++) {
                          for (int l = k; l < nBasis2; l++) {
                            for (int m = 0; m < nBasis3; m++) {
                              for (int n = m; n < nBasis3; n++) {
                                double tijklmn = basis_->getThreeBodyIntegrals(i, j, k, l, m, n, mode1, mode2, mode3);
                                if (std::abs(tijklmn) > integralTol_) {
                                  std::pair<int, int> inds = getThreeBodyIndex(i, j, k, l, m, n, mode1, mode2, mode3);
                                  tmpMat.coeffRef(inds.first, inds.second) = tijklmn;
                                  std::pair<int, int> inds2 = getThreeBodyIndex(i, j, k, l, n, m, mode1, mode2, mode3);
                                  tmpMat.coeffRef(inds2.first, inds2.second) = tijklmn;
                                  std::pair<int, int> inds3 = getThreeBodyIndex(i, j, l, k, m, n, mode1, mode2, mode3);
                                  tmpMat.coeffRef(inds3.first, inds3.second) = tijklmn;
                                  std::pair<int, int> inds4 = getThreeBodyIndex(i, j, l, k, n, m, mode1, mode2, mode3);
                                  tmpMat.coeffRef(inds4.first, inds4.second) = tijklmn;
                                  std::pair<int, int> inds5 = getThreeBodyIndex(j, i, k, l, m, n, mode1, mode2, mode3);
                                  tmpMat.coeffRef(inds5.first, inds5.second) = tijklmn;
                                  std::pair<int, int> inds6 = getThreeBodyIndex(j, i, k, l, n, m, mode1, mode2, mode3);
                                  tmpMat.coeffRef(inds6.first, inds6.second) = tijklmn;
                                  std::pair<int, int> inds7 = getThreeBodyIndex(j, i, l, k, m, n, mode1, mode2, mode3);
                                  tmpMat.coeffRef(inds7.first, inds7.second) = tijklmn;
                                  std::pair<int, int> inds8 = getThreeBodyIndex(j, i, l, k, n, m, mode1, mode2, mode3);
                                  tmpMat.coeffRef(inds8.first, inds8.second) = tijklmn;
                                  // probably there is a smarter way than
                                  // storing this 2^3 * 3! times...
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                    tmpMat.makeCompressed();
                    std::tuple<int, int, int> tup = std::make_tuple(mode1, mode2, mode3);
                    threeBodyIntegrals_[tup] = tmpMat;
                  }
                }
              }
            }
          }
        }
      }
    }
  } else {
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    // Prepare so all threads have everything they need
    int numCalcs = 0;
    if (id == 0) { // Master thread is distributing work
      // Generate queue of tasks. Each calculation is specified by {calcIdx,
      // mode1, mode2, mode3, modali, modalj, modalk, modall, modalm, modaln}
      // modali & modalj are of mode1, modalk & modall of mode2, modalm & modaln
      // of mode3
      std::queue<std::tuple<int, int, int, int, int, int, int, int, int, int>> calcsToDo;
      int calcIdx = 0;
      for (int mode1 = 0; mode1 < nModes; mode1++) {
        if (std::find(coupledModes_.begin(), coupledModes_.end(), mode1) != coupledModes_.end()) {
          for (int mode2 = mode1 + 1; mode2 < nModes; mode2++) {
            if (std::find(coupledModes_.begin(), coupledModes_.end(), mode2) != coupledModes_.end()) {
              for (int mode3 = mode2 + 1; mode3 < nModes; mode3++) {
                if (std::find(coupledModes_.begin(), coupledModes_.end(), mode3) != coupledModes_.end()) {
                  std::array<int, 3> key = {mode1, mode2, mode3};
                  std::sort(key.begin(), key.end());
                  auto find = set3B_.find(key);
  #ifdef OTF
                  if (std::dynamic_pointer_cast<OTFPes>(basis_->pes_) || std::dynamic_pointer_cast<PesFromSPs>(basis_->pes_) ||
                      std::dynamic_pointer_cast<HybridPes>(basis_->pes_) || find != set3B_.end()) {
  #else
                  if (std::dynamic_pointer_cast<PesFromSPs>(basis_->pes_) || find != set3B_.end()) {
  #endif
                    int nBasis1 = basis_->getBasisSize(mode1);
                    int nBasis2 = basis_->getBasisSize(mode2);
                    int nBasis3 = basis_->getBasisSize(mode3);
                    for (int i = 0; i < nBasis1; i++) {
                      for (int j = i; j < nBasis1; j++) {
                        for (int k = 0; k < nBasis2; k++) {
                          for (int l = k; l < nBasis2; l++) {
                            for (int m = 0; m < nBasis3; m++) {
                              for (int n = m; n < nBasis3; n++) {
                                std::tuple<int, int, int, int, int, int, int, int, int, int> calc =
                                    std::make_tuple(calcIdx, mode1, mode2, mode3, i, j, k, l, m, n);
                                calcsToDo.push(calc);
                                ++calcIdx;
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
      numCalcs = calcIdx;
      // Preliminary operator
      std::vector<bool> finished(size, false);
      std::vector<int> numCalcsPerThread(size, 0);
      finished[0] = true;
      int threadToSendCalc = -1;
      // Actual work distribution - loop until all threads have receive
      // termination signal
      while (!std::all_of(finished.begin(), finished.end(), [](bool i) { return i; })) {
        MPI_Status status;
        MPI_Recv(&threadToSendCalc, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        if (threadToSendCalc < 0) {
          finished[status.MPI_SOURCE] = true;
        } else {
          assert(threadToSendCalc == status.MPI_SOURCE);
          std::tuple<int, int, int, int, int, int, int, int, int, int> calc;
          if (!calcsToDo.empty()) {
            calc = calcsToDo.front();
            calcsToDo.pop();
            numCalcsPerThread[threadToSendCalc] += 1;
          } else {
            calc = std::make_tuple(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
          }
          MPI_Send(&calc, 10, MPI_INT, threadToSendCalc, 1, MPI_COMM_WORLD);
        }
      }
      // Now all calcs should be done
      MPI_Barrier(MPI_COMM_WORLD);
      // Now retreive all the data -- loop over threads
      std::vector<int> allCalcIdxs;
      std::vector<double> allCalcResults;

      for (int worker = 1; worker < size; worker++) {
        int numCalcsOfWorker = numCalcsPerThread[worker];
        std::vector<int> calcIdxs;
        std::vector<double> calcResults;
        calcIdxs.resize(numCalcsOfWorker);
        calcResults.resize(numCalcsOfWorker);
        MPI_Recv(&calcIdxs[0], numCalcsOfWorker, MPI_INT, worker, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&calcResults[0], numCalcsOfWorker, MPI_DOUBLE, worker, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        allCalcIdxs.insert(allCalcIdxs.end(), calcIdxs.begin(), calcIdxs.end());
        allCalcResults.insert(allCalcResults.end(), calcResults.begin(), calcResults.end());
      }

      MPI_Barrier(MPI_COMM_WORLD);

      std::vector<int> indexOrdering(allCalcIdxs.size());
      std::iota(indexOrdering.begin(), indexOrdering.end(), 0);
      std::sort(indexOrdering.begin(), indexOrdering.end(),
                [&](int A, int B) -> bool { return allCalcIdxs[A] < allCalcIdxs[B]; });

      calcIdx = 0;
      for (int mode1 = 0; mode1 < nModes; mode1++) {
        if (std::find(coupledModes_.begin(), coupledModes_.end(), mode1) != coupledModes_.end()) {
          for (int mode2 = mode1 + 1; mode2 < nModes; mode2++) {
            if (std::find(coupledModes_.begin(), coupledModes_.end(), mode2) != coupledModes_.end()) {
              for (int mode3 = mode2 + 1; mode3 < nModes; mode3++) {
                if (std::find(coupledModes_.begin(), coupledModes_.end(), mode3) != coupledModes_.end()) {
                  std::array<int, 3> key = {mode1, mode2, mode3};
                  std::sort(key.begin(), key.end());
                  auto find = set3B_.find(key);
  #ifdef OTF
                  if (std::dynamic_pointer_cast<OTFPes>(basis_->pes_) || std::dynamic_pointer_cast<PesFromSPs>(basis_->pes_) ||
                      std::dynamic_pointer_cast<HybridPes>(basis_->pes_) || find != set3B_.end()) {
  #else
                  if (std::dynamic_pointer_cast<PesFromSPs>(basis_->pes_) || find != set3B_.end()) {
  #endif
                    int nBasis1 = basis_->getBasisSize(mode1);
                    int nBasis2 = basis_->getBasisSize(mode2);
                    int nBasis3 = basis_->getBasisSize(mode3);
                    Eigen::SparseMatrix<double> tmpMat(nBasis1 * nBasis1, nBasis2 * nBasis2 * nBasis3 * nBasis3);
                    Eigen::SparseMatrix<double> tmpMat1(nBasis1 * nBasis1, nBasis2 * nBasis2 * nBasis3 * nBasis3);
                    Eigen::SparseMatrix<double> tmpMat2(nBasis2 * nBasis2, nBasis1 * nBasis1 * nBasis3 * nBasis3);
                    Eigen::SparseMatrix<double> tmpMat3(nBasis2 * nBasis2, nBasis1 * nBasis1 * nBasis3 * nBasis3);
                    Eigen::SparseMatrix<double> tmpMat4(nBasis3 * nBasis3, nBasis2 * nBasis2 * nBasis1 * nBasis1);
                    Eigen::SparseMatrix<double> tmpMat5(nBasis3 * nBasis3, nBasis2 * nBasis2 * nBasis1 * nBasis1);
                    tmpMat.reserve(nBasis1 * nBasis1 * nBasis2 * nBasis2 * nBasis3 * nBasis3);
                    tmpMat1.reserve(nBasis1 * nBasis1 * nBasis2 * nBasis2 * nBasis3 * nBasis3);
                    tmpMat2.reserve(nBasis1 * nBasis1 * nBasis2 * nBasis2 * nBasis3 * nBasis3);
                    tmpMat3.reserve(nBasis1 * nBasis1 * nBasis2 * nBasis2 * nBasis3 * nBasis3);
                    tmpMat4.reserve(nBasis1 * nBasis1 * nBasis2 * nBasis2 * nBasis3 * nBasis3);
                    tmpMat5.reserve(nBasis1 * nBasis1 * nBasis2 * nBasis2 * nBasis3 * nBasis3);
                    for (int i = 0; i < nBasis1; i++) {
                      for (int j = i; j < nBasis1; j++) {
                        for (int k = 0; k < nBasis2; k++) {
                          for (int l = k; l < nBasis2; l++) {
                            for (int m = 0; m < nBasis3; m++) {
                              for (int n = m; n < nBasis3; n++) {
                                double tijklmn = allCalcResults[indexOrdering[calcIdx]];
                                if (std::abs(tijklmn) > integralTol_) {
                                  // this is for mode1 - mode2 - mode3
                                  std::pair<int, int> inds1 = getThreeBodyIndex(i, j, k, l, m, n, mode1, mode2, mode3);
                                  tmpMat.coeffRef(inds1.first, inds1.second) = tijklmn;
                                  std::pair<int, int> inds2 = getThreeBodyIndex(i, j, k, l, n, m, mode1, mode2, mode3);
                                  tmpMat.coeffRef(inds2.first, inds2.second) = tijklmn;
                                  std::pair<int, int> inds3 = getThreeBodyIndex(i, j, l, k, m, n, mode1, mode2, mode3);
                                  tmpMat.coeffRef(inds3.first, inds3.second) = tijklmn;
                                  std::pair<int, int> inds4 = getThreeBodyIndex(i, j, l, k, n, m, mode1, mode2, mode3);
                                  tmpMat.coeffRef(inds4.first, inds4.second) = tijklmn;
                                  std::pair<int, int> inds5 = getThreeBodyIndex(j, i, k, l, m, n, mode1, mode2, mode3);
                                  tmpMat.coeffRef(inds5.first, inds5.second) = tijklmn;
                                  std::pair<int, int> inds6 = getThreeBodyIndex(j, i, k, l, n, m, mode1, mode2, mode3);
                                  tmpMat.coeffRef(inds6.first, inds6.second) = tijklmn;
                                  std::pair<int, int> inds7 = getThreeBodyIndex(j, i, l, k, m, n, mode1, mode2, mode3);
                                  tmpMat.coeffRef(inds7.first, inds7.second) = tijklmn;
                                  std::pair<int, int> inds8 = getThreeBodyIndex(j, i, l, k, n, m, mode1, mode2, mode3);
                                  tmpMat.coeffRef(inds8.first, inds8.second) = tijklmn;
                                  // this is for mode1 - mode3 - mode2
                                  std::pair<int, int> inds11 = getThreeBodyIndex(i, j, m, n, k, l, mode1, mode3, mode2);
                                  tmpMat1.coeffRef(inds11.first, inds11.second) = tijklmn;
                                  std::pair<int, int> inds12 = getThreeBodyIndex(i, j, n, m, k, l, mode1, mode3, mode2);
                                  tmpMat1.coeffRef(inds12.first, inds12.second) = tijklmn;
                                  std::pair<int, int> inds13 = getThreeBodyIndex(i, j, m, n, l, k, mode1, mode3, mode2);
                                  tmpMat1.coeffRef(inds13.first, inds13.second) = tijklmn;
                                  std::pair<int, int> inds14 = getThreeBodyIndex(i, j, n, m, l, k, mode1, mode3, mode2);
                                  tmpMat1.coeffRef(inds14.first, inds14.second) = tijklmn;
                                  std::pair<int, int> inds15 = getThreeBodyIndex(j, i, m, n, k, l, mode1, mode3, mode2);
                                  tmpMat1.coeffRef(inds15.first, inds15.second) = tijklmn;
                                  std::pair<int, int> inds16 = getThreeBodyIndex(j, i, n, m, k, l, mode1, mode3, mode2);
                                  tmpMat1.coeffRef(inds16.first, inds16.second) = tijklmn;
                                  std::pair<int, int> inds17 = getThreeBodyIndex(j, i, m, n, l, k, mode1, mode3, mode2);
                                  tmpMat1.coeffRef(inds17.first, inds17.second) = tijklmn;
                                  std::pair<int, int> inds18 = getThreeBodyIndex(j, i, n, m, l, k, mode1, mode3, mode2);
                                  tmpMat1.coeffRef(inds18.first, inds18.second) = tijklmn;
                                  // this is for mode2 - mode1 - mode3
                                  std::pair<int, int> inds21 = getThreeBodyIndex(k, l, i, j, m, n, mode2, mode1, mode3);
                                  tmpMat2.coeffRef(inds21.first, inds21.second) = tijklmn;
                                  std::pair<int, int> inds22 = getThreeBodyIndex(k, l, i, j, n, m, mode2, mode1, mode3);
                                  tmpMat2.coeffRef(inds22.first, inds22.second) = tijklmn;
                                  std::pair<int, int> inds23 = getThreeBodyIndex(l, k, i, j, m, n, mode2, mode1, mode3);
                                  tmpMat2.coeffRef(inds23.first, inds23.second) = tijklmn;
                                  std::pair<int, int> inds24 = getThreeBodyIndex(l, k, i, j, n, m, mode2, mode1, mode3);
                                  tmpMat2.coeffRef(inds24.first, inds24.second) = tijklmn;
                                  std::pair<int, int> inds25 = getThreeBodyIndex(k, l, j, i, m, n, mode2, mode1, mode3);
                                  tmpMat2.coeffRef(inds25.first, inds25.second) = tijklmn;
                                  std::pair<int, int> inds26 = getThreeBodyIndex(k, l, j, i, n, m, mode2, mode1, mode3);
                                  tmpMat2.coeffRef(inds26.first, inds26.second) = tijklmn;
                                  std::pair<int, int> inds27 = getThreeBodyIndex(l, k, j, i, m, n, mode2, mode1, mode3);
                                  tmpMat2.coeffRef(inds27.first, inds27.second) = tijklmn;
                                  std::pair<int, int> inds28 = getThreeBodyIndex(l, k, j, i, n, m, mode2, mode1, mode3);
                                  tmpMat2.coeffRef(inds28.first, inds28.second) = tijklmn;
                                  // this is for mode2 - mode3 - mode1
                                  std::pair<int, int> inds31 = getThreeBodyIndex(k, l, m, n, i, j, mode2, mode3, mode1);
                                  tmpMat3.coeffRef(inds31.first, inds31.second) = tijklmn;
                                  std::pair<int, int> inds32 = getThreeBodyIndex(k, l, n, m, i, j, mode2, mode3, mode1);
                                  tmpMat3.coeffRef(inds32.first, inds32.second) = tijklmn;
                                  std::pair<int, int> inds33 = getThreeBodyIndex(l, k, m, n, i, j, mode2, mode3, mode1);
                                  tmpMat3.coeffRef(inds33.first, inds33.second) = tijklmn;
                                  std::pair<int, int> inds34 = getThreeBodyIndex(l, k, n, m, i, j, mode2, mode3, mode1);
                                  tmpMat3.coeffRef(inds34.first, inds34.second) = tijklmn;
                                  std::pair<int, int> inds35 = getThreeBodyIndex(k, l, m, n, j, i, mode2, mode3, mode1);
                                  tmpMat3.coeffRef(inds35.first, inds35.second) = tijklmn;
                                  std::pair<int, int> inds36 = getThreeBodyIndex(k, l, n, m, j, i, mode2, mode3, mode1);
                                  tmpMat3.coeffRef(inds36.first, inds36.second) = tijklmn;
                                  std::pair<int, int> inds37 = getThreeBodyIndex(l, k, m, n, j, i, mode2, mode3, mode1);
                                  tmpMat3.coeffRef(inds37.first, inds37.second) = tijklmn;
                                  std::pair<int, int> inds38 = getThreeBodyIndex(l, k, n, m, j, i, mode2, mode3, mode1);
                                  tmpMat3.coeffRef(inds38.first, inds38.second) = tijklmn;
                                  // this is for mode3 - mode1 - mode2
                                  std::pair<int, int> inds41 = getThreeBodyIndex(m, n, i, j, k, l, mode3, mode1, mode2);
                                  tmpMat4.coeffRef(inds41.first, inds41.second) = tijklmn;
                                  std::pair<int, int> inds42 = getThreeBodyIndex(n, m, i, j, k, l, mode3, mode1, mode2);
                                  tmpMat4.coeffRef(inds42.first, inds42.second) = tijklmn;
                                  std::pair<int, int> inds43 = getThreeBodyIndex(m, n, i, j, l, k, mode3, mode1, mode2);
                                  tmpMat4.coeffRef(inds43.first, inds43.second) = tijklmn;
                                  std::pair<int, int> inds44 = getThreeBodyIndex(n, m, i, j, l, k, mode3, mode1, mode2);
                                  tmpMat4.coeffRef(inds44.first, inds44.second) = tijklmn;
                                  std::pair<int, int> inds45 = getThreeBodyIndex(m, n, j, i, k, l, mode3, mode1, mode2);
                                  tmpMat4.coeffRef(inds45.first, inds45.second) = tijklmn;
                                  std::pair<int, int> inds46 = getThreeBodyIndex(n, m, j, i, k, l, mode3, mode1, mode2);
                                  tmpMat4.coeffRef(inds46.first, inds46.second) = tijklmn;
                                  std::pair<int, int> inds47 = getThreeBodyIndex(m, n, j, i, l, k, mode3, mode1, mode2);
                                  tmpMat4.coeffRef(inds47.first, inds47.second) = tijklmn;
                                  std::pair<int, int> inds48 = getThreeBodyIndex(n, m, j, i, l, k, mode3, mode1, mode2);
                                  tmpMat4.coeffRef(inds48.first, inds48.second) = tijklmn;
                                  // this is for mode3 - mode2 - mode1
                                  std::pair<int, int> inds51 = getThreeBodyIndex(m, n, k, l, i, j, mode3, mode2, mode1);
                                  tmpMat5.coeffRef(inds51.first, inds51.second) = tijklmn;
                                  std::pair<int, int> inds52 = getThreeBodyIndex(n, m, k, l, i, j, mode3, mode2, mode1);
                                  tmpMat5.coeffRef(inds52.first, inds52.second) = tijklmn;
                                  std::pair<int, int> inds53 = getThreeBodyIndex(m, n, l, k, i, j, mode3, mode2, mode1);
                                  tmpMat5.coeffRef(inds53.first, inds53.second) = tijklmn;
                                  std::pair<int, int> inds54 = getThreeBodyIndex(n, m, l, k, i, j, mode3, mode2, mode1);
                                  tmpMat5.coeffRef(inds54.first, inds54.second) = tijklmn;
                                  std::pair<int, int> inds55 = getThreeBodyIndex(m, n, k, l, j, i, mode3, mode2, mode1);
                                  tmpMat5.coeffRef(inds55.first, inds55.second) = tijklmn;
                                  std::pair<int, int> inds56 = getThreeBodyIndex(n, m, k, l, j, i, mode3, mode2, mode1);
                                  tmpMat5.coeffRef(inds56.first, inds56.second) = tijklmn;
                                  std::pair<int, int> inds57 = getThreeBodyIndex(m, n, l, k, j, i, mode3, mode2, mode1);
                                  tmpMat5.coeffRef(inds57.first, inds57.second) = tijklmn;
                                  std::pair<int, int> inds58 = getThreeBodyIndex(n, m, l, k, j, i, mode3, mode2, mode1);
                                  tmpMat5.coeffRef(inds58.first, inds58.second) = tijklmn;
                                  // probably there is a smarter way than
                                  // storing this 2^3 * 3! times...
                                }
                                ++calcIdx;
                              }
                            }
                          }
                        }
                      }
                      tmpMat.makeCompressed();
                      tmpMat1.makeCompressed();
                      tmpMat2.makeCompressed();
                      tmpMat3.makeCompressed();
                      tmpMat4.makeCompressed();
                      tmpMat5.makeCompressed();
                      std::tuple<int, int, int> tup = std::make_tuple(mode1, mode2, mode3);
                      std::tuple<int, int, int> tup1 = std::make_tuple(mode1, mode3, mode2);
                      std::tuple<int, int, int> tup2 = std::make_tuple(mode2, mode1, mode3);
                      std::tuple<int, int, int> tup3 = std::make_tuple(mode2, mode3, mode1);
                      std::tuple<int, int, int> tup4 = std::make_tuple(mode3, mode1, mode2);
                      std::tuple<int, int, int> tup5 = std::make_tuple(mode3, mode2, mode1);
                      threeBodyIntegrals_[tup] = tmpMat;
                      threeBodyIntegrals_[tup1] = tmpMat1;
                      threeBodyIntegrals_[tup2] = tmpMat2;
                      threeBodyIntegrals_[tup3] = tmpMat3;
                      threeBodyIntegrals_[tup4] = tmpMat4;
                      threeBodyIntegrals_[tup5] = tmpMat5;
                    }
                  }
                }
              }
            }
          }
        }
      }
    } else { // all other threads are doing calculations until they receive the
             // termination code
      std::vector<int> calcIdxs;
      std::vector<double> calcResults;
      bool doCalcs = true;
      while (doCalcs) {
        MPI_Send(&id, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        std::tuple<int, int, int, int, int, int, int, int, int, int> calc;
        MPI_Recv(&calc, 10, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        int calcIdx, mode1, mode2, mode3, i, j, k, l, m, n;
        std::tie(calcIdx, mode1, mode2, mode3, i, j, k, l, m, n) = calc;
        if (calcIdx < 0) { // Termination code! Work of that thread is done
          int done = -1;
          MPI_Send(&done, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
          doCalcs = false;
          break;
        }
        double threebody = basis_->getThreeBodyIntegrals(i, j, k, l, m, n, mode1, mode2, mode3);
        calcIdxs.push_back(calcIdx);
        calcResults.push_back(threebody);
      }
      MPI_Barrier(MPI_COMM_WORLD);

      MPI_Send(&calcIdxs[0], calcIdxs.size(), MPI_INT, 0, 2, MPI_COMM_WORLD);
      MPI_Send(&calcResults[0], calcResults.size(), MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);

      MPI_Barrier(MPI_COMM_WORLD);

    } // Now all calcs should be done
    MPI_Barrier(MPI_COMM_WORLD);
  }
  double after3B = MPI_Wtime();
  if (id == 0)
    std::cout << "Total 3B part took " << after3B - start3B << " seconds." << std::endl;
}
#endif

Eigen::MatrixXd VSCF::getOverlap(int mode) {
  return overlaps_[mode];
}

int VSCF::getBasisSize(int mode) {
  return basis_->getBasisSize(mode);
}

std::pair<int, int> VSCF::getTwoBodyIndex(int i, int j, int k, int l, int modei, int modej) const {
  int Nbasisj = basis_->getBasisSize(modej);
  int Nbasisi = basis_->getBasisSize(modei);
  return std::make_pair(i * Nbasisi + j, k * Nbasisj + l);
}

std::pair<int, int> VSCF::getThreeBodyIndex(int i, int j, int k, int l, int m, int n, int mode1, int mode2, int mode3) const {
  int Nbasis1 = basis_->getBasisSize(mode1);
  int Nbasis2 = basis_->getBasisSize(mode2);
  int Nbasis3 = basis_->getBasisSize(mode3);
  return std::make_pair(i * Nbasis1 + j, (k * Nbasis2 + l) * (Nbasis3 * Nbasis3) + (m * Nbasis3) + n);
}

int VSCF::getDensityIndex(int k, int l, int modej) const {
  int nBasisj = basis_->getBasisSize(modej);
  return (k * nBasisj + l);
}

Eigen::MatrixXd VSCF::getMeanFieldOperator(int mode) {
  return oneBodyIntegrals_[mode];
}

Eigen::MatrixXd VSCF::getMeanFieldOperator(int mode, const std::vector<Eigen::VectorXd>& coefficients) {
  return (oneBodyIntegrals_[mode] + contractTwoBodyIntegral(mode, coefficients));
}

Eigen::MatrixXd VSCF::getMeanFieldOperator3(int mode, const std::vector<Eigen::VectorXd>& coefficients) {
  return (oneBodyIntegrals_[mode] + contractTwoBodyIntegral(mode, coefficients) + contractThreeBodyIntegral(mode, coefficients));
}

Eigen::MatrixXd VSCF::getKineticOperator(int mode) {
  int nBasis = basis_->getBasisSize(mode);
  Eigen::MatrixXd onebodyKin = Eigen::MatrixXd::Zero(nBasis, nBasis);
  for (int i = 0; i < nBasis; i++) {
    for (int j = 0; j < nBasis; j++) {
      if (i >= j) {
        onebodyKin(i, j) = basis_->getKinetic(i, j, mode);
      }
    }
  }
  Eigen::MatrixXd transp = onebodyKin.triangularView<Eigen::StrictlyLower>();
  transp.transposeInPlace();
  return (onebodyKin + transp);
}

double VSCF::getEnergy(const std::vector<Eigen::VectorXd>& coeffs, const std::vector<double>& evals) {
  double energy = 0;
  for (int mode = 0; mode < nModes; mode++) {
    int nBasis = basis_->getBasisSize(mode);
    Eigen::MatrixXd tmpcorr = contractTwoBodyIntegral(mode, coeffs);
    energy += evals[mode] - 1. / 2 * coeffs[mode].transpose() * tmpcorr * coeffs[mode];
  }
  return energy;
}

double VSCF::getEnergy3(const std::vector<Eigen::VectorXd>& coeffs, const std::vector<double>& evals) {
  double energy = 0;
  for (int mode = 0; mode < nModes; mode++) {
    int nBasis = basis_->getBasisSize(mode);
    Eigen::MatrixXd tmpcorr = contractTwoBodyIntegral(mode, coeffs);
    Eigen::MatrixXd tmpcorr3 = contractThreeBodyIntegral(mode, coeffs);
    energy += evals[mode] - 1. / 2 * coeffs[mode].transpose() * tmpcorr * coeffs[mode] -
              5. / 6 * coeffs[mode].transpose() * tmpcorr3 * coeffs[mode];
  }
  return energy;
}

double VSCF::getEnergy(const std::vector<double>& evals) {
  double energy = 0;
  for (int mode = 0; mode < nModes; mode++) {
    energy += evals[mode];
  }
  return energy;
}

void VSCF::printAllEnergies(const std::vector<Eigen::VectorXd>& coeffs, std::vector<Eigen::VectorXd>& allEvals) {
  double energy = 0;
  for (int mode = 0; mode < nModes; mode++) {
    Eigen::VectorXd evals = allEvals[mode];
    std::sort(evals.data(), evals.data() + evals.size());
    std::cout << "Mode " << mode << " has evals " << Utils::Constants::invCentimeter_per_hartree * evals.transpose()
              << std::endl;
    Eigen::MatrixXd tmpcorr = contractTwoBodyIntegral(mode, coeffs);
    Eigen::MatrixXd tmpcorr3 = contractThreeBodyIntegral(mode, coeffs);
    int conf = occupation[mode];
    energy += evals[conf] - 1. / 2 * coeffs[mode].transpose() * tmpcorr * coeffs[mode] -
              5. / 6 * coeffs[mode].transpose() * tmpcorr3 * coeffs[mode];
  }
  std::cout << "Energy computed for all modes: " << Utils::Constants::invCentimeter_per_hartree * energy << std::endl;
  double energyC = 0;
  for (int mode = 0; mode < nModes; mode++) {
    if (std::find(coupledModes_.begin(), coupledModes_.end(), mode) != coupledModes_.end()) {
      Eigen::VectorXd evals = allEvals[mode];
      std::sort(evals.data(), evals.data() + evals.size());
      Eigen::MatrixXd tmpcorr = contractTwoBodyIntegral(mode, coeffs);
      Eigen::MatrixXd tmpcorr3 = contractThreeBodyIntegral(mode, coeffs);
      int conf = occupation[mode];
      energyC += evals[conf] - 1. / 2 * coeffs[mode].transpose() * tmpcorr * coeffs[mode] -
                 5. / 6 * coeffs[mode].transpose() * tmpcorr3 * coeffs[mode];
    }
  }
  std::cout << "Energy computed for dumpOnlyCoupledModes: " << Utils::Constants::invCentimeter_per_hartree * energyC
            << std::endl;
}

Eigen::VectorXd VSCF::ConstructDensity(int modej, const Eigen::VectorXd& coeffLeft, const Eigen::VectorXd& coeffRight) {
  Eigen::MatrixXd densityMat = coeffLeft * coeffRight.transpose();
  int nBasisj = basis_->getBasisSize(modej);
  Eigen::VectorXd density(nBasisj * nBasisj);
  for (int k = 0; k < nBasisj; k++) {
    for (int l = 0; l < nBasisj; l++) {
      int a = getDensityIndex(k, l, modej);
      density(a) = densityMat(k, l);
    }
  }
  return density;
}

Eigen::VectorXd VSCF::ConstructDensity3(int mode2, const Eigen::VectorXd& coeffLeft2, const Eigen::VectorXd& coeffRight2,
                                        int mode3, const Eigen::VectorXd& coeffLeft3, const Eigen::VectorXd& coeffRight3) {
  Eigen::MatrixXd densityMat2 = coeffLeft2 * coeffRight2.transpose();
  Eigen::MatrixXd densityMat3 = coeffLeft3 * coeffRight3.transpose();
  int nBasis2 = basis_->getBasisSize(mode2);
  int nBasis3 = basis_->getBasisSize(mode3);
  Eigen::VectorXd density(nBasis2 * nBasis2 * nBasis3 * nBasis3);
  for (int k = 0; k < nBasis2; k++) {
    for (int l = 0; l < nBasis2; l++) {
      for (int m = 0; m < nBasis3; m++) {
        for (int n = 0; n < nBasis3; n++) {
          int a = getDensityIndex(k, l, mode2);
          int b = getDensityIndex(m, n, mode3);
          density(a * nBasis3 * nBasis3 + b) = densityMat2(k, l) * densityMat3(m, n);
        }
      }
    }
  }
  return density;
}

Eigen::MatrixXd VSCF::contractTwoBodyIntegral(int mode, const std::vector<Eigen::VectorXd>& coeff) {
  int nBasis = basis_->getBasisSize(mode);
  Eigen::MatrixXd twoBodyIntegral = Eigen::MatrixXd::Zero(nBasis, nBasis);
#pragma omp parallel for
  for (int modej = 0; modej < nModes; modej++) {
    if (modej != mode) {
      Eigen::VectorXd density = this->ConstructDensity(modej, coeff[modej], coeff[modej]);
      std::pair<int, int> ij = std::make_pair(mode, modej);
      if (twoBodyIntegrals_.find(ij) != twoBodyIntegrals_.end()) {
        Eigen::VectorXd tmpTwob = twoBodyIntegrals_[ij] * density;
        Eigen::MatrixXd tmpTwoB = Eigen::MatrixXd::Zero(nBasis, nBasis);
        for (int i = 0; i < nBasis; i++) {
          for (int j = 0; j < nBasis; j++) {
            int a = getDensityIndex(i, j, mode);
            tmpTwoB(i, j) = tmpTwob(a);
          }
        }
#pragma omp critical
        { twoBodyIntegral += tmpTwoB; }
      }
    }
  }
  return twoBodyIntegral;
}

Eigen::MatrixXd VSCF::contractThreeBodyIntegral(int mode, const std::vector<Eigen::VectorXd>& coeff) {
  int nBasis = basis_->getBasisSize(mode);
  Eigen::MatrixXd threeBodyIntegral = Eigen::MatrixXd::Zero(nBasis, nBasis);
#pragma omp parallel for collapse(2)
  for (int modei = 0; modei < nModes; modei++) {
    for (int modej = 0; modej < nModes; modej++) {
      if ((modei != mode) && (modej != mode) && (modej != modei)) {
        std::array<int, 3> key = {mode, modei, modej};
        std::sort(key.begin(), key.end());
        auto find = set3B_.find(key);
        std::tuple<int, int, int> tup = std::make_tuple(mode, modei, modej);
#ifdef OTF
        if (((std::dynamic_pointer_cast<OTFPes>(basis_->pes_) || std::dynamic_pointer_cast<HybridPes>(basis_->pes_) ||
              std::dynamic_pointer_cast<PesFromSPs>(basis_->pes_)) &&
             (threeBodyIntegrals_.find(tup) != threeBodyIntegrals_.end())) ||
            find != set3B_.end()) {
#else
        if ((std::dynamic_pointer_cast<PesFromSPs>(basis_->pes_) &&
             (threeBodyIntegrals_.find(tup) != threeBodyIntegrals_.end())) ||
            find != set3B_.end()) {
#endif
          Eigen::VectorXd density =
              this->ConstructDensity3(modei, coeff[modei], coeff[modei], modej, coeff[modej], coeff[modej]);
          Eigen::VectorXd tmpThreeb = threeBodyIntegrals_[tup] * density;
          Eigen::MatrixXd tmpThreeB = Eigen::MatrixXd::Zero(nBasis, nBasis);
          for (int i = 0; i < nBasis; i++) {
            for (int j = 0; j < nBasis; j++) {
              int a = getDensityIndex(i, j, mode);
              tmpThreeB(i, j) = tmpThreeb(a);
            }
          }
#pragma omp critical
          { threeBodyIntegral += tmpThreeB; }
        }
      }
    }
  }
  return threeBodyIntegral;
}

std::vector<Eigen::VectorXd> VSCF::getStartingGuess(bool randomize) {
  std::vector<Eigen::VectorXd> guesses(nModes);
  if (randomize) {
    // this does not seem to be properly implemented, as nothing is done with
    // the random number
    std::random_device rd;  // to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> distrib(0, 1);
    for (int mode = 0; mode < nModes; mode++) {
      Eigen::MatrixXd F = getMeanFieldOperator(mode);
      Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> guess;
      guess.compute(F, overlaps_[mode]);
      Eigen::MatrixXd vtemp = guess.eigenvectors();
      Eigen::VectorXd etemp = guess.eigenvalues();

      std::vector<std::pair<double, Eigen::VectorXd>> guessPairs(etemp.size());
      for (int i = 0; i < etemp.size(); i++) {
        guessPairs[i].first = etemp(i);
        guessPairs[i].second = vtemp.col(i);
      }

      std::sort(guessPairs.begin(), guessPairs.end(), [](auto a, auto b) { return a.first < b.first; });

      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> s_evals;
      s_evals.compute(overlaps_[mode]);
      Eigen::VectorXd stemp = s_evals.eigenvalues().cwiseInverse().cwiseSqrt();
      Eigen::MatrixXd eigenVectors = s_evals.eigenvectors();
      Eigen::MatrixXd X = eigenVectors * stemp.asDiagonal() * eigenVectors.transpose();
      Eigen::VectorXd tmp = guessPairs[0].second;
      for (int iElement = 1; iElement < guessPairs.size(); iElement++) {
        tmp += guessPairs[iElement].second * distrib(gen);
      }
      Eigen::VectorXd vecGuess = 1 / sqrt(guessPairs[0].second.transpose() * guessPairs[0].second) * guessPairs[0].second;
      guesses[mode] = X * vecGuess;
    }
  } else {
    for (int mode = 0; mode < nModes; mode++) {
      Eigen::MatrixXd F = getMeanFieldOperator(mode);
      Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> guess;
      guess.compute(F, overlaps_[mode]);
      Eigen::MatrixXd vtemp = guess.eigenvectors();
      Eigen::VectorXd etemp = guess.eigenvalues();

      std::vector<std::pair<double, Eigen::VectorXd>> guessPairs(etemp.size());
      for (int i = 0; i < etemp.size(); i++) {
        guessPairs[i].first = etemp(i);
        guessPairs[i].second = vtemp.col(i);
      }

      std::sort(guessPairs.begin(), guessPairs.end(), [](auto a, auto b) { return a.first < b.first; });

      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> s_evals;
      s_evals.compute(overlaps_[mode]);
      Eigen::VectorXd stemp = s_evals.eigenvalues().cwiseInverse().cwiseSqrt();
      Eigen::MatrixXd eigenVectors = s_evals.eigenvectors();
      Eigen::MatrixXd X = eigenVectors * stemp.asDiagonal() * eigenVectors.transpose();

      int conf = occupation[mode]; // for state-specific initial guess

      Eigen::VectorXd vecGuess =
          1 / sqrt(guessPairs[conf].second.transpose() * guessPairs[conf].second) * guessPairs[conf].second;

      guesses[mode] = X * vecGuess;
    }
  }
  return guesses;
}

void VSCF::dumpIntegrals(const std::vector<std::vector<std::pair<double, Eigen::VectorXd>>>& evalPairs) {
  for (int mode = 0; mode < nModes; mode++) {
    for (int modal1 = 0; modal1 < nMax_[mode]; modal1++) {
      for (int modal2 = 0; modal2 < nMax_[mode]; modal2++) {
        std::array<int, 3> inds = {mode, modal1, modal2};
        Eigen::VectorXd coeffModal1 = evalPairs[mode][modal1].second;
        Eigen::VectorXd coeffModal2 = evalPairs[mode][modal2].second;
        double resOneBody = coeffModal1.transpose() * oneBodyIntegrals_[mode] * coeffModal2;
        oneBodyModals_[inds] = resOneBody;
      }
    }
  }
  if (potentialOrder == 1) {
    IntegralDumper::writeIntegrals(oneBodyModals_, nMax_, nModes, integralDumpFile, dumpOnlyCoupledModes_, coupledModes_);
  } else if (potentialOrder == 2 || potentialOrder == 5 || potentialOrder == 6 || potentialOrder == 7) {
#pragma omp parallel for collapse(2)
    for (int modei = 0; modei < nModes; modei++) {
      for (int modej = 0; modej < nModes; modej++) {
        if (modej > modei) {
          for (int modali1 = 0; modali1 < nMax_[modei]; modali1++) {
            for (int modali2 = 0; modali2 < nMax_[modei]; modali2++) {
              for (int modalj1 = 0; modalj1 < nMax_[modej]; modalj1++) {
                for (int modalj2 = 0; modalj2 < nMax_[modej]; modalj2++) {
                  Eigen::VectorXd coeffModali1 = evalPairs[modei][modali1].second;
                  Eigen::VectorXd coeffModali2 = evalPairs[modei][modali2].second;
                  Eigen::VectorXd coeffModalj1 = evalPairs[modej][modalj1].second;
                  Eigen::VectorXd coeffModalj2 = evalPairs[modej][modalj2].second;
                  Eigen::VectorXd densityi = this->ConstructDensity(modei, coeffModali1, coeffModali2);
                  Eigen::VectorXd densityj = this->ConstructDensity(modej, coeffModalj1, coeffModalj2);
#pragma omp critical
                  {
                    if (twoBodyIntegrals_.find(std::make_pair(modei, modej)) != twoBodyIntegrals_.end()) {
                      Eigen::VectorXd tmpTwoB = twoBodyIntegrals_.at(std::make_pair(modei, modej)) * densityj;

                      double resTwoBody = densityi.transpose() * tmpTwoB;

                      twoBodyModals_[{modei, modali1, modali2, modej, modalj1, modalj2}] = resTwoBody;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    IntegralDumper::writeIntegrals(oneBodyModals_, twoBodyModals_, nMax_, nModes, integralDumpFile,
                                   dumpOnlyCoupledModes_, coupledModes_);
  } else if (potentialOrder == 3) {
#pragma omp parallel for collapse(2)
    for (int modei = 0; modei < nModes; modei++) {
      for (int modej = 0; modej < nModes; modej++) {
        if (modej > modei) {
          for (int modali1 = 0; modali1 < nMax_[modei]; modali1++) {
            for (int modali2 = 0; modali2 < nMax_[modei]; modali2++) {
              for (int modalj1 = 0; modalj1 < nMax_[modej]; modalj1++) {
                for (int modalj2 = 0; modalj2 < nMax_[modej]; modalj2++) {
                  Eigen::VectorXd coeffModali1 = evalPairs[modei][modali1].second;
                  Eigen::VectorXd coeffModali2 = evalPairs[modei][modali2].second;
                  Eigen::VectorXd coeffModalj1 = evalPairs[modej][modalj1].second;
                  Eigen::VectorXd coeffModalj2 = evalPairs[modej][modalj2].second;
                  Eigen::VectorXd densityi = this->ConstructDensity(modei, coeffModali1, coeffModali2);
                  Eigen::VectorXd densityj = this->ConstructDensity(modej, coeffModalj1, coeffModalj2);
#pragma omp critical
                  {
                    if (twoBodyIntegrals_.find(std::make_pair(modei, modej)) != twoBodyIntegrals_.end()) {
                      Eigen::VectorXd tmpTwoB = twoBodyIntegrals_.at(std::make_pair(modei, modej)) * densityj;
                      double resTwoBody = densityi.transpose() * tmpTwoB;

                      twoBodyModals_[{modei, modali1, modali2, modej, modalj1, modalj2}] = resTwoBody;
                    }
                  }
                  // we have three densities and the threeBodyIntegrals are a
                  // sparse matrix with special indexing
                  for (int modek = 0; modek < nModes; modek++) {
                    if (modek > modej) {
                      for (int modalk1 = 0; modalk1 < nMax_[modek]; modalk1++) {
                        for (int modalk2 = 0; modalk2 < nMax_[modek]; modalk2++) {
                          Eigen::VectorXd coeffModalk1 = evalPairs[modek][modalk1].second;
                          Eigen::VectorXd coeffModalk2 = evalPairs[modek][modalk2].second;
                          Eigen::VectorXd density3 =
                              this->ConstructDensity3(modej, coeffModalj1, coeffModalj2, modek, coeffModalk1, coeffModalk2);
#pragma omp critical
                          {
                            if (threeBodyIntegrals_.find(std::make_tuple(modei, modej, modek)) != threeBodyIntegrals_.end()) {
                              Eigen::VectorXd tmpThreeB =
                                  threeBodyIntegrals_.at(std::make_tuple(modei, modej, modek)) * density3;
                              double resThreeBody = densityi.transpose() * tmpThreeB;

                              threeBodyModals_[{modei, modali1, modali2, modej, modalj1, modalj2, modek, modalk1, modalk2}] =
                                  resThreeBody;
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
    IntegralDumper::writeIntegrals(oneBodyModals_, twoBodyModals_, threeBodyModals_, nMax_, nModes, integralDumpFile,
                                   dumpOnlyCoupledModes_, coupledModes_);
  }
}

tensor<3> VSCF::getOneBodyModals() {
  return oneBodyModals_;
}

tensor<6> VSCF::getTwoBodyModals() {
  return twoBodyModals_;
}

tensor<9> VSCF::getThreeBodyModals() {
  return threeBodyModals_;
}

} // namespace Scine::Colibri
