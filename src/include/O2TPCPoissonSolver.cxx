// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file O2TPCPoissonSolver.cxx
/// \brief This class provides implementation of Poisson Eq
/// solver by MultiGrid Method
///
///
///
/// \author  Matthias Kleiner <matthias.kleiner@cern.ch>
/// \date Aug 21, 2020

#include <TMath.h>
#include "O2TPCPoissonSolver.h"

#include <iostream>
#include <numeric>

using namespace o2::tpc;

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nz, Nr, Nphi>::poissonSolver3D(DataContainer& matricesV, const DataContainer& matricesCharge, const int symmetry)
{
  // poissonMultiGrid3D2D(matricesV, matricesCharge, symmetry);
  poissonMultiGrid3D(matricesV, matricesCharge, symmetry);
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nz, Nr, Nphi>::poissonMultiGrid3D2D(DataContainer& matricesV, const DataContainer& matricesCharge, const int symmetry)
{
  Info("PoissonMultiGrid3D2D", "%s", Form("in Poisson Solver 3D multiGrid semi coarsening Nr=%lu, cols=%lu, Nphi=%lu \n", Nz, Nr, Nphi));

  // Check that the number of Nr and Nz is suitable for a binary expansion
  if (!isPowerOfTwo((Nr - 1))) {
    Error("PoissonMultiGrid3D2D", "Poisson3DMultiGrid - Error in the number of Nr. Must be 2**M + 1");
    return;
  }
  if (!isPowerOfTwo((Nz - 1))) {
    Error("PoissonMultiGrid3D2D", "Poisson3DMultiGrid - Error in the number of Nz. Must be 2**N - 1");
    return;
  }
  if (Nphi <= 3) {
    Error("PoissonMultiGrid3D2D", "Poisson3DMultiGrid - Error in the number of Nphi. Must be larger than 3");
    return;
  }
  if (Nphi > 1000) {
    Error("PoissonMultiGrid3D2D", "Poisson3D  Nphi > 1000 is not allowed (nor wise) ");
    return;
  }

  const DataT gridSpacingR = getSpacingR();
  const DataT gridSpacingZ = getSpacingZ();
  const DataT gridSpacingPhi = getSpacingPhi();
  const DataT ratioPhi = gridSpacingR * gridSpacingR / (gridSpacingPhi * gridSpacingPhi); // ratio_{phi} = gridSize_{r} / gridSize_{phi}
  const DataT ratioZ = gridSpacingR * gridSpacingR / (gridSpacingZ * gridSpacingZ);       // ratio_{Z} = gridSize_{r} / gridSize_{z}

  // Solve Poisson's equation in cylindrical coordinates by multiGrid technique
  // Allow for different size grid spacing in R and Z directions
  int nGridRow = 0; // number grid
  int nGridCol = 0; // number grid
  int nnRow = Nr;
  int nnCol = Nz;

  while (nnRow >>= 1) {
    nGridRow++;
  }
  while (nnCol >>= 1) {
    nGridCol++;
  }

  const int maxVal = std::max(nGridRow, nGridCol); // Calculate the number of nLoop for the binary expansion
  const size_t nLoop = (maxVal > o2::tpc::MGParameters::maxLoop) ? o2::tpc::MGParameters::maxLoop : maxVal;
  unsigned int iOne = 1; // index i in gridSize r (original)
  unsigned int jOne = 1; // index j in gridSize z (original)

  std::vector<Matrix3D> tvArrayV(nLoop);     // potential <--> error
  std::vector<Matrix3D> tvChargeFMG(nLoop);  // charge is restricted in full multiGrid
  std::vector<Matrix3D> tvCharge(nLoop);     // charge <--> residue
  std::vector<Matrix3D> tvPrevArrayV(nLoop); // error calculation
  std::vector<Matrix3D> tvResidue(nLoop);    // residue calculation

  for (unsigned int count = 1; count <= nLoop; count++) {
    const unsigned int tnRRow = iOne == 1 ? Nr : Nr / iOne + 1;
    const unsigned int tnZColumn = jOne == 1 ? Nz : Nz / jOne + 1;
    tvResidue[count - 1].resize(tnRRow, tnZColumn, Nphi);
    tvPrevArrayV[count - 1].resize(tnRRow, tnZColumn, Nphi);

    // memory for the finest grid is from parameters
    tvChargeFMG[count - 1].resize(tnRRow, tnZColumn, Nphi);
    tvArrayV[count - 1].resize(tnRRow, tnZColumn, Nphi);
    tvCharge[count - 1].resize(tnRRow, tnZColumn, Nphi);

    if (count == 1) {
      std::move(matricesCharge.begin(), matricesCharge.end(), tvChargeFMG[count - 1].data());
      std::move(matricesV.begin(), matricesV.end(), tvArrayV[count - 1].data());
    } else {
      restrict3D(tvChargeFMG[count - 1], tvChargeFMG[count - 2], tnRRow, tnZColumn, Nphi, Nphi);
      restrictBoundary3D(tvArrayV[count - 1], tvArrayV[count - 2], tnRRow, tnZColumn, Nphi, Nphi);
    }
    iOne = 2 * iOne; // doubling
    jOne = 2 * jOne; // doubling
  }

  std::array<DataT, Nr> coefficient1{};        // coefficient1(Nr) for storing (1 + h_{r}/2r_{i}) from central differences in r direction
  std::array<DataT, Nr> coefficient2{};        // coefficient2(Nr) for storing (1 + h_{r}/2r_{i}) from central differences in r direction
  std::array<DataT, Nr> coefficient3{};        // coefficient3(Nr) for storing (1/r_{i}^2) from central differences in phi direction
  std::array<DataT, Nr> coefficient4{};        // coefficient4(Nr) for storing  1/2
  std::array<DataT, Nr> inverseCoefficient4{}; // inverse of coefficient4(Nr)

  // Case full multi grid (FMG)
  if (o2::tpc::MGParameters::cycleType == FCycle) {
    // 1) Relax on the coarsest grid
    iOne = iOne / 2;
    jOne = jOne / 2;
    int tnRRow = iOne == 1 ? Nr : Nr / iOne + 1;
    int tnZColumn = jOne == 1 ? Nz : Nz / jOne + 1;

    const DataT h = getSpacingR() * iOne;
    const DataT h2 = h * h;
    const DataT iOne2 = iOne * iOne;
    const DataT tempRatioPhi = ratioPhi * iOne2; // Used tobe divided by ( m_one * m_one ) when m_one was != 1
    const DataT tempRatioZ = ratioZ * iOne2 / (jOne * jOne);

    calcCoefficients(1, tnRRow - 1, h, tempRatioZ, tempRatioPhi, coefficient1, coefficient2, coefficient3, coefficient4);

    // relax on the coarsest level
    relax3D(tvArrayV[nLoop - 1], tvChargeFMG[nLoop - 1], tnRRow, tnZColumn, Nphi, symmetry, h2, tempRatioZ, coefficient1, coefficient2, coefficient3, coefficient4);

    // 2) Do multiGrid v-cycle from coarsest to finest
    for (int count = nLoop - 2; count >= 0; count--) {
      // move to finer grid
      iOne = iOne * 0.5;
      jOne = jOne * 0.5;
      tnRRow = iOne == 1 ? Nr : Nr / iOne + 1;
      tnZColumn = jOne == 1 ? Nz : Nz / jOne + 1;

      // 2) a) Interpolate potential for h -> 2h (coarse -> fine)
      interp3D(tvArrayV[count], tvArrayV[count + 1], tnRRow, tnZColumn, Nphi, Nphi);

      // 2) c) Copy the restricted charge to charge for calculation
      tvCharge[count] = tvChargeFMG[count]; //copy

      // 2) c) Do V cycle o2::tpc::MGParameters::nMGCycle times at most
      for (int mgCycle = 0; mgCycle < o2::tpc::MGParameters::nMGCycle; mgCycle++) {
        // Copy the potential to temp array for convergence calculation
        tvPrevArrayV[count] = tvArrayV[count];

        // 2) c) i) Call V cycle from grid count+1 (current fine level) to nLoop (coarsest)
        vCycle3D2D(symmetry, count + 1, nLoop, o2::tpc::MGParameters::nPre, o2::tpc::MGParameters::nPost, ratioZ, ratioPhi, tvArrayV, tvCharge, tvResidue, coefficient1, coefficient2, coefficient3, coefficient4, inverseCoefficient4);

        const DataT convergenceError = getConvergenceError(tvArrayV[count], tvPrevArrayV[count]);

        /// if already converge just break move to finer grid
        if (convergenceError <= sConvergenceError) {
          break;
        }
      }
    }
  } // Case V multi grid (VMG)

  std::move(tvArrayV[0].storage.begin(), tvArrayV[0].storage.end(), matricesV.data());
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nz, Nr, Nphi>::poissonMultiGrid3D(DataContainer& matricesV, const DataContainer& matricesCharge, const int symmetry)
{

  // const Float_t gridSpacingR = (AliTPCPoissonSolver<DataT>::fgkOFCRadius - TPCParameters<DataT>::IFCRADIUS) / (Nr - 1); // h_{r}
  // const Float_t gridSizeZ = AliTPCPoissonSolver<DataT>::fgkTPCZ0 / (Nz - 1);                // h_{z}
  // const Float_t ratioZ = gridSpacingR * gridSpacingR / (gridSizeZ * gridSizeZ);                  // ratio_{Z} = gridSize_{r} / gridSize_{z}

  const DataT gridSpacingR = getSpacingR();
  const DataT gridSpacingZ = getSpacingZ();
  const DataT gridSpacingPhi = getSpacingPhi();
  const DataT ratioPhi = gridSpacingR * gridSpacingR / (gridSpacingPhi * gridSpacingPhi); // ratio_{phi} = gridSize_{r} / gridSize_{phi}
  const DataT ratioZ = gridSpacingR * gridSpacingR / (gridSpacingZ * gridSpacingZ);       // ratio_{Z} = gridSize_{r} / gridSize_{z}

  Float_t gridSizePhi = TMath::TwoPi() / Nphi; // h_{phi}
  Float_t h, h2, radius;
  Float_t tempRatioPhi, tempRatioZ;

  Float_t convergenceError; // Convergence error

  Info("PoissonMultiGrid3D", "%s", Form("in Poisson Solver 3D multi grid full coarsening  Nr=%d, cols=%d, Nphi=%d \n", Nr, Nz, Nphi));

  // Check that the number of Nr and Nz is suitable for a binary expansion
  if (!isPowerOfTwo((Nr - 1))) {
    Error("PoissonMultiGrid3D", "Poisson3DMultiGrid - Error in the number of Nr. Must be 2**M + 1");
    return;
  }
  if (!isPowerOfTwo((Nz - 1))) {
    Error("PoissonMultiGrid3D", "Poisson3DMultiGrid - Error in the number of Nz. Must be 2**N - 1");
    return;
  }
  if (Nphi <= 3) {
    Error("PoissonMultiGrid3D", "Poisson3DMultiGrid - Error in the number of Nphi. Must be larger than 3");
    return;
  }
  if (Nphi > 1000) {
    Error("PoissonMultiGrid3D", "Poisson3D  Nphi > 1000 is not allowed (nor wise) ");
    return;
  }

  // Solve Poisson's equation in cylindrical coordinates by multi grid technique
  // Allow for different size grid spacing in R and Z directions

  Int_t nGridRow = 0; // number grid
  Int_t nGridCol = 0; // number grid
  Int_t nGridPhi = 0;

  Int_t nnRow;
  Int_t nnCol;
  Int_t nnPhi;

  nnRow = Nr;
  while (nnRow >>= 1) {
    nGridRow++;
  }

  nnCol = Nz;
  while (nnCol >>= 1) {
    nGridCol++;
  }

  nnPhi = Nphi;

  while (nnPhi % 2 == 0) {
    nGridPhi++;
    nnPhi /= 2;
  }

  Info("PoissonMultiGrid3D", "%s", Form("nGridRow=%d, nGridCol=%d, nGridPhi=%d", nGridRow, nGridCol, nGridPhi));
  Int_t nLoop = TMath::Max(nGridRow, nGridCol); // Calculate the number of nLoop for the binary expansion
  nLoop = TMath::Max(nLoop, nGridPhi);

  // Vector for storing multi grid array
  Int_t iOne = 1; // index i in gridSize r (original)
  Int_t jOne = 1; // index j in gridSize z (original)
  Int_t kOne = 1; // index k in gridSize phi
  Int_t tnRRow = Nr, tnZColumn = Nz, tPhiSlice = Nphi, otPhiSlice;

  // 1)	Memory allocation for multi grid
  std::vector<Matrix3D> tvArrayV(nLoop);     // potential <--> error
  std::vector<Matrix3D> tvChargeFMG(nLoop);  // charge is restricted in full multiGrid
  std::vector<Matrix3D> tvCharge(nLoop);     // charge <--> residue
  std::vector<Matrix3D> tvPrevArrayV(nLoop); // error calculation
  std::vector<Matrix3D> tvResidue(nLoop);    // residue calculation

  std::array<DataT, Nr> coefficient1{};        // coefficient1(Nr) for storing (1 + h_{r}/2r_{i}) from central differences in r direction
  std::array<DataT, Nr> coefficient2{};        // coefficient2(Nr) for storing (1 + h_{r}/2r_{i}) from central differences in r direction
  std::array<DataT, Nr> coefficient3{};        // coefficient3(Nr) for storing (1/r_{i}^2) from central differences in phi direction
  std::array<DataT, Nr> coefficient4{};        // coefficient4(Nr) for storing  1/2
  std::array<DataT, Nr> inverseCoefficient4{}; // inverse of coefficient4(Nr)

  for (Int_t count = 1; count <= nLoop; count++) {
    // tnRRow,tnZColumn in new grid
    tnRRow = iOne == 1 ? Nr : Nr / iOne + 1;
    tnZColumn = jOne == 1 ? Nz : Nz / jOne + 1;
    tPhiSlice = kOne == 1 ? Nphi : Nphi / kOne;
    tPhiSlice = tPhiSlice < nnPhi ? nnPhi : tPhiSlice;

    // allocate memory for residue
    tvResidue[count - 1].resize(tnRRow, tnZColumn, tPhiSlice);
    tvPrevArrayV[count - 1].resize(tnRRow, tnZColumn, tPhiSlice);

    tvChargeFMG[count - 1].resize(tnRRow, tnZColumn, tPhiSlice);
    tvArrayV[count - 1].resize(tnRRow, tnZColumn, tPhiSlice);
    tvCharge[count - 1].resize(tnRRow, tnZColumn, tPhiSlice);

    // memory for the finest grid is from parameters
    if (count == 1) {
      // tvChargeFMG[count - 1] = matricesCharge;
      // tvArrayV[count - 1] = matricesV;
      // tvCharge[count - 1] = matricesCharge;
      std::move(matricesCharge.begin(), matricesCharge.end(), tvChargeFMG[count - 1].data());
      std::move(matricesV.begin(), matricesV.end(), tvArrayV[count - 1].data());
      tvCharge[count - 1] = tvChargeFMG[count - 1];
    }
    // } else {
      // allocate for coarser grid
      // tvChargeFMG[count - 1] = new TMatrixD*[tPhiSlice];
      // tvArrayV[count - 1] = new TMatrixD*[tPhiSlice];
      // tvCharge[count - 1] = new TMatrixD*[tPhiSlice];
      // for (Int_t k = 0; k < tPhiSlice; k++) {
      //   tvArrayV[count - 1][k] = new TMatrixD(tnRRow, tnZColumn);
      //   tvCharge[count - 1][k] = new TMatrixD(tnRRow, tnZColumn);
      //   tvChargeFMG[count - 1][k] = new TMatrixD(tnRRow, tnZColumn);
      // }
    // }
    iOne = 2 * iOne; // doubling
    jOne = 2 * jOne; // doubling
    kOne = 2 * kOne;
  }

  // Case full multi grid (FMG)

  if (o2::tpc::MGParameters::cycleType == FCycle) {
    // Restrict the charge to coarser grid
    iOne = 2;
    jOne = 2;
    kOne = 2;
    otPhiSlice = Nphi;

    // 1) Restrict Charge and Boundary to coarser grid
    for (Int_t count = 2; count <= nLoop; count++) {
      // tnRRow,tnZColumn in new grid
      tnRRow = iOne == 1 ? Nr : Nr / iOne + 1;
      tnZColumn = jOne == 1 ? Nz : Nz / jOne + 1;
      tPhiSlice = kOne == 1 ? Nphi : Nphi / kOne;
      tPhiSlice = tPhiSlice < nnPhi ? nnPhi : tPhiSlice;

      Info("PoissonMultiGrid3D", "%s", Form("Restrict3D, tnRRow=%d, tnZColumn=%d, newPhiSlice=%d, oldPhiSlice=%d\n", tnRRow, tnZColumn, tPhiSlice, otPhiSlice));
      restrict3D(tvChargeFMG[count - 1], tvChargeFMG[count - 2], tnRRow, tnZColumn, tPhiSlice, otPhiSlice);
      // copy boundary values of V
      restrictBoundary3D(tvArrayV[count - 1], tvArrayV[count - 2], tnRRow, tnZColumn, tPhiSlice, otPhiSlice);
      otPhiSlice = tPhiSlice;

      iOne = 2 * iOne; // doubling
      jOne = 2 * jOne; // doubling
      kOne = 2 * kOne;
    }

    // Relax on the coarsest grid
    // FMG
    // 2) Relax on the coarsest grid

    // move to the coarsest + 1
    iOne = iOne / 2;
    jOne = jOne / 2;
    kOne = kOne / 2;
    tnRRow = iOne == 1 ? Nr : Nr / iOne + 1;
    tnZColumn = jOne == 1 ? Nz : Nz / jOne + 1;
    tPhiSlice = kOne == 1 ? Nphi : Nphi / kOne;
    tPhiSlice = tPhiSlice < nnPhi ? nnPhi : tPhiSlice;
    otPhiSlice = tPhiSlice;

    h = gridSpacingR * iOne;
    h2 = h * h;
    gridSizePhi = TMath::TwoPi() / tPhiSlice;           // h_{phi}
    tempRatioPhi = h * h / (gridSizePhi * gridSizePhi); // ratio_{phi} = gridSize_{r} / gridSize_{phi}
    tempRatioZ = ratioZ * iOne * iOne / (jOne * jOne);
    for (Int_t i = 1; i < tnRRow - 1; i++) {
      radius = TPCParameters<DataT>::IFCRADIUS + i * h;
      coefficient1[i] = 1.0 + h / (2 * radius);
      coefficient2[i] = 1.0 - h / (2 * radius);
      coefficient3[i] = tempRatioPhi / (radius * radius);
      coefficient4[i] = 0.5 / (1.0 + tempRatioZ + coefficient3[i]);
    }
    // 3) Relax on the coarsest grid
    // std::cout<<"A"<<std::endl;
    relax3D(tvArrayV[nLoop - 1], tvChargeFMG[nLoop - 1], tnRRow, tnZColumn, tPhiSlice, symmetry, h2, tempRatioZ, coefficient1, coefficient2, coefficient3, coefficient4);
    // std::cout<<"B"<<std::endl;

    // 4) V Cycle from coarsest to finest
    for (Int_t count = nLoop - 2; count >= 0; count--) {
      // move to finer grid
      // coefficient1.clear();
      // coefficient2.clear();
      // coefficient3.clear();
      // coefficient4.clear();
      // inverseCoefficient4.clear();
      std::fill(std::begin(coefficient1), std::end(coefficient1), 0);
      std::fill(std::begin(coefficient2), std::end(coefficient2), 0);
      std::fill(std::begin(coefficient3), std::end(coefficient3), 0);
      std::fill(std::begin(coefficient4), std::end(coefficient4), 0);
      std::fill(std::begin(inverseCoefficient4), std::end(inverseCoefficient4), 0);

      iOne = iOne / 2;
      jOne = jOne / 2;
      kOne = kOne / 2;

      tnRRow = iOne == 1 ? Nr : Nr / iOne + 1;
      tnZColumn = jOne == 1 ? Nz : Nz / jOne + 1;
      tPhiSlice = kOne == 1 ? Nphi : Nphi / kOne;
      tPhiSlice = tPhiSlice < nnPhi ? nnPhi : tPhiSlice;
      // 4) a) interpolate from 2h --> h grid
      // std::cout<<"C"<<std::endl;
      interp3D(tvArrayV[count], tvArrayV[count + 1], tnRRow, tnZColumn, tPhiSlice, otPhiSlice);
      // std::cout<<"D"<<std::endl;

      // Copy the relax charge to the tvCharge
      if (count > 0) {
        // for (Int_t m = 0; m < tPhiSlice; m++) {
          // *tvCharge[count][m] = *tvChargeFMG[count][m]; //copy
        // }
        tvCharge[count] = tvChargeFMG[count];
      }
      for (Int_t mgCycle = 0; mgCycle < o2::tpc::MGParameters::nMGCycle; mgCycle++) {
        // copy to store previous potential
        // for (Int_t m = 0; m < tPhiSlice; m++) {
          // *tvPrevArrayV[count][m] = *tvArrayV[count][m]; //copy
        // }
        tvPrevArrayV[count] = tvArrayV[count];
        // std::cout<<"E"<<std::endl;
        vCycle3D(symmetry, count + 1, nLoop, o2::tpc::MGParameters::nPre, o2::tpc::MGParameters::nPost, ratioZ, tvArrayV, tvCharge, tvResidue, coefficient1, coefficient2, coefficient3, coefficient4, inverseCoefficient4);

        /// converge error
        // std::cout<<"F"<<std::endl;
        convergenceError = getConvergenceError(tvArrayV[count], tvPrevArrayV[count]);
        // std::cout<<"G"<<std::endl;
        //// error counting /////
        if (count == 0) {
          // (*fErrorConvergenceNormInf)(mgCycle) = convergenceError;
          // (*fError)(mgCycle) = GetExactError(matricesV, tvPrevArrayV[count], Nphi);
        }
        /// if already converge just break move to finer grid
        if (convergenceError <= sConvergenceError) {
          // fIterations = mgCycle + 1;
          break;
        }
      }
      // keep old slice information
      otPhiSlice = tPhiSlice;
    }

  } else if (o2::tpc::MGParameters::cycleType == VCycle) {
    // V-cycle
    Int_t gridFrom = 1;
    Int_t gridTo = nLoop;

    for (Int_t mgCycle = 0; mgCycle < o2::tpc::MGParameters::nMGCycle; mgCycle++) {
      // copy to store previous potential
      // for (Int_t m = 0; m < Nphi; m++) {
        // *tvPrevArrayV[0][m] = *tvArrayV[0][m]; //copy
      // }
      tvPrevArrayV[0] = tvArrayV[0];
      // Do V Cycle from the coarsest to finest grid
      vCycle3D(symmetry, gridFrom, gridTo, o2::tpc::MGParameters::nPre, o2::tpc::MGParameters::nPost, ratioZ, tvArrayV, tvCharge, tvResidue, coefficient1, coefficient2, coefficient3, coefficient4, inverseCoefficient4);
      // convergence error
      convergenceError = getConvergenceError(tvArrayV[0], tvPrevArrayV[0]);
      // (*fErrorConvergenceNormInf)(mgCycle) = convergenceError;
      // (*fError)(mgCycle) = GetExactError(matricesV, tvPrevArrayV[0], Nphi);
      // if error already achieved then stop mg iteration
      if (convergenceError <= sConvergenceError) {
        //Info("PoissonMultiGrid3D",Form("Exact Err: %f, MG Iteration : %d", (*fError)(mgCycle), mgCycle));
        // fIterations = mgCycle + 1;
        break;
      }
    }
  }
  std::move(tvArrayV[0].storage.begin(), tvArrayV[0].storage.end(), matricesV.data());
}


template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nz, Nr, Nphi>::vCycle3D2D(const int symmetry, const int gridFrom, const int gridTo, const int nPre, const int nPost, const DataT ratioZ, const DataT ratioPhi,
                                                         std::vector<Matrix3D>& tvArrayV, std::vector<Matrix3D>& tvCharge, std::vector<Matrix3D>& tvResidue, std::array<DataT, Nr>& coefficient1,
                                                         std::array<DataT, Nr>& coefficient2, std::array<DataT, Nr>& coefficient3, std::array<DataT, Nr>& coefficient4, std::array<DataT, Nr>& inverseCoefficient4) const
{
  int iOne = 1 << (gridFrom - 1);
  int jOne = 1 << (gridFrom - 1);
  int tnRRow = iOne == 1 ? Nr : Nr / iOne + 1;
  int tnZColumn = jOne == 1 ? Nz : Nz / jOne + 1;

  for (int count = gridFrom; count <= gridTo - 1; ++count) {
    const DataT h = getSpacingR() * iOne;
    const DataT h2 = h * h;
    const DataT ih2 = 1.0 / h2;
    const int iOne2 = iOne * iOne;
    const DataT tempRatioPhi = ratioPhi * iOne2; // Used tobe divided by ( m_one * m_one ) when m_one was != 1
    const DataT tempRatioZ = ratioZ * iOne2 / (jOne * jOne);

    calcCoefficients(1, tnRRow - 1, h, tempRatioZ, tempRatioPhi, coefficient1, coefficient2, coefficient3, coefficient4);
    for (unsigned int i = 1; i < tnRRow - 1; ++i) {
      inverseCoefficient4[i] = 1.0 / coefficient4[i];
    }

    //Info("VCycle3D2D","Before Pre-smoothing");
    // 1) Pre-Smoothing: Gauss-Seidel Relaxation or Jacobi
    for (int jPre = 1; jPre <= nPre; ++jPre) {
      relax3D(tvArrayV[count - 1], tvCharge[count - 1], tnRRow, tnZColumn, Nphi, symmetry, h2, tempRatioZ, coefficient1, coefficient2, coefficient3, coefficient4);
    } // end pre smoothing

    // 2) Residue calculation
    residue3D(tvResidue[count - 1], tvArrayV[count - 1], tvCharge[count - 1], tnRRow, tnZColumn, Nphi, symmetry, ih2, tempRatioZ, coefficient1, coefficient2, coefficient3, inverseCoefficient4);

    iOne = 2 * iOne;
    jOne = 2 * jOne;
    tnRRow = iOne == 1 ? Nr : Nr / iOne + 1;
    tnZColumn = jOne == 1 ? Nz : Nz / jOne + 1;

    //3) Restriction
    restrict3D(tvCharge[count], tvResidue[count - 1], tnRRow, tnZColumn, Nphi, Nphi);

    //4) Zeroing coarser V
    std::fill(tvArrayV[count].storage.begin(), tvArrayV[count].storage.end(), 0);
  }

  // coarsest grid
  const DataT h = getSpacingR() * iOne;
  const DataT h2 = h * h;

  const int iOne2 = iOne * iOne;
  const DataT tempRatioPhi = ratioPhi * iOne2;
  const DataT tempRatioZ = ratioZ * iOne2 / (jOne * jOne);

  calcCoefficients(1, tnRRow - 1, h, tempRatioZ, tempRatioPhi, coefficient1, coefficient2, coefficient3, coefficient4);

  // 3) Relax on the coarsest grid
  relax3D(tvArrayV[gridTo - 1], tvCharge[gridTo - 1], tnRRow, tnZColumn, Nphi, symmetry, h2, tempRatioZ, coefficient1, coefficient2, coefficient3, coefficient4);

  // back to fine
  for (int count = gridTo - 1; count >= gridFrom; --count) {
    iOne = iOne / 2;
    jOne = jOne / 2;

    tnRRow = iOne == 1 ? Nr : Nr / iOne + 1;
    tnZColumn = jOne == 1 ? Nz : Nz / jOne + 1;

    const DataT h = getSpacingR() * iOne;
    const DataT h2 = h * h;
    const int iOne2 = iOne * iOne;
    const DataT tempRatioPhi = ratioPhi * iOne2; // Used tobe divided by ( m_one * m_one ) when m_one was != 1
    const DataT tempRatioZ = ratioZ * iOne2 / (jOne * jOne);

    // 4) Interpolation/Prolongation
    addInterp3D(tvArrayV[count - 1], tvArrayV[count], tnRRow, tnZColumn, Nphi, Nphi);

    calcCoefficients(1, tnRRow - 1, h, tempRatioZ, tempRatioPhi, coefficient1, coefficient2, coefficient3, coefficient4);

    // 5) Post-Smoothing: Gauss-Seidel Relaxation
    for (int jPost = 1; jPost <= nPost; ++jPost) {
      relax3D(tvArrayV[count - 1], tvCharge[count - 1], tnRRow, tnZColumn, Nphi, symmetry, h2, tempRatioZ, coefficient1, coefficient2, coefficient3, coefficient4);
    } // end post smoothing
  }
}


template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nz, Nr, Nphi>::vCycle3D(const int symmetry, const int gridFrom, const int gridTo, const int nPre, const int nPost, const DataT ratioZ, std::vector<Matrix3D>& tvArrayV,
                std::vector<Matrix3D>& tvCharge, std::vector<Matrix3D>& tvResidue, std::array<DataT, Nr>& coefficient1, std::array<DataT, Nr>& coefficient2, std::array<DataT, Nr>& coefficient3,
                std::array<DataT, Nr>& coefficient4, std::array<DataT, Nr>& inverseCoefficient4) const
{

  Float_t h, h2, ih2, tempRatioZ, tempRatioPhi, radius, tempGridSizePhi;
  Int_t iOne, jOne, kOne, tnRRow, tnZColumn, tPhiSlice, otPhiSlice, count, nnPhi;
  const DataT gridSpacingR = getSpacingR();

  iOne = 1 << (gridFrom - 1);
  jOne = 1 << (gridFrom - 1);
  kOne = 1 << (gridFrom - 1);

  nnPhi = Nphi;

  while (nnPhi % 2 == 0) {
    nnPhi /= 2;
  }

  tnRRow = iOne == 1 ? Nr : Nr / iOne + 1;
  tnZColumn = jOne == 1 ? Nz : Nz / jOne + 1;
  tPhiSlice = kOne == 1 ? Nphi : Nphi / kOne;
  tPhiSlice = tPhiSlice < nnPhi ? nnPhi : tPhiSlice;

  //Info("VCycle3D",Form("Grid information: tnRRow=%d, tcols=%d, tPhiSlice=%d\n", tnRRow,tnZColumn,tPhiSlice));

  for (count = gridFrom; count <= gridTo - 1; ++count) {
    otPhiSlice = tPhiSlice;

    h = gridSpacingR * iOne;
    h2 = h * h;
    ih2 = 1.0 / h2;
    tempGridSizePhi = TMath::TwoPi() / tPhiSlice; // phi now is multiGrid

    tempRatioPhi = h * h / (tempGridSizePhi * tempGridSizePhi); // ratio_{phi} = gridSize_{r} / gridSize_{phi}

    tempRatioZ = ratioZ * iOne * iOne / (jOne * jOne);

    for (Int_t i = 1; i < tnRRow - 1; i++) {
      radius = TPCParameters<DataT>::IFCRADIUS + i * h;
      // const DataT radiusInv = 1. / (TPCParameters<DataT>::IFCRADIUS + i * h);
      coefficient1[i] = 1.0 + h / (2 * radius);
      coefficient2[i] = 1.0 - h / (2 * radius);
      coefficient3[i] = tempRatioPhi / (radius * radius);
      coefficient4[i] = 0.5 / (1.0 + tempRatioZ + coefficient3[i]);
      inverseCoefficient4[i] = 1.0 / coefficient4[i];
    }

    // matricesCurrentV = tvArrayV[count - 1];
    // matricesCurrentCharge = tvCharge[count - 1];
    // residue = tvResidue[count - 1];

    //Info("VCycle3D","Before Pre-smoothing");
    //matricesCurrentV->Print();

    // 1) Pre-Smoothing: Gauss-Seidel Relaxation or Jacobi
    for (Int_t jPre = 1; jPre <= nPre; jPre++) {
      relax3D(tvArrayV[count - 1], tvCharge[count - 1], tnRRow, tnZColumn, tPhiSlice, symmetry, h2, tempRatioZ, coefficient1, coefficient2, coefficient3, coefficient4);
    } // end pre smoothing

    // 2) Residue calculation
// std::cout<<"1:"<<std::endl;
    residue3D(tvResidue[count - 1], tvArrayV[count - 1], tvCharge[count - 1], tnRRow, tnZColumn, tPhiSlice, symmetry, ih2, tempRatioZ, coefficient1, coefficient2, coefficient3, inverseCoefficient4);

    iOne = 2 * iOne;
    jOne = 2 * jOne;
    kOne = 2 * kOne;
    tnRRow = iOne == 1 ? Nr : Nr / iOne + 1;
    tnZColumn = jOne == 1 ? Nz : Nz / jOne + 1;
    tPhiSlice = Nphi / kOne;
    tPhiSlice = tPhiSlice < nnPhi ? nnPhi : tPhiSlice;

    // matricesCurrentCharge = tvCharge[count];
    // matricesCurrentV = tvArrayV[count];
    //3) Restriction
    restrict3D(tvCharge[count], tvResidue[count - 1], tnRRow, tnZColumn, tPhiSlice, otPhiSlice);
    // std::cout<<"2:"<<std::endl;

    //4) Zeroing coarser V
    // for (Int_t m = 0; m < tPhiSlice; m++) {
      // matricesCurrentV[m]->Zero();
    // }
    // std::cout<<"tvArAA: "<< tvArrayV[count].nR << std::endl;
    // std::cout<<"count: "<< count << std::endl;
    std::fill(tvArrayV[count].storage.begin(), tvArrayV[count].storage.end(), 0);
    // std::cout<<"coun2t: "<< count << std::endl;
  }

  // std::cout<<"count END: "<< count << std::endl;
  // std::cout<<"gridTo - 1: "<< gridTo - 1 << std::endl;
  count = gridTo - 1;

  // coarsest grid
  h = gridSpacingR * iOne;
  h2 = h * h;
  tempGridSizePhi = TMath::TwoPi() / tPhiSlice; // phi now is multiGrid

  tempRatioPhi = h * h / (tempGridSizePhi * tempGridSizePhi); // ratio_{phi} = gridSize_{r} / gridSize_{phi}
  tempRatioZ = ratioZ * iOne * iOne / (jOne * jOne);

  for (Int_t i = 1; i < tnRRow - 1; i++) {
    radius = TPCParameters<DataT>::IFCRADIUS + i * h;
    coefficient1[i] = 1.0 + h / (2 * radius);
    coefficient2[i] = 1.0 - h / (2 * radius);
    coefficient3[i] = tempRatioPhi / (radius * radius);
    coefficient4[i] = 0.5 / (1.0 + tempRatioZ + coefficient3[i]);
  }
  // std::cout<<"3:"<<std::endl;
  // std::cout<<"tvArrayV[count]: "<< tvArrayV[count].nR << std::endl;
  // std::cout<<"count: "<< count << std::endl;
  // 3) Relax on the coarsest grid
  relax3D(tvArrayV[count], tvCharge[count], tnRRow, tnZColumn, tPhiSlice, symmetry, h2, tempRatioZ, coefficient1, coefficient2, coefficient3, coefficient4);
  // std::cout<<"3a:"<<std::endl;
  // back to fine
  for (count = gridTo - 1; count >= gridFrom; count--) {
    otPhiSlice = tPhiSlice;

    iOne = iOne / 2;
    jOne = jOne / 2;
    kOne = kOne / 2;

    tnRRow = iOne == 1 ? Nr : Nr / iOne + 1;
    tnZColumn = jOne == 1 ? Nz : Nz / jOne + 1;
    tPhiSlice = kOne == 1 ? Nphi : Nphi / kOne;
    tPhiSlice = tPhiSlice < nnPhi ? nnPhi : tPhiSlice;

    h = gridSpacingR * iOne;
    h2 = h * h;
    tempGridSizePhi = TMath::TwoPi() / tPhiSlice; // phi now is multiGrid

    tempRatioPhi = h * h / (tempGridSizePhi * tempGridSizePhi); // ratio_{phi} = gridSize_{r} / gridSize_{phi}

    tempRatioZ = ratioZ * iOne * iOne / (jOne * jOne);
  // std::cout<<"4:"<<std::endl;
    // matricesCurrentCharge = tvCharge[count - 1];
    // matricesCurrentV = tvArrayV[count - 1];
    // matricesCurrentVC = tvArrayV[count];

    // 4) Interpolation/Prolongation

    addInterp3D(tvArrayV[count - 1], tvArrayV[count], tnRRow, tnZColumn, tPhiSlice, otPhiSlice);

    for (Int_t i = 1; i < tnRRow - 1; i++) {
      radius = TPCParameters<DataT>::IFCRADIUS + i * h;
      coefficient1[i] = 1.0 + h / (2 * radius);
      coefficient2[i] = 1.0 - h / (2 * radius);
      coefficient3[i] = tempRatioPhi / (radius * radius);
      coefficient4[i] = 0.5 / (1.0 + tempRatioZ + coefficient3[i]);
    }
  // std::cout<<"5:"<<std::endl;
    // 5) Post-Smoothing: Gauss-Seidel Relaxation
    for (Int_t jPost = 1; jPost <= nPost; jPost++) {
      relax3D(tvArrayV[count - 1], tvCharge[count - 1], tnRRow, tnZColumn, tPhiSlice, symmetry, h2, tempRatioZ, coefficient1, coefficient2, coefficient3, coefficient4);
    }
  }
}


template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nz, Nr, Nphi>::residue3D(Matrix3D& residue, const Matrix3D& matricesCurrentV, const Matrix3D& matricesCurrentCharge, const int tnRRow, const int tnZColumn, const int iPhi, const int symmetry,
                                                        const DataT ih2, const DataT tempRatioZ, const std::array<DataT, Nr>& coefficient1, const std::array<DataT, Nr>& coefficient2, const std::array<DataT, Nr>& coefficient3, const std::array<DataT, Nr>& inverseCoefficient4) const
{
#pragma omp parallel for // parallising this loop is possible - but using more than 2 cores makes it slower -
  for (int m = 0; m < iPhi; ++m) {
    int mp1 = m + 1;
    int signPlus = 1;
    int mm1 = m - 1;
    int signMinus = 1;

    // Reflection symmetry in phi (e.g. symmetry at sector boundaries, or half sectors, etc.)
    if (symmetry == 1) {
      if (mp1 > iPhi - 1) {
        mp1 = iPhi - 2;
      }
      if (mm1 < 0) {
        mm1 = 1;
      }
    }
    // Anti-symmetry in phi
    else if (symmetry == -1) {
      if (mp1 > iPhi - 1) {
        mp1 = iPhi - 2;
        signPlus = -1;
      }
      if (mm1 < 0) {
        mm1 = 1;
        signMinus = -1;
      }
    } else { // No Symmetries in phi, no boundaries, the calculation is continuous across all phi
      if (mp1 > iPhi - 1) {
        mp1 = m + 1 - iPhi;
      }
      if (mm1 < 0) {
        mm1 = m - 1 + iPhi;
      }
    }

    for (int j = 1; j < tnZColumn - 1; ++j) {
      for (int i = 1; i < tnRRow - 1; ++i) {
        residue(i, j, m) = ih2 * (coefficient2[i] * matricesCurrentV(i - 1, j, m) + tempRatioZ * (matricesCurrentV(i, j - 1, m) + matricesCurrentV(i, j + 1, m)) + coefficient1[i] * matricesCurrentV(i + 1, j, m) +
                                  coefficient3[i] * (signPlus * matricesCurrentV(i, j, mp1) + signMinus * matricesCurrentV(i, j, mm1)) - inverseCoefficient4[i] * matricesCurrentV(i, j, m)) +
                           matricesCurrentCharge(i, j, m);
      } // end cols
    }   // end Nr
  }
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nz, Nr, Nphi>::interp3D(Matrix3D& matricesCurrentV, const Matrix3D& matricesCurrentVC, const int tnRRow, const int tnZColumn, const int newPhiSlice, const int oldPhiSlice) const
{
  // Do restrict 2 D for each slice
  if (newPhiSlice == 2 * oldPhiSlice) {
    for (int m = 0; m < newPhiSlice; m += 2) {
      // assuming no symmetry
      int mm = m * 0.5;
      int mmPlus = mm + 1;
      int mp1 = m + 1;

      // round
      if (mmPlus > oldPhiSlice - 1) {
        mmPlus = mm + 1 - oldPhiSlice;
      }
      if (mp1 > newPhiSlice - 1) {
        mp1 = m + 1 - newPhiSlice;
      }

      for (int j = 2; j < tnZColumn - 1; j += 2) {
        for (int i = 2; i < tnRRow - 1; i += 2) {
          const int iHalf = i * 0.5;
          const int jHalf = j * 0.5;
          matricesCurrentV(i, j, m) = matricesCurrentVC(iHalf, jHalf, mm);
          // point on corner lines at phi direction
          matricesCurrentV(i, j, mp1) = 0.5 * (matricesCurrentVC(iHalf, jHalf, mm) + matricesCurrentVC(iHalf, jHalf, mmPlus));
        }
      }

      for (int j = 1; j < tnZColumn - 1; j += 2) {
        for (int i = 2; i < tnRRow - 1; i += 2) {
          const int iHalf = i * 0.5;
          const int jHalf = j * 0.5;
          matricesCurrentV(i, j, m) = 0.5 * (matricesCurrentVC(iHalf, jHalf, mm) + matricesCurrentVC(iHalf, jHalf + 1, mm));
          // point on corner lines at phi direction
          matricesCurrentV(i, j, mp1) = 0.25 * (matricesCurrentVC(iHalf, jHalf, mm) + matricesCurrentVC(iHalf, jHalf + 1, mm) + matricesCurrentVC(iHalf, jHalf, mmPlus) + matricesCurrentVC(iHalf, jHalf + 1, mmPlus));
        }
      }

      for (int j = 2; j < tnZColumn - 1; j += 2) {
        for (int i = 1; i < tnRRow - 1; i += 2) {
          const int iHalf = i * 0.5;
          const int jHalf = j * 0.5;
          matricesCurrentV(i, j, m) = 0.5 * (matricesCurrentVC(iHalf, jHalf, mm) + matricesCurrentVC(iHalf + 1, jHalf, mm));
          // point on line at phi direction
          matricesCurrentV(i, j, mp1) = 0.25 * ((matricesCurrentVC(iHalf, jHalf, mm) + matricesCurrentVC(iHalf, jHalf, mmPlus)) + (matricesCurrentVC(iHalf + 1, jHalf, mmPlus) + matricesCurrentVC(iHalf + 1, jHalf, mm)));
        }
      }

      for (int j = 1; j < tnZColumn - 1; j += 2) {
        for (int i = 1; i < tnRRow - 1; i += 2) {
          const int iHalf = i * 0.5;
          const int jHalf = j * 0.5;
          matricesCurrentV(i, j, m) = 0.25 * ((matricesCurrentVC(iHalf, jHalf, mm) + matricesCurrentVC(iHalf, jHalf + 1, mm)) + (matricesCurrentVC(iHalf + 1, jHalf, mm) + matricesCurrentVC(iHalf + 1, jHalf + 1, mm)));
          // point at the center at phi direction
          matricesCurrentV(i, j, mp1) = 0.125 * ((matricesCurrentVC(iHalf, jHalf, mm) + matricesCurrentVC(iHalf, jHalf + 1, mm) + matricesCurrentVC(iHalf, jHalf, mmPlus) + matricesCurrentVC(iHalf, jHalf + 1, mmPlus)) +
                                                 (matricesCurrentVC(iHalf + 1, jHalf, mm) + matricesCurrentVC(iHalf + 1, jHalf + 1, mm) + matricesCurrentVC(iHalf + 1, jHalf, mmPlus) + matricesCurrentVC(iHalf + 1, jHalf + 1, mmPlus)));
        }
      }
    }

  } else {
#pragma omp parallel for // no change
    for (int m = 0; m < newPhiSlice; m++) {
      interp2D(matricesCurrentV, matricesCurrentVC, tnRRow, tnZColumn, m);
    }
  }
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nz, Nr, Nphi>::interp2D(Matrix3D& matricesCurrentV, const Matrix3D& matricesCurrentVC, const int tnRRow, const int tnZColumn, const int iphi) const
{
  for (int j = 2; j < tnZColumn - 1; j += 2) {
    for (int i = 2; i < tnRRow - 1; i += 2) {
      matricesCurrentV(i, j, iphi) = matricesCurrentVC(i / 2, j / 2, iphi);
    }
  }

  for (int j = 1; j < tnZColumn - 1; j += 2) {
    for (int i = 2; i < tnRRow - 1; i += 2) {
      const int iHalf = i * 0.5;
      const int jHalf = j * 0.5;
      matricesCurrentV(i, j, iphi) = 0.5 * (matricesCurrentVC(iHalf, jHalf, iphi) + matricesCurrentVC(iHalf, jHalf + 1, iphi));
    }
  }

  for (int j = 2; j < tnZColumn - 1; j += 2) {
    for (int i = 1; i < tnRRow - 1; i += 2) {
      const int iHalf = i * 0.5;
      const int jHalf = j * 0.5;
      matricesCurrentV(i, j, iphi) = 0.5 * (matricesCurrentVC(iHalf, jHalf, iphi) + matricesCurrentVC(iHalf + 1, jHalf, iphi));
    }
  }

  // only if full
  if (o2::tpc::MGParameters::gtType == Full) {
    for (int j = 1; j < tnZColumn - 1; j += 2) {
      for (int i = 1; i < tnRRow - 1; i += 2) {
        const int iHalf = i * 0.5;
        const int jHalf = j * 0.5;
        matricesCurrentV(i, j, iphi) = 0.25 * (matricesCurrentVC(iHalf, jHalf, iphi) + matricesCurrentVC(iHalf, jHalf + 1, iphi) + matricesCurrentVC(iHalf + 1, jHalf, iphi) + matricesCurrentVC(iHalf + 1, jHalf + 1, iphi));
      }
    }
  }
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nz, Nr, Nphi>::addInterp3D(Matrix3D& matricesCurrentV, const Matrix3D& matricesCurrentVC, const int tnRRow, const int tnZColumn, const int newPhiSlice, const int oldPhiSlice) const
{
  // Do restrict 2 D for each slice
  if (newPhiSlice == 2 * oldPhiSlice) {
    for (int m = 0; m < newPhiSlice; m += 2) {
      // assuming no symmetry
      int mm = m * 0.5;
      int mmPlus = mm + 1;
      int mp1 = m + 1;

      // round
      if (mmPlus > (oldPhiSlice)-1) {
        mmPlus = mm + 1 - (oldPhiSlice);
      }
      if (mp1 > (newPhiSlice)-1) {
        mp1 = m + 1 - (newPhiSlice);
      }

      for (int j = 2; j < tnZColumn - 1; j += 2) {
        for (int i = 2; i < tnRRow - 1; i += 2) {
          const int iHalf = i * 0.5;
          const int jHalf = j * 0.5;
          matricesCurrentV(i, j, m) += matricesCurrentVC(iHalf, jHalf, mm);
          // point on corner lines at phi direction
          matricesCurrentV(i, j, mp1) += 0.5 * (matricesCurrentVC(iHalf, jHalf, mm) + matricesCurrentVC(iHalf, jHalf, mmPlus));
        }
      }

      for (int j = 1; j < tnZColumn - 1; j += 2) {
        for (int i = 2; i < tnRRow - 1; i += 2) {
          const int iHalf = i * 0.5;
          const int jHalf = j * 0.5;
          matricesCurrentV(i, j, m) += 0.5 * (matricesCurrentVC(iHalf, jHalf, mm) + matricesCurrentVC(iHalf, jHalf + 1, mm));
          // point on corner lines at phi direction
          matricesCurrentV(i, j, mp1) += 0.25 * (matricesCurrentVC(iHalf, jHalf, mm) + matricesCurrentVC(iHalf, jHalf + 1, mm) + matricesCurrentVC(iHalf, jHalf, mmPlus) + matricesCurrentVC(iHalf, jHalf + 1, mmPlus));
        }
      }

      for (int j = 2; j < tnZColumn - 1; j += 2) {
        for (int i = 1; i < tnRRow - 1; i += 2) {
          const int iHalf = i * 0.5;
          const int jHalf = j * 0.5;
          matricesCurrentV(i, j, m) += 0.5 * (matricesCurrentVC(iHalf, jHalf, mm) + matricesCurrentVC(iHalf + 1, jHalf, mm));
          // point on line at phi direction
          matricesCurrentV(i, j, mp1) += 0.25 * ((matricesCurrentVC(iHalf, jHalf, mm) + matricesCurrentVC(iHalf, jHalf, mmPlus)) + (matricesCurrentVC(iHalf + 1, jHalf, mmPlus) + matricesCurrentVC(iHalf + 1, jHalf, mm)));
        }
      }

      for (int j = 1; j < tnZColumn - 1; j += 2) {
        for (int i = 1; i < tnRRow - 1; i += 2) {
          const int iHalf = i * 0.5;
          const int jHalf = j * 0.5;
          matricesCurrentV(i, j, m) += 0.25 * ((matricesCurrentVC(iHalf, jHalf, mm) + matricesCurrentVC(iHalf, jHalf + 1, mm)) + (matricesCurrentVC(iHalf + 1, jHalf, mm) + matricesCurrentVC(iHalf + 1, jHalf + 1, mm)));
          // point at the center at phi direction
          matricesCurrentV(i, j, mp1) += 0.125 * ((matricesCurrentVC(iHalf, jHalf, mm) + matricesCurrentVC(iHalf, jHalf + 1, mm) + matricesCurrentVC(iHalf, jHalf, mmPlus) + matricesCurrentVC(iHalf, jHalf + 1, mmPlus)) +
                                                  (matricesCurrentVC(iHalf + 1, jHalf, mm) + matricesCurrentVC(iHalf + 1, jHalf + 1, mm) + matricesCurrentVC(iHalf + 1, jHalf, mmPlus) + matricesCurrentVC(iHalf + 1, jHalf + 1, mmPlus)));
        }
      }
    }

  } else {
#pragma omp parallel for // no change
    for (int m = 0; m < newPhiSlice; m++) {
      addInterp2D(matricesCurrentV, matricesCurrentVC, tnRRow, tnZColumn, m);
    }
  }
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nz, Nr, Nphi>::addInterp2D(Matrix3D& matricesCurrentV, const Matrix3D& matricesCurrentVC, const int tnRRow, const int tnZColumn, const int iphi) const
{
  for (int j = 2; j < tnZColumn - 1; j += 2) {
    for (int i = 2; i < tnRRow - 1; i += 2) {
      matricesCurrentV(i, j, iphi) = matricesCurrentV(i, j, iphi) + matricesCurrentVC(i * 0.5, j * 0.5, iphi);
    }
  }

  for (int j = 1; j < tnZColumn - 1; j += 2) {
    for (int i = 2; i < tnRRow - 1; i += 2) {
      const int iHalf = 0.5 * i;
      const int jHalf = 0.5 * j;
      matricesCurrentV(i, j, iphi) = matricesCurrentV(i, j, iphi) + 0.5 * (matricesCurrentVC(iHalf, jHalf, iphi) + matricesCurrentVC(iHalf, jHalf + 1, iphi));
    }
  }

  for (int j = 2; j < tnZColumn - 1; j += 2) {
    for (int i = 1; i < tnRRow - 1; i += 2) {
      const int iHalf = 0.5 * i;
      const int jHalf = 0.5 * j;
      matricesCurrentV(i, j, iphi) = matricesCurrentV(i, j, iphi) + 0.5 * (matricesCurrentVC(iHalf, jHalf, iphi) + matricesCurrentVC(iHalf + 1, jHalf, iphi));
    }
  }

  // only if full
  if (o2::tpc::MGParameters::gtType == Full) {
    for (int j = 1; j < tnZColumn - 1; j += 2) {
      for (int i = 1; i < tnRRow - 1; i += 2) {
        const int iHalf = 0.5 * i;
        const int jHalf = 0.5 * j;
        matricesCurrentV(i, j, iphi) = matricesCurrentV(i, j, iphi) + 0.25 * (matricesCurrentVC(iHalf, jHalf, iphi) + matricesCurrentVC(iHalf, jHalf + 1, iphi) + matricesCurrentVC(iHalf + 1, jHalf, iphi) + matricesCurrentVC(iHalf + 1, jHalf + 1, iphi));
      }
    }
  }
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nz, Nr, Nphi>::relax3D(Matrix3D& matricesCurrentV, const Matrix3D& matricesCurrentCharge, const int tnRRow, const int tnZColumn, const int iPhi, const int symmetry, const DataT h2,
                                                      const DataT tempRatioZ, const std::array<DataT, Nr>& coefficient1, const std::array<DataT, Nr>& coefficient2, const std::array<DataT, Nr>& coefficient3, const std::array<DataT, Nr>& coefficient4) const
{
  // Gauss-Seidel (Read Black}
  if (o2::tpc::MGParameters::relaxType == GaussSeidel) {
    // for each slice
    // std::cout<<"X"<<std::endl;
    for (int iPass = 1; iPass <= 2; ++iPass) {
      const int msw = (iPass % 2) ? 1 : 2;
      for (int m = 0; m < iPhi; ++m) {
        const int jsw = ((msw + m) % 2) ? 1 : 2;

        int mp1 = m + 1;
        int signPlus = 1;
        int mm1 = m - 1;
        int signMinus = 1;
        // Reflection symmetry in phi (e.g. symmetry at sector boundaries, or half sectors, etc.)
        if (symmetry == 1) {
          if (mp1 > iPhi - 1) {
            mp1 = iPhi - 2;
          }
          if (mm1 < 0) {
            mm1 = 1;
          }
        }
        // Anti-symmetry in phi
        else if (symmetry == -1) {
          if (mp1 > iPhi - 1) {
            mp1 = iPhi - 2;
            signPlus = -1;
          }
          if (mm1 < 0) {
            mm1 = 1;
            signMinus = -1;
          }
        } else { // No Symmetries in phi, no boundaries, the calculation is continuous across all phi
          if (mp1 > iPhi - 1) {
            mp1 = m + 1 - iPhi;
          }
          if (mm1 < 0) {
            mm1 = m - 1 + iPhi;
          }
        }
        // std::cout<<"Y"<<std::endl;
        int isw = jsw;
        for (int j = 1; j < tnZColumn - 1; j++, isw = 3 - isw) {
          for (int i = isw; i < tnRRow - 1; i += 2) {
            // std::cout<<"i: "<< i << std::endl;
            // std::cout<<"j: "<< j << std::endl;
            // std::cout<<"m: "<< m << std::endl;
            // std::cout<<"mm1: "<< mm1 << std::endl;
            // std::cout<<"mp1: "<< mp1 << std::endl;
            // (matricesCurrentV)(i, j, m) = 11;
            // std::cout<<"matricesCurrentV.nR"<< matricesCurrentV.nR << std::endl;
            // std::cout<<"matricesCurrentV.nR"<< matricesCurrentV.nZ << std::endl;
            // std::cout<<"matricesCurrentV.nR"<< matricesCurrentV.nPhi << std::endl;

            (matricesCurrentV)(i, j, m) = (coefficient2[i] * (matricesCurrentV)(i - 1, j, m) + tempRatioZ * ((matricesCurrentV)(i, j - 1, m) + (matricesCurrentV)(i, j + 1, m)) + coefficient1[i] * (matricesCurrentV)(i + 1, j, m) + coefficient3[i] * (signPlus * (matricesCurrentV)(i, j, mp1) + signMinus * (matricesCurrentV)(i, j, mm1)) + (h2 * (matricesCurrentCharge)(i, j, m))) * coefficient4[i];
            // std::cout<<"--------"<< std::endl;
          } // end cols
        }   // end Nr
      }     // end phi
    }       // end sweep
  } else if (o2::tpc::MGParameters::relaxType == Jacobi) {
    // for each slice
    for (int m = 0; m < iPhi; m++) {
      int mp1 = m + 1;
      int signPlus = 1;
      int mm1 = m - 1;
      int signMinus = 1;

      // Reflection symmetry in phi (e.g. symmetry at sector boundaries, or half sectors, etc.)
      if (symmetry == 1) {
        if (mp1 > iPhi - 1) {
          mp1 = iPhi - 2;
        }
        if (mm1 < 0) {
          mm1 = 1;
        }
      }
      // Anti-symmetry in phi
      else if (symmetry == -1) {
        if (mp1 > iPhi - 1) {
          mp1 = iPhi - 2;
          signPlus = -1;
        }
        if (mm1 < 0) {
          mm1 = 1;
          signMinus = -1;
        }
      } else { // No Symmetries in phi, no boundaries, the calculation is continuous across all phi
        if (mp1 > iPhi - 1) {
          mp1 = m + 1 - iPhi;
        }
        if (mm1 < 0) {
          mm1 = m - 1 + iPhi;
        }
      }
      // Jacobian
      for (int j = 1; j < tnZColumn - 1; j++) {
        for (int i = 1; i < tnRRow - 1; i++) {
          (matricesCurrentV)(i, j, m) = (coefficient2[i] * (matricesCurrentV)(i - 1, j, m) + tempRatioZ * ((matricesCurrentV)(i, j - 1, m) + (matricesCurrentV)(i, j + 1, m)) + coefficient1[i] * (matricesCurrentV)(i + 1, j, m) + coefficient3[i] * (signPlus * (matricesCurrentV)(i, j, mp1) + signMinus * (matricesCurrentV)(i, j, mm1)) + (h2 * (matricesCurrentCharge)(i, j, m))) * coefficient4[i];
        } // end cols
      }   // end Nr
    }     // end phi
  } else {
    // Case weighted Jacobi
    // TODO
  }
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nz, Nr, Nphi>::restrictBoundary3D(Matrix3D& matricesCurrentCharge, const Matrix3D& residue, const int tnRRow, const int tnZColumn, const int newPhiSlice, const int oldPhiSlice) const
{
  // in case of full 3d and the Nphi is also coarsening
  if (2 * newPhiSlice == oldPhiSlice) {
    for (int m = 0, mm = 0; m < newPhiSlice; ++m, mm += 2) {
      // for boundary
      for (int j = 0, jj = 0; j < tnZColumn; ++j, jj += 2) {
        matricesCurrentCharge(0, j, m) = residue(0, jj, mm);
        matricesCurrentCharge(tnRRow - 1, j, m) = residue((tnRRow - 1) * 2, jj, mm);
      }

      // for boundary
      for (int i = 0, ii = 0; i < tnRRow; ++i, ii += 2) {
        matricesCurrentCharge(i, 0, m) = residue(ii, 0, mm);
        matricesCurrentCharge(i, tnZColumn - 1, m) = residue(ii, (tnZColumn - 1) * 2, mm);
      }
    } // end phis
  } else {
    for (int m = 0; m < newPhiSlice; ++m) {
      restrictBoundary2D(matricesCurrentCharge, residue, tnRRow, tnZColumn, m);
    }
  }
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nz, Nr, Nphi>::restrictBoundary2D(Matrix3D& matricesCurrentCharge, const Matrix3D& residue, const int tnRRow, const int tnZColumn, const int iphi) const
{
  // for boundary
  for (int j = 0, jj = 0; j < tnZColumn; ++j, jj += 2) {
    matricesCurrentCharge(0, j, iphi) = residue(0, jj, iphi);
    matricesCurrentCharge(tnRRow - 1, j, iphi) = residue((tnRRow - 1) * 2, jj, iphi);
  }

  // for boundary
  for (int i = 0, ii = 0; i < tnRRow; ++i, ii += 2) {
    matricesCurrentCharge(i, 0, iphi) = residue(ii, 0, iphi);
    matricesCurrentCharge(i, tnZColumn - 1, iphi) = residue(ii, (tnZColumn - 1) * 2, iphi);
  }
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nz, Nr, Nphi>::restrict3D(Matrix3D& matricesCurrentCharge, const Matrix3D& residue, const int tnRRow, const int tnZColumn, const int newPhiSlice, const int oldPhiSlice) const
{
  if (2 * newPhiSlice == oldPhiSlice) {
    int mm = 0;
    for (int m = 0; m < newPhiSlice; m++, mm += 2) {

      // assuming no symmetry
      int mp1 = mm + 1;
      int mm1 = mm - 1;

      if (mp1 > (oldPhiSlice)-1) {
        mp1 = mm + 1 - (oldPhiSlice);
      }
      if (mm1 < 0) {
        mm1 = mm - 1 + (oldPhiSlice);
      }

      for (int i = 1, ii = 2; i < tnRRow - 1; i++, ii += 2) {
        for (int j = 1, jj = 2; j < tnZColumn - 1; j++, jj += 2) {

          // at the same plane
          const int iip1 = ii + 1;
          const int iim1 = ii - 1;
          const int jjp1 = jj + 1;
          const int jjm1 = jj - 1;
          const DataT s1 = residue(iip1, jj, mm) + residue(iim1, jj, mm) + residue(ii, jjp1, mm) + residue(ii, jjm1, mm) + residue(ii, jj, mp1) + residue(ii, jj, mm1);

          const DataT s2 = (residue(iip1, jjp1, mm) + residue(iip1, jjm1, mm) + residue(iip1, jj, mp1) + residue(iip1, jj, mm1)) +
                           (residue(iim1, jjm1, mm) + residue(iim1, jjp1, mm) + residue(iim1, jj, mp1) + residue(iim1, jj, mm1)) +
                           residue(ii, jjm1, mp1) + residue(ii, jjp1, mm1) + residue(ii, jjm1, mm1) + residue(ii, jjp1, mp1);

          const DataT s3 = (residue(iip1, jjp1, mp1) + residue(iip1, jjm1, mp1) + residue(iip1, jjp1, mm1) + residue(iip1, jjm1, mm1)) +
                           (residue(iim1, jjm1, mm1) + residue(iim1, jjp1, mm1) + residue(iim1, jjm1, mp1) + residue(iim1, jjp1, mp1));

          matricesCurrentCharge(i, j, m) = 0.125 * residue(ii, jj, mm) + 0.0625 * s1 + 0.03125 * s2 + 0.015625 * s3;
        } // end cols
      }   // end Nr

      // for boundary
      for (int j = 0, jj = 0; j < tnZColumn; j++, jj += 2) {
        matricesCurrentCharge(0, j, m) = residue(0, jj, mm);
        matricesCurrentCharge(tnRRow - 1, j, m) = residue((tnRRow - 1) * 2, jj, mm);
      }

      // for boundary
      for (int i = 0, ii = 0; i < tnRRow; i++, ii += 2) {
        matricesCurrentCharge(i, 0, m) = residue(ii, 0, mm);
        matricesCurrentCharge(i, tnZColumn - 1, m) = residue(ii, (tnZColumn - 1) * 2, mm);
      }
    } // end phis

  } else {
    for (int m = 0; m < newPhiSlice; m++) {
      restrict2D(matricesCurrentCharge, residue, tnRRow, tnZColumn, m);
    }
  }
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nz, Nr, Nphi>::restrict2D(Matrix3D& matricesCurrentCharge, const Matrix3D& residue, const int tnRRow, const int tnZColumn, const int iphi) const
{
  for (int i = 1, ii = 2; i < tnRRow - 1; i++, ii += 2) {
    for (int j = 1, jj = 2; j < tnZColumn - 1; j++, jj += 2) {
      const int iip1 = ii + 1;
      const int iim1 = ii - 1;
      const int jjp1 = jj + 1;
      const int jjm1 = jj - 1;
      if (o2::tpc::MGParameters::gtType == Half) {
        // half
        matricesCurrentCharge(i, j, iphi) = 0.5 * residue(ii, jj, iphi) + 0.125 * (residue(iip1, jj, iphi) + residue(iim1, jj, iphi) + residue(ii, jjp1, iphi) + residue(ii, jjm1, iphi));
      } else if (o2::tpc::MGParameters::gtType == Full) {
        matricesCurrentCharge(i, j, iphi) = 0.25 * residue(ii, jj, iphi) + 0.125 * (residue(iip1, jj, iphi) + residue(iim1, jj, iphi) + residue(ii, jjp1, iphi) + residue(ii, jjm1, iphi)) +
                                            0.0625 * (residue(iip1, jjp1, iphi) + residue(iim1, jjp1, iphi) + residue(iip1, jjm1, iphi) + residue(iim1, jjm1, iphi));
      }
    } // end cols
  }   // end Nr

  // boundary
  // for boundary
  for (int j = 0, jj = 0; j < tnZColumn; j++, jj += 2) {
    matricesCurrentCharge(0, j, iphi) = residue(0, jj, iphi);
    matricesCurrentCharge(tnRRow - 1, j, iphi) = residue((tnRRow - 1) * 2, jj, iphi);
  }

  // for boundary
  for (int i = 0, ii = 0; i < tnRRow; i++, ii += 2) {
    matricesCurrentCharge(i, 0, iphi) = residue(ii, 0, iphi);
    matricesCurrentCharge(i, tnZColumn - 1, iphi) = residue(ii, (tnZColumn - 1) * 2, iphi);
  }
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
bool O2TPCPoissonSolver<DataT, Nz, Nr, Nphi>::isPowerOfTwo(int i) const
{
  int j = 0;
  while (i > 0) {
    j += (i & 1);
    i = (i >> 1);
  }
  if (j == 1) {
    return true;
  }
  return false;
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
DataT O2TPCPoissonSolver<DataT, Nz, Nr, Nphi>::getConvergenceError(const Matrix3D& matricesCurrentV, Matrix3D& prevArrayV) const
{
  std::vector<DataT> errorArr(prevArrayV.nPhi);

  // subtract the two matrices
  std::transform(prevArrayV.storage.begin(), prevArrayV.storage.end(), matricesCurrentV.storage.begin(), prevArrayV.storage.begin(), std::minus<DataT>());

#pragma omp parallel for // parallising this loop is possible - but using more than 2 cores makes it slower -
  for (unsigned int m = 0; m < prevArrayV.nPhi; m++) {
    // square each entry in the vector and sum them up
    const auto phiStep = prevArrayV.nR * prevArrayV.nZ; // number of points in one phi slice
    const auto start = prevArrayV.storage.begin() + m * phiStep;
    const auto end = start + phiStep;
    errorArr[m] = std::inner_product(start, end, start, 0.); // inner product "Sum (matrix[a]*matrix[a])"
  }
  // return largest error
  return *std::max_element(std::begin(errorArr), std::end(errorArr));
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nz, Nr, Nphi>::calcCoefficients(unsigned int from, unsigned int to, const DataT h, const DataT tempRatioZ, const DataT tempRatioPhi, std::array<DataT, Nr>& coefficient1, std::array<DataT, Nr>& coefficient2, std::array<DataT, Nr>& coefficient3, std::array<DataT, Nr>& coefficient4) const
{
  for (unsigned int i = from; i < to; ++i) {
    const DataT radiusInv = 1. / (TPCParameters<DataT>::IFCRADIUS + i * h);
    const DataT hRadiusTmp = h * 0.5 * radiusInv;
    coefficient1[i] = 1.0 + hRadiusTmp;
    coefficient2[i] = 1.0 - hRadiusTmp;
    coefficient3[i] = tempRatioPhi * radiusInv * radiusInv;
    coefficient4[i] = 0.5 / (1.0 + tempRatioZ + coefficient3[i]);
  }
}

template class o2::tpc::O2TPCPoissonSolver<float, 17, 17, 90>;
template class o2::tpc::O2TPCPoissonSolver<double, 17, 17, 90>;
template class o2::tpc::O2TPCPoissonSolver<float, 65, 65, 90>;
template class o2::tpc::O2TPCPoissonSolver<float, 129, 129, 180>;
template class o2::tpc::O2TPCPoissonSolver<double, 129, 129, 180>;
