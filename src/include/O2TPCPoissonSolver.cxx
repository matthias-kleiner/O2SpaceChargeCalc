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
  poissonMultiGrid3D2D(matricesV, matricesCharge, symmetry);
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
    relax3D(tvArrayV[nLoop - 1], tvChargeFMG[nLoop - 1], tnRRow, tnZColumn, symmetry, h2, tempRatioZ, coefficient1, coefficient2, coefficient3, coefficient4);

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
      relax3D(tvArrayV[count - 1], tvCharge[count - 1], tnRRow, tnZColumn, symmetry, h2, tempRatioZ, coefficient1, coefficient2, coefficient3, coefficient4);
    } // end pre smoothing

    // 2) Residue calculation
    residue3D(tvResidue[count - 1], tvArrayV[count - 1], tvCharge[count - 1], tnRRow, tnZColumn, symmetry, ih2, tempRatioZ, coefficient1, coefficient2, coefficient3, inverseCoefficient4);

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
  relax3D(tvArrayV[gridTo - 1], tvCharge[gridTo - 1], tnRRow, tnZColumn, symmetry, h2, tempRatioZ, coefficient1, coefficient2, coefficient3, coefficient4);

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
      relax3D(tvArrayV[count - 1], tvCharge[count - 1], tnRRow, tnZColumn, symmetry, h2, tempRatioZ, coefficient1, coefficient2, coefficient3, coefficient4);
    } // end post smoothing
  }
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nz, Nr, Nphi>::residue3D(Matrix3D& residue, const Matrix3D& matricesCurrentV, const Matrix3D& matricesCurrentCharge, const int tnRRow, const int tnZColumn, const int symmetry,
                                                        const DataT ih2, const DataT tempRatioZ, const std::array<DataT, Nr>& coefficient1, const std::array<DataT, Nr>& coefficient2, const std::array<DataT, Nr>& coefficient3, const std::array<DataT, Nr>& inverseCoefficient4) const
{
#pragma omp parallel for // parallising this loop is possible - but using more than 2 cores makes it slower -
  for (int m = 0; m < Nphi; ++m) {
    int mp1 = m + 1;
    int signPlus = 1;
    int mm1 = m - 1;
    int signMinus = 1;

    // Reflection symmetry in phi (e.g. symmetry at sector boundaries, or half sectors, etc.)
    if (symmetry == 1) {
      if (mp1 > Nphi - 1) {
        mp1 = Nphi - 2;
      }
      if (mm1 < 0) {
        mm1 = 1;
      }
    }
    // Anti-symmetry in phi
    else if (symmetry == -1) {
      if (mp1 > Nphi - 1) {
        mp1 = Nphi - 2;
        signPlus = -1;
      }
      if (mm1 < 0) {
        mm1 = 1;
        signMinus = -1;
      }
    } else { // No Symmetries in phi, no boundaries, the calculation is continuous across all phi
      if (mp1 > Nphi - 1) {
        mp1 = m + 1 - Nphi;
      }
      if (mm1 < 0) {
        mm1 = m - 1 + Nphi;
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
void O2TPCPoissonSolver<DataT, Nz, Nr, Nphi>::relax3D(Matrix3D& matricesCurrentV, const Matrix3D& matricesCurrentCharge, const int tnRRow, const int tnZColumn, const int symmetry, const DataT h2,
                                                      const DataT tempRatioZ, const std::array<DataT, Nr>& coefficient1, const std::array<DataT, Nr>& coefficient2, const std::array<DataT, Nr>& coefficient3, const std::array<DataT, Nr>& coefficient4) const
{
  // Gauss-Seidel (Read Black}
  if (o2::tpc::MGParameters::relaxType == GaussSeidel) {
    // for each slice
    for (int iPass = 1; iPass <= 2; ++iPass) {
      const int msw = (iPass % 2) ? 1 : 2;
      for (int m = 0; m < Nphi; ++m) {
        const int jsw = ((msw + m) % 2) ? 1 : 2;

        int mp1 = m + 1;
        int signPlus = 1;
        int mm1 = m - 1;
        int signMinus = 1;
        // Reflection symmetry in phi (e.g. symmetry at sector boundaries, or half sectors, etc.)
        if (symmetry == 1) {
          if (mp1 > Nphi - 1) {
            mp1 = Nphi - 2;
          }
          if (mm1 < 0) {
            mm1 = 1;
          }
        }
        // Anti-symmetry in phi
        else if (symmetry == -1) {
          if (mp1 > Nphi - 1) {
            mp1 = Nphi - 2;
            signPlus = -1;
          }
          if (mm1 < 0) {
            mm1 = 1;
            signMinus = -1;
          }
        } else { // No Symmetries in phi, no boundaries, the calculation is continuous across all phi
          if (mp1 > Nphi - 1) {
            mp1 = m + 1 - Nphi;
          }
          if (mm1 < 0) {
            mm1 = m - 1 + Nphi;
          }
        }

        int isw = jsw;
        for (int j = 1; j < tnZColumn - 1; j++, isw = 3 - isw) {
          for (int i = isw; i < tnRRow - 1; i += 2) {
            (matricesCurrentV)(i, j, m) = (coefficient2[i] * (matricesCurrentV)(i - 1, j, m) + tempRatioZ * ((matricesCurrentV)(i, j - 1, m) + (matricesCurrentV)(i, j + 1, m)) + coefficient1[i] * (matricesCurrentV)(i + 1, j, m) + coefficient3[i] * (signPlus * (matricesCurrentV)(i, j, mp1) + signMinus * (matricesCurrentV)(i, j, mm1)) + (h2 * (matricesCurrentCharge)(i, j, m))) * coefficient4[i];
          } // end cols
        }   // end Nr
      }     // end phi
    }       // end sweep
  } else if (o2::tpc::MGParameters::relaxType == Jacobi) {
    // for each slice
    for (int m = 0; m < Nphi; m++) {
      int mp1 = m + 1;
      int signPlus = 1;
      int mm1 = m - 1;
      int signMinus = 1;

      // Reflection symmetry in phi (e.g. symmetry at sector boundaries, or half sectors, etc.)
      if (symmetry == 1) {
        if (mp1 > Nphi - 1) {
          mp1 = Nphi - 2;
        }
        if (mm1 < 0) {
          mm1 = 1;
        }
      }
      // Anti-symmetry in phi
      else if (symmetry == -1) {
        if (mp1 > Nphi - 1) {
          mp1 = Nphi - 2;
          signPlus = -1;
        }
        if (mm1 < 0) {
          mm1 = 1;
          signMinus = -1;
        }
      } else { // No Symmetries in phi, no boundaries, the calculation is continuous across all phi
        if (mp1 > Nphi - 1) {
          mp1 = m + 1 - Nphi;
        }
        if (mm1 < 0) {
          mm1 = m - 1 + Nphi;
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
  DataT errorArr[Nphi]{};

  // subtract the two matrices
  std::transform(prevArrayV.storage.begin(), prevArrayV.storage.end(), matricesCurrentV.storage.begin(), prevArrayV.storage.begin(), std::minus<DataT>());

#pragma omp parallel for // parallising this loop is possible - but using more than 2 cores makes it slower -
  for (unsigned int m = 0; m < Nphi; m++) {
    // square each entry in the vector and sum them up
    const auto phiStep = prevArrayV.nX * prevArrayV.nY; // number of points in one phi slice
    const auto start = prevArrayV.storage.begin() + m * phiStep;
    const auto end = start + phiStep;
    errorArr[m] = std::inner_product(start, end, start, 0.); // inner product "Sum (matrix[a]*matrix[a])"
  }
  // return largest error
  return *std::max_element(errorArr, errorArr + Nphi);
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
