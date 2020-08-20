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
/// \author Rifki Sadikin <rifki.sadikin@cern.ch>, Indonesian Institute of Sciences
/// \date Nov 20, 2017

#include <TMath.h>
#include "O2TPCPoissonSolver.h"

#include <iostream>
#include <numeric>

// #include "RegularGrid3D.h"

#include "Rtypes.h"
/// \cond CLASSIMP
// ClassImp(O2TPCPoissonSolver);
// templateClassImp(O2TPCPoissonSolver);
/// \endcond

/// Provides poisson solver in Cylindrical 3D (TPC geometry)
///
/// Strategy based on parameter settings (fStrategy and fMgParameters)provided
/// * Cascaded multi grid with S.O.R
/// * Geometric MultiGrid
///		* Cycles: V, W, Full
///		* Relaxation: Jacobi, Weighted-Jacobi, Gauss-Seidel
///		* Grid transfer operators: Full, Half
/// * Spectral Methods (TODO)
///
/// \param matricesV TMatrixD** potential in 3D matrix
/// \param matricesCharge TMatrixD** charge density in 3D matrix (side effect)
/// \param Nr int number of Nr in the r direction of TPC
/// \param Nz int number of Nz in z direction of TPC
/// \param Nphi int number of Nphi in phi direction of T{C
/// \param maxIteration int maximum iteration for relaxation method
/// \param symmetry int symmetry or not
///
/// \pre Charge density distribution in **matricesCharge** is known and boundary values for **matricesV** are set
/// \post Numerical solution for potential distribution is calculated and stored in each rod at **matricesV**
template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nr, Nz, Nphi>::PoissonSolver3D(DataContainer& matricesV, DataContainer& matricesCharge,
                                                              int maxIteration,
                                                              int symmetry)
{
  PoissonMultiGrid3D2D(matricesV, matricesCharge, symmetry);
}

/// 3D - Solve Poisson's Equation in 3D by MultiGrid with constant phi slices
///
///    NOTE: In order for this algorithm to work, the number of Nr and Nz must be a power of 2 plus one.
///    The number of Nr and Z Column can be different.
///
///    R Row       ==  2**M + 1
///    Z Column  ==  2**N + 1
///    Phi Slice  ==  Arbitrary but greater than 3
///
///		 Solving: \f$  \nabla^{2}V(r,\phi,z) = - f(r,\phi,z) \f$
///
/// Algorithm for MultiGrid Full Cycle (FMG)
/// - Relax on the coarsest grid
/// - Do from coarsest to finest
///     - Interpolate potential from coarse -> fine
///   - Do V-Cycle to the current coarse level to the coarsest
///   - Stop if converged
///
/// DeltaPhi in Radians
/// \param matricesV TMatrixD** potential in 3D matrix \f$ V(r,\phi,z) \f$
/// \param matricesCharge TMatrixD** charge density in 3D matrix (side effect) \f$ - f(r,\phi,z) \f$
/// \param Nr int number of Nr in the r direction of TPC
/// \param Nz int number of Nz in z direction of TPC
/// \param Nphi int number of Nphi in phi direction of T{C
/// \param maxIteration int maximum iteration for relaxation method (NOT USED)
/// \param symmetry int symmetry (TODO for symmetry = 1)
//
///    SYMMETRY = 0 if no phi symmetries, and no phi boundary condition
///    = 1 if we have reflection symmetry at the boundaries (eg. sector symmetry or half sector symmetries).
///
template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nr, Nz, Nphi>::PoissonMultiGrid3D2D(DataContainer& matricesV, DataContainer& matricesCharge, int symmetry)
{
  const DataT ratioPhi = GRIDSPACINGR * GRIDSPACINGR / (GRIDSPACINGPHI * GRIDSPACINGPHI); // ratio_{phi} = gridSize_{r} / gridSize_{phi}
  const DataT ratioZ = GRIDSPACINGR * GRIDSPACINGR / (GRIDSPACINGZ * GRIDSPACINGZ);       // ratio_{Z} = gridSize_{r} / gridSize_{z}

  Info("PoissonMultiGrid3D2D", "%s", Form("in Poisson Solver 3D multiGrid semi coarsening Nr=%lu, cols=%lu, Nphi=%lu \n", Nr, Nz, Nphi));

  // Check that the number of Nr and Nz is suitable for a binary expansion
  if (!IsPowerOfTwo((Nr - 1))) {
    Error("PoissonMultiGrid3D2D", "Poisson3DMultiGrid - Error in the number of Nr. Must be 2**M + 1");
    return;
  }
  if (!IsPowerOfTwo((Nz - 1))) {
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

  int nLoop = std::max(nGridRow, nGridCol); // Calculate the number of nLoop for the binary expansion
  nLoop = (nLoop > fMgParameters.maxLoop) ? fMgParameters.maxLoop : nLoop;
  int iOne = 1; // index i in gridSize r (original)
  int jOne = 1; // index j in gridSize z (original)

  std::vector<Matrix3D> tvArrayV(nLoop);     // potential <--> error
  std::vector<Matrix3D> tvChargeFMG(nLoop);  // charge is restricted in full multiGrid
  std::vector<Matrix3D> tvCharge(nLoop);     // charge <--> residue
  std::vector<Matrix3D> tvPrevArrayV(nLoop); // error calculation
  std::vector<Matrix3D> tvResidue(nLoop);    // residue calculation

  for (int count = 1; count <= nLoop; count++) {
    const int tnRRow = iOne == 1 ? Nr : Nr / iOne + 1;
    const int tnZColumn = jOne == 1 ? Nz : Nz / jOne + 1;
    tvResidue[count - 1].resize(tnRRow, tnZColumn, Nphi);
    tvPrevArrayV[count - 1].resize(tnRRow, tnZColumn, Nphi);

    // memory for the finest grid is from parameters
    tvChargeFMG[count - 1].resize(tnRRow, tnZColumn, Nphi);
    tvArrayV[count - 1].resize(tnRRow, tnZColumn, Nphi);
    tvCharge[count - 1].resize(tnRRow, tnZColumn, Nphi);

    if (count == 1) {
      std::copy(matricesCharge.begin(), matricesCharge.end(), tvChargeFMG[count - 1].mData.data());
      // std::copy(matricesCharge.begin(), matricesCharge.end(), tvCharge[count - 1].mData.data());
      std::copy(matricesV.begin(), matricesV.end(), tvArrayV[count - 1].mData.data());
    } else {
      Restrict3D(tvChargeFMG[count - 1], tvChargeFMG[count - 2], tnRRow, tnZColumn, Nphi, Nphi);
      RestrictBoundary3D(tvArrayV[count - 1], tvArrayV[count - 2], tnRRow, tnZColumn, Nphi, Nphi);
    }
    iOne = 2 * iOne; // doubling
    jOne = 2 * jOne; // doubling
  }

  std::vector<float> coefficient1(Nr);        // coefficient1(Nr) for storing (1 + h_{r}/2r_{i}) from central differences in r direction
  std::vector<float> coefficient2(Nr);        // coefficient2(Nr) for storing (1 + h_{r}/2r_{i}) from central differences in r direction
  std::vector<float> coefficient3(Nr);        // coefficient3(Nr) for storing (1/r_{i}^2) from central differences in phi direction
  std::vector<float> coefficient4(Nr);        // coefficient4(Nr) for storing  1/2
  std::vector<float> inverseCoefficient4(Nr); // inverse of coefficient4(Nr)

  // Case full multi grid (FMG)
  if (fMgParameters.cycleType == kFCycle) {
    // 1) Relax on the coarsest grid
    iOne = iOne / 2;
    jOne = jOne / 2;
    int tnRRow = iOne == 1 ? Nr : Nr / iOne + 1;
    int tnZColumn = jOne == 1 ? Nz : Nz / jOne + 1;

    const DataT h = GRIDSPACINGR * iOne;
    const DataT h2 = h * h;

    const DataT iOne2 = iOne * iOne;
    const DataT tempRatioPhi = ratioPhi * iOne2; // Used tobe divided by ( m_one * m_one ) when m_one was != 1
    const DataT tempRatioZ = ratioZ * iOne2 / (jOne * jOne);

    for (int i = 1; i < tnRRow - 1; i++) {
      const DataT radius = O2TPCPoissonSolver<DataT, Nr, Nz, Nphi>::fgkIFCRadius + i * h;
      coefficient1[i] = 1.0 + h / (2 * radius);
      coefficient2[i] = 1.0 - h / (2 * radius);
      coefficient3[i] = tempRatioPhi / (radius * radius);
      coefficient4[i] = 0.5 / (1.0 + tempRatioZ + coefficient3[i]);
    }
    // relax on the coarsest level
    Relax3D(tvArrayV[nLoop - 1], tvChargeFMG[nLoop - 1], tnRRow, tnZColumn, symmetry, h2, tempRatioZ, coefficient1, coefficient2, coefficient3, coefficient4);

    // 2) Do multiGrid v-cycle from coarsest to finest
    for (int count = nLoop - 2; count >= 0; count--) {
      // move to finer grid
      iOne = iOne / 2;
      jOne = jOne / 2;
      tnRRow = iOne == 1 ? Nr : Nr / iOne + 1;
      tnZColumn = jOne == 1 ? Nz : Nz / jOne + 1;
      // 2) a) Interpolate potential for h -> 2h (coarse -> fine)
      Interp3D(tvArrayV[count], tvArrayV[count + 1], tnRRow, tnZColumn, Nphi, Nphi);

      // 2) c) Copy the restricted charge to charge for calculation
      tvCharge[count] = tvChargeFMG[count]; //copy
      // 2) c) Do V cycle fMgParameters.nMGCycle times at most

      for (int mgCycle = 0; mgCycle < fMgParameters.nMGCycle; mgCycle++) {
        // Copy the potential to temp array for convergence calculation
        tvPrevArrayV[count] = tvArrayV[count]; //copy

        // 2) c) i) Call V cycle from grid count+1 (current fine level) to nLoop (coarsest)
        VCycle3D2D(symmetry, count + 1, nLoop, fMgParameters.nPre, fMgParameters.nPost,
                   ratioZ, ratioPhi, tvArrayV, tvCharge, tvResidue, coefficient1, coefficient2, coefficient3,
                   coefficient4, inverseCoefficient4);

        const DataT convergenceError = GetConvergenceError(tvArrayV[count], tvPrevArrayV[count]);

        if (count == 0) {
          (*fErrorConvergenceNormInf)(mgCycle) = convergenceError;
          // (*fError)(mgCycle) = GetExactError(matricesV, tvPrevArrayV[count], Nphi);
        }

        /// if already converge just break move to finer grid
        if (convergenceError <= fgConvergenceError) {
          fIterations = mgCycle + 1;
          break;
        }
      }
    }
  } // Case V multi grid (VMG)

  std::move(tvArrayV[0].mData.begin(), tvArrayV[0].mData.end(), matricesV.data());
}

/// VCycle 3D2D, V Cycle 3D in multiGrid with constant Nphi
/// fine-->coarsest-->fine, propagating the residue to correct initial guess of V
///
/// Algorithm:
///
///    NOTE: In order for this algorithm to work, the number of Nr and Nz must be a power of 2 plus one.
///    The number of Nr and Z Column can be different.
///
///    R Row       ==  2**M + 1
///    Z Column    ==  2**N + 1
///    Phi Slice  ==  Arbitrary but greater than 3
///
///    DeltaPhi in Radians
/// \param Nr int number of grid in Nr (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
/// \param Nz int number of grid in Nz (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1
/// \param gridFrom const int finest level of grid
/// \param gridTo const int coarsest level of grid
/// \param nPre const int number of smoothing before coarsening
/// \param nPost const int number of smoothing after coarsening
/// \param gridSizeR const Float_t grid size in r direction (OPTION,  recalculate)
/// \param ratio const Float_t ratio between square of grid r and grid z (OPTION,  recalculate)
/// \param tvArrayV vector<TMatrixD *> vector of V potential in different grids
/// \param tvCharge vector<TMatrixD *> vector of charge distribution in different grids
/// \param tvResidue vector<TMatrixD *> vector of residue calculation in different grids
/// \param coefficient1 std::vector<float>& coefficient for relaxation (r direction)
/// \param coefficient2 std::vector<float>& coefficient for relaxation (r direction)
/// \param coefficient3 std::vector<float>& coefficient for relaxation (ratio r/z)
/// \param coefficient4 std::vector<float>& coefficient for relaxation (ratio for grid_r)
/// \param inverseCoefficient4 std::vector<float>& coefficient for relaxation (inverse coefficient4)
///
template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nr, Nz, Nphi>::VCycle3D2D(const int symmetry,
                                                         const int gridFrom, const int gridTo, const int nPre, const int nPost,
                                                         const DataT ratioZ, const DataT ratioPhi,
                                                         std::vector<Matrix3D>& tvArrayV, std::vector<Matrix3D>& tvCharge,
                                                         std::vector<Matrix3D>& tvResidue, std::vector<float>& coefficient1,
                                                         std::vector<float>& coefficient2, std::vector<float>& coefficient3,
                                                         std::vector<float>& coefficient4,
                                                         std::vector<float>& inverseCoefficient4)
{
  int iOne = 1 << (gridFrom - 1);
  int jOne = 1 << (gridFrom - 1);
  int tnRRow = iOne == 1 ? Nr : Nr / iOne + 1;
  int tnZColumn = jOne == 1 ? Nz : Nz / jOne + 1;

  for (int count = gridFrom; count <= gridTo - 1; ++count) {
    const DataT h = GRIDSPACINGR * iOne;
    const DataT h2 = h * h;
    const DataT ih2 = 1.0 / h2;
    const int iOne2 = iOne * iOne;
    const DataT tempRatioPhi = ratioPhi * iOne2; // Used tobe divided by ( m_one * m_one ) when m_one was != 1
    const DataT tempRatioZ = ratioZ * iOne2 / (jOne * jOne);

    for (int i = 1; i < tnRRow - 1; ++i) {
      const DataT radius = O2TPCPoissonSolver<DataT, Nr, Nz, Nphi>::fgkIFCRadius + i * h;
      coefficient1[i] = 1.0 + h / (2 * radius);
      coefficient2[i] = 1.0 - h / (2 * radius);
      coefficient3[i] = tempRatioPhi / (radius * radius);
      coefficient4[i] = 0.5 / (1.0 + tempRatioZ + coefficient3[i]);
      inverseCoefficient4[i] = 1.0 / coefficient4[i];
    }

    //Info("VCycle3D2D","Before Pre-smoothing");
    // 1) Pre-Smoothing: Gauss-Seidel Relaxation or Jacobi
    for (int jPre = 1; jPre <= nPre; ++jPre) {
      Relax3D(tvArrayV[count - 1], tvCharge[count - 1], tnRRow, tnZColumn, symmetry, h2, tempRatioZ, coefficient1, coefficient2, coefficient3, coefficient4);
    } // end pre smoothing

    // 2) Residue calculation
    Residue3D(tvResidue[count - 1], tvArrayV[count - 1], tvCharge[count - 1], tnRRow, tnZColumn, symmetry, ih2, tempRatioZ, coefficient1, coefficient2, coefficient3, inverseCoefficient4);

    iOne = 2 * iOne;
    jOne = 2 * jOne;
    tnRRow = iOne == 1 ? Nr : Nr / iOne + 1;
    tnZColumn = jOne == 1 ? Nz : Nz / jOne + 1;

    //3) Restriction
    Restrict3D(tvCharge[count], tvResidue[count - 1], tnRRow, tnZColumn, Nphi, Nphi);

    //4) Zeroing coarser V
    std::fill(tvArrayV[count].mData.begin(), tvArrayV[count].mData.end(), 0);
  }

  // coarsest grid
  const DataT h = GRIDSPACINGR * iOne;
  const DataT h2 = h * h;

  const int iOne2 = iOne * iOne;
  const DataT tempRatioPhi = ratioPhi * iOne2;
  const DataT tempRatioZ = ratioZ * iOne2 / (jOne * jOne);

  for (int i = 1; i < tnRRow - 1; ++i) {
    const DataT radius = fgkIFCRadius + i * h;
    coefficient1[i] = 1.0 + h / (2 * radius);
    coefficient2[i] = 1.0 - h / (2 * radius);
    coefficient3[i] = tempRatioPhi / (radius * radius);
    coefficient4[i] = 0.5 / (1.0 + tempRatioZ + coefficient3[i]);
  }

  // 3) Relax on the coarsest grid
  Relax3D(tvArrayV[gridTo - 1], tvCharge[gridTo - 1], tnRRow, tnZColumn, symmetry, h2, tempRatioZ, coefficient1, coefficient2, coefficient3, coefficient4);

  // back to fine
  for (int count = gridTo - 1; count >= gridFrom; --count) {
    iOne = iOne / 2;
    jOne = jOne / 2;

    tnRRow = iOne == 1 ? Nr : Nr / iOne + 1;
    tnZColumn = jOne == 1 ? Nz : Nz / jOne + 1;

    const DataT h = GRIDSPACINGR * iOne;
    const DataT h2 = h * h;

    const int iOne2 = iOne * iOne;
    const DataT tempRatioPhi = ratioPhi * iOne2; // Used tobe divided by ( m_one * m_one ) when m_one was != 1
    const DataT tempRatioZ = ratioZ * iOne2 / (jOne * jOne);

    // 4) Interpolation/Prolongation
    AddInterp3D(tvArrayV[count - 1], tvArrayV[count], tnRRow, tnZColumn, Nphi, Nphi);

    for (int i = 1; i < tnRRow - 1; ++i) {
      const DataT radius = O2TPCPoissonSolver<DataT, Nr, Nz, Nphi>::fgkIFCRadius + i * h;
      coefficient1[i] = 1.0 + h / (2 * radius);
      coefficient2[i] = 1.0 - h / (2 * radius);
      coefficient3[i] = tempRatioPhi / (radius * radius);
      coefficient4[i] = 0.5 / (1.0 + tempRatioZ + coefficient3[i]);
    }

    // 5) Post-Smoothing: Gauss-Seidel Relaxation
    for (int jPost = 1; jPost <= nPost; ++jPost) {
      Relax3D(tvArrayV[count - 1], tvCharge[count - 1], tnRRow, tnZColumn, symmetry, h2, tempRatioZ, coefficient1, coefficient2, coefficient3, coefficient4);
    } // end post smoothing
  }
}

/// Residue3D
///
///    Compute residue from V(.) where V(.) is numerical potential and f(.).
///		 residue used 7 stencil in cylindrical coordinate
///
/// Using the following equations
/// \f$ U_{i,j,k} = (1 + \frac{1}{r_{i}h_{r}}) U_{i+1,j,k}  + (1 - \frac{1}{r_{i}h_{r}}) U_{i+1,j,k}  \f$
///
/// \param residue TMatrixD** residue in 3D (matrices of matrix)
/// \param matricesCurrentV TMatrixD** potential in 3D (matrices of matrix)
/// \param matricesCurrentCharge TMatrixD** charge in 3D
/// \param Nr const int number of Nr in the r direction of TPC
/// \param Nz const int number of Nz in z direction of TPC
/// \param Nphi const int number of Nphi in phi direction of TPC
/// \param symmetry const int is the cylinder has symmetry
/// \param ih2 const Float_t \f$ 1/ h_{r}^{2} \f$
/// \param tempRatioZ const Float_t ration between grid size in z-direction and r-direction
/// \param coefficient1 std::vector<float> coefficient for \f$  V_{x+1,y,z} \f$
/// \param coefficient2 std::vector<float> coefficient for \f$  V_{x-1,y,z} \f$
/// \param coefficient3 std::vector<float> coefficient for z
/// \param inverseCoefficient4 std::vector<float> inverse coefficient for f(r,\phi,z)
///
template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nr, Nz, Nphi>::Residue3D(Matrix3D& residue, Matrix3D& matricesCurrentV, Matrix3D& matricesCurrentCharge, const int tnRRow, const int tnZColumn, const int symmetry,
                                                        const DataT ih2, const DataT tempRatioZ, std::vector<float>& coefficient1, std::vector<float>& coefficient2,
                                                        std::vector<float>& coefficient3, std::vector<float>& inverseCoefficient4)
{
  for (int m = 0; m < Nphi; ++m) {
    int mPlus = m + 1;
    int signPlus = 1;
    int mMinus = m - 1;
    int signMinus = 1;

    // Reflection symmetry in phi (e.g. symmetry at sector boundaries, or half sectors, etc.)
    if (symmetry == 1) {
      if (mPlus > Nphi - 1) {
        mPlus = Nphi - 2;
      }
      if (mMinus < 0) {
        mMinus = 1;
      }
    }
    // Anti-symmetry in phi
    else if (symmetry == -1) {
      if (mPlus > Nphi - 1) {
        mPlus = Nphi - 2;
        signPlus = -1;
      }
      if (mMinus < 0) {
        mMinus = 1;
        signMinus = -1;
      }
    } else { // No Symmetries in phi, no boundaries, the calculation is continuous across all phi
      if (mPlus > Nphi - 1) {
        mPlus = m + 1 - Nphi;
      }
      if (mMinus < 0) {
        mMinus = m - 1 + Nphi;
      }
    }

    for (int j = 1; j < tnZColumn - 1; ++j) {
      for (int i = 1; i < tnRRow - 1; ++i) {
        residue(i, j, m) =
          ih2 * (coefficient2[i] * matricesCurrentV(i - 1, j, m) + tempRatioZ * (matricesCurrentV(i, j - 1, m) + matricesCurrentV(i, j + 1, m)) + coefficient1[i] * matricesCurrentV(i + 1, j, m) +
          coefficient3[i] * (signPlus * matricesCurrentV(i, j, mPlus) + signMinus * matricesCurrentV(i, j, mMinus)) - inverseCoefficient4[i] * matricesCurrentV(i, j, m)) + matricesCurrentCharge(i, j, m);
      } // end cols
    }   // end Nr
  }
}

/// Interpolation/Prolongation in 3D
///
/// Interpolation is a map from coarse grid (h) to fine grid (2h)
///
/// In case of 3D
/// Full weighting:
/// \f[ (R u)_{i,j,k} = \frac{1}{2} u_{2i,2j,2k} + \frac{1}{4} S_{1} + \frac{1}{8} S_{2} + \frac{1}{16} S_{3}\f]
///
///
/// Restriction in all direction r-phi-z
/// restriction in phi only if oldPhi == 2*newPhi
/// \param matricesCurrentV TMatrixD** finer grid h
/// \param curArrayCV TMatrixD ** coarse grid 2h
/// \param tnRRow int number of grid in Nr (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
/// \param tnZColumn int number of grid in Nz (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1
/// \param newPhiSlice int number of Nphi (in phi-direction) for coarser grid
/// \param oldPhiSlice int number of Nphi (in phi-direction) for finer grid
///
template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nr, Nz, Nphi>::Interp3D(Matrix3D& matricesCurrentV, Matrix3D& matricesCurrentVC, const int tnRRow,
                                                       const int tnZColumn, const int newPhiSlice, const int oldPhiSlice)
{

  // Do restrict 2 D for each slice
  if (newPhiSlice == 2 * oldPhiSlice) {
    int mPlus, mmPlus;
    int mm = 0;

    for (int m = 0; m < newPhiSlice; m += 2) {
      // assuming no symmetry
      mm = m / 2;
      mmPlus = mm + 1;
      mPlus = m + 1;

      // round
      if (mmPlus > (oldPhiSlice)-1) {
        mmPlus = mm + 1 - (oldPhiSlice);
      }
      if (mPlus > (newPhiSlice)-1) {
        mPlus = m + 1 - (newPhiSlice);
      }

      for (int j = 2; j < tnZColumn - 1; j += 2) {
        for (int i = 2; i < tnRRow - 1; i += 2) {
          matricesCurrentV(i, j, m) = matricesCurrentVC(i / 2, j / 2, mm);
          // point on corner lines at phi direction
          matricesCurrentV(i, j, mPlus) = 0.5 * (matricesCurrentVC(i / 2, j / 2, mm) + matricesCurrentVC(i / 2, j / 2, mmPlus));
        }
      }

      for (int j = 1; j < tnZColumn - 1; j += 2) {
        for (int i = 2; i < tnRRow - 1; i += 2) {
          matricesCurrentV(i, j, m) = 0.5 * (matricesCurrentVC(i / 2, j / 2, mm) + matricesCurrentVC(i / 2, j / 2 + 1, mm));
          // point on corner lines at phi direction
          matricesCurrentV(i, j, mPlus) = 0.25 * (matricesCurrentVC(i / 2, j / 2, mm) + matricesCurrentVC(i / 2, j / 2 + 1, mm) + matricesCurrentVC(i / 2, j / 2, mmPlus) +
                                                  matricesCurrentVC(i / 2, j / 2 + 1, mmPlus));
        }
      }

      for (int j = 2; j < tnZColumn - 1; j += 2) {
        for (int i = 1; i < tnRRow - 1; i += 2) {
          matricesCurrentV(i, j, m) = 0.5 * (matricesCurrentVC(i / 2, j / 2, mm) + matricesCurrentVC(i / 2 + 1, j / 2, mm));
          // point on line at phi direction
          matricesCurrentV(i, j, mPlus) = 0.25 * ((matricesCurrentVC(i / 2, j / 2, mm) + matricesCurrentVC(i / 2, j / 2, mmPlus)) +
                                                  (matricesCurrentVC(i / 2 + 1, j / 2, mmPlus) + matricesCurrentVC(i / 2 + 1, j / 2, mm)));
        }
      }

      for (int j = 1; j < tnZColumn - 1; j += 2) {
        for (int i = 1; i < tnRRow - 1; i += 2) {
          matricesCurrentV(i, j, m) = 0.25 * ((matricesCurrentVC(i / 2, j / 2, mm) + matricesCurrentVC(i / 2, j / 2 + 1, mm)) +
                                              (matricesCurrentVC(i / 2 + 1, j / 2, mm) + matricesCurrentVC(i / 2 + 1, j / 2 + 1, mm)));
          // point at the center at phi direction
          matricesCurrentV(i, j, mPlus) = 0.125 * ((matricesCurrentVC(i / 2, j / 2, mm) + matricesCurrentVC(i / 2, j / 2 + 1, mm) + matricesCurrentVC(i / 2, j / 2, mmPlus) +
                                                    matricesCurrentVC(i / 2, j / 2 + 1, mmPlus)) +
                                                   (matricesCurrentVC(i / 2 + 1, j / 2, mm) + matricesCurrentVC(i / 2 + 1, j / 2 + 1, mm) +
                                                    matricesCurrentVC(i / 2 + 1, j / 2, mmPlus) + matricesCurrentVC(i / 2 + 1, j / 2 + 1, mmPlus)));
        }
      }
    }

  } else {
    for (int m = 0; m < newPhiSlice; m++) {
      Interp2D(matricesCurrentV, matricesCurrentVC, tnRRow, tnZColumn, m);
    }
  }
}

/// Interpolation/Prolongation in 2D
///
/// Interpolation is a map from coarse grid (h) to fine grid (2h)
///
/// In case of 2D
/// Full weighting:
/// \f[ (R u)_{i,j,k} = \frac{1}{2} u_{2i,2j,2k} + \frac{1}{4} S_{1} + \frac{1}{8} S_{2} + \frac{1}{16} S_{3}\f]
///
///
/// Restriction in all direction r-phi-z
/// \param matricesCurrentV TMatrixD** finer grid h
/// \param curArrayCV TMatrixD ** coarse grid 2h
/// \param tnRRow int number of grid in Nr (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
/// \param tnZColumn int number of grid in Nz (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1
///
template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nr, Nz, Nphi>::Interp2D(Matrix3D& matricesCurrentV, Matrix3D& matricesCurrentVC, const int tnRRow,
                                                       const int tnZColumn, const int iphi)
{
  for (int j = 2; j < tnZColumn - 1; j += 2) {
    for (int i = 2; i < tnRRow - 1; i += 2) {
      matricesCurrentV(i, j, iphi) = matricesCurrentVC(i / 2, j / 2, iphi);
    }
  }

  for (int j = 1; j < tnZColumn - 1; j += 2) {
    for (int i = 2; i < tnRRow - 1; i += 2) {
      matricesCurrentV(i, j, iphi) = 0.5 * (matricesCurrentVC(i / 2, j / 2, iphi) + matricesCurrentVC(i / 2, j / 2 + 1, iphi));
    }
  }

  for (int j = 2; j < tnZColumn - 1; j += 2) {
    for (int i = 1; i < tnRRow - 1; i += 2) {
      matricesCurrentV(i, j, iphi) = 0.5 * (matricesCurrentVC(i / 2, j / 2, iphi) + matricesCurrentVC(i / 2 + 1, j / 2, iphi));
    }
  }

  // only if full
  if (fMgParameters.gtType == kFull) {
    for (int j = 1; j < tnZColumn - 1; j += 2) {
      for (int i = 1; i < tnRRow - 1; i += 2) {
        matricesCurrentV(i, j, iphi) = 0.25 * (matricesCurrentVC(i / 2, j / 2, iphi) + matricesCurrentVC(i / 2, j / 2 + 1, iphi) + matricesCurrentVC(i / 2 + 1, j / 2, iphi) + matricesCurrentVC(i / 2 + 1, j / 2 + 1, iphi));
      }
    }
  }
}

/// Prolongation with Addition for 3D
///
/// Interpolation with addition from coarse level (2h) -->  fine level (h)
///
/// Interpolation in all direction r-phi-z
/// Interpolation in phi only if oldPhi == 2*newPhi
/// \param matricesCurrentV TMatrixD& fine grid h
/// \param matricesCurrentVC TMatrixD& coarse grid 2h
/// \param tnRRow int number of grid in Nr (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
/// \param tnZColumn int number of grid in Nz (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1a
/// \param newPhiSlice int number of Nphi (in phi-direction) for coarser grid
/// \param oldPhiSlice int number of Nphi (in phi-direction) for finer grid
///
template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nr, Nz, Nphi>::AddInterp3D(Matrix3D& matricesCurrentV, Matrix3D& matricesCurrentVC, const int tnRRow,
                                                          const int tnZColumn,
                                                          const int newPhiSlice, const int oldPhiSlice)
{
  // Do restrict 2 D for each slice

  //const Float_t  h   =  (O2TPCPoissonSolver<DataT,Nr,Nz,Nphi>::fgkOFCRadius-O2TPCPoissonSolver<DataT,Nr,Nz,Nphi>::fgkIFCRadius) / ((tnRRow-1)/2); // h_{r}
  //Float_t radius,ratio;
  //std::vector<float> coefficient1((tnRRow-1) / 2 );  // coefficient1(Nr) for storing (1 + h_{r}/2r_{i}) from central differences in r direction
  //std::vector<float> coefficient2((tnRRow-1) / 2);  // coefficient2(Nr) for storing (1 + h_{r}/2r_{i}) from central differences in r direction

  if (newPhiSlice == 2 * oldPhiSlice) {
    int mPlus, mmPlus;
    int mm = 0;

    for (int m = 0; m < newPhiSlice; m += 2) {

      // assuming no symmetry
      mm = m / 2;
      mmPlus = mm + 1;
      mPlus = m + 1;

      // round
      if (mmPlus > (oldPhiSlice)-1) {
        mmPlus = mm + 1 - (oldPhiSlice);
      }
      if (mPlus > (newPhiSlice)-1) {
        mPlus = m + 1 - (newPhiSlice);
      }

      for (int j = 2; j < tnZColumn - 1; j += 2) {
        for (int i = 2; i < tnRRow - 1; i += 2) {
          matricesCurrentV(i, j, m) += matricesCurrentVC(i / 2, j / 2, mm);
          // point on corner lines at phi direction
          matricesCurrentV(i, j, mPlus) += 0.5 * (matricesCurrentVC(i / 2, j / 2, mm) + matricesCurrentVC(i / 2, j / 2, mmPlus));
        }
      }

      for (int j = 1; j < tnZColumn - 1; j += 2) {
        for (int i = 2; i < tnRRow - 1; i += 2) {
          matricesCurrentV(i, j, m) += 0.5 * (matricesCurrentVC(i / 2, j / 2, mm) + matricesCurrentVC(i / 2, j / 2 + 1, mm));
          // point on corner lines at phi direction
          matricesCurrentV(i, j, mPlus) += 0.25 * (matricesCurrentVC(i / 2, j / 2, mm) + matricesCurrentVC(i / 2, j / 2 + 1, mm) + matricesCurrentVC(i / 2, j / 2, mmPlus) +
                                                   matricesCurrentVC(i / 2, j / 2 + 1, mmPlus));
        }
      }

      for (int j = 2; j < tnZColumn - 1; j += 2) {
        for (int i = 1; i < tnRRow - 1; i += 2) {
          matricesCurrentV(i, j, m) += 0.5 * (matricesCurrentVC(i / 2, j / 2, mm) + matricesCurrentVC(i / 2 + 1, j / 2, mm));

          // point on line at phi direction
          matricesCurrentV(i, j, mPlus) += 0.25 * ((matricesCurrentVC(i / 2, j / 2, mm) + matricesCurrentVC(i / 2, j / 2, mmPlus)) +
                                                   (matricesCurrentVC(i / 2 + 1, j / 2, mmPlus) + matricesCurrentVC(i / 2 + 1, j / 2, mm)));
        }
      }

      for (int j = 1; j < tnZColumn - 1; j += 2) {
        for (int i = 1; i < tnRRow - 1; i += 2) {
          matricesCurrentV(i, j, m) += 0.25 * ((matricesCurrentVC(i / 2, j / 2, mm) + matricesCurrentVC(i / 2, j / 2 + 1, mm)) +
                                               (matricesCurrentVC(i / 2 + 1, j / 2, mm) + matricesCurrentVC(i / 2 + 1, j / 2 + 1, mm)));

          // point at the center at phi direction
          matricesCurrentV(i, j, mPlus) += 0.125 * ((matricesCurrentVC(i / 2, j / 2, mm) + matricesCurrentVC(i / 2, j / 2 + 1, mm) + matricesCurrentVC(i / 2, j / 2, mmPlus) +
                                                     matricesCurrentVC(i / 2, j / 2 + 1, mmPlus)) +
                                                    (matricesCurrentVC(i / 2 + 1, j / 2, mm) + matricesCurrentVC(i / 2 + 1, j / 2 + 1, mm) +
                                                     matricesCurrentVC(i / 2 + 1, j / 2, mmPlus) + matricesCurrentVC(i / 2 + 1, j / 2 + 1, mmPlus)));
        }
      }
    }

  } else {
    for (int m = 0; m < newPhiSlice; m++) {
      AddInterp2D(matricesCurrentV, matricesCurrentVC, tnRRow, tnZColumn, m);
    }
  }
}

/// Prolongation with Addition for 2D
///
/// Interpolation with addition from coarse level (2h) -->  fine level (h)
///
/// Interpolation in all direction r-phi-z
/// \param matricesCurrentV TMatrixD& fine grid h
/// \param matricesCurrentVC TMatrixD& coarse grid 2h
/// \param tnRRow int number of grid in Nr (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
/// \param tnZColumn int number of grid in Nz (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1a
///
template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nr, Nz, Nphi>::AddInterp2D(Matrix3D& matricesCurrentV, Matrix3D& matricesCurrentVC, const int tnRRow,
                                                          const int tnZColumn, const int iphi)
{
  for (int j = 2; j < tnZColumn - 1; j += 2) {
    for (int i = 2; i < tnRRow - 1; i += 2) {
      matricesCurrentV(i, j, iphi) = matricesCurrentV(i, j, iphi) + matricesCurrentVC(i / 2, j / 2, iphi);
    }
  }

  for (int j = 1; j < tnZColumn - 1; j += 2) {
    for (int i = 2; i < tnRRow - 1; i += 2) {
      matricesCurrentV(i, j, iphi) =
        matricesCurrentV(i, j, iphi) + 0.5 * (matricesCurrentVC(i / 2, j / 2, iphi) + matricesCurrentVC(i / 2, j / 2 + 1, iphi));
    }
  }

  for (int j = 2; j < tnZColumn - 1; j += 2) {
    for (int i = 1; i < tnRRow - 1; i += 2) {
      matricesCurrentV(i, j, iphi) =
        matricesCurrentV(i, j, iphi) + 0.5 * (matricesCurrentVC(i / 2, j / 2, iphi) + matricesCurrentVC(i / 2 + 1, j / 2, iphi));
    }
  }

  // only if full
  if (fMgParameters.gtType == kFull) {
    for (int j = 1; j < tnZColumn - 1; j += 2) {
      for (int i = 1; i < tnRRow - 1; i += 2) {
        matricesCurrentV(i, j, iphi) =
          matricesCurrentV(i, j, iphi) + 0.25 * (matricesCurrentVC(i / 2, j / 2, iphi) + matricesCurrentVC(i / 2, j / 2 + 1, iphi) + matricesCurrentVC(i / 2 + 1, j / 2, iphi) + matricesCurrentVC(i / 2 + 1, j / 2 + 1, iphi));
      }
    }
  }
}

/// Relax3D
///
///    Relaxation operation for multiGrid
///		 relaxation used 7 stencil in cylindrical coordinate
///
/// Using the following equations
/// \f$ U_{i,j,k} = (1 + \frac{1}{r_{i}h_{r}}) U_{i+1,j,k}  + (1 - \frac{1}{r_{i}h_{r}}) U_{i+1,j,k}  \f$
///
/// \param matricesCurrentV TMatrixD** potential in 3D (matrices of matrix)
/// \param matricesCurrentCharge TMatrixD** charge in 3D
/// \param Nr const int number of Nr in the r direction of TPC
/// \param Nz const int number of Nz in z direction of TPC
/// \param Nphi const int number of Nphi in phi direction of TPC
/// \param symmetry const int is the cylinder has symmetry
/// \param h2 const Float_t \f$  h_{r}^{2} \f$
/// \param tempRatioZ const Float_t ration between grid size in z-direction and r-direction
/// \param coefficient1 std::vector<float> coefficient for \f$  V_{x+1,y,z} \f$
/// \param coefficient2 std::vector<float> coefficient for \f$  V_{x-1,y,z} \f$
/// \param coefficient3 std::vector<float> coefficient for z
/// \param coefficient4 std::vector<float> coefficient for f(r,\phi,z)
///
template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nr, Nz, Nphi>::Relax3D(Matrix3D& matricesCurrentV, Matrix3D& matricesCurrentCharge, const int tnRRow,
                                                      const int tnZColumn,
                                                      const int symmetry, const Float_t h2,
                                                      const Float_t tempRatioZ, std::vector<float>& coefficient1,
                                                      std::vector<float>& coefficient2,
                                                      std::vector<float>& coefficient3, std::vector<float>& coefficient4)
{

  int mPlus, mMinus, signPlus, signMinus;

  // Gauss-Seidel (Read Black}
  if (fMgParameters.relaxType == kGaussSeidel) {
    // for each slice
    for (int iPass = 1; iPass <= 2; ++iPass) {
      const int msw = (iPass % 2) ? 1 : 2;
      for (int m = 0; m < Nphi; ++m) {
        const int jsw = ((msw + m) % 2) ? 1 : 2;

        mPlus = m + 1;
        signPlus = 1;
        mMinus = m - 1;
        signMinus = 1;
        // Reflection symmetry in phi (e.g. symmetry at sector boundaries, or half sectors, etc.)
        if (symmetry == 1) {
          if (mPlus > Nphi - 1) {
            mPlus = Nphi - 2;
          }
          if (mMinus < 0) {
            mMinus = 1;
          }
        }
        // Anti-symmetry in phi
        else if (symmetry == -1) {
          if (mPlus > Nphi - 1) {
            mPlus = Nphi - 2;
            signPlus = -1;
          }
          if (mMinus < 0) {
            mMinus = 1;
            signMinus = -1;
          }
        } else { // No Symmetries in phi, no boundaries, the calculation is continuous across all phi
          if (mPlus > Nphi - 1) {
            mPlus = m + 1 - Nphi;
          }
          if (mMinus < 0) {
            mMinus = m - 1 + Nphi;
          }
        }

        int isw = jsw;
        for (int j = 1; j < tnZColumn - 1; j++, isw = 3 - isw) {
          for (int i = isw; i < tnRRow - 1; i += 2) {
            (matricesCurrentV)(i, j, m) = (coefficient2[i] * (matricesCurrentV)(i - 1, j, m) + tempRatioZ * ((matricesCurrentV)(i, j - 1, m) + (matricesCurrentV)(i, j + 1, m)) + coefficient1[i] * (matricesCurrentV)(i + 1, j, m) + coefficient3[i] * (signPlus * (matricesCurrentV)(i, j, mPlus) + signMinus * (matricesCurrentV)(i, j, mMinus)) + (h2 * (matricesCurrentCharge)(i, j, m))) * coefficient4[i];
          } // end cols
        }   // end Nr
      }     // end phi
    }       // end sweep
  } else if (fMgParameters.relaxType == kJacobi) {
    // for each slice
    for (int m = 0; m < Nphi; m++) {

      mPlus = m + 1;
      signPlus = 1;
      mMinus = m - 1;
      signMinus = 1;

      // Reflection symmetry in phi (e.g. symmetry at sector boundaries, or half sectors, etc.)
      if (symmetry == 1) {
        if (mPlus > Nphi - 1) {
          mPlus = Nphi - 2;
        }
        if (mMinus < 0) {
          mMinus = 1;
        }
      }
      // Anti-symmetry in phi
      else if (symmetry == -1) {
        if (mPlus > Nphi - 1) {
          mPlus = Nphi - 2;
          signPlus = -1;
        }
        if (mMinus < 0) {
          mMinus = 1;
          signMinus = -1;
        }
      } else { // No Symmetries in phi, no boundaries, the calculation is continuous across all phi
        if (mPlus > Nphi - 1) {
          mPlus = m + 1 - Nphi;
        }
        if (mMinus < 0) {
          mMinus = m - 1 + Nphi;
        }
      }
      // Jacobian
      for (int j = 1; j < tnZColumn - 1; j++) {
        for (int i = 1; i < tnRRow - 1; i++) {
          (matricesCurrentV)(i, j, m) = (coefficient2[i] * (matricesCurrentV)(i - 1, j, m) + tempRatioZ * ((matricesCurrentV)(i, j - 1, m) + (matricesCurrentV)(i, j + 1, m)) + coefficient1[i] * (matricesCurrentV)(i + 1, j, m) + coefficient3[i] * (signPlus * (matricesCurrentV)(i, j, mPlus) + signMinus * (matricesCurrentV)(i, j, mMinus)) + (h2 * (matricesCurrentCharge)(i, j, m))) * coefficient4[i];
        } // end cols
      }   // end Nr

    } // end phi

  } else {
    // Case weighted Jacobi
    // TODO
  }
}

/// Restrict Boundary in 3D
///
/// Pass boundary information to coarse grid
///
/// \param matricesCurrentCharge TMatrixD** coarser grid 2h
/// \param residue TMatrixD ** fine grid h
/// \param tnRRow int number of grid in Nr (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
/// \param tnZColumn int number of grid in Nz (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1
/// \param newPhiSlice int number of Nphi (in phi-direction) for coarser grid
/// \param oldPhiSlice int number of Nphi (in phi-direction) for finer grid
///
template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nr, Nz, Nphi>::RestrictBoundary3D(Matrix3D& matricesCurrentCharge, Matrix3D& residue, const int tnRRow,
                                                                 const int tnZColumn, const int newPhiSlice, const int oldPhiSlice)
{
  // in case of full 3d and the Nphi is also coarsening

  if (2 * newPhiSlice == oldPhiSlice) {

    for (int m = 0, mm = 0; m < newPhiSlice; m++, mm += 2) {

      // TMatrixD& arrayResidue = *residue[mm];
      // TMatrixD& arrayCharge = *matricesCurrentCharge[m];
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
      RestrictBoundary2D(matricesCurrentCharge, residue, tnRRow, tnZColumn, m);
    }
  }
}

/// RestrictBoundary2D
///
///    Boundary transfer  restrict from fine -> coarse grid
///
/// \param matricesCurrentCharge TMatrixD& coarse grid (2h)
/// \param residue TMatrixD& fine grid  (h)
/// \param Nr const int number of Nr in the r direction of TPC
/// \param Nz const int number of Nz in z direction of TPC
///
template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nr, Nz, Nphi>::RestrictBoundary2D(Matrix3D& matricesCurrentCharge, Matrix3D& residue, const int tnRRow,
                                                                 const int tnZColumn, const int iphi)
{
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

/// Restriction in 3D
///
/// Restriction is a map from fine grid (h) to coarse grid (2h)
///
/// In case of 3D
/// Full weighting:
/// \f[ (R u)_{i,j,k} = \frac{1}{2} u_{2i,2j,2k} + \frac{1}{4} S_{1} + \frac{1}{8} S_{2} + \frac{1}{16} S_{3}\f]
///
///
/// Restriction in all direction r-phi-z
/// restriction in phi only if oldPhi == 2*newPhi
/// \param matricesCurrentCharge TMatrixD** coarser grid 2h
/// \param residue TMatrixD ** fine grid h
/// \param tnRRow int number of grid in Nr (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
/// \param tnZColumn int number of grid in Nz (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1
/// \param newPhiSlice int number of Nphi (in phi-direction) for coarser grid
/// \param oldPhiSlice int number of Nphi (in phi-direction) for finer grid
///
template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nr, Nz, Nphi>::Restrict3D(Matrix3D& matricesCurrentCharge, Matrix3D& residue, const int tnRRow,
                                                         const int tnZColumn,
                                                         const int newPhiSlice, const int oldPhiSlice)
{

  Double_t s1, s2, s3;

  if (2 * newPhiSlice == oldPhiSlice) {

    int mPlus, mMinus;
    int mm = 0;

    for (int m = 0; m < newPhiSlice; m++, mm += 2) {

      // assuming no symmetry
      mPlus = mm + 1;
      mMinus = mm - 1;

      if (mPlus > (oldPhiSlice)-1) {
        mPlus = mm + 1 - (oldPhiSlice);
      }
      if (mMinus < 0) {
        mMinus = mm - 1 + (oldPhiSlice);
      }

      // TMatrixD& arrayResidue = *residue[mm];
      // TMatrixD& arrayResidueP = *residue[mPlus];
      // TMatrixD& arrayResidueM = *residue[mMinus]; // slice
      // TMatrixD& arrayCharge = *matricesCurrentCharge[m];

      for (int i = 1, ii = 2; i < tnRRow - 1; i++, ii += 2) {
        for (int j = 1, jj = 2; j < tnZColumn - 1; j++, jj += 2) {

          // at the same plane
          s1 = residue(ii + 1, jj, mm) + residue(ii - 1, jj, mm) + residue(ii, jj + 1, mm) +
               residue(ii, jj - 1, mm) + residue(ii, jj, mPlus) + residue(ii, jj, mMinus);
          s2 = (residue(ii + 1, jj + 1, mm) + residue(ii + 1, jj - 1, mm) + residue(ii + 1, jj, mPlus) +
                residue(ii + 1, jj, mMinus)) +
               (residue(ii - 1, jj - 1, mm) + residue(ii - 1, jj + 1, mm) + residue(ii - 1, jj, mPlus) +
                residue(ii - 1, jj, mMinus)) +
               residue(ii, jj - 1, mPlus) + residue(ii, jj + 1, mMinus) + residue(ii, jj - 1, mMinus) +
               residue(ii, jj + 1, mPlus);

          s3 = (residue(ii + 1, jj + 1, mPlus) + residue(ii + 1, jj - 1, mPlus) + residue(ii + 1, jj + 1, mMinus) +
                residue(ii + 1, jj - 1, mMinus)) +
               (residue(ii - 1, jj - 1, mMinus) + residue(ii - 1, jj + 1, mMinus) + residue(ii - 1, jj - 1, mPlus) +
                residue(ii - 1, jj + 1, mPlus));

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
      Restrict2D(matricesCurrentCharge, residue, tnRRow, tnZColumn, m);
    }
  }
}

/// Restrict2D
///
///    Grid transfer operator, restrict from fine -> coarse grid
///		 provide full-half weighting
///
///		 \[ \frac{1}{16}\left( \begin{array}{ccc}
///      1 & 2 & 1 \\
///      2 & 4 & 2 \\
///      1 & 2 & 1 \end{array} \right) \]
///
/// \param matricesCurrentCharge TMatrixD& coarse grid (2h)
/// \param residue TMatrixD& fine grid  (h)
/// \param Nr const int number of Nr in the r direction of TPC
/// \param Nz const int number of Nz in z direction of TPC
///
template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCPoissonSolver<DataT, Nr, Nz, Nphi>::Restrict2D(Matrix3D& matricesCurrentCharge, Matrix3D& residue, const int tnRRow,
                                                         const int tnZColumn, const int iphi)
{

  for (int i = 1, ii = 2; i < tnRRow - 1; i++, ii += 2) {
    for (int j = 1, jj = 2; j < tnZColumn - 1; j++, jj += 2) {
      if (fMgParameters.gtType == kHalf) {
        // half
        matricesCurrentCharge(i, j, iphi) = 0.5 * residue(ii, jj, iphi) +
                                            0.125 *
                                              (residue(ii + 1, jj, iphi) + residue(ii - 1, jj, iphi) + residue(ii, jj + 1, iphi) +
                                               residue(ii, jj - 1, iphi));

      } else
        // full
        if (fMgParameters.gtType == kFull) {
        matricesCurrentCharge(i, j, iphi) = 0.25 * residue(ii, jj, iphi) +
                                            0.125 *
                                              (residue(ii + 1, jj, iphi) + residue(ii - 1, jj, iphi) + residue(ii, jj + 1, iphi) +
                                               residue(ii, jj - 1, iphi)) +
                                            0.0625 *
                                              (residue(ii + 1, jj + 1, iphi) + residue(ii - 1, jj + 1, iphi) + residue(ii + 1, jj - 1, iphi) +
                                               residue(ii - 1, jj - 1, iphi));
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

/// Helper function to check if the integer is equal to a power of two
/// \param i int the number
/// \return 1 if it is a power of two, else 0
template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
int O2TPCPoissonSolver<DataT, Nr, Nz, Nphi>::IsPowerOfTwo(int i) const
{
  int j = 0;
  while (i > 0) {
    j += (i & 1);
    i = (i >> 1);
  }
  if (j == 1) {
    return (1); // True
  }
  return (0); // False
}

///
/// Relative error calculation: comparison with exact solution
///
/// \param matricesCurrentV TMatrixD** current potential (numerical solution)
/// \param tempArrayV TMatrixD** temporary matrix for calculating error
/// \param Nr int number of grid in Nr (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
/// \param Nz int number of grid in Nz (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1
/// \param Nphi const int phi slices
///
template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
Double_t
  O2TPCPoissonSolver<DataT, Nr, Nz, Nphi>::GetConvergenceError(Matrix3D& matricesCurrentV, Matrix3D& prevArrayV)
{
  Double_t error = 0.0;

  std::transform(prevArrayV.mData.begin(), prevArrayV.mData.end(), matricesCurrentV.mData.begin(), prevArrayV.mData.begin(), std::minus<DataT>());

  for (int m = 0; m < Nphi; m++) {

    // square each entry in the vector and sum them up
    const auto phiStep = prevArrayV.nX * prevArrayV.nY;
    const auto start = prevArrayV.mData.begin() + m * phiStep;
    const auto end = start + phiStep;

    const DataT inner_product = std::inner_product(start, end, start, 0.);

    //get maximum error
    if (inner_product > error) {
      error = inner_product;
    }
  }
  return error;
}

///
/// Relative error calculation: comparison with exact solution
///
/// \param matricesCurrentV TMatrixD** current potential (numerical solution)
/// \param tempArrayV TMatrixD** temporary matrix for calculating error
/// \param Nr int number of grid in Nr (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
/// \param Nz int number of grid in Nz (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1
/// \param Nphi const int phi slices
///
// template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
// Double_t O2TPCPoissonSolver<DataT, Nr, Nz, Nphi>::GetExactError(DataContainer& matricesCurrentV, Matrix3D& tempArrayV, const int Nphi)
// {
//   Double_t error = 0.0;
//
//   if (fExactPresent == false) {
//
//     std::transform(fExactSolution.mData.begin(), fExactSolution.mData.end(), matricesCurrentV.begin(), tempArrayV.mData.begin(), std::minus<DataT>());
//     const DataT tmpExact = 1./fMaxExact;
//     std::transform(tempArrayV.mData.begin(), tempArrayV.mData.end(), tempArrayV.mData.begin(), [&tmpExact](auto& c){return c*tmpExact;});
//
//     for (int m = 0; m < Nphi; m++) {
//       // (*tempArrayV[m]) = (*fExactSolution[m]) - (*matricesCurrentV[m]);
//       // (*tempArrayV[m]) *= 1.0 / GetMaxExact();
//
//
//       // const DataT inner_product = std::inner_product( tempArrayV.mData.begin(), tempArrayV.mData.end(), tempArrayV.mData.begin(), 0 );
//       // square each entry in the vector and sum them up
//       const auto phiStep = tempArrayV.nX*tempArrayV.nY;
//       const auto start = tempArrayV.mData.begin() + m*phiStep;
//       const auto end = start + phiStep;
//
//       const DataT inner_product = std::inner_product( start, end, start, 0. );
//       std::cout<<"inner_product: "<<inner_product<<std::endl;
//
//       if (inner_product > error) {
//         error = inner_product;
//       }
//
//       // if (tempArrayV[m]->E2Norm() > error) {
//         // error = tempArrayV[m]->E2Norm();
//       // }
//       //printf("%f\n",tempArrayV[m]->E2Norm();
//     }
//   }
//   return error;
// }

template class O2TPCPoissonSolver<float, 17, 17, 90>;
template class O2TPCPoissonSolver<float, 65, 65, 90>;
template class O2TPCPoissonSolver<float, 129, 129, 180>;
template class O2TPCPoissonSolver<double, 129, 129, 180>;
