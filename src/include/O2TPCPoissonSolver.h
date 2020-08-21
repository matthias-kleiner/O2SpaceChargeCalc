// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file O2TPCPoissonSolver.h
/// \brief This class provides implementation of Poisson Eq
/// solver by MultiGrid Method
/// Old version of this class can be found in the AliTPCPoissonSolver.h
///
///
/// \author  Matthias Kleiner <matthias.kleiner@cern.ch>
/// \date Aug 21, 2020

#ifndef O2TPCPOISSONSOLVER_H
#define O2TPCPOISSONSOLVER_H

#include "DataContainer3D.h"
#include "O2TPCPoissonSolverStructs.h"

/// \tparam DataT the type of data which is used during the calculations
/// \tparam Nr number of vertices in r direction
/// \tparam Nz number of vertices in z direction
/// \tparam Nphi number of vertices in phi direction
template <typename DataT = float, size_t Nr = 129, size_t Nz = 129, size_t Nphi = 180>
class O2TPCPoissonSolver
{
 public:
  using DataContainer = DataContainer3D<DataT, Nr, Nz, Nphi>;
  using Matrix3D = Matrix3D<DataT>;

  /// constructor
  ///
  O2TPCPoissonSolver() : mErrorConvergenceNormInf(mMgParameters.nMGCycle) {}

  /// Provides poisson solver in Cylindrical 3D (TPC geometry)
  ///
  /// Strategy based on parameter settings (fStrategy and mMgParameters)provided
  /// * Cascaded multi grid with S.O.R
  /// * Geometric MultiGrid
  ///		* Cycles: V, W, Full
  ///		* Relaxation: Jacobi, Weighted-Jacobi, Gauss-Seidel
  ///		* Grid transfer operators: Full, Half
  /// * Spectral Methods (TODO)
  ///
  /// \param matricesV potential in 3D matrix
  /// \param matricesCharge charge density in 3D matrix (side effect)
  /// \param symmetry symmetry or not
  ///
  /// \pre Charge density distribution in **matricesCharge** is known and boundary values for **matricesV** are set
  /// \post Numerical solution for potential distribution is calculated and stored in each rod at **matricesV**
  void PoissonSolver3D(DataContainer& matricesV, const DataContainer& matricesChargeDensities, const int symmetry);

  static constexpr DataT TPCZ0{249.7};                          ///< nominal gating grid position
  static constexpr DataT IFCRadius{83.5};                       ///< Mean Radius of the Inner Field Cage ( 82.43 min,  83.70 max) (cm)
  static constexpr DataT OFCRadius{254.5};                      ///< Mean Radius of the Outer Field Cage (252.55 min, 256.45 max) (cm)
  static constexpr DataT ZOffSet{0.2};                          ///< Offset from CE: calculate all distortions closer to CE as if at this point
  static constexpr DataT CATHODEV{-100000.0};                   ///< Cathode Voltage (volts)
  static constexpr DataT GG{-70.0};                             ///< Gating Grid voltage (volts)
  static constexpr DataT DVDE{0.0024};                          ///< [cm/V] drift velocity dependency on the E field (from Magboltz for NeCO2N2 at standard environment)
  static constexpr DataT EM{-1.602176487e-19 / 9.10938215e-31}; ///< charge/mass in [C/kg]
  static constexpr DataT E0{8.854187817e-12};                   ///< vacuum permittivity [A·s/(V·m)]

  static constexpr DataT RMIN = IFCRadius;                        ///< min radius
  static constexpr DataT ZMIN = 0;                                ///< min z coordinate
  static constexpr DataT PHIMIN = 0;                              ///< min phi coordinate
  static constexpr DataT RMAX = OFCRadius;                        ///< max radius
  static constexpr DataT ZMAX = TPCZ0;                            ///< max z coordinate
  static constexpr DataT PHIMAX = 2 * M_PI;                       ///< max phi coordinate // TODO CHANGE TO O2
  static constexpr DataT GRIDSPACINGR = (RMAX - RMIN) / (Nr - 1); ///< grid spacing in r direction
  static constexpr DataT GRIDSPACINGZ = (ZMAX - ZMIN) / (Nz - 1); ///< grid spacing in z direction
  static constexpr DataT GRIDSPACINGPHI = PHIMAX / Nphi;          ///< grid spacing in phi direction

  inline static DataT sConvergenceError{1e-3}; ///< Error tolerated

 private:
  // static constexpr DataT EXACTERR{1e-4};      ///< Error tolerated
  int mIterations{};                             ///< number of maximum iteration
  MGParameters mMgParameters{};                  ///< parameters multi grid
  std::vector<DataT> mErrorConvergenceNormInf{}; ///< for storing convergence error normInf

  ///
  /// Relative error calculation: comparison with exact solution
  ///
  /// \param matricesCurrentV TMatrixD** current potential (numerical solution)
  /// \param tempArrayV TMatrixD** temporary matrix for calculating error
  /// \param Nr int number of grid in Nr (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
  /// \param Nz int number of grid in Nz (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1
  /// \param Nphi const int phi slices
  ///
  DataT GetConvergenceError(const Matrix3D& currentMatricesV, Matrix3D& prevArrayV) const;

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
  /// \param matricesV potential in 3D matrix \f$ V(r,\phi,z) \f$
  /// \param matricesCharge charge density in 3D matrix (side effect) \f$ - f(r,\phi,z) \f$
  /// \param Nr number of bins in r direction of TPC
  /// \param Nz number of bins in z direction of TPC
  /// \param Nphi number of bins in phi direction of TPC
  /// \param symmetry symmetry (TODO for symmetry = 1)
  //
  ///    SYMMETRY = 0 if no phi symmetries, and no phi boundary condition
  ///    = 1 if we have reflection symmetry at the boundaries (eg. sector symmetry or half sector symmetries).
  ///
  void PoissonMultiGrid3D2D(DataContainer& matricesV, const DataContainer& matricesChargeDensities, const int symmetry);

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
  /// \param matricesCurrentCharge coarse grid (2h)
  /// \param residue fine grid  (h)
  /// \param tnRRow number of bins in r direction of TPC
  /// \param tnZColumn number of bins in z direction of TPC
  /// \param iphi phi bin
  ///
  void Restrict2D(Matrix3D& matrixCharge, const Matrix3D& residue, const int tnRRow, const int tnZColumn, const int iphi) const;

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
  void Restrict3D(Matrix3D& matricesCharge, const Matrix3D& residue, const int tnRRow, const int tnZColumn, const int newPhiSlice, const int oldPhiSlice) const;

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
  void RestrictBoundary3D(Matrix3D& matricesCharge, const Matrix3D& residue, const int tnRRow, const int tnZColumn, const int newPhiSlice, const int oldPhiSlice) const;

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
  void Relax3D(Matrix3D& currentMatricesV, const Matrix3D& matricesCharge, const int tnRRow, const int tnZColumn, const int symmetry, const DataT h2, const DataT tempRatioZ,
               const std::array<DataT, Nr>& vectorCoefficient1, const std::array<DataT, Nr>& vectorCoefficient2, const std::array<DataT, Nr>& vectorCoefficient3, const std::array<DataT, Nr>& vectorCoefficient4) const;

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
  void Interp2D(Matrix3D& matrixV, const Matrix3D& matrixVC, const int tnRRow, const int tnZColumn, const int iphi) const;

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
  void Interp3D(Matrix3D& currentMatricesV, const Matrix3D& currentMatricesVC, const int tnRRow, const int tnZColumn, const int newPhiSlice, const int oldPhiSlice) const;

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
  void AddInterp3D(Matrix3D& currentMatricesV, const Matrix3D& currentMatricesVC, const int tnRRow, const int tnZColumn, const int newPhiSlice, const int oldPhiSlice) const;

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
  void AddInterp2D(Matrix3D& matrixV, const Matrix3D& matrixVC, const int tnRRow, const int tnZColumn, const int iphi) const;

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
  void VCycle3D2D(const int symmetry, const int gridFrom, const int gridTo, const int nPre, const int nPost, const DataT ratioZ, const DataT ratioPhi, std::vector<Matrix3D>& tvArrayV,
                  std::vector<Matrix3D>& tvCharge, std::vector<Matrix3D>& tvResidue, std::array<DataT, Nr>& vectorCoefficient1, std::array<DataT, Nr>& vectorCoefficient2, std::array<DataT, Nr>& vectorCoefficient3,
                  std::array<DataT, Nr>& vectorCoefficient4, std::array<DataT, Nr>& vectorInverseCoefficient4) const;

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
  void Residue3D(Matrix3D& residue, const Matrix3D& currentMatricesV, const Matrix3D& matricesCharge, const int tnRRow, const int tnZColumn, const int symmetry, const DataT ih2, const DataT tempRatio,
                 const std::array<DataT, Nr>& vectorCoefficient1, const std::array<DataT, Nr>& vectorCoefficient2, const std::array<DataT, Nr>& vectorCoefficient3, const std::array<DataT, Nr>& vectorInverseCoefficient4) const;

  /// RestrictBoundary2D
  ///
  ///    Boundary transfer  restrict from fine -> coarse grid
  ///
  /// \param matricesCurrentCharge TMatrixD& coarse grid (2h)
  /// \param residue TMatrixD& fine grid  (h)
  /// \param tnRRow const int number of tnRRow in the r direction of TPC
  /// \param tnZColumn const int number of tnZColumn in z direction of TPC
  /// \param iphi phi bin
  ///
  void RestrictBoundary2D(Matrix3D& matrixCharge, const Matrix3D& residue, const int tnRRow, const int tnZColumn, const int iphi) const;

  /// Helper function to check if the integer is equal to a power of two
  /// \param i int the number
  /// \return 1 if it is a power of two, else 0
  bool IsPowerOfTwo(int i) const;
};

#endif
