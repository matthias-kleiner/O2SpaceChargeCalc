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
#include "RegularGrid3D.h"

namespace o2
{
namespace tpc
{

/// The O2TPCPoissonSolver class represents methods to solve the poisson equation.
/// Original version with more methods can be found in AliTPCPoissonSolver.
/// only the default method is implemented in this class.

/// \tparam DataT the type of data which is used during the calculations
/// \tparam Nr number of vertices in r direction (2^N + 1)
/// \tparam Nz number of vertices in z direction (2^M + 1)
/// \tparam Nphi number of vertices in phi direction
template <typename DataT = float, size_t Nz = 129, size_t Nr = 129, size_t Nphi = 180>
class O2TPCPoissonSolver
{
 public:
  using RegularGrid = RegularGrid3D<DataT, Nz, Nr, Nphi>;
  using DataContainer = DataContainer3D<DataT, Nz, Nr, Nphi>;
  using Matrix3D = Matrix3D<DataT>;

  /// default constructor
  O2TPCPoissonSolver(const RegularGrid& gridProperties) : mGrid3D{gridProperties} {};

  /// Provides poisson solver in Cylindrical 3D (TPC geometry)
  ///
  /// Strategy based on parameter settings (mMgParameters)provided
  /// * Cascaded multi grid with S.O.R
  /// * Geometric MultiGrid
  ///		* Cycles: V, W, Full
  ///		* Relaxation: Jacobi, Weighted-Jacobi, Gauss-Seidel
  ///		* Grid transfer operators: Full, Half
  /// * Spectral Methods (TODO)
  ///
  /// \param matricesV potential in 3D
  /// \param matricesCharge charge density in 3D (side effect)
  /// \param symmetry symmetry or not
  ///
  /// \pre Charge density distribution in **matricesCharge** is known and boundary values for **matricesV** are set
  /// \post Numerical solution for potential distribution is calculated and stored in each rod at **matricesV**
  void poissonSolver3D(DataContainer& matricesV, const DataContainer& matricesCharge, const int symmetry);

  inline static DataT sConvergenceError{1e-3}; ///< Error tolerated

  DataT getSpacingZ() const { return mGrid3D.getSpacingX(); }
  DataT getSpacingR() const { return mGrid3D.getSpacingY(); }
  DataT getSpacingPhi() const { return mGrid3D.getSpacingZ(); }

 private:
  const RegularGrid& mGrid3D{}; ///< grid properties. member is set in O2TPCSpaceCharge3DCalc

  /// Relative error calculation: comparison with exact solution
  ///
  /// \param matricesCurrentV current potential (numerical solution)
  /// \param prevArrayV content from matricesCurrentV from previous iteration
  DataT getConvergenceError(const Matrix3D& matricesCurrentV, Matrix3D& prevArrayV) const;

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
  /// \param symmetry symmetry (TODO for symmetry = 1)
  //
  ///    SYMMETRY = 0 if no phi symmetries, and no phi boundary condition
  ///    = 1 if we have reflection symmetry at the boundaries (eg. sector symmetry or half sector symmetries).
  void poissonMultiGrid3D2D(DataContainer& matricesV, const DataContainer& matricesCharge, const int symmetry);

  void poissonMultiGrid3D(DataContainer& matricesV, const DataContainer& matricesCharge, const int symmetry);

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
  void restrict2D(Matrix3D& matricesCurrentCharge, const Matrix3D& residue, const int tnRRow, const int tnZColumn, const int iphi) const;

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
  /// \param matricesCurrentCharge coarser grid 2h
  /// \param residue fine grid h
  /// \param tnRRow number of grid in Nr (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
  /// \param tnZColumn number of grid in Nz (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1
  /// \param newPhiSlice number of Nphi (in phi-direction) for coarser grid
  /// \param oldPhiSlice number of Nphi (in phi-direction) for finer grid
  void restrict3D(Matrix3D& matricesCurrentCharge, const Matrix3D& residue, const int tnRRow, const int tnZColumn, const int newPhiSlice, const int oldPhiSlice) const;

  /// Restrict Boundary in 3D
  ///
  /// Pass boundary information to coarse grid
  ///
  /// \param matricesCurrentCharge coarser grid 2h
  /// \param residue fine grid h
  /// \param tnRRow number of grid in Nr (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
  /// \param tnZColumn number of grid in Nz (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1
  /// \param newPhiSlice number of Nphi (in phi-direction) for coarser grid
  /// \param oldPhiSlice number of Nphi (in phi-direction) for finer grid
  void restrictBoundary3D(Matrix3D& matricesCurrentCharge, const Matrix3D& residue, const int tnRRow, const int tnZColumn, const int newPhiSlice, const int oldPhiSlice) const;

  ///    Relaxation operation for multiGrid
  ///		 relaxation used 7 stencil in cylindrical coordinate
  ///
  /// Using the following equations
  /// \f$ U_{i,j,k} = (1 + \frac{1}{r_{i}h_{r}}) U_{i+1,j,k}  + (1 - \frac{1}{r_{i}h_{r}}) U_{i+1,j,k}  \f$
  ///
  /// \param matricesCurrentV potential in 3D (matrices of matrix)
  /// \param matricesCurrentCharge charge in 3D
  /// \param tnRRow number of grid in Nr (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
  /// \param tnZColumn number of grid in Nz (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1
  /// \param symmetry is the cylinder has symmetry
  /// \param h2 \f$  h_{r}^{2} \f$
  /// \param tempRatioZration between grid size in z-direction and r-direction
  /// \param coefficient1 coefficients for \f$  V_{x+1,y,z} \f$
  /// \param coefficient2 coefficients for \f$  V_{x-1,y,z} \f$
  /// \param coefficient3 coefficients for z
  /// \param coefficient4 coefficients for f(r,\phi,z)
  void relax3D(Matrix3D& matricesCurrentV, const Matrix3D& matricesCurrentCharge, const int tnRRow, const int tnZColumn, const int iPhi, const int symmetry, const DataT h2, const DataT tempRatioZ,
               const std::array<DataT, Nr>& coefficient1, const std::array<DataT, Nr>& coefficient2, const std::array<DataT, Nr>& coefficient3, const std::array<DataT, Nr>& coefficient4) const;

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
  /// \param matricesCurrentV finer grid h
  /// \param matricesCurrentVC coarse grid 2h
  /// \param tnRRow number of grid in Nr (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
  /// \param tnZColumn number of grid in Nz (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1
  /// \param iphi phi bin
  void interp2D(Matrix3D& matricesCurrentV, const Matrix3D& matricesCurrentVC, const int tnRRow, const int tnZColumn, const int iphi) const;

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
  /// \param matricesCurrentV finer grid h
  /// \param matricesCurrentVC coarse grid 2h
  /// \param tnRRow number of grid in Nr (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
  /// \param tnZColumn number of grid in Nz (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1
  /// \param newPhiSlice number of Nphi (in phi-direction) for coarser grid
  /// \param oldPhiSlice number of Nphi (in phi-direction) for finer grid
  void interp3D(Matrix3D& matricesCurrentV, const Matrix3D& matricesCurrentVC, const int tnRRow, const int tnZColumn, const int newPhiSlice, const int oldPhiSlice) const;

  /// Prolongation with Addition for 3D
  ///
  /// Interpolation with addition from coarse level (2h) -->  fine level (h)
  ///
  /// Interpolation in all direction r-phi-z
  /// Interpolation in phi only if oldPhi == 2*newPhi
  /// \param matricesCurrentV fine grid h
  /// \param matricesCurrentVC coarse grid 2h
  /// \param tnRRow number of grid in Nr (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
  /// \param tnZColumn number of grid in Nz (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1a
  /// \param newPhiSlice number of Nphi (in phi-direction) for coarser grid
  /// \param oldPhiSlice number of Nphi (in phi-direction) for finer grid
  void addInterp3D(Matrix3D& matricesCurrentV, const Matrix3D& matricesCurrentVC, const int tnRRow, const int tnZColumn, const int newPhiSlice, const int oldPhiSlice) const;

  /// Prolongation with Addition for 2D
  ///
  /// Interpolation with addition from coarse level (2h) -->  fine level (h)
  ///
  /// Interpolation in all direction r-phi-z
  /// \param matricesCurrentV fine grid h
  /// \param matricesCurrentVC coarse grid 2h
  /// \param tnRRow number of grid in Nr (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
  /// \param tnZColumn number of grid in Nz (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1a
  /// \param iphi phi bin
  void addInterp2D(Matrix3D& matricesCurrentV, const Matrix3D& matricesCurrentVC, const int tnRRow, const int tnZColumn, const int iphi) const;

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
  /// \param gridFrom finest level of grid
  /// \param gridTo coarsest level of grid
  /// \param nPre number of smoothing before coarsening
  /// \param nPost number of smoothing after coarsening
  /// \param ratioZ ratio between square of grid r and grid z (OPTION,  recalculate)
  /// \param ratioPhi ratio between square of grid r and grid phi (OPTION,  recalculate)
  /// \param tvArrayV vector of V potential in different grids
  /// \param tvCharge vector of charge distribution in different grids
  /// \param tvResidue vector of residue calculation in different grids
  /// \param coefficient1 coefficient for relaxation (r direction)
  /// \param coefficient2 coefficient for relaxation (r direction)
  /// \param coefficient3 coefficient for relaxation (ratio r/z)
  /// \param coefficient4 coefficient for relaxation (ratio for grid_r)
  /// \param inverseCoefficient4 coefficient for relaxation (inverse coefficient4)
  void vCycle3D2D(const int symmetry, const int gridFrom, const int gridTo, const int nPre, const int nPost, const DataT ratioZ, const DataT ratioPhi, std::vector<Matrix3D>& tvArrayV,
                  std::vector<Matrix3D>& tvCharge, std::vector<Matrix3D>& tvResidue, std::array<DataT, Nr>& coefficient1, std::array<DataT, Nr>& coefficient2, std::array<DataT, Nr>& coefficient3,
                  std::array<DataT, Nr>& coefficient4, std::array<DataT, Nr>& inverseCoefficient4) const;

  void vCycle3D(const int symmetry, const int gridFrom, const int gridTo, const int nPre, const int nPost, const DataT ratioZ, std::vector<Matrix3D>& tvArrayV,
                  std::vector<Matrix3D>& tvCharge, std::vector<Matrix3D>& tvResidue, std::array<DataT, Nr>& coefficient1, std::array<DataT, Nr>& coefficient2, std::array<DataT, Nr>& coefficient3,
                  std::array<DataT, Nr>& coefficient4, std::array<DataT, Nr>& inverseCoefficient4) const;

  /// Residue3D
  ///
  ///    Compute residue from V(.) where V(.) is numerical potential and f(.).
  ///		 residue used 7 stencil in cylindrical coordinate
  ///
  /// Using the following equations
  /// \f$ U_{i,j,k} = (1 + \frac{1}{r_{i}h_{r}}) U_{i+1,j,k}  + (1 - \frac{1}{r_{i}h_{r}}) U_{i+1,j,k}  \f$
  ///
  /// \param residue residue in 3D (matrices of matrix)
  /// \param matricesCurrentV potential in 3D (matrices of matrix)
  /// \param matricesCurrentCharge charge in 3D
  /// \param Nphi number of Nphi in phi direction of TPC
  /// \param symmetry is the cylinder has symmetry
  /// \param ih2 \f$ 1/ h_{r}^{2} \f$
  /// \param tempRatioZ ration between grid size in z-direction and r-direction
  /// \param coefficient1 coefficient for \f$  V_{x+1,y,z} \f$
  /// \param coefficient2 coefficient for \f$  V_{x-1,y,z} \f$
  /// \param coefficient3 coefficient for z
  /// \param inverseCoefficient4 inverse coefficient for f(r,\phi,z)
  void residue3D(Matrix3D& residue, const Matrix3D& matricesCurrentV, const Matrix3D& matricesCurrentCharge, const int tnRRow, const int tnZColumn, const int iPhi, const int symmetry, const DataT ih2, const DataT tempRatioZ,
                 const std::array<DataT, Nr>& coefficient1, const std::array<DataT, Nr>& coefficient2, const std::array<DataT, Nr>& coefficient3, const std::array<DataT, Nr>& inverseCoefficient4) const;

  ///    Boundary transfer  restrict from fine -> coarse grid
  ///
  /// \param matricesCurrentCharge coarse grid (2h)
  /// \param residue fine grid  (h)
  /// \param tnRRow number of tnRRow in the r direction of TPC
  /// \param tnZColumn number of tnZColumn in z direction of TPC
  /// \param iphi phi bin
  ///
  void restrictBoundary2D(Matrix3D& matricesCurrentCharge, const Matrix3D& residue, const int tnRRow, const int tnZColumn, const int iphi) const;

  // calculate coefficients
  void calcCoefficients(unsigned int from, unsigned int to, const DataT h, const DataT tempRatioZ, const DataT tempRatioPhi, std::array<DataT, Nr>& coefficient1,
                        std::array<DataT, Nr>& coefficient2, std::array<DataT, Nr>& coefficient3, std::array<DataT, Nr>& coefficient4) const;

  /// Helper function to check if the integer is equal to a power of two
  /// \param i the number
  /// \return 1 if it is a power of two, else 0
  bool isPowerOfTwo(int i) const;
};

} // namespace tpc
} // namespace o2

#endif
