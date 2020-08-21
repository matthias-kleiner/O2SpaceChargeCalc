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

#include <iostream>
#include <iomanip>

template <typename DataT = float, size_t Nr = 129, size_t Nz = 129, size_t Nphi = 180>
class O2TPCPoissonSolver : public TNamed
{
 public:
  using DataContainer = DataContainer3D<DataT, Nr, Nz, Nphi>;

  struct Matrix3D {

    Matrix3D(const unsigned int nXTmp, const unsigned int nYTmp, const unsigned int nZTmp) : nX{nYTmp}, nY{nXTmp}, nZ{nZTmp}
    {
      mData.resize(nXTmp * nYTmp * nZTmp);
    };

    Matrix3D(){};

    // ix and iy must be swapped since in the original version TMatrixD where used! TMatrixD(row,column)
    DataT& operator()(const unsigned int iy, const unsigned int ix, const unsigned int iz)
    {
      return mData[ix + nX * (iy + nY * iz)];
    }

    const DataT& operator()(const unsigned int iy, const unsigned int ix, const unsigned int iz) const
    {
      return mData[ix + nX * (iy + nY * iz)];
    }

    void resize(const unsigned int nXTmp, const unsigned int nYTmp, const unsigned int nZTmp)
    {
      nX = nXTmp;
      nY = nYTmp;
      nZ = nZTmp;
      mData.resize(nXTmp * nYTmp * nZTmp);
    }

    void print(const int nZStart = 0, int nZMax = -1)
    {
      nZMax = nZMax == -1 ? nZ - 1 : nZMax;
      std::ostream& out = std::cout;
      out.precision(3);
      auto&& w = std::setw(9);
      out << std::endl;

      for (unsigned int iz = nZStart; iz <= nZMax; ++iz) {
        out << "z layer: " << iz << std::endl;
        // print top x row
        out << "⎡" << w << mData[0 + nX * (0 + nY * iz)];
        for (unsigned int ix = 1; ix < nX; ++ix) {
          out << ", " << w << mData[ix + nX * (0 + nY * iz)];
        }
        out << " ⎤" << std::endl;

        for (unsigned int iy = 1; iy < nY - 1; ++iy) {
          out << "⎢" << w << mData[0 + nX * (iy + nY * iz)];
          for (unsigned int ix = 1; ix < nX; ++ix) {
            out << ", " << w << mData[ix + nX * (iy + nY * iz)];
          }
          out << " ⎥" << std::endl;
        }

        out << "⎣" << w << mData[0 + nX * ((nY - 1) + nY * iz)];
        for (unsigned int ix = 1; ix < nX; ++ix) {
          out << ", " << w << mData[ix + nX * ((nY - 1) + nY * iz)];
        }
        out << " ⎦" << std::endl;
        out << std::endl;
        out << std::endl;
      }
    }

    unsigned int nX{};
    unsigned int nY{};
    unsigned int nZ{};
    std::vector<DataT> mData{};
  };

  ///< Enumeration of Cycles Type
  enum CycleType {
    kVCycle = 0, ///< V Cycle
    kWCycle = 1, ///< W Cycle (TODO)
    kFCycle = 2  ///< Full Cycle
  };

  ///< Fine -> Coarse Grid transfer operator types
  enum GridTransferType {
    kHalf = 0, ///< Half weighting
    kFull = 1, ///< Full weighting
  };

  ///< Smoothing (Relax) operator types
  enum RelaxType {
    kJacobi = 0,         ///< Jacobi (5 Stencil 2D, 7 Stencil 3D_
    kWeightedJacobi = 1, ///< (TODO)
    kGaussSeidel = 2     ///< Gauss Seidel 2D (2 Color, 5 Stencil), 3D (7 Stencil)
  };

  ///< Parameters choice for MultiGrid    algorithm
  struct MGParameters {
    bool isFull3D = false;              ///<  TRUE: full coarsening, FALSE: semi coarsening
    CycleType cycleType = kFCycle;      ///< cycleType follow  CycleType
    GridTransferType gtType = kFull;    ///< gtType grid transfer type follow GridTransferType
    RelaxType relaxType = kGaussSeidel; ///< relaxType follow RelaxType
    int gamma;                          ///< number of iteration at coarsest level
    int nPre = 2;                       ///< number of iteration for pre smoothing
    int nPost = 2;                      ///< number of iteration for post smoothing
    int nMGCycle = 200;                 ///< number of multi grid cycle (V type)
    int maxLoop = 6;                    ///< the number of tree-deep of multi grid
  };

  /// constructor
  ///
  O2TPCPoissonSolver() : fErrorConvergenceNormInf(fMgParameters.nMGCycle) {}

  void PoissonSolver3D(DataContainer& matricesV, const DataContainer& matricesChargeDensities, const int symmetry);

  void PoissonMultiGrid3D2D(DataContainer& matricesV, const DataContainer& matricesChargeDensities, const int symmetry);
  void Restrict2D(Matrix3D& matrixCharge, const Matrix3D& residue, const int tnRRow, const int tnZColumn, const int iphi) const;
  void Restrict3D(Matrix3D& matricesCharge, const Matrix3D& residue, const int tnRRow, const int tnZColumn, const int newPhiSlice, const int oldPhiSlice) const;
  void RestrictBoundary3D(Matrix3D& matricesCharge, const Matrix3D& residue, const int tnRRow, const int tnZColumn, const int newPhiSlice, const int oldPhiSlice) const;
  void Relax3D(Matrix3D& currentMatricesV, const Matrix3D& matricesCharge, const int tnRRow, const int tnZColumn, const int symmetry, const DataT h2, const DataT tempRatioZ,
               const std::vector<DataT>& vectorCoefficient1, const std::vector<DataT>& vectorCoefficient2, const std::vector<DataT>& vectorCoefficient3, const std::vector<DataT>& vectorCoefficient4) const;

  void Interp2D(Matrix3D& matrixV, const Matrix3D& matrixVC, const int tnRRow, const int tnZColumn, const int iphi) const;
  void Interp3D(Matrix3D& currentMatricesV, const Matrix3D& currentMatricesVC, const int tnRRow, const int tnZColumn, const int newPhiSlice, const int oldPhiSlice) const;
  void AddInterp3D(Matrix3D& currentMatricesV, const Matrix3D& currentMatricesVC, const int tnRRow, const int tnZColumn, const int newPhiSlice, const int oldPhiSlice) const;
  void AddInterp2D(Matrix3D& matrixV, const Matrix3D& matrixVC, const int tnRRow, const int tnZColumn, const int iphi) const;

  void VCycle3D2D(const int symmetry, const int gridFrom, const int gridTo, const int nPre, const int nPost, const DataT ratioZ, const DataT ratioPhi, std::vector<Matrix3D>& tvArrayV,
                  std::vector<Matrix3D>& tvCharge, std::vector<Matrix3D>& tvResidue, std::vector<DataT>& vectorCoefficient1, std::vector<DataT>& vectorCoefficient2, std::vector<DataT>& vectorCoefficient3,
                  std::vector<DataT>& vectorCoefficient4, std::vector<DataT>& vectorInverseCoefficient4) const;

  void Residue3D(Matrix3D& residue, const Matrix3D& currentMatricesV, const Matrix3D& matricesCharge, const int tnRRow, const int tnZColumn, const int symmetry, const DataT ih2, const DataT tempRatio,
                 const std::vector<DataT>& vectorCoefficient1, const std::vector<DataT>& vectorCoefficient2, const std::vector<DataT>& vectorCoefficient3, const std::vector<DataT>& vectorInverseCoefficient4) const;

  void RestrictBoundary2D(Matrix3D& matrixCharge, const Matrix3D& residue, const int tnRRow, const int tnZColumn, const int iphi) const;

  bool IsPowerOfTwo(int i) const;
  DataT GetConvergenceError(const Matrix3D& currentMatricesV, Matrix3D& prevArrayV) const;

  static constexpr DataT fgkTPCZ0{249.7};                          ///< nominal gating grid position
  static constexpr DataT fgkIFCRadius{83.5};                       ///< Mean Radius of the Inner Field Cage ( 82.43 min,  83.70 max) (cm)
  static constexpr DataT fgkOFCRadius{254.5};                      ///< Mean Radius of the Outer Field Cage (252.55 min, 256.45 max) (cm)
  static constexpr DataT fgkZOffSet{0.2};                          ///< Offset from CE: calculate all distortions closer to CE as if at this point
  static constexpr DataT fgkCathodeV{-100000.0};                   ///< Cathode Voltage (volts)
  static constexpr DataT fgkGG{-70.0};                             ///< Gating Grid voltage (volts)
  static constexpr DataT fgkdvdE{0.0024};                          ///< [cm/V] drift velocity dependency on the E field (from Magboltz for NeCO2N2 at standard environment)
  static constexpr DataT fgkEM{-1.602176487e-19 / 9.10938215e-31}; ///< charge/mass in [C/kg]
  static constexpr DataT fgke0{8.854187817e-12};                   ///< vacuum permittivity [A·s/(V·m)]

  static constexpr DataT RMIN = fgkIFCRadius;                     ///< min radius
  static constexpr DataT ZMIN = 0;                                ///< min z coordinate
  static constexpr DataT PHIMIN = 0;                              ///< min phi coordinate
  static constexpr DataT RMAX = fgkOFCRadius;                     ///< max radius
  static constexpr DataT ZMAX = fgkTPCZ0;                         ///< max z coordinate
  static constexpr DataT PHIMAX = 2 * M_PI;                       ///< max phi coordinate // TODO CHANGE TO O2
  static constexpr DataT GRIDSPACINGR = (RMAX - RMIN) / (Nr - 1); ///< grid spacing in r direction
  static constexpr DataT GRIDSPACINGZ = (ZMAX - ZMIN) / (Nz - 1); ///< grid spacing in z direction
  static constexpr DataT GRIDSPACINGPHI = PHIMAX / Nphi;          ///< grid spacing in phi direction

  static constexpr DataT fgExactErr{1e-4};      ///< Error tolerated
  inline static DataT fgConvergenceError{1e-3}; ///< Error tolerated
  int fIterations;                              ///< number of maximum iteration
  MGParameters fMgParameters;                   ///< parameters multi grid
  std::vector<DataT> fErrorConvergenceNormInf{}; ///< for storing convergence error normInf

 private:
};

#endif
