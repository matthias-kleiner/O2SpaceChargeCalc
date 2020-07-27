#ifndef O2TPCSpaceCharge3DCalc_H
#define O2TPCSpaceCharge3DCalc_H

#include "TriCubic.h"
#include "AliRoot/AliTPCPoissonSolver.h"

template <typename DataT = float, unsigned int Nr = 4, unsigned int Nz = 4, unsigned int Nphi = 4>
class O2TPCSpaceCharge3DCalc
{

 public:
  O2TPCSpaceCharge3DCalc() = default;

  constexpr DataT getGridSpacingR() const { return mGridSpacingR; }
  constexpr DataT getGridSpacingZ() const { return mGridSpacingZ; }
  constexpr DataT getGridSpacingPhi() const { return mGridSpacingPhi; } // TODO CHANGE TO O2
  constexpr DataT getEzField() const { return (AliTPCPoissonSolver::fgkCathodeV - AliTPCPoissonSolver::fgkGG) / AliTPCPoissonSolver::fgkTPCZ0; }
  constexpr DataT getRMin() const { return AliTPCPoissonSolver::fgkIFCRadius; }
  RegularGrid3D<DataT, Nr, Nz, Nphi> getGrid3D() const { return mGrid3D; }

 private:
  static constexpr DataT mGridSpacingR = 3;// (AliTPCPoissonSolver::fgkOFCRadius - AliTPCPoissonSolver::fgkIFCRadius) / (Nr - 1) ;
  static constexpr DataT mGridSpacingZ = 3;//AliTPCPoissonSolver::fgkTPCZ0 / (Nz - 1);
  static constexpr DataT mGridSpacingPhi = 2 * M_PI / Nphi;
  static constexpr DataT mRMin = 0;//AliTPCPoissonSolver::fgkIFCRadius;
  static constexpr DataT mZMin = 0;
  static constexpr DataT mPhiMin = 0;

  RegularGrid3D<DataT, Nr, Nz, Nphi> mGrid3D{mRMin, mZMin, mPhiMin, mGridSpacingR, mGridSpacingZ, mGridSpacingPhi}; ///< this grid contains the values for the local distortions/corrections, electric field etc.
};

#endif
