#ifndef O2TPCSpaceCharge3DCalc_H
#define O2TPCSpaceCharge3DCalc_H

#include "TriCubic.h"
#include "AliRoot/AliTPCPoissonSolver.h"

/// \tparam DataT the type of data which is used during the calculations
/// \tparam Nr number of vertices in r direction
/// \tparam Nz number of vertices in z direction
/// \tparam Nphi number of vertices in phi direction
template <typename DataT = float, unsigned int Nr = 4, unsigned int Nz = 4, unsigned int Nphi = 4>
class O2TPCSpaceCharge3DCalc
{

 public:
  O2TPCSpaceCharge3DCalc() = default;


  // stepp 0:
  // Info("AliTPCSpaceCharge3DCalc::InitSpaceCharge3DPoissonIntegralDz", "%s", Form("Step = 0: Fill Boundary and Charge Densities"));
  void fillBoundaryAndChargeDensities();

  // stepp 1:
  //Info("AliTPCSpaceCharge3DCalc::InitSpaceCharge3DPoissonIntegralDz", "%s", Form("Step 1: Poisson solver: %f\n", w.CpuTime()));
  void poissonSolver();

  // stepp 2:
  // Info("AliTPCSpaceCharge3DCalc::InitSpaceCharge3DPoissonIntegralDz", "%s", Form("Step 2: Electric Field Calculation: %f\n", w.CpuTime()));
  void calcEField();

  // stepp 3:
  // Info("AliTPCSpaceCharge3DCalc::InitSpaceCharge3DPoissonIntegralDz", "%s", Form("Step 3: Local distortion and correction cpu time: %f\n", w.CpuTime()));
  void calcLocalDistortionsCorrections();

  //step 4:
  // Info("AliTPCSpaceCharge3DCalc::InitSpaceCharge3DPoissonIntegralDz", "%s", Form("Step 4: Global correction/distortion cpu time: %f\n", w.CpuTime()));
  void calcGlobalDistortionsCorrections();


  constexpr DataT getGridSpacingR() const { return mGridSpacingR; }
  constexpr DataT getGridSpacingZ() const { return mGridSpacingZ; }
  constexpr DataT getGridSpacingPhi() const { return mGridSpacingPhi; }
  constexpr DataT getEzField() const { return (ASolv::fgkCathodeV - ASolv::fgkGG) / ASolv::fgkTPCZ0; }
  constexpr DataT getRMin() const { return ASolv::fgkIFCRadius; }
  RegularGrid3D<DataT, Nr, Nz, Nphi> getGrid3D() const { return mGrid3D; }

 private:
  using ASolv = AliTPCPoissonSolver<DataT>;

  static constexpr DataT mGridSpacingR = (ASolv::fgkOFCRadius - ASolv::fgkIFCRadius) / (Nr - 1); ///< grid spacing in r direction
  static constexpr DataT mGridSpacingZ = ASolv::fgkTPCZ0 / (Nz - 1); ///< grid spacing in z direction
  static constexpr DataT mGridSpacingPhi = 2 * M_PI / Nphi; ///< grid spacing in phi direction // TODO CHANGE TO O2
  static constexpr DataT mRMin = ASolv::fgkIFCRadius; ///< min radius
  static constexpr DataT mZMin = 0; ///< min z coordinate
  static constexpr DataT mPhiMin = 0; ///< min phi coordinate

  RegularGrid3D<DataT, Nr, Nz, Nphi> mGrid3D{mRMin, mZMin, mPhiMin, mGridSpacingR, mGridSpacingZ, mGridSpacingPhi}; ///< this grid contains the values for the local distortions/corrections, electric field etc.
};

#endif
