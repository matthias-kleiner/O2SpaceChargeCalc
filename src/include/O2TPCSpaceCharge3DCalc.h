// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  O2TPCSpaceCharge3DCalc.h
/// \brief Definition of O2TPCSpaceCharge3DCalc class
///
/// \author  Matthias Kleiner <matthias.kleiner@cern.ch>

#ifndef O2TPCSpaceCharge3DCalc_H
#define O2TPCSpaceCharge3DCalc_H

#include "TriCubic.h"
#include "AliRoot/AliTPCPoissonSolver.h"
#include "SpaceChargeStructs.h"

// Root includes
#include "TF1.h" /// for numerical intergration only
#include "TTree.h"

#ifdef WITH_OPENMP
#include <omp.h>
#endif

// temporary for nearest neighbour search
#include "GeometricalTools/NearestNeighborQuery.h"

/// \tparam DataT the type of data which is used during the calculations
/// \tparam Nr number of vertices in r direction
/// \tparam Nz number of vertices in z direction
/// \tparam Nphi number of vertices in phi direction
template <typename DataT = float, size_t Nr = 17, size_t Nz = 17, size_t Nphi = 90>
class O2TPCSpaceCharge3DCalc
{
  using RegularGrid = RegularGrid3D<DataT, Nz, Nr, Nphi>;

 public:
  O2TPCSpaceCharge3DCalc() = default;

  // stepp 0:
  // TODO change this to accept histogram as input for density
  void fillBoundaryAndChargeDensities(const AnalyticalFields<DataT>& formulaStruct);

  // stepp 1:
  void poissonSolver(const int maxIteration = 300, const DataT stoppingConvergence = 1e-8);

  // stepp 2:
  void calcEField();

  // stepp 3:
  /// lcorrections=false -> distortions, lcorrections=true->corrections
  template <typename ElectricFields = AnalyticalFields<DataT>>
  void calcLocalDistortionsCorrections(const bool lcorrections, ElectricFields& formulaStruct);

  //step 4:
  // Info("AliTPCSpaceCharge3DCalc::InitSpaceCharge3DPoissonIntegralDz", "%s", Form("Step 4: Global correction/distortion cpu time: %f\n", w.CpuTime()));
  void calcGlobalDistortionsCorrections();

  template <typename ElectricFields = AnalyticalFields<DataT>>
  void calcGlobalDistortions(const ElectricFields& formulaStruct);

  template <typename ElectricFields = AnalyticalFields<DataT>>
  void calcGlobalCorrections(const ElectricFields& formulaStruct);

  // calculate global distortions using the global corrections RENAME NAME
  void calcGlobalDistWithGlobalCorrIterative(const DistCorrInterpolator<DataT, RegularGrid>& globCorr, const int maxIter = 200, const DataT convZ = 0.05, const DataT convR = 0.05, const DataT convPhi = 0.05);

  static constexpr DataT getGridSpacingR() { return mGridSpacingR; }
  static constexpr DataT getGridSpacingZ() { return mGridSpacingZ; }
  static constexpr DataT getGridSpacingPhi() { return mGridSpacingPhi; }
  static constexpr DataT getEzField() { return (ASolv::fgkCathodeV - ASolv::fgkGG) / ASolv::fgkTPCZ0; }

  static constexpr DataT getRMin() { return mRMin; }
  static constexpr DataT getZMin() { return mZMin; }
  static constexpr DataT getRMax() { return mGridSpacingR * (Nr - 1) + mRMin; };
  static constexpr DataT getZMax() { return mGridSpacingZ * (Nz - 1) + mZMin; }

  const RegularGrid& getGrid3D() const { return mGrid3D; }

  NumericalFields<DataT, RegularGrid> getNumericalFieldsInterpolator() const
  {
    if (!mIsEfieldSet) {
      std::cout << "============== E-Fields are not set! returning ==============" << std::endl;
    }
    NumericalFields<DataT, RegularGrid> numFields(mElectricFieldEr, mElectricFieldEz, mElectricFieldEphi);
    return numFields;
  }

  DistCorrInterpolator<DataT, RegularGrid> getLocalDistInterpolator() const
  {
    if (!mIsLocalDistSet) {
      std::cout << "============== local distortions not set! returning ==============" << std::endl;
    }
    DistCorrInterpolator<DataT, RegularGrid> numFields(mLocalDistdR, mLocalDistdZ, mLocalDistdRPhi);
    return numFields;
  }

  DistCorrInterpolator<DataT, RegularGrid> getLocalCorrInterpolator() const
  {
    if (!mIsLocalCorrSet) {
      std::cout << "============== local corrections not set! returning ==============" << std::endl;
    }
    DistCorrInterpolator<DataT, RegularGrid> numFields(mLocalCorrdR, mLocalCorrdZ, mLocalCorrdRPhi);
    return numFields;
  }

  DistCorrInterpolator<DataT, RegularGrid> getGlobalDistInterpolator() const
  {
    if (!mIsGlobalDistSet) {
      std::cout << "============== global distortions not set! returning ==============" << std::endl;
    }
    DistCorrInterpolator<DataT, RegularGrid> numFields(mGlobalDistdR, mGlobalDistdZ, mGlobalDistdRPhi);
    return numFields;
  }

  DistCorrInterpolator<DataT, RegularGrid> getGlobalCorrInterpolator() const
  {
    if (!mIsGlobalCorrSet) {
      std::cout << "============== global corrections not set! returning ==============" << std::endl;
    }
    DistCorrInterpolator<DataT, RegularGrid> numFields(mGlobalCorrdR, mGlobalCorrdZ, mGlobalCorrdRPhi);
    return numFields;
  }

  DataT getLocalDistR(size_t iz, size_t ir, size_t iphi) const
  {
    return mLocalDistdR(iz, ir, iphi);
  }

  DataT getLocalDistZ(size_t iz, size_t ir, size_t iphi) const
  {
    return mLocalDistdZ(iz, ir, iphi);
  }

  DataT getLocalDistRPhi(size_t iz, size_t ir, size_t iphi) const
  {
    return mLocalDistdRPhi(iz, ir, iphi);
  }

  DataT getLocalCorrR(size_t iz, size_t ir, size_t iphi) const
  {
    return mLocalCorrdR(iz, ir, iphi);
  }

  DataT getLocalCorrZ(size_t iz, size_t ir, size_t iphi) const
  {
    return mLocalCorrdZ(iz, ir, iphi);
  }

  DataT getLocalCorrRPhi(size_t iz, size_t ir, size_t iphi) const
  {
    return mLocalCorrdRPhi(iz, ir, iphi);
  }

  DataT getGlobalDistR(size_t iz, size_t ir, size_t iphi) const
  {
    return mGlobalDistdR(iz, ir, iphi);
  }

  DataT getGlobalDistZ(size_t iz, size_t ir, size_t iphi) const
  {
    return mGlobalDistdZ(iz, ir, iphi);
  }

  DataT getGlobalDistRPhi(size_t iz, size_t ir, size_t iphi) const
  {
    return mGlobalDistdRPhi(iz, ir, iphi);
  }

  DataT getGlobalCorrR(size_t iz, size_t ir, size_t iphi) const
  {
    return mGlobalCorrdR(iz, ir, iphi);
  }

  DataT getGlobalCorrZ(size_t iz, size_t ir, size_t iphi) const
  {
    return mGlobalCorrdZ(iz, ir, iphi);
  }

  DataT getGlobalCorrRPhi(size_t iz, size_t ir, size_t iphi) const
  {
    return mGlobalCorrdRPhi(iz, ir, iphi);
  }

  DataT getEr(size_t iz, size_t ir, size_t iphi) const
  {
    return mElectricFieldEr(iz, ir, iphi);
  }

  /// return the numerically calculated electric field. To get the analytical field use the struct eval function
  DataT getEz(size_t iz, size_t ir, size_t iphi) const
  {
    return mElectricFieldEz(iz, ir, iphi);
  }

  DataT getEphi(size_t iz, size_t ir, size_t iphi) const
  {
    return mElectricFieldEphi(iz, ir, iphi);
  }

  int getIntegrationSteps() const { return mStepWidth; }

  DataT getPhiVertex(size_t indexPhi) const
  {
    return mGrid3D.getZVertex(indexPhi);
  }

  DataT getRVertex(size_t indexR) const
  {
    return mGrid3D.getYVertex(indexR);
  }

  DataT getZVertex(size_t indexZ) const
  {
    return mGrid3D.getXVertex(indexZ);
  }

  const RegularGrid& getGridEr() const
  {
    return mElectricFieldEr;
  }

  const RegularGrid& getGridEz() const
  {
    return mElectricFieldEz;
  }

  const RegularGrid& getGridEphi() const
  {
    return mElectricFieldEphi;
  }

  void setOmegaTauT1T2(DataT omegaTau, DataT t1, DataT t2)
  {
    const DataT wt0 = t2 * omegaTau;
    fC0 = 1 / (1 + wt0 * wt0);
    const DataT wt1 = t1 * omegaTau;
    fC1 = wt1 / (1 + wt1 * wt1);
  };

  void setC0C1(DataT c0, DataT c1)
  {
    fC0 = c0;
    fC1 = c1;
  }

  void setIntegrationSteps(const int nSteps) { mStepWidth = nSteps; }

  /// numerical integration strategys
  enum IntegrationStrategy { Trapezoidal = 0,
                             Simpson = 1,
                             Root = 2,
                             SimpsonExperimental = 3 };

  // calculat the global distortions with interpolation of local distortions or by taking the electric field
  enum GlobalDistortionType {
    LocalDistortions = 0,
    ElectricField = 1
  };

  int dumpElectricFields(TFile& outf) const
  {
    if (!mIsEfieldSet) {
      std::cout << "============== E-Fields are not set! returning ==============" << std::endl;
      return 0;
    }
    const int er = mElectricFieldEr.storeValuesToFile(outf, "fieldEr");
    const int ez = mElectricFieldEz.storeValuesToFile(outf, "fieldEz");
    const int ephi = mElectricFieldEphi.storeValuesToFile(outf, "fieldEphi");
    return er + ez + ephi;
  }

  void setElectricFieldsFromFile(TFile& inpf)
  {
    mElectricFieldEr.initFromFile(inpf, "fieldEr");
    mElectricFieldEz.initFromFile(inpf, "fieldEz");
    mElectricFieldEphi.initFromFile(inpf, "fieldEphi");
    mIsEfieldSet = true;
  }

  int dumpPotential(TFile& outf) const
  {
    return mPotential.storeValuesToFile(outf, "potential");
  }

  void setPotentialFromFile(TFile& inpf)
  {
    mPotential.initFromFile(inpf, "potential");
  }

  int dumpGlobalDistortions(TFile& outf) const
  {
    if (!mIsGlobalDistSet) {
      std::cout << "============== global distortions are not set! returning ==============" << std::endl;
      return 0;
    }
    const int er = mGlobalDistdR.storeValuesToFile(outf, "distR");
    const int ez = mGlobalDistdZ.storeValuesToFile(outf, "distZ");
    const int ephi = mGlobalDistdRPhi.storeValuesToFile(outf, "distRPhi");
    return er + ez + ephi;
  }

  void setGlobalDistortionsFromFile(TFile& inpf)
  {
    mIsGlobalDistSet = true;
    mGlobalDistdR.initFromFile(inpf, "distR");
    mGlobalDistdZ.initFromFile(inpf, "distZ");
    mGlobalDistdRPhi.initFromFile(inpf, "distRPhi");
  }

  int dumpGlobalCorrections(TFile& outf) const
  {
    if (!mIsGlobalCorrSet) {
      std::cout << "============== global corrections are not set! returning ==============" << std::endl;
      return 0;
    }
    const int er = mGlobalCorrdR.storeValuesToFile(outf, "corrR");
    const int ez = mGlobalCorrdZ.storeValuesToFile(outf, "corrZ");
    const int ephi = mGlobalCorrdRPhi.storeValuesToFile(outf, "corrRPhi");
    return er + ez + ephi;
  }

  void setGlobalCorrectionsFromFile(TFile& inpf)
  {
    mIsGlobalCorrSet = true;
    mGlobalCorrdR.initFromFile(inpf, "corrR");
    mGlobalCorrdZ.initFromFile(inpf, "corrZ");
    mGlobalCorrdRPhi.initFromFile(inpf, "corrRPhi");
  }

  int dumpLocalCorrections(TFile& outf) const
  {
    if (!mIsLocalCorrSet) {
      std::cout << "============== local corrections are not set! returning ==============" << std::endl;
      return 0;
    }
    const int lCorrdR = mLocalCorrdR.storeValuesToFile(outf, "lcorrR");
    const int lCorrdZ = mLocalCorrdZ.storeValuesToFile(outf, "lcorrZ");
    const int lCorrdRPhi = mLocalCorrdRPhi.storeValuesToFile(outf, "lcorrRPhi");
    return lCorrdR + lCorrdZ + lCorrdRPhi;
  }

  void setLocalCorrectionsFromFile(TFile& inpf)
  {
    mIsLocalCorrSet = true;
    mLocalCorrdR.initFromFile(inpf, "lcorrR");
    mLocalCorrdZ.initFromFile(inpf, "lcorrZ");
    mLocalCorrdRPhi.initFromFile(inpf, "lcorrRPhi");
  }

  int dumpLocalDistortions(TFile& outf) const
  {
    if (!mIsLocalDistSet) {
      std::cout << "============== local distortions are not set! returning ==============" << std::endl;
      return 0;
    }
    const int lCorrdR = mLocalDistdR.storeValuesToFile(outf, "ldistR");
    const int lCorrdZ = mLocalDistdZ.storeValuesToFile(outf, "ldistZ");
    const int lCorrdRPhi = mLocalDistdRPhi.storeValuesToFile(outf, "ldistRPhi");
    return lCorrdR + lCorrdZ + lCorrdRPhi;
  }

  void setLocalDistortionsFromFile(TFile& inpf)
  {
    mIsLocalDistSet = true;
    mLocalDistdR.initFromFile(inpf, "ldistR");
    mLocalDistdZ.initFromFile(inpf, "ldistZ");
    mLocalDistdRPhi.initFromFile(inpf, "ldistRPhi");
  }

  static DataT regulatePhi(const DataT phi)
  {
    const DataT twoPi = 2 * M_PI;
    DataT phiTmp = phi;
    while (phiTmp < 0.0) {
      phiTmp += twoPi; // TODO USE O2 for twoPi
    }
    while (phiTmp > twoPi) {
      phiTmp -= twoPi;
    }
    return phiTmp;
  }

  DataT regulateZ(const DataT pos) const
  {
    return mGrid3D.clampToGrid(pos, 0);
  }

  DataT regulateR(const DataT pos) const
  {
    return mGrid3D.clampToGrid(pos, 1);
  }

 private:
  using ASolv = AliTPCPoissonSolver<DataT>;

  static constexpr DataT mRMin = ASolv::fgkIFCRadius; ///< min radius
  static constexpr DataT mZMin = ASolv::fgkZOffSet;   ///< min z coordinate
  static constexpr DataT mPhiMin = 0;
  static constexpr DataT mGridSpacingR = (ASolv::fgkOFCRadius - mRMin) / (Nr - 1); ///< grid spacing in r direction
  static constexpr DataT mGridSpacingZ = (ASolv::fgkTPCZ0 - mZMin) / (Nz - 1);     ///< grid spacing in z direction
  static constexpr DataT mGridSpacingPhi = 2 * M_PI / Nphi;                        ///< grid spacing in phi direction // TODO CHANGE TO O2                                                           ///< min phi coordinate
  int mNumericalIntegrationStrategy = SimpsonExperimental;                         ///< numerical integration strategy of integration of the E-Field: 0: trapezoidal, 1: Simpson, 2: Root (only for analytical formula case)
  unsigned int mNumericalIntegrationSteps = 1;                                     ///< number of intervalls during numerical integration are taken.
  int mStepWidth = 1;                                                              ///< during the calculation of the corrections/distortions it is assumed that the electron drifts on a line from deltaZ = z0 -> z1. The value sets the deltaZ width: 1: deltaZ=zBin/1, 5: deltaZ=zBin/5

  DataT fC0 = 0; ///< coefficient C0 (compare Jim Thomas's notes for definitions)
  DataT fC1 = 0; ///< coefficient C1 (compare Jim Thomas's notes for definitions)

  bool mIsEfieldSet = false;
  bool mIsLocalCorrSet = false;
  bool mIsLocalDistSet = false;
  bool mIsGlobalCorrSet = false;
  bool mIsGlobalDistSet = false;

  RegularGrid mGrid3D{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi}; ///< this grid contains the values for the local distortions/corrections, electric field etc.

  RegularGrid mLocalDistdR{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};
  RegularGrid mLocalDistdZ{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};
  RegularGrid mLocalDistdRPhi{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};

  RegularGrid mLocalCorrdR{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};
  RegularGrid mLocalCorrdZ{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};
  RegularGrid mLocalCorrdRPhi{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};

  RegularGrid mGlobalDistdR{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};
  RegularGrid mGlobalDistdZ{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};
  RegularGrid mGlobalDistdRPhi{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};

  RegularGrid mGlobalCorrdR{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};
  RegularGrid mGlobalCorrdZ{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};
  RegularGrid mGlobalCorrdRPhi{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};

  RegularGrid mDensity{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};
  RegularGrid mPotential{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};
  RegularGrid mElectricFieldEr{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};
  RegularGrid mElectricFieldEz{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};
  RegularGrid mElectricFieldEphi{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};

  DataT getInvSpacingR() const
  {
    return mGrid3D.getInvSpacingY();
  }

  DataT getInvSpacingZ() const
  {
    return mGrid3D.getInvSpacingX();
  }

  DataT getInvSpacingPhi() const
  {
    return mGrid3D.getInvSpacingZ();
  }

  /// calculate distortions or corrections analytical with electric fields
  template <typename ElectricFields = AnalyticalFields<DataT>>
  void calcDistortions(const DataT p1r, const DataT p1phi, const DataT p1z, const DataT p2z, DataT& ddR, DataT& ddRPhi, DataT& ddZ, ElectricFields& formulaStruct) const;

  template <typename ElectricFields = AnalyticalFields<DataT>>
  void integrateEFieldsRoot(const DataT p1r, const DataT p1phi, const DataT p1z, const DataT p2z, DataT& localIntErOverEz, DataT& localIntEPhiOverEz, DataT& localIntDeltaEz, ElectricFields& formulaStruct) const;

  template <typename ElectricFields = AnalyticalFields<DataT>>
  void integrateEFieldsTrapezoidal(const DataT p1r, const DataT p1phi, const DataT p1z, const DataT p2z, DataT& localIntErOverEz, DataT& localIntEPhiOverEz, DataT& localIntDeltaEz, ElectricFields& formulaStruct) const;

  template <typename ElectricFields = AnalyticalFields<DataT>>
  void integrateEFieldsSimpson(const DataT p1r, const DataT p1phi, const DataT p1z, const DataT p2z, DataT& localIntErOverEz, DataT& localIntEPhiOverEz, DataT& localIntDeltaEz, ElectricFields& formulaStruct) const;

  template <typename ElectricFields = AnalyticalFields<DataT>>
  void integrateEFieldsSimpsonExperimental(const DataT p1r, const DataT p2r, const DataT p1phi, const DataT p2phi, const DataT p1z, const DataT p2z, DataT& localIntErOverEz, DataT& localIntEPhiOverEz, DataT& localIntDeltaEz, ElectricFields& formulaStruct) const;

  void processGlobalDistCorr(const DataT radius, const DataT phi, const DataT z0Tmp, const DataT z1Tmp, DataT& ddR, DataT& ddPhi, DataT& ddZ, const AnalyticalFields<DataT>& formulaStruct) const
  {
    calcDistortions(radius, phi, z0Tmp, z1Tmp, ddR, ddPhi, ddZ, formulaStruct);
  }

  void processGlobalDistCorr(const DataT radius, const DataT phi, const DataT z0Tmp, const DataT z1Tmp, DataT& ddR, DataT& ddPhi, DataT& ddZ, const NumericalFields<DataT, RegularGrid>& formulaStruct) const
  {
    calcDistortions(radius, phi, z0Tmp, z1Tmp, ddR, ddPhi, ddZ, formulaStruct);
  }

  void processGlobalDistCorr(const DataT radius, const DataT phi, const DataT z0Tmp, const DataT z1Tmp, DataT& ddR, DataT& ddRPhi, DataT& ddZ, const DistCorrInterpolator<DataT, RegularGrid>& formulaStruct) const
  {
    ddR = formulaStruct.evaldR(z0Tmp, radius, phi);
    ddZ = formulaStruct.evaldZ(z0Tmp, radius, phi);
    ddRPhi = formulaStruct.evaldRPhi(z0Tmp, radius, phi) / radius;
  }
};

///
/// ========================================================================================================
///       Inline implementations of some methods
/// ========================================================================================================
///

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
template <typename ElectricFields>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::integrateEFieldsRoot(const DataT p1r, const DataT p1phi, const DataT p1z, const DataT p2z, DataT& localIntErOverEz, DataT& localIntEPhiOverEz, DataT& localIntDeltaEz, ElectricFields& formulaStruct) const
{
  const DataT ezField = getEzField();
  TF1 fErOverEz("fErOverEz", [&](double* x, double* p) { (void)p; return static_cast<double>(formulaStruct.evalEr(p1r, p1phi, static_cast<DataT>(x[0])) / (formulaStruct.evalEz(p1r, p1phi, static_cast<DataT>(x[0])) + ezField)); }, p1z, p2z, 1);
  localIntErOverEz = static_cast<DataT>(fErOverEz.Integral(p1z, p2z));

  TF1 fEphiOverEz("fEPhiOverEz", [&](double* x, double* p) { (void)p; return static_cast<double>(formulaStruct.evalEphi(p1r, p1phi, static_cast<DataT>(x[0])) / (formulaStruct.evalEz(p1r, p1phi, static_cast<DataT>(x[0])) + ezField)); }, p1z, p2z, 1);
  localIntEPhiOverEz = static_cast<DataT>(fEphiOverEz.Integral(p1z, p2z));

  TF1 fEz("fEZOverEz", [&](double* x, double* p) { (void)p; return static_cast<double>(formulaStruct.evalEz(p1r, p1phi, static_cast<DataT>(x[0])) - ezField); }, p1z, p2z, 1);
  localIntDeltaEz = static_cast<DataT>(fEz.Integral(p1z, p2z));
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
template <typename ElectricFields>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::integrateEFieldsTrapezoidal(const DataT p1r, const DataT p1phi, const DataT p1z, const DataT p2z, DataT& localIntErOverEz, DataT& localIntEPhiOverEz, DataT& localIntDeltaEz, ElectricFields& formulaStruct) const
{
  //========trapezoidal rule see: https://en.wikipedia.org/wiki/Trapezoidal_rule ==============
  const DataT ezField = getEzField();

  const DataT fielder0 = formulaStruct.evalEr(p1z, p1r, p1phi);
  const DataT fieldez0 = formulaStruct.evalEz(p1z, p1r, p1phi);
  const DataT fieldephi0 = formulaStruct.evalEphi(p1z, p1r, p1phi);

  const DataT fielder1 = formulaStruct.evalEr(p2z, p1r, p1phi);
  const DataT fieldez1 = formulaStruct.evalEz(p2z, p1r, p1phi);
  const DataT fieldephi1 = formulaStruct.evalEphi(p2z, p1r, p1phi);

  const DataT eZ0 = ezField + fieldez0;
  const DataT eZ1 = ezField + fieldez1;

  const unsigned int nSteps = mNumericalIntegrationSteps;
  const DataT deltaX = (p2z - p1z) / nSteps;

  DataT fieldSumEr = 0;
  DataT fieldSumEphi = 0;
  DataT fieldSumEz = 0;

  for (unsigned int i = 1; i < nSteps; ++i) {
    const DataT xk1Tmp = p1z + i * deltaX;
    const DataT ezField1 = formulaStruct.evalEz(xk1Tmp, p1r, p1phi);
    const DataT ezField1Denominator = 1 / ezField + ezField1;

    fieldSumEr += formulaStruct.evalEr(xk1Tmp, p1r, p1phi) * ezField1Denominator;
    fieldSumEphi += formulaStruct.evalEphi(xk1Tmp, p1r, p1phi) * ezField1Denominator;
    fieldSumEz += ezField1;
  }

  localIntErOverEz = deltaX * (fieldSumEr + static_cast<DataT>(0.5) * (fielder0 / eZ0 + fielder1 / eZ1));
  localIntEPhiOverEz = deltaX * (fieldSumEphi + static_cast<DataT>(0.5) * (fieldephi0 / eZ0 + fieldephi1 / eZ1));
  localIntDeltaEz = deltaX * (fieldSumEz + static_cast<DataT>(0.5) * (fieldez0 + fieldez1));
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
template <typename ElectricFields>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::integrateEFieldsSimpson(const DataT p1r, const DataT p1phi, const DataT p1z, const DataT p2z, DataT& localIntErOverEz, DataT& localIntEPhiOverEz, DataT& localIntDeltaEz, ElectricFields& formulaStruct) const
{
  //==========simpsons rule see: https://en.wikipedia.org/wiki/Simpson%27s_rule =============================
  const DataT ezField = getEzField();

  const DataT fielder0 = formulaStruct.evalEr(p1z, p1r, p1phi);
  const DataT fieldez0 = formulaStruct.evalEz(p1z, p1r, p1phi);
  const DataT fieldephi0 = formulaStruct.evalEphi(p1z, p1r, p1phi);

  const DataT fielder1 = formulaStruct.evalEr(p2z, p1r, p1phi);
  const DataT fieldez1 = formulaStruct.evalEz(p2z, p1r, p1phi);
  const DataT fieldephi1 = formulaStruct.evalEphi(p2z, p1r, p1phi);

  const DataT eZ0 = ezField + fieldez0;
  const DataT eZ1 = ezField + fieldez1;

  const unsigned int nSteps = mNumericalIntegrationSteps;
  const DataT deltaX = (p2z - p1z) / nSteps;

  DataT fieldSum1ErOverEz = 0;
  DataT fieldSum2ErOverEz = 0;
  DataT fieldSum1EphiOverEz = 0;
  DataT fieldSum2EphiOverEz = 0;
  DataT fieldSum1Ez = 0;
  DataT fieldSum2Ez = 0;

  for (unsigned int i = 1; i < nSteps; ++i) {
    const DataT xk1Tmp = p1z + i * deltaX;
    const DataT xk2 = xk1Tmp - static_cast<DataT>(0.5) * deltaX;

    const DataT ezField1 = formulaStruct.evalEz(xk1Tmp, p1r, p1phi);
    const DataT ezField2 = formulaStruct.evalEz(xk2, p1r, p1phi);
    const DataT ezField1Denominator = 1 / (ezField + ezField1);
    const DataT ezField2Denominator = 1 / (ezField + ezField2);

    fieldSum1ErOverEz += formulaStruct.evalEr(xk1Tmp, p1r, p1phi) * ezField1Denominator;
    fieldSum2ErOverEz += formulaStruct.evalEr(xk2, p1r, p1phi) * ezField2Denominator;

    fieldSum1EphiOverEz += formulaStruct.evalEphi(xk1Tmp, p1r, p1phi) * ezField1Denominator;
    fieldSum2EphiOverEz += formulaStruct.evalEphi(xk2, p1r, p1phi) * ezField2Denominator;

    fieldSum1Ez += ezField1;
    fieldSum2Ez += ezField2;
  }
  const DataT xk2N = (p2z - static_cast<DataT>(0.5) * deltaX);
  const DataT ezField2 = formulaStruct.evalEz(xk2N, p1r, p1phi);
  const DataT ezField2Denominator = 1 / (ezField + ezField2);
  fieldSum2ErOverEz += formulaStruct.evalEr(xk2N, p1r, p1phi) * ezField2Denominator;
  fieldSum2EphiOverEz += formulaStruct.evalEphi(xk2N, p1r, p1phi) * ezField2Denominator;
  fieldSum2Ez += ezField2;

  const DataT deltaXSimpsonSixth = deltaX / 6;
  localIntErOverEz = deltaXSimpsonSixth * (2 * fieldSum1ErOverEz + 4 * fieldSum2ErOverEz + fielder0 / eZ0 + fielder1 / eZ1);
  localIntEPhiOverEz = deltaXSimpsonSixth * (2 * fieldSum1EphiOverEz + 4 * fieldSum2EphiOverEz + fieldephi0 / eZ0 + fieldephi1 / eZ1);
  localIntDeltaEz = deltaXSimpsonSixth * (2 * fieldSum1Ez + 4 * fieldSum2Ez + fieldez0 + fieldez1);
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
template <typename ElectricFields>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::integrateEFieldsSimpsonExperimental(const DataT p1r, const DataT p2r, const DataT p1phi, const DataT p2phi, const DataT p1z, const DataT p2z, DataT& localIntErOverEz, DataT& localIntEPhiOverEz, DataT& localIntDeltaEz, ElectricFields& formulaStruct) const
{
  //==========simpsons rule see: https://en.wikipedia.org/wiki/Simpson%27s_rule =============================
  const DataT ezField = getEzField();

  const DataT fielder0 = formulaStruct.evalEr(p1z, p1r, p1phi);
  const DataT fieldez0 = formulaStruct.evalEz(p1z, p1r, p1phi);
  const DataT fieldephi0 = formulaStruct.evalEphi(p1z, p1r, p1phi);

  const DataT fielder1 = formulaStruct.evalEr(p2z, p2r, p2phi);
  const DataT fieldez1 = formulaStruct.evalEz(p2z, p2r, p2phi);
  const DataT fieldephi1 = formulaStruct.evalEphi(p2z, p2r, p2phi);

  const DataT eZ0 = ezField + fieldez0;
  const DataT eZ1 = ezField + fieldez1;

  const unsigned int nSteps = mNumericalIntegrationSteps;
  const DataT deltaX = (p2z - p1z) / nSteps;

  DataT fieldSum1ErOverEz = 0;
  DataT fieldSum2ErOverEz = 0;
  DataT fieldSum1EphiOverEz = 0;
  DataT fieldSum2EphiOverEz = 0;
  DataT fieldSum1Ez = 0;
  DataT fieldSum2Ez = 0;

  const DataT xk2N = (p2z - static_cast<DataT>(0.5) * deltaX);
  const DataT ezField2 = formulaStruct.evalEz(xk2N, 0.5 * (p1r + p2r), 0.5 * (p1phi + p2phi));
  const DataT ezField2Denominator = 1 / (ezField + ezField2);
  fieldSum2ErOverEz += formulaStruct.evalEr(xk2N, 0.5 * (p1r + p2r), 0.5 * (p1phi + p2phi)) * ezField2Denominator;
  fieldSum2EphiOverEz += formulaStruct.evalEphi(xk2N, 0.5 * (p1r + p2r), 0.5 * (p1phi + p2phi)) * ezField2Denominator;
  fieldSum2Ez += ezField2;

  const DataT deltaXSimpsonSixth = deltaX / 6;
  localIntErOverEz = deltaXSimpsonSixth * (2 * fieldSum1ErOverEz + 4 * fieldSum2ErOverEz + fielder0 / eZ0 + fielder1 / eZ1);
  localIntEPhiOverEz = deltaXSimpsonSixth * (2 * fieldSum1EphiOverEz + 4 * fieldSum2EphiOverEz + fieldephi0 / eZ0 + fieldephi1 / eZ1);
  localIntDeltaEz = deltaXSimpsonSixth * (2 * fieldSum1Ez + 4 * fieldSum2Ez + fieldez0 + fieldez1);
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
template <typename ElectricFields>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::calcLocalDistortionsCorrections(const bool lcorrections, ElectricFields& formulaStruct)
{
#pragma omp parallel for
  for (size_t iPhi = 0; iPhi < Nphi; iPhi++) {
    const DataT phi = getPhiVertex(iPhi);
    for (size_t iR = 0; iR < Nr; iR++) {
      const DataT radius = getRVertex(iR);
      for (size_t iZ = 0; iZ < Nz - 1; ++iZ) {

        // set z coordinated depending on distortions or correction calculation
        const DataT z0 = lcorrections ? getZVertex(iZ + 1) : getZVertex(iZ);
        const DataT z1 = lcorrections ? getZVertex(iZ) : getZVertex(iZ + 1);

        const int iSteps = getIntegrationSteps();
        const DataT stepSize = (z1 - z0) / iSteps;

        DataT drTmp = 0;
        DataT dRPhiTmp = 0;
        DataT dPhiTmp = 0;
        DataT dzTmp = 0;

        for (int iter = 0; iter < iSteps; ++iter) {
          const DataT z0Tmp = regulateZ(z0 + iter * stepSize + dzTmp);
          const DataT z1Tmp = regulateZ(z0Tmp + stepSize);

          DataT ddR = 0;
          DataT ddPhi = 0;
          DataT ddZ = 0;
          const DataT radiusTmp = regulateR(radius + drTmp);
          const DataT phiTmp = regulatePhi(phi + dPhiTmp);

          calcDistortions(radiusTmp, phiTmp, z0Tmp, z1Tmp, ddR, ddPhi, ddZ, formulaStruct);
          drTmp += ddR;
          dRPhiTmp += ddPhi * radiusTmp;
          dPhiTmp += ddPhi;
          dzTmp += ddZ;
        }
        if (lcorrections == true) {
          mLocalCorrdR(iZ + 1, iR, iPhi) = drTmp;
          mLocalCorrdRPhi(iZ + 1, iR, iPhi) = dRPhiTmp;
          mLocalCorrdZ(iZ + 1, iR, iPhi) = dzTmp;
        } else {
          mLocalDistdR(iZ, iR, iPhi) = drTmp;
          mLocalDistdRPhi(iZ, iR, iPhi) = dRPhiTmp;
          mLocalDistdZ(iZ, iR, iPhi) = dzTmp;
        }
      }
    }
  }
  if (lcorrections) {
    mIsLocalCorrSet = true;
  } else {
    mIsLocalDistSet = true;
  }
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
template <typename ElectricFields>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::calcDistortions(const DataT p1r, const DataT p1phi, const DataT p1z, const DataT p2z, DataT& ddR, DataT& ddRPhi, DataT& ddZ, ElectricFields& formulaStruct) const
{
  DataT localIntErOverEz = 0;
  DataT localIntEPhiOverEz = 0;
  DataT localIntDeltaEz = 0;

  const int iN = mNumericalIntegrationStrategy == SimpsonExperimental ? 10 : 1;
  for (int i = 0; i < iN; ++i) {
    switch (mNumericalIntegrationStrategy) {
      case Simpson:
        integrateEFieldsSimpson(p1r, p1phi, p1z, p2z, localIntErOverEz, localIntEPhiOverEz, localIntDeltaEz, formulaStruct);
        break;
      case SimpsonExperimental:
        integrateEFieldsSimpsonExperimental(p1r, regulateR(p1r + ddR), p1phi, p1phi + ddRPhi, p1z, regulateZ(p2z + ddZ), localIntErOverEz, localIntEPhiOverEz, localIntDeltaEz, formulaStruct);
        break;
      case Trapezoidal:
        integrateEFieldsTrapezoidal(p1r, p1phi, p1z, p2z, localIntErOverEz, localIntEPhiOverEz, localIntDeltaEz, formulaStruct);
        break;
      case Root:
        integrateEFieldsRoot(p1r, p1phi, p1z, p2z, localIntErOverEz, localIntEPhiOverEz, localIntDeltaEz, formulaStruct);
        break;
      default:
        std::cout << "no matching case: Using Simpson" << std::endl;
        integrateEFieldsSimpson(p1r, p1phi, p1z, p2z, localIntErOverEz, localIntEPhiOverEz, localIntDeltaEz, formulaStruct);
    }
    ddR = fC0 * localIntErOverEz + fC1 * localIntEPhiOverEz;
    ddRPhi = (fC0 * localIntEPhiOverEz - fC1 * localIntErOverEz) / p1r;
    ddZ = -1 * localIntDeltaEz * ASolv::fgkdvdE;
  }
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
template <typename ElectricFields>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::calcGlobalDistortions(const ElectricFields& formulaStruct)
{
// loop over tpc volume and let the electron drift from each vertex tothe readout of the tpc
#pragma omp parallel for
  for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
    for (size_t iR = 0; iR < Nr; ++iR) {
      for (size_t iZ = 0; iZ < Nz - 1; ++iZ) {

        DataT drDist = 0.0;
        DataT dPhiDist = 0.0;
        DataT dzDist = 0.0;

        int iter = 0;
        for (;;) {
          // const DataT z0 = getZVertex(iZ); // the electron starts at phi, radius, z0
          // const DataT z1 = getZVertex(iZ + 1);
          const DataT stepSize = formulaStruct.id == 2 ? mGridSpacingZ : mGridSpacingZ / getIntegrationSteps();

          const DataT z0Tmp = getZVertex(iZ) + dzDist + iter * stepSize;

          const DataT z1Tmp = regulateZ(z0Tmp + stepSize);
          const DataT radius = regulateR(getRVertex(iR) + drDist);
          const DataT phi = regulatePhi(getPhiVertex(iPhi) + dPhiDist);

          DataT ddR = 0;
          DataT ddPhi = 0;
          DataT ddZ = 0;

          // get the distortion from interpolation of local distortions or electric field
          processGlobalDistCorr(radius, phi, z0Tmp, z1Tmp, ddR, ddPhi, ddZ, formulaStruct);

          // add local distortion
          drDist += ddR;
          dPhiDist += ddPhi;
          dzDist += ddZ;

          if (z0Tmp + stepSize >= ASolv::fgkTPCZ0) { // set loop to exit if teh readout is reached
            break;
          }
          ++iter;
        }
        mGlobalDistdR(iZ, iR, iPhi) = drDist;
        mGlobalDistdRPhi(iZ, iR, iPhi) = dPhiDist * getRVertex(iR);
        mGlobalDistdZ(iZ, iR, iPhi) = dzDist;
      }
    }
  }
  mIsGlobalDistSet = true;
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
template <typename ElectricFields>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::calcGlobalCorrections(const ElectricFields& formulaStruct)
{
#pragma omp parallel for
  for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
    for (size_t iR = 0; iR < Nr; ++iR) {
      DataT drCorr = 0;
      DataT dPhiCorr = 0;
      DataT dzCorr = 0;
      for (size_t iZ = Nz - 1; iZ >= 1; --iZ) {

        const int iSteps = formulaStruct.id == 2 ? 1 : getIntegrationSteps();
        for (int iter = 0; iter < iSteps; ++iter) {
          const DataT z0 = getZVertex(iZ); // the electron starts at phi, radius, z0
          // const DataT z1 = getZVertex(iZ - 1);

          const DataT stepSize = -mGridSpacingZ / iSteps;
          const DataT radius = regulateR(getRVertex(iR) + drCorr);
          const DataT phi = regulatePhi(getPhiVertex(iPhi) + dPhiCorr);
          const DataT z0Tmp = regulateZ(z0 + dzCorr + iter * stepSize);
          const DataT z1Tmp = regulateZ(z0Tmp + stepSize);

          DataT ddR = 0;
          DataT ddPhi = 0;
          DataT ddZ = 0;

          // get the distortion from interpolation of local distortions or electric field
          processGlobalDistCorr(radius, phi, z0Tmp, z1Tmp, ddR, ddPhi, ddZ, formulaStruct);

          drCorr += ddR;
          dPhiCorr += ddPhi;
          dzCorr += ddZ;
        }

        mGlobalCorrdR(iZ - 1, iR, iPhi) = drCorr;
        mGlobalCorrdRPhi(iZ - 1, iR, iPhi) = dPhiCorr * getRVertex(iR);
        mGlobalCorrdZ(iZ - 1, iR, iPhi) = dzCorr;
      }
    }
  }
  mIsGlobalCorrSet = true;
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::calcGlobalDistWithGlobalCorrIterative(const DistCorrInterpolator<DataT, RegularGrid>& globCorr, const int maxIter, const DataT convZ, const DataT convR, const DataT convPhi)
{
  // store all values here for kdtree
  const int nPoints = Nz * Nr * Nphi;
  std::vector<gte::PositionSite<3, float>> sites;
  sites.reserve(nPoints);

  for (unsigned int iPhi = 0; iPhi < Nphi; ++iPhi) {
    for (unsigned int iR = 0; iR < Nr; ++iR) {
      for (unsigned int iZ = 0; iZ < Nz; ++iZ) {
        const DataT radius = getRVertex(iR);
        const DataT z = getZVertex(iZ);
        const DataT phi = getPhiVertex(iPhi);

        const DataT globalCorrR = mGlobalCorrdR(iZ, iR, iPhi);
        const DataT globalCorrRPhi = mGlobalCorrdRPhi(iZ, iR, iPhi);
        const DataT globalCorrZ = mGlobalCorrdZ(iZ, iR, iPhi);

        const DataT posRCorr = radius + globalCorrR; // position of global correction
        const DataT posPhiCorr = regulatePhi(phi + globalCorrRPhi / radius);
        const DataT posZCorr = z + globalCorrZ;

        if (posRCorr >= mRMin && posRCorr <= ASolv::fgkOFCRadius && posZCorr >= mZMin && posZCorr <= ASolv::fgkTPCZ0) {
          const std::array<float, 3> position{static_cast<float>((posZCorr - mZMin) * getInvSpacingZ()), static_cast<float>((posRCorr - mRMin) * getInvSpacingR()), static_cast<float>((posPhiCorr - mPhiMin) * getInvSpacingPhi())};
          const std::array<unsigned int, 3> positionIndex{iZ, iR, iPhi};
          const gte::PositionSite<3, float> siteTmp(position, positionIndex);
          sites.emplace_back(siteTmp);
        }
      }
    }
  }

  const int maxLeafSize = 10;
  const int maxLevel = 10;
  gte::NearestNeighborQuery<3, float, gte::PositionSite<3, float>> kdTree(sites, maxLeafSize, maxLevel);

#pragma omp parallel for
  for (unsigned int iZ = 0; iZ < Nz; ++iZ) {
    for (unsigned int iR = 0; iR < Nr; ++iR) {
      for (unsigned int iPhi = 0; iPhi < Nphi; ++iPhi) {
        // find nearest neighbour
        const DataT radius = getRVertex(iR);
        const DataT z = getZVertex(iZ);
        const DataT phi = getPhiVertex(iPhi);

        const std::array<float, 3> valuesQuery = {static_cast<float>((z - mZMin) * getInvSpacingZ()), static_cast<float>((radius - mRMin) * getInvSpacingR()), static_cast<float>((phi - mPhiMin) * getInvSpacingPhi())};
        const gte::Vector<3, float> point(valuesQuery);
        const float radiusSearch = 3; // larger radius -> more cpu time

        const int MaxNeighbors = 1;
        std::array<int, MaxNeighbors> neighbors{};
        const int nNeighbours = kdTree.FindNeighbors<MaxNeighbors>(point, radiusSearch, neighbors);

        if (nNeighbours == 0) {
          // TODO make automatic search radius
          std::cout << "no Neighbour found :( use larger search radius!" << std::endl;
          continue;
        }

        const unsigned int index = neighbors[0];
        const DataT nearestZ = sites[index].position[0] * mGridSpacingZ + mZMin;
        const DataT nearestR = sites[index].position[1] * mGridSpacingR + mRMin;
        const DataT nearestPhi = sites[index].position[2] * mGridSpacingPhi + mPhiMin;

        const unsigned int nearestiZ = sites[index].index[0];
        const unsigned int nearestiR = sites[index].index[1];
        const unsigned int nearestiPhi = sites[index].index[2];

        //start algorithm: use tricubic upsampling to numerically approach the query point
        // 1. calculate difference from nearest point to query point with stepwidth factor x
        const DataT rStepWidth = static_cast<DataT>(0.1);
        DataT stepR = (radius - nearestR) * rStepWidth;
        const DataT zStepWidth = static_cast<DataT>(0.1);
        DataT stepZ = (z - nearestZ) * zStepWidth;
        const DataT phiStepWidth = static_cast<DataT>(0.1);
        DataT stepPhi = (phi - nearestPhi) * phiStepWidth;

        // needed to check for convergence
        DataT lastDistanceR = std::numeric_limits<DataT>::max();
        DataT lastDistanceZ = std::numeric_limits<DataT>::max();
        DataT lastDistancePhi = std::numeric_limits<DataT>::max();

        int count = 0;

        DataT corrdR = 0;
        DataT corrdRPhi = 0;
        DataT corrdZ = 0;

        const bool safe = false;

        for (int iter = 0; iter < maxIter; ++iter) {
          // 2. get new points coordinates
          const DataT rPos = regulateR(getRVertex(nearestiR) + stepR);
          const DataT zPos = regulateZ(getZVertex(nearestiZ) + stepZ);
          const DataT phiPosUnreg = getPhiVertex(nearestiPhi) + stepPhi;
          const DataT phiPos = regulatePhi(phiPosUnreg);

          corrdR = globCorr.evaldR(zPos, rPos, phiPos, safe);
          const DataT rpos = rPos + corrdR;

          const DataT corrRPhi = globCorr.evaldRPhi(zPos, rPos, phiPos, safe);
          corrdRPhi = corrRPhi / rPos * rpos;
          const DataT phipos = phiPosUnreg + corrRPhi / rPos;

          corrdZ = globCorr.evaldZ(zPos, rPos, phiPos, safe);
          const DataT zpos = zPos + corrdZ;

          const DataT distanceR = radius - rpos;
          const DataT distanceZ = z - zpos;
          const DataT distancePhi = phi - phipos;
          stepR += distanceR * rStepWidth;
          stepZ += distanceZ * zStepWidth;
          stepPhi += distancePhi * phiStepWidth;

          const DataT totaldistRDiv = lastDistanceR == 0 ? 0 : std::abs(1 - std::abs(distanceR / lastDistanceR)); // should be larger than 0
          const bool checkR = totaldistRDiv <= convR;                                                             // if the improvemnt in distance is smaller than changeRFactor set the flag

          const DataT totaldistZDiv = lastDistanceZ == 0 ? 0 : std::abs(1 - std::abs(distanceZ / lastDistanceZ)); // should be larger than 0
          const bool checkZ = totaldistZDiv <= convZ;                                                             // if the improvemnt in distance is smaller than changeRFactor set the flag

          const DataT totaldistPhiDiv = lastDistancePhi == 0 ? 0 : std::abs(1 - std::abs(distancePhi / lastDistancePhi)); // should be larger than 0
          const bool checkPhi = totaldistPhiDiv <= convPhi;                                                               // if the improvemnt in distance is smaller than changeRFactor set the flag

          if (checkR && checkZ && checkPhi) {
            break;
          }

          lastDistanceR = distanceR;
          lastDistanceZ = distanceZ;
          lastDistancePhi = distancePhi;
          ++count;
        }
        mGlobalDistdR(iZ, iR, iPhi) = -corrdR;
        mGlobalDistdRPhi(iZ, iR, iPhi) = -corrdRPhi;
        mGlobalDistdZ(iZ, iR, iPhi) = -corrdZ;
      }
    }
  }
}

#endif
