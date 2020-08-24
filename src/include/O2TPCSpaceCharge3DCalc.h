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
#include "O2TPCPoissonSolver.h"

#include "SpaceChargeStructs.h"
#include "O2/Defs.h"
#include "RegularGrid3D.h"

// Root includes
#include "TF1.h"   /// for numerical intergration only
#include "TTree.h" /// for debugging
#include "TH3.h"
#include <chrono>

#ifdef WITH_OPENMP
#include <omp.h>
#endif

// for nearest neighbour search
// #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Search_traits_3.h>
// #include <CGAL/Search_traits_adapter.h>
// #include <CGAL/Orthogonal_k_neighbor_search.h>
// #include <CGAL/property_map.h>
// #include <boost/iterator/zip_iterator.hpp>

namespace o2
{
namespace tpc
{

/// \tparam DataT the type of data which is used during the calculations
/// \tparam Nr number of vertices in r direction
/// \tparam Nz number of vertices in z direction
/// \tparam Nphi number of vertices in phi direction
template <typename DataT = float, size_t Nz = 17, size_t Nr = 17, size_t Nphi = 90>
class O2TPCSpaceCharge3DCalc
{
  using RegularGrid = RegularGrid3D<DataT, Nz, Nr, Nphi>;
  using DataContainer = DataContainer3D<DataT, Nz, Nr, Nphi>;

 public:
  O2TPCSpaceCharge3DCalc() = default;

  // rebin the input space charge density histogram to desired binning
  void rebinDensityHisto(const TH3& hOrig, TH3& hRebined) const;

  // set the charge density from TH3 histogram containing the space charge density
  void fillChargeDensityFromHisto(TFile& fInp, const char* name);

  // fill boundary and ChargeDensities
  // posson solver
  // local distortions and corrections
  // global distortions and correction
  // \param mode mode=0: full analytical, mode=1 full numerical
  void performFullRun(const AnalyticalFields<DataT>& formulaStruct, const int mode, const bool electricFieldGlobCorrDist, TFile& file, const o2::tpc::Side side);
  void performGlobalCorrDist(TFile& file, const o2::tpc::Side side);

  void setFromFile(TFile& file, const o2::tpc::Side side);

  // stepp 0: this function fills the internal storage for density and boundary conditions for potential
  /// \param formulaStruct struct containing a method to evaluate the density and potential
  void fillBoundaryAndChargeDensities(const AnalyticalFields<DataT>& formulaStruct);

  // stepp 1: use the AliTPCPoissonSolver class to numerically calculate the potential with space charge density and boundary conditions from potential
  void poissonSolver(const o2::tpc::Side side, const int maxIteration = 300, const DataT stoppingConvergence = 1e-8);

  // stepp 2: calculate numerically the electric field from the potential
  void calcEField(const o2::tpc::Side side);

  /// set electric field from analytic formula
  /// \param formulaStruct struct containing a method to evaluate the density and potential
  void setEField(const AnalyticalFields<DataT>& formulaStruct);

  // stepp 3:
  /// \param type calculate local corrections or local distortions: type=0->distortions, type=1->corrections
  /// \param formulaStruct struct containing a method to evaluate the electric field Er, Ez, Ephi
  template <typename ElectricFields = AnalyticalFields<DataT>>
  void calcLocalDistortionsCorrections(const int type, ElectricFields& formulaStruct);

  //step 4: calculate global distortions by using the electric field or the local distortions
  /// \param formulaStruct struct containing a method to evaluate the electric field Er, Ez, Ephi or the local distortions
  template <typename Formulas = AnalyticalFields<DataT>>
  void calcGlobalDistortions(const Formulas& formulaStruct);

  //step 4: calculate global corrections by using the electric field or the local corrections
  /// \param formulaStruct struct containing a method to evaluate the electric field Er, Ez, Ephi or the local corrections
  template <typename Formulas = AnalyticalFields<DataT>>
  void calcGlobalCorrections(const Formulas& formulaStruct);

  // calculate global distortions using the global corrections
  /// \param globCorr interpolator of global corrections
  /// \param maxIter maximum iterations per global distortion
  /// \param convZ convergence criteria for z direction: if the ratio from the position from last iteration zPosLast compared to positon from current iteration zPosCurr is smaller than this value set converged to true abs(1-abs(zPosLast/zPosCurr))<convZ
  /// \param convR convergence criteria for r direction: if the ratio from the position from last iteration rPosLast compared to positon from current iteration rPosCurr is smaller than this value set converged to true abs(1-abs(rPosLast/rPosCurr))<convR
  /// \param convPhi convergence criteria for phi direction: if the ratio from the position from last iteration phiPosLast compared to positon from current iteration phiPosCurr is smaller than this value set converged to true abs(1-abs(PhiPosLast/PhiPosCurr))<convPhi
  /// \param approachZ when the difference between the desired z coordinate and the position of the global correction is deltaZ, approach the desired z coordinate by deltaZ * \p approachZ.
  /// \param approachR when the difference between the desired r coordinate and the position of the global correction is deltaR, approach the desired r coordinate by deltaR * \p approachR.
  /// \param approachPhi when the difference between the desired phi coordinate and the position of the global correction is deltaPhi, approach the desired phi coordinate by deltaPhi * \p approachPhi.
  void calcGlobalDistWithGlobalCorrIterative(const DistCorrInterpolator<DataT, Nz, Nr, Nphi>& globCorr, const int maxIter = 200, const DataT convZ = 0.05, const DataT convR = 0.05, const DataT convPhi = 0.05, const DataT approachZ = 0.1, const DataT approachR = 0.1, const DataT approachPhi = 0.1);

  /// Get grid spacing in r direction
  static constexpr DataT getGridSpacingR() { return GRIDSPACINGR; }

  /// Get grid spacing in z direction
  static constexpr DataT getGridSpacingZ() { return GRIDSPACINGZ; }

  /// Get grid spacing in phi direction
  static constexpr DataT getGridSpacingPhi() { return GRIDSPACINGPHI; }

  /// Get constant electric field
  static constexpr DataT getEzField() { return (TPCParameters<DataT>::CATHODEV - TPCParameters<DataT>::GG) / TPCParameters<DataT>::TPCZ0; }

  /// Get inner radius of tpc
  static constexpr DataT getRMin() { return RMIN; }

  /// Get min z position which is used during the calaculations
  static constexpr DataT getZMin() { return ZMIN; }

  /// Get min phi
  static constexpr DataT getPhiMin() { return PHIMIN; }

  /// Get max r
  static constexpr DataT getRMax() { return RMAX; };

  /// Get max z
  static constexpr DataT getZMax() { return ZMAX; }

  /// Get max phi
  static constexpr DataT getPhiMax() { return PHIMAX; }

  // get side of TPC for z coordinate TODO rewrite this
  int getSide(const DataT z) const { return z <= 0 ? 0 : 1; }

  /// Get the grid object
  const RegularGrid& getGrid3D() const { return mGrid3D; }

  /// Get struct containing interpolators for electrical fields
  NumericalFields<DataT, Nz, Nr, Nphi> getElectricFieldsInterpolator(const o2::tpc::Side side) const;

  /// Get struct containing interpolators for local distortions dR, dZ, dPhi
  DistCorrInterpolator<DataT, Nz, Nr, Nphi> getLocalDistInterpolator(const o2::tpc::Side side) const;

  /// Get struct containing interpolators for local corrections dR, dZ, dPhi
  DistCorrInterpolator<DataT, Nz, Nr, Nphi> getLocalCorrInterpolator(const o2::tpc::Side side) const;

  /// Get struct containing interpolators for global distortions dR, dZ, dPhi
  DistCorrInterpolator<DataT, Nz, Nr, Nphi> getGlobalDistInterpolator(const o2::tpc::Side side) const;

  /// Get struct containing interpolators for global corrections dR, dZ, dPhi
  DistCorrInterpolator<DataT, Nz, Nr, Nphi> getGlobalCorrInterpolator(const o2::tpc::Side side) const;

  /// \param vertex in iz dimension
  /// \param vertex in ir dimension
  /// \param vertex in iphi dimension
  /// \return returns local distortion dR for given vertex
  DataT getLocalDistR(const size_t iz, const size_t ir, const size_t iphi, const o2::tpc::Side side) const { return mLocalDistdR[side](iz, ir, iphi); }

  /// \param vertex in iz dimension
  /// \param vertex in ir dimension
  /// \param vertex in iphi dimension
  /// \return returns local distortion dZ for given vertex
  DataT getLocalDistZ(const size_t iz, const size_t ir, const size_t iphi, const o2::tpc::Side side) const { return mLocalDistdZ[side](iz, ir, iphi); }

  /// \param vertex in iz dimension
  /// \param vertex in ir dimension
  /// \param vertex in iphi dimension
  /// \return returns local distortion dRPhi for given vertex
  DataT getLocalDistRPhi(const size_t iz, const size_t ir, const size_t iphi, const o2::tpc::Side side) const { return mLocalDistdRPhi[side](iz, ir, iphi); }

  /// \param vertex in iz dimension
  /// \param vertex in ir dimension
  /// \param vertex in iphi dimension
  /// \return returns local correction dR for given vertex
  DataT getLocalCorrR(const size_t iz, const size_t ir, const size_t iphi, const o2::tpc::Side side) const { return mLocalCorrdR[side](iz, ir, iphi); }

  /// \param vertex in iz dimension
  /// \param vertex in ir dimension
  /// \param vertex in iphi dimension
  /// \return returns local correction dZ for given vertex
  DataT getLocalCorrZ(const size_t iz, const size_t ir, const size_t iphi, const o2::tpc::Side side) const { return mLocalCorrdZ[side](iz, ir, iphi); }

  /// \param vertex in iz dimension
  /// \param vertex in ir dimension
  /// \param vertex in iphi dimension
  /// \return returns local correction dRPhi for given vertex
  DataT getLocalCorrRPhi(const size_t iz, const size_t ir, const size_t iphi, const o2::tpc::Side side) const { return mLocalCorrdRPhi[side](iz, ir, iphi); }

  /// Get global distortion dR for vertex
  DataT getGlobalDistR(const size_t iz, const size_t ir, const size_t iphi, const o2::tpc::Side side) const { return mGlobalDistdR[side](iz, ir, iphi); }

  /// Get global distortion dZ for vertex
  DataT getGlobalDistZ(const size_t iz, const size_t ir, const size_t iphi, const o2::tpc::Side side) const { return mGlobalDistdZ[side](iz, ir, iphi); }

  /// Get global distortion dRPhi for vertex
  DataT getGlobalDistRPhi(const size_t iz, const size_t ir, const size_t iphi, const o2::tpc::Side side) const { return mGlobalDistdRPhi[side](iz, ir, iphi); }

  /// Get global correction dR for vertex
  DataT getGlobalCorrR(const size_t iz, const size_t ir, const size_t iphi, const o2::tpc::Side side) const { return mGlobalCorrdR[side](iz, ir, iphi); }

  /// Get global correction dZ for vertex
  DataT getGlobalCorrZ(const size_t iz, const size_t ir, const size_t iphi, const o2::tpc::Side side) const { return mGlobalCorrdZ[side](iz, ir, iphi); }

  /// Get global correction dRPhi for vertex
  DataT getGlobalCorrRPhi(const size_t iz, const size_t ir, const size_t iphi, const o2::tpc::Side side) const { return mGlobalCorrdRPhi[side](iz, ir, iphi); }

  /// Get global electric Field Er for vertex
  DataT getEr(const size_t iz, const size_t ir, const size_t iphi, const o2::tpc::Side side) const { return mElectricFieldEr[side](iz, ir, iphi); }

  /// Get global electric Field Ez for vertex
  DataT getEz(const size_t iz, const size_t ir, const size_t iphi, const o2::tpc::Side side) const { return mElectricFieldEz[side](iz, ir, iphi); }

  /// Get global electric Field Ephi for vertex
  DataT getEphi(const size_t iz, const size_t ir, const size_t iphi, const o2::tpc::Side side) const { return mElectricFieldEphi[side](iz, ir, iphi); }

  /// Get density for vertex
  DataT getDensity(const size_t iz, const size_t ir, const size_t iphi, const o2::tpc::Side side) const { return mDensity[side](iz, ir, iphi); }

  /// Get the step width which is used for the calculation of the correction/distortions in units of the z-bin
  int getStepWidth() const { return 1 / mStepWidth; }

  /// Get phi vertex psotion for index in phi direction
  DataT getPhiVertex(const size_t indexPhi) const { return mGrid3D.getZVertex(indexPhi); }

  /// Get r vertex psotion for index in r direction
  DataT getRVertex(const size_t indexR) const { return mGrid3D.getYVertex(indexR); }

  /// Get z vertex psotion for index in z direction
  DataT getZVertex(const size_t indexZ) const { return mGrid3D.getXVertex(indexZ); }

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

  void setNStep(const int nSteps) { mStepWidth = nSteps; }

  void setNumericalIntegrationStrategy(int strategy)
  {
    mNumericalIntegrationStrategy = strategy;
  }

  /// numerical integration strategys
  enum IntegrationStrategy { Trapezoidal = 0,        ///< trapezoidal integration (https://en.wikipedia.org/wiki/Trapezoidal_rule). straight electron drift line assumed: z0->z1, r0->r0, phi0->phi0
                             Simpson = 1,            ///< simpon integration. see: https://en.wikipedia.org/wiki/Simpson%27s_rule. straight electron drift line assumed: z0->z1, r0->r0, phi0->phi0
                             Root = 2,               ///< Root integration. straight electron drift line assumed: z0->z1, r0->r0, phi0->phi0
                             SimpsonExperimental = 3 ///< simpon integration, but using an iterative method to approximate the drift path. No straight electron drift line assumed: z0->z1, r0->r1, phi0->phi1
  };

  /// write electric field to root file
  int dumpElectricFields(TFile& outf, const o2::tpc::Side side) const;

  /// set electric field from root file
  void setElectricFieldsFromFile(TFile& inpf, const o2::tpc::Side side);

  /// write potential to root file
  int dumpPotential(TFile& outf, const o2::tpc::Side side) const { return mPotential[side].writeToFile(outf, Form("potential_side%s", getSideName(side).data())); }

  /// set potential from root file
  void setPotentialFromFile(TFile& inpf, const o2::tpc::Side side) { mPotential[side].initFromFile(inpf, Form("potential_side%s", getSideName(side).data())); }

  /// write potential to root file
  int dumpDensity(TFile& outf, const o2::tpc::Side side) const { return mDensity[side].writeToFile(outf, Form("density_side%s", getSideName(side).data())); }

  /// set potential from root file
  void setDensityFromFile(TFile& inpf, const o2::tpc::Side side) { mDensity[side].initFromFile(inpf, Form("density_side%s", getSideName(side).data())); }

  /// write global distortions to root file
  int dumpGlobalDistortions(TFile& outf, const o2::tpc::Side side) const;

  /// set global distortions from root file
  void setGlobalDistortionsFromFile(TFile& inpf, const o2::tpc::Side side);

  /// write global corrections to root file
  int dumpGlobalCorrections(TFile& outf, const o2::tpc::Side side) const;

  /// set global corrections from root file
  void setGlobalCorrectionsFromFile(TFile& inpf, const o2::tpc::Side side);

  /// write local corrections to root file
  int dumpLocalCorrections(TFile& outf, const o2::tpc::Side side) const;

  /// set local corrections from root file
  void setLocalCorrectionsFromFile(TFile& inpf, const o2::tpc::Side side);

  /// write local distortions to root file
  int dumpLocalDistortions(TFile& outf, const o2::tpc::Side side) const;

  /// set local distortions from root file
  void setLocalDistortionsFromFile(TFile& inpf, const o2::tpc::Side side);

  /// set z coordinate between min z max z
  DataT regulateZ(const DataT posZ) const { return mGrid3D.clampToGrid(posZ, 0); }

  /// set r coordinate between min r max r
  DataT regulateR(const DataT posR) const { return mGrid3D.clampToGrid(posR, 1); }

  /// set phi coordinate between min phi max phi
  DataT regulatePhi(const DataT posPhi) const { return mGrid3D.clampToGridCircular(posPhi, 2); }

 private:
  using ASolvAli = AliTPCPoissonSolver<DataT>;
  using ASolv = o2::tpc::O2TPCPoissonSolver<DataT, Nz, Nr, Nphi>;

  static constexpr DataT RMIN = TPCParameters<DataT>::IFCRADIUS;  ///< min radius
  static constexpr DataT ZMIN = 0;                                ///< min z coordinate
  static constexpr DataT PHIMIN = 0;                              ///< min phi coordinate
  static constexpr DataT RMAX = TPCParameters<DataT>::OFCRADIUS;  ///< max radius
  static constexpr DataT ZMAX = TPCParameters<DataT>::TPCZ0;      ///< max z coordinate
  static constexpr DataT PHIMAX = 2 * M_PI;                       ///< max phi coordinate // TODO CHANGE TO O2
  static constexpr DataT GRIDSPACINGR = (RMAX - RMIN) / (Nr - 1); ///< grid spacing in r direction
  static constexpr DataT GRIDSPACINGZ = (ZMAX - ZMIN) / (Nz - 1); ///< grid spacing in z direction
  static constexpr DataT GRIDSPACINGPHI = PHIMAX / Nphi;          ///< grid spacing in phi direction

  int mNumericalIntegrationStrategy = SimpsonExperimental; ///< numerical integration strategy of integration of the E-Field: 0: trapezoidal, 1: Simpson, 2: Root (only for analytical formula case)
  unsigned int mNumericalIntegrationSteps = 1;             ///< number of intervalls during numerical integration are taken.
  int mStepWidth = 1;                                      ///< during the calculation of the corrections/distortions it is assumed that the electron drifts on a line from deltaZ = z0 -> z1. The value sets the deltaZ width: 1: deltaZ=zBin/1, 5: deltaZ=zBin/5

  DataT fC0 = 0; ///< coefficient C0 (compare Jim Thomas's notes for definitions)
  DataT fC1 = 0; ///< coefficient C1 (compare Jim Thomas's notes for definitions)

  static constexpr int FNSIDES = o2::tpc::SIDES;
  bool mIsEfieldSet[FNSIDES]{};     ///< flag if E-fields are set
  bool mIsLocalCorrSet[FNSIDES]{};  ///< flag if local corrections are set
  bool mIsLocalDistSet[FNSIDES]{};  ///< flag if local distortions are set
  bool mIsGlobalCorrSet[FNSIDES]{}; ///< flag if global corrections are set
  bool mIsGlobalDistSet[FNSIDES]{}; ///< flag if global distortions are set

  RegularGrid mGrid3D{ZMIN, RMIN, PHIMIN, GRIDSPACINGZ, GRIDSPACINGR, GRIDSPACINGPHI}; ///< grid properties

  DataContainer mLocalDistdR[FNSIDES]{};    ///< data storage for local distortions dR
  DataContainer mLocalDistdZ[FNSIDES]{};    ///< data storage for local distortions dZ
  DataContainer mLocalDistdRPhi[FNSIDES]{}; ///< data storage for local distortions dRPhi

  DataContainer mLocalCorrdR[FNSIDES]{};    ///< data storage for local corrections dR
  DataContainer mLocalCorrdZ[FNSIDES]{};    ///< data storage for local corrections dZ
  DataContainer mLocalCorrdRPhi[FNSIDES]{}; ///< data storage for local corrections dRPhi

  DataContainer mGlobalDistdR[FNSIDES]{};    ///< data storage for global distortions dR
  DataContainer mGlobalDistdZ[FNSIDES]{};    ///< data storage for global distortions dZ
  DataContainer mGlobalDistdRPhi[FNSIDES]{}; ///< data storage for global distortions dRPhi

  DataContainer mGlobalCorrdR[FNSIDES]{};    ///< data storage for global corrections dR
  DataContainer mGlobalCorrdZ[FNSIDES]{};    ///< data storage for global corrections dZ
  DataContainer mGlobalCorrdRPhi[FNSIDES]{}; ///< data storage for global corrections dRPhi

  DataContainer mDensity[FNSIDES]{};           ///< data storage for space charge density
  DataContainer mPotential[FNSIDES]{};         ///< data storage for the potential
  DataContainer mElectricFieldEr[FNSIDES]{};   ///< data storage for the electric field Er
  DataContainer mElectricFieldEz[FNSIDES]{};   ///< data storage for the electric field Ez
  DataContainer mElectricFieldEphi[FNSIDES]{}; ///< data storage for the electric field Ephi

  DataT getInvSpacingZ() const { return mGrid3D.getInvSpacingX(); }

  DataT getInvSpacingR() const { return mGrid3D.getInvSpacingY(); }

  DataT getInvSpacingPhi() const { return mGrid3D.getInvSpacingZ(); }

  std::string getSideName(const o2::tpc::Side side) const
  {
    return side == o2::tpc::Side::A ? "A" : "C";
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

  void processGlobalDistCorr(const DataT radius, const DataT phi, const DataT z0Tmp, const DataT z1Tmp, DataT& ddR, DataT& ddPhi, DataT& ddZ, const NumericalFields<DataT, Nz, Nr, Nphi>& formulaStruct) const
  {
    calcDistortions(radius, phi, z0Tmp, z1Tmp, ddR, ddPhi, ddZ, formulaStruct);
  }

  void processGlobalDistCorr(const DataT radius, const DataT phi, const DataT z0Tmp, const DataT z1Tmp, DataT& ddR, DataT& ddRPhi, DataT& ddZ, const DistCorrInterpolator<DataT, Nz, Nr, Nphi>& formulaStruct) const
  {
    ddR = formulaStruct.evaldR(z0Tmp, radius, phi);
    ddZ = formulaStruct.evaldZ(z0Tmp, radius, phi);
    ddRPhi = formulaStruct.evaldRPhi(z0Tmp, radius, phi) / radius;
  }
};

///
/// ========================================================================================================
///                                Inline implementations of some methods
/// ========================================================================================================
///

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
template <typename ElectricFields>
void O2TPCSpaceCharge3DCalc<DataT, Nz, Nr, Nphi>::integrateEFieldsRoot(const DataT p1r, const DataT p1phi, const DataT p1z, const DataT p2z, DataT& localIntErOverEz, DataT& localIntEPhiOverEz, DataT& localIntDeltaEz, ElectricFields& formulaStruct) const
{
  const DataT ezField = getEzField();
  TF1 fErOverEz("fErOverEz", [&](double* x, double* p) { (void)p; return static_cast<double>(formulaStruct.evalEr(p1r, p1phi, static_cast<DataT>(x[0])) / (formulaStruct.evalEz(p1r, p1phi, static_cast<DataT>(x[0])) + ezField)); }, p1z, p2z, 1);
  localIntErOverEz = static_cast<DataT>(fErOverEz.Integral(p1z, p2z));

  TF1 fEphiOverEz("fEPhiOverEz", [&](double* x, double* p) { (void)p; return static_cast<double>(formulaStruct.evalEphi(p1r, p1phi, static_cast<DataT>(x[0])) / (formulaStruct.evalEz(p1r, p1phi, static_cast<DataT>(x[0])) + ezField)); }, p1z, p2z, 1);
  localIntEPhiOverEz = static_cast<DataT>(fEphiOverEz.Integral(p1z, p2z));

  TF1 fEz("fEZOverEz", [&](double* x, double* p) { (void)p; return static_cast<double>(formulaStruct.evalEz(p1r, p1phi, static_cast<DataT>(x[0])) - ezField); }, p1z, p2z, 1);
  localIntDeltaEz = static_cast<DataT>(fEz.Integral(p1z, p2z));
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
template <typename ElectricFields>
void O2TPCSpaceCharge3DCalc<DataT, Nz, Nr, Nphi>::integrateEFieldsTrapezoidal(const DataT p1r, const DataT p1phi, const DataT p1z, const DataT p2z, DataT& localIntErOverEz, DataT& localIntEPhiOverEz, DataT& localIntDeltaEz, ElectricFields& formulaStruct) const
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

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
template <typename ElectricFields>
void O2TPCSpaceCharge3DCalc<DataT, Nz, Nr, Nphi>::integrateEFieldsSimpson(const DataT p1r, const DataT p1phi, const DataT p1z, const DataT p2z, DataT& localIntErOverEz, DataT& localIntEPhiOverEz, DataT& localIntDeltaEz, ElectricFields& formulaStruct) const
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

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
template <typename ElectricFields>
void O2TPCSpaceCharge3DCalc<DataT, Nz, Nr, Nphi>::integrateEFieldsSimpsonExperimental(const DataT p1r, const DataT p2r, const DataT p1phi, const DataT p2phi, const DataT p1z, const DataT p2z, DataT& localIntErOverEz, DataT& localIntEPhiOverEz, DataT& localIntDeltaEz, ElectricFields& formulaStruct) const
{
  //==========simpsons rule see: https://en.wikipedia.org/wiki/Simpson%27s_rule =============================
  const DataT ezField = getEzField();
  const DataT p2phiSave = regulatePhi(p2phi);

  const DataT fielder0 = formulaStruct.evalEr(p1z, p1r, p1phi);
  const DataT fieldez0 = formulaStruct.evalEz(p1z, p1r, p1phi);
  const DataT fieldephi0 = formulaStruct.evalEphi(p1z, p1r, p1phi);

  const DataT fielder1 = formulaStruct.evalEr(p2z, p2r, p2phiSave);
  const DataT fieldez1 = formulaStruct.evalEz(p2z, p2r, p2phiSave);
  const DataT fieldephi1 = formulaStruct.evalEphi(p2z, p2r, p2phiSave);

  const DataT eZ0 = ezField + fieldez0;
  const DataT eZ1 = ezField + fieldez1;

  DataT fieldSum1ErOverEz = 0;
  DataT fieldSum2ErOverEz = 0;
  DataT fieldSum1EphiOverEz = 0;
  DataT fieldSum2EphiOverEz = 0;
  DataT fieldSum1Ez = 0;
  DataT fieldSum2Ez = 0;

  const DataT pHalfZ = 0.5 * (p1z + p2z);                        // dont needs to be regulated since p1z and p2z are already regulated
  const DataT pHalfPhiSave = regulatePhi(0.5 * (p1phi + p2phi)); // needs to be regulated since p2phi is not regulated
  const DataT pHalfR = 0.5 * (p1r + p2r);

  const DataT ezField2 = formulaStruct.evalEz(pHalfZ, pHalfR, pHalfPhiSave);
  const DataT ezField2Denominator = 1 / (ezField + ezField2);
  fieldSum2ErOverEz += formulaStruct.evalEr(pHalfZ, pHalfR, pHalfPhiSave) * ezField2Denominator;
  fieldSum2EphiOverEz += formulaStruct.evalEphi(pHalfZ, pHalfR, pHalfPhiSave) * ezField2Denominator;
  fieldSum2Ez += ezField2;

  const DataT deltaXSimpsonSixth = (p2z - p1z) / 6;
  localIntErOverEz = deltaXSimpsonSixth * (2 * fieldSum1ErOverEz + 4 * fieldSum2ErOverEz + fielder0 / eZ0 + fielder1 / eZ1);
  localIntEPhiOverEz = deltaXSimpsonSixth * (2 * fieldSum1EphiOverEz + 4 * fieldSum2EphiOverEz + fieldephi0 / eZ0 + fieldephi1 / eZ1);
  localIntDeltaEz = deltaXSimpsonSixth * (2 * fieldSum1Ez + 4 * fieldSum2Ez + fieldez0 + fieldez1);
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
template <typename ElectricFields>
void O2TPCSpaceCharge3DCalc<DataT, Nz, Nr, Nphi>::calcLocalDistortionsCorrections(const int type, ElectricFields& formulaStruct)
{
  const o2::tpc::Side side = formulaStruct.getSide();
#pragma omp parallel for
  for (size_t iPhi = 0; iPhi < Nphi; iPhi++) {
    const DataT phi = getPhiVertex(iPhi);
    for (size_t iR = 0; iR < Nr; iR++) {
      const DataT radius = getRVertex(iR);
      for (size_t iZ = 0; iZ < Nz - 1; ++iZ) {

        // set z coordinated depending on distortions or correction calculation
        const DataT z0 = type == 1 ? getZVertex(iZ + 1) : getZVertex(iZ);
        const DataT z1 = type == 1 ? getZVertex(iZ) : getZVertex(iZ + 1);

        const int iSteps = mStepWidth;
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

          if (z0Tmp + stepSize >= ZMAX || z0Tmp + stepSize < ZMIN) { // set loop to exit if the readout or central electrode is reached
            break;
          }
        }

        if (type == 1) {
          mLocalCorrdR[side](iZ + 1, iR, iPhi) = drTmp;
          mLocalCorrdRPhi[side](iZ + 1, iR, iPhi) = dRPhiTmp;
          mLocalCorrdZ[side](iZ + 1, iR, iPhi) = dzTmp;
        } else {
          mLocalDistdR[side](iZ, iR, iPhi) = drTmp;
          mLocalDistdRPhi[side](iZ, iR, iPhi) = dRPhiTmp;
          mLocalDistdZ[side](iZ, iR, iPhi) = dzTmp;
        }
      }
    }
  }
  if (type == 1) {
    mIsLocalCorrSet[side] = true;
  } else {
    mIsLocalDistSet[side] = true;
  }
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
template <typename ElectricFields>
void O2TPCSpaceCharge3DCalc<DataT, Nz, Nr, Nphi>::calcDistortions(const DataT p1r, const DataT p1phi, const DataT p1z, const DataT p2z, DataT& ddR, DataT& ddRPhi, DataT& ddZ, ElectricFields& formulaStruct) const
{
  DataT localIntErOverEz = 0;
  DataT localIntEPhiOverEz = 0;
  DataT localIntDeltaEz = 0;

  const int iN = mNumericalIntegrationStrategy == SimpsonExperimental ? 5 : 1;
  for (int i = 0; i < iN; ++i) {
    switch (mNumericalIntegrationStrategy) {
      case Simpson:
        integrateEFieldsSimpson(p1r, p1phi, p1z, p2z, localIntErOverEz, localIntEPhiOverEz, localIntDeltaEz, formulaStruct);
        break;
      case SimpsonExperimental:
        integrateEFieldsSimpsonExperimental(p1r, regulateR(p1r + ddR), p1phi, (p1phi + ddRPhi), p1z, regulateZ(p2z + ddZ), localIntErOverEz, localIntEPhiOverEz, localIntDeltaEz, formulaStruct);
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
    ddZ = -1 * localIntDeltaEz * TPCParameters<DataT>::DVDE;
  }
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
template <typename Formulas>
void O2TPCSpaceCharge3DCalc<DataT, Nz, Nr, Nphi>::calcGlobalDistortions(const Formulas& formulaStruct)
{
  const o2::tpc::Side side = formulaStruct.getSide();
  // loop over tpc volume and let the electron drift from each vertex to the readout of the tpc
#pragma omp parallel for
  for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
    const DataT phi0 = getPhiVertex(iPhi);
    for (size_t iR = 0; iR < Nr; ++iR) {
      const DataT r0 = getRVertex(iR);
      for (size_t iZ = 0; iZ < Nz - 1; ++iZ) {

        DataT drDist = 0.0;
        DataT dPhiDist = 0.0;
        DataT dzDist = 0.0;

        const DataT z0 = getZVertex(iZ); // the electron starts at z0, r0, phi0
        const DataT stepSize = formulaStruct.ID == 2 ? GRIDSPACINGZ : GRIDSPACINGZ / mStepWidth;

        int iter = 0;
        for (;;) {
          const DataT z0Tmp = z0 + dzDist + iter * stepSize;
          const DataT z1Tmp = regulateZ(z0Tmp + stepSize);
          const DataT radius = regulateR(r0 + drDist);
          const DataT phi = regulatePhi(phi0 + dPhiDist);

          DataT ddR = 0;
          DataT ddPhi = 0;
          DataT ddZ = 0;

          // get the distortion from interpolation of local distortions or electric field
          processGlobalDistCorr(radius, phi, z0Tmp, z1Tmp, ddR, ddPhi, ddZ, formulaStruct);

          // add local distortion
          drDist += ddR;
          dPhiDist += ddPhi;
          dzDist += ddZ;

          if (z0Tmp + stepSize >= ZMAX) { // set loop to exit if the readout is reached
            break;
          }
          ++iter;
        }
        mGlobalDistdR[side](iZ, iR, iPhi) = drDist;
        mGlobalDistdRPhi[side](iZ, iR, iPhi) = dPhiDist * getRVertex(iR);
        mGlobalDistdZ[side](iZ, iR, iPhi) = dzDist;
      }
    }
  }
  mIsGlobalDistSet[side] = true;
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
template <typename Formulas>
void O2TPCSpaceCharge3DCalc<DataT, Nz, Nr, Nphi>::calcGlobalCorrections(const Formulas& formulaStruct)
{
  const o2::tpc::Side side = formulaStruct.getSide();
#pragma omp parallel for
  for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
    const DataT phi0 = getPhiVertex(iPhi);
    for (size_t iR = 0; iR < Nr; ++iR) {
      const DataT r0 = getRVertex(iR);
      DataT drCorr = 0;
      DataT dPhiCorr = 0;
      DataT dzCorr = 0;
      for (size_t iZ = Nz - 1; iZ >= 1; --iZ) {
        const int iSteps = formulaStruct.ID == 2 ? 1 : mStepWidth; // if one used local corrections no step width is nneded. since it is already used for calculation of the local corrections
        const DataT z0 = getZVertex(iZ);                           // the electron starts at z0, r0, phi0
        const DataT stepSize = -GRIDSPACINGZ / iSteps;

        for (int iter = 0; iter < iSteps; ++iter) {
          const DataT radius = regulateR(r0 + drCorr);
          const DataT phi = regulatePhi(phi0 + dPhiCorr);
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
        mGlobalCorrdR[side](iZ - 1, iR, iPhi) = drCorr;
        mGlobalCorrdRPhi[side](iZ - 1, iR, iPhi) = dPhiCorr * getRVertex(iR);
        mGlobalCorrdZ[side](iZ - 1, iR, iPhi) = dzCorr;
      }
    }
  }
  mIsGlobalCorrSet[side] = true;
}

} // namespace tpc
} // namespace o2

#endif
