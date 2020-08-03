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

// Root includes
#include "TF1.h" /// for numerical intergration only

#ifdef WITH_OPENMP
#include <omp.h>
#endif

template <typename DataT = float>
struct AnalyticalFields {

  DataT parA{1e-5}; ///< parameter [0] of functions
  DataT parB{0.5};  ///< parameter [1] of functions
  DataT parC{1e-4}; ///< parameter [2] of functions
  static constexpr unsigned int id = 0; ///< needed to distinguish between the differrent structs

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for electric field Er for given coordinate
  DataT evalEr(DataT z, DataT r, DataT phi) const
  {
    return erFunc(z, r, phi);
  }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for electric field Ez for given coordinate
  DataT evalEz(DataT z, DataT r, DataT phi) const
  {
    return ezFunc(z, r, phi);
  }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for electric field Ephi for given coordinate
  DataT evalEphi(DataT z, DataT r, DataT phi) const
  {
    return ephiFunc(z, r, phi);
  }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for the potential for given coordinate
  DataT evalPotential(DataT z, DataT r, DataT phi) const
  {
    return potentialFunc(z, r, phi);
  }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for the space charge density for given coordinate
  DataT evalDensity(DataT z, DataT r, DataT phi) const
  {
    return densityFunc(z, r, phi);
  }

  /// analytical potential
  std::function<DataT(DataT, DataT, DataT)> potentialFunc = [& parA = parA, &parB = parB, &parC = parC](DataT z, DataT r, DataT phi) {
    return -parA * (std::pow((-r + 254.5 + 83.5), 4) - 338.0 * std::pow((-r + 254.5 + 83.5), 3) + 21250.75 * std::pow((-r + 254.5 + 83.5), 2)) * std::cos(parB * phi) * std::cos(parB * phi) * std::exp(-1 * parC * (z - 125) * (z - 125));
  };

  /// analytical space charge - NOTE: if the space charge density is calculated analytical there would be a - sign in the formula (-parA)  - however since its an e- the sign is flipped (IS THIS CORRECT??? see for minus sign: AliTPCSpaceCharge3DCalc::SetPotentialBoundaryAndChargeFormula)-
  std::function<DataT(DataT, DataT, DataT)> densityFunc = [& parA = parA, &parB = parB, &parC = parC](DataT z, DataT r, DataT phi) {
    return parA * ((1 / r * 16 * (-3311250 + 90995.5 * r - 570.375 * r * r + r * r * r)) * std::cos(parB * phi) * std::cos(parB * phi) * std::exp(-1 * parC * (z - 125) * (z - 125)) +
                   (std::pow(-r + 254.5 + 83.5, 4) - 338.0 * std::pow(-r + 254.5 + 83.5, 3) + 21250.75 * std::pow(-r + 254.5 + 83.5, 2)) / (r * r) * std::exp(-1 * parC * (z - 125) * (z - 125)) * -2 * parB * parB * std::cos(2 * parB * phi) +
                   (std::pow(-r + 254.5 + 83.5, 4) - 338.0 * std::pow(-r + 254.5 + 83.5, 3) + 21250.75 * std::pow(-r + 254.5 + 83.5, 2)) * std::cos(parB * phi) * std::cos(parB * phi) * 2 * parC * std::exp(-1 * parC * (z - 125) * (z - 125)) * (2 * parC * (z - 125) * (z - 125) - 1));
  };

  /// analytical electric field Er
  std::function<DataT(DataT, DataT, DataT)> erFunc = [& parA = parA, &parB = parB, &parC = parC](DataT z, DataT r, DataT phi) {
    return parA * 4 * (r * r * r - 760.5 * r * r + 181991 * r - 1.3245 * std::pow(10, 7)) * std::cos(parB * phi) * std::cos(parB * phi) * std::exp(-1 * parC * (z - 125) * (z - 125));
  };

  /// analytical electric field Ephi
  std::function<DataT(DataT, DataT, DataT)> ephiFunc = [& parA = parA, &parB = parB, &parC = parC](DataT z, DataT r, DataT phi) {
    return parA * (std::pow(-r + 254.5 + 83.5, 4) - 338.0 * std::pow(-r + 254.5 + 83.5, 3) + 21250.75 * (-r + 254.5 + 83.5) * (-r + 254.5 + 83.5)) / r * std::exp(-1 * parC * (z - 125) * (z - 125)) * -parB * std::sin(2 * parB * phi);
  };

  /// analytical electric field Ez
  std::function<DataT(DataT, DataT, DataT)> ezFunc = [& parA = parA, &parB = parB, &parC = parC](DataT z, DataT r, DataT phi) {
    return parA * (std::pow(-r + 254.5 + 83.5, 4) - 338.0 * std::pow(-r + 254.5 + 83.5, 3) + 21250.75 * (-r + 254.5 + 83.5) * (-r + 254.5 + 83.5)) * std::cos(parB * phi) * std::cos(parB * phi) * -2 * parC * (z - 125) * std::exp(-1 * parC * (z - 125) * (z - 125));
  };
};

template <typename DataT = float, typename Grid3D = RegularGrid3D<>>
struct NumericalFields {
  // using RegularGrid = RegularGrid3D<DataT, Nz, Nr, Nphi>;

  NumericalFields(const Grid3D& gridErTmp, const Grid3D& gridEzTmp, const Grid3D& gridEphiTmp) : gridEr{gridErTmp}, gridEz{gridEzTmp}, gridEphi{gridEphiTmp} {};

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for electric field Er for given coordinate
  DataT evalEr(DataT z, DataT r, DataT phi) const
  {
    return interpolatorEr(z, r, phi);
  }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for electric field Ez for given coordinate
  DataT evalEz(DataT z, DataT r, DataT phi) const
  {
    return interpolatorEz(z, r, phi);
  }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for electric field Ephi for given coordinate
  DataT evalEphi(DataT z, DataT r, DataT phi) const
  {
    return interpolatorEphi(z, r, phi);
  }

  const Grid3D& gridEr{};   // adress to the data container of the grid
  const Grid3D& gridEz{};   // adress to the data container of the grid
  const Grid3D& gridEphi{}; // adress to the data container of the grid
  const bool circularZ = false;
  const bool circularR = false;
  const bool circularPhi = true;
  TriCubicInterpolator<DataT, Grid3D> interpolatorEr{gridEr, circularZ, circularR, circularPhi};
  TriCubicInterpolator<DataT, Grid3D> interpolatorEz{gridEz, circularZ, circularR, circularPhi};
  TriCubicInterpolator<DataT, Grid3D> interpolatorEphi{gridEphi, circularZ, circularR, circularPhi};
  static constexpr unsigned int id = 1; ///< needed to distinguish between the differrent structs
};

template <typename DataT = float, typename Grid3D = RegularGrid3D<>>
struct LocalDistCorr {

  LocalDistCorr(const Grid3D& lDistCorrdR, const Grid3D& lDistCorrdZ, const Grid3D& lDistCorrdRPhi) : mDistCorrdR{lDistCorrdR}, mDistCorrdZ{lDistCorrdZ}, mDistCorrdRPhi{lDistCorrdRPhi} {};

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for the local distortion or correction dR for given coordinate
  DataT evaldR(DataT z, DataT r, DataT phi) const
  {
    return mInterpolatorDistCorrdR(z, r, phi);
  }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for the local distortion or correction dZ for given coordinate
  DataT evaldZ(DataT z, DataT r, DataT phi) const
  {
    return mInterpolatorDistCorrdZ(z, r, phi);
  }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for the local distortion or correction dRPhi for given coordinate
  DataT evaldRPhi(DataT z, DataT r, DataT phi) const
  {
    return mInterpolatorDistCorrdRPhi(z, r, phi);
  }

  const Grid3D& mDistCorrdR{};    // adress to the data container of the grid
  const Grid3D& mDistCorrdZ{};    // adress to the data container of the grid
  const Grid3D& mDistCorrdRPhi{}; // adress to the data container of the grid
  const bool circularZ = false;
  const bool circularR = false;
  const bool circularPhi = true;
  TriCubicInterpolator<DataT, Grid3D> mInterpolatorDistCorrdR{mDistCorrdR, circularZ, circularR, circularPhi};
  TriCubicInterpolator<DataT, Grid3D> mInterpolatorDistCorrdZ{mDistCorrdZ, circularZ, circularR, circularPhi};
  TriCubicInterpolator<DataT, Grid3D> mInterpolatorDistCorrdRPhi{mDistCorrdRPhi, circularZ, circularR, circularPhi};
  static constexpr unsigned int id = 2; ///< needed to distinguish between the differrent structs
};

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

  static constexpr DataT getGridSpacingR() { return mGridSpacingR; }
  static constexpr DataT getGridSpacingZ() { return mGridSpacingZ; }
  static constexpr DataT getGridSpacingPhi() { return mGridSpacingPhi; }
  static constexpr DataT getEzField() { return (ASolv::fgkCathodeV - ASolv::fgkGG) / ASolv::fgkTPCZ0; }
  // constexpr DataT getRMin() const { return ASolv::fgkIFCRadius; }
  const RegularGrid& getGrid3D() const { return mGrid3D; }

  NumericalFields<DataT, RegularGrid> getNumericalFieldsInterpolator() const
  {
    NumericalFields<DataT, RegularGrid> numFields(mElectricFieldEr, mElectricFieldEz, mElectricFieldEphi);
    return numFields;
  }

  LocalDistCorr<DataT, RegularGrid> getLocalDistInterpolator() const
  {
    LocalDistCorr<DataT, RegularGrid> numFields(mLocalDistdR, mLocalDistdZ, mLocalDistdRPhi);
    return numFields;
  }

  LocalDistCorr<DataT, RegularGrid> getLocalCorrInterpolator() const
  {
    LocalDistCorr<DataT, RegularGrid> numFields(mLocalCorrdR, mLocalCorrdZ, mLocalCorrdRPhi);
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
    const int er = mGlobalDistdR.storeValuesToFile(outf, "distR");
    const int ez = mGlobalDistdZ.storeValuesToFile(outf, "distZ");
    const int ephi = mGlobalDistdRPhi.storeValuesToFile(outf, "distRPhi");
    return er + ez + ephi;
  }

  void setGlobalDistortionsFromFile(TFile& inpf)
  {
    mGlobalDistdR.initFromFile(inpf, "distR");
    mGlobalDistdZ.initFromFile(inpf, "distZ");
    mGlobalDistdRPhi.initFromFile(inpf, "distRPhi");
  }

  void printPotential() const
  {
    std::cout << mPotential;
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

  DataT regulateZ(const DataT pos)
  {
    return mGrid3D.clampToGrid(pos, 0);
  }

  DataT regulateR(const DataT pos)
  {
    return mGrid3D.clampToGrid(pos, 1);
  }

  // TODO make this circular
  DataT regulatePhi(const DataT pos)
  {
    return mGrid3D.clampToGrid(pos, 2);
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

  void processGlobalDist(const DataT radius, const DataT phi, const DataT z0Tmp, const DataT z1Tmp, DataT& ddR, DataT& ddPhi, DataT& ddZ, const AnalyticalFields<DataT>& formulaStruct) const{
    calcDistortions(radius, phi, z0Tmp, z1Tmp, ddR, ddPhi, ddZ, formulaStruct);
  }

  void processGlobalDist(const DataT radius, const DataT phi, const DataT z0Tmp, const DataT z1Tmp, DataT& ddR, DataT& ddPhi, DataT& ddZ, const NumericalFields<DataT, RegularGrid>& formulaStruct) const{
    calcDistortions(radius, phi, z0Tmp, z1Tmp, ddR, ddPhi, ddZ, formulaStruct);
  }

  void processGlobalDist(const DataT radius, const DataT phi, const DataT z0Tmp, const DataT z1Tmp, DataT& ddR, DataT& ddRPhi, DataT& ddZ, const LocalDistCorr<DataT, RegularGrid>& formulaStruct) const{
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
          const DataT z0Tmp = z0 + iter * stepSize + dzTmp;
          const DataT z1Tmp = z0Tmp + stepSize;

          DataT ddR = 0;
          DataT ddPhi = 0;
          DataT ddZ = 0;
          const DataT radiusTmp = radius + drTmp;
          const DataT phiTmp = phi + dPhiTmp;

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
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
template <typename ElectricFields>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::calcDistortions(const DataT p1r, const DataT p1phi, const DataT p1z, const DataT p2z, DataT& ddR, DataT& ddRPhi, DataT& ddZ, ElectricFields& formulaStruct) const
{
  DataT localIntErOverEz = 0;
  DataT localIntEPhiOverEz = 0;
  DataT localIntDeltaEz = 0;

  const int iN = mNumericalIntegrationStrategy == SimpsonExperimental ? 4 : 1;
  for (int i = 0; i < iN; ++i) {
    switch (mNumericalIntegrationStrategy) {
      case Simpson:
        integrateEFieldsSimpson(p1r, p1phi, p1z, p2z, localIntErOverEz, localIntEPhiOverEz, localIntDeltaEz, formulaStruct);
        break;
      case SimpsonExperimental:
        integrateEFieldsSimpsonExperimental(p1r, p1r + ddR, p1phi, p1phi + ddRPhi, p1z, p2z + ddZ, localIntErOverEz, localIntEPhiOverEz, localIntDeltaEz, formulaStruct);
        break;
      case Trapezoidal:;
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
    std::cout << "iPhi: " << iPhi << std::endl;
    for (size_t iR = 0; iR < Nr; ++iR) {
      for (size_t iZ = 0; iZ < Nz - 1; ++iZ) {

        DataT drDist = 0.0;
        DataT dPhi = 0.0;
        DataT dzDist = 0.0;

        bool readoutNotreached = true;
        int iter = 0;
        while (readoutNotreached) {
          const DataT z0 = getZVertex(iZ); // the electron starts at phi, radius, z0
          const DataT z1 = getZVertex(iZ + 1);
          const DataT stepSize = formulaStruct.id == 2 ? (z1 - z0) : (z1 - z0) / getIntegrationSteps();

          const DataT z0Tmp = z0 + dzDist + iter * stepSize;

          if (z0Tmp + stepSize >= ASolv::fgkTPCZ0) { // set loop to exit if teh readout is reached
            readoutNotreached = false;
          }

          const DataT z1Tmp = regulateZ(z0Tmp + stepSize);
          const DataT radius = regulateR(getRVertex(iR) + drDist);
          const DataT phi = regulatePhi(getPhiVertex(iPhi) + dPhi); // TODO phi doenst have to be regulated (circular)

          DataT ddR = 0;
          DataT ddPhi = 0;
          DataT ddZ = 0;

          // get the distortion from interpolation of local distortions or electric field
          processGlobalDist(radius, phi, z0Tmp, z1Tmp, ddR, ddPhi, ddZ, formulaStruct);

          // add local distortion
          drDist += ddR;
          dPhi += ddPhi;
          dzDist += ddZ;

          ++iter;
        }
        mGlobalDistdR(iZ, iR, iPhi) = drDist;
        mGlobalDistdRPhi(iZ, iR, iPhi) = dPhi * getRVertex(iR);
        mGlobalDistdZ(iZ, iR, iPhi) = dzDist;
      }
    }
  }
}
#endif
