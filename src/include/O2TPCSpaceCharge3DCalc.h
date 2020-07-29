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

template <typename DataT = float, size_t Nr = 17, size_t Nz = 17, size_t Nphi = 90>
struct NumericalFields {
  using RegularGrid = RegularGrid3D<DataT, Nz, Nr, Nphi>;

  NumericalFields(const RegularGrid& gridEr, const RegularGrid& gridEz, const RegularGrid& gridEphi) : mGridEr{gridEr}, mGridEz{gridEz}, mGridEphi{gridEphi} {};

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for electric field Er for given coordinate
  DataT evalEr(DataT z, DataT r, DataT phi) const
  {
    return mInterpolatorEr(z, r, phi);
  }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for electric field Ez for given coordinate
  DataT evalEz(DataT z, DataT r, DataT phi) const
  {
    return mInterpolatorEz(z, r, phi);
  }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for electric field Ephi for given coordinate
  DataT evalEphi(DataT z, DataT r, DataT phi) const
  {
    return mInterpolatorEphi(z, r, phi);
  }

  const RegularGrid& mGridEr{};   // adress to the data container of the grid
  const RegularGrid& mGridEz{};   // adress to the data container of the grid
  const RegularGrid& mGridEphi{}; // adress to the data container of the grid
  const bool mCircularZ = false;
  const bool mCircularR = false;
  const bool mCircularPhi = true;
  TriCubicInterpolator<DataT, RegularGrid> mInterpolatorEr{mGridEr, mCircularZ, mCircularR, mCircularPhi};
  TriCubicInterpolator<DataT, RegularGrid> mInterpolatorEz{mGridEz, mCircularZ, mCircularR, mCircularPhi};
  TriCubicInterpolator<DataT, RegularGrid> mInterpolatorEphi{mGridEphi, mCircularZ, mCircularR, mCircularPhi};
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
  // Info("AliTPCSpaceCharge3DCalc::InitSpaceCharge3DPoissonIntegralDz", "%s", Form("Step = 0: Fill Boundary and Charge Densities"));
  void fillBoundaryAndChargeDensities(AnalyticalFields<DataT>& formulaStruct, const int maxIteration, const DataT stoppingConvergence);

  // stepp 1:
  //Info("AliTPCSpaceCharge3DCalc::InitSpaceCharge3DPoissonIntegralDz", "%s", Form("Step 1: Poisson solver: %f\n", w.CpuTime()));
  void poissonSolver();

  // stepp 2:
  // Info("AliTPCSpaceCharge3DCalc::InitSpaceCharge3DPoissonIntegralDz", "%s", Form("Step 2: Electric Field Calculation: %f\n", w.CpuTime()));
  void calcEField();

  // stepp 3:
  // Info("AliTPCSpaceCharge3DCalc::InitSpaceCharge3DPoissonIntegralDz", "%s", Form("Step 3: Local distortion and correction cpu time: %f\n", w.CpuTime()));
  /// lcorrections=false -> distortions, lcorrections=true->corrections
  template <typename ElectricFields = AnalyticalFields<DataT>>
  void calcLocalDistortionsCorrections(const bool lcorrections, ElectricFields& formulaStruct);

  /// calculate distortions or corrections analytical with given funcion
  template <typename ElectricFields = AnalyticalFields<DataT>>
  void getDistortionsAnalytical(const DataT p1r, const DataT p1phi, const DataT p1z, const DataT p2z, DataT& ddR, DataT& ddRPhi, DataT& ddZ, ElectricFields& formulaStruct) const;

  //step 4:
  // Info("AliTPCSpaceCharge3DCalc::InitSpaceCharge3DPoissonIntegralDz", "%s", Form("Step 4: Global correction/distortion cpu time: %f\n", w.CpuTime()));
  void calcGlobalDistortionsCorrections();

  constexpr DataT getGridSpacingR() const { return mGridSpacingR; }
  constexpr DataT getGridSpacingZ() const { return mGridSpacingZ; }
  constexpr DataT getGridSpacingPhi() const { return mGridSpacingPhi; }
  constexpr DataT getEzField() const { return (ASolv::fgkCathodeV - ASolv::fgkGG) / ASolv::fgkTPCZ0; }
  // constexpr DataT getRMin() const { return ASolv::fgkIFCRadius; }
  const RegularGrid3D<DataT, Nr, Nz, Nphi>& getGrid3D() const { return mGrid3D; }

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

  int getIntegrationSteps() const { return mIntegrationSteps; }

  DataT getPhiVertex(size_t index) const
  {
    return mGrid3D.getZVertex(index);
  }

  DataT getRVertex(size_t index) const
  {
    return mGrid3D.getYVertex(index);
  }

  DataT getZVertex(size_t index) const
  {
    return mGrid3D.getXVertex(index);
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

  void setIntegrationSteps(const int nSteps) { mIntegrationSteps = nSteps; }

  /// numerical integration strategys
  enum IntegrationStrategy { Trapezoidal = 0,
                             Simpson = 1,
                             Root = 2 };

 private:
  using ASolv = AliTPCPoissonSolver<DataT>;

  static constexpr DataT mGridSpacingR = (ASolv::fgkOFCRadius - ASolv::fgkIFCRadius) / (Nr - 1); ///< grid spacing in r direction
  static constexpr DataT mGridSpacingZ = ASolv::fgkTPCZ0 / (Nz - 1);                             ///< grid spacing in z direction
  static constexpr DataT mGridSpacingPhi = 2 * M_PI / Nphi;                                      ///< grid spacing in phi direction // TODO CHANGE TO O2
  static constexpr DataT mRMin = ASolv::fgkIFCRadius;                                            ///< min radius
  static constexpr DataT mZMin = 0;                                                              ///< min z coordinate
  static constexpr DataT mPhiMin = 0;                                                            ///< min phi coordinate
  int mNumericalIntegrationStrategy = Simpson;                                                   ///< numerical integration strategy of integration of the E-Field: 0: trapezoidal, 1: Simpson, 2: Root (only for analytical formula case)
  int mIntegrationSteps = 1;                                                                     ///< number of integration steps performed between each bin in Z. e.g.: 1: direct integration from z[i] -> z[i+1]    2: z[i] -> z[i] + (z[i+1] - z[i])/2 > z[i+1]

  DataT fC0 = 0.f; ///< coefficient C0 (compare Jim Thomas's notes for definitions)
  DataT fC1 = 0.f; ///< coefficient C1 (compare Jim Thomas's notes for definitions)

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

  // using RegularGrid = RegularGrid3D<DataT, Nz, Nr, Nphi>;
  RegularGrid mGrid3D{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi}; ///< this grid contains the values for the local distortions/corrections, electric field etc.

  RegularGrid mLocalDistdR{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};
  RegularGrid mLocalDistdZ{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};
  RegularGrid mLocalDistdRPhi{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};

  RegularGrid mLocalCorrdR{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};
  RegularGrid mLocalCorrdZ{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};
  RegularGrid mLocalCorrdRPhi{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};

  RegularGrid mPotential{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};
  RegularGrid mElectricFieldEr{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};
  RegularGrid mElectricFieldEz{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};
  RegularGrid mElectricFieldEphi{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};
};

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
template <typename ElectricFields>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::calcLocalDistortionsCorrections(const bool lcorrections, ElectricFields& formulaStruct)
{
  // #pragma omp parallel for num_threads(nTHREADS)
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

          getDistortionsAnalytical(radiusTmp, phiTmp, z0Tmp, z1Tmp, ddR, ddPhi, ddZ, formulaStruct);
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
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::getDistortionsAnalytical(const DataT p1r, const DataT p1phi, const DataT p1z, const DataT p2z, DataT& ddR, DataT& ddRPhi, DataT& ddZ, ElectricFields& formulaStruct) const
{
  DataT localIntErOverEz = 0;
  DataT localIntEPhiOverEz = 0;
  DataT localIntDeltaEz = 0;
  const DataT ezField = getEzField();

  if (mNumericalIntegrationStrategy == Root) {
    // (void)p; -> supress warning of unused parameter
    TF1 fErOverEz("fErOverEz", [&](double* x, double* p) { (void)p; return static_cast<double>(formulaStruct.evalEr(p1r, p1phi, static_cast<DataT>(x[0])) / (formulaStruct.evalEz(p1r, p1phi, static_cast<DataT>(x[0])) + ezField)); }, p1z, p2z, 1);
    localIntErOverEz = static_cast<DataT>(fErOverEz.Integral(p1z, p2z));

    TF1 fEphiOverEz("fEPhiOverEz", [&](double* x, double* p) { (void)p; return static_cast<double>(formulaStruct.evalEphi(p1r, p1phi, static_cast<DataT>(x[0])) / (formulaStruct.evalEz(p1r, p1phi, static_cast<DataT>(x[0])) + ezField)); }, p1z, p2z, 1);
    localIntEPhiOverEz = static_cast<DataT>(fEphiOverEz.Integral(p1z, p2z));

    TF1 fEz("fEZOverEz", [&](double* x, double* p) { (void)p; return static_cast<double>(formulaStruct.evalEz(p1r, p1phi, static_cast<DataT>(x[0])) - ezField); }, p1z, p2z, 1);
    localIntDeltaEz = static_cast<DataT>(fEz.Integral(p1z, p2z));
  } else {
    const DataT fielder0 = formulaStruct.evalEr(p1z, p1r, p1phi);
    const DataT fieldez0 = formulaStruct.evalEz(p1z, p1r, p1phi);
    const DataT fieldephi0 = formulaStruct.evalEphi(p1z, p1r, p1phi);

    const DataT fielder1 = formulaStruct.evalEr(p2z, p1r, p1phi);
    const DataT fieldez1 = formulaStruct.evalEz(p2z, p1r, p1phi);
    const DataT fieldephi1 = formulaStruct.evalEphi(p2z, p1r, p1phi);

    const DataT eZ0 = ezField + fieldez0;
    const DataT eZ1 = ezField + fieldez1;

    const int nSteps = 1; //getIntegrationSteps(); //mNumericalIntegrationSteps;
    const DataT deltaX = (p2z - p1z) / nSteps;

    // trapezoidal integration
    if (mNumericalIntegrationStrategy == Trapezoidal) {
      //========trapezoidal rule==============
      DataT fieldSumEr = 0;
      DataT fieldSumEphi = 0;
      DataT fieldSumEz = 0;
      for (int i = 1; i < nSteps; ++i) {
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
    } else if (mNumericalIntegrationStrategy == Simpson) {
      //==========simpsons rule see: https://en.wikipedia.org/wiki/Simpson%27s_rule =============================
      DataT fieldSum1ErOverEz = 0;
      DataT fieldSum2ErOverEz = 0;
      DataT fieldSum1EphiOverEz = 0;
      DataT fieldSum2EphiOverEz = 0;
      DataT fieldSum1Ez = 0;
      DataT fieldSum2Ez = 0;

      for (int i = 1; i < nSteps; ++i) {
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
  }

  ddR = fC0 * localIntErOverEz + fC1 * localIntEPhiOverEz;
  ddRPhi = (fC0 * localIntEPhiOverEz - fC1 * localIntErOverEz) / p1r;
  ddZ = -1 * localIntDeltaEz * ASolv::fgkdvdE;
}

#endif
