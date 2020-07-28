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
#include "TFormula.h"

#ifdef WITH_OPENMP
#include <omp.h>
#endif

template <typename DataT = float>
struct Formulas {

  DataT parA{1e-5}; ///< parameter [0] of functions
  DataT parB{0.5}; ///< parameter [1] of functions
  DataT parC{1e-4}; ///< parameter [2] of functions

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for electric field Er for given coordinate
  DataT evalEr(DataT r, DataT phi, DataT z) const
  {
    return erFunc(r, phi, z);
  }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for electric field Ez for given coordinate
  DataT evalEz(DataT r, DataT phi, DataT z) const
  {
    return ezFunc(r, phi, z);
  }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for electric field Ephi for given coordinate
  DataT evalEphi(DataT r, DataT phi, DataT z) const
  {
    return ephiFunc(r, phi, z);
  }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for the potential for given coordinate
  DataT evalPotential(DataT r, DataT phi, DataT z) const
  {
    return potentialFunc(r, phi, z);
  }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for the space charge density for given coordinate
  DataT evalDensity(DataT r, DataT phi, DataT z) const
  {
    return densityFunc(r, phi, z);
  }

  /// analytical potential
  std::function<DataT(DataT, DataT, DataT)> potentialFunc = [& parA = parA, &parB = parB, &parC = parC](DataT r, DataT phi, DataT z) {
    return -parA * (std::pow((-r + 254.5 + 83.5), 4) - 338.0 * std::pow((-r + 254.5 + 83.5), 3) + 21250.75 * std::pow((-r + 254.5 + 83.5), 2)) * std::cos(parB * phi) * std::cos(parB * phi) * std::exp(-1 * parC * (z - 125) * (z - 125));
  };

  /// analytical space charge density
  std::function<DataT(DataT, DataT, DataT)> densityFunc = [& parA = parA, &parB = parB, &parC = parC](DataT r, DataT phi, DataT z) {
    return -parA * ((1 / r * 16 * (-3311250 + 90995.5 * r - 570.375 * r * r + r * r * r)) * std::cos(parB * phi) * std::cos(parB * phi) * std::exp(-1 * parC * (z - 125) * (z - 125)) +
                    (std::pow(-r + 254.5 + 83.5, 4) - 338.0 * std::pow(-r + 254.5 + 83.5, 3) + 21250.75 * std::pow(-r + 254.5 + 83.5, 2)) / (r * r) * std::exp(-1 * parC * (z - 125) * (z - 125)) * -2 * parB * parB * std::cos(2 * parB * phi) +
                    (std::pow(-r + 254.5 + 83.5, 4) - 338.0 * std::pow(-r + 254.5 + 83.5, 3) + 21250.75 * std::pow(-r + 254.5 + 83.5, 2)) * std::cos(parB * phi) * std::cos(parB * phi) * 2 * parC * std::exp(-1 * parC * (z - 125) * (z - 125)) * (2 * parC * (z - 125) * (z - 125) - 1));
  };

  /// analytical electric field Er
  std::function<DataT(DataT, DataT, DataT)> erFunc = [& parA = parA, &parB = parB, &parC = parC](DataT r, DataT phi, DataT z) {
    return parA * 4 * (r * r * r - 760.5 * r * r + 181991 * r - 1.3245 * std::pow(10, 7)) * std::cos(parB * phi) * std::cos(parB * phi) * std::exp(-1 * parC * (z - 125) * (z - 125));
  };

  /// analytical electric field Ephi
  std::function<DataT(DataT, DataT, DataT)> ephiFunc = [& parA = parA, &parB = parB, &parC = parC](DataT r, DataT phi, DataT z) {
    return parA * (std::pow(-r + 254.5 + 83.5, 4) - 338.0 * std::pow(-r + 254.5 + 83.5, 3) + 21250.75 * (-r + 254.5 + 83.5) * (-r + 254.5 + 83.5)) / r * std::exp(-1 * parC * (z - 125) * (z - 125)) * -parB * std::sin(2 * parB * phi);
  };

  /// analytical electric field Ez
  std::function<DataT(DataT, DataT, DataT)> ezFunc = [& parA = parA, &parB = parB, &parC = parC](DataT r, DataT phi, DataT z) {
    return parA * (std::pow(-r + 254.5 + 83.5, 4) - 338.0 * std::pow(-r + 254.5 + 83.5, 3) + 21250.75 * (-r + 254.5 + 83.5) * (-r + 254.5 + 83.5)) * std::cos(parB * phi) * std::cos(parB * phi) * -2 * parC * (z - 125) * std::exp(-1 * parC * (z - 125) * (z - 125));
  };
};

/// \tparam DataT the type of data which is used during the calculations
/// \tparam Nr number of vertices in r direction
/// \tparam Nz number of vertices in z direction
/// \tparam Nphi number of vertices in phi direction
template <typename DataT = float, size_t Nr = 17, size_t Nz = 17, size_t Nphi = 90>
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
  // template<typename structFields>
  // void calcLocalDistortionsCorrections(structFields formulaStruct);

  /// type=0 -> distortions, type=1->corrections
  template<typename ElectricFields = Formulas<DataT>>
  void calcLocalDistortionsCorrections(const bool lcorrections, ElectricFields& formulaStruct);

  /// calculate distortions or corrections analytical with given funcion
  void getDistortionsAnalytical(const DataT p1r, const DataT p1phi, const DataT p1z, const DataT p2z, DataT& ddR, DataT& ddRPhi, DataT& ddZ, Formulas<DataT> formulaStruct) const;

  //step 4:
  // Info("AliTPCSpaceCharge3DCalc::InitSpaceCharge3DPoissonIntegralDz", "%s", Form("Step 4: Global correction/distortion cpu time: %f\n", w.CpuTime()));
  void calcGlobalDistortionsCorrections();

  // constexpr DataT getGridSpacingR() const { return mGridSpacingR; }
  // constexpr DataT getGridSpacingZ() const { return mGridSpacingZ; }
  // constexpr DataT getGridSpacingPhi() const { return mGridSpacingPhi; }
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

  int getIntegrationSteps() const { return mIntegrationSteps; }

  void setOmegaTauT1T2(DataT omegaTau, DataT t1, DataT t2)
  {
    const DataT wt0 = t2 * omegaTau;
    fC0 = 1. / (1. + wt0 * wt0);
    const DataT wt1 = t1 * omegaTau;
    fC1 = wt1 / (1. + wt1 * wt1);
  };

  void setC0C1(DataT c0, DataT c1)
  {
    fC0 = c0;
    fC1 = c1;
  }

  void setIntegrationSteps(const int nSteps) { mIntegrationSteps = nSteps; }

 private:
  using ASolv = AliTPCPoissonSolver<DataT>;

  static constexpr DataT mGridSpacingR = (ASolv::fgkOFCRadius - ASolv::fgkIFCRadius) / (Nr - 1); ///< grid spacing in r direction
  static constexpr DataT mGridSpacingZ = ASolv::fgkTPCZ0 / (Nz - 1);                             ///< grid spacing in z direction
  static constexpr DataT mGridSpacingPhi = 2 * M_PI / Nphi;                                      ///< grid spacing in phi direction // TODO CHANGE TO O2
  static constexpr DataT mRMin = ASolv::fgkIFCRadius;                                            ///< min radius
  static constexpr DataT mZMin = 0;                                                              ///< min z coordinate
  static constexpr DataT mPhiMin = 0;                                                            ///< min phi coordinate
  int mNumericalIntegrationStrategy = 1;                                                         ///< numerical integration strategy of integration of the E-Field: 0: trapezoidal, 1: Simpson, 2: Root (only for analytical formula case)
  int mIntegrationSteps = 1;                                                                     ///< number of integration steps performed between each bin in Z. e.g.: 1: direct integration from z[i] -> z[i+1]    2: z[i] -> z[i] + (z[i+1] - z[i])/2 > z[i+1]

  DataT fC0 = 0.f; ///< coefficient C0 (compare Jim Thomas's notes for definitions)
  DataT fC1 = 0.f; ///< coefficient C1 (compare Jim Thomas's notes for definitions)

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

  RegularGrid3D<DataT, Nz, Nr, Nphi> mGrid3D{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi}; ///< this grid contains the values for the local distortions/corrections, electric field etc.

  RegularGrid3D<DataT, Nz, Nr, Nphi> mLocalDistdR{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};
  RegularGrid3D<DataT, Nz, Nr, Nphi> mLocalDistdZ{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};
  RegularGrid3D<DataT, Nz, Nr, Nphi> mLocalDistdRPhi{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};

  RegularGrid3D<DataT, Nz, Nr, Nphi> mLocalCorrdR{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};
  RegularGrid3D<DataT, Nz, Nr, Nphi> mLocalCorrdZ{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};
  RegularGrid3D<DataT, Nz, Nr, Nphi> mLocalCorrdRPhi{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};
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


#endif
