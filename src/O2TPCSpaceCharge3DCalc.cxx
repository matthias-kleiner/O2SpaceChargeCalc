#include "include/O2TPCSpaceCharge3DCalc.h"

#include "TF1.h" /// for numerical intergration only

templateClassImp(O2TPCSpaceCharge3DCalc);
// const int nTHREADS = 8;

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::getDistortionsAnalytical(const DataT p1r, const DataT p1phi, const DataT p1z, const DataT p2z, DataT& ddR, DataT& ddRPhi, DataT& ddZ, AnalyticalFields<DataT>& formulaStruct) const
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

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::fillBoundaryAndChargeDensities(AnalyticalFields<DataT>& formulaStruct, const int maxIteration, const DataT stoppingConvergence)
{
  RegularGrid3D<DataT, Nz, Nr, Nphi> gridPotentialTmp{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};
  RegularGrid3D<DataT, Nz, Nr, Nphi> gridDensityTmp{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};

  for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
    const DataT phi = getPhiVertex(iPhi);
    for (size_t iR = 0; iR < Nr; ++iR) {
      const DataT radius = getRVertex(iR);
      for (size_t iZ = 0; iZ < Nz; ++iZ) {
        const DataT z = getZVertex(iZ);
        gridDensityTmp(iZ, iR, iPhi) = formulaStruct.evalDensity(z, radius, phi);

        if ((iR == 0) || (iR == (Nr - 1)) || (iZ == 0) || (iZ == (Nz - 1))) {
          gridPotentialTmp(iZ, iR, iPhi) = formulaStruct.evalPotential(z, radius, phi);
        }
      }
    }
  }

  // TODO MODIFY AliTPCPoissonSolver class to accept grid instead TMATRIXD
  TMatrixD* matricesPotential[Nphi];
  TMatrixD* matricesDensity[Nphi];
  for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
    matricesPotential[iPhi] = new TMatrixD(Nr, Nz);
    matricesDensity[iPhi] = new TMatrixD(Nr, Nz);
    for (size_t iR = 0; iR < Nr; ++iR) {
      for (size_t iZ = 0; iZ < Nz; ++iZ) {
        (*matricesPotential[iPhi])(iR, iZ) = gridPotentialTmp(iZ, iR, iPhi);
        (*matricesDensity[iPhi])(iR, iZ) = gridDensityTmp(iZ, iR, iPhi);
      }
    }
  }

  ASolv::fgConvergenceError = stoppingConvergence;
  ASolv poissonSolver;
  const int symmetry = 0;
  poissonSolver.PoissonSolver3D(matricesPotential, matricesDensity, Nr, Nz, Nphi, maxIteration, symmetry);

  //convert potential back to regular grid
  for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
    for (size_t iR = 0; iR < Nr; ++iR) {
      for (size_t iZ = 0; iZ < Nz; ++iZ) {
        mPotential(iZ, iR, iPhi) = static_cast<DataT>((*matricesPotential[iPhi])(iR, iZ));
      }
    }
  }

  for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
    delete matricesPotential[iPhi];
    delete matricesDensity[iPhi];
  }
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::calcEField()
{
  const int symmetry = 0;

  for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
    size_t tmpPlus = iPhi + 1;
    int signPlus = 1;
    int tmpMinus = static_cast<int>(iPhi - 1);
    int signMinus = 1;
    if (symmetry == 1 || symmetry == -1) { // Reflection symmetry in phi (e.g. symmetry at sector boundaries, or half sectors, etc.)
      if (tmpPlus > Nphi - 1) {
        if (symmetry == -1) {
          signPlus = -1;
        }
        tmpPlus = Nphi - 2;
      }
      if (tmpMinus < 0) {
        tmpMinus = 1; // SHOULD IT BE =0?
        if (symmetry == -1) {
          signMinus = -1;
        }
      }
    } else { // No Symmetries in phi, no boundaries, the calculations is continuous across all phi
      if (tmpPlus > Nphi - 1) {
        tmpPlus = iPhi + 1 - Nphi;
      }
      if (tmpMinus < 0) {
        tmpMinus = static_cast<int>(iPhi - 1 + Nphi);
      }
    }

    const size_t tmpMinusS = static_cast<size_t>(tmpMinus);

    // for non-boundary V
    for (size_t iR = 1; iR < Nr - 1; iR++) {
      const DataT radius = getRVertex(iR);
      for (size_t iZ = 1; iZ < Nz - 1; iZ++) {
        mElectricFieldEr(iZ, iR, iPhi) = -1 * (mPotential(iZ, iR + 1, iPhi) - mPotential(iZ, iR - 1, iPhi)) * static_cast<DataT>(0.5) * getInvSpacingR();                                     // r direction
        mElectricFieldEz(iZ, iR, iPhi) = -1 * (mPotential(iZ + 1, iR, iPhi) - mPotential(iZ - 1, iR, iPhi)) * static_cast<DataT>(0.5) * getInvSpacingZ();                                     // z direction
        mElectricFieldEphi(iZ, iR, iPhi) = -1 * (signPlus * mPotential(iZ, iR, tmpPlus) - signMinus * mPotential(iZ, iR, tmpMinusS)) * static_cast<DataT>(0.5) * getInvSpacingPhi() / radius; // phi direction
      }
    }

    // for boundary-r
    for (size_t iZ = 0; iZ < Nz; iZ++) {
      mElectricFieldEr(iZ, 0, iPhi) = -1 * (-static_cast<DataT>(0.5) * mPotential(iZ, 2, iPhi) + 2 * mPotential(iZ, 1, iPhi) - static_cast<DataT>(1.5) * mPotential(iZ, 0, iPhi)) * getInvSpacingR();                    // forward difference
      mElectricFieldEr(iZ, Nr - 1, iPhi) = -1 * (static_cast<DataT>(1.5) * mPotential(iZ, Nr - 1, iPhi) - 2 * mPotential(iZ, Nr - 2, iPhi) + static_cast<DataT>(0.5) * mPotential(iZ, Nr - 3, iPhi)) * getInvSpacingR(); // backward difference
    }

    for (size_t iR = 0; iR < Nr; iR += Nr - 1) {
      const DataT radius = getRVertex(iR);
      for (size_t iZ = 1; iZ < Nz - 1; iZ++) {
        mElectricFieldEz(iZ, iR, iPhi) = -1 * (mPotential(iZ + 1, iR, iPhi) - mPotential(iZ - 1, iR, iPhi)) * static_cast<DataT>(0.5) * getInvSpacingZ();                                     // z direction
        mElectricFieldEphi(iZ, iR, iPhi) = -1 * (signPlus * mPotential(iZ, iR, tmpPlus) - signMinus * mPotential(iZ, iR, tmpMinusS)) * static_cast<DataT>(0.5) * getInvSpacingPhi() / radius; // phi direction
      }
    }

    // for boundary-z
    for (size_t iR = 0; iR < Nr; ++iR) {
      mElectricFieldEz(0, iR, iPhi) = -1 * (-static_cast<DataT>(0.5) * mPotential(2, iR, iPhi) + 2 * mPotential(1, iR, iPhi) - static_cast<DataT>(1.5) * mPotential(0, iR, iPhi)) * getInvSpacingZ();
      mElectricFieldEz(Nz - 1, iR, iPhi) = -1 * (static_cast<DataT>(1.5) * mPotential(Nz - 1, iR, iPhi) - 2 * mPotential(Nz - 2, iR, iPhi) + static_cast<DataT>(0.5) * mPotential(Nz - 3, iR, iPhi)) * getInvSpacingZ();
    }

    for (size_t iR = 0; iR < Nr - 1; ++iR) {
      const DataT radius = getRVertex(iR);
      for (size_t iZ = 0; iZ < Nz; iZ += Nz - 1) {
        mElectricFieldEr(iZ, iR, iPhi) = -1 * (mPotential(iZ, iR + 1, iPhi) - mPotential(iZ, iR - 1, iPhi)) * static_cast<DataT>(0.5) * getInvSpacingR();                                     // r direction
        mElectricFieldEphi(iZ, iR, iPhi) = -1 * (signPlus * mPotential(iZ, iR, tmpPlus) - signMinus * mPotential(iZ, iR, tmpMinusS)) * static_cast<DataT>(0.5) * getInvSpacingPhi() / radius; // phi direction
      }
    }

    // corner points for EPhi
    for (size_t iR = 0; iR < Nr; iR += Nr - 1) {
      const DataT radius = getRVertex(iR);
      for (size_t iZ = 0; iZ < Nz; iZ += Nz - 1) {
        mElectricFieldEphi(iZ, iR, iPhi) = -1 * (signPlus * mPotential(iZ, iR, tmpPlus) - signMinus * mPotential(iZ, iR, tmpMinusS)) * static_cast<DataT>(0.5) * getInvSpacingPhi() / radius; // phi direction
      }
    }
  }
}

template class O2TPCSpaceCharge3DCalc<float, 17, 17, 90>;
template class O2TPCSpaceCharge3DCalc<double, 17, 17, 90>;
