#include "include/O2TPCSpaceCharge3DCalc.h"

templateClassImp(O2TPCSpaceCharge3DCalc);
// const int nTHREADS = 8;

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::fillBoundaryAndChargeDensities(AnalyticalFields<DataT>& formulaStruct, const int maxIteration, const DataT stoppingConvergence)
{
  RegularGrid gridPotentialTmp{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};
  RegularGrid gridDensityTmp{mZMin, mRMin, mPhiMin, mGridSpacingZ, mGridSpacingR, mGridSpacingPhi};

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

    for (size_t iR = 1; iR < Nr - 1; ++iR) {
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
template class O2TPCSpaceCharge3DCalc<float, 129, 129, 180>;
// template class O2TPCSpaceCharge3DCalc<double, 17, 17, 90>;
