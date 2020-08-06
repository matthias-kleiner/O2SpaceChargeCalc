#include "include/O2TPCSpaceCharge3DCalc.h"

templateClassImp(O2TPCSpaceCharge3DCalc);
// const int nTHREADS = 8;

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::performFullRun(const AnalyticalFields<DataT>& formulas, const int mode, TFile& file)
{

  using timer = std::chrono::high_resolution_clock;

  if (mode == 0) {
    std::cout << std::endl;
    auto startTotal = timer::now();

    auto start = timer::now();
    calcLocalDistortionsCorrections(0, formulas); // local distortion calculation
    auto stop = timer::now();
    std::chrono::duration<float> time = stop - start;
    std::cout << "local distortions analytical: " << time.count() << std::endl;
    dumpLocalDistortions(file);

    start = timer::now();
    calcLocalDistortionsCorrections(1, formulas); // local correction calculation
    stop = timer::now();
    time = stop - start;
    std::cout << "local corrections analytical: " << time.count() << std::endl;
    dumpLocalCorrections(file);

    start = timer::now();
    calcGlobalDistortions(formulas);
    stop = timer::now();
    time = stop - start;
    std::cout << "global distortions analytical: " << time.count() << std::endl;
    dumpGlobalDistortions(file);

    start = timer::now();
    calcGlobalCorrections(formulas);
    stop = timer::now();
    time = stop - start;
    std::cout << "global corrections analytical: " << time.count() << std::endl;
    dumpGlobalCorrections(file);

    auto stopTotal = timer::now();
    time = stopTotal - startTotal;
    std::cout << "everything is done. Total Time: " << time.count() << std::endl;
    std::cout << std::endl;
  }
  if (mode == 1) {
    std::cout << std::endl;
    auto startTotal = timer::now();

    auto start = timer::now();
    fillBoundaryAndChargeDensities(formulas);
    auto stop = timer::now();
    std::chrono::duration<float> time = stop - start;
    std::cout << "filling boundary and charge density: " << time.count() << std::endl;

    start = timer::now();
    poissonSolver(300, 1e-8);
    stop = timer::now();
    time = stop - start;
    std::cout << "poissonSolver: " << time.count() << std::endl;
    dumpPotential(file);

    start = timer::now();
    calcEField();
    stop = timer::now();
    time = stop - start;
    std::cout << "electric field calculation: " << time.count() << std::endl;
    dumpElectricFields(file);

    const auto numEFields = getElectricFieldsInterpolator();
    start = timer::now();
    calcLocalDistortionsCorrections(0, numEFields); // local distortion calculation
    stop = timer::now();
    time = stop - start;
    std::cout << "local distortions numerical: " << time.count() << std::endl;
    dumpLocalDistortions(file);

    start = timer::now();
    calcLocalDistortionsCorrections(1, numEFields); // local correction calculation
    stop = timer::now();
    time = stop - start;
    std::cout << "local corrections numerical: " << time.count() << std::endl;
    dumpLocalCorrections(file);

    start = timer::now();
    const auto lDistInterpolator = getLocalDistInterpolator();
    calcGlobalDistortions(lDistInterpolator);
    stop = timer::now();
    time = stop - start;
    std::cout << "global distortions with local distortions: " << time.count() << std::endl;
    dumpGlobalDistortions(file);

    start = timer::now();
    const auto lCorrInterpolator = getLocalCorrInterpolator();
    calcGlobalCorrections(lCorrInterpolator);
    stop = timer::now();
    time = stop - start;
    std::cout << "global corrections with local corrections: " << time.count() << std::endl;
    dumpGlobalCorrections(file);

    auto stopTotal = timer::now();
    time = stopTotal - startTotal;
    std::cout << "everything is done. Total Time: " << time.count() << std::endl;
    std::cout << std::endl;
  }
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::setFromFile(const int mode, TFile& file)
{
  setLocalDistortionsFromFile(file);
  setLocalCorrectionsFromFile(file);
  setGlobalDistortionsFromFile(file);
  setGlobalCorrectionsFromFile(file);
  if (mode == 1) {
    setPotentialFromFile(file);
    setElectricFieldsFromFile(file);
  }
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::fillBoundaryAndChargeDensities(const AnalyticalFields<DataT>& formulaStruct)
{
  for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
    const DataT phi = getPhiVertex(iPhi);
    for (size_t iR = 0; iR < Nr; ++iR) {
      const DataT radius = getRVertex(iR);
      for (size_t iZ = 0; iZ < Nz; ++iZ) {
        const DataT z = getZVertex(iZ);
        mDensity(iZ, iR, iPhi) = formulaStruct.evalDensity(z, radius, phi);

        if ((iR == 0) || (iR == (Nr - 1)) || (iZ == 0) || (iZ == (Nz - 1))) {
          mPotential(iZ, iR, iPhi) = formulaStruct.evalPotential(z, radius, phi);
        }
      }
    }
  }
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::poissonSolver(const int maxIteration, const DataT stoppingConvergence)
{
  // TODO MODIFY AliTPCPoissonSolver class to accept grid instead TMATRIXD
  TMatrixD* matricesPotential[Nphi];
  TMatrixD* matricesDensity[Nphi];
  for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
    matricesPotential[iPhi] = new TMatrixD(Nr, Nz);
    matricesDensity[iPhi] = new TMatrixD(Nr, Nz);
    for (size_t iR = 0; iR < Nr; ++iR) {
      for (size_t iZ = 0; iZ < Nz; ++iZ) {
        (*matricesPotential[iPhi])(iR, iZ) = mPotential(iZ, iR, iPhi);
        (*matricesDensity[iPhi])(iR, iZ) = mDensity(iZ, iR, iPhi);
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
  mIsEfieldSet = true;
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::calcGlobalDistWithGlobalCorrIterative(const DistCorrInterpolator<DataT, Nr, Nz, Nphi>& globCorr, const int maxIter, const DataT convZ, const DataT convR, const DataT convPhi)
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

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
NumericalFields<DataT, Nr, Nz, Nphi> O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::getElectricFieldsInterpolator() const
{
  if (!mIsEfieldSet) {
    std::cout << "============== E-Fields are not set! returning ==============" << std::endl;
  }
  NumericalFields<DataT, Nr, Nz, Nphi> numFields(mElectricFieldEr, mElectricFieldEz, mElectricFieldEphi, mGrid3D);
  return numFields;
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
DistCorrInterpolator<DataT, Nr, Nz, Nphi> O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::getLocalDistInterpolator() const
{
  if (!mIsLocalDistSet) {
    std::cout << "============== local distortions not set! returning ==============" << std::endl;
  }

  DistCorrInterpolator<DataT, Nr, Nz, Nphi> numFields(mLocalDistdR, mLocalDistdZ, mLocalDistdRPhi, mGrid3D);
  return numFields;
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
DistCorrInterpolator<DataT, Nr, Nz, Nphi> O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::getLocalCorrInterpolator() const
{
  if (!mIsLocalCorrSet) {
    std::cout << "============== local corrections not set! returning ==============" << std::endl;
  }
  DistCorrInterpolator<DataT, Nr, Nz, Nphi> numFields(mLocalCorrdR, mLocalCorrdZ, mLocalCorrdRPhi, mGrid3D);
  return numFields;
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
DistCorrInterpolator<DataT, Nr, Nz, Nphi> O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::getGlobalDistInterpolator() const
{
  if (!mIsGlobalDistSet) {
    std::cout << "============== global distortions not set! returning ==============" << std::endl;
  }
  DistCorrInterpolator<DataT, Nr, Nz, Nphi> numFields(mGlobalDistdR, mGlobalDistdZ, mGlobalDistdRPhi, mGrid3D);
  return numFields;
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
DistCorrInterpolator<DataT, Nr, Nz, Nphi> O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::getGlobalCorrInterpolator() const
{
  if (!mIsGlobalCorrSet) {
    std::cout << "============== global corrections not set! returning ==============" << std::endl;
  }
  DistCorrInterpolator<DataT, Nr, Nz, Nphi> numFields(mGlobalCorrdR, mGlobalCorrdZ, mGlobalCorrdRPhi, mGrid3D);
  return numFields;
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
int O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::dumpElectricFields(TFile& outf) const
{
  if (!mIsEfieldSet) {
    std::cout << "============== E-Fields are not set! returning ==============" << std::endl;
    return 0;
  }
  const int er = mElectricFieldEr.writeToFile(outf, "fieldEr");
  const int ez = mElectricFieldEz.writeToFile(outf, "fieldEz");
  const int ephi = mElectricFieldEphi.writeToFile(outf, "fieldEphi");
  return er + ez + ephi;
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::setElectricFieldsFromFile(TFile& inpf)
{
  mElectricFieldEr.initFromFile(inpf, "fieldEr");
  mElectricFieldEz.initFromFile(inpf, "fieldEz");
  mElectricFieldEphi.initFromFile(inpf, "fieldEphi");
  mIsEfieldSet = true;
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
int O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::dumpGlobalDistortions(TFile& outf) const
{
  if (!mIsGlobalDistSet) {
    std::cout << "============== global distortions are not set! returning ==============" << std::endl;
    return 0;
  }
  const int er = mGlobalDistdR.writeToFile(outf, "distR");
  const int ez = mGlobalDistdZ.writeToFile(outf, "distZ");
  const int ephi = mGlobalDistdRPhi.writeToFile(outf, "distRPhi");
  return er + ez + ephi;
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::setGlobalDistortionsFromFile(TFile& inpf)
{
  mIsGlobalDistSet = true;
  mGlobalDistdR.initFromFile(inpf, "distR");
  mGlobalDistdZ.initFromFile(inpf, "distZ");
  mGlobalDistdRPhi.initFromFile(inpf, "distRPhi");
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
int O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::dumpGlobalCorrections(TFile& outf) const
{
  if (!mIsGlobalCorrSet) {
    std::cout << "============== global corrections are not set! returning ==============" << std::endl;
    return 0;
  }
  const int er = mGlobalCorrdR.writeToFile(outf, "corrR");
  const int ez = mGlobalCorrdZ.writeToFile(outf, "corrZ");
  const int ephi = mGlobalCorrdRPhi.writeToFile(outf, "corrRPhi");
  return er + ez + ephi;
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::setGlobalCorrectionsFromFile(TFile& inpf)
{
  mIsGlobalCorrSet = true;
  mGlobalCorrdR.initFromFile(inpf, "corrR");
  mGlobalCorrdZ.initFromFile(inpf, "corrZ");
  mGlobalCorrdRPhi.initFromFile(inpf, "corrRPhi");
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
int O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::dumpLocalCorrections(TFile& outf) const
{
  if (!mIsLocalCorrSet) {
    std::cout << "============== local corrections are not set! returning ==============" << std::endl;
    return 0;
  }
  const int lCorrdR = mLocalCorrdR.writeToFile(outf, "lcorrR");
  const int lCorrdZ = mLocalCorrdZ.writeToFile(outf, "lcorrZ");
  const int lCorrdRPhi = mLocalCorrdRPhi.writeToFile(outf, "lcorrRPhi");
  return lCorrdR + lCorrdZ + lCorrdRPhi;
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::setLocalCorrectionsFromFile(TFile& inpf)
{
  mIsLocalCorrSet = true;
  mLocalCorrdR.initFromFile(inpf, "lcorrR");
  mLocalCorrdZ.initFromFile(inpf, "lcorrZ");
  mLocalCorrdRPhi.initFromFile(inpf, "lcorrRPhi");
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
int O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::dumpLocalDistortions(TFile& outf) const
{
  if (!mIsLocalDistSet) {
    std::cout << "============== local distortions are not set! returning ==============" << std::endl;
    return 0;
  }
  const int lCorrdR = mLocalDistdR.writeToFile(outf, "ldistR");
  const int lCorrdZ = mLocalDistdZ.writeToFile(outf, "ldistZ");
  const int lCorrdRPhi = mLocalDistdRPhi.writeToFile(outf, "ldistRPhi");
  return lCorrdR + lCorrdZ + lCorrdRPhi;
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::setLocalDistortionsFromFile(TFile& inpf)
{
  mIsLocalDistSet = true;
  mLocalDistdR.initFromFile(inpf, "ldistR");
  mLocalDistdZ.initFromFile(inpf, "ldistZ");
  mLocalDistdRPhi.initFromFile(inpf, "ldistRPhi");
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
DataT O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::regulatePhi(const DataT phi)
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

template class O2TPCSpaceCharge3DCalc<float, 17, 17, 90>;
template class O2TPCSpaceCharge3DCalc<float, 129, 129, 180>;

// using DataTF = float;
// using O2TPCSpaceCharge3DCalc17 = O2TPCSpaceCharge3DCalc<DataTF, 17, 17, 90>;
// using O2TPCSpaceCharge3DCalc129 = O2TPCSpaceCharge3DCalc<DataTF, 129, 129, 180>;
//
// template class O2TPCSpaceCharge3DCalc<DataTF, 17, 17, 90>;
// template class O2TPCSpaceCharge3DCalc<DataTF, 129, 129, 180>;
//
// using NumFields17 = NumericalFields<DataTF, 17, 17, 90>;
// using NumFields129 = NumericalFields<DataTF, 129, 129, 180>;
// using AnaFields17 = AnalyticalFields<DataTF>;
// using AnaFields129 = AnalyticalFields<DataTF>;
//
// template void O2TPCSpaceCharge3DCalc17::calcDistortions(DataTF, DataTF, DataTF, DataTF, DataTF&, DataTF&, DataTF&, NumFields17&) const;
// template void O2TPCSpaceCharge3DCalc129::calcDistortions(DataTF, DataTF, DataTF, DataTF, DataTF&, DataTF&, DataTF&, NumFields129&) const;
// template void O2TPCSpaceCharge3DCalc17::calcDistortions(DataTF, DataTF, DataTF, DataTF, DataTF&, DataTF&, DataTF&, AnaFields17&) const;
// template void O2TPCSpaceCharge3DCalc129::calcDistortions(DataTF, DataTF, DataTF, DataTF, DataTF&, DataTF&, DataTF&, AnaFields129&) const;
