#include "include/O2TPCSpaceCharge3DCalc.h"
// #include "utils/TreeStreamRedirector.h"

templateClassImp(O2TPCSpaceCharge3DCalc);

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::performFullRun(const AnalyticalFields<DataT>& formulas, const int mode, const bool electricFieldGlobCorrDist, TFile& file, const o2::tpc::Side side)
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
    dumpLocalDistortions(file, side);

    start = timer::now();
    calcLocalDistortionsCorrections(1, formulas); // local correction calculation
    stop = timer::now();
    time = stop - start;
    std::cout << "local corrections analytical: " << time.count() << std::endl;
    dumpLocalCorrections(file, side);

    start = timer::now();
    calcGlobalDistortions(formulas);
    stop = timer::now();
    time = stop - start;
    std::cout << "global distortions analytical: " << time.count() << std::endl;
    dumpGlobalDistortions(file, side);

    start = timer::now();
    calcGlobalCorrections(formulas);
    stop = timer::now();
    time = stop - start;
    std::cout << "global corrections analytical: " << time.count() << std::endl;
    dumpGlobalCorrections(file, side);

    auto stopTotal = timer::now();
    time = stopTotal - startTotal;
    std::cout << "everything is done. Total Time: " << time.count() << std::endl;
    std::cout << std::endl;
  }
  if (mode == 1 || mode == 2) {
    std::cout << std::endl;
    auto startTotal = timer::now();

    if (mode == 1) {
      auto start = timer::now();
      fillBoundaryAndChargeDensities(formulas);
      auto stop = timer::now();
      std::chrono::duration<float> time = stop - start;
      std::cout << "filling boundary and charge density: " << time.count() << std::endl;
    } else {
      const std::string pathToRootFile{"/Users/matthias/alice/o2_macros/SpaceCharges/ATO-498/InputSCDensityHistograms_8000events.root"};
      TFile fileSCDensity(pathToRootFile.data(), "READ");
      const char* histoName = "inputSCDensity3D_8000_avg";
      fillChargeDensityFromHisto(fileSCDensity, histoName);
      fileSCDensity.Close();
    }
    dumpDensity(file, side);

    auto start = timer::now();
    poissonSolver(side, 300, 1e-8);
    auto stop = timer::now();
    std::chrono::duration<float> time = stop - start;
    std::cout << "poissonSolver: " << time.count() << std::endl;
    dumpPotential(file, side);

    start = timer::now();
    calcEField(side);
    stop = timer::now();
    time = stop - start;
    std::cout << "electric field calculation: " << time.count() << std::endl;
    dumpElectricFields(file, side);

    const auto numEFields = getElectricFieldsInterpolator(side);
    start = timer::now();
    calcLocalDistortionsCorrections(0, numEFields); // local distortion calculation
    stop = timer::now();
    time = stop - start;
    std::cout << "local distortions numerical: " << time.count() << std::endl;
    dumpLocalDistortions(file, side);

    start = timer::now();
    calcLocalDistortionsCorrections(1, numEFields); // local correction calculation
    stop = timer::now();
    time = stop - start;
    std::cout << "local corrections numerical: " << time.count() << std::endl;
    dumpLocalCorrections(file, side);

    start = timer::now();
    const auto lDistInterpolator = getLocalDistInterpolator(side);
    electricFieldGlobCorrDist ? calcGlobalDistortions(numEFields) : calcGlobalDistortions(lDistInterpolator);
    stop = timer::now();
    time = stop - start;
    std::cout << "global distortions with local distortions: " << time.count() << std::endl;
    dumpGlobalDistortions(file, side);

    start = timer::now();
    const auto lCorrInterpolator = getLocalCorrInterpolator(side);
    electricFieldGlobCorrDist ? calcGlobalCorrections(numEFields) : calcGlobalCorrections(lCorrInterpolator);
    stop = timer::now();
    time = stop - start;
    std::cout << "global corrections with local corrections: " << time.count() << std::endl;
    dumpGlobalCorrections(file, side);

    auto stopTotal = timer::now();
    time = stopTotal - startTotal;
    std::cout << "everything is done. Total Time: " << time.count() << std::endl;
    std::cout << std::endl;
  }
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::performGlobalCorrDist(TFile& file, const o2::tpc::Side side)
{

  using timer = std::chrono::high_resolution_clock;
  auto startTotal = timer::now();
  std::cout << std::endl;

  setElectricFieldsFromFile(file, side);

  auto start = timer::now();
  const auto numEFields = getElectricFieldsInterpolator(side);
  calcGlobalDistortions(numEFields);
  auto stop = timer::now();
  std::chrono::duration<float> time = stop - start;
  std::cout << "global distortions with local distortions: " << time.count() << std::endl;
  dumpGlobalDistortions(file, side);

  start = timer::now();
  calcGlobalCorrections(numEFields);
  stop = timer::now();
  time = stop - start;
  std::cout << "global corrections with local corrections: " << time.count() << std::endl;
  dumpGlobalCorrections(file, side);

  auto stopTotal = timer::now();
  time = stopTotal - startTotal;
  std::cout << "everything is done. Total Time: " << time.count() << std::endl;
  std::cout << std::endl;
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::setFromFile(TFile& file, const o2::tpc::Side side)
{
  setDensityFromFile(file, side);
  setPotentialFromFile(file, side);
  setElectricFieldsFromFile(file, side);
  setLocalDistortionsFromFile(file, side);
  setLocalCorrectionsFromFile(file, side);
  setGlobalDistortionsFromFile(file, side);
  setGlobalCorrectionsFromFile(file, side);
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::fillBoundaryAndChargeDensities(const AnalyticalFields<DataT>& formulaStruct)
{
  const o2::tpc::Side side = formulaStruct.getSide();
  for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
    const DataT phi = getPhiVertex(iPhi);
    for (size_t iR = 0; iR < Nr; ++iR) {
      const DataT radius = getRVertex(iR);
      for (size_t iZ = 0; iZ < Nz; ++iZ) {
        const DataT z = getZVertex(iZ);
        mDensity[side](iZ, iR, iPhi) = formulaStruct.evalDensity(z, radius, phi);

        if ((iR == 0) || (iR == (Nr - 1)) || (iZ == 0) || (iZ == (Nz - 1))) {
          mPotential[side](iZ, iR, iPhi) = formulaStruct.evalPotential(z, radius, phi);
        }
      }
    }
  }
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::poissonSolver(const o2::tpc::Side side, const int maxIteration, const DataT stoppingConvergence)
{
  // TODO MODIFY AliTPCPoissonSolver class to accept grid instead TMATRIXD
  TMatrixD* matricesPotential[Nphi];
  TMatrixD* matricesDensity[Nphi];
  for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
    matricesPotential[iPhi] = new TMatrixD(Nr, Nz);
    matricesDensity[iPhi] = new TMatrixD(Nr, Nz);
    for (size_t iR = 0; iR < Nr; ++iR) {
      for (size_t iZ = 0; iZ < Nz; ++iZ) {
        (*matricesPotential[iPhi])(iR, iZ) = mPotential[side](iZ, iR, iPhi);
        (*matricesDensity[iPhi])(iR, iZ) = mDensity[side](iZ, iR, iPhi);
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
        mPotential[side](iZ, iR, iPhi) = static_cast<DataT>((*matricesPotential[iPhi])(iR, iZ));
      }
    }
  }

  for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
    delete matricesPotential[iPhi];
    delete matricesDensity[iPhi];
  }
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::setEField(const AnalyticalFields<DataT>& formulaStruct)
{
  const o2::tpc::Side side = formulaStruct.getSide();
  for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
    for (size_t iR = 0; iR < Nr; ++iR) {
      for (size_t iZ = 0; iZ < Nz; ++iZ) {
        const DataT radius = getRVertex(iR);
        const DataT z = getZVertex(iZ);
        const DataT phi = getPhiVertex(iPhi);
        mElectricFieldEr[side](iZ, iR, iPhi) = formulaStruct.evalEr(z, radius, phi);
        mElectricFieldEz[side](iZ, iR, iPhi) = formulaStruct.evalEz(z, radius, phi);
        mElectricFieldEphi[side](iZ, iR, iPhi) = formulaStruct.evalEphi(z, radius, phi);
      }
    }
  }
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::calcEField(const o2::tpc::Side side)
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
        mElectricFieldEr[side](iZ, iR, iPhi) = -1 * (mPotential[side](iZ, iR + 1, iPhi) - mPotential[side](iZ, iR - 1, iPhi)) * static_cast<DataT>(0.5) * getInvSpacingR();                                     // r direction
        mElectricFieldEz[side](iZ, iR, iPhi) = -1 * (mPotential[side](iZ + 1, iR, iPhi) - mPotential[side](iZ - 1, iR, iPhi)) * static_cast<DataT>(0.5) * getInvSpacingZ();                                     // z direction
        mElectricFieldEphi[side](iZ, iR, iPhi) = -1 * (signPlus * mPotential[side](iZ, iR, tmpPlus) - signMinus * mPotential[side](iZ, iR, tmpMinusS)) * static_cast<DataT>(0.5) * getInvSpacingPhi() / radius; // phi direction
      }
    }

    // for boundary-r
    for (size_t iZ = 0; iZ < Nz; iZ++) {
      mElectricFieldEr[side](iZ, 0, iPhi) = -1 * (-static_cast<DataT>(0.5) * mPotential[side](iZ, 2, iPhi) + 2 * mPotential[side](iZ, 1, iPhi) - static_cast<DataT>(1.5) * mPotential[side](iZ, 0, iPhi)) * getInvSpacingR();                    // forward difference
      mElectricFieldEr[side](iZ, Nr - 1, iPhi) = -1 * (static_cast<DataT>(1.5) * mPotential[side](iZ, Nr - 1, iPhi) - 2 * mPotential[side](iZ, Nr - 2, iPhi) + static_cast<DataT>(0.5) * mPotential[side](iZ, Nr - 3, iPhi)) * getInvSpacingR(); // backward difference
    }

    for (size_t iR = 0; iR < Nr; iR += Nr - 1) {
      const DataT radius = getRVertex(iR);
      for (size_t iZ = 1; iZ < Nz - 1; iZ++) {
        mElectricFieldEz[side](iZ, iR, iPhi) = -1 * (mPotential[side](iZ + 1, iR, iPhi) - mPotential[side](iZ - 1, iR, iPhi)) * static_cast<DataT>(0.5) * getInvSpacingZ();                                     // z direction
        mElectricFieldEphi[side](iZ, iR, iPhi) = -1 * (signPlus * mPotential[side](iZ, iR, tmpPlus) - signMinus * mPotential[side](iZ, iR, tmpMinusS)) * static_cast<DataT>(0.5) * getInvSpacingPhi() / radius; // phi direction
      }
    }

    // for boundary-z
    for (size_t iR = 0; iR < Nr; ++iR) {
      mElectricFieldEz[side](0, iR, iPhi) = -1 * (-static_cast<DataT>(0.5) * mPotential[side](2, iR, iPhi) + 2 * mPotential[side](1, iR, iPhi) - static_cast<DataT>(1.5) * mPotential[side](0, iR, iPhi)) * getInvSpacingZ();
      mElectricFieldEz[side](Nz - 1, iR, iPhi) = -1 * (static_cast<DataT>(1.5) * mPotential[side](Nz - 1, iR, iPhi) - 2 * mPotential[side](Nz - 2, iR, iPhi) + static_cast<DataT>(0.5) * mPotential[side](Nz - 3, iR, iPhi)) * getInvSpacingZ();
    }

    for (size_t iR = 1; iR < Nr - 1; ++iR) {
      const DataT radius = getRVertex(iR);
      for (size_t iZ = 0; iZ < Nz; iZ += Nz - 1) {
        mElectricFieldEr[side](iZ, iR, iPhi) = -1 * (mPotential[side](iZ, iR + 1, iPhi) - mPotential[side](iZ, iR - 1, iPhi)) * static_cast<DataT>(0.5) * getInvSpacingR();                                     // r direction
        mElectricFieldEphi[side](iZ, iR, iPhi) = -1 * (signPlus * mPotential[side](iZ, iR, tmpPlus) - signMinus * mPotential[side](iZ, iR, tmpMinusS)) * static_cast<DataT>(0.5) * getInvSpacingPhi() / radius; // phi direction
      }
    }

    // corner points for EPhi
    for (size_t iR = 0; iR < Nr; iR += Nr - 1) {
      const DataT radius = getRVertex(iR);
      for (size_t iZ = 0; iZ < Nz; iZ += Nz - 1) {
        mElectricFieldEphi[side](iZ, iR, iPhi) = -1 * (signPlus * mPotential[side](iZ, iR, tmpPlus) - signMinus * mPotential[side](iZ, iR, tmpMinusS)) * static_cast<DataT>(0.5) * getInvSpacingPhi() / radius; // phi direction
      }
    }
  }
  mIsEfieldSet[side] = true;
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::calcGlobalDistWithGlobalCorrIterative(const DistCorrInterpolator<DataT, Nr, Nz, Nphi>& globCorr, const int maxIter, const DataT convZ, const DataT convR, const DataT convPhi, const DataT approachZ, const DataT approachR, const DataT approachPhi)
{
  // for nearest neighbour search see: https://doc.cgal.org/latest/Spatial_searching/Spatial_searching_2searching_with_point_with_info_8cpp-example.html
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef Kernel::Point_3 Point_3;
  typedef boost::tuple<Point_3, std::array<unsigned int, 3>> Point_Index;
  typedef CGAL::Search_traits_3<Kernel> Traits_base;
  typedef CGAL::Search_traits_adapter<Point_Index, CGAL::Nth_of_tuple_property_map<0, Point_Index>, Traits_base> Traits;
  typedef CGAL::Orthogonal_k_neighbor_search<Traits> K_neighbor_search;
  typedef K_neighbor_search::Tree Tree;

  const o2::tpc::Side side = globCorr.getSide();
  const int nPoints = Nz * Nr * Nphi; // maximum number of points

  std::vector<Point_3> points;
  points.reserve(nPoints);
  std::vector<std::array<unsigned int, 3>> indices;
  indices.reserve(nPoints);

  // loop over global corrections and calculate the positions of the global correction
  for (unsigned int iPhi = 0; iPhi < Nphi; ++iPhi) {
    const DataT phi = getPhiVertex(iPhi);
    for (unsigned int iR = 0; iR < Nr; ++iR) {
      const DataT radius = getRVertex(iR);
      for (unsigned int iZ = 0; iZ < Nz; ++iZ) {
        const DataT z = getZVertex(iZ);

        const DataT globalCorrR = mGlobalCorrdR[side](iZ, iR, iPhi);
        const DataT globalCorrRPhi = mGlobalCorrdRPhi[side](iZ, iR, iPhi);
        const DataT globalCorrZ = mGlobalCorrdZ[side](iZ, iR, iPhi);

        const DataT posRCorr = radius + globalCorrR; // position of global correction
        const DataT posPhiCorr = regulatePhi(phi + globalCorrRPhi / radius);
        const DataT posZCorr = z + globalCorrZ;

        // check if the postion lies in the TPC volume
        if (posRCorr >= RMIN && posRCorr <= RMAX && posZCorr >= ZMIN && posZCorr <= ZMAX) {
          // normalize coordinates to gridspacing for searching the nearest neighbour. Otherwise one would have to convert the coordinates to x,y,z
          points.emplace_back(Point_3(posZCorr * getInvSpacingZ(), posRCorr * getInvSpacingR(), posPhiCorr * getInvSpacingPhi()));
          indices.emplace_back(std::array<unsigned int, 3>{iZ, iR, iPhi});
        }
      }
    }
  }

  // Insert data points in the tree
  Tree tree(boost::make_zip_iterator(boost::make_tuple(points.begin(), indices.begin())), boost::make_zip_iterator(boost::make_tuple(points.end(), indices.end())));

#pragma omp parallel for
  for (unsigned int iZ = 0; iZ < Nz; ++iZ) {
    const DataT z = getZVertex(iZ);
    for (unsigned int iR = 0; iR < Nr; ++iR) {
      const DataT radius = getRVertex(iR);
      for (unsigned int iPhi = 0; iPhi < Nphi; ++iPhi) {
        const DataT phi = getPhiVertex(iPhi);

        // find nearest neighbour
        const Point_3 query(z * getInvSpacingZ(), radius * getInvSpacingR(), phi * getInvSpacingPhi());
        const int neighbors = 1; // we need only one nearest neighbour
        const K_neighbor_search search(tree, query, neighbors);

        const K_neighbor_search::iterator it = search.begin();
        const DataT nearestZ = static_cast<DataT>(boost::get<0>(it->first)[0]) * GRIDSPACINGZ;
        const DataT nearestR = static_cast<DataT>(boost::get<0>(it->first)[1]) * GRIDSPACINGR;
        const DataT nearestPhi = static_cast<DataT>(boost::get<0>(it->first)[2]) * GRIDSPACINGPHI;

        const unsigned int nearestiZ = boost::get<1>(it->first)[0];
        const unsigned int nearestiR = boost::get<1>(it->first)[1];
        const unsigned int nearestiPhi = boost::get<1>(it->first)[2];

        //
        //==========================================================================================
        //==== start algorithm: use tricubic upsampling to numerically approach the query point ====
        //==========================================================================================
        //
        // 1. calculate difference from nearest point to query point with stepwidth factor x
        // and approach the new point
        //
        DataT stepR = (radius - nearestR) * approachR;
        DataT stepZ = (z - nearestZ) * approachZ;
        DataT stepPhi = (phi - nearestPhi) * approachPhi;

        // needed to check for convergence
        DataT lastDistanceR = std::numeric_limits<DataT>::max();
        DataT lastDistanceZ = std::numeric_limits<DataT>::max();
        DataT lastDistancePhi = std::numeric_limits<DataT>::max();

        DataT corrdR = 0;
        DataT corrdRPhi = 0;
        DataT corrdZ = 0;
        const bool safe = false;

        for (int iter = 0; iter < maxIter; ++iter) {
          // 2. get new point coordinates
          const DataT rPos = regulateR(getRVertex(nearestiR) + stepR);
          const DataT zPos = regulateZ(getZVertex(nearestiZ) + stepZ);
          const DataT phiPosUnreg = getPhiVertex(nearestiPhi) + stepPhi;
          const DataT phiPos = regulatePhi(phiPosUnreg);

          // interpolate global correction at new point and calculate position of global correction
          corrdR = globCorr.evaldR(zPos, rPos, phiPos, safe);
          const DataT rpos = rPos + corrdR;

          const DataT corrPhi = globCorr.evaldRPhi(zPos, rPos, phiPos, safe) / rPos;
          corrdRPhi = corrPhi * rpos; // normalize to new r coordinate
          const DataT phipos = phiPosUnreg + corrPhi;

          corrdZ = globCorr.evaldZ(zPos, rPos, phiPos, safe);
          const DataT zpos = zPos + corrdZ;

          // calculate distance to new point
          const DataT distanceR = radius - rpos;
          const DataT distanceZ = z - zpos;
          const DataT distancePhi = phi - phipos;

          // approach desired coordinate
          stepR += distanceR * approachR;
          stepZ += distanceZ * approachZ;
          stepPhi += distancePhi * approachPhi;

          // check for convergence
          const DataT totaldistRDiv = lastDistanceR == 0 ? 0 : std::abs(1 - std::abs(distanceR / lastDistanceR));         // should be larger than 0
          const DataT totaldistZDiv = lastDistanceZ == 0 ? 0 : std::abs(1 - std::abs(distanceZ / lastDistanceZ));         // should be larger than 0
          const DataT totaldistPhiDiv = lastDistancePhi == 0 ? 0 : std::abs(1 - std::abs(distancePhi / lastDistancePhi)); // should be larger than 0

          // stop algorithm if converged
          if (totaldistRDiv <= convR && totaldistZDiv <= convZ && totaldistPhiDiv <= convPhi) {
            break;
          }

          // save current distance to desired coordinate for convergency check
          lastDistanceR = distanceR;
          lastDistanceZ = distanceZ;
          lastDistancePhi = distancePhi;
        }
        // set global distortions if algorithm converged or iterations exceed max numbers of iterations
        mGlobalDistdR[side](iZ, iR, iPhi) = -corrdR;
        mGlobalDistdRPhi[side](iZ, iR, iPhi) = -corrdRPhi;
        mGlobalDistdZ[side](iZ, iR, iPhi) = -corrdZ;
      }
    }
  }
}


template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
NumericalFields<DataT, Nr, Nz, Nphi> O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::getElectricFieldsInterpolator(const o2::tpc::Side side) const
{
  if (!mIsEfieldSet[side]) {
    std::cout << "============== E-Fields are not set! returning ==============" << std::endl;
  }
  NumericalFields<DataT, Nr, Nz, Nphi> numFields(mElectricFieldEr[side], mElectricFieldEz[side], mElectricFieldEphi[side], mGrid3D, side);
  return numFields;
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
DistCorrInterpolator<DataT, Nr, Nz, Nphi> O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::getLocalDistInterpolator(const o2::tpc::Side side) const
{
  if (!mIsLocalDistSet[side]) {
    std::cout << "============== local distortions not set! returning ==============" << std::endl;
  }

  DistCorrInterpolator<DataT, Nr, Nz, Nphi> numFields(mLocalDistdR[side], mLocalDistdZ[side], mLocalDistdRPhi[side], mGrid3D, side);
  return numFields;
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
DistCorrInterpolator<DataT, Nr, Nz, Nphi> O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::getLocalCorrInterpolator(const o2::tpc::Side side) const
{
  if (!mIsLocalCorrSet[side]) {
    std::cout << "============== local corrections not set! returning ==============" << std::endl;
  }
  DistCorrInterpolator<DataT, Nr, Nz, Nphi> numFields(mLocalCorrdR[side], mLocalCorrdZ[side], mLocalCorrdRPhi[side], mGrid3D, side);
  return numFields;
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
DistCorrInterpolator<DataT, Nr, Nz, Nphi> O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::getGlobalDistInterpolator(const o2::tpc::Side side) const
{
  if (!mIsGlobalDistSet[side]) {
    std::cout << "============== global distortions not set! returning ==============" << std::endl;
  }
  DistCorrInterpolator<DataT, Nr, Nz, Nphi> numFields(mGlobalDistdR[side], mGlobalDistdZ[side], mGlobalDistdRPhi[side], mGrid3D, side);
  return numFields;
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
DistCorrInterpolator<DataT, Nr, Nz, Nphi> O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::getGlobalCorrInterpolator(const o2::tpc::Side side) const
{
  if (!mIsGlobalCorrSet[side]) {
    std::cout << "============== global corrections not set! returning ==============" << std::endl;
  }
  DistCorrInterpolator<DataT, Nr, Nz, Nphi> numFields(mGlobalCorrdR[side], mGlobalCorrdZ[side], mGlobalCorrdRPhi[side], mGrid3D, side);
  return numFields;
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
int O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::dumpElectricFields(TFile& outf, const o2::tpc::Side side) const
{
  if (!mIsEfieldSet[side]) {
    std::cout << "============== E-Fields are not set! returning ==============" << std::endl;
    return 0;
  }
  const std::string sideName = getSideName(side);
  const int er = mElectricFieldEr[side].writeToFile(outf, Form("fieldEr_side%s", sideName.data()));
  const int ez = mElectricFieldEz[side].writeToFile(outf, Form("fieldEz_side%s", sideName.data()));
  const int ephi = mElectricFieldEphi[side].writeToFile(outf, Form("fieldEphi_side%s", sideName.data()));
  return er + ez + ephi;
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::setElectricFieldsFromFile(TFile& inpf, const o2::tpc::Side side)
{
  const std::string sideName = getSideName(side);
  mElectricFieldEr[side].initFromFile(inpf, Form("fieldEr_side%s", sideName.data()));
  mElectricFieldEz[side].initFromFile(inpf, Form("fieldEz_side%s", sideName.data()));
  mElectricFieldEphi[side].initFromFile(inpf, Form("fieldEphi_side%s", sideName.data()));
  mIsEfieldSet[side] = true;
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
int O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::dumpGlobalDistortions(TFile& outf, const o2::tpc::Side side) const
{
  if (!mIsGlobalDistSet[side]) {
    std::cout << "============== global distortions are not set! returning ==============" << std::endl;
    return 0;
  }
  const std::string sideName = getSideName(side);
  const int er = mGlobalDistdR[side].writeToFile(outf, Form("distR_side%s", sideName.data()));
  const int ez = mGlobalDistdZ[side].writeToFile(outf, Form("distZ_side%s", sideName.data()));
  const int ephi = mGlobalDistdRPhi[side].writeToFile(outf, Form("distRphi_side%s", sideName.data()));
  return er + ez + ephi;
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::setGlobalDistortionsFromFile(TFile& inpf, const o2::tpc::Side side)
{
  mIsGlobalDistSet[side] = true;
  const std::string sideName = getSideName(side);
  mGlobalDistdR[side].initFromFile(inpf, Form("distR_side%s", sideName.data()));
  mGlobalDistdZ[side].initFromFile(inpf, Form("distZ_side%s", sideName.data()));
  mGlobalDistdRPhi[side].initFromFile(inpf, Form("distRphi_side%s", sideName.data()));
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
int O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::dumpGlobalCorrections(TFile& outf, const o2::tpc::Side side) const
{
  if (!mIsGlobalCorrSet[side]) {
    std::cout << "============== global corrections are not set! returning ==============" << std::endl;
    return 0;
  }
  const std::string sideName = getSideName(side);
  const int er = mGlobalCorrdR[side].writeToFile(outf, Form("corrR_side%s", sideName.data()));
  const int ez = mGlobalCorrdZ[side].writeToFile(outf, Form("corrZ_side%s", sideName.data()));
  const int ephi = mGlobalCorrdRPhi[side].writeToFile(outf, Form("corrRPhi_side%s", sideName.data()));
  return er + ez + ephi;
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::setGlobalCorrectionsFromFile(TFile& inpf, const o2::tpc::Side side)
{
  mIsGlobalCorrSet[side] = true;
  const std::string sideName = getSideName(side);
  mGlobalCorrdR[side].initFromFile(inpf, Form("corrR_side%s", sideName.data()));
  mGlobalCorrdZ[side].initFromFile(inpf, Form("corrZ_side%s", sideName.data()));
  mGlobalCorrdRPhi[side].initFromFile(inpf, Form("corrRPhi_side%s", sideName.data()));
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
int O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::dumpLocalCorrections(TFile& outf, const o2::tpc::Side side) const
{
  if (!mIsLocalCorrSet[side]) {
    std::cout << "============== local corrections are not set! returning ==============" << std::endl;
    return 0;
  }
  const std::string sideName = getSideName(side);
  const int lCorrdR = mLocalCorrdR[side].writeToFile(outf, Form("lcorrR_side%s", sideName.data()));
  const int lCorrdZ = mLocalCorrdZ[side].writeToFile(outf, Form("lcorrZ_side%s", sideName.data()));
  const int lCorrdRPhi = mLocalCorrdRPhi[side].writeToFile(outf, Form("lcorrRPhi_side%s", sideName.data()));
  return lCorrdR + lCorrdZ + lCorrdRPhi;
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::setLocalCorrectionsFromFile(TFile& inpf, const o2::tpc::Side side)
{
  const std::string sideName = getSideName(side);
  const bool lCorrdR = mLocalCorrdR[side].initFromFile(inpf, Form("lcorrR_side%s", sideName.data()));
  const bool lCorrdZ = mLocalCorrdZ[side].initFromFile(inpf, Form("lcorrZ_side%s", sideName.data()));
  const bool lCorrdRPhi = mLocalCorrdRPhi[side].initFromFile(inpf, Form("lcorrRPhi_side%s", sideName.data()));
  if (lCorrdR && lCorrdZ && lCorrdRPhi) {
    mIsLocalCorrSet[side] = true;
  } else {
    mIsLocalCorrSet[side] = false;
  }
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
int O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::dumpLocalDistortions(TFile& outf, const o2::tpc::Side side) const
{
  if (!mIsLocalDistSet[side]) {
    std::cout << "============== local distortions are not set! returning ==============" << std::endl;
    return 0;
  }
  const std::string sideName = getSideName(side);
  const int lDistdR = mLocalDistdR[side].writeToFile(outf, Form("ldistR_side%s", sideName.data()));
  const int lDistdZ = mLocalDistdZ[side].writeToFile(outf, Form("ldistZ_side%s", sideName.data()));
  const int lDistdRPhi = mLocalDistdRPhi[side].writeToFile(outf, Form("ldistRPhi_side%s", sideName.data()));
  return lDistdR + lDistdZ + lDistdRPhi;
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::setLocalDistortionsFromFile(TFile& inpf, const o2::tpc::Side side)
{
  const std::string sideName = getSideName(side);
  const bool lDistdR = mLocalDistdR[side].initFromFile(inpf, Form("ldistR_side%s", sideName.data()));
  const bool lDistdZ = mLocalDistdZ[side].initFromFile(inpf, Form("ldistZ_side%s", sideName.data()));
  const bool lDistdRPhi = mLocalDistdRPhi[side].initFromFile(inpf, Form("ldistRPhi_side%s", sideName.data()));

  if (lDistdR && lDistdZ && lDistdRPhi) {
    mIsLocalDistSet[side] = true;
  } else {
    mIsLocalDistSet[side] = false;
  }
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::fillChargeDensityFromHisto(TFile& fInp, const char* name)
{
  TH3* hisSCDensity3D = (TH3*)fInp.Get(name);
  TH3D hRebin{};
  rebinDensityHisto(*hisSCDensity3D, hRebin);

  for (int side = o2::tpc::Side::A; side <= o2::tpc::SIDES; ++side) {
    std::cout << "side: " << side << std::endl;
    for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
      // const DataT phi = getPhiVertex(iPhi);
      for (size_t iR = 0; iR < Nr; ++iR) {
        // const DataT radius = getRVertex(iR);
        for (size_t iZ = 0; iZ < Nz; ++iZ) {
          const size_t zBin = side == o2::tpc::Side::A ? Nz + iZ + 1 : Nz - iZ; // TODO CHECK THIS!
          mDensity[side](iZ, iR, iPhi) = hRebin.GetBinContent(iPhi, iR, zBin);
        }
      }
    }
  }
}

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::rebinDensityHisto(const TH3& hOrig, TH3& hRebin) const
{

  const int nBinsPhiNew = Nphi;
  const int nBinsRNew = Nr;
  const int nBinsZNew = 2 * Nz;

  const auto phiLow = hOrig.GetXaxis()->GetBinLowEdge(1);
  const auto phiUp = hOrig.GetXaxis()->GetBinUpEdge(hOrig.GetNbinsX());
  const auto rLow = hOrig.GetYaxis()->GetBinLowEdge(1);
  const auto rUp = hOrig.GetYaxis()->GetBinUpEdge(hOrig.GetNbinsY());
  const auto zLow = hOrig.GetZaxis()->GetBinLowEdge(1);
  const auto zUp = hOrig.GetZaxis()->GetBinUpEdge(hOrig.GetNbinsZ());
  hRebin.SetBins(nBinsPhiNew, phiLow, phiUp, nBinsRNew, rLow, rUp, nBinsZNew, zLow, zUp);

  for (int iBinPhi = 1; iBinPhi <= nBinsPhiNew; ++iBinPhi) {
    const auto phiLowEdge = hRebin.GetXaxis()->GetBinLowEdge(iBinPhi);
    const auto phiUpEdge = hRebin.GetXaxis()->GetBinUpEdge(iBinPhi);

    const int phiLowBinOrig = hOrig.GetXaxis()->FindBin(phiLowEdge);
    const int phiUpBinOrig = hOrig.GetXaxis()->FindBin(phiUpEdge);

    // calculate the weights (area of original bin lies in the new bin / binwidthOrig) of the first and last bins
    const auto binWidthPhiOrig = hOrig.GetXaxis()->GetBinWidth(phiLowBinOrig);
    const auto lowerBinWeightPhi = std::abs(phiLowEdge - hOrig.GetXaxis()->GetBinUpEdge(phiLowBinOrig)) / binWidthPhiOrig;
    const auto upperBinWeightPhi = std::abs(phiUpEdge - hOrig.GetXaxis()->GetBinLowEdge(phiUpBinOrig)) / binWidthPhiOrig;

    for (int iBinR = 1; iBinR <= nBinsRNew; ++iBinR) {
      const auto rLowEdge = hRebin.GetYaxis()->GetBinLowEdge(iBinR);
      const auto rUpEdge = hRebin.GetYaxis()->GetBinUpEdge(iBinR);

      const int rLowBinOrig = hOrig.GetYaxis()->FindBin(rLowEdge);
      const int rUpBinOrig = hOrig.GetYaxis()->FindBin(rUpEdge);

      // calculate the weights (area of original bin lies in the new bin / binwidthOrig) of the first and last bins
      const auto binWidthROrig = hOrig.GetYaxis()->GetBinWidth(rLowBinOrig);
      const auto lowerBinWeightR = std::abs(rLowEdge - hOrig.GetYaxis()->GetBinUpEdge(rLowBinOrig)) / binWidthROrig;
      const auto upperBinWeightR = std::abs(rUpEdge - hOrig.GetYaxis()->GetBinLowEdge(rUpBinOrig)) / binWidthROrig;

      for (int iBinZ = 1; iBinZ <= nBinsZNew; ++iBinZ) {
        const auto zLowEdge = hRebin.GetZaxis()->GetBinLowEdge(iBinZ);
        const auto zUpEdge = hRebin.GetZaxis()->GetBinUpEdge(iBinZ);
        const auto zCenter = hRebin.GetZaxis()->GetBinCenter(iBinZ);

        int zLowBinOrig = hOrig.GetZaxis()->FindBin(zLowEdge);
        int zUpBinOrig = hOrig.GetZaxis()->FindBin(zUpEdge);
        const int currside = getSide(zCenter); // set the side of the current z-bin
        // get the side of the lowest and uppest bin from the orig histo
        const int sideLowOrig = getSide(hOrig.GetZaxis()->GetBinCenter(zLowBinOrig));
        const int sideUpOrig = getSide(hOrig.GetZaxis()->GetBinCenter(zUpBinOrig));

        // make bounds/side check of the zLowBinOrig and zUpBinOrig bins. They must be on the same side as the currside!!!
        if (currside != sideLowOrig && zLowBinOrig != zUpBinOrig) {
          // if the lower bins from the orig histo are not on the same side as the rebinned increase the binnumber until they are on the same side
          bool notequal = true;
          do {
            zLowBinOrig += 1;
            if (zLowBinOrig > zUpBinOrig) {
              std::cout << "SOMETHING WENT WRONG: SETTING BINS TO: " << zUpBinOrig << std::endl;
              zLowBinOrig = zUpBinOrig;
              notequal = false;
            }
            const int sideTmp = getSide(hOrig.GetZaxis()->GetBinCenter(zLowBinOrig));
            if (sideTmp == currside) {
              notequal = false;
            }
          } while (notequal);
        }

        if (currside != sideUpOrig && zLowBinOrig != zUpBinOrig) {
          // if the upper bins from the orig histo are not on the same side as the rebinned increase the binnumber until they are on the same side
          bool notequal = true;
          do {
            zUpBinOrig -= 1;
            if (zUpBinOrig < zLowBinOrig) {
              std::cout << "SOMETHING WENT WRONG: SETTING BINS TO: " << zLowBinOrig << std::endl;
              zUpBinOrig = zLowBinOrig;
              notequal = false;
            }
            const int sideTmp = getSide(hOrig.GetZaxis()->GetBinCenter(zUpBinOrig));
            if (sideTmp == currside) {
              notequal = false;
            }
          } while (notequal);
        }

        const auto binWidthZOrig = hOrig.GetZaxis()->GetBinWidth(zLowBinOrig);
        const auto lowerBinWeightZ = std::abs(zLowEdge - hOrig.GetZaxis()->GetBinUpEdge(zLowBinOrig)) / binWidthZOrig;
        const auto upperBinWeightZ = std::abs(zUpEdge - hOrig.GetZaxis()->GetBinLowEdge(zUpBinOrig)) / binWidthZOrig;

        // get the mean value of the original histogram of the found bin range
        DataT sum = 0;
        DataT sumW = 0;
        for (int iPhi = phiLowBinOrig; iPhi <= phiUpBinOrig; ++iPhi) {
          DataT weightPhi = 1;
          if (iPhi == phiLowBinOrig) {
            weightPhi = lowerBinWeightPhi;
          } else if (iPhi == phiUpBinOrig) {
            weightPhi = upperBinWeightPhi;
          }

          for (int iR = rLowBinOrig; iR <= rUpBinOrig; ++iR) {
            DataT weightR = 1;
            if (iR == rLowBinOrig) {
              weightR = lowerBinWeightR;
            } else if (iR == rUpBinOrig) {
              weightR = upperBinWeightR;
            }

            for (int iZ = zLowBinOrig; iZ <= zUpBinOrig; ++iZ) {
              DataT weightZ = 1;
              if (iZ == zLowBinOrig) {
                weightZ = lowerBinWeightZ;
              } else if (iZ == zUpBinOrig) {
                weightZ = upperBinWeightZ;
              }
              const auto val = hOrig.GetBinContent(iPhi, iR, iZ);
              // if(val==0){
              // what to do now???
              // }
              const auto totalWeight = weightPhi * weightR * weightZ;
              sum += val * totalWeight;
              sumW += totalWeight;
            }
          }
        }
        sum /= sumW;
        hRebin.SetBinContent(iBinPhi, iBinR, iBinZ, sum);
      }
    }
  }
  // dump histos to ttree for debugging
  // o2::utils::TreeStreamRedirector pcstream("rebinDebug.root", "UPDATE");
  // for (int iBinPhi = 1; iBinPhi <= nBinsPhiNew; ++iBinPhi) {
  //   for (int iBinR = 1; iBinR <= nBinsRNew; ++iBinR) {
  //     for (int iBinZ = 1; iBinZ <= nBinsZNew; ++iBinZ) {
  //       float binCont = hRebin.GetBinContent(iBinPhi, iBinR, iBinZ);
  //       float bincenterZ = hRebin.GetZaxis()->GetBinCenter(iBinZ);
  //       float bincenterR = hRebin.GetYaxis()->GetBinCenter(iBinR);
  //       float bincenterPhi = hRebin.GetXaxis()->GetBinCenter(iBinPhi);
  //
  //       (pcstream) << "rebin"
  //                  << "iBinPhi=" << iBinPhi << "iBinR=" << iBinR << "iBinZ=" << iBinZ << "binCont=" << binCont << "z=" << bincenterZ << "r=" << bincenterR << "phi" << bincenterPhi << "\n";
  //     }
  //   }
  // }
  // const int nBinsPhiOrig = hOrig.GetXaxis()->GetNbins();
  // const int nBinsROrig = hOrig.GetYaxis()->GetNbins();
  // const int nBinsZOrig = hOrig.GetZaxis()->GetNbins();
  //
  // for (int iBinPhi = 1; iBinPhi <= nBinsPhiOrig; ++iBinPhi) {
  //   for (int iBinR = 1; iBinR <= nBinsROrig; ++iBinR) {
  //     for (int iBinZ = 1; iBinZ <= nBinsZOrig; ++iBinZ) {
  //       float binCont = hOrig.GetBinContent(iBinPhi, iBinR, iBinZ);
  //       float bincenterZ = hOrig.GetZaxis()->GetBinCenter(iBinZ);
  //       float bincenterR = hOrig.GetYaxis()->GetBinCenter(iBinR);
  //       float bincenterPhi = hOrig.GetXaxis()->GetBinCenter(iBinPhi);
  //
  //       (pcstream) << "orig"
  //                   << "iBinPhi=" << iBinPhi << "iBinR=" << iBinR << "iBinZ=" << iBinZ << "binCont=" << binCont << "z=" << bincenterZ << "r=" << bincenterR << "phi" << bincenterPhi << "\n";
  //     }
  //   }
  // }
  // pcstream.GetFile()->cd();
  // pcstream.Close();
}

template class O2TPCSpaceCharge3DCalc<float, 17, 17, 90>;
template class O2TPCSpaceCharge3DCalc<float, 65, 65, 180>;
template class O2TPCSpaceCharge3DCalc<float, 129, 129, 180>;
template class O2TPCSpaceCharge3DCalc<float, 129, 129, 360>;
template class O2TPCSpaceCharge3DCalc<float, 257, 257, 360>;

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
