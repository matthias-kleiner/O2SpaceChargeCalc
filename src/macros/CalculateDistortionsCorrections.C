#include "../include/O2TPCSpaceCharge3DCalc.h"

/// calculate distortions and corrections using an input histogram containing the density. All calculated values are dumped to one file
/// \param spaceChargeCalc O2TPCSpaceCharge3DCalc object in which the calculations are performed
/// \param fileOut output file where the calculated values are stored
/// \param side side of the TPC
/// \param path path to the root file containing the input histogram
/// \param histoName name of the histogram in the root file
/// \param globalEFieldType settings for global distortions/corrections: 0: using electric field for calculation of global distortions/corrections, 1: using local dis/corr interpolator for calculation of global distortions/corrections
/// \param globalDistType settings for global distortions: 0: standard method (start calculation of global distortion at each voxel in the tpc and follow electron drift to readout -slow-), 1: interpolation of global corrections (use the global corrections to apply an iterative approach to obtain the global distortions -fast-)
/// \param nThreads number of threads which are used (poisson solver uses only 2 cores)
template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void calculateDistortionsCorrections(o2::tpc::O2TPCSpaceCharge3DCalc<DataT, Nz, Nr, Nphi>& spaceChargeCalc, TFile& fileOut, const o2::tpc::Side side, const char* path, const char* histoName, const int globalEFieldType, const int globalDistType, const int nThreads = 1)
{
  using timer = std::chrono::high_resolution_clock;

  std::cout << "====== STARTING CALCULATIION OF DISTORTIONS AND CORRECTIONS BY USING A HISTOGRAM AS INPUT ======" << std::endl;
  const std::array sGlobalType{"Electric fields", "local distortion/correction interpolator"};
  std::cout << "calculation of global distortions and corrections are performed by using: " << sGlobalType[globalEFieldType] << std::endl;

  const std::array sGlobalDistType{"Standard method", "interpolation of global corrections"};
  std::cout << "calculation of global distortions performed by following method: " << sGlobalDistType[globalDistType] << std::endl;
  std::cout << std::endl;

  auto startTotal = timer::now();

  TFile fileOutSCDensity(path, "READ");
  spaceChargeCalc.fillChargeDensityFromHisto(fileOutSCDensity, histoName);
  fileOutSCDensity.Close();
  spaceChargeCalc.dumpDensity(fileOut, side);

  auto start = timer::now();
  omp_set_num_threads(2);
  spaceChargeCalc.poissonSolver(side, 300, 1e-8);
  auto stop = timer::now();
  std::chrono::duration<float> time = stop - start;
  std::cout << "poissonSolver: " << time.count() << std::endl;
  spaceChargeCalc.dumpPotential(fileOut, side);

  omp_set_num_threads(nThreads);
  start = timer::now();
  spaceChargeCalc.calcEField(side);
  stop = timer::now();
  time = stop - start;
  std::cout << "electric field calculation: " << time.count() << std::endl;
  spaceChargeCalc.dumpElectricFields(fileOut, side);

  const auto numEFields = spaceChargeCalc.getElectricFieldsInterpolator(side);
  start = timer::now();
  int type = o2::tpc::O2TPCSpaceCharge3DCalc<>::Type::Distortions;
  spaceChargeCalc.calcLocalDistortionsCorrections(type, numEFields); // local distortion calculation
  stop = timer::now();
  time = stop - start;
  std::cout << "local distortions numerical: " << time.count() << std::endl;
  spaceChargeCalc.dumpLocalDistortions(fileOut, side);

  start = timer::now();
  type = o2::tpc::O2TPCSpaceCharge3DCalc<>::Type::Corrections;
  spaceChargeCalc.calcLocalDistortionsCorrections(type, numEFields); // local correction calculation
  stop = timer::now();
  time = stop - start;
  std::cout << "local corrections numerical: " << time.count() << std::endl;
  spaceChargeCalc.dumpLocalCorrections(fileOut, side);

  start = timer::now();
  const auto lCorrInterpolator = spaceChargeCalc.getLocalCorrInterpolator(side);
  (globalEFieldType == 1) ? spaceChargeCalc.calcGlobalCorrections(lCorrInterpolator) : spaceChargeCalc.calcGlobalCorrections(numEFields);
  stop = timer::now();
  time = stop - start;
  std::cout << "global corrections: " << time.count() << std::endl;
  spaceChargeCalc.dumpGlobalCorrections(fileOut, side);

  start = timer::now();
  const auto lDistInterpolator = spaceChargeCalc.getLocalDistInterpolator(side);
  if (globalDistType == 0) {
    (globalEFieldType == 1) ? spaceChargeCalc.calcGlobalDistortions(lDistInterpolator) : spaceChargeCalc.calcGlobalDistortions(numEFields);
  } else {
    const auto globalCorrInterpolator = spaceChargeCalc.getGlobalCorrInterpolator(side);
    spaceChargeCalc.calcGlobalDistWithGlobalCorrIterative(globalCorrInterpolator);
  }
  stop = timer::now();
  time = stop - start;
  std::cout << "global distortions: " << time.count() << std::endl;
  spaceChargeCalc.dumpGlobalDistortions(fileOut, side);

  auto stopTotal = timer::now();
  time = stopTotal - startTotal;
  std::cout << "====== everything is done. Total Time: " << time.count() << " number of threads: " << omp_get_max_threads() << " ======" << std::endl;
  std::cout << std::endl;
}

/// \param spaceChargeCalc O2TPCSpaceCharge3DCalc object in which the calculations are performed
/// \param anaFields struct containing the analytical electric fields potential, space charge
/// \param fileOut output file where the calculated values are stored
/// \param side side of the TPC
/// \param globalEFieldType settings for global distortions/corrections: 0: using electric field for calculation of global distortions/corrections, 1: using local dis/corr interpolator for calculation of global distortions/corrections
/// \param globalDistType settings for global distortions: 0: standard method (start calculation of global distortion at each voxel in the tpc and follow electron drift to readout -slow-), 1: interpolation of global corrections (use the global corrections to apply an iterative approach to obtain the global distortions -fast-)
/// \param eFieldType setting for the electric field: 0: use analytical formula for the eletrical field for all calculations, 1: use the tricubic interpolator for the electric field
/// \param usePoissonSolver use poisson solver to calculate the potential or get the potential from the analytical formula 0: use analytical formula, 1: use poisson solver to calculate the potential (also calculates Efields using the obtained potential)
/// \param nThreads number of threads which are used
template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void calculateDistortionsCorrectionsAnalytical(o2::tpc::O2TPCSpaceCharge3DCalc<DataT, Nz, Nr, Nphi>& spaceChargeCalc, AnalyticalFields<DataT> anaFields, TFile& fileOut, const int globalEFieldType, const int eFieldType, const int globalDistType, const int usePoissonSolver, const int nThreads)
{
  using timer = std::chrono::high_resolution_clock;
  omp_set_num_threads(nThreads);

  std::cout << "====== STARTING CALCULATIION OF DISTORTIONS AND CORRECTIONS BY USING A ANALYTICAL FORMULA AS INPUT ======" << std::endl;
  const std::array sGlobalType{"Electric fields", "local distortion/correction interpolator"};
  std::cout << "calculation of global distortions and corrections are performed by using: " << sGlobalType[globalEFieldType] << std::endl;

  const std::array sGlobalDistType{"Standard method", "interpolation of global corrections"};
  std::cout << "calculation of global distortions performed by following method: " << sGlobalDistType[globalDistType] << std::endl;

  const std::array sEfieldType{"analytical formula", "tricubic interpolator"};
  std::cout << "Using the E-fields from: " << sEfieldType[eFieldType] << std::endl;

  const std::array sUsePoissonSolver{"analytical formula", "poisson solver"};
  std::cout << "Using the Potential from: " << sUsePoissonSolver[usePoissonSolver] << std::endl;
  std::cout << std::endl;

  const o2::tpc::Side side = anaFields.getSide();

  spaceChargeCalc.setChargeDensity(anaFields);
  spaceChargeCalc.dumpDensity(fileOut, side);

  if (usePoissonSolver == 1) {
    spaceChargeCalc.setPotentialBoundary(anaFields);
    auto start = timer::now();
    omp_set_num_threads(2);
    spaceChargeCalc.poissonSolver(side, 300, 1e-8);
    auto stop = timer::now();
    std::chrono::duration<float> time = stop - start;
    std::cout << "poissonSolver: " << time.count() << std::endl;
    omp_set_num_threads(nThreads);
  } else {
    spaceChargeCalc.setPotential(anaFields);
  }
  spaceChargeCalc.dumpPotential(fileOut, side);

  if (usePoissonSolver == 1) {
    auto start = timer::now();
    spaceChargeCalc.calcEField(side);
    auto stop = timer::now();
    std::chrono::duration<float> time = stop - start;
    std::cout << "electric field calculation: " << time.count() << std::endl;
  }
  else{
    spaceChargeCalc.setEField(anaFields);
  }
  spaceChargeCalc.dumpElectricFields(fileOut, side);

  auto startTotal = timer::now();
  const auto numEFields = spaceChargeCalc.getElectricFieldsInterpolator(side);

  auto start = timer::now();
  int type = o2::tpc::O2TPCSpaceCharge3DCalc<>::Type::Distortions;

  (eFieldType == 1) ? spaceChargeCalc.calcLocalDistortionsCorrections(type, numEFields) : spaceChargeCalc.calcLocalDistortionsCorrections(type, anaFields); // local distortion calculation

  auto stop = timer::now();
  std::chrono::duration<float> time = stop - start;
  std::cout << "local distortions analytical: " << time.count() << std::endl;
  spaceChargeCalc.dumpLocalDistortions(fileOut, side);

  start = timer::now();
  type = o2::tpc::O2TPCSpaceCharge3DCalc<>::Type::Corrections;
  (eFieldType == 1) ? spaceChargeCalc.calcLocalDistortionsCorrections(type, numEFields) : spaceChargeCalc.calcLocalDistortionsCorrections(type, anaFields); // local correction calculation
  stop = timer::now();
  time = stop - start;
  std::cout << "local corrections analytical: " << time.count() << std::endl;
  spaceChargeCalc.dumpLocalCorrections(fileOut, side);

  start = timer::now();
  const auto lCorrInterpolator = spaceChargeCalc.getLocalCorrInterpolator(side);
  if (globalEFieldType == 1) {
    spaceChargeCalc.calcGlobalCorrections(lCorrInterpolator);
  } else if (eFieldType == 1) {
    spaceChargeCalc.calcGlobalCorrections(numEFields);
  } else {
    spaceChargeCalc.calcGlobalCorrections(anaFields);
  }

  stop = timer::now();
  time = stop - start;
  std::cout << "global corrections analytical: " << time.count() << std::endl;
  spaceChargeCalc.dumpGlobalCorrections(fileOut, side);

  start = timer::now();
  const auto lDistInterpolator = spaceChargeCalc.getLocalDistInterpolator(side);

  if (globalDistType == 0) {
    if (globalEFieldType == 1) {
      spaceChargeCalc.calcGlobalDistortions(lDistInterpolator);
    } else if (eFieldType == 1) {
      spaceChargeCalc.calcGlobalDistortions(numEFields);
    } else {
      spaceChargeCalc.calcGlobalDistortions(anaFields);
    }
  } else {
    const auto globalCorrInterpolator = spaceChargeCalc.getGlobalCorrInterpolator(side);
    spaceChargeCalc.calcGlobalDistWithGlobalCorrIterative(globalCorrInterpolator);
  }

  stop = timer::now();
  time = stop - start;
  std::cout << "global distortions analytical: " << time.count() << std::endl;
  spaceChargeCalc.dumpGlobalDistortions(fileOut, side);

  auto stopTotal = timer::now();
  time = stopTotal - startTotal;
  std::cout << "====== everything is done. Total Time: " << time.count() << " number of threads: " << omp_get_max_threads() << " ======" << std::endl;
  std::cout << std::endl;
}
