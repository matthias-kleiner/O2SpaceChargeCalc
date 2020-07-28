#ifndef O2TPCSpaceCharge3DCalc_H
#define O2TPCSpaceCharge3DCalc_H

#include "TriCubic.h"
#include "AliRoot/AliTPCPoissonSolver.h"

// Root includes
#include "TFormula.h"

#ifdef WITH_OPENMP
#include <omp.h>
#endif

/// \tparam DataT the type of data which is used during the calculations
/// \tparam Nr number of vertices in r direction
/// \tparam Nz number of vertices in z direction
/// \tparam Nphi number of vertices in phi direction
template <typename DataT = float, unsigned int Nr = 4, unsigned int Nz = 4, unsigned int Nphi = 4>
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
  void calcLocalDistortionsCorrections();

  //step 4:
  // Info("AliTPCSpaceCharge3DCalc::InitSpaceCharge3DPoissonIntegralDz", "%s", Form("Step 4: Global correction/distortion cpu time: %f\n", w.CpuTime()));
  void calcGlobalDistortionsCorrections();

  constexpr DataT getGridSpacingR() const { return mGridSpacingR; }
  constexpr DataT getGridSpacingZ() const { return mGridSpacingZ; }
  constexpr DataT getGridSpacingPhi() const { return mGridSpacingPhi; }
  constexpr DataT getEzField() const { return (ASolv::fgkCathodeV - ASolv::fgkGG) / ASolv::fgkTPCZ0; }
  constexpr DataT getRMin() const { return ASolv::fgkIFCRadius; }
  RegularGrid3D<DataT, Nr, Nz, Nphi> getGrid3D() const { return mGrid3D; }

 private:
  using ASolv = AliTPCPoissonSolver<DataT>;

  static constexpr DataT mGridSpacingR = (ASolv::fgkOFCRadius - ASolv::fgkIFCRadius) / (Nr - 1); ///< grid spacing in r direction
  static constexpr DataT mGridSpacingZ = ASolv::fgkTPCZ0 / (Nz - 1);                             ///< grid spacing in z direction
  static constexpr DataT mGridSpacingPhi = 2 * M_PI / Nphi;                                      ///< grid spacing in phi direction // TODO CHANGE TO O2
  static constexpr DataT mRMin = ASolv::fgkIFCRadius;                                            ///< min radius
  static constexpr DataT mZMin = 0;                                                              ///< min z coordinate
  static constexpr DataT mPhiMin = 0;                                                            ///< min phi coordinate

  RegularGrid3D<DataT, Nr, Nz, Nphi> mGrid3D{mRMin, mZMin, mPhiMin, mGridSpacingR, mGridSpacingZ, mGridSpacingPhi}; ///< this grid contains the values for the local distortions/corrections, electric field etc.
};

template <typename DataT>
struct Formulas {

  DataT parA = -1e-5;
  DataT parB = 0.5;
  DataT parC = 1e-4;

  DataT evalEr(DataT r, DataT phi, DataT z) const
  {
    return erFunc(r, phi, z);
  }

  DataT evalEz(DataT r, DataT phi, DataT z) const
  {
    return ezFunc(r, phi, z);
  }

  DataT evalEphi(DataT r, DataT phi, DataT z) const
  {
    return ephiFunc(r, phi, z);
  }

  DataT evalPotential(DataT r, DataT phi, DataT z) const
  {
    return potentialFunc(r, phi, z);
  }

  DataT evalDensity(DataT r, DataT phi, DataT z) const
  {
    return densityFunc(r, phi, z);
  }

  /// analytical potential
  // std::string potential{"[0]*( (-x+254.5+83.5)^4 - 338.0 *(-x+254.5+83.5)^3 + 21250.75 * (-x+254.5+83.5)^2)*cos([1]*y)^2*exp(-1* [2] * (z-125)^2)"};
  std::function<DataT(DataT, DataT, DataT)> potentialFunc = [& parA = parA, &parB = parB, &parC = parC](DataT r, DataT phi, DataT z) {
    return -parA * (std::pow((-r + 254.5 + 83.5), 4) - 338.0 * std::pow((-r + 254.5 + 83.5), 3) + 21250.75 * std::pow((-r + 254.5 + 83.5), 2)) * std::cos(parB * phi) * std::cos(parB * phi) * std::exp(-1 * parC * (z - 125) * (z - 125));
  };

  /// analytical space charge density
  // const std::string secDerivativeR = "[0]*( 1/x * 16 * (-3311250 + 90995.5*x - 570.375*x^2 + x^3) )*cos([1]*y)^2*exp(-1* [2] * (z-125)^2)";
  // const std::string secDerivativePhi = "[0]*( (-x+254.5+83.5)^4 - 338.0 *(-x+254.5+83.5)^3 + 21250.75 * (-x+254.5+83.5)^2)/(x*x)*exp(-1* [2] * (z-125)^2) * - 2 * [1] * [1] * cos(2 * [1] * y) ";
  // const std::string secDerivativeZ = "[0]*( (-x+254.5+83.5)^4 - 338.0 *(-x+254.5+83.5)^3 + 21250.75 * (-x+254.5+83.5)^2)*cos([1]*y)^2 * 2*[2]*exp(-1* [2] * (z-125)^2) * (2*[2]* (z-125)^2 - 1  )";
  // const std::string rho = secDerivativeR + "+" + secDerivativePhi + "+" + secDerivativeZ;
  std::function<DataT(DataT, DataT, DataT)> densityFunc = [& parA = parA, &parB = parB, &parC = parC](DataT r, DataT phi, DataT z) {
    return -parA * ((1 / r * 16 * (-3311250 + 90995.5 * r - 570.375 * r * r + r * r * r)) * std::cos(parB * phi) * std::cos(parB * phi) * std::exp(-1 * parC * (z - 125) * (z - 125)) +
                    (std::pow(-r + 254.5 + 83.5, 4) - 338.0 * std::pow(-r + 254.5 + 83.5, 3) + 21250.75 * std::pow(-r + 254.5 + 83.5, 2)) / (r * r) * std::exp(-1 * parC * (z - 125) * (z - 125)) * -2 * parB * parB * std::cos(2 * parB * phi) +
                    (std::pow(-r + 254.5 + 83.5, 4) - 338.0 * std::pow(-r + 254.5 + 83.5, 3) + 21250.75 * std::pow(-r + 254.5 + 83.5, 2)) * std::cos(parB * phi) * std::cos(parB * phi) * 2 * parC * std::exp(-1 * parC * (z - 125) * (z - 125)) * (2 * parC * (z - 125) * (z - 125) - 1));
  };

  // std::string er{"[0]*4*(x^3 - 760.5*x^2 + 181991*x - 1.3245*10^7 )*cos([1]*y)^2*exp(-1* [2] * (z-125)^2)"};
  std::function<DataT(DataT, DataT, DataT)> erFunc = [& parA = parA, &parB = parB, &parC = parC](DataT r, DataT phi, DataT z) {
    return parA * 4 * (r * r * r - 760.5 * r * r + 181991 * r - 1.3245 * std::pow(10, 7)) * std::cos(parB * phi) * std::cos(parB * phi) * std::exp(-1 * parC * (z - 125) * (z - 125));
  };

  // std::string ephi{"[0]*( (-x+254.5+83.5)^4 - 338.0 *(-x+254.5+83.5)^3 + 21250.75 * (-x+254.5+83.5)^2)/x*exp(-1* [2] * (z-125)^2) * - [1] * sin(2*[1]*y)"};
  std::function<DataT(DataT, DataT, DataT)> ephiFunc = [& parA = parA, &parB = parB, &parC = parC](DataT r, DataT phi, DataT z) {
    return parA * (std::pow(-r + 254.5 + 83.5, 4) - 338.0 * std::pow(-r + 254.5 + 83.5, 3) + 21250.75 * (-r + 254.5 + 83.5) * (-r + 254.5 + 83.5)) / r * std::exp(-1 * parC * (z - 125) * (z - 125)) * -parB * std::sin(2 * parB * phi);
  };

  // std::string ez{"[0]*( (-x+254.5+83.5)^4 - 338.0 *(-x+254.5+83.5)^3 + 21250.75 * (-x+254.5+83.5)^2)*cos([1]*y)^2* -2 *[2]*(z-125) * exp(-1* [2] * (z-125)^2)"};
  std::function<DataT(DataT, DataT, DataT)> ezFunc = [& parA = parA, &parB = parB, &parC = parC](DataT r, DataT phi, DataT z) {
    return parA * (std::pow(-r + 254.5 + 83.5, 4) - 338.0 * std::pow(-r + 254.5 + 83.5, 3) + 21250.75 * (-r + 254.5 + 83.5) * (-r + 254.5 + 83.5)) * std::cos(parB * phi) * std::cos(parB * phi) * -2 * parC * (z - 125) * std::exp(-1 * parC * (z - 125) * (z - 125));
  };
};

#endif
