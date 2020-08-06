#ifndef AnalyticalFields_H
#define AnalyticalFields_H

#include <functional>
#include <cmath>
#include "TriCubic.h"

template <typename DataT = float>
struct AnalyticalFields {

  DataT parA{1e-5};                     ///< parameter [0] of functions
  DataT parB{0.5};                      ///< parameter [1] of functions
  DataT parC{1e-4};                     ///< parameter [2] of functions
  static constexpr unsigned int ID = 0; ///< needed to distinguish between the differrent structs

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
  using RegularGrid = RegularGrid3D<DataT, Nr, Nz, Nphi>;
  using DataContainer = DataContainer3D<DataT, Nr, Nz, Nphi>;
  using TriCubic = TriCubicInterpolator<DataT, Nr, Nz, Nphi>;
  NumericalFields(const DataContainer& gridErTmp, const DataContainer& gridEzTmp, const DataContainer& gridEphiTmp, const RegularGrid& gridProperties) : gridEr{gridErTmp}, gridEz{gridEzTmp}, gridEphi{gridEphiTmp}, gridInf{gridProperties} {};

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

  const DataContainer& gridEr{};   // adress to the data container of the grid
  const DataContainer& gridEz{};   // adress to the data container of the grid
  const DataContainer& gridEphi{}; // adress to the data container of the grid
  const RegularGrid& gridInf{};

  const bool circularZ = false;
  const bool circularR = false;
  const bool circularPhi = true;
  TriCubic interpolatorEr{gridEr, gridInf, circularZ, circularR, circularPhi};
  TriCubic interpolatorEz{gridEz, gridInf, circularZ, circularR, circularPhi};
  TriCubic interpolatorEphi{gridEphi, gridInf, circularZ, circularR, circularPhi};
  static constexpr unsigned int ID = 1; ///< needed to distinguish between the differrent structs
};

template <typename DataT = float, size_t Nr = 17, size_t Nz = 17, size_t Nphi = 90>
struct DistCorrInterpolator {
  using RegularGrid = RegularGrid3D<DataT, Nr, Nz, Nphi>;
  using DataContainer = DataContainer3D<DataT, Nr, Nz, Nphi>;
  using TriCubic = TriCubicInterpolator<DataT, Nr, Nz, Nphi>;

  DistCorrInterpolator(const DataContainer& distCorrdR, const DataContainer& distCorrdZ, const DataContainer& distCorrdRPhi, const RegularGrid& gridProperties) : distCorrdR{distCorrdR}, distCorrdZ{distCorrdZ}, distCorrdRPhi{distCorrdRPhi}, gridInf{gridProperties} {};

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for the local distortion or correction dR for given coordinate
  DataT evaldR(const DataT z, const DataT r, const DataT phi, const bool safe = true) const
  {
    return interpolatorDistCorrdR(z, r, phi, safe);
  }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for the local distortion or correction dZ for given coordinate
  DataT evaldZ(const DataT z, const DataT r, const DataT phi, const bool safe = true) const
  {
    return interpolatorDistCorrdZ(z, r, phi, safe);
  }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for the local distortion or correction dRPhi for given coordinate
  DataT evaldRPhi(const DataT z, const DataT r, const DataT phi, const bool safe = true) const
  {
    return interpolatorDistCorrdRPhi(z, r, phi, safe);
  }

  const DataContainer& distCorrdR{};    // adress to the data container of the grid
  const DataContainer& distCorrdZ{};    // adress to the data container of the grid
  const DataContainer& distCorrdRPhi{}; // adress to the data container of the grid
  const RegularGrid& gridInf{};
  const bool circularZ = false;
  const bool circularR = false;
  const bool circularPhi = true;
  TriCubic interpolatorDistCorrdR{distCorrdR, gridInf, circularZ, circularR, circularPhi};
  TriCubic interpolatorDistCorrdZ{distCorrdZ, gridInf, circularZ, circularR, circularPhi};
  TriCubic interpolatorDistCorrdRPhi{distCorrdRPhi, gridInf, circularZ, circularR, circularPhi};
  static constexpr unsigned int ID = 2; ///< needed to distinguish between the differrent structs
};

#endif
