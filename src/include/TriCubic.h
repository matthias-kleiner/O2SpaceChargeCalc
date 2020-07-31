// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  TriCubic.h
/// \brief Definition of TriCubic class
///
/// \author  Matthias Kleiner <matthias.kleiner@cern.ch>

#ifndef TRI_H
#define TRI_H

#include "matrix.h"
#include "RegularGrid3D.h"
#include <omp.h>

template <typename DataT = float, typename Grid3D = RegularGrid3D<>>
class TriCubicInterpolator
{
  using VDataT = Vc::Vector<DataT>;

 public:
  /// Constructor for a tricubic interpolator
  /// \param gridData memory adress of the 3D-Grid struct
  TriCubicInterpolator(const Grid3D& gridData, const bool circularX = false, const bool circularY = false, const bool circularZ = false) : mGridData{gridData}, mCircularX(circularX), mCircularY(circularY), mCircularZ(circularZ){};

  // interpolate value at given coordinate
  /// \param x x coordinate
  /// \param y y coordinate
  /// \param z z coordinate
  /// \param safe checks if the given coordinates lie in the grid. if a value is out of the grid, the coordinate will be set to the border of the grid
  /// \return returns the interpolated value at given coordinate
  DataT operator()(const DataT x, const DataT y, const DataT z, const bool safe = true) const
  {
    const Vector<DataT, FDim> coordinates{{x, y, z}};  // vector holding the coordinates
    const auto relPos = processInp(coordinates, safe); // vector containing the relative position to
    const DataT res = interpolate(relPos);
    return res;
  }

  // interpolate derivative at given coordinate
  /// \param x x coordinate
  /// \param y y coordinate
  /// \param z z coordinate
  /// \param derx order of derivative d/dx: derx=1 -> d/dx f(x,y,z), derx=2 -> d^2/dx^2 f(x,y,z), derx=3 -> d^3/dx^3 f(x,y,z)
  /// \param derz order of derivative d/dy: dery=1 -> d/dy f(x,y,z), dery=2 -> d^2/dy^2 f(x,y,z), dery=3 -> d^3/dy^3 f(x,y,z)
  /// \param derz order of derivative d/dz: derz=1 -> d/dz f(x,y,z), derz=2 -> d^2/dz^2 f(x,y,z), derz=3 -> d^3/dz^3 f(x,y,z)
  /// derx=1 and dery=2 -> d/dx * d^2/dy^2 * f(x,y,z)
  /// \param safe checks if the given coordinates lie in the grid. if a value is out of the grid, the coordinate will be set to the border of the grid
  /// \return returns the interpolated derivative at given coordinate
  DataT operator()(const DataT x, const DataT y, const DataT z, const size_t derx, const size_t dery, const size_t derz, const bool safe) const
  {
    const Vector<DataT, FDim> coordinates{{x, y, z}}; // vector holding the coordinates
    const auto relPos = processInp(coordinates, safe);
    return evalDerivative(relPos[0], relPos[1], relPos[2], derx, dery, derz);
  }

  /// \return returns the maximum number of threads the tricubic interpolator can be used with (value should be = omp_get_max_threads())
  int getNThreads() const
  {
    return mNThreads;
  }

  /// \return returns the number of the thread. Each thread should have an individual thread number
  int getThreadNum() const
  {
    return mThreadnum;
  }

  /// \return performs a check if the interpolator can be used with current number of threads
  bool checkThreadSafety() const
  {
    return mNThreads >= omp_get_num_threads();
  }

 private:
  // matrix needed to compute the coefficients
  inline static Vc::Memory<VDataT, 64> mA[64]{
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-3, 3, 0, 0, 0, 0, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {2, -2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {9, -9, -9, 9, 0, 0, 0, 0, 6, 3, -6, -3, 0, 0, 0, 0, 6, -6, 3, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-6, 6, 6, -6, 0, 0, 0, 0, -3, -3, 3, 3, 0, 0, 0, 0, -4, 4, -2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -2, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {2, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-6, 6, 6, -6, 0, 0, 0, 0, -4, -2, 4, 2, 0, 0, 0, 0, -3, 3, -3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -1, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {4, -4, -4, 4, 0, 0, 0, 0, 2, 2, -2, -2, 0, 0, 0, 0, 2, -2, 2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, 0, 0, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, -1, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, -9, -9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 3, -6, -3, 0, 0, 0, 0, 6, -6, 3, -3, 0, 0, 0, 0, 4, 2, 2, 1, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -6, 6, 6, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, -3, 3, 3, 0, 0, 0, 0, -4, 4, -2, 2, 0, 0, 0, 0, -2, -2, -1, -1, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -6, 6, 6, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, -2, 4, 2, 0, 0, 0, 0, -3, 3, -3, 3, 0, 0, 0, 0, -2, -1, -2, -1, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -4, -4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, -2, -2, 0, 0, 0, 0, 2, -2, 2, -2, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0},
    {-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {9, -9, 0, 0, -9, 9, 0, 0, 6, 3, 0, 0, -6, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, -6, 0, 0, 3, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 2, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-6, 6, 0, 0, 6, -6, 0, 0, -3, -3, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, 4, 0, 0, -2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -2, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, -1, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, -9, 0, 0, -9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 3, 0, 0, -6, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, -6, 0, 0, 3, -3, 0, 0, 4, 2, 0, 0, 2, 1, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -6, 6, 0, 0, 6, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, -3, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, 4, 0, 0, -2, 2, 0, 0, -2, -2, 0, 0, -1, -1, 0, 0},
    {9, 0, -9, 0, -9, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 3, 0, -6, 0, -3, 0, 6, 0, -6, 0, 3, 0, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 9, 0, -9, 0, -9, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 3, 0, -6, 0, -3, 0, 6, 0, -6, 0, 3, 0, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 2, 0, 1, 0},
    {-27, 27, 27, -27, 27, -27, -27, 27, -18, -9, 18, 9, 18, 9, -18, -9, -18, 18, -9, 9, 18, -18, 9, -9, -18, 18, 18, -18, -9, 9, 9, -9, -12, -6, -6, -3, 12, 6, 6, 3, -12, -6, 12, 6, -6, -3, 6, 3, -12, 12, -6, 6, -6, 6, -3, 3, -8, -4, -4, -2, -4, -2, -2, -1},
    {18, -18, -18, 18, -18, 18, 18, -18, 9, 9, -9, -9, -9, -9, 9, 9, 12, -12, 6, -6, -12, 12, -6, 6, 12, -12, -12, 12, 6, -6, -6, 6, 6, 6, 3, 3, -6, -6, -3, -3, 6, 6, -6, -6, 3, 3, -3, -3, 8, -8, 4, -4, 4, -4, 2, -2, 4, 4, 2, 2, 2, 2, 1, 1},
    {-6, 0, 6, 0, 6, 0, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, -3, 0, 3, 0, 3, 0, -4, 0, 4, 0, -2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, -2, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, -6, 0, 6, 0, 6, 0, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, -3, 0, 3, 0, 3, 0, -4, 0, 4, 0, -2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, -2, 0, -1, 0, -1, 0},
    {18, -18, -18, 18, -18, 18, 18, -18, 12, 6, -12, -6, -12, -6, 12, 6, 9, -9, 9, -9, -9, 9, -9, 9, 12, -12, -12, 12, 6, -6, -6, 6, 6, 3, 6, 3, -6, -3, -6, -3, 8, 4, -8, -4, 4, 2, -4, -2, 6, -6, 6, -6, 3, -3, 3, -3, 4, 2, 4, 2, 2, 1, 2, 1},
    {-12, 12, 12, -12, 12, -12, -12, 12, -6, -6, 6, 6, 6, 6, -6, -6, -6, 6, -6, 6, 6, -6, 6, -6, -8, 8, 8, -8, -4, 4, 4, -4, -3, -3, -3, -3, 3, 3, 3, 3, -4, -4, 4, 4, -2, -2, 2, 2, -4, 4, -4, 4, -2, 2, -2, 2, -2, -2, -2, -2, -1, -1, -1, -1},
    {2, 0, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-6, 6, 0, 0, 6, -6, 0, 0, -4, -2, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, -3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -1, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {4, -4, 0, 0, -4, 4, 0, 0, 2, 2, 0, 0, -2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -6, 6, 0, 0, 6, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, -2, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, -3, 3, 0, 0, -2, -1, 0, 0, -2, -1, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -4, 0, 0, -4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, -2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 2, -2, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0},
    {-6, 0, 6, 0, 6, 0, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, 0, -2, 0, 4, 0, 2, 0, -3, 0, 3, 0, -3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, -1, 0, -2, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, -6, 0, 6, 0, 6, 0, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, 0, -2, 0, 4, 0, 2, 0, -3, 0, 3, 0, -3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, -1, 0, -2, 0, -1, 0},
    {18, -18, -18, 18, -18, 18, 18, -18, 12, 6, -12, -6, -12, -6, 12, 6, 12, -12, 6, -6, -12, 12, -6, 6, 9, -9, -9, 9, 9, -9, -9, 9, 8, 4, 4, 2, -8, -4, -4, -2, 6, 3, -6, -3, 6, 3, -6, -3, 6, -6, 3, -3, 6, -6, 3, -3, 4, 2, 2, 1, 4, 2, 2, 1},
    {-12, 12, 12, -12, 12, -12, -12, 12, -6, -6, 6, 6, 6, 6, -6, -6, -8, 8, -4, 4, 8, -8, 4, -4, -6, 6, 6, -6, -6, 6, 6, -6, -4, -4, -2, -2, 4, 4, 2, 2, -3, -3, 3, 3, -3, -3, 3, 3, -4, 4, -2, 2, -4, 4, -2, 2, -2, -2, -1, -1, -2, -2, -1, -1},
    {4, 0, -4, 0, -4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, -2, 0, -2, 0, 2, 0, -2, 0, 2, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 4, 0, -4, 0, -4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, -2, 0, -2, 0, 2, 0, -2, 0, 2, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0},
    {-12, 12, 12, -12, 12, -12, -12, 12, -8, -4, 8, 4, 8, 4, -8, -4, -6, 6, -6, 6, 6, -6, 6, -6, -6, 6, 6, -6, -6, 6, 6, -6, -4, -2, -4, -2, 4, 2, 4, 2, -4, -2, 4, 2, -4, -2, 4, 2, -3, 3, -3, 3, -3, 3, -3, 3, -2, -1, -2, -1, -2, -1, -2, -1},
    {8, -8, -8, 8, -8, 8, 8, -8, 4, 4, -4, -4, -4, -4, 4, 4, 4, -4, 4, -4, -4, 4, -4, 4, 4, -4, -4, 4, 4, -4, -4, 4, 2, 2, 2, 2, -2, -2, -2, -2, 2, 2, -2, -2, 2, 2, -2, -2, 2, -2, 2, -2, 2, -2, 2, -2, 1, 1, 1, 1, 1, 1, 1, 1}};

  inline static Matrix<DataT, 64> mMatrixA{mA};

  static constexpr unsigned int FDim = Grid3D::getDim(); // dimensions of the grid (only 3 supported)
  static constexpr unsigned int FX = Grid3D::getFX();    // index for x coordinate
  static constexpr unsigned int FY = Grid3D::getFY();    // index for y coordinate
  static constexpr unsigned int FZ = Grid3D::getFZ();    // index for z coordinate

  const Grid3D& mGridData{}; // adress to the data container of the grid

  const bool mCircularX{}; ///< sets circular padding in x dimension
  const bool mCircularY{}; ///< sets circular padding in y dimension
  const bool mCircularZ{}; ///< sets circular padding in z dimension

  inline static thread_local const int mThreadnum{omp_get_thread_num()}; ///< save for each thread the thread number to get fast access to the correct array
  const int mNThreads{omp_get_max_threads()};
  
  std::unique_ptr<Vector<DataT, 64>[]> mCoefficients = std::make_unique<Vector<DataT, 64>[]>(mNThreads);              ///< coefficients needed to interpolate a value
  std::unique_ptr<Vector<unsigned int, FDim>[]> mLastInd = std::make_unique<Vector<unsigned int, FDim>[]>(mNThreads); ///< stores the index for the cell, where the coefficients are already evaluated (only the coefficients for one-the last cell are stored)
  std::unique_ptr<bool[]> mInitialized = std::make_unique<bool[]>(mNThreads);                                         ///< sets the flag if the coefficients are evaluated at least once

  // use std::pow?
  DataT uiPow(const DataT base, unsigned int exponent) const;

  const Vector<DataT, 3> processInp(const Vector<DataT, 3>& coordinates, const bool safe) const;

  void calcCoefficients(const unsigned int ix, const unsigned int iy, const unsigned int iz) const;

  DataT interpolate(const Vector<DataT, 3>& pos) const;

  DataT evalDerivative(const DataT dx, const DataT dy, const DataT dz, const size_t derx, const size_t dery, const size_t derz) const;

  // for circular padding
  void getDataIndexCircularArray(const int index0, const int dim, int arr[]) const;

  // for non circular padding
  void getDataIndexNonCircularArray(const int index0, const int dim, int arr[]) const;

  // this helps to get circular and non circular padding indices
  int getRegulatedDelta(const int index0, const int delta, const int dim, const int offs) const
  {
    const int regulatedDelta = mGridData.isIndexInGrid(index0 + delta, dim) ? delta : offs;
    return regulatedDelta;
  }

  void initInterpolator(const unsigned int ix, const unsigned int iy, const unsigned int iz) const
  {
    calcCoefficients(ix, iy, iz);

    // store current cell
    mInitialized[mThreadnum] = true;
    mLastInd[mThreadnum][FX] = ix;
    mLastInd[mThreadnum][FY] = iy;
    mLastInd[mThreadnum][FZ] = iz;
  }
};

///
/// ========================================================================================================
///       Inline implementations of some methods
/// ========================================================================================================
///


template <typename DataT, typename Grid3D>
DataT TriCubicInterpolator<DataT, Grid3D>::evalDerivative(const DataT dx, const DataT dy, const DataT dz, const size_t derx, const size_t dery, const size_t derz) const
{

  DataT ret{};
  DataT cont{};

  for (size_t i = derx; i < 4; i++) {
    for (size_t j = dery; j < 4; j++) {
      for (size_t k = derz; k < 4; k++) {

        const size_t index = i + j * 4 + 16 * k;
        cont = mCoefficients[index] * uiPow(dx, i - derx) * uiPow(dy, j - dery) * uiPow(dz, k - derz);
        for (size_t w = 0; w < derx; w++) {
          cont *= (i - w);
        }
        for (size_t w = 0; w < dery; w++) {
          cont *= (j - w);
        }
        for (size_t w = 0; w < derz; w++) {
          cont *= (k - w);
        }
        ret += cont;
      }
    }
  }
  const auto invSpacing = mGridData.getInvSpacing();
  const DataT norm = uiPow(invSpacing[FX], derx) * uiPow(invSpacing[FY], dery) * uiPow(invSpacing[FZ], derz);
  return (ret * norm);
}

template <typename DataT, typename Grid3D>
void TriCubicInterpolator<DataT, Grid3D>::calcCoefficients(const unsigned int ix, const unsigned int iy, const unsigned int iz) const
{
  int deltaX[4]{};
  int deltaY[4]{};
  int deltaZ[4]{};

  // set padding type: circular or standard
  mCircularX ? getDataIndexCircularArray(ix, FX, deltaX) : getDataIndexNonCircularArray(ix, FX, deltaX);
  mCircularY ? getDataIndexCircularArray(iy, FY, deltaY) : getDataIndexNonCircularArray(iy, FY, deltaY);
  mCircularZ ? getDataIndexNonCircularArray(iz, FZ, deltaZ) : getDataIndexNonCircularArray(iz, FZ, deltaZ);

  // indices to datat storage
  const size_t i_x_y_z = mGridData.getDataIndex(ix, iy, iz);
  const size_t i_xp1_y_z = i_x_y_z + deltaX[2];
  const size_t i_x_yp1_z = i_x_y_z + deltaY[2];
  const size_t i_xp1_yp1_z = i_x_y_z + deltaX[2] + deltaY[2];
  const size_t i_x_y_zp1 = i_x_y_z + deltaZ[2];
  const size_t i_xp1_y_zp1 = i_x_y_z + deltaX[2] + deltaZ[2];
  const size_t i_x_yp1_zp1 = i_x_y_z + deltaY[2] + deltaZ[2];
  const size_t i_xp1_yp1_zp1 = i_x_y_z + deltaX[2] + deltaY[2] + deltaZ[2];

  const size_t i_xp2_y_z = i_x_y_z + deltaX[3];
  const size_t i_xm1_y_z = i_x_y_z + deltaX[1];
  const size_t i_xm1_yp1_z = i_x_y_z + deltaX[1] + deltaY[2];
  const size_t i_xp2_yp1_z = i_x_y_z + deltaX[3] + deltaY[2];
  const size_t i_xm1_y_zp1 = i_x_y_z + deltaX[1] + deltaZ[2];
  const size_t i_xp2_y_zp1 = i_x_y_z + deltaX[3] + deltaZ[2];
  const size_t i_xm1_yp1_zp1 = i_x_y_z + deltaX[1] + deltaY[2] + deltaZ[2];
  const size_t i_xp2_yp1_zp1 = i_x_y_z + deltaX[3] + deltaY[2] + deltaZ[2];

  const size_t i_x_ym1_z = i_x_y_z + deltaY[1];
  const size_t i_xp1_ym1_z = i_x_y_z + deltaX[2] + deltaY[1];
  const size_t i_x_yp2_z = i_x_y_z + deltaY[3];
  const size_t i_xp1_yp2_z = i_x_y_z + deltaX[2] + deltaY[3];
  const size_t i_x_ym1_zp1 = i_x_y_z + deltaY[1] + deltaZ[2];
  const size_t i_xp1_ym1_zp1 = i_x_y_z + deltaX[2] + deltaY[1] + deltaZ[2];
  const size_t i_x_yp2_zp1 = i_x_y_z + deltaY[3] + deltaZ[2];
  const size_t i_xp1_yp2_zp1 = i_x_y_z + deltaX[2] + deltaY[3] + deltaZ[2];

  const size_t i_x_y_zm1 = i_x_y_z + deltaZ[1];
  const size_t i_xp1_y_zm1 = i_x_y_z + deltaX[2] + deltaZ[1];
  const size_t i_x_yp1_zm1 = i_x_y_z + deltaY[2] + deltaZ[1];
  const size_t i_xp1_yp1_zm1 = i_x_y_z + deltaX[2] + deltaY[2] + deltaZ[1];
  const size_t i_x_y_zp2 = i_x_y_z + deltaZ[3];
  const size_t i_xp1_y_zp2 = i_x_y_z + deltaX[2] + deltaZ[3];
  const size_t i_x_yp1_zp2 = i_x_y_z + deltaY[2] + deltaZ[3];
  const size_t i_xp1_yp1_zp2 = i_x_y_z + deltaX[2] + deltaY[2] + deltaZ[3];

  const size_t i_xm1_ym1_z = i_x_y_z + deltaX[1] + deltaY[1];
  const size_t i_xp2_ym1_z = i_x_y_z + deltaX[3] + deltaY[1];
  const size_t i_xm1_yp2_z = i_x_y_z + deltaX[1] + deltaY[3];
  const size_t i_xp2_yp2_z = i_x_y_z + deltaX[3] + deltaY[3];
  const size_t i_xm1_ym1_zp1 = i_x_y_z + deltaX[1] + deltaY[1] + deltaZ[2];
  const size_t i_xp2_ym1_zp1 = i_x_y_z + deltaX[3] + deltaY[1] + deltaZ[2];
  const size_t i_xm1_yp2_zp1 = i_x_y_z + deltaX[1] + deltaY[3] + deltaZ[2];
  const size_t i_xp2_yp2_zp1 = i_x_y_z + deltaX[3] + deltaY[3] + deltaZ[2];

  const size_t i_xm1_y_zm1 = i_x_y_z + deltaX[1] + deltaZ[1];
  const size_t i_xp2_y_zm1 = i_x_y_z + deltaX[3] + deltaZ[1];
  const size_t i_xm1_yp1_zm1 = i_x_y_z + deltaX[1] + deltaY[2] + deltaZ[1];
  const size_t i_xp2_yp1_zm1 = i_x_y_z + deltaX[3] + deltaY[2] + deltaZ[1];
  const size_t i_xm1_y_zp2 = i_x_y_z + deltaX[1] + deltaZ[3];
  const size_t i_xp2_y_zp2 = i_x_y_z + deltaX[3] + deltaZ[3];
  const size_t i_xm1_yp1_zp2 = i_x_y_z + deltaX[1] + deltaY[2] + deltaZ[3];
  const size_t i_xp2_yp1_zp2 = i_x_y_z + deltaX[3] + deltaY[2] + deltaZ[3];

  const size_t i_x_ym1_zm1 = i_x_y_z + deltaY[1] + deltaZ[1];
  const size_t i_xp1_ym1_zm1 = i_x_y_z + deltaX[2] + deltaY[1] + deltaZ[1];
  const size_t i_x_yp2_zm1 = i_x_y_z + deltaY[3] + deltaZ[1];
  const size_t i_xp1_yp2_zm1 = i_x_y_z + deltaX[2] + deltaY[3] + deltaZ[1];
  const size_t i_x_ym1_zp2 = i_x_y_z + deltaY[1] + deltaZ[3];
  const size_t i_xp1_ym1_zp2 = i_x_y_z + deltaX[2] + deltaY[1] + deltaZ[3];
  const size_t i_x_yp2_zp2 = i_x_y_z + deltaY[3] + deltaZ[3];
  const size_t i_xp1_yp2_zp2 = i_x_y_z + deltaX[2] + deltaY[3] + deltaZ[3];

  const size_t i_xm1_ym1_zm1 = i_x_y_z + deltaX[1] + deltaY[1] + deltaZ[1];
  const size_t i_xp2_ym1_zm1 = i_x_y_z + deltaX[3] + deltaY[1] + deltaZ[1];
  const size_t i_xm1_yp2_zm1 = i_x_y_z + deltaX[1] + deltaY[3] + deltaZ[1];
  const size_t i_xp2_yp2_zm1 = i_x_y_z + deltaX[3] + deltaY[3] + deltaZ[1];
  const size_t i_xm1_ym1_zp2 = i_x_y_z + deltaX[1] + deltaY[1] + deltaZ[3];
  const size_t i_xp2_ym1_zp2 = i_x_y_z + deltaX[3] + deltaY[1] + deltaZ[3];
  const size_t i_xm1_yp2_zp2 = i_x_y_z + deltaX[1] + deltaY[3] + deltaZ[3];
  const size_t i_xp2_yp2_zp2 = i_x_y_z + deltaX[3] + deltaY[3] + deltaZ[3];

  // load values to tmp Vc
  const Vector<DataT, 24> vecDeriv1A{{mGridData[i_xp1_y_z], mGridData[i_xp2_y_z], mGridData[i_xp1_yp1_z], mGridData[i_xp2_yp1_z], mGridData[i_xp1_y_zp1], mGridData[i_xp2_y_zp1], mGridData[i_xp1_yp1_zp1], mGridData[i_xp2_yp1_zp1],
                                      mGridData[i_x_yp1_z], mGridData[i_xp1_yp1_z], mGridData[i_x_yp2_z], mGridData[i_xp1_yp2_z], mGridData[i_x_yp1_zp1], mGridData[i_xp1_yp1_zp1], mGridData[i_x_yp2_zp1], mGridData[i_xp1_yp2_zp1],
                                      mGridData[i_x_y_zp1], mGridData[i_xp1_y_zp1], mGridData[i_x_yp1_zp1], mGridData[i_xp1_yp1_zp1], mGridData[i_x_y_zp2], mGridData[i_xp1_y_zp2], mGridData[i_x_yp1_zp2], mGridData[i_xp1_yp1_zp2]}};

  const Vector<DataT, 24> vecDeriv1B{{mGridData[i_xm1_y_z], mGridData[i_x_y_z], mGridData[i_xm1_yp1_z], mGridData[i_x_yp1_z], mGridData[i_xm1_y_zp1], mGridData[i_x_y_zp1], mGridData[i_xm1_yp1_zp1], mGridData[i_x_yp1_zp1],
                                      mGridData[i_x_ym1_z], mGridData[i_xp1_ym1_z], mGridData[i_x_y_z], mGridData[i_xp1_y_z], mGridData[i_x_ym1_zp1], mGridData[i_xp1_ym1_zp1], mGridData[i_x_y_zp1], mGridData[i_xp1_y_zp1],
                                      mGridData[i_x_y_zm1], mGridData[i_xp1_y_zm1], mGridData[i_x_yp1_zm1], mGridData[i_xp1_yp1_zm1], mGridData[i_x_y_z], mGridData[i_xp1_y_z], mGridData[i_x_yp1_z], mGridData[i_xp1_yp1_z]}};

  const Vector<DataT, 24> vecDeriv2A{{mGridData[i_xp1_yp1_z], mGridData[i_xp2_yp1_z], mGridData[i_xp1_yp2_z], mGridData[i_xp2_yp2_z], mGridData[i_xp1_yp1_zp1], mGridData[i_xp2_yp1_zp1], mGridData[i_xp1_yp2_zp1], mGridData[i_xp2_yp2_zp1],
                                      mGridData[i_xp1_y_zp1], mGridData[i_xp2_y_zp1], mGridData[i_xp1_yp1_zp1], mGridData[i_xp2_yp1_zp1], mGridData[i_xp1_y_zp2], mGridData[i_xp2_y_zp2], mGridData[i_xp1_yp1_zp2], mGridData[i_xp2_yp1_zp2],
                                      mGridData[i_x_yp1_zp1], mGridData[i_xp1_yp1_zp1], mGridData[i_x_yp2_zp1], mGridData[i_xp1_yp2_zp1], mGridData[i_x_yp1_zp2], mGridData[i_xp1_yp1_zp2], mGridData[i_x_yp2_zp2], mGridData[i_xp1_yp2_zp2]}};

  const Vector<DataT, 24> vecDeriv2B{{mGridData[i_xm1_yp1_z], mGridData[i_x_yp1_z], mGridData[i_xm1_yp2_z], mGridData[i_x_yp2_z], mGridData[i_xm1_yp1_zp1], mGridData[i_x_yp1_zp1], mGridData[i_xm1_yp2_zp1], mGridData[i_x_yp2_zp1],
                                      mGridData[i_xm1_y_zp1], mGridData[i_x_y_zp1], mGridData[i_xm1_yp1_zp1], mGridData[i_x_yp1_zp1], mGridData[i_xm1_y_zp2], mGridData[i_x_y_zp2], mGridData[i_xm1_yp1_zp2], mGridData[i_x_yp1_zp2],
                                      mGridData[i_x_ym1_zp1], mGridData[i_xp1_ym1_zp1], mGridData[i_x_y_zp1], mGridData[i_xp1_y_zp1], mGridData[i_x_ym1_zp2], mGridData[i_xp1_ym1_zp2], mGridData[i_x_y_zp2], mGridData[i_xp1_y_zp2]}};

  const Vector<DataT, 24> vecDeriv2C{{mGridData[i_xp1_ym1_z], mGridData[i_xp2_ym1_z], mGridData[i_xp1_y_z], mGridData[i_xp2_y_z], mGridData[i_xp1_ym1_zp1], mGridData[i_xp2_ym1_zp1], mGridData[i_xp1_y_zp1], mGridData[i_xp2_y_zp1],
                                      mGridData[i_xp1_y_zm1], mGridData[i_xp2_y_zm1], mGridData[i_xp1_yp1_zm1], mGridData[i_xp2_yp1_zm1], mGridData[i_xp1_y_z], mGridData[i_xp2_y_z], mGridData[i_xp1_yp1_z], mGridData[i_xp2_yp1_z],
                                      mGridData[i_x_yp1_zm1], mGridData[i_xp1_yp1_zm1], mGridData[i_x_yp2_zm1], mGridData[i_xp1_yp2_zm1], mGridData[i_x_yp1_z], mGridData[i_xp1_yp1_z], mGridData[i_x_yp2_z], mGridData[i_xp1_yp2_z]}};

  const Vector<DataT, 24> vecDeriv2D{{mGridData[i_xm1_ym1_z], mGridData[i_x_ym1_z], mGridData[i_xm1_y_z], mGridData[i_x_y_z], mGridData[i_xm1_ym1_zp1], mGridData[i_x_ym1_zp1], mGridData[i_xm1_y_zp1], mGridData[i_x_y_zp1],
                                      mGridData[i_xm1_y_zm1], mGridData[i_x_y_zm1], mGridData[i_xm1_yp1_zm1], mGridData[i_x_yp1_zm1], mGridData[i_xm1_y_z], mGridData[i_x_y_z], mGridData[i_xm1_yp1_z], mGridData[i_x_yp1_z],
                                      mGridData[i_x_ym1_zm1], mGridData[i_xp1_ym1_zm1], mGridData[i_x_y_zm1], mGridData[i_xp1_y_zm1], mGridData[i_x_ym1_z], mGridData[i_xp1_ym1_z], mGridData[i_x_y_z], mGridData[i_xp1_y_z]}};

  const Vector<DataT, 8> vecDeriv3A{{mGridData[i_xp1_yp1_zp1], mGridData[i_xp2_yp1_zp1], mGridData[i_xp1_yp2_zp1], mGridData[i_xp2_yp2_zp1], mGridData[i_xp1_yp1_zp2], mGridData[i_xp2_yp1_zp2], mGridData[i_xp1_yp2_zp2], mGridData[i_xp2_yp2_zp2]}};

  const Vector<DataT, 8> vecDeriv3B{{mGridData[i_xm1_yp1_zp1], mGridData[i_x_yp1_zp1], mGridData[i_xm1_yp2_zp1], mGridData[i_x_yp2_zp1], mGridData[i_xm1_yp1_zp2], mGridData[i_x_yp1_zp2], mGridData[i_xm1_yp2_zp2], mGridData[i_x_yp2_zp2]}};

  const Vector<DataT, 8> vecDeriv3C{{mGridData[i_xp1_ym1_zp1], mGridData[i_xp2_ym1_zp1], mGridData[i_xp1_y_zp1], mGridData[i_xp2_y_zp1], mGridData[i_xp1_ym1_zp2], mGridData[i_xp2_ym1_zp2], mGridData[i_xp1_y_zp2], mGridData[i_xp2_y_zp2]}};

  const Vector<DataT, 8> vecDeriv3D{{mGridData[i_xm1_ym1_zp1], mGridData[i_x_ym1_zp1], mGridData[i_xm1_y_zp1], mGridData[i_x_y_zp1], mGridData[i_xm1_ym1_zp2], mGridData[i_x_ym1_zp2], mGridData[i_xm1_y_zp2], mGridData[i_x_y_zp2]}};

  const Vector<DataT, 8> vecDeriv3E{{mGridData[i_xp1_yp1_zm1], mGridData[i_xp2_yp1_zm1], mGridData[i_xp1_yp2_zm1], mGridData[i_xp2_yp2_zm1], mGridData[i_xp1_yp1_z], mGridData[i_xp2_yp1_z], mGridData[i_xp1_yp2_z], mGridData[i_xp2_yp2_z]}};

  const Vector<DataT, 8> vecDeriv3F{{mGridData[i_xm1_yp1_zm1], mGridData[i_x_yp1_zm1], mGridData[i_xm1_yp2_zm1], mGridData[i_x_yp2_zm1], mGridData[i_xm1_yp1_z], mGridData[i_x_yp1_z], mGridData[i_xm1_yp2_z], mGridData[i_x_yp2_z]}};

  const Vector<DataT, 8> vecDeriv3G{{mGridData[i_xp1_ym1_zm1], mGridData[i_xp2_ym1_zm1], mGridData[i_xp1_y_zm1], mGridData[i_xp2_y_zm1], mGridData[i_xp1_ym1_z], mGridData[i_xp2_ym1_z], mGridData[i_xp1_y_z], mGridData[i_xp2_y_z]}};

  const Vector<DataT, 8> vecDeriv3H{{mGridData[i_xm1_ym1_zm1], mGridData[i_x_ym1_zm1], mGridData[i_xm1_y_zm1], mGridData[i_x_y_zm1], mGridData[i_xm1_ym1_z], mGridData[i_x_ym1_z], mGridData[i_xm1_y_z], mGridData[i_x_y_z]}};

  const DataT fac1{0.5};
  const Vector<DataT, 24> vfac1{
    {fac1, fac1, fac1, fac1, fac1, fac1, fac1, fac1,
     fac1, fac1, fac1, fac1, fac1, fac1, fac1, fac1,
     fac1, fac1, fac1, fac1, fac1, fac1, fac1, fac1}};

  const DataT fac2{0.25};
  const Vector<DataT, 24> vfac2{
    {fac2, fac2, fac2, fac2, fac2, fac2, fac2, fac2,
     fac2, fac2, fac2, fac2, fac2, fac2, fac2, fac2,
     fac2, fac2, fac2, fac2, fac2, fac2, fac2, fac2}};

  const DataT fac3{0.125};
  const Vector<DataT, 8> vfac3{{fac3, fac3, fac3, fac3, fac3, fac3, fac3, fac3}};

  const Vector<DataT, 24> vecDeriv1Res{vfac1 * (vecDeriv1A - vecDeriv1B)};
  const Vector<DataT, 24> vecDeriv2Res{vfac2 * (vecDeriv2A - vecDeriv2B - vecDeriv2C + vecDeriv2D)};
  const Vector<DataT, 8> vecDeriv3Res{vfac3 * (vecDeriv3A - vecDeriv3B - vecDeriv3C + vecDeriv3D - vecDeriv3E + vecDeriv3F + vecDeriv3G - vecDeriv3H)};

  const Vector<DataT, 64> matrixPar{{mGridData[i_x_y_z], mGridData[i_xp1_y_z], mGridData[i_x_yp1_z], mGridData[i_xp1_yp1_z], mGridData[i_x_y_zp1], mGridData[i_xp1_y_zp1], mGridData[i_x_yp1_zp1], mGridData[i_xp1_yp1_zp1],
                                     vecDeriv1Res[0], vecDeriv1Res[1], vecDeriv1Res[2], vecDeriv1Res[3], vecDeriv1Res[4], vecDeriv1Res[5], vecDeriv1Res[6], vecDeriv1Res[7], vecDeriv1Res[8], vecDeriv1Res[9], vecDeriv1Res[10],
                                     vecDeriv1Res[11], vecDeriv1Res[12], vecDeriv1Res[13], vecDeriv1Res[14], vecDeriv1Res[15], vecDeriv1Res[16], vecDeriv1Res[17], vecDeriv1Res[18], vecDeriv1Res[19], vecDeriv1Res[20], vecDeriv1Res[21],
                                     vecDeriv1Res[22], vecDeriv1Res[23], vecDeriv2Res[0], vecDeriv2Res[1], vecDeriv2Res[2], vecDeriv2Res[3], vecDeriv2Res[4], vecDeriv2Res[5], vecDeriv2Res[6], vecDeriv2Res[7], vecDeriv2Res[8], vecDeriv2Res[9],
                                     vecDeriv2Res[10], vecDeriv2Res[11], vecDeriv2Res[12], vecDeriv2Res[13], vecDeriv2Res[14], vecDeriv2Res[15], vecDeriv2Res[16], vecDeriv2Res[17], vecDeriv2Res[18], vecDeriv2Res[19], vecDeriv2Res[20],
                                     vecDeriv2Res[21], vecDeriv2Res[22], vecDeriv2Res[23], vecDeriv3Res[0], vecDeriv3Res[1], vecDeriv3Res[2], vecDeriv3Res[3], vecDeriv3Res[4], vecDeriv3Res[5], vecDeriv3Res[6], vecDeriv3Res[7]}};

  // calc coeffiecients
  mCoefficients[mThreadnum] = mMatrixA * matrixPar;
}

template <typename DataT, typename Grid3D>
DataT TriCubicInterpolator<DataT, Grid3D>::interpolate(const Vector<DataT, 3>& pos) const
{
  // the formula for evaluating the interpolation is as follows:
  // f(x,y,z) = \sum_{i,j,k=0}^3 a_{ijk} * x^{i} * y^{j} * z^{k}
  // memValX contains the x^{i} values

  const Vector<DataT, FDim> vals0{{1, 1, 1}};
  const Vector<DataT, FDim> vals2{pos * pos};
  const Vector<DataT, FDim> vals3{vals2 * pos};

  const DataT valX[4]{vals0[FX], pos[FX], vals2[FX], vals3[FX]};
  const Vector<DataT, 64> vecValX{{valX[0], valX[1], valX[2], valX[3],
                                   valX[0], valX[1], valX[2], valX[3],
                                   valX[0], valX[1], valX[2], valX[3],
                                   valX[0], valX[1], valX[2], valX[3],
                                   valX[0], valX[1], valX[2], valX[3],
                                   valX[0], valX[1], valX[2], valX[3],
                                   valX[0], valX[1], valX[2], valX[3],
                                   valX[0], valX[1], valX[2], valX[3],
                                   valX[0], valX[1], valX[2], valX[3],
                                   valX[0], valX[1], valX[2], valX[3],
                                   valX[0], valX[1], valX[2], valX[3],
                                   valX[0], valX[1], valX[2], valX[3],
                                   valX[0], valX[1], valX[2], valX[3],
                                   valX[0], valX[1], valX[2], valX[3],
                                   valX[0], valX[1], valX[2], valX[3],
                                   valX[0], valX[1], valX[2], valX[3]}};

  const DataT valY[4]{vals0[FY], pos[FY], vals2[FY], vals3[FY]};
  const Vector<DataT, 64> vecValY{{valY[0], valY[0], valY[0], valY[0],
                                   valY[1], valY[1], valY[1], valY[1],
                                   valY[2], valY[2], valY[2], valY[2],
                                   valY[3], valY[3], valY[3], valY[3],
                                   valY[0], valY[0], valY[0], valY[0],
                                   valY[1], valY[1], valY[1], valY[1],
                                   valY[2], valY[2], valY[2], valY[2],
                                   valY[3], valY[3], valY[3], valY[3],
                                   valY[0], valY[0], valY[0], valY[0],
                                   valY[1], valY[1], valY[1], valY[1],
                                   valY[2], valY[2], valY[2], valY[2],
                                   valY[3], valY[3], valY[3], valY[3],
                                   valY[0], valY[0], valY[0], valY[0],
                                   valY[1], valY[1], valY[1], valY[1],
                                   valY[2], valY[2], valY[2], valY[2],
                                   valY[3], valY[3], valY[3], valY[3]}};

  const DataT valZ[4]{vals0[FZ], pos[FZ], vals2[FZ], vals3[FZ]};
  const Vector<DataT, 64> vecValZ{{valZ[0], valZ[0], valZ[0], valZ[0],
                                   valZ[0], valZ[0], valZ[0], valZ[0],
                                   valZ[0], valZ[0], valZ[0], valZ[0],
                                   valZ[0], valZ[0], valZ[0], valZ[0],
                                   valZ[1], valZ[1], valZ[1], valZ[1],
                                   valZ[1], valZ[1], valZ[1], valZ[1],
                                   valZ[1], valZ[1], valZ[1], valZ[1],
                                   valZ[1], valZ[1], valZ[1], valZ[1],
                                   valZ[2], valZ[2], valZ[2], valZ[2],
                                   valZ[2], valZ[2], valZ[2], valZ[2],
                                   valZ[2], valZ[2], valZ[2], valZ[2],
                                   valZ[2], valZ[2], valZ[2], valZ[2],
                                   valZ[3], valZ[3], valZ[3], valZ[3],
                                   valZ[3], valZ[3], valZ[3], valZ[3],
                                   valZ[3], valZ[3], valZ[3], valZ[3],
                                   valZ[3], valZ[3], valZ[3], valZ[3]}};

  const DataT result = sum(mCoefficients[mThreadnum] * vecValX * vecValY * vecValZ);
  return result;
}

template <typename DataT, typename Grid3D>
const Vector<DataT, 3> TriCubicInterpolator<DataT, Grid3D>::processInp(const Vector<DataT, 3>& coordinates, const bool safe) const
{
  const auto minGrid = mGridData.getGridMin();                      // vector containing the min x,y,z value of the grid
  const auto invSpacing = mGridData.getInvSpacing();                // vector containing the grid spacing for dim x,y,z
  Vector<DataT, FDim> posRel{(coordinates - minGrid) * invSpacing}; // needed for the grid index

  if (safe) {
    mGridData.clampToGrid(posRel);
  }

  const unsigned int ix = static_cast<unsigned int>(posRel[FX]);
  const unsigned int iy = static_cast<unsigned int>(posRel[FY]);
  const unsigned int iz = static_cast<unsigned int>(posRel[FZ]);

  const Vector<unsigned int, FDim> index{{ix, iy, iz}};

  if (!mInitialized[mThreadnum] || !(mLastInd[mThreadnum] == index)) {
    initInterpolator(index[FX], index[FY], index[FZ]);
  }

  const Vector<DataT, FDim> indexTmp{{static_cast<DataT>(ix), static_cast<DataT>(iy), static_cast<DataT>(iz)}};
  const Vector<DataT, FDim> relPos{posRel - indexTmp};
  return relPos;
}

// for circular padding
template <typename DataT, typename Grid3D>
void TriCubicInterpolator<DataT, Grid3D>::getDataIndexCircularArray(const int index0, const int dim, int arr[]) const
{
  const int delta_min1 = getRegulatedDelta(index0, -1, dim, mGridData.getN(dim) - 1);
  const int delta_min2 = getRegulatedDelta(index0, -2, dim, mGridData.getN(dim) - 2);
  const int delta_plus1 = getRegulatedDelta(index0, +1, dim, 1 - mGridData.getN(dim));
  const int delta_plus2 = getRegulatedDelta(index0, +2, dim, 2 - mGridData.getN(dim));

  arr[0] = mGridData.getDeltaDataIndex(delta_min2, dim);
  arr[1] = mGridData.getDeltaDataIndex(delta_min1, dim);
  arr[2] = mGridData.getDeltaDataIndex(delta_plus1, dim);
  arr[3] = mGridData.getDeltaDataIndex(delta_plus2, dim);
}

template <typename DataT, typename Grid3D>
void TriCubicInterpolator<DataT, Grid3D>::getDataIndexNonCircularArray(const int index0, const int dim, int arr[]) const
{
  const int delta_min1 = getRegulatedDelta(index0, -1, dim, 0);
  const int delta_min2 = getRegulatedDelta(index0, -2, dim, delta_min1);
  const int delta_plus1 = getRegulatedDelta(index0, +1, dim, 0);
  const int delta_plus2 = getRegulatedDelta(index0, +2, dim, delta_plus1);

  arr[0] = mGridData.getDeltaDataIndex(delta_min2, dim);
  arr[1] = mGridData.getDeltaDataIndex(delta_min1, dim);
  arr[2] = mGridData.getDeltaDataIndex(delta_plus1, dim);
  arr[3] = mGridData.getDeltaDataIndex(delta_plus2, dim);
}

template <typename DataT, typename Grid3D>
DataT TriCubicInterpolator<DataT, Grid3D>::uiPow(const DataT base, unsigned int exponent) const
{
  DataT result = 1;
  // infinite for loop
  for (;;) {
    // check if x is uneven number
    if (exponent & 1) {
      result *= base;
    }
    exponent >>= 1;
    if (!exponent) {
      break;
    }
    base *= base;
  }
  return result;
}

#endif
