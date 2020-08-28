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

#ifndef TRICUBIC_H
#define TRICUBIC_H

#include "Vector.h"
#include "RegularGrid3D.h"
#include "DataContainer3D.h"
#include <omp.h>

namespace o2
{
namespace tpc
{
///
/// The TriCubic class represents tricubic interpolation on a regular 3-Dim grid.
/// The algorithm which is used is based on the method developed by F. Lekien and J. Marsden and is described
/// in 'Tricubic Interpolation in Three Dimensions (2005)'  http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.89.7835
/// In this method in a first step 64 coefficients are computed by using a predefined 64*64 matrix.
/// These coefficients have to be computed for each cell in the grid, but are only computed when querying a point in a given cell.
/// The calculated coefficient is then stored for only the last cell and will be reused if the next query point lies in the same cell.
///
///
///  ----- creation of a tricubic interpolator -----
///
///  // define the grid:
///  // define number of vertices per dimension
///  const int xvertices = 40;
///  const int yvertices = 40;
///  const int zvertices = 40;
///
///  // define min range
///  float xmin = 0;
///  float ymin = 0;
///  float zmin = 0;
///
///  // define spacing between grid vertices
///  float xSpacing = 0.25;
///  float ySpacing = 0.25;
///  float zSpacing = 2 * M_PI / zvertices;
///
///  // create grid and datacontainer object
///  o2::tpc::RegularGrid3D<float, xvertices, yvertices, zvertices> grid3D(xmin, ymin, zmin, xSpacing, ySpacing, zSpacing);
///  DataContainer3D<float,xvertices,yvertices,zvertices> data3D;
///
///  // fill the DataContainer3D with some values
///  for (int ix = 0; ix < xvertices; ++ix) {
///   for (int iy = 0; iy < yvertices; ++iy) {
///     for (int iz = 0; iz < zvertices; ++iz) {
///       const float ixPos = xSpacing * ix + xmin;
///       const float iyPos = ySpacing * iy + ymin;
///       const float izPos = zSpacing * iz + zmin;
///       data3D(ix, iy, iz) = std::sin(iyPos * ixPos / 10.) + std::cos(izPos); // some arbitrary function is used here
///     }
///   }
/// }
/// // define periodicity of dimension
/// bool periodicX = false;
/// bool periodicY = false;
/// bool periodicZ = true;
///
/// // create tricubic interpolator
/// o2::tpc::TriCubicInterpolator<float, xvertices, yvertices, zvertices> interpolator(data3D, grid3D, periodicX, periodicY, periodicZ);
///
/// // query some values
/// for (float ix = grid3D.getGridMinX(); ix < grid3D.getGridMaxX(); ix += xSpacing / 3.) {
///   for (float iy = grid3D.getGridMinY(); iy < grid3D.getGridMaxY(); iy += ySpacing / 3.) {
///     for (float iz = grid3D.getGridMinZ() - 2 * zSpacing; iz < grid3D.getGridMaxZ() + 2 * zSpacing; iz += zSpacing / 3.) {
///       const float xQuery = ix;
///       const float yQuery = iy;
///       const float zQuery = iz;
///
///       const float interpolatedValue = interpolator(xQuery, yQuery, zQuery);
///       const float trueValue = std::sin(yQuery * xQuery / 10.) + std::cos(zQuery);
///
///       const float interpolatedDerivative = interpolator(xQuery, yQuery, zQuery, 1, 1, 0);
///       const float trueDerivative = 1 / 10. * std::cos(yQuery * xQuery / 10.) - yQuery / 10. * std::sin(yQuery * xQuery / 10.) * xQuery / 10.;
///     }
///   }
/// }
///

/// \tparam DataT the type of data which is used during the calculations
/// \tparam Nx number of vertices in r direction
/// \tparam Ny number of vertices in z direction
/// \tparam Nz number of vertices in phi direction
template <typename DataT = float, size_t Nx = 17, size_t Ny = 17, size_t Nz = 90>
class TriCubicInterpolator
{
  using Grid3D = RegularGrid3D<DataT, Nx, Ny, Nz>;
  using DataContainer = DataContainer3D<DataT, Nx, Ny, Nz>;
  using VDataT = Vc::Vector<DataT>;

 public:
  /// Constructor for a tricubic interpolator
  /// \param gridData struct containing access to the values of the grid
  /// \param gridProperties properties of the 3D grid
  /// \param circularX if set to true periodic boundary conditions are used in x direction
  /// \param circularY if set to true periodic boundary conditions are used in y direction
  /// \param circularZ if set to true periodic boundary conditions are used in z direction
  TriCubicInterpolator(const DataContainer& gridData, const Grid3D& gridProperties, const bool circularX = false, const bool circularY = false, const bool circularZ = false) : mGridData{gridData}, mGridProperties{gridProperties}, mCircular{{static_cast<int>(circularX), static_cast<int>(circularY), static_cast<int>(circularZ)}} {};

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

  /// interpolate derivative at given coordinate
  /// \param x x coordinate
  /// \param y y coordinate
  /// \param z z coordinate
  /// \param derx order of derivative d/dx: derx=1 -> d/dx f(x,y,z), derx=2 -> d^2/dx^2 f(x,y,z), derx=3 -> d^3/dx^3 f(x,y,z)
  /// \param derz order of derivative d/dy: dery=1 -> d/dy f(x,y,z), dery=2 -> d^2/dy^2 f(x,y,z), dery=3 -> d^3/dy^3 f(x,y,z)
  /// \param derz order of derivative d/dz: derz=1 -> d/dz f(x,y,z), derz=2 -> d^2/dz^2 f(x,y,z), derz=3 -> d^3/dz^3 f(x,y,z)
  /// derx=1 and dery=2 -> d/dx * d^2/dy^2 * f(x,y,z)
  /// \param safe checks if the given coordinates lie in the grid. if a value is out of the grid, the coordinate will be set to the border of the grid
  /// \return returns the interpolated derivative at given coordinate
  DataT operator()(const DataT x, const DataT y, const DataT z, const size_t derx, const size_t dery, const size_t derz, const bool safe = true) const
  {
    const Vector<DataT, FDim> coordinates{{x, y, z}}; // vector holding the coordinates
    const auto relPos = processInp(coordinates, safe);
    return evalDerivative(relPos[0], relPos[1], relPos[2], derx, dery, derz);
  }

  /// set which type of extrapolation is used at the grid boundaries (linear or parabol can be used with periodic z axis and non periodic x and y axis).
  /// \param extrapolationType sets type of extrapolation. See enum ExtrapolationType for different types
  void setExtrapolationType(const int extrapolationType)
  {
    mExtrapolationType = extrapolationType;
  }

  /// \return returns the extrapolation technique for missing boundary values/
  int getExtrapolationType() const
  {
    return mExtrapolationType;
  }

  /// \return returns the maximum number of threads the tricubic interpolator can be used with (value should be = omp_get_max_threads())
  int getNThreads() const
  {
    return mNThreads;
  }

  /// \return returns the number of the thread. Each thread should have an individual thread number
  int getThreadNum() const
  {
    return sThreadnum;
  }

  /// \return performs a check if the interpolator can be used with maximum number of threads
  bool checkThreadSafety() const
  {
    return mNThreads >= omp_get_max_threads();
  }

  enum ExtrapolationType {
    Linear = 0,
    Parabola = 1,
    None = 2
  };

 private:
  // matrix containing the 'relationship between the derivatives at the corners of the elements and the coefficients'
  inline static Vc::Memory<VDataT, 64> sMat[64]{
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
    {8, -8, -8, 8, -8, 8, 8, -8, 4, 4, -4, -4, -4, -4, 4, 4, 4, -4, 4, -4, -4, 4, -4, 4, 4, -4, -4, 4, 4, -4, -4, 4, 2, 2, 2, 2, -2, -2, -2, -2, 2, 2, -2, -2, 2, 2, -2, -2, 2, -2, 2, -2, 2, -2, 2, -2, 1, 1, 1, 1, 1, 1, 1, 1}}; ///< matrix containing the 'relationship between the derivatives at the corners of the elements and the coefficients'

  inline static Matrix<DataT, 64> sMatrixA{sMat}; ///< this matrix is used for vectorized operations with the 64*64 matrix

  static constexpr unsigned int FDim = Grid3D::getDim(); ///< dimensions of the grid
  static constexpr unsigned int FX = Grid3D::getFX();    ///< index for x coordinate
  static constexpr unsigned int FY = Grid3D::getFY();    ///< index for y coordinate
  static constexpr unsigned int FZ = Grid3D::getFZ();    ///< index for z coordinate

  const DataContainer& mGridData{}; ///< adress to the data container of the grid
  const Grid3D& mGridProperties{};  ///< adress to the properties of the grid
  const Vector<int, 3> mCircular{}; ///< information of periodic boundary conditions are used for x,y,z

  inline static thread_local const size_t sThreadnum{static_cast<size_t>(omp_get_thread_num())}; ///< save for each thread the thread number to get fast access to the correct array
  const int mNThreads{omp_get_max_threads()};                                                    ///< number of threads the tricubic interpolator can be used with

  std::unique_ptr<Vector<DataT, 64>[]> mCoefficients = std::make_unique<Vector<DataT, 64>[]>(mNThreads);              ///< coefficients needed to interpolate a value
  std::unique_ptr<Vector<unsigned int, FDim>[]> mLastInd = std::make_unique<Vector<unsigned int, FDim>[]>(mNThreads); ///< stores the index for the cell, where the coefficients are already evaluated (only the coefficients for the last cell are stored)
  std::unique_ptr<bool[]> mInitialized = std::make_unique<bool[]>(mNThreads);                                         ///< sets the flag if the coefficients are evaluated at least once

  int mExtrapolationType = ExtrapolationType::Parabola; ///< sets which type of extrapolation for missing points at boundary is used. Linear and Parabola is only supported for perdiodic z axis and non periodic x and y axis;

  //                 DEFINITION OF enum GridPos
  //========================================================
  //              Y
  //              |             6------F---7
  //              |           / |        / |
  //              |         K   G YR   L   H
  //              |       /     |    /     |
  //              |      2---B------3      |
  //              |      |      |   |      |
  //              |      |      4---|---E--5
  //              |      C XL  /    D XR  /
  //              |      |   I  YL  |   J
  //              |      | /        | /
  //              |      0---A------1
  //              |------------------------------- X
  //            /
  //          /
  //        /
  //      Z
  //========================================================

  enum GridPos {
    InnerVolume = 26,
    Edge0 = 0,
    Edge1 = 1,
    Edge2 = 2,
    Edge3 = 3,
    Edge4 = 4,
    Edge5 = 5,
    Edge6 = 6,
    Edge7 = 7,
    LineA = 8,
    LineB = 9,
    LineC = 10,
    LineD = 11,
    LineE = 12,
    LineF = 13,
    LineG = 14,
    LineH = 15,
    LineI = 16,
    LineJ = 17,
    LineK = 18,
    LineL = 19,
    SideXRight = 20,
    SideXLeft = 21,
    SideYRight = 22,
    SideYLeft = 23,
    SideZRight = 24,
    SideZLeft = 25
  };

  // use std::pow?
  DataT uiPow(DataT base, unsigned int exponent) const;

  const Vector<DataT, 3> processInp(const Vector<DataT, 3>& coordinates, const bool safe) const;

  // calculate the coefficients needed for the interpolation using the 64*64 matrix.
  // this is the 'slow' part of the code and might be optimized
  void calcCoefficients(const unsigned int ix, const unsigned int iy, const unsigned int iz) const;

  void calcCoefficientsExtrapolation(const unsigned int ix, const unsigned int iy, const unsigned int iz) const;

  DataT interpolate(const Vector<DataT, 3>& pos) const;

  DataT evalDerivative(const DataT dx, const DataT dy, const DataT dz, const size_t derx, const size_t dery, const size_t derz) const;

  // for periodic boundary conditions
  void getDataIndexCircularArray(const int index0, const int dim, int arr[]) const;

  // for non periodic boundary conditions
  void getDataIndexNonCircularArray(const int index0, const int dim, int arr[]) const;

  // this helps to get circular and non circular padding indices
  int getRegulatedDelta(const int index0, const int delta, const unsigned int dim, const int offs) const
  {
    return mGridProperties.isIndexInGrid(index0 + delta, dim) ? delta : offs;
  }

  void initInterpolator(const unsigned int ix, const unsigned int iy, const unsigned int iz) const;

  DataT extrapolation(const DataT valk, const DataT valk1, const int valk2) const;

  DataT linearExtrapolation(const DataT valk, const DataT valk1) const;

  DataT parabolExtrapolation(const DataT valk, const DataT valk1, const DataT valk2) const;

  int findPos(const int ix, const int iy, const int iz) const;

  bool isInInnerVolume(const int ix, const int iy, const int iz, int& posType) const
  {
    if (ix >= 1 && ix < Nx - 2 && iy >= 1 && iy < Ny - 2 && iz >= 1 && iz < Nz - 2) {
      posType = InnerVolume;
      return true;
    }
    return false;
  }

  bool findEdge(const int ix, const int iy, const int iz, int& posType) const;

  bool findLine(const int ix, const int iy, const int iz, int& posType) const;

  bool findSide(const int ix, const int iy, const int iz, int& posType) const;

  bool isSideRight(const int ind, const int dim) const
  {
    if (ind == mGridProperties.getN(dim) - 2) {
      return true;
    }
    return false;
  }

  bool isSideLeft(const int ind) const
  {
    if (ind == 0) {
      return true;
    }
    return false;
  }
};

///
/// ========================================================================================================
///                                         Inline implementations
/// ========================================================================================================
///

template <typename DataT, size_t Nx, size_t Ny, size_t Nz>
void TriCubicInterpolator<DataT, Nx, Ny, Nz>::initInterpolator(const unsigned int ix, const unsigned int iy, const unsigned int iz) const
{
  switch (mExtrapolationType) {
    case ExtrapolationType::None:
      calcCoefficients(ix, iy, iz);
      break;
    default:
      calcCoefficientsExtrapolation(ix, iy, iz);
      break;
  }

  // store current cell
  mInitialized[sThreadnum] = true;
  mLastInd[sThreadnum][FX] = ix;
  mLastInd[sThreadnum][FY] = iy;
  mLastInd[sThreadnum][FZ] = iz;
}

template <typename DataT, size_t Nx, size_t Ny, size_t Nz>
DataT TriCubicInterpolator<DataT, Nx, Ny, Nz>::evalDerivative(const DataT dx, const DataT dy, const DataT dz, const size_t derx, const size_t dery, const size_t derz) const
{
  //TODO optimize this
  DataT ret{};
  for (size_t i = derx; i < 4; i++) {
    for (size_t j = dery; j < 4; j++) {
      for (size_t k = derz; k < 4; k++) {

        const size_t index = i + j * 4 + 16 * k;
        DataT cont = mCoefficients[sThreadnum][index] * uiPow(dx, i - derx) * uiPow(dy, j - dery) * uiPow(dz, k - derz);
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
  const auto invSpacing = mGridProperties.getInvSpacing();
  const DataT norm = uiPow(invSpacing[FX], derx) * uiPow(invSpacing[FY], dery) * uiPow(invSpacing[FZ], derz);
  return (ret * norm);
}

template <typename DataT, size_t Nx, size_t Ny, size_t Nz>
void TriCubicInterpolator<DataT, Nx, Ny, Nz>::calcCoefficients(const unsigned int ix, const unsigned int iy, const unsigned int iz) const
{
  int deltaX[3]{};
  int deltaY[3]{};
  int deltaZ[3]{};

  // set padding type: circular or standard
  mCircular[0] ? getDataIndexCircularArray(ix, FX, deltaX) : getDataIndexNonCircularArray(ix, FX, deltaX);
  mCircular[1] ? getDataIndexCircularArray(iy, FY, deltaY) : getDataIndexNonCircularArray(iy, FY, deltaY);
  mCircular[2] ? getDataIndexCircularArray(iz, FZ, deltaZ) : getDataIndexNonCircularArray(iz, FZ, deltaZ);

  // indices to datat storage
  const int i0 = 0;
  const int i1 = 1;
  const int i2 = 2;
  const size_t ii_x_y_z = mGridData.getDataIndex(ix, iy, iz);
  const DataT i_x_y_z = mGridData[ii_x_y_z];
  const DataT i_xp1_y_z = mGridData[ii_x_y_z + deltaX[i1]];
  const DataT i_x_yp1_z = mGridData[ii_x_y_z + deltaY[i1]];
  const DataT i_xp1_yp1_z = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1]];
  const DataT i_x_y_zp1 = mGridData[ii_x_y_z + deltaZ[i1]];
  const DataT i_xp1_y_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i1]];
  const DataT i_x_yp1_zp1 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i1]];
  const DataT i_xp1_yp1_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i1]];

  const DataT i_xp2_y_z = mGridData[ii_x_y_z + deltaX[i2]];
  const DataT i_xm1_y_z = mGridData[ii_x_y_z + deltaX[i0]];
  const DataT i_xm1_yp1_z = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1]];
  const DataT i_xp2_yp1_z = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1]];
  const DataT i_xm1_y_zp1 = mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i1]];
  const DataT i_xp2_y_zp1 = mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i1]];
  const DataT i_xm1_yp1_zp1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i1]];
  const DataT i_xp2_yp1_zp1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1] + deltaZ[i1]];

  const DataT i_x_ym1_z = mGridData[ii_x_y_z + deltaY[i0]];
  const DataT i_xp1_ym1_z = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0]];
  const DataT i_x_yp2_z = mGridData[ii_x_y_z + deltaY[i2]];
  const DataT i_xp1_yp2_z = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2]];
  const DataT i_x_ym1_zp1 = mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i1]];
  const DataT i_xp1_ym1_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i1]];
  const DataT i_x_yp2_zp1 = mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i1]];
  const DataT i_xp1_yp2_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2] + deltaZ[i1]];

  const DataT i_x_y_zm1 = mGridData[ii_x_y_z + deltaZ[i0]];
  const DataT i_xp1_y_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i0]];
  const DataT i_x_yp1_zm1 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i0]];
  const DataT i_xp1_yp1_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i0]];
  const DataT i_x_y_zp2 = mGridData[ii_x_y_z + deltaZ[i2]];
  const DataT i_xp1_y_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i2]];
  const DataT i_x_yp1_zp2 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i2]];
  const DataT i_xp1_yp1_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i2]];

  const DataT i_xm1_ym1_z = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i0]];
  const DataT i_xp2_ym1_z = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0]];
  const DataT i_xm1_yp2_z = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2]];
  const DataT i_xp2_yp2_z = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i2]];
  const DataT i_xm1_ym1_zp1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i0] + deltaZ[i1]];
  const DataT i_xp2_ym1_zp1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0] + deltaZ[i1]];
  const DataT i_xm1_yp2_zp1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2] + deltaZ[i1]];
  const DataT i_xp2_yp2_zp1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i2] + deltaZ[i1]];

  const DataT i_xm1_y_zm1 = mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i0]];
  const DataT i_xp2_y_zm1 = mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i0]];
  const DataT i_xm1_yp1_zm1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i0]];
  const DataT i_xp2_yp1_zm1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1] + deltaZ[i0]];
  const DataT i_xm1_y_zp2 = mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i2]];
  const DataT i_xp2_y_zp2 = mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i2]];
  const DataT i_xm1_yp1_zp2 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i2]];
  const DataT i_xp2_yp1_zp2 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1] + deltaZ[i2]];

  const DataT i_x_ym1_zm1 = mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i0]];
  const DataT i_xp1_ym1_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i0]];
  const DataT i_x_yp2_zm1 = mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i0]];
  const DataT i_xp1_yp2_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2] + deltaZ[i0]];
  const DataT i_x_ym1_zp2 = mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i2]];
  const DataT i_xp1_ym1_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i2]];
  const DataT i_x_yp2_zp2 = mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i2]];
  const DataT i_xp1_yp2_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2] + deltaZ[i2]];

  const DataT i_xm1_ym1_zm1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i0] + deltaZ[i0]];
  const DataT i_xp2_ym1_zm1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0] + deltaZ[i0]];
  const DataT i_xm1_yp2_zm1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2] + deltaZ[i0]];
  const DataT i_xp2_yp2_zm1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i2] + deltaZ[i0]];
  const DataT i_xm1_ym1_zp2 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i0] + deltaZ[i2]];
  const DataT i_xp2_ym1_zp2 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0] + deltaZ[i2]];
  const DataT i_xm1_yp2_zp2 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2] + deltaZ[i2]];
  const DataT i_xp2_yp2_zp2 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i2] + deltaZ[i2]];

  // load values to tmp Vc
  // needed for first derivative
  const Vector<DataT, 24> vecDeriv1A{
    {i_xp1_y_z, i_xp2_y_z, i_xp1_yp1_z, i_xp2_yp1_z, i_xp1_y_zp1, i_xp2_y_zp1, i_xp1_yp1_zp1, i_xp2_yp1_zp1,
     i_x_yp1_z, i_xp1_yp1_z, i_x_yp2_z, i_xp1_yp2_z, i_x_yp1_zp1, i_xp1_yp1_zp1, i_x_yp2_zp1, i_xp1_yp2_zp1,
     i_x_y_zp1, i_xp1_y_zp1, i_x_yp1_zp1, i_xp1_yp1_zp1, i_x_y_zp2, i_xp1_y_zp2, i_x_yp1_zp2, i_xp1_yp1_zp2}};

  const Vector<DataT, 24> vecDeriv1B{
    {i_xm1_y_z, i_x_y_z, i_xm1_yp1_z, i_x_yp1_z, i_xm1_y_zp1, i_x_y_zp1, i_xm1_yp1_zp1, i_x_yp1_zp1,
     i_x_ym1_z, i_xp1_ym1_z, i_x_y_z, i_xp1_y_z, i_x_ym1_zp1, i_xp1_ym1_zp1, i_x_y_zp1, i_xp1_y_zp1,
     i_x_y_zm1, i_xp1_y_zm1, i_x_yp1_zm1, i_xp1_yp1_zm1, i_x_y_z, i_xp1_y_z, i_x_yp1_z, i_xp1_yp1_z}};

  // needed for second derivative
  const Vector<DataT, 24> vecDeriv2A{
    {i_xp1_yp1_z, i_xp2_yp1_z, i_xp1_yp2_z, i_xp2_yp2_z, i_xp1_yp1_zp1, i_xp2_yp1_zp1, i_xp1_yp2_zp1, i_xp2_yp2_zp1,
     i_xp1_y_zp1, i_xp2_y_zp1, i_xp1_yp1_zp1, i_xp2_yp1_zp1, i_xp1_y_zp2, i_xp2_y_zp2, i_xp1_yp1_zp2, i_xp2_yp1_zp2,
     i_x_yp1_zp1, i_xp1_yp1_zp1, i_x_yp2_zp1, i_xp1_yp2_zp1, i_x_yp1_zp2, i_xp1_yp1_zp2, i_x_yp2_zp2, i_xp1_yp2_zp2}};

  const Vector<DataT, 24> vecDeriv2B{
    {i_xm1_yp1_z, i_x_yp1_z, i_xm1_yp2_z, i_x_yp2_z, i_xm1_yp1_zp1, i_x_yp1_zp1, i_xm1_yp2_zp1, i_x_yp2_zp1,
     i_xm1_y_zp1, i_x_y_zp1, i_xm1_yp1_zp1, i_x_yp1_zp1, i_xm1_y_zp2, i_x_y_zp2, i_xm1_yp1_zp2, i_x_yp1_zp2,
     i_x_ym1_zp1, i_xp1_ym1_zp1, i_x_y_zp1, i_xp1_y_zp1, i_x_ym1_zp2, i_xp1_ym1_zp2, i_x_y_zp2, i_xp1_y_zp2}};

  const Vector<DataT, 24> vecDeriv2C{
    {i_xp1_ym1_z, i_xp2_ym1_z, i_xp1_y_z, i_xp2_y_z, i_xp1_ym1_zp1, i_xp2_ym1_zp1, i_xp1_y_zp1, i_xp2_y_zp1,
     i_xp1_y_zm1, i_xp2_y_zm1, i_xp1_yp1_zm1, i_xp2_yp1_zm1, i_xp1_y_z, i_xp2_y_z, i_xp1_yp1_z, i_xp2_yp1_z,
     i_x_yp1_zm1, i_xp1_yp1_zm1, i_x_yp2_zm1, i_xp1_yp2_zm1, i_x_yp1_z, i_xp1_yp1_z, i_x_yp2_z, i_xp1_yp2_z}};

  const Vector<DataT, 24> vecDeriv2D{
    {i_xm1_ym1_z, i_x_ym1_z, i_xm1_y_z, i_x_y_z, i_xm1_ym1_zp1, i_x_ym1_zp1, i_xm1_y_zp1, i_x_y_zp1,
     i_xm1_y_zm1, i_x_y_zm1, i_xm1_yp1_zm1, i_x_yp1_zm1, i_xm1_y_z, i_x_y_z, i_xm1_yp1_z, i_x_yp1_z,
     i_x_ym1_zm1, i_xp1_ym1_zm1, i_x_y_zm1, i_xp1_y_zm1, i_x_ym1_z, i_xp1_ym1_z, i_x_y_z, i_xp1_y_z}};

  // needed for third derivative
  const Vector<DataT, 8> vecDeriv3A{{i_xp1_yp1_zp1, i_xp2_yp1_zp1, i_xp1_yp2_zp1, i_xp2_yp2_zp1, i_xp1_yp1_zp2, i_xp2_yp1_zp2, i_xp1_yp2_zp2, i_xp2_yp2_zp2}};
  const Vector<DataT, 8> vecDeriv3B{{i_xm1_yp1_zp1, i_x_yp1_zp1, i_xm1_yp2_zp1, i_x_yp2_zp1, i_xm1_yp1_zp2, i_x_yp1_zp2, i_xm1_yp2_zp2, i_x_yp2_zp2}};
  const Vector<DataT, 8> vecDeriv3C{{i_xp1_ym1_zp1, i_xp2_ym1_zp1, i_xp1_y_zp1, i_xp2_y_zp1, i_xp1_ym1_zp2, i_xp2_ym1_zp2, i_xp1_y_zp2, i_xp2_y_zp2}};
  const Vector<DataT, 8> vecDeriv3D{{i_xm1_ym1_zp1, i_x_ym1_zp1, i_xm1_y_zp1, i_x_y_zp1, i_xm1_ym1_zp2, i_x_ym1_zp2, i_xm1_y_zp2, i_x_y_zp2}};
  const Vector<DataT, 8> vecDeriv3E{{i_xp1_yp1_zm1, i_xp2_yp1_zm1, i_xp1_yp2_zm1, i_xp2_yp2_zm1, i_xp1_yp1_z, i_xp2_yp1_z, i_xp1_yp2_z, i_xp2_yp2_z}};
  const Vector<DataT, 8> vecDeriv3F{{i_xm1_yp1_zm1, i_x_yp1_zm1, i_xm1_yp2_zm1, i_x_yp2_zm1, i_xm1_yp1_z, i_x_yp1_z, i_xm1_yp2_z, i_x_yp2_z}};
  const Vector<DataT, 8> vecDeriv3G{{i_xp1_ym1_zm1, i_xp2_ym1_zm1, i_xp1_y_zm1, i_xp2_y_zm1, i_xp1_ym1_z, i_xp2_ym1_z, i_xp1_y_z, i_xp2_y_z}};
  const Vector<DataT, 8> vecDeriv3H{{i_xm1_ym1_zm1, i_x_ym1_zm1, i_xm1_y_zm1, i_x_y_zm1, i_xm1_ym1_z, i_x_ym1_z, i_xm1_y_z, i_x_y_z}};

  // factor for first derivative
  const DataT fac1{0.5};
  const Vector<DataT, 24> vfac1{
    {fac1, fac1, fac1, fac1, fac1, fac1, fac1, fac1,
     fac1, fac1, fac1, fac1, fac1, fac1, fac1, fac1,
     fac1, fac1, fac1, fac1, fac1, fac1, fac1, fac1}};

  // factor for second derivative
  const DataT fac2{0.25};
  const Vector<DataT, 24> vfac2{
    {fac2, fac2, fac2, fac2, fac2, fac2, fac2, fac2,
     fac2, fac2, fac2, fac2, fac2, fac2, fac2, fac2,
     fac2, fac2, fac2, fac2, fac2, fac2, fac2, fac2}};

  // factor for third derivative
  const DataT fac3{0.125};
  const Vector<DataT, 8> vfac3{{fac3, fac3, fac3, fac3, fac3, fac3, fac3, fac3}};

  // compute the derivatives
  const Vector<DataT, 24> vecDeriv1Res{vfac1 * (vecDeriv1A - vecDeriv1B)};
  const Vector<DataT, 24> vecDeriv2Res{vfac2 * (vecDeriv2A - vecDeriv2B - vecDeriv2C + vecDeriv2D)};
  const Vector<DataT, 8> vecDeriv3Res{vfac3 * (vecDeriv3A - vecDeriv3B - vecDeriv3C + vecDeriv3D - vecDeriv3E + vecDeriv3F + vecDeriv3G - vecDeriv3H)};

  const Vector<DataT, 64> matrixPar{
    {i_x_y_z, i_xp1_y_z, i_x_yp1_z, i_xp1_yp1_z, i_x_y_zp1, i_xp1_y_zp1, i_x_yp1_zp1, i_xp1_yp1_zp1,
     vecDeriv1Res[0], vecDeriv1Res[1], vecDeriv1Res[2], vecDeriv1Res[3], vecDeriv1Res[4], vecDeriv1Res[5], vecDeriv1Res[6], vecDeriv1Res[7], vecDeriv1Res[8], vecDeriv1Res[9], vecDeriv1Res[10],
     vecDeriv1Res[11], vecDeriv1Res[12], vecDeriv1Res[13], vecDeriv1Res[14], vecDeriv1Res[15], vecDeriv1Res[16], vecDeriv1Res[17], vecDeriv1Res[18], vecDeriv1Res[19], vecDeriv1Res[20], vecDeriv1Res[21],
     vecDeriv1Res[22], vecDeriv1Res[23], vecDeriv2Res[0], vecDeriv2Res[1], vecDeriv2Res[2], vecDeriv2Res[3], vecDeriv2Res[4], vecDeriv2Res[5], vecDeriv2Res[6], vecDeriv2Res[7], vecDeriv2Res[8], vecDeriv2Res[9],
     vecDeriv2Res[10], vecDeriv2Res[11], vecDeriv2Res[12], vecDeriv2Res[13], vecDeriv2Res[14], vecDeriv2Res[15], vecDeriv2Res[16], vecDeriv2Res[17], vecDeriv2Res[18], vecDeriv2Res[19], vecDeriv2Res[20],
     vecDeriv2Res[21], vecDeriv2Res[22], vecDeriv2Res[23], vecDeriv3Res[0], vecDeriv3Res[1], vecDeriv3Res[2], vecDeriv3Res[3], vecDeriv3Res[4], vecDeriv3Res[5], vecDeriv3Res[6], vecDeriv3Res[7]}};

  // calc coeffiecients
  mCoefficients[sThreadnum] = sMatrixA * matrixPar;
}

template <typename DataT, size_t Nx, size_t Ny, size_t Nz>
DataT TriCubicInterpolator<DataT, Nx, Ny, Nz>::extrapolation(const DataT valk, const DataT valk1, const int valk2) const
{
  switch (mExtrapolationType) {
    case ExtrapolationType::Linear:
    default:
      return linearExtrapolation(valk, valk1);
      break;
    case ExtrapolationType::Parabola:
      return parabolExtrapolation(valk, valk1, valk2);
      break;
  }
}

template <typename DataT, size_t Nx, size_t Ny, size_t Nz>
DataT TriCubicInterpolator<DataT, Nx, Ny, Nz>::linearExtrapolation(const DataT valk, const DataT valk1) const
{
  const DataT val = 2 * valk - valk1;
  return val;
}

template <typename DataT, size_t Nx, size_t Ny, size_t Nz>
DataT TriCubicInterpolator<DataT, Nx, Ny, Nz>::parabolExtrapolation(const DataT valk, const DataT valk1, const DataT valk2) const
{
  const DataT val = 3 * (valk - valk1) + valk2; // legendre polynom with x0=0, x1=1, x2=2 and x=-1
  return val;
}

template <typename DataT, size_t Nx, size_t Ny, size_t Nz>
void TriCubicInterpolator<DataT, Nx, Ny, Nz>::calcCoefficientsExtrapolation(const unsigned int ix, const unsigned int iy, const unsigned int iz) const
{
  int deltaX[3]{mGridProperties.getDeltaDataIndex(-1, 0), mGridProperties.getDeltaDataIndex(1, 0), mGridProperties.getDeltaDataIndex(2, 0)};
  int deltaY[3]{mGridProperties.getDeltaDataIndex(-1, 1), mGridProperties.getDeltaDataIndex(1, 1), mGridProperties.getDeltaDataIndex(2, 1)};
  int deltaZ[3]{};
  getDataIndexCircularArray(iz, FZ, deltaZ);
  const int location = findPos(ix, iy, iz);

  // indices to data storage
  const int i0 = 0;
  const int i1 = 1;
  const int i2 = 2;
  const int ii_x_y_z = mGridData.getDataIndex(ix, iy, iz);
  const DataT i_x_y_z = mGridData[ii_x_y_z];
  DataT i_xp1_y_z{};
  DataT i_x_yp1_z{};
  DataT i_xp1_yp1_z{};
  DataT i_x_y_zp1{};
  DataT i_xp1_y_zp1{};
  DataT i_x_yp1_zp1{};
  DataT i_xm1_y_z{};
  DataT i_xm1_yp1_z{};
  DataT i_xm1_y_zp1{};
  DataT i_xm1_yp1_zp1{};
  DataT i_xm1_ym1_z{};
  DataT i_xm1_yp2_z{};
  DataT i_xm1_ym1_zp1{};
  DataT i_xm1_yp2_zp1{};
  DataT i_xm1_y_zm1{};
  DataT i_xm1_yp1_zm1{};
  DataT i_xm1_y_zp2{};
  DataT i_xm1_yp1_zp2{};
  DataT i_xm1_ym1_zm1{};
  DataT i_xm1_yp2_zm1{};
  DataT i_xm1_ym1_zp2{};
  DataT i_xm1_yp2_zp2{};
  DataT i_xp1_yp1_zp1{};
  DataT i_xp2_y_z{};
  DataT i_xp2_yp1_z{};
  DataT i_xp2_y_zp1{};
  DataT i_xp2_yp1_zp1{};
  DataT i_x_ym1_z{};
  DataT i_xp1_ym1_z{};
  DataT i_x_yp2_z{};
  DataT i_xp1_yp2_z{};
  DataT i_x_ym1_zp1{};
  DataT i_xp1_ym1_zp1{};
  DataT i_x_yp2_zp1{};
  DataT i_xp1_yp2_zp1{};
  DataT i_x_y_zm1{};
  DataT i_xp1_y_zm1{};
  DataT i_x_yp1_zm1{};
  DataT i_xp1_yp1_zm1{};
  DataT i_x_y_zp2{};
  DataT i_xp1_y_zp2{};
  DataT i_x_yp1_zp2{};
  DataT i_xp1_yp1_zp2{};
  DataT i_xp2_ym1_z{};
  DataT i_xp2_yp2_z{};
  DataT i_xp2_ym1_zp1{};
  DataT i_xp2_yp2_zp1{};
  DataT i_xp2_y_zm1{};
  DataT i_xp2_yp1_zm1{};
  DataT i_xp2_y_zp2{};
  DataT i_xp2_yp1_zp2{};
  DataT i_x_ym1_zm1{};
  DataT i_xp1_ym1_zm1{};
  DataT i_x_yp2_zm1{};
  DataT i_xp1_yp2_zm1{};
  DataT i_x_ym1_zp2{};
  DataT i_xp1_ym1_zp2{};
  DataT i_x_yp2_zp2{};
  DataT i_xp1_yp2_zp2{};
  DataT i_xp2_ym1_zm1{};
  DataT i_xp2_yp2_zm1{};
  DataT i_xp2_ym1_zp2{};
  DataT i_xp2_yp2_zp2{};

  // DataContainer3D<int, 4, 4 ,4> arrB;

  // /*
  const int ind[4][4][4]{
  {
    {ii_x_y_z + deltaZ[i0] + deltaY[i0] + deltaX[i0],        ind[0][0][0] - deltaX[i0],  ind[0][0][1] - deltaX[i0],  ind[0][0][2] - deltaX[i0]},
    { ind[0][0][0] - deltaY[i0],                              ind[0][1][0] - deltaX[i0],  ind[0][1][1] - deltaX[i0],  ind[0][1][2] - deltaX[i0]},
    { ind[0][1][0] - deltaY[i0],                              ind[0][2][0] - deltaX[i0],  ind[0][2][1] - deltaX[i0],  ind[0][2][2] - deltaX[i0]},
    { ind[0][2][0] - deltaY[i0],                              ind[0][3][0] - deltaX[i0],  ind[0][3][1] - deltaX[i0],  ind[0][3][2] - deltaX[i0]}
  },
  {
    { ind[0][0][0] - deltaZ[i0],                              ind[1][0][0] - deltaX[i0],  ind[1][0][1] - deltaX[i0],  ind[1][0][2] - deltaX[i0]},
    { ind[1][0][0] - deltaY[i0],                              ind[1][1][0] - deltaX[i0],  ind[1][1][1] - deltaX[i0],  ind[1][1][2] - deltaX[i0]},
    { ind[1][1][0] - deltaY[i0],                              ind[1][2][0] - deltaX[i0],  ind[1][2][1] - deltaX[i0],  ind[1][2][2] - deltaX[i0]},
    { ind[1][2][0] - deltaY[i0],                              ind[1][3][0] - deltaX[i0],  ind[1][3][1] - deltaX[i0],  ind[1][3][2] - deltaX[i0]}
  },
  {
    { ind[1][0][0] - deltaZ[i0],                              ind[2][0][0] - deltaX[i0],  ind[2][0][1] - deltaX[i0],  ind[2][0][2] - deltaX[i0]},
    { ind[2][0][0] - deltaY[i0],                              ind[2][1][0] - deltaX[i0],  ind[2][1][1] - deltaX[i0],  ind[2][1][2] - deltaX[i0]},
    { ind[2][1][0] - deltaY[i0],                              ind[2][2][0] - deltaX[i0],  ind[2][2][1] - deltaX[i0],  ind[2][2][2] - deltaX[i0]},
    { ind[2][2][0] - deltaY[i0],                              ind[2][3][0] - deltaX[i0],  ind[2][3][1] - deltaX[i0],  ind[2][3][2] - deltaX[i0]}
  },
  {
      { ind[2][0][0] - deltaZ[i0],                              ind[3][0][0] - deltaX[i0],  ind[3][0][1] - deltaX[i0],  ind[3][0][2] - deltaX[i0]},
      { ind[3][0][0] - deltaY[i0],                              ind[3][1][0] - deltaX[i0],  ind[3][1][1] - deltaX[i0],  ind[3][1][2] - deltaX[i0]},
      { ind[3][1][0] - deltaY[i0],                              ind[3][2][0] - deltaX[i0],  ind[3][2][1] - deltaX[i0],  ind[3][2][2] - deltaX[i0]},
      { ind[3][2][0] - deltaY[i0],                              ind[3][3][0] - deltaX[i0],  ind[3][3][1] - deltaX[i0],  ind[3][3][2] - deltaX[i0]}
    }
 };

//   std::cout<<"-------"<<std::endl;
//   for (int k = 0; k < 4; ++k) {
//     std::cout<<"---k--- "<< k << std::endl;
//   for (int i = 0; i < 4; ++i) {
//      for (int j = 0; j < 4; ++j) {
//          std::cout<< ind[k][i][j] << " ";
//      }
//      std::cout<<std::endl;
//  }
//  std::cout<<std::endl;
//  std::cout<<std::endl;
// }
//  std::cout<<"-------"<<std::endl;
// */

  switch (location) {
    case InnerVolume:
    case SideZRight:
    case SideZLeft:
    default:
      i_xp1_y_z = mGridData[ii_x_y_z + deltaX[i1]];

      // ind[0]
      i_x_yp1_z = mGridData[ii_x_y_z + deltaY[i1]];
      i_xp1_yp1_z = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1]];
      i_x_y_zp1 = mGridData[ii_x_y_z + deltaZ[i1]];
      i_xp1_y_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i1]];
      i_x_yp1_zp1 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i1]];
      i_xm1_y_z = mGridData[ii_x_y_z + deltaX[i0]];
      i_xm1_yp1_z = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1]];
      i_xm1_y_zp1 = mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i1]];
      i_xm1_yp1_zp1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i1]];
      i_xm1_ym1_z = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i0]];
      i_xm1_yp2_z = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2]];
      i_xm1_ym1_zp1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i0] + deltaZ[i1]];
      i_xm1_yp2_zp1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2] + deltaZ[i1]];
      i_xm1_y_zm1 = mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i0]];
      i_xm1_yp1_zm1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i0]];
      i_xm1_y_zp2 = mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i2]];
      i_xm1_yp1_zp2 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i2]];
      i_xm1_ym1_zm1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i0] + deltaZ[i0]];
      i_xm1_yp2_zm1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2] + deltaZ[i0]];
      i_xm1_ym1_zp2 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i0] + deltaZ[i2]];
      i_xm1_yp2_zp2 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2] + deltaZ[i2]];
      i_xp1_yp1_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i1]];
      i_xp2_y_z = mGridData[ii_x_y_z + deltaX[i2]];
      i_xp2_yp1_z = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1]];
      i_xp2_y_zp1 = mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i1]];
      i_xp2_yp1_zp1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1] + deltaZ[i1]];
      i_x_ym1_z = mGridData[ii_x_y_z + deltaY[i0]];
      i_xp1_ym1_z = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0]];
      i_x_yp2_z = mGridData[ii_x_y_z + deltaY[i2]];
      i_xp1_yp2_z = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2]];
      i_x_ym1_zp1 = mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i1]];
      i_xp1_ym1_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i1]];
      i_x_yp2_zp1 = mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i1]];
      i_xp1_yp2_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2] + deltaZ[i1]];
      i_x_y_zm1 = mGridData[ii_x_y_z + deltaZ[i0]];
      i_xp1_y_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i0]];
      i_x_yp1_zm1 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i0]];
      i_xp1_yp1_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i0]];
      i_x_y_zp2 = mGridData[ii_x_y_z + deltaZ[i2]];
      i_xp1_y_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i2]];
      i_x_yp1_zp2 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i2]];
      i_xp1_yp1_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i2]];
      i_xp2_ym1_z = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0]];
      i_xp2_yp2_z = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i2]];
      i_xp2_ym1_zp1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0] + deltaZ[i1]];
      i_xp2_yp2_zp1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i2] + deltaZ[i1]];
      i_xp2_y_zm1 = mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i0]];
      i_xp2_yp1_zm1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1] + deltaZ[i0]];
      i_xp2_y_zp2 = mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i2]];
      i_xp2_yp1_zp2 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1] + deltaZ[i2]];
      i_x_ym1_zm1 = mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i0]];
      i_xp1_ym1_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i0]];
      i_x_yp2_zm1 = mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i0]];
      i_xp1_yp2_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2] + deltaZ[i0]];
      i_x_ym1_zp2 = mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i2]];
      i_xp1_ym1_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i2]];
      i_x_yp2_zp2 = mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i2]];
      i_xp1_yp2_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2] + deltaZ[i2]];
      i_xp2_ym1_zm1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0] + deltaZ[i0]];
      i_xp2_yp2_zm1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i2] + deltaZ[i0]];
      i_xp2_ym1_zp2 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0] + deltaZ[i2]];
      i_xp2_yp2_zp2 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i2] + deltaZ[i2]];
      break;

    case SideXRight:
    case LineD:
    case LineH:
      i_xp1_y_z = mGridData[ii_x_y_z + deltaX[i1]];
      i_x_yp1_z = mGridData[ii_x_y_z + deltaY[i1]];
      i_xp1_yp1_z = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1]];
      i_x_y_zp1 = mGridData[ii_x_y_z + deltaZ[i1]];
      i_xp1_y_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i1]];
      i_x_yp1_zp1 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i1]];
      i_xm1_y_z = mGridData[ii_x_y_z + deltaX[i0]];
      i_xm1_yp1_z = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1]];
      i_xm1_y_zp1 = mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i1]];
      i_xm1_yp1_zp1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i1]];
      i_xm1_ym1_z = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i0]];
      i_xm1_yp2_z = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2]];
      i_xm1_ym1_zp1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i0] + deltaZ[i1]];
      i_xm1_yp2_zp1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2] + deltaZ[i1]];
      i_xm1_y_zm1 = mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i0]];
      i_xm1_yp1_zm1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i0]];
      i_xm1_y_zp2 = mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i2]];
      i_xm1_yp1_zp2 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i2]];
      i_xm1_ym1_zm1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i0] + deltaZ[i0]];
      i_xm1_yp2_zm1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2] + deltaZ[i0]];
      i_xm1_ym1_zp2 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i0] + deltaZ[i2]];
      i_xm1_yp2_zp2 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2] + deltaZ[i2]];
      i_xp1_yp1_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i1]];
      i_x_ym1_z = mGridData[ii_x_y_z + deltaY[i0]];
      i_xp1_ym1_z = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0]];
      i_x_yp2_z = mGridData[ii_x_y_z + deltaY[i2]];
      i_xp1_yp2_z = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2]];
      i_x_ym1_zp1 = mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i1]];
      i_xp1_ym1_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i1]];
      i_x_yp2_zp1 = mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i1]];
      i_xp1_yp2_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2] + deltaZ[i1]];
      i_x_y_zm1 = mGridData[ii_x_y_z + deltaZ[i0]];
      i_xp1_y_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i0]];
      i_x_yp1_zm1 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i0]];
      i_xp1_yp1_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i0]];
      i_x_y_zp2 = mGridData[ii_x_y_z + deltaZ[i2]];
      i_xp1_y_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i2]];
      i_x_yp1_zp2 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i2]];
      i_xp1_yp1_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i2]];
      i_x_ym1_zm1 = mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i0]];
      i_xp1_ym1_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i0]];
      i_x_yp2_zm1 = mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i0]];
      i_xp1_yp2_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2] + deltaZ[i0]];
      i_x_ym1_zp2 = mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i2]];
      i_xp1_ym1_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i2]];
      i_x_yp2_zp2 = mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i2]];
      i_xp1_yp2_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2] + deltaZ[i2]];

      i_xp2_y_z = extrapolation(mGridData[ii_x_y_z + deltaX[i1]], mGridData[ii_x_y_z], mGridData[ii_x_y_z + deltaX[i0]]);
      i_xp2_yp1_z = extrapolation(mGridData[ii_x_y_z + deltaY[i1] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i1]], mGridData[ii_x_y_z + deltaY[i1] + deltaX[i0]]);
      i_xp2_y_zp1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i1] + deltaX[i1]], mGridData[ii_x_y_z + deltaZ[i1]], mGridData[ii_x_y_z + deltaZ[i1] + deltaX[i0]]);
      i_xp2_yp1_zp1 = extrapolation(mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i1] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i1]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i1] + deltaX[i0]]);
      i_xp2_ym1_z = extrapolation(mGridData[ii_x_y_z + deltaY[i0] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i0]], mGridData[ii_x_y_z + deltaY[i0] + deltaX[i0]]);
      i_xp2_yp2_z = extrapolation(mGridData[ii_x_y_z + deltaY[i2] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i2]], mGridData[ii_x_y_z + deltaY[i2] + deltaX[i0]]);
      i_xp2_ym1_zp1 = extrapolation(mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i1] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i1]], mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i1] + deltaX[i0]]);
      i_xp2_yp2_zp1 = extrapolation(mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i1] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i1]], mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i1] + deltaX[i0]]);
      i_xp2_y_zm1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i0] + deltaX[i1]], mGridData[ii_x_y_z + deltaZ[i0]], mGridData[ii_x_y_z + deltaZ[i0] + deltaX[i0]]);
      i_xp2_yp1_zm1 = extrapolation(mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i0] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i0]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i0] + deltaX[i0]]);
      i_xp2_y_zp2 = extrapolation(mGridData[ii_x_y_z + deltaZ[i2] + deltaX[i1]], mGridData[ii_x_y_z + deltaZ[i2]], mGridData[ii_x_y_z + deltaZ[i2] + deltaX[i0]]);
      i_xp2_yp1_zp2 = extrapolation(mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i2] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i2]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i2] + deltaX[i0]]);
      i_xp2_ym1_zm1 = extrapolation(mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i0] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i0]], mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i0] + deltaX[i0]]);
      i_xp2_yp2_zm1 = extrapolation(mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i0] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i0]], mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i0] + deltaX[i0]]);
      i_xp2_ym1_zp2 = extrapolation(mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i2] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i2]], mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i2] + deltaX[i0]]);
      i_xp2_yp2_zp2 = extrapolation(mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i2] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i2]], mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i2] + deltaX[i0]]);
      break;

    case SideYRight:
    case LineB:
    case LineF:
      i_xm1_yp2_z = extrapolation(mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i0]], mGridData[ii_x_y_z + deltaX[i0] + deltaY[i0]]);
      i_xm1_yp2_zp1 = extrapolation(mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i1] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i1] + deltaY[i0]]);
      i_xm1_yp2_zm1 = extrapolation(mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i0] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i0] + deltaY[i0]]);
      i_xm1_yp2_zp2 = extrapolation(mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i2] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i2] + deltaY[i0]]);
      i_x_yp2_z = extrapolation(mGridData[ii_x_y_z + deltaY[i1]], mGridData[ii_x_y_z], mGridData[ii_x_y_z + deltaY[i0]]);
      i_xp1_yp2_z = extrapolation(mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0]]);
      i_x_yp2_zp1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i1] + deltaY[i1]], mGridData[ii_x_y_z + deltaZ[i1]], mGridData[ii_x_y_z + deltaZ[i1] + deltaY[i0]]);
      i_xp1_yp2_zp1 = extrapolation(mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i1] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i1] + deltaY[i0]]);
      i_xp2_yp2_z = extrapolation(mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i2]], mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0]]);
      i_xp2_yp2_zp1 = extrapolation(mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i1] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i1] + deltaY[i0]]);
      i_x_yp2_zm1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i0] + deltaY[i1]], mGridData[ii_x_y_z + deltaZ[i0]], mGridData[ii_x_y_z + deltaZ[i0] + deltaY[i0]]);
      i_xp1_yp2_zm1 = extrapolation(mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i0] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i0] + deltaY[i0]]);
      i_x_yp2_zp2 = extrapolation(mGridData[ii_x_y_z + deltaZ[i2] + deltaY[i1]], mGridData[ii_x_y_z + deltaZ[i2]], mGridData[ii_x_y_z + deltaZ[i2] + deltaY[i0]]);
      i_xp1_yp2_zp2 = extrapolation(mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i2] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i2] + deltaY[i0]]);
      i_xp2_yp2_zm1 = extrapolation(mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i0] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i0] + deltaY[i0]]);
      i_xp2_yp2_zp2 = extrapolation(mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i2] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i2] + deltaY[i0]]);

      i_xp1_y_z = mGridData[ii_x_y_z + deltaX[i1]];
      i_x_yp1_z = mGridData[ii_x_y_z + deltaY[i1]];
      i_xp1_yp1_z = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1]];
      i_x_y_zp1 = mGridData[ii_x_y_z + deltaZ[i1]];
      i_xp1_y_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i1]];
      i_x_yp1_zp1 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i1]];
      i_xm1_y_z = mGridData[ii_x_y_z + deltaX[i0]];
      i_xm1_yp1_z = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1]];
      i_xm1_y_zp1 = mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i1]];
      i_xm1_yp1_zp1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i1]];
      i_xm1_ym1_z = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i0]];
      i_xm1_ym1_zp1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i0] + deltaZ[i1]];
      i_xm1_y_zm1 = mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i0]];
      i_xm1_yp1_zm1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i0]];
      i_xm1_y_zp2 = mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i2]];
      i_xm1_yp1_zp2 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i2]];
      i_xm1_ym1_zm1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i0] + deltaZ[i0]];
      i_xm1_ym1_zp2 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i0] + deltaZ[i2]];
      i_xp1_yp1_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i1]];
      i_xp2_y_z = mGridData[ii_x_y_z + deltaX[i2]];
      i_xp2_yp1_z = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1]];
      i_xp2_y_zp1 = mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i1]];
      i_xp2_yp1_zp1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1] + deltaZ[i1]];
      i_x_ym1_z = mGridData[ii_x_y_z + deltaY[i0]];
      i_xp1_ym1_z = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0]];
      i_x_ym1_zp1 = mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i1]];
      i_xp1_ym1_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i1]];
      i_x_y_zm1 = mGridData[ii_x_y_z + deltaZ[i0]];
      i_xp1_y_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i0]];
      i_x_yp1_zm1 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i0]];
      i_xp1_yp1_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i0]];
      i_x_y_zp2 = mGridData[ii_x_y_z + deltaZ[i2]];
      i_xp1_y_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i2]];
      i_x_yp1_zp2 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i2]];
      i_xp1_yp1_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i2]];
      i_xp2_ym1_z = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0]];
      i_xp2_ym1_zp1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0] + deltaZ[i1]];
      i_xp2_y_zm1 = mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i0]];
      i_xp2_yp1_zm1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1] + deltaZ[i0]];
      i_xp2_y_zp2 = mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i2]];
      i_xp2_yp1_zp2 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1] + deltaZ[i2]];
      i_x_ym1_zm1 = mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i0]];
      i_xp1_ym1_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i0]];
      i_x_ym1_zp2 = mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i2]];
      i_xp1_ym1_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i2]];
      i_xp2_ym1_zm1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0] + deltaZ[i0]];
      i_xp2_ym1_zp2 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0] + deltaZ[i2]];
      break;

    case SideYLeft:
    case LineA:
    case LineE:
      i_xm1_ym1_z = extrapolation(mGridData[ii_x_y_z + deltaX[i0]], mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i0] + 2 * deltaY[i1]]);
      i_xm1_ym1_zp1 = extrapolation(mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i0] + 2 * deltaY[i1] + deltaZ[i1]]);
      i_xm1_ym1_zm1 = extrapolation(mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i0] + 2 * deltaY[i1] + deltaZ[i0]]);
      i_xm1_ym1_zp2 = extrapolation(mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i0] + 2 * deltaY[i1] + deltaZ[i2]]);
      i_x_ym1_z = extrapolation(mGridData[ii_x_y_z], mGridData[ii_x_y_z + deltaY[i1]], mGridData[ii_x_y_z + 2 * deltaY[i1]]);
      i_xp1_ym1_z = extrapolation(mGridData[ii_x_y_z + deltaX[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i1] + 2 * deltaY[i1]]);
      i_x_ym1_zp1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i1]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i1]], mGridData[ii_x_y_z + 2 * deltaY[i1] + deltaZ[i1]]);
      i_xp1_ym1_zp1 = extrapolation(mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i1] + 2 * deltaY[i1] + deltaZ[i1]]);
      i_xp2_ym1_z = extrapolation(mGridData[ii_x_y_z + deltaX[i2]], mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i2] + 2 * deltaY[i1]]);
      i_xp2_ym1_zp1 = extrapolation(mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i2] + 2 * deltaY[i1] + deltaZ[i1]]);
      i_x_ym1_zm1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i0]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i0]], mGridData[ii_x_y_z + 2 * deltaY[i1] + deltaZ[i0]]);
      i_xp1_ym1_zm1 = extrapolation(mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i1] + 2 * deltaY[i1] + deltaZ[i0]]);
      i_x_ym1_zp2 = extrapolation(mGridData[ii_x_y_z + deltaZ[i2]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i2]], mGridData[ii_x_y_z + 2 * deltaY[i1] + deltaZ[i2]]);
      i_xp1_ym1_zp2 = extrapolation(mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i1] + 2 * deltaY[i1] + deltaZ[i2]]);
      i_xp2_ym1_zm1 = extrapolation(mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i2] + 2 * deltaY[i1] + deltaZ[i0]]);
      i_xp2_ym1_zp2 = extrapolation(mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i2] + 2 * deltaY[i1] + deltaZ[i2]]);

      i_xp1_y_z = mGridData[ii_x_y_z + deltaX[i1]];
      i_x_yp1_z = mGridData[ii_x_y_z + deltaY[i1]];
      i_xp1_yp1_z = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1]];
      i_x_y_zp1 = mGridData[ii_x_y_z + deltaZ[i1]];
      i_xp1_y_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i1]];
      i_x_yp1_zp1 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i1]];
      i_xm1_y_z = mGridData[ii_x_y_z + deltaX[i0]];
      i_xm1_yp1_z = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1]];
      i_xm1_y_zp1 = mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i1]];
      i_xm1_yp1_zp1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i1]];
      i_xm1_yp2_z = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2]];
      i_xm1_yp2_zp1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2] + deltaZ[i1]];
      i_xm1_y_zm1 = mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i0]];
      i_xm1_yp1_zm1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i0]];
      i_xm1_y_zp2 = mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i2]];
      i_xm1_yp1_zp2 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i2]];
      i_xm1_yp2_zm1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2] + deltaZ[i0]];
      i_xm1_yp2_zp2 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2] + deltaZ[i2]];
      i_xp1_yp1_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i1]];
      i_xp2_y_z = mGridData[ii_x_y_z + deltaX[i2]];
      i_xp2_yp1_z = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1]];
      i_xp2_y_zp1 = mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i1]];
      i_xp2_yp1_zp1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1] + deltaZ[i1]];
      i_x_yp2_z = mGridData[ii_x_y_z + deltaY[i2]];
      i_xp1_yp2_z = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2]];
      i_x_yp2_zp1 = mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i1]];
      i_xp1_yp2_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2] + deltaZ[i1]];
      i_x_y_zm1 = mGridData[ii_x_y_z + deltaZ[i0]];
      i_xp1_y_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i0]];
      i_x_yp1_zm1 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i0]];
      i_xp1_yp1_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i0]];
      i_x_y_zp2 = mGridData[ii_x_y_z + deltaZ[i2]];
      i_xp1_y_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i2]];
      i_x_yp1_zp2 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i2]];
      i_xp1_yp1_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i2]];
      i_xp2_yp2_z = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i2]];
      i_xp2_yp2_zp1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i2] + deltaZ[i1]];
      i_xp2_y_zm1 = mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i0]];
      i_xp2_yp1_zm1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1] + deltaZ[i0]];
      i_xp2_y_zp2 = mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i2]];
      i_xp2_yp1_zp2 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1] + deltaZ[i2]];
      i_x_yp2_zm1 = mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i0]];
      i_xp1_yp2_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2] + deltaZ[i0]];
      i_x_yp2_zp2 = mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i2]];
      i_xp1_yp2_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2] + deltaZ[i2]];
      i_xp2_yp2_zm1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i2] + deltaZ[i0]];
      i_xp2_yp2_zp2 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i2] + deltaZ[i2]];
      break;

    case SideXLeft:
    case LineC:
    case LineG:
      // case extrapolationType::Parabola:
      i_xm1_y_z = extrapolation(mGridData[ii_x_y_z], mGridData[ii_x_y_z + deltaX[i1]], mGridData[ii_x_y_z + 2 * deltaX[i1]]);
      i_xm1_yp1_z = extrapolation(mGridData[ii_x_y_z + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1]], mGridData[ii_x_y_z + 2 * deltaX[i1] + deltaY[i1]]);
      i_xm1_y_zp1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i1]], mGridData[ii_x_y_z + 2 * deltaX[i1] + deltaZ[i1]]);
      i_xm1_yp1_zp1 = extrapolation(mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i1]], mGridData[ii_x_y_z + 2 * deltaX[i1] + deltaY[i1] + deltaZ[i1]]);
      i_xm1_ym1_z = extrapolation(mGridData[ii_x_y_z + deltaY[i0]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0]], mGridData[ii_x_y_z + 2 * deltaX[i1] + deltaY[i0]]);
      i_xm1_yp2_z = extrapolation(mGridData[ii_x_y_z + deltaY[i2]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2]], mGridData[ii_x_y_z + 2 * deltaX[i1] + deltaY[i2]]);
      i_xm1_ym1_zp1 = extrapolation(mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i1]], mGridData[ii_x_y_z + 2 * deltaX[i1] + deltaY[i0] + deltaZ[i1]]);
      i_xm1_yp2_zp1 = extrapolation(mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2] + deltaZ[i1]], mGridData[ii_x_y_z + 2 * deltaX[i1] + deltaY[i2] + deltaZ[i1]]);
      i_xm1_y_zm1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i0]], mGridData[ii_x_y_z + 2 * deltaX[i1] + deltaZ[i0]]);
      i_xm1_yp1_zm1 = extrapolation(mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i0]], mGridData[ii_x_y_z + 2 * deltaX[i1] + deltaY[i1] + deltaZ[i0]]);
      i_xm1_y_zp2 = extrapolation(mGridData[ii_x_y_z + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i2]], mGridData[ii_x_y_z + 2 * deltaX[i1] + deltaZ[i2]]);
      i_xm1_yp1_zp2 = extrapolation(mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i2]], mGridData[ii_x_y_z + 2 * deltaX[i1] + deltaY[i1] + deltaZ[i2]]);
      i_xm1_ym1_zm1 = extrapolation(mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i0]], mGridData[ii_x_y_z + 2 * deltaX[i1] + deltaY[i0] + deltaZ[i0]]);
      i_xm1_yp2_zm1 = extrapolation(mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2] + deltaZ[i0]], mGridData[ii_x_y_z + 2 * deltaX[i1] + deltaY[i2] + deltaZ[i0]]);
      i_xm1_ym1_zp2 = extrapolation(mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i2]], mGridData[ii_x_y_z + 2 * deltaX[i1] + deltaY[i0] + deltaZ[i2]]);
      i_xm1_yp2_zp2 = extrapolation(mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2] + deltaZ[i2]], mGridData[ii_x_y_z + 2 * deltaX[i1] + deltaY[i2] + deltaZ[i2]]);

      i_xp1_y_z = mGridData[ii_x_y_z + deltaX[i1]];
      i_x_yp1_z = mGridData[ii_x_y_z + deltaY[i1]];
      i_xp1_yp1_z = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1]];
      i_x_y_zp1 = mGridData[ii_x_y_z + deltaZ[i1]];
      i_xp1_y_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i1]];
      i_x_yp1_zp1 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i1]];
      i_xp1_yp1_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i1]];
      i_xp2_y_z = mGridData[ii_x_y_z + deltaX[i2]];
      i_xp2_yp1_z = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1]];
      i_xp2_y_zp1 = mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i1]];
      i_xp2_yp1_zp1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1] + deltaZ[i1]];
      i_x_ym1_z = mGridData[ii_x_y_z + deltaY[i0]];
      i_xp1_ym1_z = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0]];
      i_x_yp2_z = mGridData[ii_x_y_z + deltaY[i2]];
      i_xp1_yp2_z = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2]];
      i_x_ym1_zp1 = mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i1]];
      i_xp1_ym1_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i1]];
      i_x_yp2_zp1 = mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i1]];
      i_xp1_yp2_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2] + deltaZ[i1]];
      i_x_y_zm1 = mGridData[ii_x_y_z + deltaZ[i0]];
      i_xp1_y_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i0]];
      i_x_yp1_zm1 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i0]];
      i_xp1_yp1_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i0]];
      i_x_y_zp2 = mGridData[ii_x_y_z + deltaZ[i2]];
      i_xp1_y_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i2]];
      i_x_yp1_zp2 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i2]];
      i_xp1_yp1_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i2]];
      i_xp2_ym1_z = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0]];
      i_xp2_yp2_z = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i2]];
      i_xp2_ym1_zp1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0] + deltaZ[i1]];
      i_xp2_yp2_zp1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i2] + deltaZ[i1]];
      i_xp2_y_zm1 = mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i0]];
      i_xp2_yp1_zm1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1] + deltaZ[i0]];
      i_xp2_y_zp2 = mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i2]];
      i_xp2_yp1_zp2 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1] + deltaZ[i2]];
      i_x_ym1_zm1 = mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i0]];
      i_xp1_ym1_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i0]];
      i_x_yp2_zm1 = mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i0]];
      i_xp1_yp2_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2] + deltaZ[i0]];
      i_x_ym1_zp2 = mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i2]];
      i_xp1_ym1_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i2]];
      i_x_yp2_zp2 = mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i2]];
      i_xp1_yp2_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2] + deltaZ[i2]];
      i_xp2_ym1_zm1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0] + deltaZ[i0]];
      i_xp2_yp2_zm1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i2] + deltaZ[i0]];
      i_xp2_ym1_zp2 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0] + deltaZ[i2]];
      i_xp2_yp2_zp2 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i2] + deltaZ[i2]];
      break;

    case Edge0:
    case Edge4:
    case LineI:
      // case extrapolationType::Parabola:
      i_xm1_y_z = extrapolation(mGridData[ii_x_y_z], mGridData[ii_x_y_z + deltaX[i1]], mGridData[ii_x_y_z + 2 * deltaX[i1]]);
      i_xm1_yp1_z = extrapolation(mGridData[ii_x_y_z + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1]], mGridData[ii_x_y_z + 2 * deltaX[i1] + deltaY[i1]]);
      i_xm1_y_zp1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i1]], mGridData[ii_x_y_z + 2 * deltaX[i1] + deltaZ[i1]]);
      i_xm1_yp1_zp1 = extrapolation(mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i1]], mGridData[ii_x_y_z + 2 * deltaX[i1] + deltaY[i1] + deltaZ[i1]]);
      i_xm1_yp2_z = extrapolation(mGridData[ii_x_y_z + deltaY[i2]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2]], mGridData[ii_x_y_z + 2 * deltaX[i1] + deltaY[i2]]);
      i_xm1_yp2_zp1 = extrapolation(mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2] + deltaZ[i1]], mGridData[ii_x_y_z + 2 * deltaX[i1] + deltaY[i2] + deltaZ[i1]]);
      i_xm1_y_zp2 = extrapolation(mGridData[ii_x_y_z + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i2]], mGridData[ii_x_y_z + 2 * deltaX[i1] + deltaZ[i2]]);
      i_xm1_yp1_zp2 = extrapolation(mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i2]], mGridData[ii_x_y_z + 2 * deltaX[i1] + deltaY[i1] + deltaZ[i2]]);
      i_xm1_yp2_zp2 = extrapolation(mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2] + deltaZ[i2]], mGridData[ii_x_y_z + 2 * deltaX[i1] + deltaY[i2] + deltaZ[i2]]);
      i_xm1_yp1_zm1 = extrapolation(mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i0]], mGridData[ii_x_y_z + 2 * deltaX[i1] + deltaY[i1] + deltaZ[i0]]);
      i_xm1_yp2_zm1 = extrapolation(mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2] + deltaZ[i0]], mGridData[ii_x_y_z + 2 * deltaX[i1] + deltaY[i2] + deltaZ[i0]]);
      i_xm1_y_zm1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i0]], mGridData[ii_x_y_z + 2 * deltaX[i1] + deltaZ[i0]]);

      // these can be extrapolated directly
      i_x_ym1_z = extrapolation(mGridData[ii_x_y_z], mGridData[ii_x_y_z + deltaY[i1]], mGridData[ii_x_y_z + 2 * deltaY[i1]]);
      i_xp1_ym1_z = extrapolation(mGridData[ii_x_y_z + deltaX[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i1] + 2 * deltaY[i1]]);
      i_x_ym1_zp1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i1]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i1]], mGridData[ii_x_y_z + 2 * deltaY[i1] + deltaZ[i1]]);
      i_xp1_ym1_zp1 = extrapolation(mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i1] + 2 * deltaY[i1] + deltaZ[i1]]);
      i_xp2_ym1_z = extrapolation(mGridData[ii_x_y_z + deltaX[i2]], mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i2] + 2 * deltaY[i1]]);
      i_xp2_ym1_zp1 = extrapolation(mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i2] + 2 * deltaY[i1] + deltaZ[i1]]);
      i_x_ym1_zm1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i0]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i0]], mGridData[ii_x_y_z + 2 * deltaY[i1] + deltaZ[i0]]);
      i_xp1_ym1_zm1 = extrapolation(mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i1] + 2 * deltaY[i1] + deltaZ[i0]]);
      i_x_ym1_zp2 = extrapolation(mGridData[ii_x_y_z + deltaZ[i2]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i2]], mGridData[ii_x_y_z + 2 * deltaY[i1] + deltaZ[i2]]);
      i_xp1_ym1_zp2 = extrapolation(mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i1] + 2 * deltaY[i1] + deltaZ[i2]]);
      i_xp2_ym1_zm1 = extrapolation(mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i2] + 2 * deltaY[i1] + deltaZ[i0]]);
      i_xp2_ym1_zp2 = extrapolation(mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i2] + 2 * deltaY[i1] + deltaZ[i2]]);

      // these need some steps first
      i_xm1_ym1_z = extrapolation(mGridData[ii_x_y_z], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1]], mGridData[ii_x_y_z + 2 * deltaX[i1] + 2 * deltaY[i1]]);
      i_xm1_ym1_zp1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i1]], mGridData[ii_x_y_z + 2 * deltaX[i1] + 2 * deltaY[i1] + deltaZ[i1]]);
      i_xm1_ym1_zm1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i0]], mGridData[ii_x_y_z + 2 * deltaX[i1] + 2 * deltaY[i1] + deltaZ[i0]]);
      i_xm1_ym1_zp2 = extrapolation(mGridData[ii_x_y_z + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i2]], mGridData[ii_x_y_z + 2 * deltaX[i1] + 2 * deltaY[i1] + deltaZ[i2]]);

      i_xp1_y_z = mGridData[ii_x_y_z + deltaX[i1]];
      i_x_yp1_z = mGridData[ii_x_y_z + deltaY[i1]];
      i_xp1_yp1_z = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1]];
      i_x_y_zp1 = mGridData[ii_x_y_z + deltaZ[i1]];
      i_xp1_y_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i1]];
      i_x_yp1_zp1 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i1]];
      i_xp1_yp1_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i1]];
      i_xp2_y_z = mGridData[ii_x_y_z + deltaX[i2]];
      i_xp2_yp1_z = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1]];
      i_xp2_y_zp1 = mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i1]];
      i_xp2_yp1_zp1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1] + deltaZ[i1]];
      i_x_yp2_z = mGridData[ii_x_y_z + deltaY[i2]];
      i_xp1_yp2_z = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2]];
      i_x_yp2_zp1 = mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i1]];
      i_xp1_yp2_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2] + deltaZ[i1]];
      i_x_y_zm1 = mGridData[ii_x_y_z + deltaZ[i0]];
      i_xp1_y_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i0]];
      i_x_yp1_zm1 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i0]];
      i_xp1_yp1_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i0]];
      i_x_y_zp2 = mGridData[ii_x_y_z + deltaZ[i2]];
      i_xp1_y_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i2]];
      i_x_yp1_zp2 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i2]];
      i_xp1_yp1_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i2]];
      i_xp2_yp2_z = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i2]];
      i_xp2_yp2_zp1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i2] + deltaZ[i1]];
      i_xp2_y_zm1 = mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i0]];
      i_xp2_yp1_zm1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1] + deltaZ[i0]];
      i_xp2_y_zp2 = mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i2]];
      i_xp2_yp1_zp2 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1] + deltaZ[i2]];
      i_x_yp2_zm1 = mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i0]];
      i_xp1_yp2_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2] + deltaZ[i0]];
      i_x_yp2_zp2 = mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i2]];
      i_xp1_yp2_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2] + deltaZ[i2]];
      i_xp2_yp2_zm1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i2] + deltaZ[i0]];
      i_xp2_yp2_zp2 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i2] + deltaZ[i2]];
      break;

    case Edge1:
    case Edge5:
    case LineJ:
      // these can be extrapolated directly
      i_xp2_y_z = extrapolation(mGridData[ii_x_y_z + deltaX[i1]], mGridData[ii_x_y_z], mGridData[ii_x_y_z + deltaX[i0]]);
      i_xp2_yp1_z = extrapolation(mGridData[ii_x_y_z + deltaY[i1] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i1]], mGridData[ii_x_y_z + deltaY[i1] + deltaX[i0]]);
      i_xp2_y_zp1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i1] + deltaX[i1]], mGridData[ii_x_y_z + deltaZ[i1]], mGridData[ii_x_y_z + deltaZ[i1] + deltaX[i0]]);
      i_xp2_yp1_zp1 = extrapolation(mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i1] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i1]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i1] + deltaX[i0]]);
      i_xp2_yp2_z = extrapolation(mGridData[ii_x_y_z + deltaY[i2] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i2]], mGridData[ii_x_y_z + deltaY[i2] + deltaX[i0]]);
      i_xp2_yp2_zp1 = extrapolation(mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i1] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i1]], mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i1] + deltaX[i0]]);
      i_xp2_y_zm1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i0] + deltaX[i1]], mGridData[ii_x_y_z + deltaZ[i0]], mGridData[ii_x_y_z + deltaZ[i0] + deltaX[i0]]);
      i_xp2_yp1_zm1 = extrapolation(mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i0] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i0]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i0] + deltaX[i0]]);
      i_xp2_y_zp2 = extrapolation(mGridData[ii_x_y_z + deltaZ[i2] + deltaX[i1]], mGridData[ii_x_y_z + deltaZ[i2]], mGridData[ii_x_y_z + deltaZ[i2] + deltaX[i0]]);
      i_xp2_yp1_zp2 = extrapolation(mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i2] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i2]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i2] + deltaX[i0]]);
      i_xp2_yp2_zm1 = extrapolation(mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i0] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i0]], mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i0] + deltaX[i0]]);
      i_xp2_yp2_zp2 = extrapolation(mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i2] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i2]], mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i2] + deltaX[i0]]);

      // these can be extrapolated directly
      i_xm1_ym1_z = extrapolation(mGridData[ii_x_y_z + deltaX[i0]], mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2]]);
      i_xm1_ym1_zp1 = extrapolation(mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2] + deltaZ[i1]]);
      i_xm1_ym1_zm1 = extrapolation(mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2] + deltaZ[i0]]);
      i_xm1_ym1_zp2 = extrapolation(mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2] + deltaZ[i2]]);
      i_x_ym1_z = extrapolation(mGridData[ii_x_y_z], mGridData[ii_x_y_z + deltaY[i1]], mGridData[ii_x_y_z + deltaY[i2]]);
      i_x_ym1_zp1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i1]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i1]], mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i1]]);
      i_x_ym1_zm1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i0]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i0]], mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i0]]);
      i_x_ym1_zp2 = extrapolation(mGridData[ii_x_y_z + deltaZ[i2]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i2]], mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i2]]);

      // these need some steps first
      i_xp1_ym1_zm1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i0]], mGridData[ii_x_y_z + 2 * deltaX[i0] + deltaY[i2] + deltaZ[i0]]);
      i_xp1_ym1_z = extrapolation(mGridData[ii_x_y_z], mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1]], mGridData[ii_x_y_z + 2 * deltaX[i0] + 2 * deltaY[i1]]);
      i_xp1_ym1_zp2 = extrapolation(mGridData[ii_x_y_z + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i2]], mGridData[ii_x_y_z + 2 * deltaX[i0] + deltaY[i2] + deltaZ[i2]]);
      i_xp1_ym1_zp1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i1]], mGridData[ii_x_y_z + 2 * deltaX[i0] + deltaY[i2] + deltaZ[i1]]);

      i_xp2_ym1_z = extrapolation(mGridData[ii_x_y_z + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2]]);
      i_xp2_ym1_zp1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i1] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2] + deltaZ[i1]]);
      i_xp2_ym1_zm1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i0] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2] + deltaZ[i0]]);
      i_xp2_ym1_zp2 = extrapolation(mGridData[ii_x_y_z + deltaZ[i2] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2] + deltaZ[i2]]);

      i_xp1_y_z = mGridData[ii_x_y_z + deltaX[i1]];
      i_x_yp1_z = mGridData[ii_x_y_z + deltaY[i1]];
      i_xp1_yp1_z = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1]];
      i_x_y_zp1 = mGridData[ii_x_y_z + deltaZ[i1]];
      i_xp1_y_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i1]];
      i_x_yp1_zp1 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i1]];
      i_xm1_y_z = mGridData[ii_x_y_z + deltaX[i0]];
      i_xm1_yp1_z = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1]];
      i_xm1_y_zp1 = mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i1]];
      i_xm1_yp1_zp1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i1]];
      i_xm1_yp2_z = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2]];
      i_xm1_yp2_zp1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2] + deltaZ[i1]];
      i_xm1_y_zm1 = mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i0]];
      i_xm1_yp1_zm1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i0]];
      i_xm1_y_zp2 = mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i2]];
      i_xm1_yp1_zp2 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i2]];
      i_xm1_yp2_zm1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2] + deltaZ[i0]];
      i_xm1_yp2_zp2 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i2] + deltaZ[i2]];
      i_xp1_yp1_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i1]];
      i_x_yp2_z = mGridData[ii_x_y_z + deltaY[i2]];
      i_xp1_yp2_z = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2]];
      i_x_yp2_zp1 = mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i1]];
      i_xp1_yp2_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2] + deltaZ[i1]];
      i_x_y_zm1 = mGridData[ii_x_y_z + deltaZ[i0]];
      i_xp1_y_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i0]];
      i_x_yp1_zm1 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i0]];
      i_xp1_yp1_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i0]];
      i_x_y_zp2 = mGridData[ii_x_y_z + deltaZ[i2]];
      i_xp1_y_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i2]];
      i_x_yp1_zp2 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i2]];
      i_xp1_yp1_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i2]];
      i_x_yp2_zm1 = mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i0]];
      i_xp1_yp2_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2] + deltaZ[i0]];
      i_x_yp2_zp2 = mGridData[ii_x_y_z + deltaY[i2] + deltaZ[i2]];
      i_xp1_yp2_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i2] + deltaZ[i2]];
      break;

    case Edge2:
    case Edge6:
    case LineK:
      i_xm1_y_z = extrapolation(mGridData[ii_x_y_z], mGridData[ii_x_y_z + deltaX[i1]], mGridData[ii_x_y_z + deltaX[i2]]);
      i_xm1_y_zp1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i1]]);
      i_xm1_ym1_z = extrapolation(mGridData[ii_x_y_z + deltaY[i0]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0]], mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0]]);
      i_xm1_ym1_zp1 = extrapolation(mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0] + deltaZ[i1]]);
      i_xm1_y_zm1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i0]]);
      i_xm1_y_zp2 = extrapolation(mGridData[ii_x_y_z + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i2]]);
      i_xm1_ym1_zm1 = extrapolation(mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0] + deltaZ[i0]]);
      i_xm1_ym1_zp2 = extrapolation(mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0] + deltaZ[i2]]);

      i_x_yp2_z = extrapolation(mGridData[ii_x_y_z + deltaY[i1]], mGridData[ii_x_y_z], mGridData[ii_x_y_z + deltaY[i0]]);
      i_xp1_yp2_z = extrapolation(mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0]]);
      i_x_yp2_zp1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i1] + deltaY[i1]], mGridData[ii_x_y_z + deltaZ[i1]], mGridData[ii_x_y_z + deltaZ[i1] + deltaY[i0]]);
      i_xp1_yp2_zp1 = extrapolation(mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i1] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i1] + deltaY[i0]]);
      i_xp2_yp2_z = extrapolation(mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i2]], mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0]]);
      i_xp2_yp2_zp1 = extrapolation(mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i1] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i1] + deltaY[i0]]);
      i_x_yp2_zm1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i0] + deltaY[i1]], mGridData[ii_x_y_z + deltaZ[i0]], mGridData[ii_x_y_z + deltaZ[i0] + deltaY[i0]]);
      i_xp1_yp2_zm1 = extrapolation(mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i0] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i0] + deltaY[i0]]);
      i_x_yp2_zp2 = extrapolation(mGridData[ii_x_y_z + deltaZ[i2] + deltaY[i1]], mGridData[ii_x_y_z + deltaZ[i2]], mGridData[ii_x_y_z + deltaZ[i2] + deltaY[i0]]);
      i_xp1_yp2_zp2 = extrapolation(mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i2] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i2] + deltaY[i0]]);
      i_xp2_yp2_zm1 = extrapolation(mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i0] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i0] + deltaY[i0]]);
      i_xp2_yp2_zp2 = extrapolation(mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i2] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i2] + deltaY[i0]]);

      i_xm1_yp1_z = extrapolation(mGridData[ii_x_y_z], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0]], mGridData[ii_x_y_z + deltaX[i2] + 2 * deltaY[i0]]);
      i_xm1_yp1_zp1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i2] + 2 * deltaY[i0] + deltaZ[i1]]);
      i_xm1_yp1_zm1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i2] + 2 * deltaY[i0] + deltaZ[i0]]);
      i_xm1_yp1_zp2 = extrapolation(mGridData[ii_x_y_z + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i2] + 2 * deltaY[i0] + deltaZ[i2]]);

      i_xm1_yp2_zp1 = extrapolation(mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0] + deltaZ[i1]]);
      i_xm1_yp2_zm1 = extrapolation(mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0] + deltaZ[i0]]);
      i_xm1_yp2_zp2 = extrapolation(mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0] + deltaZ[i2]]);
      i_xm1_yp2_z = extrapolation(mGridData[ii_x_y_z + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i0] + deltaX[i2]]);

      i_xp1_y_z = mGridData[ii_x_y_z + deltaX[i1]];
      i_x_yp1_z = mGridData[ii_x_y_z + deltaY[i1]];
      i_xp1_yp1_z = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1]];
      i_x_y_zp1 = mGridData[ii_x_y_z + deltaZ[i1]];
      i_xp1_y_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i1]];
      i_x_yp1_zp1 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i1]];
      i_xp1_yp1_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i1]];
      i_xp2_y_z = mGridData[ii_x_y_z + deltaX[i2]];
      i_xp2_yp1_z = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1]];
      i_xp2_y_zp1 = mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i1]];
      i_xp2_yp1_zp1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1] + deltaZ[i1]];
      i_x_ym1_z = mGridData[ii_x_y_z + deltaY[i0]];
      i_xp1_ym1_z = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0]];
      i_x_ym1_zp1 = mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i1]];
      i_xp1_ym1_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i1]];
      i_x_y_zm1 = mGridData[ii_x_y_z + deltaZ[i0]];
      i_xp1_y_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i0]];
      i_x_yp1_zm1 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i0]];
      i_xp1_yp1_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i0]];
      i_x_y_zp2 = mGridData[ii_x_y_z + deltaZ[i2]];
      i_xp1_y_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i2]];
      i_x_yp1_zp2 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i2]];
      i_xp1_yp1_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i2]];
      i_xp2_ym1_z = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0]];
      i_xp2_ym1_zp1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0] + deltaZ[i1]];
      i_xp2_y_zm1 = mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i0]];
      i_xp2_yp1_zm1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1] + deltaZ[i0]];
      i_xp2_y_zp2 = mGridData[ii_x_y_z + deltaX[i2] + deltaZ[i2]];
      i_xp2_yp1_zp2 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i1] + deltaZ[i2]];
      i_x_ym1_zm1 = mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i0]];
      i_xp1_ym1_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i0]];
      i_x_ym1_zp2 = mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i2]];
      i_xp1_ym1_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i2]];
      i_xp2_ym1_zm1 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0] + deltaZ[i0]];
      i_xp2_ym1_zp2 = mGridData[ii_x_y_z + deltaX[i2] + deltaY[i0] + deltaZ[i2]];
      break;

    case Edge3:
    case Edge7:
    case LineL:
      i_xm1_yp2_z = extrapolation(mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i0]], mGridData[ii_x_y_z + deltaX[i0] + deltaY[i0]]);
      i_xm1_yp2_zp1 = extrapolation(mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i1] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i1] + deltaY[i0]]);
      i_xm1_yp2_zm1 = extrapolation(mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i0] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i0] + deltaY[i0]]);
      i_xm1_yp2_zp2 = extrapolation(mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i2] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i2] + deltaY[i0]]);
      i_x_yp2_z = extrapolation(mGridData[ii_x_y_z + deltaY[i1]], mGridData[ii_x_y_z], mGridData[ii_x_y_z + deltaY[i0]]);
      i_xp1_yp2_z = extrapolation(mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0]]);
      i_x_yp2_zp1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i1] + deltaY[i1]], mGridData[ii_x_y_z + deltaZ[i1]], mGridData[ii_x_y_z + deltaZ[i1] + deltaY[i0]]);
      i_xp1_yp2_zp1 = extrapolation(mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i1] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i1] + deltaY[i0]]);
      i_x_yp2_zm1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i0] + deltaY[i1]], mGridData[ii_x_y_z + deltaZ[i0]], mGridData[ii_x_y_z + deltaZ[i0] + deltaY[i0]]);
      i_xp1_yp2_zm1 = extrapolation(mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i0] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i0]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i0] + deltaY[i0]]);
      i_x_yp2_zp2 = extrapolation(mGridData[ii_x_y_z + deltaZ[i2] + deltaY[i1]], mGridData[ii_x_y_z + deltaZ[i2]], mGridData[ii_x_y_z + deltaZ[i2] + deltaY[i0]]);
      i_xp1_yp2_zp2 = extrapolation(mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i2] + deltaY[i1]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i2]], mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i2] + deltaY[i0]]);

      i_xp2_y_z = extrapolation(mGridData[ii_x_y_z + deltaX[i1]], mGridData[ii_x_y_z], mGridData[ii_x_y_z + deltaX[i0]]);
      i_xp2_yp1_z = extrapolation(mGridData[ii_x_y_z + deltaY[i1] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i1]], mGridData[ii_x_y_z + deltaY[i1] + deltaX[i0]]);
      i_xp2_y_zp1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i1] + deltaX[i1]], mGridData[ii_x_y_z + deltaZ[i1]], mGridData[ii_x_y_z + deltaZ[i1] + deltaX[i0]]);
      i_xp2_yp1_zp1 = extrapolation(mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i1] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i1]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i1] + deltaX[i0]]);
      i_xp2_ym1_z = extrapolation(mGridData[ii_x_y_z + deltaY[i0] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i0]], mGridData[ii_x_y_z + deltaY[i0] + deltaX[i0]]);
      i_xp2_ym1_zp1 = extrapolation(mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i1] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i1]], mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i1] + deltaX[i0]]);
      i_xp2_y_zm1 = extrapolation(mGridData[ii_x_y_z + deltaZ[i0] + deltaX[i1]], mGridData[ii_x_y_z + deltaZ[i0]], mGridData[ii_x_y_z + deltaZ[i0] + deltaX[i0]]);
      i_xp2_yp1_zm1 = extrapolation(mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i0] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i0]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i0] + deltaX[i0]]);
      i_xp2_y_zp2 = extrapolation(mGridData[ii_x_y_z + deltaZ[i2] + deltaX[i1]], mGridData[ii_x_y_z + deltaZ[i2]], mGridData[ii_x_y_z + deltaZ[i2] + deltaX[i0]]);
      i_xp2_yp1_zp2 = extrapolation(mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i2] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i2]], mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i2] + deltaX[i0]]);
      i_xp2_ym1_zm1 = extrapolation(mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i0] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i0]], mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i0] + deltaX[i0]]);
      i_xp2_ym1_zp2 = extrapolation(mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i2] + deltaX[i1]], mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i2]], mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i2] + deltaX[i0]]);

      i_xp2_yp2_zp2 = extrapolation(mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i2] + deltaX[i1]], mGridData[ii_x_y_z + deltaZ[i2]], mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i2] + deltaX[i0]]);
      i_xp2_yp2_zp1 = extrapolation(mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i1] + deltaX[i1]], mGridData[ii_x_y_z + deltaZ[i1]], mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i1] + deltaX[i0]]);
      i_xp2_yp2_zm1 = extrapolation(mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i0] + deltaX[i1]], mGridData[ii_x_y_z + deltaZ[i0]], mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i0] + deltaX[i0]]);
      i_xp2_yp2_z = extrapolation(mGridData[ii_x_y_z + deltaY[i1] + deltaX[i1]], mGridData[ii_x_y_z], mGridData[ii_x_y_z + deltaY[i0] + deltaX[i0]]);

      i_xp1_y_z = mGridData[ii_x_y_z + deltaX[i1]];
      i_x_yp1_z = mGridData[ii_x_y_z + deltaY[i1]];
      i_xp1_yp1_z = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1]];
      i_x_y_zp1 = mGridData[ii_x_y_z + deltaZ[i1]];
      i_xp1_y_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i1]];
      i_x_yp1_zp1 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i1]];
      i_xm1_y_z = mGridData[ii_x_y_z + deltaX[i0]];
      i_xm1_yp1_z = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1]];
      i_xm1_y_zp1 = mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i1]];
      i_xm1_yp1_zp1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i1]];
      i_xm1_ym1_z = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i0]];
      i_xm1_ym1_zp1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i0] + deltaZ[i1]];
      i_xm1_y_zm1 = mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i0]];
      i_xm1_yp1_zm1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i0]];
      i_xm1_y_zp2 = mGridData[ii_x_y_z + deltaX[i0] + deltaZ[i2]];
      i_xm1_yp1_zp2 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i1] + deltaZ[i2]];
      i_xm1_ym1_zm1 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i0] + deltaZ[i0]];
      i_xm1_ym1_zp2 = mGridData[ii_x_y_z + deltaX[i0] + deltaY[i0] + deltaZ[i2]];
      i_xp1_yp1_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i1]];
      i_x_ym1_z = mGridData[ii_x_y_z + deltaY[i0]];
      i_xp1_ym1_z = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0]];
      i_x_ym1_zp1 = mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i1]];
      i_xp1_ym1_zp1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i1]];
      i_x_y_zm1 = mGridData[ii_x_y_z + deltaZ[i0]];
      i_xp1_y_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i0]];
      i_x_yp1_zm1 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i0]];
      i_xp1_yp1_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i0]];
      i_x_y_zp2 = mGridData[ii_x_y_z + deltaZ[i2]];
      i_xp1_y_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaZ[i2]];
      i_x_yp1_zp2 = mGridData[ii_x_y_z + deltaY[i1] + deltaZ[i2]];
      i_xp1_yp1_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i1] + deltaZ[i2]];
      i_x_ym1_zm1 = mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i0]];
      i_xp1_ym1_zm1 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i0]];
      i_x_ym1_zp2 = mGridData[ii_x_y_z + deltaY[i0] + deltaZ[i2]];
      i_xp1_ym1_zp2 = mGridData[ii_x_y_z + deltaX[i1] + deltaY[i0] + deltaZ[i2]];
      break;
  }

  /*
std::cout<< "i_xp1_y_z: " << i_xp1_y_z<<std::endl;
std::cout<< "i_x_yp1_z: " << i_x_yp1_z<<std::endl;
std::cout<< "i_xp1_yp1_z: " <<i_xp1_yp1_z<<std::endl;
std::cout<< "i_x_y_zp1: " <<i_x_y_zp1<<std::endl;
std::cout<< "i_xp1_y_zp1: " <<i_xp1_y_zp1<<std::endl;
std::cout<< "i_x_yp1_zp1: " <<i_x_yp1_zp1<<std::endl;
std::cout<< "i_xm1_y_z: " <<i_xm1_y_z<<std::endl;
std::cout<< "i_xm1_yp1_z: " <<i_xm1_yp1_z<<std::endl;
std::cout<< "i_xm1_y_zp1: " <<i_xm1_y_zp1<<std::endl;
std::cout<< "i_xm1_yp1_zp1: " <<i_xm1_yp1_zp1<<std::endl;
std::cout<< "i_xm1_ym1_z: " <<i_xm1_ym1_z<<std::endl;
std::cout<< "i_xm1_yp2_z: " <<i_xm1_yp2_z<<std::endl;
std::cout<< "i_xm1_ym1_zp1: " <<i_xm1_ym1_zp1<<std::endl;
std::cout<< "i_xm1_yp2_zp1: " <<i_xm1_yp2_zp1<<std::endl;
std::cout<< "i_xm1_y_zm1: " <<i_xm1_y_zm1<<std::endl;
std::cout<< "i_xm1_yp1_zm1: " <<i_xm1_yp1_zm1 << std::endl;
std::cout<< "i_xm1_y_zp2: " <<i_xm1_y_zp2 << std::endl;
std::cout<< "i_xm1_yp1_zp2: " <<i_xm1_yp1_zp2 << std::endl;
std::cout<< "i_xm1_ym1_zm1: " <<i_xm1_ym1_zm1 << std::endl;
std::cout<< "i_xm1_yp2_zm1: " <<i_xm1_yp2_zm1 << std::endl;
std::cout<< "i_xm1_ym1_zp2: " <<i_xm1_ym1_zp2 << std::endl;
std::cout<< "i_xm1_yp2_zp2: " <<i_xm1_yp2_zp2 << std::endl;
std::cout<< "i_xp1_yp1_zp1: " <<i_xp1_yp1_zp1 << std::endl;
std::cout<< "i_xp2_y_z: " <<i_xp2_y_z << std::endl;
std::cout<< "i_xp2_yp1_z: " <<i_xp2_yp1_z << std::endl;
std::cout<< "i_xp2_y_zp1: " <<i_xp2_y_zp1 << std::endl;
std::cout<< "i_xp2_yp1_zp1: " <<i_xp2_yp1_zp1 << std::endl;
std::cout<< "i_x_ym1_z: " <<i_x_ym1_z << std::endl;
std::cout<< "i_xp1_ym1_z: " <<i_xp1_ym1_z << std::endl;
std::cout<< "i_x_yp2_z: " <<i_x_yp2_z << std::endl;
std::cout<< "i_xp1_yp2_z: " <<i_xp1_yp2_z << std::endl;
std::cout<< "i_x_ym1_zp1: " <<i_x_ym1_zp1 << std::endl;
std::cout<< "i_xp1_ym1_zp1: " <<i_xp1_ym1_zp1 << std::endl;
std::cout<< "i_x_yp2_zp1: " <<i_x_yp2_zp1 << std::endl;
std::cout<< "i_xp1_yp2_zp1: " <<i_xp1_yp2_zp1 << std::endl;
std::cout<< "i_x_y_zm1: " <<i_x_y_zm1 << std::endl;
std::cout<< "i_xp1_y_zm1: " <<i_xp1_y_zm1 << std::endl;
std::cout<< "i_x_yp1_zm1: " <<i_x_yp1_zm1 << std::endl;
std::cout<< "i_xp1_yp1_zm1: " <<i_xp1_yp1_zm1 << std::endl;
std::cout<< "i_x_y_zp2: " <<i_x_y_zp2 << std::endl;
std::cout<< "i_xp1_y_zp2: " <<i_xp1_y_zp2 << std::endl;
std::cout<< "i_x_yp1_zp2: " <<i_x_yp1_zp2 << std::endl;
std::cout<< "i_xp1_yp1_zp2: " <<i_xp1_yp1_zp2 << std::endl;
std::cout<< "i_xp2_ym1_z: "<< i_xp2_ym1_z << std::endl;
std::cout<< "i_xp2_yp2_z: " <<i_xp2_yp2_z << std::endl;
std::cout<< "i_xp2_ym1_zp1: " <<i_xp2_ym1_zp1 << std::endl;
std::cout<< "i_xp2_yp2_zp1: " <<i_xp2_yp2_zp1 << std::endl;
std::cout<< "i_xp2_y_zm1: " <<i_xp2_y_zm1 << std::endl;
std::cout<< "i_xp2_yp1_zm1: " <<i_xp2_yp1_zm1 << std::endl;
std::cout<< "i_xp2_y_zp2: " <<i_xp2_y_zp2 << std::endl;
std::cout<< "i_xp2_yp1_zp2: " <<i_xp2_yp1_zp2 << std::endl;
std::cout<< "i_x_ym1_zm1: " <<i_x_ym1_zm1 << std::endl;
std::cout<< "i_xp1_ym1_zm1: " <<i_xp1_ym1_zm1 << std::endl;
std::cout<< "i_x_yp2_zm1: " <<i_x_yp2_zm1 << std::endl;
std::cout<< "i_xp1_yp2_zm1: " <<i_xp1_yp2_zm1 << std::endl;
std::cout<< "i_x_ym1_zp2: " <<i_x_ym1_zp2 << std::endl;
std::cout<< "i_xp1_ym1_zp2: " << i_xp1_ym1_zp2 << std::endl;
std::cout<< "i_x_yp2_zp2: "<<i_x_yp2_zp2 << std::endl;
std::cout<< "i_xp1_yp2_zp2: "<<i_xp1_yp2_zp2 << std::endl;
std::cout<< "i_xp2_ym1_zm1: "<<i_xp2_ym1_zm1 << std::endl;
std::cout<< "i_xp2_yp2_zm1: "<<i_xp2_yp2_zm1 << std::endl;
std::cout<< "i_xp2_ym1_zp2: "<<i_xp2_ym1_zp2 << std::endl;
std::cout<< "i_xp2_yp2_zp2: "<<i_xp2_yp2_zp2 << std::endl;
*/

  // load values to tmp Vc
  // needed for first derivative
  const Vector<DataT, 24> vecDeriv1A{
    {i_xp1_y_z, i_xp2_y_z, i_xp1_yp1_z, i_xp2_yp1_z, i_xp1_y_zp1, i_xp2_y_zp1, i_xp1_yp1_zp1, i_xp2_yp1_zp1,
     i_x_yp1_z, i_xp1_yp1_z, i_x_yp2_z, i_xp1_yp2_z, i_x_yp1_zp1, i_xp1_yp1_zp1, i_x_yp2_zp1, i_xp1_yp2_zp1,
     i_x_y_zp1, i_xp1_y_zp1, i_x_yp1_zp1, i_xp1_yp1_zp1, i_x_y_zp2, i_xp1_y_zp2, i_x_yp1_zp2, i_xp1_yp1_zp2}};

  const Vector<DataT, 24> vecDeriv1B{
    {i_xm1_y_z, i_x_y_z, i_xm1_yp1_z, i_x_yp1_z, i_xm1_y_zp1, i_x_y_zp1, i_xm1_yp1_zp1, i_x_yp1_zp1,
     i_x_ym1_z, i_xp1_ym1_z, i_x_y_z, i_xp1_y_z, i_x_ym1_zp1, i_xp1_ym1_zp1, i_x_y_zp1, i_xp1_y_zp1,
     i_x_y_zm1, i_xp1_y_zm1, i_x_yp1_zm1, i_xp1_yp1_zm1, i_x_y_z, i_xp1_y_z, i_x_yp1_z, i_xp1_yp1_z}};

  // needed for second derivative
  const Vector<DataT, 24> vecDeriv2A{
    {i_xp1_yp1_z, i_xp2_yp1_z, i_xp1_yp2_z, i_xp2_yp2_z, i_xp1_yp1_zp1, i_xp2_yp1_zp1, i_xp1_yp2_zp1, i_xp2_yp2_zp1,
     i_xp1_y_zp1, i_xp2_y_zp1, i_xp1_yp1_zp1, i_xp2_yp1_zp1, i_xp1_y_zp2, i_xp2_y_zp2, i_xp1_yp1_zp2, i_xp2_yp1_zp2,
     i_x_yp1_zp1, i_xp1_yp1_zp1, i_x_yp2_zp1, i_xp1_yp2_zp1, i_x_yp1_zp2, i_xp1_yp1_zp2, i_x_yp2_zp2, i_xp1_yp2_zp2}};

  const Vector<DataT, 24> vecDeriv2B{
    {i_xm1_yp1_z, i_x_yp1_z, i_xm1_yp2_z, i_x_yp2_z, i_xm1_yp1_zp1, i_x_yp1_zp1, i_xm1_yp2_zp1, i_x_yp2_zp1,
     i_xm1_y_zp1, i_x_y_zp1, i_xm1_yp1_zp1, i_x_yp1_zp1, i_xm1_y_zp2, i_x_y_zp2, i_xm1_yp1_zp2, i_x_yp1_zp2,
     i_x_ym1_zp1, i_xp1_ym1_zp1, i_x_y_zp1, i_xp1_y_zp1, i_x_ym1_zp2, i_xp1_ym1_zp2, i_x_y_zp2, i_xp1_y_zp2}};

  const Vector<DataT, 24> vecDeriv2C{
    {i_xp1_ym1_z, i_xp2_ym1_z, i_xp1_y_z, i_xp2_y_z, i_xp1_ym1_zp1, i_xp2_ym1_zp1, i_xp1_y_zp1, i_xp2_y_zp1,
     i_xp1_y_zm1, i_xp2_y_zm1, i_xp1_yp1_zm1, i_xp2_yp1_zm1, i_xp1_y_z, i_xp2_y_z, i_xp1_yp1_z, i_xp2_yp1_z,
     i_x_yp1_zm1, i_xp1_yp1_zm1, i_x_yp2_zm1, i_xp1_yp2_zm1, i_x_yp1_z, i_xp1_yp1_z, i_x_yp2_z, i_xp1_yp2_z}};

  const Vector<DataT, 24> vecDeriv2D{
    {i_xm1_ym1_z, i_x_ym1_z, i_xm1_y_z, i_x_y_z, i_xm1_ym1_zp1, i_x_ym1_zp1, i_xm1_y_zp1, i_x_y_zp1,
     i_xm1_y_zm1, i_x_y_zm1, i_xm1_yp1_zm1, i_x_yp1_zm1, i_xm1_y_z, i_x_y_z, i_xm1_yp1_z, i_x_yp1_z,
     i_x_ym1_zm1, i_xp1_ym1_zm1, i_x_y_zm1, i_xp1_y_zm1, i_x_ym1_z, i_xp1_ym1_z, i_x_y_z, i_xp1_y_z}};

  // needed for third derivative
  const Vector<DataT, 8> vecDeriv3A{{i_xp1_yp1_zp1, i_xp2_yp1_zp1, i_xp1_yp2_zp1, i_xp2_yp2_zp1, i_xp1_yp1_zp2, i_xp2_yp1_zp2, i_xp1_yp2_zp2, i_xp2_yp2_zp2}};
  const Vector<DataT, 8> vecDeriv3B{{i_xm1_yp1_zp1, i_x_yp1_zp1, i_xm1_yp2_zp1, i_x_yp2_zp1, i_xm1_yp1_zp2, i_x_yp1_zp2, i_xm1_yp2_zp2, i_x_yp2_zp2}};
  const Vector<DataT, 8> vecDeriv3C{{i_xp1_ym1_zp1, i_xp2_ym1_zp1, i_xp1_y_zp1, i_xp2_y_zp1, i_xp1_ym1_zp2, i_xp2_ym1_zp2, i_xp1_y_zp2, i_xp2_y_zp2}};
  const Vector<DataT, 8> vecDeriv3D{{i_xm1_ym1_zp1, i_x_ym1_zp1, i_xm1_y_zp1, i_x_y_zp1, i_xm1_ym1_zp2, i_x_ym1_zp2, i_xm1_y_zp2, i_x_y_zp2}};
  const Vector<DataT, 8> vecDeriv3E{{i_xp1_yp1_zm1, i_xp2_yp1_zm1, i_xp1_yp2_zm1, i_xp2_yp2_zm1, i_xp1_yp1_z, i_xp2_yp1_z, i_xp1_yp2_z, i_xp2_yp2_z}};
  const Vector<DataT, 8> vecDeriv3F{{i_xm1_yp1_zm1, i_x_yp1_zm1, i_xm1_yp2_zm1, i_x_yp2_zm1, i_xm1_yp1_z, i_x_yp1_z, i_xm1_yp2_z, i_x_yp2_z}};
  const Vector<DataT, 8> vecDeriv3G{{i_xp1_ym1_zm1, i_xp2_ym1_zm1, i_xp1_y_zm1, i_xp2_y_zm1, i_xp1_ym1_z, i_xp2_ym1_z, i_xp1_y_z, i_xp2_y_z}};
  const Vector<DataT, 8> vecDeriv3H{{i_xm1_ym1_zm1, i_x_ym1_zm1, i_xm1_y_zm1, i_x_y_zm1, i_xm1_ym1_z, i_x_ym1_z, i_xm1_y_z, i_x_y_z}};

  // factor for first derivative
  const DataT fac1{0.5};
  const Vector<DataT, 24> vfac1{
    {fac1, fac1, fac1, fac1, fac1, fac1, fac1, fac1,
     fac1, fac1, fac1, fac1, fac1, fac1, fac1, fac1,
     fac1, fac1, fac1, fac1, fac1, fac1, fac1, fac1}};

  // factor for second derivative
  const DataT fac2{0.25};
  const Vector<DataT, 24> vfac2{
    {fac2, fac2, fac2, fac2, fac2, fac2, fac2, fac2,
     fac2, fac2, fac2, fac2, fac2, fac2, fac2, fac2,
     fac2, fac2, fac2, fac2, fac2, fac2, fac2, fac2}};

  // factor for third derivative
  const DataT fac3{0.125};
  const Vector<DataT, 8> vfac3{{fac3, fac3, fac3, fac3, fac3, fac3, fac3, fac3}};

  // compute the derivatives
  const Vector<DataT, 24> vecDeriv1Res{vfac1 * (vecDeriv1A - vecDeriv1B)};
  const Vector<DataT, 24> vecDeriv2Res{vfac2 * (vecDeriv2A - vecDeriv2B - vecDeriv2C + vecDeriv2D)};
  const Vector<DataT, 8> vecDeriv3Res{vfac3 * (vecDeriv3A - vecDeriv3B - vecDeriv3C + vecDeriv3D - vecDeriv3E + vecDeriv3F + vecDeriv3G - vecDeriv3H)};

  const Vector<DataT, 64> matrixPar{
    {i_x_y_z, i_xp1_y_z, i_x_yp1_z, i_xp1_yp1_z, i_x_y_zp1, i_xp1_y_zp1, i_x_yp1_zp1, i_xp1_yp1_zp1,
     vecDeriv1Res[0], vecDeriv1Res[1], vecDeriv1Res[2], vecDeriv1Res[3], vecDeriv1Res[4], vecDeriv1Res[5], vecDeriv1Res[6], vecDeriv1Res[7], vecDeriv1Res[8], vecDeriv1Res[9], vecDeriv1Res[10],
     vecDeriv1Res[11], vecDeriv1Res[12], vecDeriv1Res[13], vecDeriv1Res[14], vecDeriv1Res[15], vecDeriv1Res[16], vecDeriv1Res[17], vecDeriv1Res[18], vecDeriv1Res[19], vecDeriv1Res[20], vecDeriv1Res[21],
     vecDeriv1Res[22], vecDeriv1Res[23], vecDeriv2Res[0], vecDeriv2Res[1], vecDeriv2Res[2], vecDeriv2Res[3], vecDeriv2Res[4], vecDeriv2Res[5], vecDeriv2Res[6], vecDeriv2Res[7], vecDeriv2Res[8], vecDeriv2Res[9],
     vecDeriv2Res[10], vecDeriv2Res[11], vecDeriv2Res[12], vecDeriv2Res[13], vecDeriv2Res[14], vecDeriv2Res[15], vecDeriv2Res[16], vecDeriv2Res[17], vecDeriv2Res[18], vecDeriv2Res[19], vecDeriv2Res[20],
     vecDeriv2Res[21], vecDeriv2Res[22], vecDeriv2Res[23], vecDeriv3Res[0], vecDeriv3Res[1], vecDeriv3Res[2], vecDeriv3Res[3], vecDeriv3Res[4], vecDeriv3Res[5], vecDeriv3Res[6], vecDeriv3Res[7]}};

  // calc coeffiecients
  mCoefficients[sThreadnum] = sMatrixA * matrixPar;
}

template <typename DataT, size_t Nx, size_t Ny, size_t Nz>
DataT TriCubicInterpolator<DataT, Nx, Ny, Nz>::interpolate(const Vector<DataT, 3>& pos) const
{
  // the formula for evaluating the interpolation is as follows:
  // f(x,y,z) = \sum_{i,j,k=0}^3 a_{ijk} * x^{i} * y^{j} * z^{k}
  // a_{ijk} is stored in  mCoefficients[]
  //

  const Vector<DataT, FDim> vals0{{1, 1, 1}};   // x^0, y^0, z^0
  const Vector<DataT, FDim> vals2{pos * pos};   // x^2, y^2, z^2
  const Vector<DataT, FDim> vals3{vals2 * pos}; // x^3, y^3, z^3

  const DataT valX[4]{vals0[FX], pos[FX], vals2[FX], vals3[FX]};
  const Vector<DataT, 64> vecValX{
    {valX[0], valX[1], valX[2], valX[3], valX[0], valX[1], valX[2], valX[3], valX[0], valX[1], valX[2], valX[3], valX[0], valX[1], valX[2], valX[3],
     valX[0], valX[1], valX[2], valX[3], valX[0], valX[1], valX[2], valX[3], valX[0], valX[1], valX[2], valX[3], valX[0], valX[1], valX[2], valX[3],
     valX[0], valX[1], valX[2], valX[3], valX[0], valX[1], valX[2], valX[3], valX[0], valX[1], valX[2], valX[3], valX[0], valX[1], valX[2], valX[3],
     valX[0], valX[1], valX[2], valX[3], valX[0], valX[1], valX[2], valX[3], valX[0], valX[1], valX[2], valX[3], valX[0], valX[1], valX[2], valX[3]}};

  const DataT valY[4]{vals0[FY], pos[FY], vals2[FY], vals3[FY]};
  const Vector<DataT, 64> vecValY{
    {valY[0], valY[0], valY[0], valY[0], valY[1], valY[1], valY[1], valY[1], valY[2], valY[2], valY[2], valY[2], valY[3], valY[3], valY[3], valY[3],
     valY[0], valY[0], valY[0], valY[0], valY[1], valY[1], valY[1], valY[1], valY[2], valY[2], valY[2], valY[2], valY[3], valY[3], valY[3], valY[3],
     valY[0], valY[0], valY[0], valY[0], valY[1], valY[1], valY[1], valY[1], valY[2], valY[2], valY[2], valY[2], valY[3], valY[3], valY[3], valY[3],
     valY[0], valY[0], valY[0], valY[0], valY[1], valY[1], valY[1], valY[1], valY[2], valY[2], valY[2], valY[2], valY[3], valY[3], valY[3], valY[3]}};

  const DataT valZ[4]{vals0[FZ], pos[FZ], vals2[FZ], vals3[FZ]};
  const Vector<DataT, 64> vecValZ{
    {valZ[0], valZ[0], valZ[0], valZ[0], valZ[0], valZ[0], valZ[0], valZ[0], valZ[0], valZ[0], valZ[0], valZ[0], valZ[0], valZ[0], valZ[0], valZ[0],
     valZ[1], valZ[1], valZ[1], valZ[1], valZ[1], valZ[1], valZ[1], valZ[1], valZ[1], valZ[1], valZ[1], valZ[1], valZ[1], valZ[1], valZ[1], valZ[1],
     valZ[2], valZ[2], valZ[2], valZ[2], valZ[2], valZ[2], valZ[2], valZ[2], valZ[2], valZ[2], valZ[2], valZ[2], valZ[2], valZ[2], valZ[2], valZ[2],
     valZ[3], valZ[3], valZ[3], valZ[3], valZ[3], valZ[3], valZ[3], valZ[3], valZ[3], valZ[3], valZ[3], valZ[3], valZ[3], valZ[3], valZ[3], valZ[3]}};

  // result = f(x,y,z) = \sum_{i,j,k=0}^3 a_{ijk}    * x^{i}   * y^{j}   * z^{k}
  const DataT result = sum(mCoefficients[sThreadnum] * vecValX * vecValY * vecValZ);
  return result;
}

template <typename DataT, size_t Nx, size_t Ny, size_t Nz>
const Vector<DataT, 3> TriCubicInterpolator<DataT, Nx, Ny, Nz>::processInp(const Vector<DataT, 3>& coordinates, const bool safe) const
{
  const Vector<DataT, FDim> epsilon{{1e-5, 1e-5, 1e-5}};                                                                // add epsilon due to rounding errors
  Vector<DataT, FDim> posRel{(coordinates - mGridProperties.getGridMin()) * mGridProperties.getInvSpacing() + epsilon}; // needed for the grid index

  if (safe) {
    mGridProperties.clampToGridCircularRel(posRel, mCircular);
    mGridProperties.clampToGridRel(posRel, mCircular);
  } else {
    mGridProperties.checkStability(posRel, mCircular);
  }

  // set last row to index last row -1 otherwise two points are missing on the right side of the grid
  const unsigned int ixTmp = static_cast<unsigned int>(posRel[FX]);
  const unsigned int iyTmp = static_cast<unsigned int>(posRel[FY]);
  const unsigned int izTmp = static_cast<unsigned int>(posRel[FZ]);
  const unsigned int ix = (!mCircular[FX] && ixTmp == Nx - 1) ? Nx - 2 : ixTmp;
  const unsigned int iy = (!mCircular[FY] && iyTmp == Ny - 1) ? Ny - 2 : iyTmp;
  const unsigned int iz = (!mCircular[FZ] && izTmp == Nz - 1) ? Nz - 2 : izTmp;

  // std::setprecision(17);
  // std::cout.precision(17);
  // std::cout<<std::endl;
  // std::cout<<"posRel[FZ]: "<<posRel[FZ]<<std::endl;
  // std::cout<<"ix: "<<ix<<std::endl;
  // std::cout<<"izTmp: "<<izTmp<<std::endl;
  // std::cout<<"iz: "<<iz<<std::endl;

  const Vector<unsigned int, FDim> index{{ix, iy, iz}};
  if (!mInitialized[sThreadnum] || !(mLastInd[sThreadnum] == index)) {
    initInterpolator(index[FX], index[FY], index[FZ]);
  }

  const Vector<DataT, FDim> indexTmp{{static_cast<DataT>(ix), static_cast<DataT>(iy), static_cast<DataT>(iz)}};
  const Vector<DataT, FDim> relPos{posRel - indexTmp};
  return relPos;
}

// for perdiodic boundary condition
template <typename DataT, size_t Nx, size_t Ny, size_t Nz>
void TriCubicInterpolator<DataT, Nx, Ny, Nz>::getDataIndexCircularArray(const int index0, const int dim, int arr[]) const
{
  const int delta_min1 = getRegulatedDelta(index0, -1, dim, mGridProperties.getN(dim) - 1);
  const int delta_plus1 = getRegulatedDelta(index0, +1, dim, 1 - mGridProperties.getN(dim));
  const int delta_plus2 = getRegulatedDelta(index0, +2, dim, 2 - mGridProperties.getN(dim));

  arr[0] = mGridProperties.getDeltaDataIndex(delta_min1, dim);
  arr[1] = mGridProperties.getDeltaDataIndex(delta_plus1, dim);
  arr[2] = mGridProperties.getDeltaDataIndex(delta_plus2, dim);
}

// for non perdiodic boundary condition
template <typename DataT, size_t Nx, size_t Ny, size_t Nz>
void TriCubicInterpolator<DataT, Nx, Ny, Nz>::getDataIndexNonCircularArray(const int index0, const int dim, int arr[]) const
{
  const int delta_min1 = getRegulatedDelta(index0, -1, dim, 0);
  const int delta_plus1 = getRegulatedDelta(index0, +1, dim, 0);
  const int delta_plus2 = getRegulatedDelta(index0, +2, dim, delta_plus1);

  arr[0] = mGridProperties.getDeltaDataIndex(delta_min1, dim);
  arr[1] = mGridProperties.getDeltaDataIndex(delta_plus1, dim);
  arr[2] = mGridProperties.getDeltaDataIndex(delta_plus2, dim);
}

template <typename DataT, size_t Nx, size_t Ny, size_t Nz>
DataT TriCubicInterpolator<DataT, Nx, Ny, Nz>::uiPow(DataT base, unsigned int exponent) const
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

template <typename DataT, size_t Nx, size_t Ny, size_t Nz>
int TriCubicInterpolator<DataT, Nx, Ny, Nz>::findPos(const int ix, const int iy, const int iz) const
{
  int pos = -1;
  if (isInInnerVolume(ix, iy, iz, pos)) {
    return pos;
  }

  if (findEdge(ix, iy, iz, pos)) {
    return pos;
  }

  if (findLine(ix, iy, iz, pos)) {
    return pos;
  }

  if (findSide(ix, iy, iz, pos)) {
    return pos;
  }
  return -1;
}

template <typename DataT, size_t Nx, size_t Ny, size_t Nz>
bool TriCubicInterpolator<DataT, Nx, Ny, Nz>::findEdge(const int ix, const int iy, const int iz, int& posType) const
{
  const int iR = 2;
  if (ix == 0 && iy == 0) {
    if (iz == 0) {
      // edge 0
      posType = GridPos::Edge0;
      return true;
    } else if (iz == Nz - iR) {
      //edge 4
      posType = GridPos::Edge4;
      return true;
    }
  } else if (ix == Nx - iR && iy == 0) {
    if (iz == 0) {
      // edge 1
      posType = GridPos::Edge1;
      return true;
    } else if (iz == Nz - iR) {
      //edge 5
      posType = GridPos::Edge5;
      return true;
    }
  } else if (ix == 0 && iy == Ny - iR) {
    if (iz == 0) {
      // edge 2
      posType = GridPos::Edge2;
      return true;
    } else if (iz == Nz - iR) {
      //edge 6
      posType = GridPos::Edge6;
      return true;
    }
  } else if (ix == Nx - iR && iy == Ny - iR) {
    if (iz == 0) {
      // edge 3
      posType = GridPos::Edge3;
      return true;
    } else if (iz == Nz - iR) {
      //edge 7
      posType = GridPos::Edge7;
      return true;
    }
  }
  return false;
}

template <typename DataT, size_t Nx, size_t Ny, size_t Nz>
bool TriCubicInterpolator<DataT, Nx, Ny, Nz>::findLine(const int ix, const int iy, const int iz, int& posType) const
{
  const int iR = 2;
  //check line
  if (iy == 0) {
    if (iz == 0) {
      // line a
      posType = GridPos::LineA;
      return true;
    } else if (iz == Nz - iR) {
      // line e
      posType = GridPos::LineE;
      return true;
    }
    if (ix == 0) {
      posType = GridPos::LineI;
      return true;
    } else if (ix == Nx - iR) {
      posType = GridPos::LineJ;
      return true;
    }
  } else if (iy == Ny - iR) {
    if (iz == 0) {
      // line b
      posType = GridPos::LineB;
      return true;
    } else if (iz == Nz - iR) {
      // line f
      posType = GridPos::LineF;
      return true;
    }
    if (ix == 0) {
      posType = GridPos::LineK;
      return true;
    } else if (ix == Nx - iR) {
      posType = GridPos::LineL;
      return true;
    }
  } else if (ix == 0) {
    if (iz == 0) {
      //line c
      posType = GridPos::LineC;
      return true;
    } else if (iz == Nz - iR) {
      // line g
      posType = GridPos::LineG;
      return true;
    }
  } else if (ix == Nx - iR) {
    if (iz == 0) {
      //line d
      posType = GridPos::LineD;
      return true;
    } else if (iz == Nz - iR) {
      // line h
      posType = GridPos::LineH;
      return true;
    }
  }
  return false;
}

template <typename DataT, size_t Nx, size_t Ny, size_t Nz>
bool TriCubicInterpolator<DataT, Nx, Ny, Nz>::findSide(const int ix, const int iy, const int iz, int& posType) const
{
  if (isSideRight(ix, FX)) {
    posType = GridPos::SideXRight;
    return true;
  } else if (isSideLeft(ix)) {
    posType = GridPos::SideXLeft;
    return true;
  }
  if (isSideRight(iy, FY)) {
    posType = GridPos::SideYRight;
    return true;
  } else if (isSideLeft(iy)) {
    posType = GridPos::SideYLeft;
    return true;
  }
  if (isSideRight(iz, FZ)) {
    posType = GridPos::SideZRight;
    return true;
  } else if (isSideLeft(iz)) {
    posType = GridPos::SideZLeft;
    return true;
  }
  return false;
}

} // namespace tpc
} // namespace o2

#endif
