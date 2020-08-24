// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  RegularGrid3D.h
/// \brief Definition of RegularGrid3D class
///
/// \author  Matthias Kleiner <matthias.kleiner@cern.ch>

#ifndef RegularGrid3D_H
#define RegularGrid3D_H

#include <iostream>
#include <iomanip>
#include "Vector.h"
#include "Rtypes.h" // for ClassDefNV
#include "DataContainer3D.h"

namespace o2
{
namespace tpc
{

/// \tparam DataT the type of data which is used during the calculations
/// \tparam Nx number of vertices in x direction
/// \tparam Ny number of vertices in y direction
/// \tparam Nz number of vertices in z direction
template <typename DataT = float, unsigned int Nx = 4, unsigned int Ny = 4, unsigned int Nz = 4>
struct RegularGrid3D {

 public:
  RegularGrid3D(const DataT xmin, const DataT ymin, const DataT zmin, const DataT spacingX, const DataT spacingY, const DataT spacingZ) : mMin{{xmin, ymin, zmin}}, mMax{{xmin + (Nx - 1)  * spacingX, ymin + (Ny - 1) * spacingY, zmin + (Nz - 1)  * spacingZ }}, mSpacing{{spacingX, spacingY, spacingZ}}, mInvSpacing{{static_cast<DataT>(1 / spacingX), static_cast<DataT>(1 / spacingY), static_cast<DataT>(1 / spacingZ)}}
  {
    initLists();
  }

  /// \param deltaX delta x index
  /// \return returns the delta index (where the data is stored) for given deltaX
  int getDeltaXDataIndex(const int deltaX) const { return deltaX; }

  /// \param deltaY delta y index
  /// \return returns the delta index (where the data is stored) for given deltaY
  int getDeltaYDataIndex(const int deltaY) const { return Nx * deltaY; }

  /// \param deltaZ delta z index
  /// \return returns the delta index (where the data is stored) for given deltaZ
  int getDeltaZDataIndex(const int deltaZ) const { return deltaZ * Ny * Nx; }

  // same as above
  /// \param delta delta index
  /// \param dim dimension of interest
  /// \return returns the delta index (where the data is stored) for given delta and dim
  int getDeltaDataIndex(const int delta, const int dim) const;

  // check if the specified index for given dimension lies in the grid
  /// \param index query index
  /// \return returns if the index lies in the grid
  bool isIndexInGrid(const int index, const unsigned int dim) const { return index < 0 ? false : (index > (sNdim[dim] - 1) ? false : true); }

  /// \param dim dimension of interest
  /// \return returns the number of vertices for given dimension for the grid
  static constexpr size_t getN(unsigned int dim) { return sNdim[dim]; }
  static constexpr size_t getNX() { return sNdim[FX]; }
  static constexpr size_t getNY() { return sNdim[FY]; }
  static constexpr size_t getNZ() { return sNdim[FZ]; }

  static constexpr unsigned int getDim() { return FDIM; } /// \return returns number of dimensions of the grid (3)
  static constexpr unsigned int getFX() { return FX; }    /// \return returns the index for dimension x (0)
  static constexpr unsigned int getFY() { return FY; }    /// \return returns the index for dimension y (1)
  static constexpr unsigned int getFZ() { return FZ; }    /// \return returns the index for dimension z (2)

  const Vector<DataT, 3>& getGridMin() const { return mMin; } /// \return returns the minimum coordinates of the grid in all dimensions
  DataT getGridMinX() const { return mMin[FX]; }              /// \return returns the minimum coordinate of the grid in x dimension
  DataT getGridMinY() const { return mMin[FY]; }              /// \return returns the minimum coordinate of the grid in y dimension
  DataT getGridMinZ() const { return mMin[FZ]; }              /// \return returns the minimum coordinate of the grid in z dimension

  DataT getGridMaxX() const { return mMax[FX]; }
  DataT getGridMaxY() const { return mMax[FY]; }
  DataT getGridMaxZ() const { return mMax[FZ]; }

  /// \return returns the inversed spacing of the grid for all dimensions
  const Vector<DataT, 3>& getInvSpacing() const { return mInvSpacing; }
  DataT getInvSpacingX() const { return mInvSpacing[FX]; }
  DataT getInvSpacingY() const { return mInvSpacing[FY]; }
  DataT getInvSpacingZ() const { return mInvSpacing[FZ]; }

  DataT getSpacingX() const { return mSpacing[FX]; }
  DataT getSpacingY() const { return mSpacing[FY]; }
  DataT getSpacingZ() const { return mSpacing[FZ]; }

  // clamp coordinates to the grid (not circular)
  /// \param pos query position which will be clamped
  /// \return returns clamped coordinate coordinate
  DataT clampToGrid(const DataT pos, const unsigned int dim) const;

  // clamp coordinates to the grid (not circular)
  /// \param pos query position which will be clamped
  /// \return returns clamped coordinate coordinate
  DataT clampToGridCircular(DataT pos, const unsigned int dim) const;

  // clamp coordinates to the grid (not circular) vectorized
  /// \param relPos query position which will be clamped in relative grid cooridnates
  /// \param circular Vector containing if a dimension is circular. e.g x=not circular, y=not circular, z=circular -> circular(0,0,1)
  void clampToGridRel(Vector<DataT, 3>& relPos, const Vector<int, 3>& circular) const;

  // clamp coordinates to the grid (circular)
  /// \param relPos query position which will be clamped in relative grid cooridnates
  /// \param circular Vector containing if a dimension is circular. e.g x=not circular, y=not circular, z=circular -> circular(0,0,1)
  void clampToGridCircularRel(Vector<DataT, 3>& relPos, const Vector<int, 3>& circular) const;

  void checkStability(Vector<DataT, 3>& relPos, const Vector<int, 3>& circular) const;

  /// \param vertex in x dimension
  /// \return returns the x positon for given vertex
  DataT getXVertex(const size_t vertex) const { return mXVertices[vertex]; }

  /// \param vertex in y dimension
  /// \return returns the y positon for given vertex
  DataT getYVertex(const size_t index) const { return mYVertices[index]; }

  /// \param vertex in z dimension
  /// \return returns the z positon for given vertex
  DataT getZVertex(const size_t index) const { return mZVertices[index]; }

 private:
  static constexpr unsigned int FDIM = 3; // dimensions of the grid (only 3 supported)
  static constexpr unsigned int FX = 0;   // index for x coordinate
  static constexpr unsigned int FY = 1;   // index for y coordinate
  static constexpr unsigned int FZ = 2;   // index for z coordinate

  const Vector<DataT, FDIM> mMin{};                                            // min vertices positions of the grid
  const Vector<DataT, FDIM> mMax{};                                            // max vertices positions of the grid
  const Vector<DataT, FDIM> mSpacing{};                                        //  spacing of the grid
  const Vector<DataT, FDIM> mInvSpacing{};                                     // inverse spacing of grid
  const inline static Vector<DataT, FDIM> sMaxIndex{{Nx - 1, Ny - 1, Nz - 1}}; // max index which is on the grid in all dimensions
  inline static Vector<int, FDIM> sNdim{{Nx, Ny, Nz}};

  DataT mXVertices[Nx]{};
  DataT mYVertices[Ny]{};
  DataT mZVertices[Nz]{};

  void initLists();

  ClassDefNV(RegularGrid3D, 1)
};

///
/// ========================================================================================================
///       Inline implementations of some methods
/// ========================================================================================================
///

template <typename DataT, unsigned int Nx, unsigned int Ny, unsigned int Nz>
DataT RegularGrid3D<DataT, Nx, Ny, Nz>::clampToGrid(const DataT pos, const unsigned int dim) const
{
  if (pos < mMin[dim]) {
    return mMin[dim];
  } else if (pos > mMax[dim]) {
    return mMax[dim];
  }
  return pos;
}

template <typename DataT, unsigned int Nx, unsigned int Ny, unsigned int Nz>
DataT RegularGrid3D<DataT, Nx, Ny, Nz>::clampToGridCircular(DataT pos, const unsigned int dim) const
{
  while (pos < mMin[dim]) {
    pos += mMax[dim] - mMin[dim] + mSpacing[dim];
  }
  while (pos >= mMax[dim] + mSpacing[dim]) {
    pos -= mMax[dim] + mSpacing[dim] - mMin[dim];
  }
  return pos;
}


template <typename DataT, unsigned int Nx, unsigned int Ny, unsigned int Nz>
void RegularGrid3D<DataT, Nx, Ny, Nz>::clampToGridRel(Vector<DataT, 3>& relPos, const Vector<int, 3>& circular) const
{
  for (size_t j = 0; j < relPos.getvectorsCount(); ++j) {
    auto vecTmp = relPos.getVector(j);
    const auto maskLeft = vecTmp < 0;
    const auto maskRight = vecTmp > (sMaxIndex.getVector(j) + simd_cast<Vc::Vector<DataT>>(circular.getVector(j)));
    vecTmp(maskLeft) = Vc::Vector<DataT>::Zero();
    vecTmp(maskRight) = sMaxIndex.getVector(j);
    relPos.setVector(j, vecTmp);
  }
}

template <typename DataT, unsigned int Nx, unsigned int Ny, unsigned int Nz>
void RegularGrid3D<DataT, Nx, Ny, Nz>::clampToGridCircularRel(Vector<DataT, 3>& relPos, const Vector<int, 3>& circular) const
{
  for (size_t j = 0; j < relPos.getvectorsCount(); ++j) {
    auto vecTmp = (relPos.getVector(j));
    auto vecIsCircular = simd_cast<Vc::Vector<DataT>>(circular.getVector(j));
    auto vecTmp3 = simd_cast<Vc::Vector<DataT>>(sNdim.getVector(j));
    while ((vecTmp * vecIsCircular < 0).count() > 0) {
      vecTmp = vecTmp + vecTmp3 * vecIsCircular;
    }
    vecTmp(vecTmp * vecIsCircular == vecTmp3) = Vc::Vector<DataT>::Zero();
    while ((vecTmp * vecIsCircular > vecTmp3).count() > 0) {
      vecTmp = vecTmp - vecTmp3 * vecIsCircular;
    }
    relPos.setVector(j, vecTmp);
  }
}

template <typename DataT, unsigned int Nx, unsigned int Ny, unsigned int Nz>
void RegularGrid3D<DataT, Nx, Ny, Nz>::checkStability(Vector<DataT, 3>& relPos, const Vector<int, 3>& circular) const
{
  for (size_t j = 0; j < relPos.getvectorsCount(); ++j) {
    auto vecTmp = (relPos.getVector(j));
    auto vecIsCircular = simd_cast<Vc::Vector<DataT>>(circular.getVector(j));
    auto vecTmp3 = simd_cast<Vc::Vector<DataT>>(sNdim.getVector(j));
    vecTmp(vecTmp * vecIsCircular == vecTmp3) = Vc::Vector<DataT>::Zero();
    relPos.setVector(j, vecTmp);
  }
}

template <typename DataT, unsigned int Nx, unsigned int Ny, unsigned int Nz>
void RegularGrid3D<DataT, Nx, Ny, Nz>::initLists()
{
  for (size_t i = 0; i < Nx; ++i) {
    mXVertices[i] = mMin[FX] + i * mSpacing[FX];
  }
  for (size_t i = 0; i < Ny; ++i) {
    mYVertices[i] = mMin[FY] + i * mSpacing[FY];
  }
  for (size_t i = 0; i < Nz; ++i) {
    mZVertices[i] = mMin[FZ] + i * mSpacing[FZ];
  }
}

template <typename DataT, unsigned int Nx, unsigned int Ny, unsigned int Nz>
int RegularGrid3D<DataT, Nx, Ny, Nz>::getDeltaDataIndex(const int delta, const int dim) const
{
  const unsigned int offset[FDIM]{1, Nx, Ny * Nx};
  const int deltaIndex = delta * offset[dim];
  return deltaIndex;
}

template <template <typename, unsigned int, unsigned int, unsigned int> typename MClass, typename DataT, unsigned int Nx, unsigned int Ny, unsigned int Nz>
std::ostream& operator<<(std::ostream& out, const MClass<DataT, Nx, Ny, Nz>& cont)
{
  out.precision(3);
  auto&& w = std::setw(9);
  out << std::endl;

  for (unsigned int iz = 0; iz < Nz; ++iz) {
    out << "z layer: " << iz << std::endl;
    // print top x row
    out << "⎡" << w << cont(0, 0, iz);
    for (unsigned int ix = 1; ix < Nx; ++ix) {
      out << ", " << w << cont(ix, 0, iz);
    }
    out << " ⎤" << std::endl;

    for (unsigned int iy = 1; iy < Ny - 1; ++iy) {
      out << "⎢" << w << cont(0, iy, iz);
      for (unsigned int ix = 1; ix < Nx; ++ix) {
        out << ", " << w << cont(ix, iy, iz);
      }
      out << " ⎥" << std::endl;
    }

    out << "⎣" << w << cont(0, Ny - 1, iz);
    for (unsigned int ix = 1; ix < Nx; ++ix) {
      out << ", " << w << cont(ix, Ny - 1, iz);
    }
    out << " ⎦" << std::endl;
    out << std::endl;
    out << std::endl;
  }
  return out;
}

}
}

#endif
