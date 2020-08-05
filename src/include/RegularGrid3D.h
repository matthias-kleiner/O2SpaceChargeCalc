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

template <typename DataT = float, unsigned int Nx = 4, unsigned int Ny = 4, unsigned int Nz = 4>
struct RegularGrid3D {

 public:
  RegularGrid3D(const DataT* gridData, const DataT xmin, const DataT ymin, const DataT zmin, const DataT spacingX, const DataT spacingY, const DataT spacingZ) : mMin{{xmin, ymin, zmin}}, mSpacing{{spacingX, spacingY, spacingZ}}, mInvSpacing{{static_cast<DataT>(1 / spacingX), static_cast<DataT>(1 / spacingY), static_cast<DataT>(1 / spacingZ)}}
  {
    initArray(gridData);
    initLists();
  }

  RegularGrid3D(const DataT xmin, const DataT ymin, const DataT zmin, const DataT spacingX, const DataT spacingY, const DataT spacingZ) : mMin{{xmin, ymin, zmin}}, mSpacing{{spacingX, spacingY, spacingZ}}, mInvSpacing{{static_cast<DataT>(1 / spacingX), static_cast<DataT>(1 / spacingY), static_cast<DataT>(1 / spacingZ)}}
  {
    initLists();
  }

  RegularGrid3D(const DataT xmin, const DataT ymin, const DataT zmin, const DataT spacingX, const DataT spacingY, const DataT spacingZ, TFile& inpf, const char* name = "data") : mMin{{xmin, ymin, zmin}}, mSpacing{{spacingX, spacingY, spacingZ}}, mInvSpacing{{static_cast<DataT>(1 / spacingX), static_cast<DataT>(1 / spacingY), static_cast<DataT>(1 / spacingZ)}}
  {
    initLists();
    initFromFile(inpf, name);
  }

  void initArray(const DataT* gridData)
  {
    memcpy(mGridData.mData.get(), gridData, Nx * Ny * Nz * sizeof(DataT));
  }

  /// \param ix x vertex
  /// \param ix y vertex
  /// \param ix z vertex
  /// \return returns the index to vertex
  size_t getDataIndex(const size_t ix, const size_t iy, const size_t iz) const
  {
    return mGridData.getDataIndex(ix, iy, iz);
  }

  /// \param deltaX delta x index
  /// \return returns the delta index (where the data is stored) for given deltaX
  int getDeltaXDataIndex(const int deltaX) const
  {
    return deltaX;
  }

  /// \param deltaY delta y index
  /// \return returns the delta index (where the data is stored) for given deltaY
  int getDeltaYDataIndex(const int deltaY) const
  {
    return Nx * deltaY;
  }

  /// \param deltaZ delta z index
  /// \return returns the delta index (where the data is stored) for given deltaZ
  int getDeltaZDataIndex(const int deltaZ) const
  {
    return deltaZ * Ny * Nx;
  }

  // same as above
  /// \param delta delta index
  /// \param dim dimension of interest
  /// \return returns the delta index (where the data is stored) for given delta and dim
  int getDeltaDataIndex(const int delta, const int dim) const
  {
    const unsigned int offset[FDim]{1, Nx, Ny * Nx};
    const int deltaIndex = delta * offset[dim];
    return deltaIndex;
  }

  // check if the specified index for given dimension lies in the grid
  /// \param index query index
  /// \return returns if the index lies in the grid
  bool isIndexInGrid(const int index, const unsigned int dim) const
  {
    return index < 0 ? false : (index > (FNdim[dim] - 1) ? false : true);
  }

  // check if the specified position for given dimension lies in the grid
  /// \param pos query position
  /// \return returns if the position lies in the grid
  bool isInGrid(const DataT pos, const unsigned int dim) const
  {
    if (pos < mMin[dim]) {
      return false;
    } else if (pos > mMin[dim] + getN(dim) * mSpacing[dim]) {
      return false;
    }
    return true;
  }

  const DataT& operator[](size_t i) const { return mGridData[i]; }
  DataT& operator[](size_t i) { return mGridData[i]; }

  const DataT& operator()(size_t ix, size_t iy, size_t iz) const
  {
    return mGridData(ix, iy, iz);
  }

  DataT& operator()(size_t ix, size_t iy, size_t iz)
  {
    return mGridData(ix, iy, iz);
  }

  /// \param dim dimension of interest
  /// \return returns the number of vertices for given dimension for the grid
  static constexpr size_t getN(unsigned int dim) { return FNdim[dim]; }
  static constexpr size_t getNX() { return FNdim[FX]; }
  static constexpr size_t getNY() { return FNdim[FY]; }
  static constexpr size_t getNZ() { return FNdim[FZ]; }

  static constexpr unsigned int getDim() { return FDim; } /// \return returns number of dimensions of the grid (3)
  static constexpr unsigned int getFX() { return FX; }    /// \return returns the index for dimension x (0)
  static constexpr unsigned int getFY() { return FY; }    /// \return returns the index for dimension y (1)
  static constexpr unsigned int getFZ() { return FZ; }    /// \return returns the index for dimension z (2)

  const Vector<DataT, 3>& getGridMin() const { return mMin; } /// \return returns the minimum coordinates of the grid in all dimensions
  DataT getGridMinX() const { return mMin[FX]; }              /// \return returns the minimum coordinate of the grid in x dimension
  DataT getGridMinY() const { return mMin[FY]; }              /// \return returns the minimum coordinate of the grid in y dimension
  DataT getGridMinZ() const { return mMin[FZ]; }              /// \return returns the minimum coordinate of the grid in z dimension

  DataT getGridMaxX() const { return mMin[FX] + getNX() / getInvSpacingX(); }
  DataT getGridMaxY() const { return mMin[FY] + getNY() / getInvSpacingY(); }
  DataT getGridMaxZ() const { return mMin[FZ] + getNZ() / getInvSpacingZ(); }

  /// \return returns the inversed spacing of the grid for all dimensions
  const Vector<DataT, 3>& getInvSpacing() const { return mInvSpacing; }
  DataT getInvSpacingX() const { return mInvSpacing[FX]; }
  DataT getInvSpacingY() const { return mInvSpacing[FY]; }
  DataT getInvSpacingZ() const { return mInvSpacing[FZ]; }

  // clamp coordinates to the grid (not circular)
  /// \param pos query position which will be clamped
  /// \return returns clamped coordinate coordinate
  DataT clampToGrid(const DataT pos, const unsigned int dim) const
  {
    if (pos < mMin[dim]) {
      return mMin[dim];
    } else if (pos > mMin[dim] + getN(dim) / mInvSpacing[dim]) {
      return mMin[dim] + getN(dim) / mInvSpacing[dim];
    }
    return pos;
  }

  // clamp coordinates to the grid (not circular) vectorized
  /// \param relPos query position which will be clamped in relative grid cooridnates
  /// \param circular Vector containing if a dimension is circular. e.g x=not circular, y=not circular, z=circular -> circular(0,0,1)
  template <size_t FDim>
  void clampToGrid(Vector<DataT, FDim>& relPos, const Vector<int, FDim>& circular) const
  {
    for (size_t j = 0; j < relPos.getvectorsCount(); ++j) {
      auto vecTmp = relPos.getVector(j);
      const auto maskLeft = vecTmp < 0;
      const auto maskRight = vecTmp > (mMaxIndex.getVector(j) + simd_cast<Vc::Vector<DataT>>(circular.getVector(j)));
      vecTmp(maskLeft) = Vc::Vector<DataT>::Zero();
      vecTmp(maskRight) = mMaxIndex.getVector(j);
      relPos.setVector(j, vecTmp);
    }
  }

  // clamp coordinate to the grid (circular)
  /// \param relPos query position which will be clamped in relative grid cooridnates
  /// \param dim dimension of relative position
  /// \return returns the circular coordinate
  DataT clampToGridCircular(DataT relPos, const unsigned int dim) const
  {
    while (relPos < 0) {
      relPos += FNdim[dim];
    }
    while (relPos > FNdim[dim]) {
      relPos -= FNdim[dim];
    }
    return relPos;
  }

  // clamp coordinates to the grid (circular)
  /// \param relPos query position which will be clamped in relative grid cooridnates
  /// \param circular Vector containing if a dimension is circular. e.g x=not circular, y=not circular, z=circular -> circular(0,0,1)
  template <size_t FDim>
  void clampToGridCircular(Vector<DataT, FDim>& relPos, const Vector<int, FDim>& circular) const
  {
    for (size_t j = 0; j < relPos.getvectorsCount(); ++j) {
      auto vecTmp = (relPos.getVector(j));
      auto vecIsCircular = simd_cast<Vc::Vector<DataT>>(circular.getVector(j));
      auto vecTmp3 = simd_cast<Vc::Vector<DataT>>(FNdim.getVector(j));
      while ((vecTmp * vecIsCircular < 0).count() > 0) {
        vecTmp = vecTmp + vecTmp3 * vecIsCircular;
      }
      while ((vecTmp * vecIsCircular > vecTmp3).count() > 0) {
        vecTmp = vecTmp - vecTmp3 * vecIsCircular;
      }
      // numerical stability check
      vecTmp(vecTmp * vecIsCircular == vecTmp3) = Vc::Vector<DataT>::Zero();

      relPos.setVector(j, vecTmp);
    }
  }

  template <size_t FDim>
  void checkStability(Vector<DataT, FDim>& relPos, const Vector<int, FDim>& circular) const
  {
    for (size_t j = 0; j < relPos.getvectorsCount(); ++j) {
      auto vecTmp = (relPos.getVector(j));
      auto vecIsCircular = simd_cast<Vc::Vector<DataT>>(circular.getVector(j));
      auto vecTmp3 = simd_cast<Vc::Vector<DataT>>(FNdim.getVector(j));
      vecTmp(vecTmp * vecIsCircular == vecTmp3) = Vc::Vector<DataT>::Zero();
      relPos.setVector(j, vecTmp);
    }
  }

  /// \param vertex in x dimension
  /// \return returns the x positon for given vertex
  DataT getXVertex(const size_t vertex) const
  {
    return mXVertices[vertex];
  }

  /// \param vertex in y dimension
  /// \return returns the y positon for given vertex
  DataT getYVertex(const size_t index) const
  {
    return mYVertices[index];
  }

  /// \param vertex in z dimension
  /// \return returns the z positon for given vertex
  DataT getZVertex(const size_t index) const
  {
    return mZVertices[index];
  }

  // set the values of the grid from root file
  void initFromFile(TFile& inpf, const char* name = "data")
  {
    const auto tmpContainer = DataContainer3D<DataT, Nx, Ny, Nz>::readFromFile(inpf, name);
    if (!tmpContainer) {
      std::cout << "Failed to load " << name << " from " << inpf.GetName() << std::endl;
      return;
    }
    memcpy(mGridData.mData.get(), tmpContainer->mData.get(), tmpContainer->getNDataPoints() * sizeof(DataT));
  }

  // set the values of the grid from root file
  int storeValuesToFile(TFile& outf, const char* name = "data") const
  {
    return mGridData.writeToFile(outf, name);
  }

 private:
  static constexpr unsigned int FDim = 3; // dimensions of the grid (only 3 supported)
  static constexpr unsigned int FX = 0;   // index for x coordinate
  static constexpr unsigned int FY = 1;   // index for y coordinate
  static constexpr unsigned int FZ = 2;   // index for z coordinate

  DataContainer3D<DataT, Nx, Ny, Nz> mGridData;

  const Vector<DataT, FDim> mMin{};                              // min positions of grid
  const Vector<DataT, FDim> mSpacing{};                          //  spacing of the grid
  const Vector<DataT, FDim> mInvSpacing{};                       // inverse spacing of grid
  const inline static Vector<DataT, FDim> mMaxIndex{{Nx - 1, Ny - 1, Nz - 1}}; // max index which is on the grid in all dimensions
  inline static Vector<int, FDim> FNdim{{Nx, Ny, Nz}};

  DataT mXVertices[Nx]{};
  DataT mYVertices[Ny]{};
  DataT mZVertices[Nz]{};

  void initLists()
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

  ClassDefNV(RegularGrid3D, 1)
};

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

#endif
