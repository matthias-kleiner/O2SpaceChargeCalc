#ifndef RegularGrid3D_H
#define RegularGrid3D_H

#include <iostream>
#include <iomanip>
#include "matrix.h"

template <typename DataT = float, unsigned int Nx = 4, unsigned int Ny = 4, unsigned int Nz = 4>
struct RegularGrid3D {

 public:
  RegularGrid3D(const DataT* gridData, const DataT xmin, const DataT ymin, const DataT zmin, const DataT spacingX, const DataT spacingY, const DataT spacingZ) : mMin{{xmin, ymin, zmin}}, mInvSpacing{{static_cast<DataT>(1 / spacingX), static_cast<DataT>(1 / spacingY), static_cast<DataT>(1 / spacingZ)}}
  {
    initArray(gridData);
  }

  RegularGrid3D(const DataT xmin=0, const DataT ymin=0, const DataT zmin=0, const DataT spacingX=1, const DataT spacingY=1, const DataT spacingZ=1) : mMin{{xmin, ymin, zmin}}, mInvSpacing{{static_cast<DataT>(1 / spacingX), static_cast<DataT>(1 / spacingY), static_cast<DataT>(1 / spacingZ)}} {}

  void initArray(const DataT* gridData)
  {
    memcpy(mGridData.get(), gridData, Nx * Ny * Nz * sizeof(DataT));
  }

  /// \param ix x vertex
  /// \param ix y vertex
  /// \param ix z vertex
  /// \return returns the index to vertex
  size_t getDataIndex(const size_t ix, const size_t iy, const size_t iz) const
  {
    const size_t indX = ix;
    const size_t indY = iy;
    const size_t indZ = iz;
    const size_t index = indX + Nx * (indY + indZ * Ny);
    return index;
  }

  int getDeltaXDataIndex(const int deltaX) const
  {
    const int deltaXindex = deltaX;
    return deltaXindex;
  }

  int getDeltaYDataIndex(const int deltaY) const
  {
    const int deltaYindex = Nx * deltaY;
    return deltaYindex;
  }

  int getDeltaZDataIndex(const int deltaZ) const
  {
    const int deltaZindex = deltaZ * Ny * Nx;
    return deltaZindex;
  }

  int getDeltaDataIndex(const int delta, const int dim) const
  {
    const unsigned int offset[FDim]{1, Nx, Ny * Nx};
    const int deltaIndex = delta * offset[dim];
    return deltaIndex;
  }

  // check if the specified index for given dimension lies in the grid
  bool isIndexInGrid(const int indexDim, const unsigned int dim) const
  {
    return indexDim < 0 ? false : (indexDim > static_cast<int>(FNdim[dim] - 1) ? false : true);
  }

  // check if the specified position for given dimension lies in the grid
  bool isInGrid(const DataT pos, const unsigned int dim) const
  {
    if (pos < mMin[dim]) {
      return false;
    } else if (pos > mMin[dim] + getN(dim) / mInvSpacing[dim]) {
      return false;
    }
    return true;
  }

  const DataT& operator[](size_t i) const { return mGridData.get()[i]; }
  DataT& operator[](size_t i) { return mGridData.get()[i]; }

  const DataT& operator()(size_t ix, size_t iy, size_t iz) const
  {
    const size_t ind = getDataIndex(ix, iy, iz);
    return mGridData.get()[ind];
  }

  DataT& operator()(size_t ix, size_t iy, size_t iz)
  {
    const size_t ind = getDataIndex(ix, iy, iz);
    return mGridData.get()[ind];
  }

  /// \param dim dimension of interest
  /// \return returns the number of vertices for given dimension for the new grid
  static constexpr size_t getN(unsigned int dim) { return FNdim[dim]; }
  static constexpr size_t getNX() { return FNdim[FX]; }
  static constexpr size_t getNY() { return FNdim[FY]; }
  static constexpr size_t getNZ() { return FNdim[FZ]; }

  static constexpr unsigned int getDim() { return FDim; }
  static constexpr unsigned int getFX() { return FX; }
  static constexpr unsigned int getFY() { return FY; }
  static constexpr unsigned int getFZ() { return FZ; }

  /// \return returns the minimum coordinates of the grid in all dimensions
  Vector<DataT, 3> getGridMin() const { return mMin; }
  DataT getGridMinX() const { return mMin[FX]; }
  DataT getGridMinY() const { return mMin[FY]; }
  DataT getGridMinZ() const { return mMin[FZ]; }

  DataT getGridMaxX() const { return mMin[FX] + getNX() / getInvSpacingX(); }
  DataT getGridMaxY() const { return mMin[FY] + getNY() / getInvSpacingY(); }
  DataT getGridMaxZ() const { return mMin[FZ] + getNZ() / getInvSpacingZ(); }

  /// \return returns the inversed spacing of the grid for all dimensions
  Vector<DataT, 3> getInvSpacing() const { return mInvSpacing; }
  DataT getInvSpacingX() const { return mInvSpacing[FX]; }
  DataT getInvSpacingY() const { return mInvSpacing[FY]; }
  DataT getInvSpacingZ() const { return mInvSpacing[FZ]; }

  // clamp coordinates to the grid
  /// \param pos coordinates in x,y,z dimensions
  template <size_t FDim>
  void clampToGrid(Vector<DataT, FDim>& pos) const
  {
    for (size_t j = 0; j < pos.GetvectorsCount(); ++j) {
      auto vecTmp = pos.GetVector(j);
      const auto maskLeft = vecTmp < 0;
      const auto maskRight = vecTmp > mMaxIndex.GetVector(j);
      vecTmp(maskLeft) = Vc::Vector<DataT>::Zero();
      vecTmp(maskRight) = mMaxIndex.GetVector(j);
      pos.setVector(j, vecTmp);
    }
  }

 private:
  static constexpr unsigned int FDim = 3; // dimensions of the grid (only 3 supported)
  static constexpr unsigned int FX = 0;   // index for x coordinate
  static constexpr unsigned int FY = 1;   // index for y coordinate
  static constexpr unsigned int FZ = 2;   // index for z coordinate

  const Vector<DataT, FDim> mMin{};                              // min positions of grid
  const Vector<DataT, FDim> mInvSpacing{};                       // inverse spacing of grid
  const Vector<DataT, FDim> mMaxIndex{{Nx - 1, Ny - 1, Nz - 1}}; // max index which is on the grid in all dimensions

  static constexpr size_t FNdim[FDim]{Nx, Ny, Nz};
  static constexpr size_t FNdataPoints{Nx * Ny * Nz};

  std::unique_ptr<DataT[]> mGridData = std::make_unique<DataT[]>(FNdataPoints);
};

template <typename DataT, unsigned int Nx, unsigned int Ny, unsigned int Nz>
std::ostream& operator<<(std::ostream& out, const RegularGrid3D<DataT, Nx, Ny, Nz>& cont)
{
  out.precision(3);
  auto&& w = std::setw(6);

  for (unsigned int iz = 0; iz < cont.getNZ(); ++iz) {
    // print top x row
    out << "⎡" << w << cont(0, 0, iz);
    for (unsigned int ix = 1; ix < cont.getNX(); ++ix) {
      out << ", " << w << cont(ix, 0, iz);
    }
    out << " ⎤" << std::endl;

    for (unsigned int iy = 1; iy < cont.getNY() - 1; ++iy) {
      out << "⎢" << w << cont(0, iy, iz);
      for (unsigned int ix = 1; ix < cont.getNX(); ++ix) {
        out << ", " << w << cont(ix, iy, iz);
      }
      out << " ⎥" << std::endl;
    }

    out << "⎣" << w << cont(0, cont.getNY() - 1, iz);
    for (unsigned int ix = 1; ix < cont.getNX(); ++ix) {
      out << ", " << w << cont(ix, cont.getNY() - 1, iz);
    }
    out << " ⎦" << std::endl;

    out << std::endl;
    out << std::endl;
  }
  return out;
}

#endif
