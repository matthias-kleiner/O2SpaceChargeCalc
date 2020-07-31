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
#include "matrix.h"
#include "Rtypes.h" // for ClassDefNV

// i/o
#include "TFile.h"

template <typename DataT = float, unsigned int Nx = 4, unsigned int Ny = 4, unsigned int Nz = 4>
struct DataContainer3D {

  static constexpr size_t FNdataPoints{Nx * Ny * Nz}; ///< number of values stored in the container
  std::unique_ptr<DataT[]> mData = std::make_unique<DataT[]>(FNdataPoints);

  const DataT& operator[](size_t i) const { return mData[i]; }
  DataT& operator[](size_t i) { return mData[i]; }

  const DataT& operator()(size_t ix, size_t iy, size_t iz) const
  {
    const size_t ind = getDataIndex(ix, iy, iz);

    return mData[ind];
  }

  DataT& operator()(size_t ix, size_t iy, size_t iz)
  {
    const size_t ind = getDataIndex(ix, iy, iz);
    return mData[ind];
  }

  size_t getDataIndex(const size_t ix, const size_t iy, const size_t iz) const
  {
    const size_t indX = ix;
    const size_t indY = iy;
    const size_t indZ = iz;
    const size_t index = indX + Nx * (indY + indZ * Ny);
    return index;
  }

  static constexpr size_t getNDataPoints()
  {
    return FNdataPoints;
  }

  int writeToFile(TFile& outf, const char* name = "data") const
  {

    if (outf.IsZombie()) {
      std::cout << "Failed to write to file " << outf.GetName() << std::endl;
      return -1;
    }

    setStreamer();
    outf.WriteObjectAny(this, DataContainer3D<DataT, Nx, Ny, Nz>::Class(), name);
    return 0;
  }

  inline static DataContainer3D<DataT, Nx, Ny, Nz>* readFromFile(TFile& inpf, const char* name = "data")
  {
    /// read from file
    if (inpf.IsZombie()) {
      std::cout << "Failed to read from file " << inpf.GetName() << std::endl;
      return nullptr;
    }
    DataContainer3D<DataT, Nx, Ny, Nz>* dataCont{nullptr};

    dataCont->setStreamer();
    dataCont = reinterpret_cast<DataContainer3D<DataT, Nx, Ny, Nz>*>(inpf.GetObjectChecked(name, DataContainer3D<DataT, Nx, Ny, Nz>::Class()));
    if (!dataCont) {
      std::cout << "Failed to load " << name << " from " << inpf.GetName() << std::endl;
      return nullptr;
    }

    return dataCont;
  }

 private:
  static constexpr void setStreamer()
  {
    auto* tClass = DataContainer3D<DataT, Nx, Ny, Nz>::Class();
    const char* className = DataContainer3D<DataT, Nx, Ny, Nz>::Class()->GetName();
    if (tClass) {
      tClass->SetStreamerFunc(dataStreamer);
    } else {
      const char* errMsg{Form("Streamer: dictionary for following class is not available: %s", className)};
      throw(std::runtime_error(errMsg));
    }
  }

  //custom data streamer for static array
  static constexpr void dataStreamer(TBuffer& buf, void* objPtr)
  {
    DataContainer3D<DataT, Nx, Ny, Nz>* dataCont = static_cast<DataContainer3D<DataT, Nx, Ny, Nz>*>(objPtr);
    if (buf.IsReading()) {
      // buf >> dataCont->mem;
      buf.ReadFastArray(dataCont->mData.get(), FNdataPoints);
    } else {
      // buf << dataCont->mem;
      buf.WriteFastArray(dataCont->mData.get(), FNdataPoints);
    }
  }
  ClassDefNV(DataContainer3D, 1)
};

template <typename DataT = float, unsigned int Nx = 4, unsigned int Ny = 4, unsigned int Nz = 4>
struct RegularGrid3D {

 public:
  RegularGrid3D(const DataT* gridData, const DataT xmin, const DataT ymin, const DataT zmin, const DataT spacingX, const DataT spacingY, const DataT spacingZ) : mMin{{xmin, ymin, zmin}}, mInvSpacing{{static_cast<DataT>(1 / spacingX), static_cast<DataT>(1 / spacingY), static_cast<DataT>(1 / spacingZ)}}
  {
    initArray(gridData);
    initLists(spacingX, spacingY, spacingZ);
  }

  RegularGrid3D(const DataT xmin = 0, const DataT ymin = 0, const DataT zmin = 0, const DataT spacingX = 1, const DataT spacingY = 1, const DataT spacingZ = 1) : mMin{{xmin, ymin, zmin}}, mInvSpacing{{static_cast<DataT>(1 / spacingX), static_cast<DataT>(1 / spacingY), static_cast<DataT>(1 / spacingZ)}}
  {
    initLists(spacingX, spacingY, spacingZ);
  }

  RegularGrid3D(TFile& inpf, const char* name = "data", const DataT xmin = 0, const DataT ymin = 0, const DataT zmin = 0, const DataT spacingX = 1, const DataT spacingY = 1, const DataT spacingZ = 1) : mMin{{xmin, ymin, zmin}}, mInvSpacing{{static_cast<DataT>(1 / spacingX), static_cast<DataT>(1 / spacingY), static_cast<DataT>(1 / spacingZ)}}
  {
    initLists(spacingX, spacingY, spacingZ);
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
    // const size_t indX = ix;
    // const size_t indY = iy;
    // const size_t indZ = iz;
    // const size_t index = indX + Nx * (indY + indZ * Ny);
    // return index;
    return mGridData.getDataIndex(ix, iy, iz);
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

  const DataT& operator[](size_t i) const { return mGridData[i]; }
  DataT& operator[](size_t i) { return mGridData[i]; }

  const DataT& operator()(size_t ix, size_t iy, size_t iz) const
  {
    // const size_t ind = getDataIndex(ix, iy, iz);
    return mGridData(ix, iy, iz);
  }

  DataT& operator()(size_t ix, size_t iy, size_t iz)
  {
    // const size_t ind = getDataIndex(ix, iy, iz);
    return mGridData(ix, iy, iz);
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

  DataT getXVertex(const size_t index) const
  {
    return listXVertices[index];
  }

  DataT getYVertex(const size_t index) const
  {
    return listYVertices[index];
  }

  DataT getZVertex(const size_t index) const
  {
    return listZVertices[index];
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

  const Vector<DataT, FDim> mMin{};                              //! min positions of grid
  const Vector<DataT, FDim> mInvSpacing{};                       //! inverse spacing of grid
  const Vector<DataT, FDim> mMaxIndex{{Nx - 1, Ny - 1, Nz - 1}}; //! max index which is on the grid in all dimensions
  static constexpr size_t FNdim[FDim]{Nx, Ny, Nz};

  DataT listXVertices[Nx]{};
  DataT listYVertices[Ny]{};
  DataT listZVertices[Nz]{};

  void initLists(const DataT spacingX, const DataT spacingY, const DataT spacingZ)
  {

    for (size_t i = 0; i < Nx; ++i) {
      listXVertices[i] = getGridMinX() + i * spacingX;
    }

    for (size_t i = 0; i < Ny; ++i) {
      listYVertices[i] = getGridMinY() + i * spacingY;
    }

    for (size_t i = 0; i < Nz; ++i) {
      listZVertices[i] = getGridMinZ() + i * spacingZ;
    }
  }
  ClassDefNV(RegularGrid3D, 1)
};

// template <typename DataT, unsigned int Nx, unsigned int Ny, unsigned int Nz>
template <template <typename, unsigned int, unsigned int, unsigned int> typename MClass, typename DataT, unsigned int Nx, unsigned int Ny, unsigned int Nz>
// std::ostream& operator<<(std::ostream& out, const RegularGrid3D<DataT, Nx, Ny, Nz>& cont)
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
