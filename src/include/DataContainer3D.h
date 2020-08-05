
#ifndef DataContainer3D_H
#define DataContainer3D_H

#include <memory>
#include <iostream>
#include "TFile.h"

/// \tparam DataT the type of data which is used during the calculations
/// \tparam Nx number of vertices in x direction
/// \tparam Ny number of vertices in y direction
/// \tparam Nz number of vertices in z direction
template <typename DataT = float, unsigned int Nx = 4, unsigned int Ny = 4, unsigned int Nz = 4>
struct DataContainer3D {

  static constexpr size_t FNdataPoints{Nx * Ny * Nz};                       ///< number of values stored in the container
  std::unique_ptr<DataT[]> mData = std::make_unique<DataT[]>(FNdataPoints); ///< storage for the data

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
    const size_t index = ix + Nx * (iy + iz * Ny);
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

  //custom data streamer for unique_ptr
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

#endif
