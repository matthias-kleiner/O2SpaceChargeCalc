
#ifndef DataContainer3D_H
#define DataContainer3D_H

#include <memory>
#include <iostream>
#include "TFile.h"

/// \tparam DataT the type of data which is used during the calculations
/// \tparam Nx number of values in x direction
/// \tparam Ny number of values in y direction
/// \tparam Nz number of values in z direction
template <typename DataT = float, unsigned int Nx = 4, unsigned int Ny = 4, unsigned int Nz = 4>
struct DataContainer3D {

  DataContainer3D() = default;
  
 public:
  class iterator
  {
   public:

    iterator(DataT* ptr) : ptr(ptr) {}
    iterator operator++()
    {
      ++ptr;
      return *this;
    }
    bool operator!=(const iterator& other) const { return ptr != other.ptr; }
    const DataT& operator*() const { return *ptr; }

   private:
    DataT* ptr;
  };

  iterator begin() const { return iterator(mData.get()); }
  iterator end() const { return iterator(mData.get() + FN); }

  static constexpr size_t FN{Nx * Ny * Nz}; ///< number of values stored in the container

  /// \param i index too data
  const DataT& operator[](size_t i) const { return mData[i]; }
  DataT& operator[](size_t i) { return mData[i]; }

  DataT* data() { return mData.get(); }

  /// \param ix index in x dimension
  /// \param iy index in y dimension
  /// \param iz index in z dimension
  const DataT& operator()(size_t ix, size_t iy, size_t iz) const
  {
    const size_t ind = getDataIndex(ix, iy, iz);
    return mData[ind];
  }

  /// \param ix index in x dimension
  /// \param iy index in y dimension
  /// \param iz index in z dimension
  DataT& operator()(size_t ix, size_t iy, size_t iz)
  {
    const size_t ind = getDataIndex(ix, iy, iz);
    return mData[ind];
  }

  /// \param ix index in x dimension
  /// \param iy index in y dimension
  /// \param iz index in z dimension
  /// \return returns the index to data
  size_t getDataIndex(const size_t ix, const size_t iy, const size_t iz) const
  {
    const size_t index = ix + Nx * (iy + iz * Ny);
    return index;
  }

  /// \return returns the number of values stored
  static constexpr size_t getNDataPoints()
  {
    return FN;
  }

  static constexpr size_t getNX()
  {
    return Nx;
  }

  static constexpr size_t getNY()
  {
    return Ny;
  }

  static constexpr size_t getNZ()
  {
    return Nz;
  }
  /// write object to file
  /// \param outf object is written to this file
  /// \param name object is save with this name
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

  /// set values from file
  bool initFromFile(TFile& inpf, const char* name = "data")
  {
    if (inpf.IsZombie()) {
      std::cout << "Failed to read from file " << inpf.GetName() << std::endl;
      return false;
    }
    DataContainer3D<DataT, Nx, Ny, Nz>* dataCont{nullptr};

    dataCont->setStreamer();
    dataCont = reinterpret_cast<DataContainer3D<DataT, Nx, Ny, Nz>*>(inpf.GetObjectChecked(name, DataContainer3D<DataT, Nx, Ny, Nz>::Class()));
    if (!dataCont) {
      std::cout << "Failed to load " << name << " from " << inpf.GetName() << std::endl;
      return false;
    }
    memcpy(mData.get(), dataCont->mData.get(), dataCont->getNDataPoints() * sizeof(DataT));
    return true;
  }

  /// get pointer to object from file
  inline static DataContainer3D<DataT, Nx, Ny, Nz>* loadFromFile(TFile& inpf, const char* name = "data")
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
  std::unique_ptr<DataT[]> mData = std::make_unique<DataT[]>(FN); ///< storage for the data
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
      buf.ReadFastArray(dataCont->mData.get(), FN);
    } else {
      // buf << dataCont->mem;
      buf.WriteFastArray(dataCont->mData.get(), FN);
    }
  }
  ClassDefNV(DataContainer3D, 1)
};

#endif
