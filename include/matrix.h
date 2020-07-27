
#ifndef MATR_H
#define MATR_H

#include <Vc/Vc>

template <typename DataT = float, size_t N = 64>
class Matrix
{
  using VDataT = Vc::Vector<DataT>;

 public:
  Matrix(const Vc::Memory<VDataT, N>* dataMatrix) : mDataMatrix(dataMatrix) {}

  // const overload of the above
  const Vc::Memory<VDataT, N>& operator[](size_t i) const { return mDataMatrix[i]; }

 private:
  const Vc::Memory<VDataT, N>* mDataMatrix{};
};

template <typename DataT = float, size_t N = 64>
class Vector
{
  using VDataT = Vc::Vector<DataT>;

 public:
  Vector(const Vc::Memory<VDataT, N>& dataVector) : mDataVector(dataVector) {}

  Vector(){};

  Vector(const DataT dataArr, size_t size)
  {
    for (int i = 0; i < size; ++i) {
      mDataVector.scalar(i) = dataArr[i];
      // mDataVector.vector(i) += V(&dataArr[i], Vc::Aligned);
      // std::cout<<'test'<<std::endl;
    }
  };

  Vector(const DataT val)
  {
    for (size_t i = 0; i < N; ++i) {
      mDataVector.scalar(i) = val;
    }
  };

  const DataT operator[](size_t i) const { return mDataVector.scalar(i); }
  DataT& operator[](size_t i) { return mDataVector.scalar(i); }

  void setVector(const size_t j, const VDataT& vector)
  {
    mDataVector.vector(j) = vector;
  }

  const Vc::Vector<DataT> GetVector(const size_t i) const { return mDataVector.vector(i); }

  Vc::Memory<VDataT, N> GetMemory() const { return mDataVector; }

  size_t GetvectorsCount() const
  {
    return mDataVector.vectorsCount();
  }

  size_t GetentriesCount() const
  {
    return mDataVector.entriesCount();
  }

  Vc::Memory<VDataT, N>& GetDataStorage()
  {
    return mDataVector;
  }

 private:
  // storage for the data
  Vc::Memory<VDataT, N> mDataVector{};
};

template <typename T, size_t N>
inline Vector<T, N> operator*(const Matrix<T, N>& a, const Vector<T, N>& b)
{
  using V = Vc::Vector<T>;

  // resulting vector c
  Vector<T, N> c;
  for (size_t i = 0; i < N; ++i) {
    V c_ij{};
    for (size_t j = 0; j < a[i].vectorsCount(); ++j) {
      c_ij += a[i].vector(j) * b.GetVector(j);
    }
    c[i] = c_ij.sum();
  }
  return c;
}

template <typename T, size_t N>
inline Vector<T, N> operator-(const Vector<T, N>& a, const Vector<T, N>& b)
{
  // resulting matrix c
  Vector<T, N> c;
  for (size_t j = 0; j < a.GetvectorsCount(); ++j) {
    c.setVector(j, a.GetVector(j) - b.GetVector(j));
  }
  return c;
}

template <typename T, size_t N>
inline Vector<T, N> operator+(const Vector<T, N>& a, const Vector<T, N>& b)
{
  // resulting matrix c
  Vector<T, N> c;
  for (size_t j = 0; j < a.GetvectorsCount(); ++j) {
    c.setVector(j, a.GetVector(j) + b.GetVector(j));
  }
  return c;
}

template <typename T, size_t N>
inline Vector<T, N> operator*(const T a, const Vector<T, N>& b)
{
  // resulting matrix c
  Vector<T, N> c;
  for (size_t j = 0; j < b.GetvectorsCount(); ++j) {
    c.setVector(j, a * b.GetVector(j));
  }
  return c;
}

template <typename DataT, size_t N>
inline DataT sum(const Vector<DataT, N>& a)
{
  // resulting matrix c
  Vc::Vector<DataT> b = a.GetVector(0);

  for (size_t j = 1; j < a.GetvectorsCount(); ++j) {
    b += a.GetVector(j);
  }
  return b.sum();
}

// multiply each erow from a vector with the row from a second vector
template <typename DataT, size_t N>
inline Vector<DataT, N> operator*(const Vector<DataT, N>& a, const Vector<DataT, N>& b)
{
  // resulting matrix c
  Vector<DataT, N> c;
  for (size_t j = 0; j < a.GetvectorsCount(); ++j) {
    c.setVector(j, a.GetVector(j) * b.GetVector(j));
  }
  return c;
}

// check if all elements are equal
template <typename DataT, size_t N>
inline bool operator==(const Vector<DataT, N>& a, const Vector<DataT, N>& b)
{
  for (size_t j = 0; j < a.GetvectorsCount(); ++j) {
    Vc::Mask<DataT> c = a.GetVector(j) != b.GetVector(j);
    if (c.count() > 0) {
      return false;
    }
  }
  return true;
}

#endif
