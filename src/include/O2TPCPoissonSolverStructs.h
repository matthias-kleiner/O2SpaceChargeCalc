#include <iostream>
#include <iomanip>

/// \tparam DataT the type of storage which is used during the calculations
template <typename DataT = float>
struct Matrix3D {

  Matrix3D(const unsigned int nXTmp, const unsigned int nYTmp, const unsigned int nZTmp) : nX{nYTmp}, nY{nXTmp}, nZ{nZTmp}
  {
    storage.resize(nXTmp * nYTmp * nZTmp);
  };

  Matrix3D(){};

  // ix and iy must be swapped since in the original version TMatrixD where used! TMatrixD(row,column)
  DataT& operator()(const unsigned int iy, const unsigned int ix, const unsigned int iz)
  {
    return storage[ix + nX * (iy + nY * iz)];
  }

  const DataT& operator()(const unsigned int iy, const unsigned int ix, const unsigned int iz) const
  {
    return storage[ix + nX * (iy + nY * iz)];
  }

  void resize(const unsigned int nXTmp, const unsigned int nYTmp, const unsigned int nZTmp)
  {
    nX = nXTmp;
    nY = nYTmp;
    nZ = nZTmp;
    storage.resize(nXTmp * nYTmp * nZTmp);
  }

  void print(const int nZStart = 0, int nZMax = -1)
  {
    nZMax = nZMax == -1 ? nZ - 1 : nZMax;
    std::ostream& out = std::cout;
    out.precision(3);
    auto&& w = std::setw(9);
    out << std::endl;

    for (unsigned int iz = nZStart; iz <= nZMax; ++iz) {
      out << "z layer: " << iz << std::endl;
      // print top x row
      out << "⎡" << w << storage[0 + nX * (0 + nY * iz)];
      for (unsigned int ix = 1; ix < nX; ++ix) {
        out << ", " << w << storage[ix + nX * (0 + nY * iz)];
      }
      out << " ⎤" << std::endl;

      for (unsigned int iy = 1; iy < nY - 1; ++iy) {
        out << "⎢" << w << storage[0 + nX * (iy + nY * iz)];
        for (unsigned int ix = 1; ix < nX; ++ix) {
          out << ", " << w << storage[ix + nX * (iy + nY * iz)];
        }
        out << " ⎥" << std::endl;
      }

      out << "⎣" << w << storage[0 + nX * ((nY - 1) + nY * iz)];
      for (unsigned int ix = 1; ix < nX; ++ix) {
        out << ", " << w << storage[ix + nX * ((nY - 1) + nY * iz)];
      }
      out << " ⎦" << std::endl;
      out << std::endl;
      out << std::endl;
    }
  }

  unsigned int nX{};
  unsigned int nY{};
  unsigned int nZ{};
  std::vector<DataT> storage{};
};

///< Enumeration of Cycles Type
enum CycleType {
  VCycle = 0, ///< V Cycle
  WCycle = 1, ///< W Cycle (TODO)
  FCycle = 2  ///< Full Cycle
};

///< Fine -> Coarse Grid transfer operator types
enum GridTransferType {
  Half = 0, ///< Half weighting
  Full = 1, ///< Full weighting
};

///< Smoothing (Relax) operator types
enum RelaxType {
  Jacobi = 0,         ///< Jacobi (5 Stencil 2D, 7 Stencil 3D_
  WeightedJacobi = 1, ///< (TODO)
  GaussSeidel = 2     ///< Gauss Seidel 2D (2 Color, 5 Stencil), 3D (7 Stencil)
};

///< Parameters choice for MultiGrid    algorithm
struct MGParameters {
  bool isFull3D = false;             ///<  TRUE: full coarsening, FALSE: semi coarsening
  CycleType cycleType = FCycle;      ///< cycleType follow  CycleType
  GridTransferType gtType = Full;    ///< gtType grid transfer type follow GridTransferType
  RelaxType relaxType = GaussSeidel; ///< relaxType follow RelaxType
  // int gamma;                          ///< number of iteration at coarsest level
  int nPre = 2;       ///< number of iteration for pre smoothing
  int nPost = 2;      ///< number of iteration for post smoothing
  int nMGCycle = 200; ///< number of multi grid cycle (V type)
  int maxLoop = 6;    ///< the number of tree-deep of multi grid
};
