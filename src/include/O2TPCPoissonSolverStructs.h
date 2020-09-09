
#include <iostream>
#include <iomanip>

namespace o2
{
namespace tpc
{

/// \tparam DataT the type of storage which is used during the calculations
template <typename DataT = float>
struct Matrix3D {

    Matrix3D(const unsigned int nRTmp, const unsigned int nZTmp, const unsigned int nPhiTmp) : nR{nZTmp}, nZ{nRTmp}, nPhi{nPhiTmp}
  {
    storage.resize(nRTmp * nZTmp * nPhiTmp);
  };

  Matrix3D(){};

  DataT& operator()(const unsigned int iR, const unsigned int iZ, const unsigned int iPhi)
  {
    return storage[iR + nR * (iZ + nZ * iPhi)];
  }

  const DataT& operator()(const unsigned int iR, const unsigned int iZ, const unsigned int iPhi) const
  {
    return storage[iR + nR * (iZ + nZ * iPhi)];
  }

  void resize(const unsigned int nRTmp, const unsigned int nZTmp, const unsigned int nPhiTmp)
  {
    nR = nRTmp;
    nZ = nZTmp;
    nPhi = nPhiTmp;
    storage.resize(nRTmp * nZTmp * nPhiTmp);
  }

  const DataT* data() const {
    return storage.data();
  }

  DataT* data(){
    return storage.data();
  }

  void print(const int nZStart = 0, int nZMax = -1)
  {
    nZMax = nZMax == -1 ? nPhi - 1 : nZMax;
    std::ostream& out = std::cout;
    out.precision(3);
    auto&& w = std::setw(9);
    out << std::endl;

    for (unsigned int iPhi = nZStart; iPhi <= nZMax; ++iPhi) {
      out << "z layer: " << iPhi << std::endl;
      // print top x row
      out << "⎡" << w << storage[0 + nR * (0 + nZ * iPhi)];
      for (unsigned int iR = 1; iR < nR; ++iR) {
        out << ", " << w << storage[iR + nR * (0 + nZ * iPhi)];
      }
      out << " ⎤" << std::endl;

      for (unsigned int iZ = 1; iZ < nZ - 1; ++iZ) {
        out << "⎢" << w << storage[0 + nR * (iZ + nZ * iPhi)];
        for (unsigned int iR = 1; iR < nR; ++iR) {
          out << ", " << w << storage[iR + nR * (iZ + nZ * iPhi)];
        }
        out << " ⎥" << std::endl;
      }

      out << "⎣" << w << storage[0 + nR * ((nZ - 1) + nZ * iPhi)];
      for (unsigned int iR = 1; iR < nR; ++iR) {
        out << ", " << w << storage[iR + nR * ((nZ - 1) + nZ * iPhi)];
      }
      out << " ⎦" << std::endl;
      out << std::endl;
      out << std::endl;
    }
  }

  unsigned int nR{};
  unsigned int nZ{};
  unsigned int nPhi{};
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
  inline static bool isFull3D = false;             ///<  TRUE: full coarsening, FALSE: semi coarsening
  inline static CycleType cycleType = FCycle;      ///< cycleType follow  CycleType
  inline static GridTransferType gtType = Full;    ///< gtType grid transfer type follow GridTransferType
  inline static RelaxType relaxType = GaussSeidel; ///< relaxType follow RelaxType
  inline static int nPre = 2;                      ///< number of iteration for pre smoothing
  inline static int nPost = 2;                     ///< number of iteration for post smoothing
  inline static int nMGCycle = 200;                ///< number of multi grid cycle (V type)
  inline static int maxLoop = 6;                   ///< the number of tree-deep of multi grid
};

template <typename DataT = float>
struct TPCParameters {
  static constexpr DataT TPCZ0{249.7};                          ///< nominal gating grid position
  static constexpr DataT IFCRADIUS{83.5};                       ///< Mean Radius of the Inner Field Cage ( 82.43 min,  83.70 max) (cm)
  static constexpr DataT OFCRADIUS{254.5};                      ///< Mean Radius of the Outer Field Cage (252.55 min, 256.45 max) (cm)
  static constexpr DataT ZOFFSET{0.2};                          ///< Offset from CE: calculate all distortions closer to CE as if at this point
  static constexpr DataT CATHODEV{-100000.0};                   ///< Cathode Voltage (volts)
  static constexpr DataT GG{-70.0};                             ///< Gating Grid voltage (volts)
  static constexpr DataT DVDE{0.0024};                          ///< [cm/V] drift velocity dependency on the E field (from Magboltz for NeCO2N2 at standard environment)
  static constexpr DataT EM{-1.602176487e-19 / 9.10938215e-31}; ///< charge/mass in [C/kg]
  static constexpr DataT E0{8.854187817e-12};                   ///< vacuum permittivity [A·s/(V·m)]
};

} // namespace tpc
} // namespace o2
