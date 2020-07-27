#ifndef TRI_H
#define TRI_H

#include "matrix.h"
#include "RegularGrid3D.h"

template <typename DataT = float, typename Grid3D = RegularGrid3D<>>
class TriCubicInterpolator
{
  using VDataT = Vc::Vector<DataT>;

 public:
  /// Constructor for a tricubic interpolator
  /// \param gridData memory adress of the 3D-Grid struct
  TriCubicInterpolator(const Grid3D& gridData, const bool circularX = false, const bool circularY = false, const bool circularZ = false) : mGridData{gridData}, mCircularX(circularX), mCircularY(circularY), mCircularZ(circularZ){};

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

  // interpolate derivative at given coordinate
  /// \param x x coordinate
  /// \param y y coordinate
  /// \param z z coordinate
  /// \param derx order of derivative d/dx: derx=1 -> d/dx f(x,y,z), derx=2 -> d^2/dx^2 f(x,y,z), derx=3 -> d^3/dx^3 f(x,y,z)
  /// \param derz order of derivative d/dy: dery=1 -> d/dy f(x,y,z), dery=2 -> d^2/dy^2 f(x,y,z), dery=3 -> d^3/dy^3 f(x,y,z)
  /// \param derz order of derivative d/dz: derz=1 -> d/dz f(x,y,z), derz=2 -> d^2/dz^2 f(x,y,z), derz=3 -> d^3/dz^3 f(x,y,z)
  /// derx=1 and dery=2 -> d/dx * d^2/dy^2 * f(x,y,z)
  /// \param safe checks if the given coordinates lie in the grid. if a value is out of the grid, the coordinate will be set to the border of the grid
  /// \return returns the interpolated derivative at given coordinate
  DataT operator()(const DataT x, const DataT y, const DataT z, const int derx, const int dery, const int derz, const bool safe) const
  {
    const Vector<DataT, FDim> coordinates{{x, y, z}}; // vector holding the coordinates
    const auto relPos = processInp(coordinates, safe);
    return EvalDerivative(relPos[0], relPos[1], relPos[2], derx, dery, derz);
  }

  void initInterpolator(const unsigned int ix, const unsigned int iy, const unsigned int iz) const
  {
    // calcCoefficientsVc(ix, iy, iz);
    calcCoefficients(ix, iy, iz);

    // store current cell
    mInitialized = true;
    mLastInd[FX] = ix;
    mLastInd[FY] = iy;
    mLastInd[FZ] = iz;
  }

 private:
  DataT uiPow(const DataT base, unsigned int exp) const
  {
    DataT result = 1;
    // infinite for loop
    for (;;) {
      // check if x is uneven number
      if (exp & 1) {
        result *= base;
      }
      exp >>= 1;
      if (!exp) {
        break;
      }
      base *= base;
    }
    return result;
  }

  template <size_t FDim>
  const Vector<DataT, FDim> processInp(const Vector<DataT, FDim>& coordinates, const bool safe) const
  {
    const auto minGrid = mGridData.getGridMin();                      // vector containing the min x,y,z value of the grid
    const auto invSpacing = mGridData.getInvSpacing();                // vector containing the grid spacing for dim x,y,z
    Vector<DataT, FDim> posRel{(coordinates - minGrid) * invSpacing}; // needed for the grid index

    if (safe) {
      mGridData.clampToGrid(posRel);
    }

    const unsigned int ix = static_cast<unsigned int>(posRel[FX]);
    const unsigned int iy = static_cast<unsigned int>(posRel[FY]);
    const unsigned int iz = static_cast<unsigned int>(posRel[FZ]);

    const Vector<unsigned int, FDim> index{{ix, iy, iz}};

    if (!mInitialized || !(mLastInd == index)) {
      initInterpolator(index[FX], index[FY], index[FZ]);
    }

    const Vector<DataT, FDim> indexTmp{{static_cast<DataT>(ix), static_cast<DataT>(iy), static_cast<DataT>(iz)}};
    // const Vector<DataT, FDim> relPos{{posRel[FX] - ix, posRel[FY] - iy, posRel[FZ] - iz}};
    const Vector<DataT, FDim> relPos{posRel - indexTmp};
    return relPos;
  }

  void calcCoefficients(const unsigned int ix, const unsigned int iy, const unsigned int iz) const
  {
    Vector<DataT, 64> matrixPar{};
    int deltaX[4]{};
    int deltaY[4]{};
    int deltaZ[4]{};

    // set padding type: circular or standard
    mCircularX ? getDataIndexCircularArray(ix, FX, deltaX) : getDataIndexNonCircularArray(ix, FX, deltaX);
    mCircularY ? getDataIndexCircularArray(iy, FY, deltaY) : getDataIndexNonCircularArray(iy, FY, deltaY);
    mCircularZ ? getDataIndexNonCircularArray(iz, FZ, deltaZ) : getDataIndexNonCircularArray(iz, FZ, deltaZ);

    const unsigned int i_x_y_z = mGridData.getDataIndex(ix, iy, iz);
    const unsigned int i_xp1_y_z = i_x_y_z + deltaX[2];
    const unsigned int i_x_yp1_z = i_x_y_z + deltaY[2];
    const unsigned int i_xp1_yp1_z = i_x_y_z + deltaX[2] + deltaY[2];
    const unsigned int i_x_y_zp1 = i_x_y_z + deltaZ[2];
    const unsigned int i_xp1_y_zp1 = i_x_y_z + deltaX[2] + deltaZ[2];
    const unsigned int i_x_yp1_zp1 = i_x_y_z + deltaY[2] + deltaZ[2];
    const unsigned int i_xp1_yp1_zp1 = i_x_y_z + deltaX[2] + deltaY[2] + deltaZ[2];
    matrixPar[0] = mGridData[i_x_y_z];
    matrixPar[1] = mGridData[i_xp1_y_z];
    matrixPar[2] = mGridData[i_x_yp1_z];
    matrixPar[3] = mGridData[i_xp1_yp1_z];
    matrixPar[4] = mGridData[i_x_y_zp1];
    matrixPar[5] = mGridData[i_xp1_y_zp1];
    matrixPar[6] = mGridData[i_x_yp1_zp1];
    matrixPar[7] = mGridData[i_xp1_yp1_zp1];

    // values of df/dx at each corner.
    const DataT fac1 = 0.5;
    const unsigned int i_xp2_y_z = i_x_y_z + deltaX[3];
    const unsigned int i_xm1_y_z = i_x_y_z + deltaX[1];
    const unsigned int i_xm1_yp1_z = i_x_y_z + deltaX[1] + deltaY[2];
    const unsigned int i_xp2_yp1_z = i_x_y_z + deltaX[3] + deltaY[2];
    const unsigned int i_xm1_y_zp1 = i_x_y_z + deltaX[1] + deltaZ[2];
    const unsigned int i_xp2_y_zp1 = i_x_y_z + deltaX[3] + deltaZ[2];
    const unsigned int i_xm1_yp1_zp1 = i_x_y_z + deltaX[1] + deltaY[2] + deltaZ[2];
    const unsigned int i_xp2_yp1_zp1 = i_x_y_z + deltaX[3] + deltaY[2] + deltaZ[2];

    // load values to tmp Vc
    Vector<DataT, 24> vecDeriv1A{};
    Vector<DataT, 24> vecDeriv1B{};
    vecDeriv1A[0] = mGridData[i_xp1_y_z];
    vecDeriv1B[0] = mGridData[i_xm1_y_z];
    vecDeriv1A[1] = mGridData[i_xp2_y_z];
    vecDeriv1B[1] = mGridData[i_x_y_z];
    vecDeriv1A[2] = mGridData[i_xp1_yp1_z];
    vecDeriv1B[2] = mGridData[i_xm1_yp1_z];
    vecDeriv1A[3] = mGridData[i_xp2_yp1_z];
    vecDeriv1B[3] = mGridData[i_x_yp1_z];
    vecDeriv1A[4] = mGridData[i_xp1_y_zp1];
    vecDeriv1B[4] = mGridData[i_xm1_y_zp1];
    vecDeriv1A[5] = mGridData[i_xp2_y_zp1];
    vecDeriv1B[5] = mGridData[i_x_y_zp1];
    vecDeriv1A[6] = mGridData[i_xp1_yp1_zp1];
    vecDeriv1B[6] = mGridData[i_xm1_yp1_zp1];
    vecDeriv1A[7] = mGridData[i_xp2_yp1_zp1];
    vecDeriv1B[7] = mGridData[i_x_yp1_zp1];

    // values of df/dy at each corner.
    const unsigned int i_x_ym1_z = i_x_y_z + deltaY[1];
    const unsigned int i_xp1_ym1_z = i_x_y_z + deltaX[2] + deltaY[1];
    const unsigned int i_x_yp2_z = i_x_y_z + deltaY[3];
    const unsigned int i_xp1_yp2_z = i_x_y_z + deltaX[2] + deltaY[3];
    const unsigned int i_x_ym1_zp1 = i_x_y_z + deltaY[1] + deltaZ[2];
    const unsigned int i_xp1_ym1_zp1 = i_x_y_z + deltaX[2] + deltaY[1] + deltaZ[2];
    const unsigned int i_x_yp2_zp1 = i_x_y_z + deltaY[3] + deltaZ[2];
    const unsigned int i_xp1_yp2_zp1 = i_x_y_z + deltaX[2] + deltaY[3] + deltaZ[2];
    vecDeriv1A[8] = mGridData[i_x_yp1_z];
    vecDeriv1B[8] = mGridData[i_x_ym1_z];
    vecDeriv1A[9] = mGridData[i_xp1_yp1_z];
    vecDeriv1B[9] = mGridData[i_xp1_ym1_z];
    vecDeriv1A[10] = mGridData[i_x_yp2_z];
    vecDeriv1B[10] = mGridData[i_x_y_z];
    vecDeriv1A[11] = mGridData[i_xp1_yp2_z];
    vecDeriv1B[11] = mGridData[i_xp1_y_z];
    vecDeriv1A[12] = mGridData[i_x_yp1_zp1];
    vecDeriv1B[12] = mGridData[i_x_ym1_zp1];
    vecDeriv1A[13] = mGridData[i_xp1_yp1_zp1];
    vecDeriv1B[13] = mGridData[i_xp1_ym1_zp1];
    vecDeriv1A[14] = mGridData[i_x_yp2_zp1];
    vecDeriv1B[14] = mGridData[i_x_y_zp1];
    vecDeriv1A[15] = mGridData[i_xp1_yp2_zp1];
    vecDeriv1B[15] = mGridData[i_xp1_y_zp1];

    // values of df/dz at each corner.
    const unsigned int i_x_y_zm1 = i_x_y_z + deltaZ[1];
    const unsigned int i_xp1_y_zm1 = i_x_y_z + deltaX[2] + deltaZ[1];
    const unsigned int i_x_yp1_zm1 = i_x_y_z + deltaY[2] + deltaZ[1];
    const unsigned int i_xp1_yp1_zm1 = i_x_y_z + deltaX[2] + deltaY[2] + deltaZ[1];
    const unsigned int i_x_y_zp2 = i_x_y_z + deltaZ[3];
    const unsigned int i_xp1_y_zp2 = i_x_y_z + deltaX[2] + deltaZ[3];
    const unsigned int i_x_yp1_zp2 = i_x_y_z + deltaY[2] + deltaZ[3];
    const unsigned int i_xp1_yp1_zp2 = i_x_y_z + deltaX[2] + deltaY[2] + deltaZ[3];
    vecDeriv1A[16] = mGridData[i_x_y_zp1];
    vecDeriv1B[16] = mGridData[i_x_y_zm1];
    vecDeriv1A[17] = mGridData[i_xp1_y_zp1];
    vecDeriv1B[17] = mGridData[i_xp1_y_zm1];
    vecDeriv1A[18] = mGridData[i_x_yp1_zp1];
    vecDeriv1B[18] = mGridData[i_x_yp1_zm1];
    vecDeriv1A[19] = mGridData[i_xp1_yp1_zp1];
    vecDeriv1B[19] = mGridData[i_xp1_yp1_zm1];
    vecDeriv1A[20] = mGridData[i_x_y_zp2];
    vecDeriv1B[20] = mGridData[i_x_y_z];
    vecDeriv1A[21] = mGridData[i_xp1_y_zp2];
    vecDeriv1B[21] = mGridData[i_xp1_y_z];
    vecDeriv1A[22] = mGridData[i_x_yp1_zp2];
    vecDeriv1B[22] = mGridData[i_x_yp1_z];
    vecDeriv1A[23] = mGridData[i_xp1_yp1_zp2];
    vecDeriv1B[23] = mGridData[i_xp1_yp1_z];

    Vector<DataT, 24> vecDeriv2A{};
    Vector<DataT, 24> vecDeriv2B{};
    Vector<DataT, 24> vecDeriv2C{};
    Vector<DataT, 24> vecDeriv2D{};

    // values of d2f/dxdy at each corner.
    const DataT fac2 = 0.25;
    const unsigned int i_xm1_ym1_z = i_x_y_z + deltaX[1] + deltaY[1];
    const unsigned int i_xp2_ym1_z = i_x_y_z + deltaX[3] + deltaY[1];
    const unsigned int i_xm1_yp2_z = i_x_y_z + deltaX[1] + deltaY[3];
    const unsigned int i_xp2_yp2_z = i_x_y_z + deltaX[3] + deltaY[3];
    const unsigned int i_xm1_ym1_zp1 = i_x_y_z + deltaX[1] + deltaY[1] + deltaZ[2];
    const unsigned int i_xp2_ym1_zp1 = i_x_y_z + deltaX[3] + deltaY[1] + deltaZ[2];
    const unsigned int i_xm1_yp2_zp1 = i_x_y_z + deltaX[1] + deltaY[3] + deltaZ[2];
    const unsigned int i_xp2_yp2_zp1 = i_x_y_z + deltaX[3] + deltaY[3] + deltaZ[2];
    vecDeriv2A[0] = mGridData[i_xp1_yp1_z];
    vecDeriv2B[0] = mGridData[i_xm1_yp1_z];
    vecDeriv2C[0] = mGridData[i_xp1_ym1_z];
    vecDeriv2D[0] = mGridData[i_xm1_ym1_z];
    vecDeriv2A[1] = mGridData[i_xp2_yp1_z];
    vecDeriv2B[1] = mGridData[i_x_yp1_z];
    vecDeriv2C[1] = mGridData[i_xp2_ym1_z];
    vecDeriv2D[1] = mGridData[i_x_ym1_z];
    vecDeriv2A[2] = mGridData[i_xp1_yp2_z];
    vecDeriv2B[2] = mGridData[i_xm1_yp2_z];
    vecDeriv2C[2] = mGridData[i_xp1_y_z];
    vecDeriv2D[2] = mGridData[i_xm1_y_z];
    vecDeriv2A[3] = mGridData[i_xp2_yp2_z];
    vecDeriv2B[3] = mGridData[i_x_yp2_z];
    vecDeriv2C[3] = mGridData[i_xp2_y_z];
    vecDeriv2D[3] = mGridData[i_x_y_z];
    vecDeriv2A[4] = mGridData[i_xp1_yp1_zp1];
    vecDeriv2B[4] = mGridData[i_xm1_yp1_zp1];
    vecDeriv2C[4] = mGridData[i_xp1_ym1_zp1];
    vecDeriv2D[4] = mGridData[i_xm1_ym1_zp1];
    vecDeriv2A[5] = mGridData[i_xp2_yp1_zp1];
    vecDeriv2B[5] = mGridData[i_x_yp1_zp1];
    vecDeriv2C[5] = mGridData[i_xp2_ym1_zp1];
    vecDeriv2D[5] = mGridData[i_x_ym1_zp1];
    vecDeriv2A[6] = mGridData[i_xp1_yp2_zp1];
    vecDeriv2B[6] = mGridData[i_xm1_yp2_zp1];
    vecDeriv2C[6] = mGridData[i_xp1_y_zp1];
    vecDeriv2D[6] = mGridData[i_xm1_y_zp1];
    vecDeriv2A[7] = mGridData[i_xp2_yp2_zp1];
    vecDeriv2B[7] = mGridData[i_x_yp2_zp1];
    vecDeriv2C[7] = mGridData[i_xp2_y_zp1];
    vecDeriv2D[7] = mGridData[i_x_y_zp1];

    // values of d2f/dxdz at each corner.
    const unsigned int i_xm1_y_zm1 = i_x_y_z + deltaX[1] + deltaZ[1];
    const unsigned int i_xp2_y_zm1 = i_x_y_z + deltaX[3] + deltaZ[1];
    const unsigned int i_xm1_yp1_zm1 = i_x_y_z + deltaX[1] + deltaY[2] + deltaZ[1];
    const unsigned int i_xp2_yp1_zm1 = i_x_y_z + deltaX[3] + deltaY[2] + deltaZ[1];
    const unsigned int i_xm1_y_zp2 = i_x_y_z + deltaX[1] + deltaZ[3];
    const unsigned int i_xp2_y_zp2 = i_x_y_z + deltaX[3] + deltaZ[3];
    const unsigned int i_xm1_yp1_zp2 = i_x_y_z + deltaX[1] + deltaY[2] + deltaZ[3];
    const unsigned int i_xp2_yp1_zp2 = i_x_y_z + deltaX[3] + deltaY[2] + deltaZ[3];

    vecDeriv2A[8] = mGridData[i_xp1_y_zp1];
    vecDeriv2B[8] = mGridData[i_xm1_y_zp1];
    vecDeriv2C[8] = mGridData[i_xp1_y_zm1];
    vecDeriv2D[8] = mGridData[i_xm1_y_zm1];
    vecDeriv2A[9] = mGridData[i_xp2_y_zp1];
    vecDeriv2B[9] = mGridData[i_x_y_zp1];
    vecDeriv2C[9] = mGridData[i_xp2_y_zm1];
    vecDeriv2D[9] = mGridData[i_x_y_zm1];
    vecDeriv2A[10] = mGridData[i_xp1_yp1_zp1];
    vecDeriv2B[10] = mGridData[i_xm1_yp1_zp1];
    vecDeriv2C[10] = mGridData[i_xp1_yp1_zm1];
    vecDeriv2D[10] = mGridData[i_xm1_yp1_zm1];
    vecDeriv2A[11] = mGridData[i_xp2_yp1_zp1];
    vecDeriv2B[11] = mGridData[i_x_yp1_zp1];
    vecDeriv2C[11] = mGridData[i_xp2_yp1_zm1];
    vecDeriv2D[11] = mGridData[i_x_yp1_zm1];
    vecDeriv2A[12] = mGridData[i_xp1_y_zp2];
    vecDeriv2B[12] = mGridData[i_xm1_y_zp2];
    vecDeriv2C[12] = mGridData[i_xp1_y_z];
    vecDeriv2D[12] = mGridData[i_xm1_y_z];
    vecDeriv2A[13] = mGridData[i_xp2_y_zp2];
    vecDeriv2B[13] = mGridData[i_x_y_zp2];
    vecDeriv2C[13] = mGridData[i_xp2_y_z];
    vecDeriv2D[13] = mGridData[i_x_y_z];
    vecDeriv2A[14] = mGridData[i_xp1_yp1_zp2];
    vecDeriv2B[14] = mGridData[i_xm1_yp1_zp2];
    vecDeriv2C[14] = mGridData[i_xp1_yp1_z];
    vecDeriv2D[14] = mGridData[i_xm1_yp1_z];
    vecDeriv2A[15] = mGridData[i_xp2_yp1_zp2];
    vecDeriv2B[15] = mGridData[i_x_yp1_zp2];
    vecDeriv2C[15] = mGridData[i_xp2_yp1_z];
    vecDeriv2D[15] = mGridData[i_x_yp1_z];

    // values of d2f/dydz at each corner.
    const unsigned int i_x_ym1_zm1 = i_x_y_z + deltaY[1] + deltaZ[1];
    const unsigned int i_xp1_ym1_zm1 = i_x_y_z + deltaX[2] + deltaY[1] + deltaZ[1];
    const unsigned int i_x_yp2_zm1 = i_x_y_z + deltaY[3] + deltaZ[1];
    const unsigned int i_xp1_yp2_zm1 = i_x_y_z + deltaX[2] + deltaY[3] + deltaZ[1];
    const unsigned int i_x_ym1_zp2 = i_x_y_z + deltaY[1] + deltaZ[3];
    const unsigned int i_xp1_ym1_zp2 = i_x_y_z + deltaX[2] + deltaY[1] + deltaZ[3];
    const unsigned int i_x_yp2_zp2 = i_x_y_z + deltaY[3] + deltaZ[3];
    const unsigned int i_xp1_yp2_zp2 = i_x_y_z + deltaX[2] + deltaY[3] + deltaZ[3];

    vecDeriv2A[16] = mGridData[i_x_yp1_zp1];
    vecDeriv2B[16] = mGridData[i_x_ym1_zp1];
    vecDeriv2C[16] = mGridData[i_x_yp1_zm1];
    vecDeriv2D[16] = mGridData[i_x_ym1_zm1];
    vecDeriv2A[17] = mGridData[i_xp1_yp1_zp1];
    vecDeriv2B[17] = mGridData[i_xp1_ym1_zp1];
    vecDeriv2C[17] = mGridData[i_xp1_yp1_zm1];
    vecDeriv2D[17] = mGridData[i_xp1_ym1_zm1];
    vecDeriv2A[18] = mGridData[i_x_yp2_zp1];
    vecDeriv2B[18] = mGridData[i_x_y_zp1];
    vecDeriv2C[18] = mGridData[i_x_yp2_zm1];
    vecDeriv2D[18] = mGridData[i_x_y_zm1];
    vecDeriv2A[19] = mGridData[i_xp1_yp2_zp1];
    vecDeriv2B[19] = mGridData[i_xp1_y_zp1];
    vecDeriv2C[19] = mGridData[i_xp1_yp2_zm1];
    vecDeriv2D[19] = mGridData[i_xp1_y_zm1];
    vecDeriv2A[20] = mGridData[i_x_yp1_zp2];
    vecDeriv2B[20] = mGridData[i_x_ym1_zp2];
    vecDeriv2C[20] = mGridData[i_x_yp1_z];
    vecDeriv2D[20] = mGridData[i_x_ym1_z];
    vecDeriv2A[21] = mGridData[i_xp1_yp1_zp2];
    vecDeriv2B[21] = mGridData[i_xp1_ym1_zp2];
    vecDeriv2C[21] = mGridData[i_xp1_yp1_z];
    vecDeriv2D[21] = mGridData[i_xp1_ym1_z];
    vecDeriv2A[22] = mGridData[i_x_yp2_zp2];
    vecDeriv2B[22] = mGridData[i_x_y_zp2];
    vecDeriv2C[22] = mGridData[i_x_yp2_z];
    vecDeriv2D[22] = mGridData[i_x_y_z];

    // values of d3f/dxdydz at each corner.
    const float fac3 = 0.125;
    Vector<DataT, 8> vecDeriv3A{};
    Vector<DataT, 8> vecDeriv3B{};
    Vector<DataT, 8> vecDeriv3C{};
    Vector<DataT, 8> vecDeriv3D{};
    Vector<DataT, 8> vecDeriv3E{};
    Vector<DataT, 8> vecDeriv3F{};
    Vector<DataT, 8> vecDeriv3G{};
    Vector<DataT, 8> vecDeriv3H{};

    const unsigned int i_xm1_ym1_zm1 = i_x_y_z + deltaX[1] + deltaY[1] + deltaZ[1];
    const unsigned int i_xp2_ym1_zm1 = i_x_y_z + deltaX[3] + deltaY[1] + deltaZ[1];
    const unsigned int i_xm1_yp2_zm1 = i_x_y_z + deltaX[1] + deltaY[3] + deltaZ[1];
    const unsigned int i_xp2_yp2_zm1 = i_x_y_z + deltaX[3] + deltaY[3] + deltaZ[1];
    const unsigned int i_xm1_ym1_zp2 = i_x_y_z + deltaX[1] + deltaY[1] + deltaZ[3];
    const unsigned int i_xp2_ym1_zp2 = i_x_y_z + deltaX[3] + deltaY[1] + deltaZ[3];
    const unsigned int i_xm1_yp2_zp2 = i_x_y_z + deltaX[1] + deltaY[3] + deltaZ[3];
    const unsigned int i_xp2_yp2_zp2 = i_x_y_z + deltaX[3] + deltaY[3] + deltaZ[3];
    vecDeriv3A[0] = mGridData[i_xp1_yp1_zp1];
    vecDeriv3B[0] = mGridData[i_xm1_yp1_zp1];
    vecDeriv3C[0] = mGridData[i_xp1_ym1_zp1];
    vecDeriv3D[0] = mGridData[i_xm1_ym1_zp1];
    vecDeriv3E[0] = mGridData[i_xp1_yp1_zm1];
    vecDeriv3F[0] = mGridData[i_xm1_yp1_zm1];
    vecDeriv3G[0] = mGridData[i_xp1_ym1_zm1];
    vecDeriv3H[0] = mGridData[i_xm1_ym1_zm1];
    vecDeriv3A[1] = mGridData[i_xp2_yp1_zp1];
    vecDeriv3B[1] = mGridData[i_x_yp1_zp1];
    vecDeriv3C[1] = mGridData[i_xp2_ym1_zp1];
    vecDeriv3D[1] = mGridData[i_x_ym1_zp1];
    vecDeriv3E[1] = mGridData[i_xp2_yp1_zm1];
    vecDeriv3F[1] = mGridData[i_x_yp1_zm1];
    vecDeriv3G[1] = mGridData[i_xp2_ym1_zm1];
    vecDeriv3H[1] = mGridData[i_x_ym1_zm1];
    vecDeriv3A[2] = mGridData[i_xp1_yp2_zp1];
    vecDeriv3B[2] = mGridData[i_xm1_yp2_zp1];
    vecDeriv3C[2] = mGridData[i_xp1_y_zp1];
    vecDeriv3D[2] = mGridData[i_xm1_y_zp1];
    vecDeriv3E[2] = mGridData[i_xp1_yp2_zm1];
    vecDeriv3F[2] = mGridData[i_xm1_yp2_zm1];
    vecDeriv3G[2] = mGridData[i_xp1_y_zm1];
    vecDeriv3H[2] = mGridData[i_xm1_y_zm1];
    vecDeriv3A[3] = mGridData[i_xp2_yp2_zp1];
    vecDeriv3B[3] = mGridData[i_x_yp2_zp1];
    vecDeriv3C[3] = mGridData[i_xp2_y_zp1];
    vecDeriv3D[3] = mGridData[i_x_y_zp1];
    vecDeriv3E[3] = mGridData[i_xp2_yp2_zm1];
    vecDeriv3F[3] = mGridData[i_x_yp2_zm1];
    vecDeriv3G[3] = mGridData[i_xp2_y_zm1];
    vecDeriv3H[3] = mGridData[i_x_y_zm1];
    vecDeriv3A[4] = mGridData[i_xp1_yp1_zp2];
    vecDeriv3B[4] = mGridData[i_xm1_yp1_zp2];
    vecDeriv3C[4] = mGridData[i_xp1_ym1_zp2];
    vecDeriv3D[4] = mGridData[i_xm1_ym1_zp2];
    vecDeriv3E[4] = mGridData[i_xp1_yp1_z];
    vecDeriv3F[4] = mGridData[i_xm1_yp1_z];
    vecDeriv3G[4] = mGridData[i_xp1_ym1_z];
    vecDeriv3H[4] = mGridData[i_xm1_ym1_z];
    vecDeriv3A[5] = mGridData[i_xp2_yp1_zp2];
    vecDeriv3B[5] = mGridData[i_x_yp1_zp2];
    vecDeriv3C[5] = mGridData[i_xp2_ym1_zp2];
    vecDeriv3D[5] = mGridData[i_x_ym1_zp2];
    vecDeriv3E[5] = mGridData[i_xp2_yp1_z];
    vecDeriv3F[5] = mGridData[i_x_yp1_z];
    vecDeriv3G[5] = mGridData[i_xp2_ym1_z];
    vecDeriv3H[5] = mGridData[i_x_ym1_z];
    vecDeriv3A[6] = mGridData[i_xp1_yp2_zp2];
    vecDeriv3B[6] = mGridData[i_xm1_yp2_zp2];
    vecDeriv3C[6] = mGridData[i_xp1_y_zp2];
    vecDeriv3D[6] = mGridData[i_xm1_y_zp2];
    vecDeriv3E[6] = mGridData[i_xp1_yp2_z];
    vecDeriv3F[6] = mGridData[i_xm1_yp2_z];
    vecDeriv3G[6] = mGridData[i_xp1_y_z];
    vecDeriv3H[6] = mGridData[i_xm1_y_z];
    vecDeriv3A[7] = mGridData[i_xp2_yp2_zp2];
    vecDeriv3B[7] = mGridData[i_x_yp2_zp2];
    vecDeriv3C[7] = mGridData[i_xp2_y_zp2];
    vecDeriv3D[7] = mGridData[i_x_y_zp2];
    vecDeriv3E[7] = mGridData[i_xp2_yp2_z];
    vecDeriv3F[7] = mGridData[i_x_yp2_z];
    vecDeriv3G[7] = mGridData[i_xp2_y_z];
    vecDeriv3H[7] = mGridData[i_x_y_z];

    const Vector<DataT, 24> vfac1{
      {fac1, fac1, fac1, fac1, fac1, fac1, fac1, fac1,
       fac1, fac1, fac1, fac1, fac1, fac1, fac1, fac1,
       fac1, fac1, fac1, fac1, fac1, fac1, fac1, fac1,
       fac1, fac1, fac1, fac1, fac1, fac1, fac1, fac1}};

    const Vector<DataT, 24> vfac2{
      {fac2, fac2, fac2, fac2, fac2, fac2, fac2, fac2,
       fac2, fac2, fac2, fac2, fac2, fac2, fac2, fac2,
       fac2, fac2, fac2, fac2, fac2, fac2, fac2, fac2}};

    const Vector<DataT, 8> vfac3{{fac3, fac3, fac3, fac3, fac3, fac3, fac3, fac3}};

    const Vector<DataT, 24> vecDeriv1Res{vfac1 * (vecDeriv1A - vecDeriv1B)};
    const Vector<DataT, 24> vecDeriv2Res{vfac2 * (vecDeriv2A - vecDeriv2B - vecDeriv2C + vecDeriv2D)};
    const Vector<DataT, 8> vecDeriv3Res{vfac3 * (vecDeriv3A - vecDeriv3B - vecDeriv3C + vecDeriv3D - vecDeriv3E + vecDeriv3F + vecDeriv3G - vecDeriv3H)};

    for (size_t i = 0; i < vecDeriv1Res.GetentriesCount(); ++i) {
      matrixPar[8 + i] = vecDeriv1Res[i];
      matrixPar[32 + i] = vecDeriv2Res[i];
    }
    for (size_t i = 0; i < vecDeriv3Res.GetentriesCount(); ++i) {
      matrixPar[56 + i] = vecDeriv3Res[i];
    }
    mCoefficients = mMatrixA * matrixPar;
  }

  DataT interpolate(const Vector<DataT, 3>& pos) const
  {
    // the formula for evaluating the interpolation is as follows:
    // f(x,y,z) = \sum_{i,j,k=0}^3 a_{ijk} * x^{i} * y^{j} * z^{k}
    // memValX contains the x^{i} values

    const Vector<DataT, FDim> vals0{{1, 1, 1}};
    const Vector<DataT, FDim> vals2{pos * pos};
    const Vector<DataT, FDim> vals3{vals2 * pos};

    const DataT valX[4]{vals0[FX], pos[FX], vals2[FX], vals3[FX]};
    const Vector<DataT, 64> vecValX{{valX[0], valX[1], valX[2], valX[3],
                                     valX[0], valX[1], valX[2], valX[3],
                                     valX[0], valX[1], valX[2], valX[3],
                                     valX[0], valX[1], valX[2], valX[3],
                                     valX[0], valX[1], valX[2], valX[3],
                                     valX[0], valX[1], valX[2], valX[3],
                                     valX[0], valX[1], valX[2], valX[3],
                                     valX[0], valX[1], valX[2], valX[3],
                                     valX[0], valX[1], valX[2], valX[3],
                                     valX[0], valX[1], valX[2], valX[3],
                                     valX[0], valX[1], valX[2], valX[3],
                                     valX[0], valX[1], valX[2], valX[3],
                                     valX[0], valX[1], valX[2], valX[3],
                                     valX[0], valX[1], valX[2], valX[3],
                                     valX[0], valX[1], valX[2], valX[3],
                                     valX[0], valX[1], valX[2], valX[3]}};

    const DataT valY[4]{vals0[FY], pos[FY], vals2[FY], vals3[FY]};
    const Vector<DataT, 64> vecValY{{valY[0], valY[0], valY[0], valY[0],
                                     valY[1], valY[1], valY[1], valY[1],
                                     valY[2], valY[2], valY[2], valY[2],
                                     valY[3], valY[3], valY[3], valY[3],
                                     valY[0], valY[0], valY[0], valY[0],
                                     valY[1], valY[1], valY[1], valY[1],
                                     valY[2], valY[2], valY[2], valY[2],
                                     valY[3], valY[3], valY[3], valY[3],
                                     valY[0], valY[0], valY[0], valY[0],
                                     valY[1], valY[1], valY[1], valY[1],
                                     valY[2], valY[2], valY[2], valY[2],
                                     valY[3], valY[3], valY[3], valY[3],
                                     valY[0], valY[0], valY[0], valY[0],
                                     valY[1], valY[1], valY[1], valY[1],
                                     valY[2], valY[2], valY[2], valY[2],
                                     valY[3], valY[3], valY[3], valY[3]}};

    const DataT valZ[4]{vals0[FZ], pos[FZ], vals2[FZ], vals3[FZ]};
    const Vector<DataT, 64> vecValZ{{valZ[0], valZ[0], valZ[0], valZ[0],
                                     valZ[0], valZ[0], valZ[0], valZ[0],
                                     valZ[0], valZ[0], valZ[0], valZ[0],
                                     valZ[0], valZ[0], valZ[0], valZ[0],
                                     valZ[1], valZ[1], valZ[1], valZ[1],
                                     valZ[1], valZ[1], valZ[1], valZ[1],
                                     valZ[1], valZ[1], valZ[1], valZ[1],
                                     valZ[1], valZ[1], valZ[1], valZ[1],
                                     valZ[2], valZ[2], valZ[2], valZ[2],
                                     valZ[2], valZ[2], valZ[2], valZ[2],
                                     valZ[2], valZ[2], valZ[2], valZ[2],
                                     valZ[2], valZ[2], valZ[2], valZ[2],
                                     valZ[3], valZ[3], valZ[3], valZ[3],
                                     valZ[3], valZ[3], valZ[3], valZ[3],
                                     valZ[3], valZ[3], valZ[3], valZ[3],
                                     valZ[3], valZ[3], valZ[3], valZ[3]}};

    const DataT result = sum(mCoefficients * vecValX * vecValY * vecValZ);
    return result;
  }

  DataT EvalDerivative(const DataT dx, const DataT dy, const DataT dz, const int derx, const int dery, const int derz) const
  {

    DataT ret{};
    DataT cont{};

    for (int i = derx; i < 4; i++) {
      for (int j = dery; j < 4; j++) {
        for (int k = derz; k < 4; k++) {

          const unsigned int index = i + j * 4 + 16 * k;
          cont = mCoefficients[index] * uiPow(dx, i - derx) * uiPow(dy, j - dery) * uiPow(dz, k - derz);
          for (int w = 0; w < derx; w++) {
            cont *= (i - w);
          }
          for (int w = 0; w < dery; w++) {
            cont *= (j - w);
          }
          for (int w = 0; w < derz; w++) {
            cont *= (k - w);
          }
          ret += cont;
        }
      }
    }
    const auto invSpacing = mGridData.getInvSpacing();
    const DataT norm = uiPow(invSpacing[FX], derx) * uiPow(invSpacing[FY], dery) * uiPow(invSpacing[FZ], derz);
    return (ret * norm);
  }

  // for circular padding
  void getDataIndexCircularArray(const int index0, const int dim, int arr[]) const
  {
    const int delta_min1 = getRegulatedDelta(index0, -1, dim, mGridData.getN(dim) - 1);
    const int delta_min2 = getRegulatedDelta(index0, -2, dim, mGridData.getN(dim) - 2);
    const int delta_plus1 = getRegulatedDelta(index0, +1, dim, 1 - mGridData.getN(dim));
    const int delta_plus2 = getRegulatedDelta(index0, +2, dim, 2 - mGridData.getN(dim));

    arr[0] = mGridData.getDeltaDataIndex(delta_min2, dim);
    arr[1] = mGridData.getDeltaDataIndex(delta_min1, dim);
    arr[2] = mGridData.getDeltaDataIndex(delta_plus1, dim);
    arr[3] = mGridData.getDeltaDataIndex(delta_plus2, dim);
  }

  // for non circular padding
  void getDataIndexNonCircularArray(const int index0, const int dim, int arr[]) const
  {
    const int delta_min1 = getRegulatedDelta(index0, -1, dim, 0);
    const int delta_min2 = getRegulatedDelta(index0, -2, dim, delta_min1);
    const int delta_plus1 = getRegulatedDelta(index0, +1, dim, 0);
    const int delta_plus2 = getRegulatedDelta(index0, +2, dim, delta_plus1);

    arr[0] = mGridData.getDeltaDataIndex(delta_min2, dim);
    arr[1] = mGridData.getDeltaDataIndex(delta_min1, dim);
    arr[2] = mGridData.getDeltaDataIndex(delta_plus1, dim);
    arr[3] = mGridData.getDeltaDataIndex(delta_plus2, dim);
  }

  // this helps to get circular and non circular padding indices
  int getRegulatedDelta(const int index0, const int delta, const int dim, const int offs) const
  {
    const int regulatedDelta = mGridData.isIndexInGrid(index0 + delta, dim) ? delta : offs;
    return regulatedDelta;
  }

  // matrix needed to compute the coefficients
  inline static Vc::Memory<VDataT, 64> mA[64]{
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
    {8, -8, -8, 8, -8, 8, 8, -8, 4, 4, -4, -4, -4, -4, 4, 4, 4, -4, 4, -4, -4, 4, -4, 4, 4, -4, -4, 4, 4, -4, -4, 4, 2, 2, 2, 2, -2, -2, -2, -2, 2, 2, -2, -2, 2, 2, -2, -2, 2, -2, 2, -2, 2, -2, 2, -2, 1, 1, 1, 1, 1, 1, 1, 1}};

  inline static Matrix<DataT, 64> mMatrixA{mA};
  mutable Vector<DataT, 64> mCoefficients{};

  static constexpr unsigned int FDim = Grid3D::getDim(); // dimensions of the grid (only 3 supported)
  static constexpr unsigned int FX = Grid3D::getFX();    // index for x coordinate
  static constexpr unsigned int FY = Grid3D::getFY();    // index for y coordinate
  static constexpr unsigned int FZ = Grid3D::getFZ();    // index for z coordinate

  // const unsigned int mNgrid[FDim]{ Grid3D::getN(FX), Grid3D::getN(FY), Grid3D::getN(FZ) }; // number of grid points in direction
  const unsigned int mNgrid[FDim]{Grid3D::getN(FX), Grid3D::getN(FY), Grid3D::getN(FZ)}; // number of grid points in direction

  const Grid3D& mGridData{}; // adress to the data container of the grid

  mutable bool mInitialized = false;             ///< sets the flag if the coefficients are evaluated at least once
  mutable Vector<unsigned int, FDim> mLastInd{}; ///< stores the index for the cell, where the coefficients are already evaluated (only the coefficients for one-the last cell are stored)
  const bool mCircularX{};
  const bool mCircularY{};
  const bool mCircularZ{};
};

#endif
