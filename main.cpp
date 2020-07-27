#include "include/O2TPCSpaceCharge3DCalc.h"
#include <iostream>

int main() {
  /* code */
  const unsigned int nGridR = 129;
  const unsigned int nGridZ = 129;
  const unsigned int nGridPhi = 180;
  O2TPCSpaceCharge3DCalc<float, nGridR, nGridZ, nGridPhi> spaceCharge3DCalc;


  std::cout<<"getGridSizeR(): " << spaceCharge3DCalc.getGridSpacingR()  << std::endl;



  return 0;
}
