#include "src/include/O2TPCSpaceCharge3DCalc.h"
#include <iostream>

// for debugging
#include "TTree.h"
#include "TFile.h"

int main()
{
  using DataT = float;

  const unsigned int nGridR = 17;
  const unsigned int nGridZ = 17;
  const unsigned int nGridPhi = 90;
  O2TPCSpaceCharge3DCalc<DataT, nGridR, nGridZ, nGridPhi> spaceCharge3DCalc;
  spaceCharge3DCalc.setOmegaTauT1T2(0.32f, 1, 1);
  spaceCharge3DCalc.setIntegrationSteps(3);

  Formulas<DataT> formulas;

  // std::cout<<"eval: "<<formulas.evalDensity(100, 1, 0)<<std::endl;
  // return 1;

  const int maxIteration = 300;
  const DataT stoppingConvergence = 1e-8;
  spaceCharge3DCalc.fillBoundaryAndChargeDensities(formulas, maxIteration, stoppingConvergence);
  spaceCharge3DCalc.calcEField();
  spaceCharge3DCalc.calcLocalDistortionsCorrections(false, formulas); // local distortion calculation
  spaceCharge3DCalc.calcLocalDistortionsCorrections(true, formulas); // local correction calculation

  // dump to disk
  TFile fDebug("debug.root", "RECREATE");
  TTree tree("tree", "tree");
  int ir{0};
  int iz{0};
  int iphi{0};
  DataT localcDistR{0};
  DataT localcDistZ{0};
  DataT localcDistRPhi{0};

  DataT localcCorrR{0};
  DataT localcCorrZ{0};
  DataT localcCorrRPhi{0};

  tree.Branch("ir", &ir);
  tree.Branch("iz", &iz);
  tree.Branch("iphi", &iphi);
  tree.Branch("ldistr", &localcDistR);
  tree.Branch("ldistz", &localcDistZ);
  tree.Branch("ldistrphi", &localcDistRPhi);

  tree.Branch("lcorrr", &localcCorrR);
  tree.Branch("lcorrz", &localcCorrZ);
  tree.Branch("lcorrrphi", &localcCorrRPhi);

  // loop over regular grid
  for (size_t indr = 0; indr < nGridR; ++indr) {
    for (size_t indphi = 0; indphi < nGridPhi; ++indphi) {
      for (size_t indz = 0; indz < nGridZ; ++indz) {
        ir = indr;
        iz = indz;
        iphi = indphi;

        localcDistR = spaceCharge3DCalc.getLocalDistR(indz,indr,indphi);
        localcDistZ = spaceCharge3DCalc.getLocalDistZ(indz,indr,indphi);
        localcDistRPhi = spaceCharge3DCalc.getLocalDistRPhi(indz,indr,indphi);

        localcCorrR = spaceCharge3DCalc.getLocalCorrR(indz,indr,indphi);
        localcCorrZ = spaceCharge3DCalc.getLocalCorrZ(indz,indr,indphi);
        localcCorrRPhi = spaceCharge3DCalc.getLocalCorrRPhi(indz,indr,indphi);

        tree.Fill();
      }
    }
  }
  tree.Write();
  fDebug.Close();

  return 0;
}
