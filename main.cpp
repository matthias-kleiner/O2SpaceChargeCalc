#include "src/include/O2TPCSpaceCharge3DCalc.h"
// #include "src/CustomStreamer.cxx"

#include <iostream>

// for debugging
#include "TTree.h"
#include "TFile.h"
#include <chrono>

int main(int argc, char const* argv[])
{
  std::stringstream strValue;
  strValue << argv[1];
  int nThreads;
  strValue >> nThreads;
  std::cout << "nThreads: " << nThreads << std::endl;
  omp_set_num_threads(nThreads);

  using DataT = float;
  const unsigned int nGridR = 17;
  const unsigned int nGridZ = nGridR;
  const unsigned int nGridPhi = 90;

  O2TPCSpaceCharge3DCalc<DataT, nGridR, nGridZ, nGridPhi> spaceCharge3DCalcAnalytical;
  spaceCharge3DCalcAnalytical.setOmegaTauT1T2(0.32f, 1, 1);
  spaceCharge3DCalcAnalytical.setNStep(1);
  int integrationStrategy = O2TPCSpaceCharge3DCalc<>::IntegrationStrategy::Simpson;
  spaceCharge3DCalcAnalytical.setNumericalIntegrationStrategy(integrationStrategy);

  O2TPCSpaceCharge3DCalc<DataT, nGridR, nGridZ, nGridPhi> spaceCharge3DCalcNumerical;
  spaceCharge3DCalcNumerical.setOmegaTauT1T2(0.32f, 1, 1);
  spaceCharge3DCalcNumerical.setNStep(1);
  spaceCharge3DCalcNumerical.setNumericalIntegrationStrategy(integrationStrategy);

  AnalyticalFields<DataT> anaFields;
  TFile fAna( Form("ANA_%i_%i_%i.root", nGridZ, nGridR, nGridPhi), "RECREATE");
  spaceCharge3DCalcAnalytical.performFullRun(anaFields,0,fAna);
  // spaceCharge3DCalcAnalytical.setFromFile(0, fAna);

  TFile fNum( Form("NUM_%i_%i_%i.root", nGridZ, nGridR, nGridPhi), "RECREATE");
  spaceCharge3DCalcNumerical.performFullRun(anaFields,1,fNum);
  // spaceCharge3DCalcNumerical.setFromFile(0, fNum);

  const auto anaGlobalDistInterpolator = spaceCharge3DCalcAnalytical.getGlobalDistInterpolator();
  const auto anaGlobalCorrInterpolator = spaceCharge3DCalcAnalytical.getGlobalCorrInterpolator();

  const auto numGlobalDistInterpolator = spaceCharge3DCalcNumerical.getGlobalDistInterpolator();
  const auto numGlobalCorrInterpolator = spaceCharge3DCalcNumerical.getGlobalCorrInterpolator();

  // //alternative approach of calculating the global distortions/correction
  // //==============================================================
  // //========== GLOBAL DISTORTIONS ANA ITERATIVE ==================
  // //==============================================================
  // start = std::chrono::high_resolution_clock::now();
  // // spaceCharge3DCalcAnalytical.calcGlobalDistWithGlobalCorrIterative(anaGlobalCorrInterpolator);
  // end = std::chrono::high_resolution_clock::now();
  // diff = end - start;
  // std::cout << "iterative algorithm analytical: " << diff.count() << std::endl;
  // // return 1;
  // //alternative approach of calculating the global distortions/correction
  // //==============================================================
  // //========== GLOBAL DISTORTIONS NUM ITERATIVE ==================
  // //==============================================================
  // start = std::chrono::high_resolution_clock::now();
  // spaceCharge3DCalcNumerical.calcGlobalDistWithGlobalCorrIterative(numGlobalCorrInterpolator);
  // end = std::chrono::high_resolution_clock::now();
  // diff = end - start;
  // std::cout << "iterative algorithm numerical: " << diff.count() << std::endl;
  // // return 1;


  TFile fDebug(Form("debug_%i_%i_%i.root", nGridR, nGridZ, nGridPhi), "RECREATE");
  TTree tree("tree", "tree");
  int ir{0};
  int iz{0};
  int iphi{0};
  DataT localcDistR{0};
  DataT localcDistZ{0};
  DataT localcDistRPhi{0};

  DataT globalDistR{0};
  DataT globalDistZ{0};
  DataT globalDistRPhi{0};

  DataT globalDistRNum{0};
  DataT globalDistZNum{0};
  DataT globalDistRPhiNum{0};

  DataT globalCorrR{0};
  DataT globalCorrZ{0};
  DataT globalCorrRPhi{0};

  DataT globalCorrRNum{0};
  DataT globalCorrZNum{0};
  DataT globalCorrRPhiNum{0};

  DataT localcCorrR{0};
  DataT localcCorrZ{0};
  DataT localcCorrRPhi{0};

  DataT localcDistRNum{0};
  DataT localcDistZNum{0};
  DataT localcDistRPhiNum{0};

  DataT localcCorrRNum{0};
  DataT localcCorrZNum{0};
  DataT localcCorrRPhiNum{0};

  DataT eRAna{0};
  DataT eZAna{0};
  DataT ePhiAna{0};
  DataT eRNum{0};
  DataT eZNum{0};
  DataT ePhiNum{0};

  DataT r{0};
  DataT phi{0};
  DataT z{0};

  tree.Branch("ir", &ir);
  tree.Branch("iz", &iz);
  tree.Branch("iphi", &iphi);
  tree.Branch("r", &r);
  tree.Branch("z", &z);
  tree.Branch("phi", &phi);
  DataT rDistortedPoint = 0;
  DataT zDistortedPoint = 0;
  DataT phiDistortedPoint = 0;

  int correctionInVolume = 0;

  tree.Branch("inVolume", &correctionInVolume);
  tree.Branch("rDistortedPoint", &rDistortedPoint);
  tree.Branch("zDistortedPoint", &zDistortedPoint);
  tree.Branch("phiDistortedPoint", &phiDistortedPoint);

  tree.Branch("ldistrana", &localcDistR);
  tree.Branch("ldistzana", &localcDistZ);
  tree.Branch("ldistrphiana", &localcDistRPhi);

  tree.Branch("distrana", &globalDistR);
  tree.Branch("distzana", &globalDistZ);
  tree.Branch("distrphiana", &globalDistRPhi);

  tree.Branch("distrnum", &globalDistRNum);
  tree.Branch("distznum", &globalDistZNum);
  tree.Branch("distrphinum", &globalDistRPhiNum);

  tree.Branch("corrrana", &globalCorrR);
  tree.Branch("corrzana", &globalCorrZ);
  tree.Branch("corrrphiana", &globalCorrRPhi);

  tree.Branch("corrrnum", &globalCorrRNum);
  tree.Branch("corrznum", &globalCorrZNum);
  tree.Branch("corrrphinum", &globalCorrRPhiNum);

  tree.Branch("lcorrrana", &localcCorrR);
  tree.Branch("lcorrzana", &localcCorrZ);
  tree.Branch("lcorrrphiana", &localcCorrRPhi);

  tree.Branch("ldistrnum", &localcDistRNum);
  tree.Branch("ldistznum", &localcDistZNum);
  tree.Branch("ldistrphinum", &localcDistRPhiNum);

  tree.Branch("lcorrrnum", &localcCorrRNum);
  tree.Branch("lcorrznum", &localcCorrZNum);
  tree.Branch("lcorrrphinum", &localcCorrRPhiNum);

  DataT corrDistPoint[3]{};
  tree.Branch("corrRDistortedPoint", &corrDistPoint[0]);
  tree.Branch("corrZDistortedPoint", &corrDistPoint[1]);
  tree.Branch("corrRPhiDistortedPoint", &corrDistPoint[2]);


  tree.Branch("erana", &eRAna);
  tree.Branch("ezana", &eZAna);
  tree.Branch("ephiana", &ePhiAna);
  tree.Branch("ernum", &eRNum);
  tree.Branch("eznum", &eZNum);
  tree.Branch("ephinum", &ePhiNum);

  // debug tree for tricubic interpolation
  TTree treeInterpolation("inter", "inter");
  treeInterpolation.Branch("r", &r);
  treeInterpolation.Branch("phi", &phi);
  treeInterpolation.Branch("z", &z);
  treeInterpolation.Branch("erana", &eRAna);
  treeInterpolation.Branch("ezana", &eZAna);
  treeInterpolation.Branch("ephiana", &ePhiAna);
  treeInterpolation.Branch("ernum", &eRNum);
  treeInterpolation.Branch("eznum", &eZNum);
  treeInterpolation.Branch("ephinum", &ePhiNum);

  // for (size_t indr = 0; indr < nGridR - 1; ++indr) {
  //   for (size_t indphi = 0; indphi < nGridPhi - 1; ++indphi) {
  //     for (size_t indz = 0; indz < nGridZ - 1; ++indz) {
  //       // interpolate at theses positions
  //       r = spaceCharge3DCalcAnalytical.getRVertex(indr) + 0.5 * spaceCharge3DCalcAnalytical.getGridSpacingR();
  //       z = spaceCharge3DCalcAnalytical.getZVertex(indz) + 0.5 * spaceCharge3DCalcAnalytical.getGridSpacingZ();
  //       phi = spaceCharge3DCalcAnalytical.getPhiVertex(indphi) + 0.5 * spaceCharge3DCalcAnalytical.getGridSpacingPhi();
  //       eRNum = numEFields.evalEr(z, r, phi);
  //       eZNum = numEFields.evalEz(z, r, phi);
  //       ePhiNum = numEFields.evalEphi(z, r, phi);
  //
  //       eZAna = anaFields.evalEz(z, r, phi);
  //       eRAna = anaFields.evalEr(z, r, phi);
  //       ePhiAna = anaFields.evalEphi(z, r, phi);
  //       treeInterpolation.Fill();
  //     }
  //   }
  // }

  // loop over regular grid
  for (size_t indr = 0; indr < nGridR; ++indr) {
    for (size_t indphi = 0; indphi < nGridPhi; ++indphi) {
      for (size_t indz = 0; indz < nGridZ; ++indz) {
        ir = indr;
        iz = indz;
        iphi = indphi;

        // get local distortions and local corrections calclulated by the the analytical Efield
        localcDistR = spaceCharge3DCalcAnalytical.getLocalDistR(indz, indr, indphi);
        localcDistZ = spaceCharge3DCalcAnalytical.getLocalDistZ(indz, indr, indphi);
        localcDistRPhi = spaceCharge3DCalcAnalytical.getLocalDistRPhi(indz, indr, indphi);

        localcCorrR = spaceCharge3DCalcAnalytical.getLocalCorrR(indz, indr, indphi);
        localcCorrZ = spaceCharge3DCalcAnalytical.getLocalCorrZ(indz, indr, indphi);
        localcCorrRPhi = spaceCharge3DCalcAnalytical.getLocalCorrRPhi(indz, indr, indphi);

        // get local distortions and local corrections calclulated by the the analytical Efield
        localcDistRNum = spaceCharge3DCalcNumerical.getLocalDistR(indz, indr, indphi);
        localcDistZNum = spaceCharge3DCalcNumerical.getLocalDistZ(indz, indr, indphi);
        localcDistRPhiNum = spaceCharge3DCalcNumerical.getLocalDistRPhi(indz, indr, indphi);

        localcCorrRNum = spaceCharge3DCalcNumerical.getLocalCorrR(indz, indr, indphi);
        localcCorrZNum = spaceCharge3DCalcNumerical.getLocalCorrZ(indz, indr, indphi);
        localcCorrRPhiNum = spaceCharge3DCalcNumerical.getLocalCorrRPhi(indz, indr, indphi);

        // get global distortions
        globalDistR = spaceCharge3DCalcAnalytical.getGlobalDistR(indz, indr, indphi);
        globalDistZ = spaceCharge3DCalcAnalytical.getGlobalDistZ(indz, indr, indphi);
        globalDistRPhi = spaceCharge3DCalcAnalytical.getGlobalDistRPhi(indz, indr, indphi);
        globalDistRNum = spaceCharge3DCalcNumerical.getGlobalDistR(indz, indr, indphi);
        globalDistZNum = spaceCharge3DCalcNumerical.getGlobalDistZ(indz, indr, indphi);
        globalDistRPhiNum = spaceCharge3DCalcNumerical.getGlobalDistRPhi(indz, indr, indphi);

        // get global corrections
        globalCorrR = spaceCharge3DCalcAnalytical.getGlobalCorrR(indz, indr, indphi);
        globalCorrZ = spaceCharge3DCalcAnalytical.getGlobalCorrZ(indz, indr, indphi);
        globalCorrRPhi = spaceCharge3DCalcAnalytical.getGlobalCorrRPhi(indz, indr, indphi);
        globalCorrRNum = spaceCharge3DCalcNumerical.getGlobalCorrR(indz, indr, indphi);
        globalCorrZNum = spaceCharge3DCalcNumerical.getGlobalCorrZ(indz, indr, indphi);
        globalCorrRPhiNum = spaceCharge3DCalcNumerical.getGlobalCorrRPhi(indz, indr, indphi);

        // get numerically calculated Efield
        eZNum = spaceCharge3DCalcNumerical.getEz(indz, indr, indphi);
        eRNum = spaceCharge3DCalcNumerical.getEr(indz, indr, indphi);
        ePhiNum = spaceCharge3DCalcNumerical.getEphi(indz, indr, indphi);

        // get analytical e field
        const DataT radius = spaceCharge3DCalcAnalytical.getRVertex(indr);
        r = radius;
        z = spaceCharge3DCalcAnalytical.getZVertex(indz);
        phi = spaceCharge3DCalcAnalytical.getPhiVertex(indphi);
        eZAna = anaFields.evalEz(z, radius, phi);
        eRAna = anaFields.evalEr(z, radius, phi);
        ePhiAna = anaFields.evalEphi(z, radius, phi);

        // distort point and then correct it
        const bool safe = true;
        const DataT posDistorted[3] = { z + globalDistZ, radius + globalDistR, phi + globalDistRPhi / radius };
        const DataT corrRDistPointTmp = anaGlobalCorrInterpolator.evaldR(posDistorted[0], posDistorted[1], posDistorted[2], safe);
        corrDistPoint[0] = corrRDistPointTmp;
        corrDistPoint[1] = anaGlobalCorrInterpolator.evaldZ(posDistorted[0], posDistorted[1], posDistorted[2], safe);
        const DataT corrRPhiDistPointTmp = anaGlobalCorrInterpolator.evaldRPhi(posDistorted[0], posDistorted[1], posDistorted[2], safe);
        corrDistPoint[2] = corrRPhiDistPointTmp  * (1 + corrRDistPointTmp/(posDistorted[1]));
        zDistortedPoint = posDistorted[0];
        rDistortedPoint = posDistorted[1];
        phiDistortedPoint = posDistorted[2];

        const DataT posCorrected[3] = { z + globalCorrZ, radius + globalCorrR, phi + globalCorrRPhi / radius };
        correctionInVolume = 1;
        if(posCorrected[0]<spaceCharge3DCalcNumerical.getZMin() || posCorrected[0]>spaceCharge3DCalcNumerical.getZMax() || posCorrected[1]<spaceCharge3DCalcNumerical.getRMin() || posCorrected[1]>spaceCharge3DCalcNumerical.getRMax()){
          correctionInVolume = 0;
        }

        tree.Fill();
      }
    }
  }
  treeInterpolation.Write();
  tree.Write();
  fDebug.Close();

  return 0;
}
