#include "src/include/O2TPCSpaceCharge3DCalc.h"
// #include "src/CustomStreamer.cxx"

#include <iostream>

// for debugging
#include "TTree.h"
#include "TFile.h"
#include <chrono>

int main()
{
  using DataT = float;

  const unsigned int nGridR = 17;
  const unsigned int nGridZ = 17;
  const unsigned int nGridPhi = 90;

  //====================== dump regular grid to file testing ====================================
  TFile fOut("grid3D.root", "RECREATE");
  // using dataContType = DataContainer3D<float,nGridR,nGridZ,nGridPhi>;
  DataContainer3D<float, nGridR, nGridZ, nGridPhi> containerOut;

  // set dummy values
  for (size_t i = 0; i < nGridR * nGridZ * nGridPhi; ++i) {
    containerOut[i] = i;
  }

  // dump to disc
  containerOut.writeToFile(fOut, "grid");
  fOut.Close();

  // create a regular grid with the dumped data
  TFile fIn("grid3D.root", "READ");
  RegularGrid3D<float, 17, 17, 90> gridIn{fIn, "grid", 0, 0, 0, 1, 1, 1};

  std::cout << gridIn;

  // modify values in regular grid
  for (size_t iPhi = 0; iPhi < nGridPhi; ++iPhi) {
    for (size_t iR = 0; iR < nGridR; ++iR) {
      for (size_t iZ = 0; iZ < nGridZ; ++iZ) {
        int ind = iZ + iR*nGridZ + nGridZ*nGridR*iPhi;
        gridIn(iZ, iR, iPhi) = ind;
      }
    }
  }

  TFile fOut2("grid3D.root", "RECREATE");
  gridIn.storeValuesToFile(fOut2, "Grid3D");
  fOut2.Close();

  TFile fIn2("grid3D.root", "READ");
  RegularGrid3D<float, 17, 17, 90> gridIn2{fIn2, "grid", 0, 0, 0, 1, 1, 1};

  std::cout << gridIn2;
  //============================================================================================
  return 1;

  O2TPCSpaceCharge3DCalc<DataT, nGridR, nGridZ, nGridPhi> spaceCharge3DCalcAnalytical;
  spaceCharge3DCalcAnalytical.setOmegaTauT1T2(0.32f, 1, 1);
  spaceCharge3DCalcAnalytical.setIntegrationSteps(10);

  O2TPCSpaceCharge3DCalc<DataT, nGridR, nGridZ, nGridPhi> spaceCharge3DCalcNumerical;
  spaceCharge3DCalcNumerical.setOmegaTauT1T2(0.32f, 1, 1);
  spaceCharge3DCalcNumerical.setIntegrationSteps(3);

  AnalyticalFields<DataT> anaFields;
  const int maxIteration = 300;
  const DataT stoppingConvergence = 1e-8;

  auto start = std::chrono::high_resolution_clock::now();
  spaceCharge3DCalcNumerical.fillBoundaryAndChargeDensities(anaFields, maxIteration, stoppingConvergence);
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<float> diff = end - start;
  std::cout << "Step 1: Poisson solver: " << diff.count() << std::endl;

  start = std::chrono::high_resolution_clock::now();
  spaceCharge3DCalcNumerical.calcEField();
  end = std::chrono::high_resolution_clock::now();
  diff = end - start;
  std::cout << "Step 2: Electric Field Calculation: " << diff.count() << std::endl;

  start = std::chrono::high_resolution_clock::now();
  // spaceCharge3DCalcAnalytical.calcLocalDistortionsCorrections(false, anaFields); // local distortion calculation
  // spaceCharge3DCalcAnalytical.calcLocalDistortionsCorrections(true, anaFields);  // local correction calculation
  end = std::chrono::high_resolution_clock::now();
  diff = end - start;
  std::cout << "Step 3: local distortions and corrections analytical: " << diff.count() << std::endl;

  NumericalFields<DataT, nGridR, nGridZ, nGridPhi> numFields(spaceCharge3DCalcNumerical.getGridEr(), spaceCharge3DCalcNumerical.getGridEz(), spaceCharge3DCalcNumerical.getGridEphi());
  start = std::chrono::high_resolution_clock::now();
  spaceCharge3DCalcNumerical.calcLocalDistortionsCorrections(false, numFields); // local distortion calculation
  spaceCharge3DCalcNumerical.calcLocalDistortionsCorrections(true, numFields);  // local correction calculation
  end = std::chrono::high_resolution_clock::now();
  diff = end - start;
  std::cout << "Step 3: local distortions and corrections numerical: " << diff.count() << std::endl;
  return 1;

  // dump to disk
  TFile fDebug(Form("debug_%i_%i_%i.root", nGridR, nGridZ, nGridPhi), "RECREATE");
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

  tree.Branch("ir", &ir);
  tree.Branch("iz", &iz);
  tree.Branch("iphi", &iphi);

  tree.Branch("ldistrana", &localcDistR);
  tree.Branch("ldistzana", &localcDistZ);
  tree.Branch("ldistrphiana", &localcDistRPhi);

  tree.Branch("lcorrrana", &localcCorrR);
  tree.Branch("lcorrzana", &localcCorrZ);
  tree.Branch("lcorrrphiana", &localcCorrRPhi);

  tree.Branch("ldistrnum", &localcDistRNum);
  tree.Branch("ldistznum", &localcDistZNum);
  tree.Branch("ldistrphinum", &localcDistRPhiNum);

  tree.Branch("lcorrrnum", &localcCorrRNum);
  tree.Branch("lcorrznum", &localcCorrZNum);
  tree.Branch("lcorrrphinum", &localcCorrRPhiNum);

  tree.Branch("erana", &eRAna);
  tree.Branch("ezana", &eZAna);
  tree.Branch("ephiana", &ePhiAna);
  tree.Branch("ernum", &eRNum);
  tree.Branch("eznum", &eZNum);
  tree.Branch("ephinum", &ePhiNum);

  // debug tree for tricubic interpolation
  TTree treeInterpolation("inter", "inter");
  DataT r{0};
  DataT phi{0};
  DataT z{0};
  treeInterpolation.Branch("r", &r);
  treeInterpolation.Branch("phi", &phi);
  treeInterpolation.Branch("z", &z);
  treeInterpolation.Branch("erana", &eRAna);
  treeInterpolation.Branch("ezana", &eZAna);
  treeInterpolation.Branch("ephiana", &ePhiAna);
  treeInterpolation.Branch("ernum", &eRNum);
  treeInterpolation.Branch("eznum", &eZNum);
  treeInterpolation.Branch("ephinum", &ePhiNum);

  for (size_t indr = 0; indr < nGridR - 1; ++indr) {
    for (size_t indphi = 0; indphi < nGridPhi - 1; ++indphi) {
      for (size_t indz = 0; indz < nGridZ - 1; ++indz) {
        // interpolate at theses positions
        r = spaceCharge3DCalcAnalytical.getRVertex(indr) + 0.5 * spaceCharge3DCalcAnalytical.getGridSpacingR();
        z = spaceCharge3DCalcAnalytical.getZVertex(indz) + 0.5 * spaceCharge3DCalcAnalytical.getGridSpacingZ();
        phi = spaceCharge3DCalcAnalytical.getPhiVertex(indphi) + 0.5 * spaceCharge3DCalcAnalytical.getGridSpacingPhi();
        eRNum = numFields.evalEr(z, r, phi);
        eZNum = numFields.evalEz(z, r, phi);
        ePhiNum = numFields.evalEphi(z, r, phi);

        eZAna = anaFields.evalEz(z, r, phi);
        eRAna = anaFields.evalEr(z, r, phi);
        ePhiAna = anaFields.evalEphi(z, r, phi);
        treeInterpolation.Fill();
      }
    }
  }

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

        // get numerically calculated Efield
        eZNum = spaceCharge3DCalcNumerical.getEz(indz, indr, indphi);
        eRNum = spaceCharge3DCalcNumerical.getEr(indz, indr, indphi);
        ePhiNum = spaceCharge3DCalcNumerical.getEphi(indz, indr, indphi);

        // get analytical e field
        const DataT radius = spaceCharge3DCalcAnalytical.getRVertex(indr);
        z = spaceCharge3DCalcAnalytical.getZVertex(indz);
        phi = spaceCharge3DCalcAnalytical.getPhiVertex(indphi);
        eZAna = anaFields.evalEz(z, radius, phi);
        eRAna = anaFields.evalEr(z, radius, phi);
        ePhiAna = anaFields.evalEphi(z, radius, phi);

        tree.Fill();
      }
    }
  }
  treeInterpolation.Write();
  tree.Write();
  fDebug.Close();

  return 0;
}
