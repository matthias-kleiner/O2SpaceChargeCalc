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
  const unsigned int nGridZ = 17;
  const unsigned int nGridPhi = 90;

  O2TPCSpaceCharge3DCalc<DataT, nGridR, nGridZ, nGridPhi> spaceCharge3DCalcAnalytical;
  spaceCharge3DCalcAnalytical.setOmegaTauT1T2(0.32f, 1, 1);
  spaceCharge3DCalcAnalytical.setIntegrationSteps(3);

  O2TPCSpaceCharge3DCalc<DataT, nGridR, nGridZ, nGridPhi> spaceCharge3DCalcNumerical;
  spaceCharge3DCalcNumerical.setOmegaTauT1T2(0.32f, 1, 1);
  spaceCharge3DCalcNumerical.setIntegrationSteps(3);

  AnalyticalFields<DataT> anaFields;

  const int maxIteration = 300;
  const DataT stoppingConvergence = 1e-8;

  //==============================================================
  //================ FILL BOUNDARY NUM ===========================
  //==============================================================
  auto start = std::chrono::high_resolution_clock::now();
  // spaceCharge3DCalcNumerical.fillBoundaryAndChargeDensities(anaFields);
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<float> diff = end - start;
  std::cout << "Step 0: Fill Boundary and Charge Densities: " << diff.count() << std::endl;

  //==============================================================
  //================ POISSON SOLVER NUM ==========================
  //==============================================================
  // start = std::chrono::high_resolution_clock::now();
  // spaceCharge3DCalcNumerical.poissonSolver(maxIteration, stoppingConvergence);
  // end = std::chrono::high_resolution_clock::now();
  // std::cout << "Step 1: Poisson solver: " << diff.count() << std::endl;

  //==============================================================
  //================ POISSON SOLVER NUM SET FROM FILE ============
  //==============================================================
  TFile fPot(Form("Potential_%i_%i_%i.root", nGridZ, nGridR, nGridPhi), "READ");
  std::cout << "========= LOADING POTENTIAL: =========" << std::endl;
  // spaceCharge3DCalcNumerical.dumpPotential(fPot);
  spaceCharge3DCalcNumerical.setPotentialFromFile(fPot);
  std::cout << "=========  POTENTIAL LOADED: ========= " << std::endl;

  //==============================================================
  //================ ELECTRIC FIELD NUMERICAL ====================
  //==============================================================
  start = std::chrono::high_resolution_clock::now();
  std::cout << "=========  CALC EFIELD: =========" << std::endl;
  spaceCharge3DCalcNumerical.calcEField(); // calculate e field
  std::cout << "=========  EFIELD CALCULATED: =========" << std::endl;
  end = std::chrono::high_resolution_clock::now();
  diff = end - start;
  std::cout << "Step 2: Electric Field Calculation: " << diff.count() << std::endl;

  // TFile fEField( Form("EField_%i_%i_%i.root", nGridZ, nGridR, nGridPhi), "READ");
  // spaceCharge3DCalcNumerical.setEFieldFromFile(fEField); // calculate e field
  // spaceCharge3DCalcNumerical.dumpElectricFields(fEField);
  // return 1;

  //==============================================================
  //================ LOCAL DISTORTIONS CORRECTIONS AND ===========
  //==============================================================
  start = std::chrono::high_resolution_clock::now();
  spaceCharge3DCalcAnalytical.calcLocalDistortionsCorrections(false, anaFields); // local distortion calculation
  spaceCharge3DCalcAnalytical.calcLocalDistortionsCorrections(true, anaFields);  // local correction calculation
  end = std::chrono::high_resolution_clock::now();
  diff = end - start;
  std::cout << "Step 3: local distortions and corrections analytical: " << diff.count() << std::endl;

  //==============================================================
  //================ LOCAL DISTORTIONS CORRECTIONS NUM ===========
  //==============================================================
  const auto numFields = spaceCharge3DCalcNumerical.getNumericalFieldsInterpolator();
  start = std::chrono::high_resolution_clock::now();
  spaceCharge3DCalcNumerical.calcLocalDistortionsCorrections(false, numFields); // local distortion calculation
  spaceCharge3DCalcNumerical.calcLocalDistortionsCorrections(true, numFields);  // local correction calculation
  end = std::chrono::high_resolution_clock::now();
  diff = end - start;
  std::cout << "Step 3: local distortions and corrections numerical: " << diff.count() << std::endl;

  const auto anaLocalDistInterpolator = spaceCharge3DCalcAnalytical.getLocalDistInterpolator();
  const auto anaLocalCorrInterpolator = spaceCharge3DCalcAnalytical.getLocalCorrInterpolator();

  const auto numLocalDistInterpolator = spaceCharge3DCalcNumerical.getLocalDistInterpolator();
  const auto numLocalCorrInterpolator = spaceCharge3DCalcNumerical.getLocalCorrInterpolator();

  //==============================================================
  //========== GLOBAL DISTORTIONS ANALYTICAL USING LOCAL DIST ====
  //==============================================================
  start = std::chrono::high_resolution_clock::now();
  // spaceCharge3DCalcAnalytical.calcGlobalDistortions(anaLocalDistInterpolator);
  spaceCharge3DCalcAnalytical.calcGlobalDistortions(anaFields);
  end = std::chrono::high_resolution_clock::now();
  TFile fGlobalDistAna(Form("GlobalDist_Ana_%i_%i_%i.root", nGridZ, nGridR, nGridPhi), "RECREATE");
  spaceCharge3DCalcAnalytical.dumpGlobalDistortions(fGlobalDistAna);
  // spaceCharge3DCalcAnalytical.setGlobalDistortionsFromFile(fGlobalDistAna);
  diff = end - start;
  std::cout << "Step 4: global distortions analytical: " << diff.count() << std::endl;

  //==============================================================
  //========== GLOBAL CORRECTIONS ANALYTICAL USING LOCAL CORR ====
  //==============================================================
  start = std::chrono::high_resolution_clock::now();
  // spaceCharge3DCalcAnalytical.calcGlobalDistortions(anaLocalCorrInterpolator);
  spaceCharge3DCalcAnalytical.calcGlobalCorrections(anaFields);
  end = std::chrono::high_resolution_clock::now();
  TFile fGlobalCorrAna(Form("GlobalCorr_Ana_%i_%i_%i.root", nGridZ, nGridR, nGridPhi), "RECREATE");
  spaceCharge3DCalcAnalytical.dumpGlobalCorrections(fGlobalCorrAna);
  // spaceCharge3DCalcAnalytical.setGlobalCorrectionsFromFile(fGlobalCorrAna);
  diff = end - start;
  std::cout << "Step 4: global corrections analytical: " << diff.count() << std::endl;

  //==============================================================
  //========== GLOBAL DISTORTIONS NUMERICAL USING LOCAL DIST ====
  //==============================================================
  start = std::chrono::high_resolution_clock::now();
  spaceCharge3DCalcNumerical.calcGlobalDistortions(numFields);
  // spaceCharge3DCalcNumerical.calcGlobalDistortions(numLocalDistInterpolator);
  end = std::chrono::high_resolution_clock::now();
  TFile fGlobalDistNum(Form("GlobalDist_Num_%i_%i_%i.root", nGridZ, nGridR, nGridPhi), "RECREATE");
  spaceCharge3DCalcNumerical.dumpGlobalDistortions(fGlobalDistNum);
  // spaceCharge3DCalcNumerical.setGlobalDistortionsFromFile(fGlobalDistNum);
  diff = end - start;
  std::cout << "Step 4: global distortions numerical: " << diff.count() << std::endl;

  //==============================================================
  //========== GLOBAL CORRECTIONS NUMERICAL ======================
  //==============================================================
  start = std::chrono::high_resolution_clock::now();
  spaceCharge3DCalcNumerical.calcGlobalCorrections(numFields);
  // spaceCharge3DCalcNumerical.calcGlobalDistortions(numLocalDistInterpolator);
  end = std::chrono::high_resolution_clock::now();
  TFile fGlobalCorrNum(Form("GlobalCorr_Num_%i_%i_%i.root", nGridZ, nGridR, nGridPhi), "RECREATE");
  spaceCharge3DCalcNumerical.dumpGlobalCorrections(fGlobalCorrNum);
  // spaceCharge3DCalcNumerical.setGlobalCorrectionsFromFile(fGlobalCorrNum);
  diff = end - start;
  std::cout << "Step 4: global corrections numerical: " << diff.count() << std::endl;

  // dump to disk
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

  tree.Branch("ir", &ir);
  tree.Branch("iz", &iz);
  tree.Branch("iphi", &iphi);

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
