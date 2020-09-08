#include "src/include/O2TPCSpaceCharge3DCalc.h"
#include "src/macros/CalculateDistortionsCorrections.C"

#include <iostream>

// for debugging
#include "TTree.h"
#include "TFile.h"
#include "src/include/O2/Defs.h"
#include <chrono>

#include "/Users/matthias/Desktop/testCube/IntpTricubic3.h"

float tricubicSpeedTest();

float FUNC(float x, float y, float z)
{
  return x * y / 10. + std::cos(z);
  // return y * x + 0.5*std::cos(z);
}

void debug();

int main(int argc, char const* argv[])
{
  // tricubicSpeedTest();
  // return 1;

  std::stringstream strValue;
  strValue << argv[1];
  int nThreads;
  strValue >> nThreads;
  std::cout << "nThreads: " << nThreads << std::endl;
  omp_set_num_threads(nThreads);



  using namespace o2::tpc;

  o2::tpc::MGParameters::maxLoop = 7; // this speeds up signifcantly the calculation of the poisson solver

  using DataT = double;
  const unsigned int nGridR = 129;
  const unsigned int nGridZ = nGridR;
  const unsigned int nGridPhi = 180;
  o2::tpc::Side side = o2::tpc::Side::C;
  AnalyticalFields<DataT> anaFields(side);
  int integrationStrategy = O2TPCSpaceCharge3DCalc<>::IntegrationStrategy::SimpsonIterative;
  // int integrationStrategy = O2TPCSpaceCharge3DCalc<>::IntegrationStrategy::Simpson;

  const int globalEFieldTypeAna = 0; // global distortions/corrections:                          0: using electric field, 1: using local dis/corr interpolator
  const int globalDistTypeAna = 0;   // global distortions:                                      0: standard method,      1: interpolation of global corrections
  const int eFieldTypeAna = 0;       // electrc field:                                           0: analytical formula,   1: tricubic interpolator
  const int usePoissonSolverAna = 0; // use poisson solver or analytical formula for potential:  0: analytical formula,   1: poisson solver

  const int globalEFieldTypeNum = 1; // global distortions/corrections:                          0: using electric field, 1: using local dis/corr interpolator
  const int globalDistTypeNum = 0;   // global distortions:                                      0: standard method,      1: interpolation of global corrections
  const int eFieldTypeNum = 1;       // electrc field:                                           0: analytical formula,   1: tricubic interpolator
  const int usePoissonSolverNum = 1; // use poisson solver or analytical formula for potential:  0: analytical formula,   1: poisson solver

  O2TPCSpaceCharge3DCalc<DataT, nGridR, nGridZ, nGridPhi> spaceCharge3DCalcAnalytical;
  spaceCharge3DCalcAnalytical.setOmegaTauT1T2(0.32f, 1, 1);
  spaceCharge3DCalcAnalytical.setNStep(10);
  spaceCharge3DCalcAnalytical.setNumericalIntegrationStrategy(integrationStrategy);

  O2TPCSpaceCharge3DCalc<DataT, nGridR, nGridZ, nGridPhi> spaceCharge3DCalcNumerical;
  spaceCharge3DCalcNumerical.setOmegaTauT1T2(0.32f, 1, 1);
  spaceCharge3DCalcNumerical.setNStep(10);
  spaceCharge3DCalcNumerical.setNumericalIntegrationStrategy(integrationStrategy);

  TFile fAna(Form("ANA_%i_%i_%i.root", nGridZ, nGridR, nGridPhi), "READ");
  spaceCharge3DCalcAnalytical.setFromFile(fAna, side);
  // calculateDistortionsCorrectionsAnalytical(spaceCharge3DCalcAnalytical, anaFields, fAna, globalEFieldTypeAna, eFieldTypeAna, globalDistTypeAna, usePoissonSolverAna, nThreads);
  fAna.Close();

  TFile fNum(Form("NUM_%i_%i_%i.root", nGridZ, nGridR, nGridPhi), "RECREATE");
  // calculateDistortionsCorrections(spaceCharge3DCalcNumerical, fNum, side, "/Users/matthias/alice/o2_macros/SpaceCharges/ATO-498/InputSCDensityHistograms_8000events.root", "inputSCDensity3D_8000_avg", globalEFieldTypeNum, globalDistTypeNum, nThreads);
  // spaceCharge3DCalcNumerical.setFromFile(fNum, side);
  calculateDistortionsCorrectionsAnalytical(spaceCharge3DCalcNumerical, anaFields, fNum, globalEFieldTypeNum, eFieldTypeNum, globalDistTypeNum, usePoissonSolverNum, nThreads);
  fNum.Close();

  // const auto numEFields = spaceCharge3DCalcNumerical.getElectricFieldsInterpolator(side);
  // int type = o2::tpc::O2TPCSpaceCharge3DCalc<>::Type::Corrections;
  // spaceCharge3DCalcNumerical.calcLocalDistortionsCorrections(type, numEFields);
  // spaceCharge3DCalcAnalytical.calcGlobalCorrections(anaFields);

  const auto anaGlobalDistInterpolator = spaceCharge3DCalcAnalytical.getGlobalDistInterpolator(side);
  const auto anaGlobalCorrInterpolator = spaceCharge3DCalcAnalytical.getGlobalCorrInterpolator(side);

  const auto numGlobalDistInterpolator = spaceCharge3DCalcNumerical.getGlobalDistInterpolator(side);
  const auto numGlobalCorrInterpolator = spaceCharge3DCalcNumerical.getGlobalCorrInterpolator(side);

  const auto anaLocalDistInterpolator = spaceCharge3DCalcAnalytical.getLocalDistInterpolator(side);
  const auto anaLocalCorrInterpolator = spaceCharge3DCalcAnalytical.getLocalCorrInterpolator(side);

  const auto numLocalDistInterpolator = spaceCharge3DCalcNumerical.getLocalDistInterpolator(side);
  const auto numLocalCorrInterpolator = spaceCharge3DCalcNumerical.getLocalCorrInterpolator(side);

  //
  // // //alternative approach of calculating the global distortions/correction
  // // //==============================================================
  // // //========== GLOBAL DISTORTIONS ANA ITERATIVE ==================
  // // //==============================================================
  // auto start = std::chrono::high_resolution_clock::now();
  // // spaceCharge3DCalcAnalytical.calcGlobalDistWithGlobalCorrIterative(anaGlobalCorrInterpolator);
  // auto end = std::chrono::high_resolution_clock::now();
  // std::chrono::duration<float> diff = end - start;
  // std::cout << "iterative algorithm analytical CGAL: " << diff.count() << std::endl;
  // // // return 1;
  // // //alternative approach of calculating the global distortions/correction
  // // //==============================================================
  // // //========== GLOBAL DISTORTIONS NUM ITERATIVE ==================
  // // //==============================================================
  // start = std::chrono::high_resolution_clock::now();
  // // spaceCharge3DCalcNumerical.calcGlobalDistWithGlobalCorrIterative(numGlobalCorrInterpolator);
  // // spaceCharge3DCalcAnalytical.calcGlobalDistWithGlobalCorrIterativeCGAL(numGlobalCorrInterpolator);
  // end = std::chrono::high_resolution_clock::now();
  // diff = end - start;
  // std::cout << "iterative algorithm numerical: " << diff.count() << std::endl;
  // // // return 1;

  TFile fDebug(Form("debug_%i_%i_%i.root", nGridR, nGridZ, nGridPhi), "RECREATE");
  TTree tree("tree", "tree");
  int ir{0};
  int iz{0};
  int iphi{0};
  DataT localDistR{0};
  DataT localDistZ{0};
  DataT localDistRPhi{0};

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

  DataT localCorrR{0};
  DataT localCorrZ{0};
  DataT localCorrRPhi{0};

  DataT localDistRNum{0};
  DataT localDistZNum{0};
  DataT localDistRPhiNum{0};

  DataT localCorrRNum{0};
  DataT localCorrZNum{0};
  DataT localCorrRPhiNum{0};

  DataT densNum{0};
  DataT eRAna{0};
  DataT eZAna{0};
  DataT ePhiAna{0};
  DataT eRNum{0};
  DataT eZNum{0};
  DataT ePhiNum{0};
  DataT potNum{0};
  DataT potAna{0};

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

  tree.Branch("ldistrana", &localDistR);
  tree.Branch("ldistzana", &localDistZ);
  tree.Branch("ldistrphiana", &localDistRPhi);

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

  tree.Branch("lcorrrana", &localCorrR);
  tree.Branch("lcorrzana", &localCorrZ);
  tree.Branch("lcorrrphiana", &localCorrRPhi);

  tree.Branch("ldistrnum", &localDistRNum);
  tree.Branch("ldistznum", &localDistZNum);
  tree.Branch("ldistrphinum", &localDistRPhiNum);

  tree.Branch("lcorrrnum", &localCorrRNum);
  tree.Branch("lcorrznum", &localCorrZNum);
  tree.Branch("lcorrrphinum", &localCorrRPhiNum);

  DataT corrDistPoint[3]{};
  tree.Branch("corrRDistortedPointnum", &corrDistPoint[0]);
  tree.Branch("corrZDistortedPointnum", &corrDistPoint[1]);
  tree.Branch("corrRPhiDistortedPointnum", &corrDistPoint[2]);

  DataT corrDistPointAna[3]{};
  tree.Branch("corrRDistortedPointana", &corrDistPointAna[0]);
  tree.Branch("corrZDistortedPointana", &corrDistPointAna[1]);
  tree.Branch("corrRPhiDistortedPointana", &corrDistPointAna[2]);

  DataT lcorrDistPoint[3]{};
  tree.Branch("lcorrRDistortedPointnum", &lcorrDistPoint[0]);
  tree.Branch("lcorrZDistortedPointnum", &lcorrDistPoint[1]);
  tree.Branch("lcorrRPhiDistortedPointnum", &lcorrDistPoint[2]);

  DataT lcorrDistPointAna[3]{};
  tree.Branch("lcorrRDistortedPointana", &lcorrDistPointAna[0]);
  tree.Branch("lcorrZDistortedPointana", &lcorrDistPointAna[1]);
  tree.Branch("lcorrRPhiDistortedPointana", &lcorrDistPointAna[2]);

  tree.Branch("densnum", &densNum);
  tree.Branch("erana", &eRAna);
  tree.Branch("ezana", &eZAna);
  tree.Branch("ephiana", &ePhiAna);
  tree.Branch("ernum", &eRNum);
  tree.Branch("eznum", &eZNum);
  tree.Branch("ephinum", &ePhiNum);
  tree.Branch("potana", &potAna);
  tree.Branch("potnum", &potNum);

  // debug tree for tricubic interpolation
  // TTree treeInterpolation("inter", "inter");
  // treeInterpolation.Branch("r", &r);
  // treeInterpolation.Branch("phi", &phi);
  // treeInterpolation.Branch("z", &z);
  // treeInterpolation.Branch("erana", &eRAna);
  // treeInterpolation.Branch("ezana", &eZAna);
  // treeInterpolation.Branch("ephiana", &ePhiAna);
  // treeInterpolation.Branch("ernum", &eRNum);
  // treeInterpolation.Branch("eznum", &eZNum);
  // treeInterpolation.Branch("ephinum", &ePhiNum);

  spaceCharge3DCalcAnalytical.setEField(anaFields);
  // const auto interpol = spaceCharge3DCalcAnalytical.getElectricFieldsInterpolator();
  // float gridSpacingZ = spaceCharge3DCalcAnalytical.getGridSpacingZ();
  // float gridSpacingPhi = spaceCharge3DCalcAnalytical.getGridSpacingPhi();

  // float zz = 0.2 + 2*gridSpacingZ + 0.001;
  // float rr = 0;
  // float phii = 0 -2*gridSpacingPhi + 0.001;
  // std::cout<< interpol.evalEr(zz, rr, phii) << std::endl;
  // std::cout<< anaFields.evalEr(zz, rr, phii) << std::endl;
  // eZNum = interpol.evalEz(zz, rr, phii);
  // ePhiNum = interpol.evalEphi(zz, rr, phii);

  // eZAna = anaFields.evalEz(zz, rr, phii);
  // eRAna = anaFields.evalEr(zz, rr, phii);
  // ePhiAna = anaFields.evalEphi(zz, rr, phii);
  // return 1;
  //
  //
  // for (size_t indr = 0; indr < nGridR - 1; ++indr) {
  //   for (size_t indphi = 0; indphi < nGridPhi - 1; ++indphi) {
  //     for (size_t indz = 0; indz < nGridZ - 1; ++indz) {
  //       // interpolate at theses positions
  //       r = spaceCharge3DCalcAnalytical.getRVertex(indr) - 1 * spaceCharge3DCalcAnalytical.getGridSpacingR();
  //       z = spaceCharge3DCalcAnalytical.getZVertex(indz) - 1 * spaceCharge3DCalcAnalytical.getGridSpacingZ();
  //       phi = spaceCharge3DCalcAnalytical.getPhiVertex(indphi) - 10 * spaceCharge3DCalcAnalytical.getGridSpacingPhi();
  //       eRNum = interpol.evalEr(z, r, phi);
  //       eZNum = interpol.evalEz(z, r, phi);
  //       ePhiNum = interpol.evalEphi(z, r, phi);
  //
  //       eZAna = anaFields.evalEz(z, r, phi);
  //       eRAna = anaFields.evalEr(z, r, phi);
  //       ePhiAna = anaFields.evalEphi(z, r, phi);
  //       treeInterpolation.Fill();
  //     }
  //   }
  // }

  // loop over regular grid
  for (size_t indphi = 0; indphi < nGridPhi; ++indphi) {
    std::cout<<"indphi: "<<indphi<<std::endl;
    for (size_t indr = 0; indr < nGridR; ++indr) {
      for (size_t indz = 0; indz < nGridZ; ++indz) {
        ir = indr;
        iz = indz;
        iphi = indphi;

        // get local distortions and local corrections calclulated by the the analytical Efield
        localDistR = spaceCharge3DCalcAnalytical.getLocalDistR(indz, indr, indphi, side);
        localDistZ = spaceCharge3DCalcAnalytical.getLocalDistZ(indz, indr, indphi, side);
        localDistRPhi = spaceCharge3DCalcAnalytical.getLocalDistRPhi(indz, indr, indphi, side);

        localCorrR = spaceCharge3DCalcAnalytical.getLocalCorrR(indz, indr, indphi, side);
        localCorrZ = spaceCharge3DCalcAnalytical.getLocalCorrZ(indz, indr, indphi, side);
        localCorrRPhi = spaceCharge3DCalcAnalytical.getLocalCorrRPhi(indz, indr, indphi, side);

        // get local distortions and local corrections calclulated by the the analytical Efield
        localDistRNum = spaceCharge3DCalcNumerical.getLocalDistR(indz, indr, indphi, side);
        localDistZNum = spaceCharge3DCalcNumerical.getLocalDistZ(indz, indr, indphi, side);
        localDistRPhiNum = spaceCharge3DCalcNumerical.getLocalDistRPhi(indz, indr, indphi, side);

        localCorrRNum = spaceCharge3DCalcNumerical.getLocalCorrR(indz, indr, indphi, side);
        localCorrZNum = spaceCharge3DCalcNumerical.getLocalCorrZ(indz, indr, indphi, side);
        localCorrRPhiNum = spaceCharge3DCalcNumerical.getLocalCorrRPhi(indz, indr, indphi, side);

        // get global distortions
        globalDistR = spaceCharge3DCalcAnalytical.getGlobalDistR(indz, indr, indphi, side);
        globalDistZ = spaceCharge3DCalcAnalytical.getGlobalDistZ(indz, indr, indphi, side);
        globalDistRPhi = spaceCharge3DCalcAnalytical.getGlobalDistRPhi(indz, indr, indphi, side);
        globalDistRNum = spaceCharge3DCalcNumerical.getGlobalDistR(indz, indr, indphi, side);
        globalDistZNum = spaceCharge3DCalcNumerical.getGlobalDistZ(indz, indr, indphi, side);
        globalDistRPhiNum = spaceCharge3DCalcNumerical.getGlobalDistRPhi(indz, indr, indphi, side);
// std::cout<<"AAA"<<std::endl;
        // get global corrections
        globalCorrR = spaceCharge3DCalcAnalytical.getGlobalCorrR(indz, indr, indphi, side);
        globalCorrZ = spaceCharge3DCalcAnalytical.getGlobalCorrZ(indz, indr, indphi, side);
        globalCorrRPhi = spaceCharge3DCalcAnalytical.getGlobalCorrRPhi(indz, indr, indphi, side);
        globalCorrRNum = spaceCharge3DCalcNumerical.getGlobalCorrR(indz, indr, indphi, side);
        globalCorrZNum = spaceCharge3DCalcNumerical.getGlobalCorrZ(indz, indr, indphi, side);
        globalCorrRPhiNum = spaceCharge3DCalcNumerical.getGlobalCorrRPhi(indz, indr, indphi, side);

        // get numerically calculated Efield
        eZNum = spaceCharge3DCalcNumerical.getEz(indz, indr, indphi, side);
        eRNum = spaceCharge3DCalcNumerical.getEr(indz, indr, indphi, side);
        ePhiNum = spaceCharge3DCalcNumerical.getEphi(indz, indr, indphi, side);

        potNum = spaceCharge3DCalcNumerical.getPotential(indz, indr, indphi, side);
        potAna = spaceCharge3DCalcAnalytical.getPotential(indz, indr, indphi, side);

        // get analytical e field
        const DataT radius = spaceCharge3DCalcAnalytical.getRVertex(indr);
        r = radius;
        z = spaceCharge3DCalcAnalytical.getZVertex(indz);
        phi = spaceCharge3DCalcAnalytical.getPhiVertex(indphi);
        eZAna = anaFields.evalEz(z, radius, phi);
        eRAna = anaFields.evalEr(z, radius, phi);
        ePhiAna = anaFields.evalEphi(z, radius, phi);
        // std::cout<<"BBB"<<std::endl;

        densNum = spaceCharge3DCalcNumerical.getDensity(indz, indr, indphi, side);
// std::cout<<"1: "<< globalDistZNum << std::endl;
        // distort point and then correct it
        // const bool safe = true;
        const DataT posDistorted[3] = {z + globalDistZNum, radius + globalDistRNum, phi + globalDistRPhiNum / radius};
// std::cout<<"1a"<<std::endl;
// std::cout<<"posDistorted[0]: "<< posDistorted[0] << std::endl;
// std::cout<<"posDistorted[1]: "<< posDistorted[1] << std::endl;
// std::cout<<"posDistorted[2]: "<< posDistorted[2] << std::endl;

        const DataT corrRDistPointTmp = numGlobalCorrInterpolator.evaldR(posDistorted[0], posDistorted[1], posDistorted[2]);
// std::cout<<"2"<<std::endl;
        corrDistPoint[0] = corrRDistPointTmp;
        corrDistPoint[1] = numGlobalCorrInterpolator.evaldZ(posDistorted[0], posDistorted[1], posDistorted[2]);
        const DataT corrRPhiDistPointTmp = numGlobalCorrInterpolator.evaldRPhi(posDistorted[0], posDistorted[1], posDistorted[2]) / posDistorted[1];
      // std::cout<<"3"<<std::endl;
        corrDistPoint[2] = corrRPhiDistPointTmp * radius;
        zDistortedPoint = posDistorted[0];
        rDistortedPoint = posDistorted[1];
        phiDistortedPoint = posDistorted[2];
        // std::cout<<"CCC"<<std::endl;

        //local
        const DataT lposDistorted[3] = {z + localDistZNum + spaceCharge3DCalcNumerical.getGridSpacingZ(), radius + localDistRNum, phi + localDistRPhiNum / radius};
        const DataT lcorrRDistPointTmp = numLocalCorrInterpolator.evaldR(lposDistorted[0], lposDistorted[1], lposDistorted[2]);
        lcorrDistPoint[0] = lcorrRDistPointTmp;
        lcorrDistPoint[1] = numLocalCorrInterpolator.evaldZ(lposDistorted[0], lposDistorted[1], lposDistorted[2]);
        const DataT lcorrRPhiDistPointTmp = numGlobalCorrInterpolator.evaldRPhi(lposDistorted[0], lposDistorted[1], lposDistorted[2]) / lposDistorted[1];
        lcorrDistPoint[2] = lcorrRPhiDistPointTmp * radius;



        const DataT posCorrected[3] = {z + globalCorrZNum, radius + globalCorrRNum, phi + globalCorrRPhiNum / radius};
        correctionInVolume = 1;
        if (posCorrected[0] < spaceCharge3DCalcNumerical.getZMin() || posCorrected[0] > spaceCharge3DCalcNumerical.getZMax() || posCorrected[1] < spaceCharge3DCalcNumerical.getRMin() || posCorrected[1] > spaceCharge3DCalcNumerical.getRMax()) {
          correctionInVolume = 0;
        }

        //=========
        // const DataT posDistortedAna[3] = {z + globalDistZ, radius + globalDistR, phi + globalDistRPhi / radius};
        // std::cout<<"NOT REG: "<< phi + globalDistRPhi / radius << std::endl;
        // std::cout<<"REG: "<< spaceCharge3DCalcNumerical.regulatePhi(phi + globalDistRPhi / radius) << std::endl;

        const DataT posDistortedAna[3] = {z + globalDistZ, radius + globalDistR, spaceCharge3DCalcNumerical.regulatePhi(phi + globalDistRPhi / radius)};
        const DataT corrRDistPointTmpAna = anaGlobalCorrInterpolator.evaldR(posDistortedAna[0], posDistortedAna[1], posDistortedAna[2]);
        corrDistPointAna[0] = corrRDistPointTmpAna;
        corrDistPointAna[1] = anaGlobalCorrInterpolator.evaldZ(posDistortedAna[0], posDistortedAna[1], posDistortedAna[2]);
        const DataT corrRPhiDistPointTmpAna = anaGlobalCorrInterpolator.evaldRPhi(posDistortedAna[0], posDistortedAna[1], posDistortedAna[2]) / posDistortedAna[1];
        corrDistPointAna[2] = corrRPhiDistPointTmpAna * radius;

        zDistortedPoint = posDistortedAna[0];
        rDistortedPoint = posDistortedAna[1];
        phiDistortedPoint = posDistortedAna[2];

        // lcorr
        const DataT lposDistortedAna[3] = {z + localDistZ + spaceCharge3DCalcNumerical.getGridSpacingZ(), radius + localDistR, spaceCharge3DCalcNumerical.regulatePhi(phi + localDistRPhi / radius)};
        const DataT lcorrRDistPointTmpAna = anaLocalCorrInterpolator.evaldR(lposDistortedAna[0], lposDistortedAna[1], lposDistortedAna[2]);
        lcorrDistPointAna[0] = lcorrRDistPointTmpAna;
        lcorrDistPointAna[1] = anaLocalCorrInterpolator.evaldZ(lposDistortedAna[0], lposDistortedAna[1], lposDistortedAna[2]);
        const DataT lcorrRPhiDistPointTmpAna = anaLocalCorrInterpolator.evaldRPhi(lposDistortedAna[0], lposDistortedAna[1], lposDistortedAna[2]) / lposDistortedAna[1];
        lcorrDistPointAna[2] = lcorrRPhiDistPointTmpAna * radius;

        tree.Fill();
      }
    }
  }
  // treeInterpolation.Write();
  tree.Write();
  fDebug.Close();

  return 0;
}

void debug()
{
  const int nGridX = 10;
  const int nGridY = 10;
  const int nGridZ = 10;
  const float xmin = 1;
  const float ymin = 0;
  const float zmin = 0;
  const float spacingX = 1;
  const float spacingY = 1;
  const float spacingZ = 1;

  o2::tpc::RegularGrid3D<float, nGridX, nGridY, nGridZ> grid3D(xmin, ymin, zmin, spacingX, spacingY, spacingZ);

  float relPosX = 10.1;

  std::cout << "getGridMinX: " << grid3D.getGridMinX() << std::endl;
  std::cout << "getGridMaxX: " << grid3D.getGridMaxX() << std::endl;
  std::cout << "------------------------------------" << std::endl;

  Vector<float, 3> relPos({relPosX, relPosX, relPosX});
  Vector<int, 3> circular({0, 0, 1});

  grid3D.clampToGridRel(relPos, circular);

  std::cout << "relPos[0]:" << relPos[0] << std::endl;
  std::cout << "relPos[2]:" << relPos[2] << std::endl;

  // const unsigned int nGridR = 257;
  // const unsigned int nGridZZ = nGridR;
  // const unsigned int nGridPhi = 360;
  // O2TPCSpaceCharge3DCalc<float, nGridR, nGridZZ, nGridPhi> sC3D;
  //
  // std::cout << "Starting:" << std::endl;
  // // float clampedX = grid3D.clampToGrid(relPosX, 0);
  // std::cout << "sC3D.regulatePhi(0): " << sC3D.regulatePhi(0) << std::endl;
  // std::cout << "sC3D.regulatePhi(0): " << sC3D.regulatePhi(-1) << std::endl;
}

float tricubicSpeedTest()
{
  // define the grid:
  // define number of vertices per dimension
  const int xvertices = 129;
  const int yvertices = 129;
  const int zvertices = 180;

  // define min range
  // float xmin = 0;
  // float xmax = 10;
  // float ymin = 0;
  // float ymax = 10;
  // float zmin = 0;
  //
  // // define spacing between grid vertices
  // float xSpacing = (xmax - xmin) / (xvertices - 1);
  // float ySpacing = (ymax - ymin) / (yvertices - 1);
  // float zSpacing = 2 * M_PI / zvertices;

  // create grid and datacontainer object
  // o2::tpc::RegularGrid3D<float, xvertices, yvertices, zvertices> grid3D(xmin, ymin, zmin, xSpacing, ySpacing, zSpacing);

  o2::tpc::O2TPCSpaceCharge3DCalc<double, xvertices, yvertices, zvertices> spaceCharge3DCalc;
  o2::tpc::RegularGrid3D<double, xvertices, yvertices, zvertices> grid3D = spaceCharge3DCalc.getGrid3D();
  DataContainer3D<double, xvertices, yvertices, zvertices> data3D;

  o2::tpc::Side side = o2::tpc::Side::C;
  AnalyticalFields<double> anaFields(side);

  // fill the DataContainer3D with some values
  for (int ix = 0; ix < xvertices; ++ix) {
    for (int iy = 0; iy < yvertices; ++iy) {
      for (int iz = 0; iz < zvertices; ++iz) {
        const double ixPos = grid3D.getXVertex(ix);
        const double iyPos = grid3D.getYVertex(iy);
        const double izPos = grid3D.getZVertex(iz);
        // data3D(ix, iy, iz) = FUNC(ixPos, iyPos, izPos);
        data3D(ix, iy, iz) = anaFields.evalEr(ixPos, iyPos, izPos);
      }
    }
  }
  // define periodicity of dimension
  bool periodicX = false;
  bool periodicY = false;
  bool periodicZ = true;

  TFile fdebug("inter.root", "RECREATE");
  TTree tree("tree", "tree");
  float IX = 0;
  float IZ = 0;
  float IY = 0;
  int IIX = 0;
  int IIZ = 0;
  int IIY = 0;
  float val = 0;
  float trueval = 0;
  tree.Branch("x", &IX);
  tree.Branch("y", &IY);
  tree.Branch("phi", &IZ);
  tree.Branch("ix", &IIX);
  tree.Branch("iy", &IIY);
  tree.Branch("iz", &IIZ);
  tree.Branch("val", &val);
  tree.Branch("trueval", &trueval);

  // create tricubic interpolator
  o2::tpc::TriCubicInterpolator<double, xvertices, yvertices, zvertices> interpolator(data3D, grid3D, periodicX, periodicY, periodicZ);
  const bool catmullRom = true;
  gte::IntpTricubic3<double> tricubicInterpolator(xvertices, yvertices, zvertices, grid3D.getGridMinX(), grid3D.getSpacingX(), grid3D.getGridMinY(), grid3D.getSpacingY(), grid3D.getGridMinZ(), grid3D.getSpacingZ(), data3D.data(), catmullRom);
  // query some values
  int a = 0;

  interpolator.setExtrapolationType(o2::tpc::TriCubicInterpolator<double, xvertices, yvertices, zvertices>::ExtrapolationType::None);

  auto start = std::chrono::high_resolution_clock::now();
// #pragma omp parallel for
  for(int irep=0; irep<1; ++irep)
  for (int ix = 0; ix < xvertices; ++ix) {
    for (int iy = -3; iy < yvertices; ++iy) {
      for (int iz = 0; iz < zvertices - 1; ++iz) {
        for (int iix = 1; iix < 3; ++iix) {
          for (int iiy = 1; iiy < 3; ++iiy) {
            for (int iiz = 1; iiz < 3; ++iiz) {



              const float xQuery = grid3D.getXVertex(ix) + iix * grid3D.getSpacingX() / 3.;

              const float yVertex = iy>=0 ? grid3D.getYVertex(iy) : grid3D.getYVertex(0) + iy * grid3D.getSpacingY();
              const float yQuery = yVertex + iiy * grid3D.getSpacingY() / 3.;
              const float zQuery = grid3D.getZVertex(iz) + iiz * grid3D.getSpacingZ() / 3.;
              if (xQuery > grid3D.getGridMaxX() || yQuery > grid3D.getGridMaxY()) {
                continue;
              }



              const float interpolatedValue = interpolator(xQuery, yQuery, zQuery);
              a += interpolatedValue;


              // /*
              IIX = ix;
              IIY = iy;
              IIZ = iz;
              IX = xQuery;
              IY = yQuery;
              IZ = zQuery;
              val = interpolatedValue;
              // trueval = FUNC(xQuery, yQuery, zQuery);
              trueval = anaFields.evalEr(xQuery, yQuery, zQuery);//FUNC(xQuery, yQuery, zQuery);
              tree.Fill();
              // */

              // const float interpolatedValueGEO = tricubicInterpolator(xQuery, yQuery, zQuery);
              // a += interpolatedValueGEO;

              // const float interpolatedDerivative = interpolator(xQuery, yQuery, zQuery, 1, 1, 0);
              // const float trueDerivative = 1 / 10. * std::cos(yQuery * xQuery / 10.) - yQuery / 10. * std::sin(yQuery * xQuery / 10.) * xQuery / 10.;
            }
          }
        }
      }
    }
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<float> diff = end - start;
  std::cout << "time interpolator: " << diff.count() << std::endl;
  std::cout << "a: " << a << std::endl;
  tree.Write();
  return a;
}
