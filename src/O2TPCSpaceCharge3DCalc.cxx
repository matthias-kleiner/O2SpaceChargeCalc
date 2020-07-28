#include "include/O2TPCSpaceCharge3DCalc.h"

#include "TF1.h" /// for numerical intergration only

templateClassImp(O2TPCSpaceCharge3DCalc);
// const int nTHREADS = 8;

template <typename DataT, size_t Nr, size_t Nz, size_t Nphi>
void O2TPCSpaceCharge3DCalc<DataT, Nr, Nz, Nphi>::getDistortionsAnalytical(const DataT p1r, const DataT p1phi, const DataT p1z, const DataT p2z, DataT& ddR, DataT& ddRPhi, DataT& ddZ, Formulas<DataT> formulaStruct) const
{

  DataT localIntErOverEz = 0;
  DataT localIntEPhiOverEz = 0;
  DataT localIntDeltaEz = 0;
  const int integrationType = mNumericalIntegrationStrategy; /// 0->Trapezoidal rule (fast), 1->Simpson rule (slower than 1, but more precise),2->Root integration (slow)
  const DataT ezField = getEzField();

  if (integrationType == 2) {
    TF1 fErOverEz("fErOverEz", [&](double* x, double* p) { return static_cast<double>(formulaStruct.evalEr(p1r, p1phi, static_cast<DataT>(x[0]) ) / (formulaStruct.evalEz(p1r, p1phi, static_cast<DataT>(x[0])) + ezField)); }, p1z, p2z, 1);
    localIntErOverEz = static_cast<DataT>(fErOverEz.Integral(p1z, p2z));

    TF1 fEphiOverEz("fEPhiOverEz", [&](double* x, double* p) { return static_cast<double>(formulaStruct.evalEphi(p1r, p1phi, static_cast<DataT>(x[0])) / (formulaStruct.evalEz(p1r, p1phi, static_cast<DataT>(x[0])) + ezField)); }, p1z, p2z, 1);
    localIntEPhiOverEz = static_cast<DataT>(fEphiOverEz.Integral(p1z, p2z));

    TF1 fEz("fEZOverEz", [&](double* x, double* p) { return static_cast<double>(formulaStruct.evalEz(p1r, p1phi, static_cast<DataT>(x[0])) - ezField); }, p1z, p2z, 1);
    localIntDeltaEz = static_cast<DataT>(fEz.Integral(p1z, p2z));
  } else {
    const DataT fielder0 = formulaStruct.evalEr(p1r, p1phi, p1z);
    const DataT fieldez0 = formulaStruct.evalEz(p1r, p1phi, p1z);
    const DataT fieldephi0 = formulaStruct.evalEphi(p1r, p1phi, p1z);

    const DataT fielder1 = formulaStruct.evalEr(p1r, p1phi, p2z);
    const DataT fieldez1 = formulaStruct.evalEz(p1r, p1phi, p2z);
    const DataT fieldephi1 = formulaStruct.evalEphi(p1r, p1phi, p2z);

    const DataT eZ0 = ezField + fieldez0;
    const DataT eZ1 = ezField + fieldez1;

    const int nSteps = 1; //getIntegrationSteps(); //mNumericalIntegrationSteps;
    const DataT deltaX = (p2z - p1z) / nSteps;

    // trapezoidal integration
    if (integrationType == 0) {
      //========trapezoidal rule==============
      DataT fieldSumEr = 0;
      DataT fieldSumEphi = 0;
      DataT fieldSumEz = 0;
      for (int i = 1; i < nSteps; ++i) {
        const DataT xk1Tmp = p1z + i * deltaX;
        const DataT ezField1 = formulaStruct.evalEz(p1r, p1phi, xk1Tmp);
        const DataT ezField1Denominator = 1. / ezField + ezField1;

        fieldSumEr += formulaStruct.evalEr(p1r, p1phi, xk1Tmp) * ezField1Denominator;
        fieldSumEphi += formulaStruct.evalEphi(p1r, p1phi, xk1Tmp) * ezField1Denominator;
        fieldSumEz += ezField1;
      }
      localIntErOverEz = deltaX * (fieldSumEr + 0.5f * (fielder0 / eZ0 + fielder1 / eZ1));
      localIntEPhiOverEz = deltaX * (fieldSumEphi + 0.5f * (fieldephi0 / eZ0 + fieldephi1 / eZ1));
      localIntDeltaEz = deltaX * (fieldSumEz + 0.5f * (fieldez0 + fieldez1));
    } else if (integrationType == 1) {
      //==========simpsons rule see: https://en.wikipedia.org/wiki/Simpson%27s_rule =============================
      DataT fieldSum1ErOverEz = 0;
      DataT fieldSum2ErOverEz = 0;
      DataT fieldSum1EphiOverEz = 0;
      DataT fieldSum2EphiOverEz = 0;
      DataT fieldSum1Ez = 0;
      DataT fieldSum2Ez = 0;

      for (int i = 1; i < nSteps; ++i) {
        const DataT xk1Tmp = p1z + i * deltaX;
        const DataT xk2 = xk1Tmp - 0.5f * deltaX;

        const DataT ezField1 = formulaStruct.evalEz(p1r, p1phi, xk1Tmp);
        const DataT ezField2 = formulaStruct.evalEz(p1r, p1phi, xk2);
        const DataT ezField1Denominator = 1.f / (ezField + ezField1);
        const DataT ezField2Denominator = 1.f / (ezField + ezField2);

        fieldSum1ErOverEz += formulaStruct.evalEr(p1r, p1phi, xk1Tmp) * ezField1Denominator;
        fieldSum2ErOverEz += formulaStruct.evalEr(p1r, p1phi, xk2) * ezField2Denominator;

        fieldSum1EphiOverEz += formulaStruct.evalEphi(p1r, p1phi, xk1Tmp) * ezField1Denominator;
        fieldSum2EphiOverEz += formulaStruct.evalEphi(p1r, p1phi, xk2) * ezField2Denominator;

        fieldSum1Ez += ezField1;
        fieldSum2Ez += ezField2;
      }
      const DataT xk2N = (p2z - 0.5f * deltaX);
      const DataT ezField2 = formulaStruct.evalEz(p1r, p1phi, xk2N);
      const DataT ezField2Denominator = 1.f / (ezField + ezField2);
      fieldSum2ErOverEz += formulaStruct.evalEr(p1r, p1phi, xk2N) * ezField2Denominator;
      fieldSum2EphiOverEz += formulaStruct.evalEphi(p1r, p1phi, xk2N) * ezField2Denominator;
      fieldSum2Ez += ezField2;

      const DataT deltaXSimpsonSixth = deltaX / 6.f;
      localIntErOverEz = deltaXSimpsonSixth * (2 * fieldSum1ErOverEz + 4 * fieldSum2ErOverEz + fielder0 / eZ0 + fielder1 / eZ1);
      localIntEPhiOverEz = deltaXSimpsonSixth * (2 * fieldSum1EphiOverEz + 4 * fieldSum2EphiOverEz + fieldephi0 / eZ0 + fieldephi1 / eZ1);
      localIntDeltaEz = deltaXSimpsonSixth * (2 * fieldSum1Ez + 4 * fieldSum2Ez + fieldez0 + fieldez1);
    }
  }

  ddR = fC0 * localIntErOverEz + fC1 * localIntEPhiOverEz;
  ddRPhi = (fC0 * localIntEPhiOverEz - fC1 * localIntErOverEz) / p1r;
  ddZ = -1 * localIntDeltaEz * ASolv::fgkdvdE;
}

template class O2TPCSpaceCharge3DCalc<float, 17, 17, 90>;
template class O2TPCSpaceCharge3DCalc<double, 17, 17, 90>;
