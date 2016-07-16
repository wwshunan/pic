
#include "TBeam.h"

using namespace std;

TBeam::TBeam() {}

TBeam::~TBeam() {}

void TBeam::TTrans(vector<double> &oX, vector<double> &oXp, double oXEmit,
                   double oXAlpha, double oXBeta) {
  double oXEmitBetaSqrt = sqrt(oXEmit / oXBeta);
  double oXEmitBeta = sqrt(oXEmit * oXBeta);
  vector<double> tXXAlpha = oX * oXAlpha;
  oXp = (oXp - tXXAlpha);
  oXp = oXp * oXEmitBetaSqrt;

  oX = oX * oXEmitBeta;
}

void TBeam::GS6D(int numPart, vector<double> &tX, vector<double> &tXp,
                 vector<double> &tY, vector<double> &tYp, vector<double> &tZ,
                 vector<double> &tZp) {

  for (auto &iNumPart : tX) {
    iNumPart = oRandom3.Gaus();
  }
  for (auto &iNumPart : tXp) {
    iNumPart = oRandom3.Gaus();
  }
  for (auto &iNumPart : tY) {
    iNumPart = oRandom3.Gaus();
  }
  for (auto &iNumPart : tYp) {
    iNumPart = oRandom3.Gaus();
  }
  for (auto &iNumPart : tZ) {
    iNumPart = oRandom3.Gaus();
  }
  for (auto &iNumPart : tZp) {
    iNumPart = oRandom3.Gaus();
  }
}

void TBeam::WB6D(int numPart, vector<double> &tX, vector<double> &tXp,
                 vector<double> &tY, vector<double> &tYp, vector<double> &tZ,
                 vector<double> &tZp) {

  int nNumPart = 0;
  double oX, oXp, oY, oYp, oZ, oZp, oS;
  while (nNumPart < numPart) {
    oX = oRandom3.Uniform(-1, 1);
    oXp = oRandom3.Uniform(-1, 1);
    oY = oRandom3.Uniform(-1, 1);
    oYp = oRandom3.Uniform(-1, 1);
    oZ = oRandom3.Uniform(-1, 1);
    oZp = oRandom3.Uniform(-1, 1);
    oS = oX * oX + oXp * oXp + oY * oY + oYp * oYp + oZ * oZ + oZp * oZp;
    if (oS > 1) {
      continue;
    }
    tX[nNumPart] = oX;
    tXp[nNumPart] = oXp;
    tY[nNumPart] = oY;
    tYp[nNumPart] = oYp;
    tZ[nNumPart] = oZ;
    tZp[nNumPart] = oZp;

    nNumPart += 1;
  }
}
void TBeam::GS4D(int numPart, vector<double> &tX, vector<double> &tXp,
                 vector<double> &tY, vector<double> &tYp) {
  for (auto &iNumPart : tX) {
    iNumPart = oRandom3.Gaus();
  }
  for (auto &iNumPart : tXp) {
    iNumPart = oRandom3.Gaus();
  }
  for (auto &iNumPart : tY) {
    iNumPart = oRandom3.Gaus();
  }
  for (auto &iNumPart : tYp) {
    iNumPart = oRandom3.Gaus();
  }
}
void TBeam::WB4D(int numPart, vector<double> &tX, vector<double> &tXp,
                 vector<double> &tY, vector<double> &tYp) {
  int nNumPart = 0;
  double oX, oXp, oY, oYp, oS;
  while (nNumPart < numPart) {
    oX = oRandom3.Uniform(-1, 1);
    oXp = oRandom3.Uniform(-1, 1);
    oY = oRandom3.Uniform(-1, 1);
    oYp = oRandom3.Uniform(-1, 1);

    oS = oX * oX + oXp * oXp + oY * oY + oYp * oYp;
    if (oS > 1) {
      continue;
    }
    tX[nNumPart] = oX;
    tXp[nNumPart] = oXp;
    tY[nNumPart] = oY;
    tYp[nNumPart] = oYp;

    nNumPart += 1;
  }
}
void TBeam::KV4D(int numPart, vector<double> &tX, vector<double> &tXp,
                 vector<double> &tY, vector<double> &tYp) {
  double oX, oXp, oY, oYp, oS1, oS2;
  int nNumPart = 0;
  while (nNumPart < numPart) {
    oX = oRandom3.Uniform(-1, 1);
    oXp = oRandom3.Uniform(-1, 1);
    oS1 = oX * oX + oXp * oXp;
    if (oS1 > 1) {
      continue;
    }
    oS2 = 2;
    while (oS2 > 1) {
      oY = oRandom3.Uniform(-1, 1);
      oYp = oRandom3.Uniform(-1, 1);
      oS2 = oY * oY + oYp * oYp;
    }
    tX[nNumPart] = oX;
    tXp[nNumPart] = oXp;
    double oArg = sqrt((1 - oS1) / oS2);
    tY[nNumPart] = oY * oArg;
    tYp[nNumPart] = oYp * oArg;
    nNumPart += 1;
  }
}
void TBeam::GS2D(int numPart, vector<double> &tR, vector<double> &tRp) {
  for (auto &iNumPart : tR) {
    iNumPart = oRandom3.Gaus();
  }
  for (auto &iNumPart : tRp) {
    iNumPart = oRandom3.Gaus();
  }
}
void TBeam::KV2D(int numPart, vector<double> &tR, vector<double> &tRp) {
  int nNumPart = 0;
  double oR, oRp, oS;
  while (nNumPart < numPart) {
    oR = oRandom3.Uniform(-1, 1);
    oRp = oRandom3.Uniform(-1, 1);

    oS = oR * oR + oRp * oRp;
    if (oS > 1) {
      continue;
    }
    tR[nNumPart] = oR;
    tRp[nNumPart] = oRp;

    nNumPart += 1;
  }
}

void TBeam::UN3D(int numPart, vector<double> &tX, vector<double> &tY,
                 vector<double> &tZ) {
  int nNumPart = 0;
  double oX, oY, oZ, oS;
  while (nNumPart < numPart) {
    oX = oRandom3.Uniform(-1, 1);
    oY = oRandom3.Uniform(-1, 1);
    oZ = oRandom3.Uniform(-1, 1);

    oS = oX * oX + oY * oY + oZ * oZ;
    if (oS > 1) {
      continue;
    }
    tX[nNumPart] = oX;
    tY[nNumPart] = oY;
    tZ[nNumPart] = oZ;

    nNumPart += 1;
  }
}

void TBeam::PartGen() {
  string oDistribution;
  fNumPart = 0;
  int numPartRec = 0;
  string numPartStr;
  string xEmitStr, yEmitStr, zEmitStr;
  string xAlphaStr, yAlphaStr, zAlphaStr;
  string xBetaStr, yBetaStr, zBetaStr;
  string tChargeStr, tMassStr;

  double numPartDbl;
  double xEmitDbl, yEmitDbl, zEmitDbl;
  double xAlphaDbl, yAlphaDbl, zAlphaDbl;
  double xBetaDbl, yBetaDbl, zBetaDbl;
  double tChargeDbl, tMassDbl;

  vector<double> tX, tY, tZ, tXp, tYp, tZp;

  // sum of all charge state particle, and malloc the space.
  for (int iMultiBeam = 0; iMultiBeam < fMultiBeam; ++iMultiBeam) {
    if (iMultiBeam < 10) {
      numPartStr = (string) "numpart_" + (char)(iMultiBeam + 49);
    } else {
      stringstream oSstream;
      oSstream.str("");
      oSstream << iMultiBeam;
      numPartStr = (string) "numpart_" + oSstream.str();
    }

    fNumPart = fNumPart + (int)fAcceDbl[numPartStr];
  }

  fXPartDistri.resize(fNumPart);
  fXpPartDistri.resize(fNumPart);
  fYPartDistri.resize(fNumPart);
  fYpPartDistri.resize(fNumPart);
  fZPartDistri.resize(fNumPart);
  fZpPartDistri.resize(fNumPart);

  fFlagPartDistri.resize(fNumPart);
  fChargePartDistri.resize(fNumPart);
  fMassPartDistri.resize(fNumPart);

  // particle generation~

  for (int iMultiBeam = 0; iMultiBeam < fMultiBeam; ++iMultiBeam) {
    if (iMultiBeam < 10) {
      oDistribution = (string) "distribution_" + (char)(iMultiBeam + 49);
      numPartStr = (string) "numpart_" + (char)(iMultiBeam + 49);
      xEmitStr = (string) "xepsilon_" + (char)(iMultiBeam + 49);
      yEmitStr = (string) "yepsilon_" + (char)(iMultiBeam + 49);
      zEmitStr = (string) "zepsilon_" + (char)(iMultiBeam + 49);
      xAlphaStr = (string) "xalpha_" + (char)(iMultiBeam + 49);
      yAlphaStr = (string) "yalpha_" + (char)(iMultiBeam + 49);
      zAlphaStr = (string) "zalpha_" + (char)(iMultiBeam + 49);
      xBetaStr = (string) "xbeta_" + (char)(iMultiBeam + 49);
      yBetaStr = (string) "ybeta_" + (char)(iMultiBeam + 49);
      zBetaStr = (string) "zbeta_" + (char)(iMultiBeam + 49);
      tChargeStr = (string) "charge_" + (char)(iMultiBeam + 49);
      tMassStr = (string) "mass_" + (char)(iMultiBeam + 49);
    } else {
      stringstream oSstream;
      oSstream.str("");
      oSstream << iMultiBeam;
      oDistribution = (string) "distribution_" + oSstream.str();
      numPartStr = (string) "numpart_" + oSstream.str();
      xEmitStr = (string) "xepsilon_" + oSstream.str();
      yEmitStr = (string) "yepsilon_" + oSstream.str();
      zEmitStr = (string) "zepsilon_" + oSstream.str();
      xAlphaStr = (string) "xalpha_" + oSstream.str();
      yAlphaStr = (string) "yalpha_" + oSstream.str();
      zAlphaStr = (string) "zalpha_" + oSstream.str();
      xBetaStr = (string) "xbeta_" + oSstream.str();
      yBetaStr = (string) "ybeta_" + oSstream.str();
      zBetaStr = (string) "zbeta_" + oSstream.str();
      tChargeStr = (string) "charge_" + oSstream.str();
      tMassStr = (string) "mass_" + oSstream.str();
    }

    numPartDbl = fAcceDbl[numPartStr];

    xEmitDbl = fAcceDbl[xEmitStr];
    yEmitDbl = fAcceDbl[yEmitStr];
    zEmitDbl = fAcceDbl[zEmitStr];
    xAlphaDbl = fAcceDbl[xAlphaStr];
    yAlphaDbl = fAcceDbl[yAlphaStr];
    zAlphaDbl = fAcceDbl[zAlphaStr];
    xBetaDbl = fAcceDbl[xBetaStr];
    yBetaDbl = fAcceDbl[yBetaStr];
    zBetaDbl = fAcceDbl[zBetaStr];

    tX.resize(numPartDbl);
    tY.resize(numPartDbl);
    tZ.resize(numPartDbl);
    tXp.resize(numPartDbl);
    tYp.resize(numPartDbl);
    tZp.resize(numPartDbl);

    if (fAcceStr[oDistribution] == "gs6d") {
      GS6D(numPartDbl, tX, tXp, tY, tYp, tZ, tZp);
    } else if (fAcceStr[oDistribution] == "wb6d") {
      WB6D(numPartDbl, tX, tXp, tY, tYp, tZ, tZp);
    } else if (fAcceStr[oDistribution] == "kv4d") {
      KV4D(numPartDbl, tX, tXp, tY, tYp);
    } else if (fAcceStr[oDistribution] == "gs4d") {
      GS4D(numPartDbl, tX, tXp, tY, tYp);
    } else if (fAcceStr[oDistribution] == "wb4d") {
      WB4D(numPartDbl, tX, tXp, tY, tYp);
    } else if (fAcceStr[oDistribution] == "kv2d") {
      KV2D(numPartDbl, tX, tXp);
      KV2D(numPartDbl, tY, tYp);
    } else if (fAcceStr[oDistribution] == "gs2d") {
      GS2D(numPartDbl, tX, tXp);
      GS2D(numPartDbl, tY, tYp);
    } else if (fAcceStr[oDistribution] == "un3d") {
      UN3D(numPartDbl, tX, tY, tZ);
    }

    TTrans(tX, tXp, xEmitDbl, xAlphaDbl, xBetaDbl);
    TTrans(tY, tYp, yEmitDbl, yAlphaDbl, yBetaDbl);
    TTrans(tZ, tZp, zEmitDbl, zAlphaDbl, zBetaDbl);

    copy(tX.begin(), tX.end(), fXPartDistri.begin() + numPartRec);
    copy(tXp.begin(), tXp.end(), fXpPartDistri.begin() + numPartRec);
    copy(tY.begin(), tY.end(), fYPartDistri.begin() + numPartRec);
    copy(tYp.begin(), tYp.end(), fYpPartDistri.begin() + numPartRec);
    copy(tZ.begin(), tZ.end(), fZPartDistri.begin() + numPartRec);
    copy(tZp.begin(), tZp.end(), fZpPartDistri.begin() + numPartRec);

    tChargeDbl = fAcceDbl[tChargeStr];
    tMassDbl = fAcceDbl[tMassStr];
    for (auto iNumPart = fChargePartDistri.begin() + numPartRec;
         iNumPart < fChargePartDistri.begin() + numPartRec + numPartDbl;
         ++iNumPart) {
      *iNumPart = tChargeDbl;
    }
    for (auto iNumPart = fMassPartDistri.begin() + numPartRec;
         iNumPart < fMassPartDistri.begin() + numPartRec + numPartDbl;
         ++iNumPart) {
      *iNumPart = tMassDbl;
    }

    numPartRec = numPartRec + (int)fAcceDbl[numPartStr];
  }
}
//***************************************************************************************
// void TBeam::MeshPoisson2D(double xMinDomain, double xMaxDomain,
//                           double yMinDomain, double yMaxDomain, int xGridLog,
//                           int yGridLog) {
//
//   xGrid = pow(2, xGridLog);
//   yGrid = pow(2, yGridLog);
//
//   qGrid.resize(xGrid);
//   for (auto &iQGrid : qGrid) {
//     iQGrid.resize(yGrid);
//   }
//
//   double xDiffDomain = xMaxDomain - xMinDomain;
//   double yDiffDomain = yMaxDomain - yMinDomain;
//   double xStepDomain = xDiffDomain / (xGrid + 1);
//   double yStepDomain = yDiffDomain / (yGrid + 1);
//
//   vector<double> mXVec = fXPartDistri - xMinDomain;
//   mXVec = mXVec / xStepDomain;
//   vector<double> nYVec = fYPartDistri - yMinDomain;
//   nYVec = nYVec / yStepDomain;
//
//   double mXDbl, nYDbl, mX1Frac, mX2Frac, nY1Frac, nY2Frac;
//   int m1XInt, n1YInt, m2XInt, n2YInt;
//
//   for (int iNumPart = 0; iNumPart < fNumPart; ++iNumPart) {
//     if (!fFlagPartDistri[iNumPart]) {
//       mXDbl = mXVec[iNumPart];
//       nYDbl = nYVec[iNumPart];
//       if ((mXDbl < 0) | (nYDbl < 0) | (mXDbl > xGrid + 1) |
//           (nYDbl > yGrid + 1)) {
//         fFlagPartDistri[iNumPart] = fZPartDistri[iNumPart];
//         continue;
//       } else {
//         m2XInt = mXDbl;
//         n2YInt = nYDbl;
//         m1XInt = m2XInt - 1;
//         n1YInt = n2YInt - 1;
//         mX1Frac = mXDbl - m2XInt;
//         nY1Frac = nYDbl - n2YInt;
//         mX2Frac = 1 - mX1Frac;
//         nY2Frac = 1 - nY1Frac;
//
//         if ((m1XInt == -1) & (n1YInt == -1)) {
//           qGrid[m2XInt][n2YInt] +=
//               fChargePartDistri[iNumPart] * mX1Frac * nY1Frac;
//         } else if ((m1XInt == -1) & (n1YInt == yGrid - 1)) {
//           qGrid[m2XInt][n1YInt] +=
//               fChargePartDistri[iNumPart] * mX1Frac * nY2Frac;
//         } else if ((m1XInt == xGrid - 1) & (n1YInt == -1)) {
//           qGrid[m1XInt][n2YInt] +=
//               fChargePartDistri[iNumPart] * mX2Frac * nY1Frac;
//         } else if ((m1XInt == xGrid - 1) & (n1YInt == yGrid - 1)) {
//           qGrid[m1XInt][n1YInt] +=
//               fChargePartDistri[iNumPart] * mX2Frac * nY2Frac;
//         } else if (m1XInt == -1) {
//           qGrid[m2XInt][n2YInt] +=
//               fChargePartDistri[iNumPart] * mX1Frac * nY1Frac;
//           qGrid[m2XInt][n1YInt] +=
//               fChargePartDistri[iNumPart] * mX1Frac * nY2Frac;
//         } else if (m1XInt == xGrid - 1) {
//           qGrid[m1XInt][n2YInt] +=
//               fChargePartDistri[iNumPart] * mX2Frac * nY1Frac;
//           qGrid[m1XInt][n1YInt] +=
//               fChargePartDistri[iNumPart] * mX2Frac * nY2Frac;
//         } else if (n1YInt == -1) {
//           qGrid[m2XInt][n2YInt] +=
//               fChargePartDistri[iNumPart] * mX1Frac * nY1Frac;
//           qGrid[m1XInt][n2YInt] +=
//               fChargePartDistri[iNumPart] * mX2Frac * nY1Frac;
//         } else if (n1YInt == yGrid - 1) {
//           qGrid[m2XInt][n1YInt] +=
//               fChargePartDistri[iNumPart] * mX1Frac * nY2Frac;
//           qGrid[m1XInt][n1YInt] +=
//               fChargePartDistri[iNumPart] * mX2Frac * nY2Frac;
//         } else {
//           qGrid[m2XInt][n1YInt] +=
//               fChargePartDistri[iNumPart] * mX1Frac * nY2Frac;
//           qGrid[m1XInt][n1YInt] +=
//               fChargePartDistri[iNumPart] * mX2Frac * nY2Frac;
//           qGrid[m2XInt][n2YInt] +=
//               fChargePartDistri[iNumPart] * mX1Frac * nY1Frac;
//           qGrid[m1XInt][n2YInt] +=
//               fChargePartDistri[iNumPart] * mX2Frac * nY1Frac;
//         }
//       }
//     }
//   }
// }

void TBeam::MeshFFT2D(double xMinDomain, double xMaxDomain, double yMinDomain,
                      double yMaxDomain, int xGridLog, int yGridLog) {

  xGrid = pow(2, xGridLog);
  yGrid = pow(2, yGridLog);

  qGrid.resize(xGrid * yGrid);

  double xDiffDomain = xMaxDomain - xMinDomain;
  double yDiffDomain = yMaxDomain - yMinDomain;
  xStepDomain = xDiffDomain / (xGrid + 1);
  yStepDomain = yDiffDomain / (yGrid + 1);

  vector<double> mXVec = fXPartDistri - xMinDomain;
  mXVec = mXVec / xStepDomain;
  vector<double> nYVec = fYPartDistri - yMinDomain;
  nYVec = nYVec / yStepDomain;

  double mXDbl, nYDbl, mX1Frac, mX2Frac, nY1Frac, nY2Frac;
  int m1XInt, n1YInt, m2XInt, n2YInt;

  for (int iNumPart = 0; iNumPart < fNumPart; ++iNumPart) {
    if (!fFlagPartDistri[iNumPart]) {
      mXDbl = mXVec[iNumPart];
      nYDbl = nYVec[iNumPart];
      if ((mXDbl < 0) | (nYDbl < 0) | (mXDbl > xGrid + 1) |
          (nYDbl > yGrid + 1)) {
        fFlagPartDistri[iNumPart] = fZPartDistri[iNumPart];
        continue;
      } else {
        m2XInt = mXDbl;
        n2YInt = nYDbl;
        m1XInt = m2XInt - 1;
        n1YInt = n2YInt - 1;
        mX1Frac = mXDbl - m2XInt;
        nY1Frac = nYDbl - n2YInt;
        mX2Frac = 1 - mX1Frac;
        nY2Frac = 1 - nY1Frac;

        if ((m1XInt == -1) & (n1YInt == -1)) {
          qGrid[m2XInt * yGrid + n2YInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY1Frac;
        } else if ((m1XInt == -1) & (n1YInt == yGrid - 1)) {
          qGrid[m2XInt * yGrid + n1YInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY2Frac;
        } else if ((m1XInt == xGrid - 1) & (n1YInt == -1)) {
          qGrid[m1XInt * yGrid + n2YInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY1Frac;
        } else if ((m1XInt == xGrid - 1) & (n1YInt == yGrid - 1)) {
          qGrid[m1XInt * yGrid + n1YInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY2Frac;
        } else if (m1XInt == -1) {
          qGrid[m2XInt * yGrid + n2YInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY1Frac;
          qGrid[m2XInt * yGrid + n1YInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY2Frac;
        } else if (m1XInt == xGrid - 1) {
          qGrid[m1XInt * yGrid + n2YInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY1Frac;
          qGrid[m1XInt * yGrid + n1YInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY2Frac;
        } else if (n1YInt == -1) {
          qGrid[m2XInt * yGrid + n2YInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY1Frac;
          qGrid[m1XInt * yGrid + n2YInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY1Frac;
        } else if (n1YInt == yGrid - 1) {
          qGrid[m2XInt * yGrid + n1YInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY2Frac;
          qGrid[m1XInt * yGrid + n1YInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY2Frac;
        } else {
          qGrid[m2XInt * yGrid + n1YInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY2Frac;
          qGrid[m1XInt * yGrid + n1YInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY2Frac;
          qGrid[m2XInt * yGrid + n2YInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY1Frac;
          qGrid[m1XInt * yGrid + n2YInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY1Frac;
        }
      }
    }
  }
  double dArea = xStepDomain * yStepDomain;
  qGrid = qGrid / dArea;
}

void TBeam::FFTPoisson2D() {
  vector<double> q2uGrid;
  q2uGrid.resize(xGrid * yGrid);
  fftw_plan plan2D;
  plan2D = fftw_plan_r2r_2d(xGrid, yGrid, &qGrid[0], &q2uGrid[0], FFTW_RODFT00,
                            FFTW_RODFT00, FFTW_MEASURE);
  fftw_execute(plan2D);

  vector<double> oKY, oKYTemp, oK;
  double dblKX;
  oKY.resize(yGrid);
  oKYTemp.resize(yGrid);
  oK.resize(xGrid * yGrid);

  for (int iKY = 0; iKY < yGrid; ++iKY) {
    oKY[iKY] = -pow(2 * sin(PI / 2 * (iKY + 1) / (yGrid + 1)) / yStepDomain, 2);
  }

  for (int iKX = 0; iKX < xGrid; ++iKX) {
    dblKX = -pow(2 * sin(PI / 2 * (iKX + 1) / (xGrid + 1)) / xStepDomain, 2);
    oKYTemp = oKY + dblKX;
    copy(oKYTemp.begin(), oKYTemp.end(), oK.begin() + iKX * yGrid);
  }

  q2uGrid = q2uGrid / oK;

  uGrid.resize(xGrid * yGrid);
  plan2D = fftw_plan_r2r_2d(xGrid, yGrid, &q2uGrid[0], &uGrid[0], FFTW_RODFT00,
                            FFTW_RODFT00, FFTW_ESTIMATE);
  fftw_execute(plan2D);
  fftw_destroy_plan(plan2D);

  double ScaledRatio = 1. / (4 * (xGrid + 1) * (yGrid + 1));
  uGrid = uGrid * ScaledRatio;

  // Plot3D(uGrid, yGrid, "uGrid", "uGrid", "AP");
  //
  // vector<double> tGrid;
  // tGrid.resize(xGrid * yGrid);
  // tGrid = uGrid / qGrid;
  // TInput oInput;
  // oInput.PrintVec<vector<double>>(tGrid);
}

void TBeam::MeshFFT3D(double xMinDomain, double xMaxDomain, double yMinDomain,
                      double yMaxDomain, double zMinDomain, double zMaxDomain,
                      int xGridLog, int yGridLog, int zGridLog) {
  //

  xGrid = pow(2, xGridLog);
  yGrid = pow(2, yGridLog);
  zGrid = pow(2, zGridLog);

  qGrid.resize(xGrid * yGrid * zGrid);
  for (auto &iQGrid : qGrid)
    iQGrid = 0;

  double xDiffDomain = xMaxDomain - xMinDomain;
  double yDiffDomain = yMaxDomain - yMinDomain;
  double zDiffDomain = zMaxDomain - zMinDomain;
  xStepDomain = xDiffDomain / (xGrid + 1);
  yStepDomain = yDiffDomain / (yGrid + 1);
  zStepDomain = zDiffDomain / (zGrid + 1);

  vector<double> mXVec = fXPartDistri - xMinDomain;
  mXVec = mXVec / xStepDomain;
  vector<double> nYVec = fYPartDistri - yMinDomain;
  nYVec = nYVec / yStepDomain;
  vector<double> kZVec = fZPartDistri - zMinDomain;
  kZVec = kZVec / zStepDomain;

  double mXDbl, nYDbl, kZDbl, mX1Frac, mX2Frac, nY1Frac, nY2Frac, kZ1Frac,
      kZ2Frac;
  int m1XInt, n1YInt, m2XInt, n2YInt, k1ZInt, k2ZInt;

  for (int iNumPart = 0; iNumPart < fNumPart; ++iNumPart) {
    if (!fFlagPartDistri[iNumPart]) {
      mXDbl = mXVec[iNumPart];
      nYDbl = nYVec[iNumPart];
      kZDbl = kZVec[iNumPart];

      if ((mXDbl < 0) | (nYDbl < 0) | (kZDbl < 0) | (mXDbl > xGrid + 1) |
          (nYDbl > yGrid + 1) | (kZDbl > zGrid + 1)) {
        fFlagPartDistri[iNumPart] = fZPartDistri[iNumPart];
        continue;
      } else {
        m2XInt = mXDbl;
        n2YInt = nYDbl;
        k2ZInt = kZDbl;
        m1XInt = m2XInt - 1;
        n1YInt = n2YInt - 1;
        k1ZInt = k2ZInt - 1;
        mX1Frac = mXDbl - m2XInt;
        nY1Frac = nYDbl - n2YInt;
        kZ1Frac = kZDbl - k2ZInt;
        mX2Frac = 1 - mX1Frac;
        nY2Frac = 1 - nY1Frac;
        kZ2Frac = 1 - kZ1Frac;

        if ((m1XInt == -1) & (n1YInt == -1) & (k1ZInt == -1)) {
          qGrid[m2XInt * (yGrid * zGrid) + n2YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY1Frac * kZ1Frac;
        } else if ((m1XInt == -1) & (n1YInt == -1) & (k1ZInt == zGrid - 1)) {
          qGrid[m2XInt * (yGrid * zGrid) + n2YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY1Frac * kZ2Frac;
        } else if ((m1XInt == -1) & (n1YInt == yGrid - 1) & (k1ZInt == -1)) {
          qGrid[m2XInt * (yGrid * zGrid) + n1YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY2Frac * kZ1Frac;
        } else if ((m1XInt == -1) & (n1YInt == yGrid - 1) &
                   (k1ZInt == zGrid - 1)) {
          qGrid[m2XInt * (yGrid * zGrid) + n1YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY2Frac * kZ2Frac;
        } else if ((m1XInt == xGrid - 1) & (n1YInt == -1) & (k1ZInt == -1)) {
          qGrid[m1XInt * (yGrid * zGrid) + n2YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY1Frac * kZ1Frac;
        } else if ((m1XInt == xGrid - 1) & (n1YInt == -1) &
                   (k1ZInt == zGrid - 1)) {
          qGrid[m1XInt * (yGrid * zGrid) + n2YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY1Frac * kZ2Frac;
        } else if ((m1XInt == xGrid - 1) & (n1YInt == yGrid - 1) &
                   (k1ZInt == -1)) {
          qGrid[m1XInt * (yGrid * zGrid) + n1YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY2Frac * kZ1Frac;
        } else if ((m1XInt == xGrid - 1) & (n1YInt == yGrid - 1) &
                   (k1ZInt == zGrid - 1)) {
          qGrid[m1XInt * (yGrid * zGrid) + n1YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY2Frac * kZ2Frac;
        } else if ((m1XInt == -1) & (n1YInt == -1)) {
          qGrid[m2XInt * (yGrid * zGrid) + n2YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY1Frac * kZ1Frac;
          qGrid[m2XInt * (yGrid * zGrid) + n2YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY1Frac * kZ2Frac;
        } else if ((m1XInt == -1) & (n1YInt == yGrid - 1)) {
          qGrid[m2XInt * (yGrid * zGrid) + n1YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY2Frac * kZ1Frac;
          qGrid[m2XInt * (yGrid * zGrid) + n1YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY2Frac * kZ2Frac;
        } else if ((m1XInt == xGrid - 1) & (n1YInt == -1)) {
          qGrid[m1XInt * (yGrid * zGrid) + n2YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY1Frac * kZ1Frac;
          qGrid[m1XInt * (yGrid * zGrid) + n2YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY1Frac * kZ2Frac;
        } else if ((m1XInt == xGrid - 1) & (n1YInt == yGrid - 1)) {
          qGrid[m1XInt * (yGrid * zGrid) + n1YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY2Frac * kZ1Frac;
          qGrid[m1XInt * (yGrid * zGrid) + n1YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY2Frac * kZ2Frac;
        } else if ((m1XInt == -1) & (k1ZInt == -1)) {
          qGrid[m2XInt * (yGrid * zGrid) + n2YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY1Frac * kZ1Frac;
          qGrid[m2XInt * (yGrid * zGrid) + n1YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY2Frac * kZ1Frac;
        } else if ((m1XInt == -1) & (k1ZInt == zGrid - 1)) {
          qGrid[m2XInt * (yGrid * zGrid) + n2YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY1Frac * kZ2Frac;
          qGrid[m2XInt * (yGrid * zGrid) + n1YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY2Frac * kZ2Frac;
        } else if ((m1XInt == xGrid - 1) & (k1ZInt == -1)) {
          qGrid[m1XInt * (yGrid * zGrid) + n2YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY1Frac * kZ1Frac;
          qGrid[m1XInt * (yGrid * zGrid) + n1YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY2Frac * kZ1Frac;
        } else if ((m1XInt == xGrid - 1) & (k1ZInt == zGrid - 1)) {
          qGrid[m1XInt * (yGrid * zGrid) + n2YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY1Frac * kZ2Frac;
          qGrid[m1XInt * (yGrid * zGrid) + n1YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY2Frac * kZ2Frac;
        } else if ((n1YInt == -1) & (k1ZInt == -1)) {
          qGrid[m2XInt * (yGrid * zGrid) + n2YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY1Frac * kZ1Frac;
          qGrid[m1XInt * (yGrid * zGrid) + n2YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY1Frac * kZ1Frac;
        } else if ((n1YInt == -1) & (k1ZInt == zGrid - 1)) {
          qGrid[m2XInt * (yGrid * zGrid) + n2YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY1Frac * kZ2Frac;
          qGrid[m1XInt * (yGrid * zGrid) + n2YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY1Frac * kZ2Frac;
        } else if ((n1YInt == yGrid - 1) & (k1ZInt == -1)) {
          qGrid[m2XInt * (yGrid * zGrid) + n1YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY2Frac * kZ1Frac;
          qGrid[m1XInt * (yGrid * zGrid) + n1YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY2Frac * kZ1Frac;
        } else if ((n1YInt == yGrid - 1) & (k1ZInt == zGrid - 1)) {
          qGrid[m2XInt * (yGrid * zGrid) + n1YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY2Frac * kZ2Frac;
          qGrid[m1XInt * (yGrid * zGrid) + n1YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY2Frac * kZ2Frac;
        } else if (m1XInt == -1) {
          qGrid[m2XInt * (yGrid * zGrid) + n2YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY1Frac * kZ1Frac;
          qGrid[m2XInt * (yGrid * zGrid) + n2YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY1Frac * kZ2Frac;
          qGrid[m2XInt * (yGrid * zGrid) + n1YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY2Frac * kZ1Frac;
          qGrid[m2XInt * (yGrid * zGrid) + n1YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY2Frac * kZ2Frac;
        } else if (m1XInt == xGrid - 1) {
          qGrid[m1XInt * (yGrid * zGrid) + n2YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY1Frac * kZ1Frac;
          qGrid[m1XInt * (yGrid * zGrid) + n2YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY1Frac * kZ2Frac;
          qGrid[m1XInt * (yGrid * zGrid) + n1YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY2Frac * kZ1Frac;
          qGrid[m1XInt * (yGrid * zGrid) + n1YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY2Frac * kZ2Frac;
        } else if (n1YInt == -1) {
          qGrid[m2XInt * (yGrid * zGrid) + n2YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY1Frac * kZ1Frac;
          qGrid[m2XInt * (yGrid * zGrid) + n2YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY1Frac * kZ2Frac;
          qGrid[m1XInt * (yGrid * zGrid) + n2YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY1Frac * kZ1Frac;
          qGrid[m1XInt * (yGrid * zGrid) + n2YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY1Frac * kZ2Frac;
        } else if (n1YInt == yGrid - 1) {
          qGrid[m2XInt * (yGrid * zGrid) + n1YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY2Frac * kZ1Frac;
          qGrid[m2XInt * (yGrid * zGrid) + n1YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY2Frac * kZ2Frac;
          qGrid[m1XInt * (yGrid * zGrid) + n1YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY2Frac * kZ1Frac;
          qGrid[m1XInt * (yGrid * zGrid) + n1YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY2Frac * kZ2Frac;
        } else if (k1ZInt == -1) {
          qGrid[m2XInt * (yGrid * zGrid) + n2YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY1Frac * kZ1Frac;
          qGrid[m2XInt * (yGrid * zGrid) + n1YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY2Frac * kZ1Frac;
          qGrid[m1XInt * (yGrid * zGrid) + n2YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY1Frac * kZ1Frac;
          qGrid[m1XInt * (yGrid * zGrid) + n1YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY2Frac * kZ1Frac;
        } else if (k1ZInt == zGrid - 1) {
          qGrid[m2XInt * (yGrid * zGrid) + n2YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY1Frac * kZ2Frac;
          qGrid[m2XInt * (yGrid * zGrid) + n1YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY2Frac * kZ2Frac;
          qGrid[m1XInt * (yGrid * zGrid) + n2YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY1Frac * kZ2Frac;
          qGrid[m1XInt * (yGrid * zGrid) + n1YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY2Frac * kZ2Frac;
        } else {
          qGrid[m2XInt * (yGrid * zGrid) + n2YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY1Frac * kZ1Frac;
          qGrid[m2XInt * (yGrid * zGrid) + n2YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY1Frac * kZ2Frac;
          qGrid[m2XInt * (yGrid * zGrid) + n1YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY2Frac * kZ1Frac;
          qGrid[m2XInt * (yGrid * zGrid) + n1YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX1Frac * nY2Frac * kZ2Frac;
          qGrid[m1XInt * (yGrid * zGrid) + n2YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY1Frac * kZ1Frac;
          qGrid[m1XInt * (yGrid * zGrid) + n2YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY1Frac * kZ2Frac;
          qGrid[m1XInt * (yGrid * zGrid) + n1YInt * zGrid + k2ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY2Frac * kZ1Frac;
          qGrid[m1XInt * (yGrid * zGrid) + n1YInt * zGrid + k1ZInt] +=
              fChargePartDistri[iNumPart] * mX2Frac * nY2Frac * kZ2Frac;
        }
      }
    }
  }

  double dArea = xStepDomain * yStepDomain * zStepDomain;
  qGrid = qGrid / dArea;
}

void TBeam::Laplace2D() {
  vector<double> qGridLaplace;
  qGridLaplace.resize(xGrid * yGrid);

  for (int iXGrid = 0; iXGrid != xGrid; ++iXGrid) {
    for (int iYGrid = 0; iYGrid != yGrid; ++iYGrid) {
      if (iXGrid == 0 && iYGrid == 0)
        qGridLaplace[iXGrid * yGrid + iYGrid] =
            (uGrid[iXGrid * yGrid + iYGrid + 1] -
             2 * uGrid[iXGrid * yGrid + iYGrid]) /
                pow(xStepDomain, 2) +
            (uGrid[(iXGrid + 1) * yGrid + iYGrid] -
             2 * uGrid[iXGrid * yGrid + iYGrid]) /
                pow(yStepDomain, 2);
      else if (iXGrid == 0 && iYGrid == yGrid - 1)
        qGridLaplace[iXGrid * yGrid + iYGrid] =
            (uGrid[iXGrid * yGrid + iYGrid + 1] -
             2 * uGrid[iXGrid * yGrid + iYGrid]) /
                pow(xStepDomain, 2) +
            (uGrid[(iXGrid - 1) * yGrid + iYGrid] -
             2 * uGrid[iXGrid * yGrid + iYGrid]) /
                pow(yStepDomain, 2);
      else if (iXGrid == xGrid - 1 && iYGrid == 0)
        qGridLaplace[iXGrid * yGrid + iYGrid] =
            (uGrid[(iXGrid - 1) * yGrid + iYGrid] -
             2 * uGrid[iXGrid * yGrid + iYGrid]) /
                pow(xStepDomain, 2) +
            (uGrid[iXGrid * yGrid + iYGrid + 1] -
             2 * uGrid[iXGrid * yGrid + iYGrid]) /
                pow(yStepDomain, 2);
      else if (iXGrid == xGrid - 1 && iYGrid == yGrid - 1)
        qGridLaplace[iXGrid * yGrid + iYGrid] =
            (uGrid[(iXGrid - 1) * yGrid + iYGrid] -
             2 * uGrid[iXGrid * yGrid + iYGrid]) /
                pow(xStepDomain, 2) +
            (uGrid[iXGrid * yGrid + iYGrid - 1] -
             2 * uGrid[iXGrid * yGrid + iYGrid]) /
                pow(yStepDomain, 2);
      else if (iYGrid == 0)
        qGridLaplace[iXGrid * yGrid + iYGrid] =
            (uGrid[(iXGrid + 1) * yGrid + iYGrid] +
             uGrid[(iXGrid - 1) * yGrid + iYGrid] -
             2 * uGrid[iXGrid * yGrid + iYGrid]) /
                pow(xStepDomain, 2) +
            (uGrid[iXGrid * yGrid + iYGrid + 1] -
             2 * uGrid[iXGrid * yGrid + iYGrid]) /
                pow(yStepDomain, 2);
      else if (iXGrid == 0)
        qGridLaplace[iXGrid * yGrid + iYGrid] =
            (uGrid[(iXGrid + 1) * yGrid + iYGrid] -
             2 * uGrid[iXGrid * yGrid + iYGrid]) /
                pow(xStepDomain, 2) +
            (uGrid[iXGrid * yGrid + iYGrid + 1] +
             uGrid[iXGrid * yGrid + iYGrid - 1] -
             2 * uGrid[iXGrid * yGrid + iYGrid]) /
                pow(yStepDomain, 2);
      else if (iYGrid == yGrid - 1)
        qGridLaplace[iXGrid * yGrid + iYGrid] =
            (uGrid[(iXGrid + 1) * yGrid + iYGrid] +
             uGrid[(iXGrid - 1) * yGrid + iYGrid] -
             2 * uGrid[iXGrid * yGrid + iYGrid]) /
                pow(xStepDomain, 2) +
            (uGrid[iXGrid * yGrid + iYGrid - 1] -
             2 * uGrid[iXGrid * yGrid + iYGrid]) /
                pow(yStepDomain, 2);
      else if (iXGrid == xGrid - 1)
        qGridLaplace[iXGrid * yGrid + iYGrid] =
            (uGrid[(iXGrid - 1) * yGrid + iYGrid] -
             2 * uGrid[iXGrid * yGrid + iYGrid]) /
                pow(xStepDomain, 2) +
            (uGrid[iXGrid * yGrid + iYGrid + 1] +
             uGrid[iXGrid * yGrid + iYGrid - 1] -
             2 * uGrid[iXGrid * yGrid + iYGrid]) /
                pow(yStepDomain, 2);
      else
        qGridLaplace[iXGrid * yGrid + iYGrid] =
            (uGrid[(iXGrid + 1) * yGrid + iYGrid] +
             uGrid[(iXGrid - 1) * yGrid + iYGrid] -
             2 * uGrid[iXGrid * yGrid + iYGrid]) /
                pow(xStepDomain, 2) +
            (uGrid[iXGrid * yGrid + iYGrid + 1] +
             uGrid[iXGrid * yGrid + iYGrid - 1] -
             2 * uGrid[iXGrid * yGrid + iYGrid]) /
                pow(yStepDomain, 2);
    }
  }

  // // Plot3D(qGridLaplace, yGrid);
  // double sumqGridLaplace = 0;
  // for (auto &iSumqGridLaplace :
  // qGridLaplace)
  //   sumqGridLaplace += iSumqGridLaplace;
  // cout << sumqGridLaplace << endl;
}

void TBeam::Gradient2D() {
  xEGrid.resize(xGrid * yGrid);
  yEGrid.resize(xGrid * yGrid);

  for (int iXGrid = 0; iXGrid != xGrid; ++iXGrid) {
    for (int iYGrid = 0; iYGrid != yGrid; ++iYGrid) {
      if (iYGrid == 0)
        yEGrid[iXGrid * yGrid + iYGrid] =
            (uGrid[iXGrid * yGrid + iYGrid + 1]) / (2 * yStepDomain);
      else if (iYGrid == yGrid - 1)
        yEGrid[iXGrid * yGrid + iYGrid] =
            -(uGrid[iXGrid * yGrid + iYGrid - 1]) / (2 * yStepDomain);
      else
        yEGrid[iXGrid * yGrid + iYGrid] = (uGrid[iXGrid * yGrid + iYGrid + 1] -
                                           uGrid[iXGrid * yGrid + iYGrid - 1]) /
                                          (2 * yStepDomain);
      if (iXGrid == 0)
        xEGrid[iXGrid * yGrid + iYGrid] =
            (uGrid[(iXGrid + 1) * yGrid + iYGrid]) / (2 * xStepDomain);
      else if (iXGrid == xGrid - 1)
        xEGrid[iXGrid * yGrid + iYGrid] =
            -(uGrid[(iXGrid - 1) * yGrid + iYGrid]) / (2 * xStepDomain);
      else
        xEGrid[iXGrid * yGrid + iYGrid] =
            (uGrid[(iXGrid + 1) * yGrid + iYGrid] -
             uGrid[(iXGrid - 1) * yGrid + iYGrid]) /
            (2 * xStepDomain);
    }
  }
  // Plot3D(xEGrid, yGrid, "xEGrid", "xEGrid");
  // oCanvas2D["xEGrid"]->WaitPrimitive();
}

void TBeam::InterPoint2D() {
  //
}

void TBeam::FFTPoisson3D() {

  vector<double> q2uGrid(xGrid * yGrid * zGrid, 0);
  uGrid.resize(xGrid * yGrid * zGrid);

  fftw_plan plan_r2hc, plan_hc2r;
  plan_r2hc =
      fftw_plan_r2r_3d(xGrid, yGrid, zGrid, &qGrid[0], &q2uGrid[0],
                       FFTW_RODFT00, FFTW_RODFT00, FFTW_R2HC, FFTW_MEASURE);
  plan_hc2r =
      fftw_plan_r2r_3d(xGrid, yGrid, zGrid, &q2uGrid[0], &uGrid[0],
                       FFTW_RODFT00, FFTW_RODFT00, FFTW_HC2R, FFTW_MEASURE);

  fftw_execute(plan_r2hc);
  int iQ2UGrid;
  double oKx, oKy, oKz, oK;
  for (int iXGrid = 0; iXGrid < xGrid; iXGrid++) {
    for (int iYGrid = 0; iYGrid < yGrid; iYGrid++) {
      for (int iZGrid = 0; iZGrid < zGrid; iZGrid++) {
        iQ2UGrid = iXGrid * (yGrid * zGrid) + iYGrid * zGrid + iZGrid;
        oKx =
            -pow(2 * sin(PI / 2 * (iXGrid + 1) / (xGrid + 1)) / xStepDomain, 2);
        oKy =
            -pow(2 * sin(PI / 2 * (iYGrid + 1) / (yGrid + 1)) / yStepDomain, 2);
        oKz = -pow(2 * sin(PI * iZGrid / zGrid) / zStepDomain, 2);
        oK = oKx + oKy + oKz;
        q2uGrid[iQ2UGrid] = q2uGrid[iQ2UGrid] / oK;
      }
    }
  }
  fftw_execute(plan_hc2r);

  double ScaledRatio = 1. / (4 * (xGrid + 1) * (yGrid + 1) * zGrid);
  uGrid = uGrid * ScaledRatio;

  fftw_destroy_plan(plan_r2hc);
  fftw_destroy_plan(plan_hc2r);

  double maxUGrid = -*min_element(uGrid.begin(), uGrid.end());

  cout << "U Max: " << maxUGrid << endl;

  vector<double> xyUGrid(xGrid * yGrid, 0);
  vector<double> xzUGrid(xGrid * zGrid, 0);
  vector<double> yzUGrid(yGrid * zGrid, 0);

  int xChoose, yChoose, zChoose;
  xChoose = xGrid / 2;
  yChoose = yGrid / 2;
  zChoose = zGrid / 2;

  int iYZ = 0, iXZ = 0, iXY = 0;
  for (int iXGrid = 0; iXGrid < xGrid; iXGrid++) {
    for (int iYGrid = 0; iYGrid < yGrid; iYGrid++) {
      for (int iZGrid = 0; iZGrid < zGrid; iZGrid++) {
        if (iXGrid == xChoose) {
          iQ2UGrid = iXGrid * (yGrid * zGrid) + iYGrid * zGrid + iZGrid;
          yzUGrid[iYZ] = uGrid[iQ2UGrid];
          ++iYZ;
        }
        if (iYGrid == yChoose) {
          iQ2UGrid = iXGrid * (yGrid * zGrid) + iYGrid * zGrid + iZGrid;
          xzUGrid[iXZ] = uGrid[iQ2UGrid];
          ++iXZ;
        }
        if (iZGrid == zChoose) {
          iQ2UGrid = iXGrid * (yGrid * zGrid) + iYGrid * zGrid + iZGrid;
          xyUGrid[iXY] = uGrid[iQ2UGrid];
          ++iXY;
        }
      }
    }
  }

  TInput oInput;
  oInput.Plot3D(xyUGrid, yGrid, "xyUGrid", "xyUGrid", "PA");
  oInput.Plot3D(xzUGrid, zGrid, "xzUGrid", "xzUGrid", "PA");
  oInput.Plot3D(yzUGrid, zGrid, "yzUGrid", "yzUGrid", "PA");

  //  //  //  //  //  //  //  //  //  //  //  //  //









}
void TBeam::FTDPoisson2D() {}
void TBeam::FTDPoisson3D() {}





void TBeam::Restriction(vector<double>& resGridCoarse, vector<double>& resGridFine, int xGrid, int yGrid, int zGrid)
{
    int xGridFine = 2 * (xGrid - 1) - 1;
    int yGridFine = 2 * (yGrid - 1) - 1;
    int zGridFine = 2 * (zGrid - 1) - 1;
    int xGridFineWrap = xGridFine + 2;
    int yGridFineWrap = yGridFine + 2;
    int zGridFineWrap = zGridFine + 2;
    int iXGridFine, iYGridFine, iZGridFine;
    int iXGridFineWrap, iYGridFineWrap, iZGridFineWrap;
    double res0XGrid, res2XGrid, res0YGrid, res2YGrid, res0ZGrid, res2ZGrid, res1Grid;
    double vertex, arris, planeCenter;
    //vector<double> res(27);
    vector<double> resGridFineWrap((xGridFine + 2) * (yGridFine + 2) * (zGridFine + 2));

    for (int iXGrid=1; iXGrid != xGridFineWrap - 1; ++iXGrid)
        for (int iYGrid=1; iYGrid != yGridFineWrap - 1; ++iYGrid)
            for (int iZGrid=1; iZGrid != zGridFineWrap - 1; ++iZGrid) { 
                resGridFineWrap[iXGrid * yGridFineWrap * zGridFineWrap + iYGrid * zGridFineWrap + iZGrid] = resGridFine[(iXGrid - 1) * yGridFine * zGridFine + (iYGrid - 1) * zGridFine + iZGrid - 1];

    for (int iXGrid=0; iXGrid != xGrid; ++iXGrid)
        for (int iYGrid=0; iYGrid != yGrid; ++iYGrid)
            for (int iZGrid=0; iZGrid != zGrid; ++iZGrid) { 
                i1XGridFineWrap = 2 * iXGrid + 1;
                i1YGridFineWrap = 2 * iYGrid + 1;
                i1ZGridFineWrap = 2 * iZGrid + 1;
                i0XGridFineWrap = i1XGridFineWrap - 1;
                i0YGridFineWrap = i1YGridFineWrap - 1;
                i0ZGridFineWrap = i1ZGridFineWrap - 1;
                i2XGridFineWrap = i2XGridFineWrap + 1;
                i2YGridFineWrap = i2YGridFineWrap + 1;
                i2ZGridFineWrap = i2ZGridFineWrap + 1;
                vertex = resGridFineWrap[i0XGridFineWrap * yGridFineWrap * zGridFineWrap + i0YGridFineWrap * zGridFineWrap + i2ZGridFineWrap] + resGridFineWrap[i2XGridFineWrap * yGridFineWrap * zGridFineWrap + i0YGridFineWrap * zGridFineWrap + i2ZGridFineWrap] + resGridFineWrap[i2XGridFineWrap * yGridFineWrap * zGridFineWrap + i2YGridFineWrap * zGridFineWrap + i2ZGridFineWrap] + resGridFineWrap[i0XGridFineWrap * yGridFineWrap * zGridFineWrap + i2YGridFineWrap * zGridFineWrap + i2ZGridFineWrap] + resGridFineWrap[i0XGridFineWrap * yGridFineWrap * zGridFineWrap + i0YGridFineWrap * zGridFineWrap + i0ZGridFineWrap] + resGridFineWrap[i2XGridFineWrap * yGridFineWrap * zGridFineWrap + i0YGridFineWrap * zGridFineWrap + i0ZGridFineWrap] + resGridFineWrap[i2XGridFineWrap * yGridFineWrap * zGridFineWrap + i2YGridFineWrap * zGridFineWrap + i0ZGridFineWrap] + resGridFineWrap[i0XGridFineWrap * yGridFineWrap * zGridFineWrap + i2YGridFineWrap * zGridFineWrap + i0ZGridFineWrap];
                
                arris = resGridFineWrap[i1XGridFineWrap * yGridFineWrap * zGridFineWrap + i0YGridFineWrap * zGridFineWrap + i2ZGridFineWrap] + resGridFineWrap[i2XGridFineWrap * yGridFineWrap * zGridFineWrap + i1YGridFineWrap * zGridFineWrap + i2ZGridFineWrap] + resGridFineWrap[i1XGridFineWrap * yGridFineWrap * zGridFineWrap + i2YGridFineWrap * zGridFineWrap + i2ZGridFineWrap] + resGridFineWrap[i0XGridFineWrap * yGridFineWrap * zGridFineWrap + i1YGridFineWrap * zGridFineWrap + i2ZGridFineWrap] + resGridFineWrap[i1XGridFineWrap * yGridFineWrap * zGridFineWrap + i0YGridFineWrap * zGridFineWrap + i0ZGridFineWrap] + resGridFineWrap[i2XGridFineWrap * yGridFineWrap * zGridFineWrap + i1YGridFineWrap * zGridFineWrap + i0ZGridFineWrap] + resGridFineWrap[i1XGridFineWrap * yGridFineWrap * zGridFineWrap + i2YGridFineWrap * zGridFineWrap + i0ZGridFineWrap] + resGridFineWrap[i0XGridFineWrap * yGridFineWrap * zGridFineWrap + i1YGridFineWrap * zGridFineWrap + i0ZGridFineWrap] + resGridFineWrap[i0XGridFineWrap * yGridFineWrap * zGridFineWrap + i0YGridFineWrap * zGridFineWrap + i1ZGridFineWrap] + resGridFineWrap[i2XGridFineWrap * yGridFineWrap * zGridFineWrap + i0YGridFineWrap * zGridFineWrap + i1ZGridFineWrap] + resGridFineWrap[i2XGridFineWrap * yGridFineWrap * zGridFineWrap + i2YGridFineWrap * zGridFineWrap + i1ZGridFineWrap] + resGridFineWrap[i0XGridFineWrap * yGridFineWrap * zGridFineWrap + i2YGridFineWrap * zGridFineWrap + i1ZGridFineWrap]; 

                planeCenter = resGridFineWrap[i1XGridFineWrap * yGridFineWrap * zGridFineWrap + i0YGridFineWrap * zGridFineWrap + i1ZGridFineWrap] + resGridFineWrap[i2XGridFineWrap * yGridFineWrap * zGridFineWrap + i1YGridFineWrap * zGridFineWrap + i1ZGridFineWrap] + resGridFineWrap[i1XGridFineWrap * yGridFineWrap * zGridFineWrap + i2YGridFineWrap * zGridFineWrap + i1ZGridFineWrap] + resGridFineWrap[i0XGridFineWrap * yGridFineWrap * zGridFineWrap + i1YGridFineWrap * zGridFineWrap + i1ZGridFineWrap] + resGridFineWrap[i1XGridFineWrap * yGridFineWrap * zGridFineWrap + i1YGridFineWrap * zGridFineWrap + i2ZGridFineWrap] + resGridFineWrap[i1XGridFineWrap * yGridFineWrap * zGridFineWrap + i1YGridFineWrap * zGridFineWrap + i0ZGridFineWrap]; 

                selfPoint = resGridFineWrap[i1XGridFineWrap * yGridFineWrap * zGridFineWrap + i1YGridFineWrap * zGridFineWrap + i1ZGridFineWrap];
                resGridCoarse[iXGrid * xGrid * yGrid + iYGrid * zGrid + iZGrid] = vertex / 8 + arris / 4 + planeCenter / 2 +                  /*
    for (int iXGrid=0; iXGrid != xGrid; ++iXGrid)
        for (int iYGrid=0; iYGrid != yGrid; ++iYGrid)
            for (int iZGrid=0; iZGrid != zGrid; ++iZGrid) { 
                iXGridFine = 2 * iXGrid;
                iYGridFine = 2 * iYGrid;
                iZGridFine = 2 * iZGrid;
                i0XGrid = iXGridFine - 1;
                i0YGrid = iYGridFine - 1;
                i0ZGrid = iZGridFine - 1;
                i2XGrid = iXGridFine + 1;
                i2YGrid = iYGridFine + 1;
                i2ZGrid = iZGridFine + 1;

                if (i0XGrid > 0)
                    res[4] = resGridFine[i0XGrid * yGridFine * zGridFine + iYGridFine * zGridFine + iZGrid];

                if (i0YGrid > 0)
                    res[10] = resGridFine[iXGrid * yGridFine * zGridFine + i0YGrid * zGridFine + iZGrid];

                if (i0XGrid > 0 && i0YGrid > 0)
                    res[1] = resGridFine[i0XGrid * yGridFine * zGridFine + i0YGrid * zGridFine + iZGrid];
                if (i0ZGrid < 0)
                    res0ZGrid =
                        resGridFine[iXGrid * (yGridFine * zGridFine) + iYGrid * zGridFine + zGridFine - 2];
                else
                    res0ZGrid =
                        resGridFine[iXGrid * (yGridFine * zGridFine) + iYGrid * zGridFine + i0ZGrid];

                if (i2XGrid >= xGridFine)
                    res2XGrid = 0;
                else
                    res2XGrid =
                        resGridFine[i2XGrid * (yGridFine * zGridFine) + iYGrid * zGridFine + iZGrid];

                if (i2YGrid >= yGridFine)
                    res2YGrid = 0;
                else
                    res2YGrid =
                        resGridFine[iXGrid * (yGridFine * zGridFine) + i2YGrid * zGridFine + iZGrid];

                if (i2ZGrid >= zGridFine)
                    res2ZGrid = resGridFine[iXGrid * (yGridFine * zGridFine) + iYGrid * zGridFine + 1];
                else
                    res2ZGrid =
                        resGridFine[iXGrid * (yGridFine * zGridFine) + iYGrid * zGridFine + i2ZGrid];

                resGridCoarse[iXGrid * yGrid * zGrid + iYGrid * zGrid + iZGrid] = 
                */

}
void TBeam::FTDMGPoisson3D()
{
    int gridNum = xGrid * yGrid * zGrid;
    int x2hGrid = (xGrid - 1) / 2 + 1;
    int y2hGrid = (yGrid - 1) / 2 + 1;
    int z2hGrid = (zGrid - 1) / 2 + 1;
    int x4hGrid = (xGrid - 1) / 4 + 1;
    int y4hGrid = (yGrid - 1) / 4 + 1;
    int z4hGrid = (zGrid - 1) / 4 + 1;
    int x8hGrid = (xGrid - 1) / 8 + 1;
    int y8hGrid = (yGrid - 1) / 8 + 1;
    int z8hGrid = (zGrid - 1) / 8 + 1;
    vector<double> uGrid(gridNum);
    vector<double> u2hGrid(gridNum / 8);
    vector<double> u4hGrid(gridNum / 64);
    vector<double> u8hGrid(gridNum / 512);
    vector<double> resGrid;

    FTDSorPoisson3D(x8hGrid, y8hGrid, z8hGrid, u8hGrid, xStepDomain * 8, yStepDomain * 8 , zStepDomain * 8, 1.5, 1000);
    Prolongation(u8hGrid, u4hGrid, x4hGrid, y4hGrid, z4hGrid);
    FTDSorPoisson3D(x4hGrid, y4hGrid, z4hGrid, u4hGrid, xStepDomain * 4, yStepDomain * 4 , zStepDomain * 4, 1, 3);

    resGrid.resize(x4hGrid * y4hGrid * z4hGrid);
    CalcResGrid(u4hGrid, resGrid, x4hGrid, y4hGrid, z4hGrid, 4);



}
void TBeam::FTDSorPoisson2D() {}
void TBeam::FieldSpaceCharge2D() {}
void TBeam::FieldSpaceCharge3D() {}

void TBeam::Laplace3D() {
  vector<double> qGridLaplace;
  qGridLaplace.resize(xGrid * yGrid * zGrid);

  int i0XGrid, i2XGrid, i0YGrid, i2YGrid, i0ZGrid, i2ZGrid;
  double d0XGrid, d2XGrid, d0YGrid, d2YGrid, d0ZGrid, d2ZGrid, d1Grid;
  for (int i1XGrid = 0; i1XGrid != xGrid; ++i1XGrid) {
    for (int i1YGrid = 0; i1YGrid != yGrid; ++i1YGrid) {
      for (int i1ZGrid = 0; i1ZGrid != zGrid; ++i1ZGrid) {

        i0XGrid = i1XGrid - 1;
        i0YGrid = i1YGrid - 1;
        i0ZGrid = i1ZGrid - 1;

        i2XGrid = i1XGrid + 1;
        i2YGrid = i1YGrid + 1;
        i2ZGrid = i1ZGrid + 1;

        d1Grid = uGrid[i1XGrid * (yGrid * zGrid) + i1YGrid * zGrid + i1ZGrid];
        if (i0XGrid < 0)
          d0XGrid = 0;
        else
          d0XGrid =
              uGrid[i0XGrid * (yGrid * zGrid) + i1YGrid * zGrid + i1ZGrid];

        if (i0YGrid < 0)
          d0YGrid = 0;
        else
          d0YGrid =
              uGrid[i1XGrid * (yGrid * zGrid) + i0YGrid * zGrid + i1ZGrid];

        if (i0ZGrid < 0)
          d0ZGrid =
              uGrid[i1XGrid * (yGrid * zGrid) + i1YGrid * zGrid + zGrid - 2];
        else
          d0ZGrid =
              uGrid[i1XGrid * (yGrid * zGrid) + i1YGrid * zGrid + i0ZGrid];

        if (i2XGrid >= xGrid)
          d2XGrid = 0;
        else
          d2XGrid =
              uGrid[i2XGrid * (yGrid * zGrid) + i1YGrid * zGrid + i1ZGrid];

        if (i2YGrid >= yGrid)
          d2YGrid = 0;
        else
          d2YGrid =
              uGrid[i1XGrid * (yGrid * zGrid) + i2YGrid * zGrid + i1ZGrid];

        if (i2ZGrid >= zGrid)
          d2ZGrid = uGrid[i1XGrid * (yGrid * zGrid) + i1YGrid * zGrid + 1];
        else
          d2ZGrid =
              uGrid[i1XGrid * (yGrid * zGrid) + i1YGrid * zGrid + i2ZGrid];

        qGridLaplace[i1XGrid * (yGrid * zGrid) + i1YGrid * zGrid + i1ZGrid] =
            (d2XGrid + d0XGrid - 2 * d1Grid) / pow(xStepDomain, 2) +
            (d2YGrid + d0YGrid - 2 * d1Grid) / pow(yStepDomain, 2) +
            (d2ZGrid + d0ZGrid - 2 * d1Grid) / pow(zStepDomain, 2);
      }
    }
  }

  // // Plot3D(qGridLaplace, yGrid);
  double sumqGridLaplace = 0;
  for (auto &iSumqGridLaplace : qGridLaplace)
    sumqGridLaplace += iSumqGridLaplace;
  cout << sumqGridLaplace << endl;
}

void TBeam::Gradient3D() {
  xEGrid.resize(xGrid * yGrid * zGrid);
  yEGrid.resize(xGrid * yGrid * zGrid);
  zEGrid.resize(xGrid * yGrid * zGrid);

  for (auto &iXEGrid : xEGrid)
    iXEGrid = 0;

  for (auto &iYEGrid : yEGrid)
    iYEGrid = 0;

  for (auto &iZEGrid : zEGrid)
    iZEGrid = 0;

  int i0XGrid, i2XGrid, i0YGrid, i2YGrid, i0ZGrid, i2ZGrid;
  double d0XGrid, d2XGrid, d0YGrid, d2YGrid, d0ZGrid, d2ZGrid;
  for (int i1XGrid = 0; i1XGrid != xGrid; ++i1XGrid) {
    for (int i1YGrid = 0; i1YGrid != yGrid; ++i1YGrid) {
      for (int i1ZGrid = 0; i1ZGrid != zGrid; ++i1ZGrid) {

        i0XGrid = i1XGrid - 1;
        i0YGrid = i1YGrid - 1;
        i0ZGrid = i1ZGrid - 1;

        i2XGrid = i1XGrid + 1;
        i2YGrid = i1YGrid + 1;
        i2ZGrid = i1ZGrid + 1;

        if (i0XGrid < 0)
          d0XGrid = 0;
        else
          d0XGrid =
              uGrid[i0XGrid * (yGrid * zGrid) + i1YGrid * zGrid + i1ZGrid];

        if (i0YGrid < 0)
          d0YGrid = 0;
        else
          d0YGrid =
              uGrid[i1XGrid * (yGrid * zGrid) + i0YGrid * zGrid + i1ZGrid];

        if (i0ZGrid < 0)
          d0ZGrid =
              uGrid[i1XGrid * (yGrid * zGrid) + i1YGrid * zGrid + zGrid - 2];
        else
          d0ZGrid =
              uGrid[i1XGrid * (yGrid * zGrid) + i1YGrid * zGrid + i0ZGrid];

        if (i2XGrid >= xGrid)
          d2XGrid = 0;
        else
          d2XGrid =
              uGrid[i2XGrid * (yGrid * zGrid) + i1YGrid * zGrid + i1ZGrid];

        if (i2YGrid >= yGrid)
          d2YGrid = 0;
        else
          d2YGrid =
              uGrid[i1XGrid * (yGrid * zGrid) + i2YGrid * zGrid + i1ZGrid];

        if (i2ZGrid >= zGrid)
          d2ZGrid = uGrid[i1XGrid * (yGrid * zGrid) + i1YGrid * zGrid + 1];
        else
          d2ZGrid =
              uGrid[i1XGrid * (yGrid * zGrid) + i1YGrid * zGrid + i2ZGrid];

        xEGrid[i1XGrid * (yGrid * zGrid) + i1YGrid * zGrid + i1ZGrid] =
            (d2XGrid - d0XGrid) / (xStepDomain * 2);
        yEGrid[i1XGrid * (yGrid * zGrid) + i1YGrid * zGrid + i1ZGrid] =
            (d2YGrid - d0YGrid) / (yStepDomain * 2);
        zEGrid[i1XGrid * (yGrid * zGrid) + i1YGrid * zGrid + i1ZGrid] =
            (d2ZGrid - d0ZGrid) / (zStepDomain * 2);
      }
    }
  }

  double maxEGrid = *max_element(xEGrid.begin(), xEGrid.end());
  double maxX = *max_element(fXPartDistri.begin(), fXPartDistri.end());

  cout << "maxEGrid  |  " << maxEGrid << "  ++  "
       << " maxX  |  " << maxX << endl;
}

void TBeam::MeshDFT3D(double xMinDomain, double xMaxDomain, double yMinDomain,
                      double yMaxDomain, double zMinDomain, double zMaxDomain,
                      int xGridLog, int yGridLog, int zGridLog) {
  //

    xGrid = pow(2, xGridLog) + 1;
    yGrid = pow(2, yGridLog) + 1;
    zGrid = pow(2, zGridLog) + 1;

    qGrid.resize(xGrid * yGrid * zGrid);
    for (auto &iQGrid : qGrid)
        iQGrid = 0;

    double xDiffDomain = xMaxDomain - xMinDomain;
    double yDiffDomain = yMaxDomain - yMinDomain;
    double zDiffDomain = zMaxDomain - zMinDomain;
    xStepDomain = xDiffDomain / (xGrid - 1);
    yStepDomain = yDiffDomain / (yGrid - 1);
    zStepDomain = zDiffDomain / (zGrid - 1);

    vector<double> mXVec, nYVec, kZVec;
    mXVec.resize(xGrid * yGrid * zGrid);
    nYVec.resize(xGrid * yGrid * zGrid);
    kZVec.resize(xGrid * yGrid * zGrid);

    int vecSz = mXVec.size();
    for (int iMXVec=0; i != vecSz; ++iMXVec) {
        mXVec[iMXVec] = (fXpPartDistri[iMXVec] - xMinDomain) / xStepDomain;
        nYVec[iMXVec] = (fYpPartDistri[iMXVec] - yMinDomain) / yStepDomain;
        kZVec[iMXVec] = (fZpPartDistri[iMXVec] - zMinDomain) / zStepDomain;
    }

    double mXDbl, nYDbl, kZDbl, mX1Frac, mX2Frac, nY1Frac, nY2Frac, kZ1Frac,
           kZ2Frac;
    int m1XInt, n1YInt, m2XInt, n2YInt, k1ZInt, k2ZInt;

    for (int iNumPart = 0; iNumPart < fNumPart; ++iNumPart) {
        if (!fFlagPartDistri[iNumPart]) {
            mXDbl = mXVec[iNumPart];
            nYDbl = nYVec[iNumPart];
            kZDbl = kZVec[iNumPart];

            if ((mXDbl < 0) | (nYDbl < 0) | (kZDbl < 0) | (mXDbl > xGrid - 1) |
                    (nYDbl > yGrid - 1) | (kZDbl > zGrid - 1)) {
                fFlagPartDistri[iNumPart] = fZPartDistri[iNumPart];
                continue;
            } else {
                m2XInt = mXDbl + 1;
                n2YInt = nYDbl + 1;
                k2ZInt = kZDbl + 1;
                m1XInt = m2XInt - 1;
                n1YInt = n2YInt - 1;
                k1ZInt = k2ZInt - 1;
                mX1Frac = mXDbl - m1XInt;
                nY1Frac = nYDbl - n1YInt;
                kZ1Frac = kZDbl - k1ZInt;
                mX2Frac = 1 - mX1Frac;
                nY2Frac = 1 - nY1Frac;
                kZ2Frac = 1 - kZ1Frac;


                qGrid[m2XInt * (yGrid * zGrid) + n2YInt * zGrid + k2ZInt] +=
                    fChargePartDistri[iNumPart] * mX1Frac * nY1Frac * kZ1Frac;
                qGrid[m2XInt * (yGrid * zGrid) + n2YInt * zGrid + k1ZInt] +=
                    fChargePartDistri[iNumPart] * mX1Frac * nY1Frac * kZ2Frac;
                qGrid[m2XInt * (yGrid * zGrid) + n1YInt * zGrid + k2ZInt] +=
                    fChargePartDistri[iNumPart] * mX1Frac * nY2Frac * kZ1Frac;
                qGrid[m2XInt * (yGrid * zGrid) + n1YInt * zGrid + k1ZInt] +=
                    fChargePartDistri[iNumPart] * mX1Frac * nY2Frac * kZ2Frac;
                qGrid[m1XInt * (yGrid * zGrid) + n2YInt * zGrid + k2ZInt] +=
                    fChargePartDistri[iNumPart] * mX2Frac * nY1Frac * kZ1Frac;
                qGrid[m1XInt * (yGrid * zGrid) + n2YInt * zGrid + k1ZInt] +=
                    fChargePartDistri[iNumPart] * mX2Frac * nY1Frac * kZ2Frac;
                qGrid[m1XInt * (yGrid * zGrid) + n1YInt * zGrid + k2ZInt] +=
                    fChargePartDistri[iNumPart] * mX2Frac * nY2Frac * kZ1Frac;
                qGrid[m1XInt * (yGrid * zGrid) + n1YInt * zGrid + k1ZInt] +=
                    fChargePartDistri[iNumPart] * mX2Frac * nY2Frac * kZ2Frac;
            }
        }
    }

    double dArea = xStepDomain * yStepDomain * zStepDomain;
    qGrid = qGrid / dArea;
}

void TBeam::FTDSorPoisson3D(int xGrid, int yGrid, int zGrid, vector<double>& uGrid, double xStep, double yStep , double zStep, double sorFactor, int iterNum)
{
    int i0XGrid, i2XGrid, i0YGrid, i2YGrid, i0ZGrid, i2ZGrid, i1Grid;
    double u0XGrid, u2XGrid, u0YGrid, u2YGrid, u0ZGrid, u2ZGrid, u1Grid;
    double q1Grid;
    double errTol = 1e-6;
    double res, resSum, resAvg;
    double xCoeff = 1 / pow(xStep, 2);
    double yCoeff = 1 / pow(yStep, 2);
    double zCoeff = 1 / pow(zStep, 2);

    do {
        resSum = 0;
        for (int i1XGrid=1; i1XGrid < xGrid - 1; ++i1XGrid)
            for (int i1YGrid=1; i1YGrid < yGrid - 1; ++i1YGrid)
                for (int i1ZGrid=0; i1ZGrid < zGrid; ++i1ZGrid) {
                    i0XGrid = i1XGrid - 1;
                    i0YGrid = i1YGrid - 1;
                    i0ZGrid = (i1ZGrid + zGrid - 2) % (zGrid - 1);
                    i2XGrid = i1XGrid + 1;
                    i2YGrid = i1YGrid + 1;
                    i2ZGrid = (i1ZGrid + zGrid) % (zGrid - 1);

                    u0XGrid = uGrid[i0XGrid*(yGrid*zGrid) + i1YGrid*zGrid + i1ZGrid];

                    u0YGrid = uGrid[i1XGrid * (yGrid * zGrid) + i0YGrid * zGrid + i1ZGrid];

                    u0ZGrid = uGrid[i1XGrid * (yGrid * zGrid) + i1YGrid * zGrid + i0ZGrid];

                    u2XGrid = uGrid[i2XGrid * (yGrid * zGrid) + i1YGrid * zGrid + i1ZGrid];

                    u2YGrid = uGrid[i1XGrid * (yGrid * zGrid) + i2YGrid * zGrid + i1ZGrid];

                    u2ZGrid = uGrid[i1XGrid * (yGrid * zGrid) + i1YGrid * zGrid + i2ZGrid];

                    i1Grid = i1XGrid * (yGrid * zGrid) + i1YGrid * zGrid + i1ZGrid;
                    u1Grid = uGrid[i1Grid];
                    q1Grid = qGrid[i1Grid];

                    res = xCoeff * (d0XGrid + d2XGrid - 2*u1Grid) + yCoeff * (d0YGrid + d2YGrid - 2*u1Grid) + zCoeff * (d0ZGrid + d2ZGrid - 2*u1Grid) - q1Grid;
                    resSum += res;
                    uGrid[i1Grid] -= 2 * res * sorFactor / (xCoeff + yCoeff + zCoeff);
                }
        resAvg = resSum / (xGrid * yGrid * zGrid);
        iterNum -= 1;
    } while (resAvg > errTol && iterNum > 0);
}

void TBeam::Prolongation(vector<double>& uGridCoarse, vector<double>& uGridFine, xGrid, yGrid, zGrid)
{
    int iGrid;
    int xGridCse = xGrid / 2 + 1;
    int yGridCse = yGrid / 2 + 1;
    int zGridCse = zGrid / 2 + 1;
    int iXGridCse, iYGridCse, iZGridCse;

    for (int iXGrid=0; iXGrid != xGrid; ++iXGrid)
        for (int iYGrid=0; iYGrid != yGrid; ++iYGrid)
            for (int iZGrid=0; iZGrid != zGrid; ++iZGrid) {
                iGrid = iXGrid * (yGrid * zGrid) + iYGrid * zGrid + iZGrid;
                iXGridCse = iXGrid / 2;
                iYGridCse = iYGrid / 2;
                iZGridCse = iZGrid / 2;
                i1XGridCse = (iXGrid - 1) / 2;
                i1YGridCse = (iYGrid - 1) / 2;
                i1ZGridCse = (iZGrid - 1) / 2;
                i2XGridCse = (iXGrid + 1) / 2;
                i2YGridCse = (iYGrid + 1) / 2;
                i2ZGridCse = (iZGrid + 1) / 2;

                if (iXGrid % 2 == 0 && iYGrid % 2 == 0 && iZGrid % 2 == 0)
                    uGridFine[iGrid] = uGridCoarse[iXGridCse * (yGridCse*zGridCse) + iYGridCse * zGridCse + iZGridCse];

                else if (iXGrid % 2 == 1 && iYGrid % 2 == 0 && iZGrid % 2 == 0)
                    uGridFine[iGrid] = 0.5 * (uGridCoarse[i1XGridCse * yGridCse * zGridCse + iYGridCse * zGridCse + iZGridCse] + uGridCoarse[i2XGridCse * yGridCse * zGridCse + iYGridCse * zGridCse + iZGridCse]);

                else if (iXGrid % 2 == 0 && iYGrid % 2 == 1 && iZGrid % 2 == 0)
                    uGridFine[iGrid] = 0.5 * (uGridCoarse[iXGridCse*yGridCse*zGridCse + i1YGridCse * zGridCse + iZGridCse] + uGridCoarse[iXGridCse*yGridCse*zGridCse + i2YGridCse * zGridCse + iZGridCse]);

                else if (iXGrid % 2 == 0 && iYGrid % 2 == 0 && iZGrid % 2 == 1)
                    uGridFine[iGrid] = 0.5 * (uGridCoarse[iXGridCse * yGridCse * zGridCse + iYGridCse * zGridCse + i1ZGridCse] + uGridCoarse[iXGridCse*yGridCse*zGridCse + iYGridCse * zGridCse + i2ZGridCse]);

                else if (iXGrid % 2 == 1 && iYGrid % 2 == 1 && iZGrid % 2 == 0)
                    uGridFine[iGrid] = 0.25 * (uGridCoarse[i1XGridCse * yGridCse * zGridCse + i1YGridCseYGrid * zGridCse + iZGridCse] + uGridCoarse[i1XGridCse * yGridCse * zGridCse + i2YGridCseYGrid * zGridCse + iZGridCse] + uGridCoarse[i2XGridCse * yGridCse * zGridCse + i1YGridCse * zGridCse + iZGridCse] + uGridCoarse[i2XGridCse * yGridCse * zGridCse + i2YGridCse * zGridCse + iZGridCse]);

                else if (iXGrid % 2 == 1 && iYGrid % 2 == 0 && iZGrid % 2 == 1)
                    uGridFine[iGrid] = 0.25 * (uGridCoarse[i1XGridCse * yGridCse * zGridCse + iYGridCse  * zGridCse  + i1ZGridCse]  + uGridCoarse[i2XGridCse * yGridCse * zGridCse + iYGridCse  * zGridCse  + i1ZGridCse] + uGridCoarse[i1XGridCse * yGridCse * zGridCse + iYGridCse  * zGridCse  + i2ZGridCse] + uGridCoarse[i2XGridCse * yGridCse * zGridCse + iYGridCse  * zGridCse  + i2ZGridCse]);

                else if (iXGrid % 2 == 0 && iYGrid % 2 == 1 && iZGrid % 2 == 1)
                    uGridFine[iGrid] = 0.25 * (uGridCoarse[iXGridCse * yGridCse * zGridCse + i1YGridCse * zGridCse + i1ZGridCse] + uGridCoarse[iXGridCse * yGridCse * zGridCse + i1YGridCse * zGridCse + i2ZGridCse] + uGridCoarse[iXGridCse * yGridCse * zGridCse + i2YGridCse * zGridCse + i1ZGridCse] + uGridCoarse[iXGridCse * yGridCse * zGridCse + i2YGridCse * zGridCse + i2ZGridCse]);

                else 
                    uGridFine[iGrid] = 0.125 * (uGridCoarse[i1XGridCse * yGridCse * zGridCse + i1YGridCse * zGridCse + i1ZGridCse] + uGridCoarse[i1XGridCse * yGridCse * zGridCse + i1YGridCse * zGridCse + i2ZGridCse] + uGridCoarse[i1XGridCse * yGridCse * zGridCse + i2YGridCse * zGridCse + i1ZGridCse] + uGridCoarse[i1XGridCse * yGridCse * zGridCse + i2YGridCse * zGridCse + i2ZGridCse] + uGridCoarse[i2XGridCse * yGridCse * zGridCse + i1YGridCse * zGridCse + i1ZGridCse] + uGridCoarse[i2XGridCse * yGridCse * zGridCse + i1YGridCse * zGridCse + i2ZGridCse] + uGridCoarse[i2XGridCse * yGridCse * zGridCse + i2YGridCse * zGridCse + i1ZGridCse] + uGridCoarse[i2XGridCse * yGridCse * zGridCse + i2YGridCse * zGridCse + i2ZGridCse]); 
            }
}

void TBeam::CalcResGrid(vector<double>& uGrid, vector<double>& resGrid, int xGrid, int yGrid, int zGrid, int level)
{
    double xCoeff = 1 / pow(xStep, 2);
    double yCoeff = 1 / pow(yStep, 2);
    double zCoeff = 1 / pow(zStep, 2);
    double u0XGrid, u2XGrid, u0YGrid, u2YGrid, u0ZGrid, u2ZGrid, u1Grid;
    int i0XGrid, i2XGrid, i0YGrid, i2YGrid, i0ZGrid, i2ZGrid, i1Grid;

    for (int i1XGrid=1; i1XGrid != xGrid - ; ++i1XGrid)
        for (int i1YGrid=1; i1YGrid != yGrid - 1; ++i1YGrid)
            for (int i1ZGrid=0; i1ZGrid != zGrid; ++i1ZGrid) {
                //resGrid[iXGrid * yGrid * zGrid + iYGrid * zGrid + iZGrid] = (qGrid[pow(level, 3) * iXGrid * yGrid * zGrid + pow(level, 2) * iYGrid * zGrid + level * iZGrid] - uGrid[iXGrid * yGrid * zGrid + iYGrid * zGrid + iZGrid]);  
                i0XGrid = i1XGrid - 1;
                i0YGrid = i1YGrid - 1;
                i0ZGrid = (i1ZGrid + zGrid - 2) % (zGrid - 1);
                i2XGrid = i1XGrid + 1;
                i2YGrid = i1YGrid + 1;
                i2ZGrid = (i1ZGrid + zGrid) % (zGrid - 1);

                u0XGrid = uGrid[i0XGrid*(yGrid*zGrid) + i1YGrid*zGrid + i1ZGrid];

                u0YGrid = uGrid[i1XGrid * (yGrid * zGrid) + i0YGrid * zGrid + i1ZGrid];

                u0ZGrid = uGrid[i1XGrid * (yGrid * zGrid) + i1YGrid * zGrid + i0ZGrid];

                u2XGrid = uGrid[i2XGrid * (yGrid * zGrid) + i1YGrid * zGrid + i1ZGrid];

                u2YGrid = uGrid[i1XGrid * (yGrid * zGrid) + i2YGrid * zGrid + i1ZGrid];

                u2ZGrid = uGrid[i1XGrid * (yGrid * zGrid) + i1YGrid * zGrid + i2ZGrid];

                i1Grid = i1XGrid * (yGrid * zGrid) + i1YGrid * zGrid + i1ZGrid;
                u1Grid = uGrid[i1Grid];
                //q1Grid = qGrid[i1Grid];
                q1Grid = qGrid[pow(level, 3) * iXGrid * yGrid * zGrid + pow(level, 2) * iYGrid * zGrid + level * iZGrid];  

                resGrid[iXGrid * yGrid * zGrid + iYGrid * zGrid + iZGrid] = xCoeff * (d0XGrid + d2XGrid - 2*u1Grid) + yCoeff * (d0YGrid + d2YGrid - 2*u1Grid) + zCoeff * (d0ZGrid + d2ZGrid - 2*u1Grid) - q1Grid;
}
