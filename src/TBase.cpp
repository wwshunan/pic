#include "TBase.h"

TBase::TBase() {}
TBase::~TBase() {}

using namespace std;

void TBase::Plot1D(vector<double> oY, string oName, string oTitle,
                   string oDraw) {
  oCanvas1D.insert(pair<string, TCanvas *>(
      oName, new TCanvas(oName.c_str(), oTitle.c_str(), 800, 600)));
  oGraph1D.insert(pair<string, TGraph *>(oName, new TGraph()));

  TGraph *oGraph = oGraph1D[oName];
  int nY = oY.size();
  for (int iY = 0; iY < nY; ++iY) {
    int iY1 = iY + 1;
    oGraph->SetPoint(iY1, iY1, oY[iY]);
  }
  oGraph1D[oName]->Draw(oDraw.c_str());
  oCanvas1D[oName]->Draw();
}

void TBase::Plot1D(vector<double> oX, vector<double> oY, string oName,
                   string oTitle, string oDraw) {
  oCanvas1D.insert(pair<string, TCanvas *>(
      oName, new TCanvas(oName.c_str(), oTitle.c_str(), 800, 600)));
  oGraph1D.insert(pair<string, TGraph *>(oName, new TGraph()));

  TGraph *oGraph = oGraph1D[oName];
  int nY = oY.size();
  for (int iY = 0; iY < nY; ++iY) {
    oGraph->SetPoint(iY, oX[iY], oY[iY]);
  }
  oGraph1D[oName]->Draw(oDraw.c_str());
  oCanvas1D[oName]->Draw();
}

void TBase::Plot3D(vector<double> oX, vector<double> oY, vector<double> oZ,
                   string oName, string oTitle, string oDraw) {
  oCanvas2D.insert(pair<string, TCanvas *>(
      oName, new TCanvas(oName.c_str(), oTitle.c_str(), 800, 600)));
  oGraph2D.insert(pair<string, TGraph2D *>(oName, new TGraph2D()));
  TGraph2D *oGraph = oGraph2D[oName];

  int nY = oY.size();
  for (int iY = 0; iY < nY; ++iY) {
    oGraph->SetPoint(iY, oX[iY], oY[iY], oZ[iY]);
  }
  oGraph2D[oName]->Draw(oDraw.c_str());
  oCanvas2D[oName]->Draw();
}

void TBase::Plot3D(vector<double> oX, int dimColume, string oName,
                   string oTitle, string oDraw) {
  oCanvas2D.insert(pair<string, TCanvas *>(
      oName, new TCanvas(oName.c_str(), oTitle.c_str(), 800, 600)));
  oGraph2D.insert(pair<string, TGraph2D *>(oName, new TGraph2D()));
  TGraph2D *oGraph = oGraph2D[oName];

  int nX = oX.size();
  int oM, oN;
  for (int iX = 0; iX < nX; ++iX) {
    oM = iX / dimColume;
    oN = iX % dimColume;
    oGraph->SetPoint(iX, oM, oN, oX[iX]);
  }
  oGraph2D[oName]->Draw(oDraw.c_str());
  oCanvas2D[oName]->Draw();
}

void TBase::Plot3D(vector<vector<double>> &oX, string oName, string oTitle,
                   string oDraw) {
  oCanvas2D.insert(pair<string, TCanvas *>(
      oName, new TCanvas(oName.c_str(), oTitle.c_str(), 800, 600)));
  oGraph2D.insert(pair<string, TGraph2D *>(oName, new TGraph2D()));
  TGraph2D *oGraph = oGraph2D[oName];

  int iCount = 0;
  int nX = oX.size();
  int nY = oX[0].size();
  for (int iX = 0; iX != nX; ++iX)
    for (int iY = 0; iY != nY; ++iY) {

      oGraph->SetPoint(iCount, iX, iY, oX[iX][iY]);
      ++iCount;
    }

  oGraph2D[oName]->Draw(oDraw.c_str());
  oCanvas2D[oName]->Draw();
}

void TBase::Plot3D(double *oX, int numX, int dimColume, string oName,
                   string oTitle, string oDraw) {

  oCanvas2D.insert(pair<string, TCanvas *>(
      oName, new TCanvas(oName.c_str(), oTitle.c_str(), 800, 600)));
  oGraph2D.insert(pair<string, TGraph2D *>(oName, new TGraph2D()));
  TGraph2D *oGraph = oGraph2D[oName];
  int nY, nX = 0;
  for (int iX = 0; iX < numX; ++iX) {
    nX = iX / dimColume;
    nY = iX % dimColume;
    oGraph->SetPoint(iX, nX, nY, oX[iX]);
  }
  oGraph2D[oName]->Draw(oDraw.c_str());
  oCanvas2D[oName]->Draw();
}

// void TBase::Plot4D(int &nX, int &nY, int &nZ, vector<double> &oData,
//                    string oName, string oTitle, string oDraw) {
//   oCanvas2D.insert(pair<string, TCanvas *>(
//       oName, new TCanvas(oName.c_str(), oTitle.c_str(), 800, 600)));
//   oGraph2D.insert(pair<string, TGraph2D *>(oName, new TGraph2D()));
//   TGraph2D *oGraph = oGraph2D[oName];
//   int nData = 0;
//   int iData;
//
//   double maxData = *(max_element(oData.begin(), oData.end()));
//
//   cout << "@@@@@@@@@@@@@@ " << maxData << endl;
//
//   for (int iX = 0; iX < nX; ++iX) {
//     for (int iY = 0; iY < nY; ++iY) {
//       for (int iZ = 0; iZ < nZ; ++iZ) {
//         iData = oData[iX * (nY * nZ) + iY * nZ + iZ];
//         oGraph->SetPoint(nData, iX, iY, iZ);
//         ++nData;
//       }
//     }
//   }
//   oGraph2D[oName]->Draw(oDraw.c_str());
//   oCanvas2D[oName]->Draw();
//
//   //
//
//   // double *xMarker = oGraph2D[oName]->GetX();
//   // double *yMarker = oGraph2D[oName]->GetY();
//   // double *zMarker = oGraph2D[oName]->GetZ();
//   //
//   // for (int iMarker = 0; iMarker < nX * nY * nZ; iMarker++) {
//   //   iData = oData[iMarker];
//   //   TMarker *marker4D = new TMarker(xMarker[iMarker], yMarker[iMarker],
//   22);
//   //   marker4D->SetMarkerStyle(20);
//   //   // marker4D->SetMarkerSize(1);
//   //   marker4D->SetMarkerColor(30 + (int)(iData / maxData) * 70);
//   //   marker4D->Draw();
//   // }
//
//   nData = 0;
//   for (int iX = 0; iX < nX; ++iX) {
//     for (int iY = 0; iY < nY; ++iY) {
//       for (int iZ = 0; iZ < nZ; ++iZ) {
//         iData = oData[nData];
//         ++nData;
//
//         TMarker *marker4D = new TMarker(iX, iY, 1);
//         marker4D->SetMarkerStyle(20);
//         marker4D->SetMarkerSize(1);
//         marker4D->SetMarkerColor(30 + (int)(iData / maxData) * 70);
//         marker4D->Draw();
//       }
//     }
//   }
// }

void TBase::Plot4D(int &nX, int &nY, int &nZ, vector<double> &oData,
                   string fileName) {
  ofstream filePlot4D(fileName.c_str());

  double iData;

  if (filePlot4D.is_open()) {
    for (int iX = 0; iX < nX; ++iX) {
      for (int iY = 0; iY < nY; ++iY) {
        for (int iZ = 0; iZ < nZ; ++iZ) {
          iData = oData[iX * (nY * nZ) + iY * nZ + iZ];
          filePlot4D << iX << " " << iY << " " << iZ << " " << iData << "\n";
          // cout << iX << " " << iY << " " << iZ << " " << iData << "\n";
        }
      }
    }
  }
  filePlot4D.close();
}
