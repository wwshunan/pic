#include "TBeam.h"

using namespace std;

// vector<double> operator+(vector<double> &aVec, vector<double> &bVec);

int main() {
  TApplication App("app", NULL, NULL);

  cout << "---------------------------------" << endl;

  TBeam oBeam;
  oBeam.Input("../input");
  oBeam.PrintMap<map<string, string>>(oBeam.fAcceStr);
  oBeam.PrintMap<map<string, double>>(oBeam.fAcceDbl);

  oBeam.PartGen();
  oBeam.MeshFFT2D(-1, 1, -1, 1, 7, 9);
  oBeam.FFTPoisson2D();
  oBeam.Gradient2D();

  // oBeam.Plot3D(oBeam.xEGrid, oBeam.yGrid);

  oBeam.Plot3D(oBeam.fXPartDistri, oBeam.fYPartDistri, oBeam.fZPartDistri,
               "Particle3D_2", "Particle3D_2", "PA*");

  double tB = 3;
  oBeam.MeshFFT3D(-tB, tB, -tB, tB, -tB, tB, 7, 6, 5);

  oBeam.Plot4D(oBeam.xGrid, oBeam.yGrid, oBeam.zGrid, oBeam.qGrid,
               "../qGrid.txt");

  oBeam.FFTPoisson3D();
  oBeam.Plot4D(oBeam.xGrid, oBeam.yGrid, oBeam.zGrid, oBeam.uGrid,
               "../uGrid.txt");

  oBeam.Laplace3D();
  oBeam.Gradient3D();

  // oBeam.oCanvas2D["Particle3D"]->WaitPrimitive();

  // MeshFFT3D();

  // oBeam.Plot1D(oBeam.fXPartDistri, oBeam.fYPartDistri);

  // *****************
  cout << "++++++++++++++++++++++++++++++++" << endl;
  App.Run();

  cout << "================================" << endl;
  return 0;
}
