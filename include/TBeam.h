#ifndef TBEAM_H
#define TBEAM_H

#include "TInput.h"

using namespace std;

class TBeam : public TInput {
public:
  TBeam();
  virtual ~TBeam();

  vector<double> fXPartDistri;
  vector<double> fXpPartDistri;
  vector<double> fYPartDistri;
  vector<double> fYpPartDistri;
  vector<double> fZPartDistri;
  vector<double> fZpPartDistri;
  vector<double> fFlagPartDistri;
  vector<double> fChargePartDistri;
  vector<double> fMassPartDistri;

  void PartGen();

  void GS6D(int, vector<double> &, vector<double> &, vector<double> &,
            vector<double> &, vector<double> &, vector<double> &);
  void WB6D(int, vector<double> &, vector<double> &, vector<double> &,
            vector<double> &, vector<double> &, vector<double> &);

  void GS4D(int, vector<double> &, vector<double> &, vector<double> &,
            vector<double> &);
  void WB4D(int, vector<double> &, vector<double> &, vector<double> &,
            vector<double> &);
  void KV4D(int, vector<double> &, vector<double> &, vector<double> &,
            vector<double> &);
  void GS2D(int, vector<double> &, vector<double> &);
  void KV2D(int, vector<double> &, vector<double> &);

  void TTrans(vector<double> &, vector<double> &, double, double, double);

  void UN3D(int, vector<double> &, vector<double> &, vector<double> &);

  //
  TInput oInputTest;

  void MeshFFT2D(double xMinDomain = -10, double xMaxDomain = 10,
                 double yMinDomain = -10, double yMaxDomain = 10,
                 int xGridLog = 5, int yGridLog = 5);
  void FFTPoisson2D();
  void Gradient2D();
  void InterPoint2D();

  void FTDPoisson2D();
  void FTDSorPoisson2D();
  void FieldSpaceCharge2D();

  void MeshFFT3D(double xMinDomain = -10, double xMaxDomain = 10,
                 double yMinDomain = -10, double yMaxDomain = 10,
                 double zMinDomain = -10, double zMaxDomain = 10,
                 int xGridLog = 5, int yGridLog = 5, int zGridLog = 5);

  void FFTPoisson3D();
  void Gradient3D();

  void FTDPoisson3D();
  void FTDSorPoisson3D();
  void FieldSpaceCharge3D();

  void Laplace2D();
  void Laplace3D();

  //  vector<vector<double>> qGrid;
  vector<double> qGrid;
  vector<double> uGrid;
  vector<double> xEGrid;
  vector<double> yEGrid;
  vector<double> zEGrid;

  int xGrid, yGrid, zGrid;

protected:
  int fNumPart;

  double xStepDomain, yStepDomain, zStepDomain;
};

#endif // TBEAM_H
