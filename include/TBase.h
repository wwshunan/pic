#ifndef TBASE_H
#define TBASE_H
#define PI 3.1415926535897932

#include "T0.h"

using std::vector;
using std::map;
using std::string;

class TBase : public T0 {
public:
  TBase();
  ~TBase();

  template <typename T> void PrintMap(T fMapPrint) {
    for (const auto &iTest : fMapPrint) {
      std::cout << iTest.first << "  |  " << iTest.second << std::endl;
    }
  }
  template <typename T> void PrintVec(T fVecPrint) {
    for (const auto &iTest : fVecPrint) {
      std::cout << iTest << std::endl;
    }
  }

  //
  // void Plot1D(vector<double>);
  // void Plot1D(int, vector<double>);
  // void Plot1D(int, double*);
  // void Plot1D(vector<double>,vector<double>);
  // void Plot1D(int,vector<double>,vector<double>);
  // void Plot1D(int,double*,double*);

  void Plot1D(vector<double>, string oName = "name", string oTitle = "title",
              string oDraw = "AP");
  // void Plot1D(int, vector<double>);
  // void Plot1D(int, double *);
  void Plot1D(vector<double>, vector<double>, string oName = "name",
              string oTitle = "title", string oDraw = "AP");
  // void Plot1D(int, vector<double>, vector<double>);
  // void Plot1D(int, double *, double *);

  map<string, TCanvas *> oCanvas1D;
  map<string, TGraph *> oGraph1D;

  map<string, TCanvas *> oCanvas2D;
  map<string, TGraph2D *> oGraph2D;
  void Plot3D(vector<double>, vector<double>, vector<double>,
              string oName = "name", string oTitle = "title",
              string oDraw = "AP");

  void Plot3D(vector<double>, int dimColume, string oName = "name",
              string oTitle = "title", string oDraw = "AP");

  void Plot3D(vector<vector<double>> &, string oName = "vector2D",
              string oTitle = "vector2D", string oDraw = "AP");
  void Plot3D(double *, int, int dimColume, string oName = "DoubleArray",
              string oTitle = "DoubleArray", string oDraw = "AP");
  // void Plot4D(int &, int &, int &, vector<double> &, string oName = "Plot4D",
  //             string oTitle = "Plot4D", string oDraw = "AP");

  void Plot4D(int &nX, int &nY, int &nZ, vector<double> &oData,
              string fileName = "Plot4D.txt");

protected:
private:
};

#endif // TBASE_H
