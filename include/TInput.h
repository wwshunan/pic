#ifndef TINPUT_H
#define TINPUT_H

#include "OPERATOR.h"
#include "TBase.h"
#include "TElement.h"
#include "TPeriod.h"

class TInput : public TBase {
public:
  TInput();
  ~TInput();

  void Input(string filename);
  void SetInput(string, double);

  string Ptrim(string &str);
  string Ltrim(string &str);
  string Rtrim(string &str);
  string Trim(string &str);

  // protected:

  map<string, double> fAcceDbl;
  map<string, string> fAcceStr;

  vector<string> fLatticeLabel;
  vector<TElement *> fLattice;

  int fMultiBeam;

private:
};

#endif
