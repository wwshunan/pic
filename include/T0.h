#ifndef T0_H
#define T0_H
// stdC++ library
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include <fftw3.h>
#include <math.h>
#include <vector>

// DIY library
#include "TVector.h"

// root library
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TRandom3.h"

class T0 {
protected:
  TRandom3 oRandom3;
};
#endif
