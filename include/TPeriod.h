#ifndef TPERIOD_H
#define TPERIOD_H
#include "TBase.h"
#include "TElement.h"

class TPeriod {
public:
  TPeriod();
  ~TPeriod();

  string fNamePeriod; // period's name

  int fNumElement; // The number of the elements in this period.
  int fCellRepeat; // Record the cell repeat.

  vector<string> fLabelPeriod;       // the element's name in this period
  vector<string> fLabelGlobal;       // the element's name in the whole lattice
  vector<TElement *> fLatticePeriod; // the element's pointer in this period.
  vector<string> fInputPreName;      // The PreName

  void Update();

protected:
private:
  string fPreName;
};

#endif
