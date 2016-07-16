#ifndef TELEMENT_H
#define TELEMENT_H
#include "TBase.h"

using namespace std;
class TElement:public TBase{
public:
  TElement(){};
  ~TElement(){};
protected:

};

//
class TDrift:public TElement {
protected:
        float l,r,ry;
public:
        TDrift(string str2) ;
        ~TDrift(){};
};

class TField_Map:public TElement {
public:

        TField_Map(string str2) ;
        ~TField_Map(){};
protected:
        int geom;
        float l,pha,r,kb,ke,ki,ka;
        string filename;

};
//

#endif
