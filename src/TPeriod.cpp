#include "TPeriod.h"

using namespace std;

TPeriod::TPeriod(){

}
TPeriod::~TPeriod(){
}



void TPeriod::Update(){
  fNumElement=fLabelPeriod.size();

  fPreName="";
  for(auto &iNumCell:fInputPreName){
    cout<<iNumCell<<"***"<<endl;
  //  fPreName=fPreName+(*iNumCell);
  }


  for(int iNumElement=0;iNumElement<fNumElement;++fNumElement){
//    fLabelGlobal[iNumElement]=fInputPreName+fLabelPeriod[iNumElement];
  }
}
