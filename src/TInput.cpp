#include "TInput.h"
using namespace std;

TInput::TInput() {}

TInput::~TInput() {}

void TInput::Input(string filename = "input") {
  string oString;
  string oStringBefore;
  string oStringAfter;
  string oOptJob = "";
  string::size_type oPosition;
  int readLatticeFlag = 0;

  //***************  add AND DEFINE
  int beamRenameFlag = 0;
  fMultiBeam = 0;

  //***************

  fstream oFile(filename, ios::in);
  if (oFile.fail()) {
    cout << "Error opening the input file!" << endl;
    exit(1);
  }

  while (!oFile.eof()) {
    getline(oFile, oString);
    transform(oString.begin(), oString.end(), oString.begin(), ::tolower);

    oString = Ptrim(oString);
    oString = Trim(oString);

    if (oString == "")
      continue;
    oPosition = oString.find("=");
    if (oPosition != string::npos) {
      oStringBefore = oString.substr(0, oPosition);
      oStringAfter = oString.substr(oPosition + 1);
      oStringBefore = Rtrim(oStringBefore);
      oStringBefore = oOptJob + oStringBefore;
      oStringAfter = Ltrim(oStringAfter);
      if (beamRenameFlag == 1 && isalpha(oStringBefore.back())) {
        oStringBefore = oStringBefore + "_1";
      }
      if (beamRenameFlag == 1 && oStringBefore.substr(0, 5) == "beam_") {
        fMultiBeam += 1;
      }

    } else if (readLatticeFlag == 0 and oString.find(":") != string::npos) {
      oOptJob = oString.substr(0, oString.find(":")) + "_";
      continue;
    } else if (oString.substr(0, 5) == "[latt") {
      readLatticeFlag = 1;

      //    oPeriod.push_back(new TPeriod());
      //    (*oPeriod.back()).fNamePeriod="lattice";

      //~~~~~~~~~~~~~~~~
    } else if (oString.substr(0, 5) == "[beam") {
      beamRenameFlag = 1;
      continue;
      //~~~~~~~~~~~~~~~~~
    } else if (oString[0] == '[') {
      oOptJob = "";
      readLatticeFlag = 0;
      beamRenameFlag = 0;
      continue;
    } else if (readLatticeFlag == 1) {
      //********************************************************   sentence....

      //*********************************************************
    } else
      continue;
    if (!oStringAfter.find_first_of("-0123456789")) {
      double oStringAfterDbl = atof(oStringAfter.c_str());

      fAcceDbl[oStringBefore] = oStringAfterDbl;
    } else {
      fAcceStr[oStringBefore] = oStringAfter;
    }
  }
  // * /  * /  * /  * /  * /  * /  * /  * /  * /  * /  * /  * /  * /  * /

  /*
  map<string,string>::iterator iMapStrStr;
  for(iMapStrStr=fAcceStr.begin();iMapStrStr!=fAcceStr.end();++iMapStrStr)
      cout<<iMapStrStr->first<<"  |   "<<iMapStrStr->second<<endl;

  cout<<"+++++++++++++++++++++++++++++"<<endl;

   map<string,double>::iterator iMapStrDbl;
   for(iMapStrDbl=fAcceDbl.begin();iMapStrDbl!=fAcceDbl.end();++iMapStrDbl)
     cout<<iMapStrDbl->first<<"  |   "<<iMapStrDbl->second<<endl;
  */
  /*
  map<string,vector<int> >::iterator iMapVectorInt;
  for(iMapVectorInt=fPeriodRecord.begin();iMapVectorInt!=fPeriodRecord.end();++iMapVectorInt){
   cout<<iMapVectorInt->first<<"  |   "<<iMapVectorInt->second[0]<<"  |
  "<<iMapVectorInt->second[1]<<endl;
  }
  */

  // * /  * /  * /  * /  * /  * /  * /  * /  * /  * /  * /  * /  * /  * /  * /
}
//###################################################
string TInput::Ltrim(string &str) {
  str.erase(str.begin(),
            find_if(str.begin(), str.end(), not1(ptr_fun(::isspace))));
  return str;
}

string TInput::Rtrim(string &str) {
  str.erase(find_if(str.rbegin(), str.rend(), not1(ptr_fun(::isspace))).base(),
            str.end());
  return str;
}

string TInput::Trim(string &str) {
  str = Ltrim(str);
  str = Rtrim(str);
  return str;
}

string TInput::Ptrim(string &str) {
  int posi = str.find("#");
  str = str.substr(0, posi);
  return str;
}
