#include "OPERATOR.h"
using namespace std;
vector<double> operator+(vector<double> &aVec, vector<double> &bVec) {
  vector<double> cVec;
  cVec.resize(aVec.size());
  transform(aVec.begin(), aVec.end(), bVec.begin(), cVec.begin(),
            plus<double>());
  return cVec;
}
vector<double> operator-(vector<double> &aVec, vector<double> &bVec) {
  vector<double> cVec;
  cVec.resize(aVec.size());
  transform(aVec.begin(), aVec.end(), bVec.begin(), cVec.begin(),
            minus<double>());
  return cVec;
}
vector<double> operator*(vector<double> &aVec, vector<double> &bVec) {
  vector<double> cVec;
  cVec.resize(aVec.size());
  transform(aVec.begin(), aVec.end(), bVec.begin(), cVec.begin(),
            multiplies<double>());
  return cVec;
}
vector<double> operator/(vector<double> &aVec, vector<double> &bVec) {
  vector<double> cVec;
  cVec.resize(aVec.size());
  transform(aVec.begin(), aVec.end(), bVec.begin(), cVec.begin(),
            divides<double>());
  return cVec;
}

vector<double> operator+(vector<double> &aVec, double &bDbl) {
  vector<double> cVec;
  cVec.resize(aVec.size());
  transform(aVec.begin(), aVec.end(), cVec.begin(),
            bind2nd(plus<double>(), bDbl));
  return cVec;
}
vector<double> operator-(vector<double> &aVec, double &bDbl) {
  vector<double> cVec;
  cVec.resize(aVec.size());
  transform(aVec.begin(), aVec.end(), cVec.begin(),
            bind2nd(minus<double>(), bDbl));
  return cVec;
}
vector<double> operator*(vector<double> &aVec, double &bDbl) {
  vector<double> cVec;
  cVec.resize(aVec.size());
  transform(aVec.begin(), aVec.end(), cVec.begin(),
            bind2nd(multiplies<double>(), bDbl));
  return cVec;
}
vector<double> operator/(vector<double> &aVec, double &bDbl) {
  vector<double> cVec;
  cVec.resize(aVec.size());
  transform(aVec.begin(), aVec.end(), cVec.begin(),
            bind2nd(divides<double>(), bDbl));
  return cVec;
}

vector<double> operator+(double &aDbl, vector<double> &bVec) {
  vector<double> cVec;
  cVec.resize(bVec.size());
  transform(bVec.begin(), bVec.end(), cVec.begin(),
            bind1st(plus<double>(), aDbl));
  return cVec;
}
vector<double> operator-(double &aDbl, vector<double> &bVec) {
  vector<double> cVec;
  cVec.resize(bVec.size());
  transform(bVec.begin(), bVec.end(), cVec.begin(),
            bind1st(minus<double>(), aDbl));
  return cVec;
}
vector<double> operator*(double &aDbl, vector<double> &bVec) {
  vector<double> cVec;
  cVec.resize(bVec.size());
  transform(bVec.begin(), bVec.end(), cVec.begin(),
            bind1st(multiplies<double>(), aDbl));
  return cVec;
}
vector<double> operator/(double &aDbl, vector<double> &bVec) {
  vector<double> cVec;
  cVec.resize(bVec.size());
  transform(bVec.begin(), bVec.end(), cVec.begin(),
            bind1st(divides<double>(), aDbl));
  return cVec;
}
