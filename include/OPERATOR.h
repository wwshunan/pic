#ifndef OPERATOR_H
#define OPERATOR_H
#include <algorithm> // std::transform
#include <cmath>
#include <functional> // std::plus
#include <vector>

using std::vector;
vector<double> operator+(vector<double> &, vector<double> &);
vector<double> operator-(vector<double> &, vector<double> &);
vector<double> operator*(vector<double> &, vector<double> &);
vector<double> operator/(vector<double> &, vector<double> &);

vector<double> operator+(vector<double> &, double &);
vector<double> operator-(vector<double> &, double &);
vector<double> operator*(vector<double> &, double &);
vector<double> operator/(vector<double> &, double &);

vector<double> operator+(double &, vector<double> &);
vector<double> operator-(double &, vector<double> &);
vector<double> operator*(double &, vector<double> &);
vector<double> operator/(double &, vector<double> &);

#endif
