#ifndef SPARENGINE_H
#define SPARENGINE_H

#include <armadillo>
#include <string.h>
#include <set>
using namespace arma;
using namespace std;



class SparEngine
{
    const double c0 = 299792458;//Speed of light (m/s)
    std::set<string> preTopoList;
public:
    SparEngine();
    cx_mat getSparams(rowvec, cx_double, cx_double, double, std::string);
    cx_mat getABCDmatrix(rowvec, double, std::string);
    cx_mat PreComputedABCD(rowvec, double, std::string);
};

#endif // SPARENGINE_H
