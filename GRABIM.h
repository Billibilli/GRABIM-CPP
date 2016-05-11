#ifndef GRABIM_H
#define GRABIM_H

#include <nlopt.hpp>
#include <armadillo>
#include "sparengine.h"
#include <queue>
#include <QString>
#include <QStringList>

using namespace arma;
using namespace std;

enum ObjectiveFunction {NINF_S11dB, NINF_POWERTRANS};

typedef struct GRABIM_Result {
    rowvec x_grid_search;
    double grid_val;
    rowvec x_nlopt;
    double nlopt_val;
    cx_vec ZS, ZL;
    vec freq;
    cx_vec S11_gridsearch, S21_gridsearch, S12_gridsearch, S22_gridsearch;
    cx_vec S11_nlopt, S21_nlopt, S12_nlopt, S22_nlopt;
    QString topology;
    QString source_path;
    QString load_path;
    QString QucsVersion;

} GRABIM_Result;


//Reference:
// [1] Broadband direct-coupled and RF matching networks. Thomas R. Cuthbert, 1999
// [2] Broadband impedance matching - Fast and simple. Thomas R. Cuthbert. RF design. 1994
class GRABIM
{
    const double c0 = 299792458;//Speed of light (m/s)

public:
    GRABIM();

    int SetSourceImpedance(cx_vec);
    int SetLoadImpedance(cx_vec);
    int SetFrequency(vec);

    int SetInitialPivot(rowvec);
    rowvec GetInitialPivot();
    int SetMatchingBand(double, double);
    int SetTopology(std::string);
    GRABIM_Result RunGRABIM();
    int SetMaxIterGridSearch(int);
    int SetThreshold(double);
    double GetThreshold();
    int SetVerbose(bool);
    double CandidateEval(rowvec);
    int SetNLoptAlg(nlopt::algorithm);
    nlopt::algorithm GetNLoptAlg();
    int SetObjectiveFunction(ObjectiveFunction);
    ObjectiveFunction GetObjectiveFunction();

private:
    rowvec GridSearch();
    rowvec LocalOptimiser(rowvec);
    int ResampleImpedances();
    mat GeneratingMatrix(int);
    double CalcInvPowerTransfer(cx_mat, cx_double, cx_double);
    rowvec InspectCandidate(rowvec);
    rowvec x_ini;

    vec freq;
    cx_vec ZS, ZL;

    bool verbose;
    std::string topology;
    unsigned int Grid_MaxIter;
    double MatchingThreshold;
    nlopt::algorithm NLoptAlgo;
    ObjectiveFunction ObjFun;

    string tolower(string str);
    string RemoveBlankSpaces(string line);

    int SearchPredefinedTopologies(rowvec &, std::string &);
    void AutoSetInitialPivot();


};

#endif // GRABIM_H