#include <iostream>
#include "matchingnetwork.h"

using namespace std;

int main()
{
    MatchingNetwork MatchingObject;

    ///// MATCHING SETTINGS ////////
    double f1 = 2e9;
    double f2 = 4e9;
    int N = 20;//Number of samples to match in the matching band

    std::string sourcefile = "./s2p_files/f531p84g0.s2p";
    std::string loadfile = "./s2p_files/TQP7M9101.s2p";

    rowvec x_ini, x_opt;


    MatchingObject.SetMatchingBand(f1, f2, N);

    // SET SOURCE AND LOAD IMPEDANCES
    //MatchingObject.SetSourceImpedance(sourcefile);
    //MatchingObject.SetLoadImpedance(loadfile);
    vec freq = linspace(5e8, 10e9, 200);

    vec ZSr = 50*ones(freq.n_rows,1);
    vec ZSi = zeros(freq.n_rows,1);
    vec ZLr = 100*ones(freq.n_rows,1);
    vec ZLi = zeros(freq.n_rows,1);

    cx_vec ZS(ZSr, ZSi);
    cx_vec ZL(ZLr, ZLi);

    MatchingObject.SetSourceImpedance(ZS, freq);
    MatchingObject.SetLoadImpedance(ZL, freq);

    MatchingObject.SetTopology("444");

    x_ini << 50 << 0.1 << 75 << 0.1 << 100 << 0.1;
    MatchingObject.SetInitialPivot(x_ini);
    MatchingObject.SetMaxIterGridSearch(1000);
    MatchingObject.SetThreshold(-30);
    GRABIM_Result R = MatchingObject.RunGRABIM();
    std::cout << "GRID SEARCH: S11_max = "<< R.grid_val << "dB <= " << R.x_grid_search << std::endl;
    std::cout << "NLOPT: S11_max = "<< R.nlopt_val << "dB <= " <<  R.x_nlopt << std::endl;

}


