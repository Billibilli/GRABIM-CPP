#ifndef IO_H
#define IO_H
#include "GRABIM.h"
#include <locale>

using namespace std;

enum terminal {SOURCE, LOAD};

class IO
{
public:
    IO();
    int exportGNUplot(GRABIM_Result, string);
    int loadS1Pdata(std::string, terminal);
    int ResampleImpedances();
    cx_mat getSourceImpedance();
    cx_mat getLoadImpedance();
    vec getFrequency();
    void set_constant_ZS_vs_freq(cx_double);
    void set_constant_ZL_vs_freq(cx_double);
    void set_matching_band(double, double);
    void setLocalOptimiser(nlopt::algorithm);
    nlopt::algorithm getLocalOptimiser();
    int ExportQucsSchematic(GRABIM_Result);
    void set_qucs_sch_path(string);
    string get_qucs_sch_path();
    void PrintNetwork_StandardOutput(GRABIM_Result);


private:
    // ZS and ZL are the source and load impedances, respectively whereas fS and fL indicates the frequencies where
    // ZS and ZL were sampled
    cx_vec ZS, ZL;
    vec fS, fL;

    vec freq;//More often than not, ZS and ZL are sampled at different frecuencies, so it is necessary to have
    // common frequency vector for pairing ZS and ZL.


    double fmatching_min, fmatching_max;
    int getFreqIndex(double);

    cx_vec ZS_matching, ZL_matching;
    vec f_matching;
    int setMatchingImpedances();
    double getS2PfreqScale(string line);
    nlopt::algorithm LocalOptAlgo;

    int SchematicParser(GRABIM_Result, int &, string &, string &, string &);
    bool CreateSchematic(string, string, string, string);

    int Nsamples;//Impedance samples within matching band


    string Num2String(int x);
    string Num2String(double x);

    string QucsSchematicPath;
};

#endif // IO_H
