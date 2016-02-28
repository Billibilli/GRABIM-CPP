#include "matchingnetwork.h"

MatchingNetwork::MatchingNetwork()
{
    verbose  = false;
}

int MatchingNetwork::SetVerbose(bool vb)
{
    verbose = vb;
}

int MatchingNetwork::SetSourceImpedance(std::string sourcepath)
{
    DeviceData DATA = LoadS2PData(sourcepath);
    if (DATA.freq.at(0) == -1) return -1;//Prevents from corrupted file format
    fS = DATA.freq;
    ZS = DATA.Z.col(3);
    if (!ZL.empty())ResampleImpedances();
    return 0;
}

int MatchingNetwork::SetSourceImpedance(cx_vec zs, vec freq)
{
    if (zs.n_rows != freq.n_rows) return -1;
    ZS = zs;
    ZS_matching = zs;
    fS = freq;
    if (!ZL.empty())ResampleImpedances();
    return 0;
}

int MatchingNetwork::SetLoadImpedance(cx_vec zl, vec freq)
{
    if (zl.n_rows != freq.n_rows) return -1;
    ZL = zl;
    ZL_matching = ZL;
    fL = freq;
    if (!ZS.empty())ResampleImpedances();
    return 0;
}

int MatchingNetwork::SetLoadImpedance(std::string loadpath)
{
    DeviceData DATA = LoadS2PData(loadpath);
    if (DATA.freq.at(0) == -1) return -1;//Prevents from corrupted file format
    fL = DATA.freq;
    ZL = DATA.Z.col(0);
    if (!ZS.empty())ResampleImpedances();
    return 0;
}


int MatchingNetwork::SetInitialPivot(rowvec x)
{
    x_ini = x;
    return 0;
}

rowvec MatchingNetwork::GetInitialPivot()
{
    return x_ini;
}

int MatchingNetwork::SetMatchingBand(double f1, double f2, int N)
{
    f_matching = linspace(f1, f2, N);
    return 0;
}

int MatchingNetwork::SetTopology(std::string s)
{
    topology = s;
    return 0;
}


DeviceData MatchingNetwork::LoadS2PData(std::string filepath)
{
    DeviceData DATA;
    std::ifstream s2pfile(filepath.c_str());//Tries to open the data file.
    if(!s2pfile.is_open())//The data file cannot be opened => error
    {
        DATA.freq.at(0) = -1;
        return DATA;
    }

    std::string line;
    double freq_scale = 1;
    double Zref = 50;

    std::getline(s2pfile, line);
    while(line.compare(0, 1, "#"))//Looking for # field
    {
        std::getline(s2pfile, line);
    }

    //Freq scale
    if (!line.compare(2,1, "G"))
    {
        freq_scale = 1e9;
    }
    else
    {
        if (!line.compare(2,1, "M"))
        {
            freq_scale = 1e6;
        }
        else
        {
            if (!line.compare(2,1, "K"))
            {
                freq_scale = 1e3;
            }
        }
    }


    //Get the impedance at which the S params were measured

    int Rindex = line.find_last_of("R");
    Rindex = line.find_first_not_of(" ", Rindex);
    Zref = atof(line.substr(Rindex+1).c_str());


    while( (!getline(s2pfile, line)) || !line.compare(0,1, "!")|| line.length() == 1 );//Looking for the start of the raw data

    //DATA beginning.
    //At this point, the number of frequency samples is not known, so it's better to
    //push data into queues and then arrange it into armadillo structures
    std::queue <double> frequency, S11M, S11A, S21M, S21A, S12M, S12A, S22M, S22A;
    unsigned int qsize=0;

    do
    {
        if (line.at(0) == '!' || line.empty() || line.length() == 1 ) break;//Data end
        //Frequency
        int index2 = line.find_first_not_of(" ");
        int index1 = line.find_first_of(" ", index2+1);
        frequency.push(atof(line.substr(0,index1).c_str()));
        index2 = line.find_first_of(" ", index1+1);
        S11M.push(atof(line.substr(index1,index2-index1).c_str()));

        index1 = line.find_first_of(" ", index2+1);
        S11A.push(atof(line.substr(index2,index1-index2).c_str()));

        index2 = line.find_first_of(" ", index1+1);
        S21M.push(atof(line.substr(index1,index2-index1).c_str()));

        index1 = line.find_first_of(" ", index2+1);
        S21A.push(atof(line.substr(index2,index1-index2).c_str()));

        index2 = line.find_first_of(" ", index1+1);
        S12M.push(atof(line.substr(index1,index2-index1).c_str()));

        index1 = line.find_first_of(" ", index2+1);
        S12A.push(atof(line.substr(index2,index1-index2).c_str()));

        index2 = line.find_first_of(" ", index1+1);
        S22M.push(atof(line.substr(index1,index2-index1).c_str()));

        index1 = line.find_first_of(" ", index2+1);
        S22A.push(atof(line.substr(index2,index1-index2).c_str()));


        qsize++;
    }while (std::getline(s2pfile, line));

    vec freq(qsize);
    cx_mat S(qsize, 4);
    cx_mat Z(qsize, 4);
    double phi, S11m, S11a, S21m, S21a, S12m, S12a, S22m, S22a;
    for (unsigned int i = 0; i < qsize; i++)
    {
        freq(i) = freq_scale*frequency.front();
        frequency.pop();

        S11m = S11M.front();
        S11a = S11A.front();
        S11M.pop();
        S11A.pop();
        phi = (datum::pi/180)*S11a;
        S(i, 0) = cx_double(S11m,0)*cx_double(cos(phi), sin(phi));

        S21m = S21M.front();
        S21a = S21A.front();
        S21M.pop();
        S21A.pop();
        phi = (datum::pi/180)*S21a;
        S(i, 1) = cx_double(S21m,0)*cx_double(cos(phi), sin(phi));


        S12m = S12M.front();
        S12a = S12A.front();
        S12M.pop();
        S12A.pop();
        phi = (datum::pi/180)*S12a;
        S(i, 2) = cx_double(S12m,0)*cx_double(cos(phi), sin(phi));



        S22m = S22M.front();
        S22a = S22A.front();
        S22M.pop();
        S22A.pop();
        phi = (datum::pi/180)*S22a;
        S(i, 3) = cx_double(S22m,0)*cx_double(cos(phi), sin(phi));


        cx_double K = cx_double(Zref,0)/((1.-S(i,0))*(1.-S(i,3))-S(i,1)*S(i,2));
        Z(i, 0) = K*(1.+S(i, 0))*(1.-S(i,3))+S(i,1)*S(i,2);
        Z(i, 1) = K*2.*S(i,2);
        Z(i, 2) = K*2.*S(i,3);
        Z(i, 3) = K*(1.-S(i, 0))*(1.+S(i,3))+S(i,1)*S(i,2);


    }

    DATA.freq = freq;
    DATA.S = S;
    DATA.Z = Z;


    return DATA;
}


int MatchingNetwork::ResampleImpedances()
{
    //Check if the inputs lie in the s2p measurements
    double fmin, fmax;
    (fS.min() > fL.min()) ? fmin = fL.min() : fmin = fS.min();
    (fS.max() > fL.max()) ? fmax = fL.max() : fmax = fS.max();

    int N = f_matching.n_rows;
    if (fmin > f_matching.min())
    {
        std::cout << "Warning: Data not found for "<< f_matching.min() << ". The lowest frequency was set to "<< fmin << std::endl;
        f_matching  = linspace(fmin,f_matching.max(), N);
    }
    if (fmax < f_matching.max())
    {
        std::cout << "Warning: Data not found for "<< f_matching.max() << ". The highest frequency was set to "<< fmax << std::endl;
        f_matching  = linspace(f_matching.min(), fmax, N);
    }

    double N_density = N/(f_matching.max() - f_matching.min());
    f_analysis = linspace(fmin, fmax, floor(N_density*(fmax-fmin)));

    //Impedance interpolation. This way both vector have the same frequency basis.
    vec ZS_r(real(ZS)), ZS_i(imag(ZS)), ZL_r(real(ZL)), ZL_i(imag(ZL));
    vec ZS_inter_R, ZS_inter_I, ZL_inter_R, ZL_inter_I;
    interp1(fS, ZS_r, f_analysis, ZS_inter_R);//Armadillo lacks of spline interpolation
    interp1(fS, ZS_i, f_analysis, ZS_inter_I);
    interp1(fL, ZL_r, f_analysis, ZL_inter_R);
    interp1(fL, ZL_i, f_analysis, ZL_inter_I);

    ZS = cx_vec(ZS_inter_R, ZS_inter_I);
    ZL = cx_vec(ZL_inter_R, ZL_inter_I);


    //Select impedance values at the matching band
    vec V = abs(f_analysis - f_matching.min());
    uvec v_find = find(V < V.min()+1e-3);
    int index1 = v_find.at(0);
    V = abs(f_analysis - f_matching.max());
    v_find = find(V < V.min()+1e-3);
    int index2 = v_find.at(0);

    ZS_matching = ZS.rows(index1, index2);
    ZL_matching = ZL.rows(index1, index2);

    //Now, it's needed to u-pdate the matching band according to f_analysis values
    f_matching = f_analysis.rows(index1, index2);
}



cx_mat MatchingNetwork::getSparams(rowvec x, cx_double zs, cx_double zl, double f)
{
    int element;
    double w = 2*datum::pi*f;
    double beta = w/c0;
    cx_double gamma = cx_double(0, beta);
    cx_mat ABCD, ABCD_t, S;
    ABCD << 1 << 0 << endr << 0 << 1 << endr;
    S << -1 << -1 << endr << -1 << -1 << endr;

    int i, k;

    for (i = 0, k=0; i < topology.length(); i++, k++)
    {
        element = atoi(topology.substr(0,1).c_str());
        switch(element)
        {
        case 0: ABCD_t << 1 << cx_double(0,w*x.at(k)) << endr << 0 << 1 << endr;
            break;
        case 1: ABCD_t << 1 << cx_double(0,-1/w*x.at(k)) << endr << 0 << 1 << endr;
            break;
        case 2: ABCD_t << 1 << 0 << endr << cx_double(0,-1/w*x.at(k)) << 1 << endr;
            break;
        case 3: ABCD_t << 1 << 0 << endr << cx_double(0, w*x.at(k)) << 1 << endr;
            break;
        case 4: ABCD_t << cosh(gamma*x.at(k+1)) << x.at(k)*sinh(gamma*x.at(k+1)) << endr << sinh(gamma*x.at(k+1))/x(k) << cosh(gamma*x.at(k+1)) << endr;
            k++;//It involves two parameters, so we need to skip the next index
            break;
        case 5: ABCD_t << 1 << 0 << endr << (tanh(gamma*x.at(k+1)))/x.at(k) << 1 << endr;
            k++;
            break;
        case 6: ABCD_t << 1 << 0 << endr << 1./(x.at(k)*tanh(gamma*x.at(k+1))) << 1 << endr;
            k++;
            break;
        default: return S;
        }
        ABCD = ABCD*ABCD_t;
    }
    //Convert ABCD to S parameters
    S(0,0) = (ABCD(0,0)*zl+ABCD(0,1)-ABCD(1,0)*conj(zs)*zl-ABCD(1,1)*conj(zs))/(ABCD(0,0)*zl+ABCD(0,1)+ABCD(1,0)*zs*zl+ABCD(1,1)*zs);
    S(0,1) = (2.*(ABCD(0,0)*ABCD(1,1)-ABCD(0,1)*ABCD(1,0))*sqrt(real(zs)*real(zl)))/(ABCD(0,0)*zl+ABCD(0,1)+ABCD(1,0)*zs*zl+ABCD(1,1)*zs);
    S(1,0) = (2.*sqrt(real(zs)*real(zl)))/(ABCD(0,0)*zl+ABCD(0,1)+ABCD(1,0)*zs*zl+ABCD(1,1)*zs);
    S(1,1) = (-ABCD(0,0)*conj(zl)+ABCD(0,1)-ABCD(1,0)*conj(zl)*zs+ABCD(1,1)*zs)/(ABCD(0,0)*zl+ABCD(0,1)+ABCD(1,0)*zs*zl+ABCD(1,1)*zs);
}

int MatchingNetwork::SetMaxIterGridSearch(int max_iter)
{
    Grid_MaxIter = max_iter;
    return 0;
}

int MatchingNetwork::SetThreshold(double th)
{
    MatchingThreshold = th;
    return 0;
}

double MatchingNetwork::GetThreshold()
{
    return MatchingThreshold;
}


GRABIM_Result MatchingNetwork::RunGRABIM()
{
    GRABIM_Result Res;
    Res.x_grid_search = GridSearch();
    Res.x_nlopt = LocalOptimiser(Res.x_grid_search);
    Res.grid_val = CandidateEval(Res.x_grid_search);
    Res.nlopt_val = CandidateEval(Res.x_nlopt);
    return Res;
}



rowvec MatchingNetwork::GridSearch()
{
    vec delta_k, best, step;
    int dim = x_ini.n_cols;
    rowvec f_xkq(pow(2,dim)+1);
    rowvec xk(x_ini);
    mat xkq = ones(pow(2, dim)+1, dim);
    mat best_pivot = ones(4, dim);
    delta_k << 1. << 1./4 << 1./16 << 1./64 << endr;
    best << 100<< 100<<100<<100<< endr;
    double fxk = CandidateEval(xk);
    double best_candidate;
    best.at(0) = fxk;
    mat C = GeneratingMatrix(dim);
    unsigned int niter = 0;

    for (int i = 0; i<4; i++, niter++)
    {
        if ((best.min() < MatchingThreshold)&&(niter > Grid_MaxIter))
        {
            return xk;
        }

        for (int q = 0; q<=pow(2, dim);q++)
        {
            step = delta_k.at(i)*C.col(q);
            xkq.row(q) = xk%(1. + step.t()*.99);
            f_xkq.at(q) = CandidateEval(xkq.row(q));
        }
        best_candidate = f_xkq.min();
        uvec v_find = find(f_xkq < best_candidate+1e-6);
        int imin = v_find.at(0);

        if (verbose)std::cout << "BEST CANDIDATE "<< i << ": " << xkq.row(imin) << " => " << best_candidate << std::endl;


        if (best_candidate < best.min())
        {
            best << best_candidate << 0 << 0 << 0 << endr;//Prepare best for new hypercube cycle
            xk = xkq.row(imin);//Set new pivot
            i = -1;//New pivot
            if (verbose) std::cout << "New pivot: " << xk << std::endl;
            continue;
        }
        else
        {
            best.at(i) = best_candidate;
            best_pivot.row(i) = xk;
        }

        if (i > 0)//Checks if the i-th cycle improved the result with regard to the (i-1)-th cycle
        {
            if ((abs(best.at(i)-best.at(i-1))<1e-3)||(best(i)>best.at(i-1)))//The current cycle did not improve the result
            {
                uvec v_find = find(best < best.min()+1e-6);

                //Suggestion for improvement. When GRABIM gets stuck, sometimes it helps to try some
                //random vectors near the current pivot. The distribution was chosen to be U(-0.5, 0.5)
                int Ncatchup  = 100;
                mat XKQ = ones(Ncatchup, xk.n_cols);
                vec FXK = ones(Ncatchup);
                for (int k = 0; k <Ncatchup; k++)
                {
                    XKQ.row(k) = best_pivot.row(0) %(randu(1, xk.n_cols)+.5);
                    FXK.at(k) =CandidateEval(XKQ.row(k));
                }
                uvec v_find_catchup = find(FXK < FXK.min()+1e-6);
                int imin = v_find_catchup.at(0);
                if (FXK.min() < best.min())
                {
                xk = XKQ.row(imin);
                if (verbose)std::cout << "=> " << xk << ": " << FXK.min()<< std::endl;
                i=-1;
                }
            }

        }
    }
    std::cout << best_candidate << std::endl;
    return xk;
}

mat MatchingNetwork::GeneratingMatrix(int dim)
{
    mat C = zeros(dim, pow(2, dim)+1);
    for (int i = 0; i <dim; i++)
    {
        mat V =  ones(1, pow(2, (i+1)));
        V.cols(0, pow(2, i)-1) = ones(1, pow(2, i));
        V.cols(pow(2,i), pow(2, i+1)-1) = -ones(1, pow(2, i));
        C.row(i).cols(0, pow(2, dim)-1) = repmat(V, 1., pow(2, dim-i-1));
    }
    return C;
}

double MatchingNetwork::CandidateEval(rowvec x)
{
    double S11_max = -1e3;
    cx_mat S;
    for (int i = 0; i < f_matching.n_rows; i++)
    {
        S = getSparams(x, ZS_matching.at(i), ZL_matching.at(i), f_matching.at(i));
        if (abs(S(0,0)) > S11_max) S11_max = abs(S(0,0));
    }
    return 20*log10(S11_max);
}


typedef struct NLoptData {
    vec f_matching;
    cx_vec ZS_matching, ZL_matching;
    std::string topology;

} NLoptData;

double myfunc(const std::vector<double> &x, std::vector<double> &grad, void *n)
{
    MatchingNetwork M;
    rowvec x_ = ones(1,x.size());
    NLoptData * N = (NLoptData *) n;
    M.SetTopology(N->topology);
    M.SetMatchingBand(N->f_matching.min(), N->f_matching.max(), N->f_matching.n_rows);
    M.SetSourceImpedance(N->ZS_matching, N->f_matching);
    M.SetLoadImpedance(N->ZL_matching, N->f_matching);
    for (int i = 0; i < x.size(); i++) x_.at(i) = x[i];

    double eval = M.CandidateEval(x_);
    std::cout<< eval << " <= " << x_<< std::endl;
    return eval;
}

rowvec MatchingNetwork::LocalOptimiser(rowvec x_grid)
{
    int dim = x_ini.n_cols;
    nlopt::opt opt(nlopt::LN_PRAXIS, dim);
    NLoptData n;
    n.f_matching = f_matching;
    n.topology = topology;
    n.ZL_matching = ZL_matching;
    n.ZS_matching = ZS_matching;
    opt.set_min_objective(myfunc, &n);
    opt.set_stopval(this->GetThreshold()-10);
    opt.set_maxeval(1e6);
    opt.set_xtol_rel(1e-4);

    //Bounds
    std::vector<double> lb(dim), ub(dim);
    int i, k;
    for (i = 0, k=0; i < topology.length(); i++, k++)
    {
        int element = atoi(topology.substr(0,1).c_str());
        if (element < 4)
        {
            lb[k] = 1e-12;
            ub[k] = 1e-6;
        }
        else
        {
            lb[k] = 5;   lb[k+1] = 1e-6;
            ub[k] = 1e4; ub[k+1] = 10;
            k++;
        }
    }


    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);

    std::vector<double> x(dim);
    for (int i = 0; i < dim; i++)x[i] = x_grid.at(i);
    double minf;
    nlopt::result result = opt.optimize(x, minf);
    std::cout << "NLopt status: " << std::endl;

    switch(result)
    {
    //Success
        case 1: std::cout << "Success!" <<std::endl;
                break;
        case 2: std::cout << "Goal achieved" <<std::endl;
                break;
        case 3: std::cout << "The algorithm reached convergence (relative)" <<std::endl;
                break;
        case 4: std::cout << "The algorithm reached convergence (absolute)" <<std::endl;
                break;
        case 5: std::cout << "Maximum number of evaluations reached" <<std::endl;
                break;
        case 6: std::cout << "Maximum time reached" <<std::endl;
            break;
     //Failure
        case -1: std::cout << "NLopt failed!" <<std::endl;
                 break;
        case -2: std::cout << "Invalid arguments" <<std::endl;
                 break;
        case -3: std::cout << "There is not enough available memory" <<std::endl;
                 break;
        case -4: std::cout << "Forced termination" <<std::endl;
                 break;
    }

    rowvec x_nlopt = ones(1, dim);
    for (int i = 0; i < dim; i++)x_nlopt.at(i) = x[i];

    return x_nlopt;

}
