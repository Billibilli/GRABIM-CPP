#include "GRABIM.h"

GRABIM::GRABIM()
{
    verbose  = false;
    Grid_MaxIter = 1000;
    MatchingThreshold = -30;//It specifies the mininum S11 [dB] required, typically,
    //S11 < 10 dB is considered as valid for common applications

    ObjFun = ObjectiveFunction::NINF_S11dB;//Sets the kind of objective functions. Everything seems to
    //suggest that NINF_S11dB gives the best results
    topology = "-1";//By default, the engine will search the optimum topology from the predefined set of circuit
    NLoptAlgo = nlopt::LN_NELDERMEAD;//By default, the local optimiser was chosen to be the Nelder-Mead algorithm
}

// It sets if the program should show detailed information about the iterations or not
int GRABIM::SetVerbose(bool vb)
{
    verbose = vb;
    return 0;
}



// It sets the source impedance vs frequency
int GRABIM::SetSourceImpedance(cx_vec zs)
{
    ZS = zs;
    return 0;
}

// It loads a s2p file and calculates the input impedance
int GRABIM::SetLoadImpedance(cx_vec zl)
{
    ZL = zl;
    return 0;
}

// It sets the load impedance vs frequency

int GRABIM::SetFrequency(vec f)
{
    freq = f;
    return 0;
}

// Sets the first point of the grid search
int GRABIM::SetInitialPivot(rowvec x)
{
    x_ini = x;
    return 0;
}

// Returns the first point of the grid search
rowvec GRABIM::GetInitialPivot()
{
    return x_ini;
}



// This function sets the ladder arrangement of the matching network
int GRABIM::SetTopology(std::string s)
{
    topology = s;
    return 0;
}


// Sets the maximum number of iterations for the grid search engine
int GRABIM::SetMaxIterGridSearch(int max_iter)
{
    Grid_MaxIter = max_iter;
    return 0;
}

// Sets the threshold of what the program consideres a good matching
int GRABIM::SetThreshold(double th)
{
    MatchingThreshold = th;
    return 0;
}

// Returns the matching threshold
double GRABIM::GetThreshold()
{
    return MatchingThreshold;
}


GRABIM_Result GRABIM::RunGRABIM()
{
    GRABIM_Result Res;
    if (!topology.compare("-1"))//The user did not entered any specific network, so it seems
    {//reasonable to try some typical wideband matching network
        rowvec Vopt;
        string candidate;

        SearchPredefinedTopologies(Vopt, candidate);

        topology = candidate;
        Res.x_grid_search = Vopt;
    }
    else//It builds an initial vector using default values 1nH for inductances, 1 pF for capacitors, 100 Ohm & lambda/4 for transmission lines
    {
        AutoSetInitialPivot();
        Res.x_grid_search = GridSearch();
    }

    Res.grid_val = CandidateEval(Res.x_grid_search);
    Res.x_nlopt = LocalOptimiser(Res.x_grid_search);
    Res.nlopt_val = CandidateEval(Res.x_nlopt);
    Res.ZS = ZS;
    Res.ZL = ZL;
    Res.freq = freq;

    //Initialize S param vectors
    Res.S11_gridsearch=cx_vec(freq,freq);
    Res.S21_gridsearch=cx_vec(freq,freq);
    Res.S12_gridsearch=cx_vec(freq,freq);
    Res.S22_gridsearch=cx_vec(freq,freq);

    Res.S11_nlopt=cx_vec(freq,freq);
    Res.S21_nlopt=cx_vec(freq,freq);
    Res.S12_nlopt=cx_vec(freq,freq);
    Res.S22_nlopt=cx_vec(freq,freq);


    //Generate S parameter results
    cx_mat S_gridsearch, S_nlopt;
    S_gridsearch << 1<<1<<endr << 1<<1<< endr;
    S_nlopt << 1<<1<<endr << 1<<1<< endr;

    SparEngine S2PEngine;

    for (unsigned int i = 0; i < freq.n_rows; i++)
    {
        // Grid search
        S_gridsearch = S2PEngine.getSparams(Res.x_grid_search, ZS.at(i), ZL.at(i), freq.at(i), topology);
        Res.S11_gridsearch.at(i) = S_gridsearch(0,0);
        Res.S21_gridsearch.at(i) = S_gridsearch(1,0);
        Res.S12_gridsearch.at(i) = S_gridsearch(0,1);
        Res.S22_gridsearch.at(i) = S_gridsearch(1,1);

        // NLopt
        S_nlopt = S2PEngine.getSparams(Res.x_nlopt, ZS.at(i), ZL.at(i), freq.at(i), topology);
        Res.S11_nlopt.at(i) = S_nlopt(0,0);
        Res.S21_nlopt.at(i) = S_nlopt(1,0);
        Res.S12_nlopt.at(i) = S_nlopt(0,1);
        Res.S22_nlopt.at(i) = S_nlopt(1,1);
    }

    Res.topology = topology.c_str();//Save network topology to create a Qucs schematic

    return Res;
}



rowvec GRABIM::GridSearch()
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

    f_xkq = 1e6*f_xkq.ones();

    /*
    //Initial guess. Sometimes, the initial point may be far away from the optimum, so it
    // seems to be reasonable to pick a handful of random vectors nearby the inital point and
    // choose the best one.

    int N_initial_guess  = 100;
    mat XKQ = ones(N_initial_guess, xk.n_cols);
    vec FXK = ones(N_initial_guess);
    for (int k = 0; k <N_initial_guess; k++)
    {
        XKQ.row(k) = xk %(randu(1, xk.n_cols)*3);
        FXK.at(k) =CandidateEval(XKQ.row(k));
    }
    uvec v_find_guess = find(FXK < FXK.min()+1e-6);
    int imin = v_find_guess.at(0);
    if (FXK.min() < fxk)
    {
        xk = XKQ.row(imin);
        fxk = CandidateEval(xk);
    }*/

    // Factorial search
    for (int i = 0; i<4; i++, niter++)
    {
        if ((best.min() < MatchingThreshold)&&(niter > Grid_MaxIter))//Stop condition
        {
            return xk;
        }

        dim = xk.n_cols;
        CheckDim(C, xkq, best_pivot, dim);//InspectCandidate() may have shrunk the network, so it is needed to update the settings

        for (int q = 0; q<=pow(2, dim);q++)//Calculates the vertices of the searching hypercube
        {
            step = delta_k.at(i)*C.col(q);
            xkq.row(q) = xk%(1. + step.t()*.99);
            f_xkq.at(q) = CandidateEval(xkq.row(q));
        }
        best_candidate = f_xkq.min();
        uvec v_find = find(f_xkq < best_candidate+1e-6);
        int imin = v_find.at(0);

        if (verbose)std::cout << "BEST CANDIDATE "<< i << ": " << xkq.row(imin) << " => " << best_candidate << std::endl;


        if (best_candidate < best.min() - 1e-3)
        {
            best << best_candidate << 0 << 0 << 0 << endr;//Prepare best for new hypercube cycle
            xk = xkq.row(imin);//Set new pivot
            xk = InspectCandidate(xk);
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
            if (best.at(i)-best.at(i-1)<-1e-3)//The current cycle did not improve the result
            {
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
                    xk = InspectCandidate(xk);
                    if (verbose)std::cout << "=> " << xk << ": " << FXK.min()<< std::endl;
                    i=-1;
                }
            }

        }
    }
    return xk;
}



//This function searches the optimum over a predefined circuit topologies
int GRABIM::SearchPredefinedTopologies(rowvec & Vopt, std::string & candidate)
{
    double gridtest, opttopo = 1e12;
    std::ifstream TopologiesFile(TopoScript_path);//Tries to open the file.
    std::string line, tag;
    rowvec Vaux;
    if(!TopologiesFile.is_open())//The data file cannot be opened => error
    {
        return -1;
    }

    cout << "Predefined networks..." << endl;

    vector<string> strlist;
    string aux;
    while (std::getline(TopologiesFile, aux))
    {
        if (aux.empty())continue;
        if (aux.substr(0,1).compare("#")==0)
        {
            tag = aux;
            continue;
        }
        topology = aux;
        std::getline(TopologiesFile, line);

        if (!line.empty())
        {//Check whether the initial pivot is given in the script or not
            istringstream f(line);
            std::string s;
            while (getline(f, s, ';'))  strlist.push_back(s);

            x_ini.resize(1, (int) strlist.size());
            for (unsigned int i = 0; i < strlist.size(); i++)
            {
                x_ini.at(i) = atof(strlist.at(i).c_str());
            }
        }
        else
        {//The script does not provide an initial pivot, so we need to use
         // the default one.
            AutoSetInitialPivot();
        }
        Vaux = GridSearch();
        gridtest = CandidateEval(Vaux);
        cout << tag << ": S11_min = " <<  gridtest << " dB" << endl;

        if (gridtest < opttopo - 0.5)//It worths to change the topology if the new one improves the result significantly
        {
            candidate = topology;
            opttopo = gridtest;
            Vopt = Vaux;
        }
    }
    return 0;
}

// Whenever the impedance of a component is not significant (e.g. low series impedance or high shunt impedance) its
// better to remove this component from the optimisation problem. This will lead to an improvement in terms of better convergence
void GRABIM::removeElement(rowvec & xk, unsigned int pos, unsigned int postopo)
{
    rowvec aux(xk.n_elem-1);
    int index = 0;
    int element = atoi(topology.substr(postopo,1).c_str());
    for (unsigned int i = 0; i < xk.n_elem; i++)//Removes component values from the working vector
    {
        if (element == 4)
        {
            if ((i == pos) || (i == pos+1)) continue;
        }
        if (i == pos) continue;
        aux.at(index) = xk.at(i);
        index++;
    }
    topology.erase(postopo, 1);//Removes the component from the topology indicator

    //Reshape xk
    xk.resize(aux.size());
    xk=aux;
}

//This function is aimed to remove reduncancies in case some network element
//is removed because of the irrelevance of its impedance
void GRABIM::CheckNetwork(rowvec & xk, unsigned int x_index, unsigned int topo_index)
{
    //The following cases must be treated separatedly
    // * Two series inductors
    // * Two series capacitors
    // * Two parallel inductors
    // * Two parallel capacitors

    if ((topo_index == 0) || (topo_index == topology.length()-1)) return;//It was removed the first or the last element. Any case below applies
    if (!strcmp(topology.substr(topo_index-1, 2).c_str(), "00"))//Two series inductors
    {
       xk.at(x_index-1) += xk.at(x_index);
       removeElement(xk, x_index, topo_index);
       return;
    }

    if (!strcmp(topology.substr(topo_index-1, 2).c_str(), "11"))//Two series capacitors
    {
       xk.at(x_index-1) = (xk.at(x_index)*xk.at(x_index-1))/(xk.at(x_index)+xk.at(x_index-1));
       removeElement(xk, x_index, topo_index);
       return;
    }

    if (!strcmp(topology.substr(topo_index-1, 2).c_str(), "00"))//Two parallel inductors
    {
       xk.at(x_index-1) = (xk.at(x_index)*xk.at(x_index-1))/(xk.at(x_index)+xk.at(x_index-1));
       removeElement(xk, x_index, topo_index);
       return;
    }

    if (!strcmp(topology.substr(topo_index-1, 2).c_str(), "11"))//Two parallel capacitors
    {
       xk.at(x_index-1) += xk.at(x_index);
       removeElement(xk, x_index, topo_index);
       return;
    }

}


// This function checks whether the candidate vector contains some
// irrelevant value or not. Unsignificant values may prevent the algorithm
// to find a better or more realistic solution. In this sense, shunt elements
// whose impedance exceeds 4kOhm are interpreted as open circuits, and conversely,
// series components whose impedance is below 1 Ohm are treated as short circuits
rowvec GRABIM::InspectCandidate(rowvec xk)
{
    double impedance, fmax = freq.max(), fmin = freq.min();
    unsigned int element;
    double wmax =2*datum::pi*fmax, wmin =2*datum::pi*fmin;
    int index = 0;
    for (unsigned int i = 0; i < topology.length(); i++)//Inspects all elements. aux is defined so as to take into account
    {                                                   // those elements which depend on more than one variable (e.g.) TL
        element = atoi(topology.substr(i,1).c_str());
        switch(element)
        {
        case 0://Series inductor
            impedance = wmax*xk.at(index);
            if (impedance < 1)
            {
                xk.at(index) = 1e-100;
                removeElement(xk, index, i);//The element seems not to have an important effect on the network, so it can
                                            // be removed
                CheckNetwork(xk, index, i);//Avoid network redundancies (parallel/series capacitors/inductances)
            }
            if (xk.at(index) < 0)//Inductance must be > 0, so it seems that a capacitor would do a better job
            {
                xk.at(index) = 1./(2e3*wmin);//High impedance capacitor
                topology[i] = '1';//Series capacitor
                cout << "Warning: The selected topology leads to L < 0" <<endl;
                cout << "Warging: Topology changed. The new topology is: " << topology << endl;
            }
            index++;
            continue;
        case 2://Shunt inductor
            impedance = wmin*xk.at(index);
            if (impedance > 4e3)
            {
                xk.at(index) = 1e10;
                removeElement(xk, index, i);
                CheckNetwork(xk, index, i);//Avoid network redundancies (parallel/series capacitors/inductances)
            }
            if (xk.at(index) < 0)//Inductance must be > 0, so it seems that a capacitor would do a better job
            {
                xk.at(index) = 1./(2e3*wmin);//High impedance capacitor
                topology[i]='3';//Shunt capacitor
                cout << "Warning: The selected topology leads to L < 0" <<endl;
                cout << "Warging: Topology changed. The new topology is: " << topology << endl;
            }
            index++;
            continue;
        case 1://Series capacitor
            impedance = 1./(wmin*xk.at(index));
            if (impedance < 1)
            {
                xk.at(index) = 100;
                removeElement(xk, index, i);
                CheckNetwork(xk, index, i);//Avoid network redundancies (parallel/series capacitors/inductances)
            }
            if (xk.at(index) < 0)//C must be > 0, so it seems that an inductor would do a better job
            {
                xk.at(index) = 5./(wmax);//Low impedance inductance
                topology[i] = '0';//Series inductance
                cout << "Warning: The selected topology leads to C < 0" <<endl;
                cout << "Warging: Topology changed. The new topology is: " << topology << endl;
            }
            index++;
            continue;
        case 3://Shunt capacitor
            impedance = 1./(wmax*xk.at(index));
            if (impedance > 4e3)
            {
                xk.at(index) = 1e-30;
                removeElement(xk, index, i);
                CheckNetwork(xk, index, i);//Avoid network redundancies (parallel/series capacitors/inductances)
            }
            if (xk.at(index) < 0)//C must be > 0, so it seems that an inductor would do a better job
            {
                xk.at(index) = 5./(wmax);//Low impedance inductance
                topology[i] = '2';//Shunt inductance
                cout << "Warning: The selected topology leads to C < 0" <<endl;
                cout << "Warging: Topology changed. The new topology is: " << topology << endl;
            }
            index++;
            continue;
        default://It is a transmission line or a stub
            double zmax = max(abs(ZS));
            if (zmax < max(abs(ZL))) zmax = max(abs(ZL));
            if ((xk.at(index) > 5*zmax) || (xk.at(index+1) < 0))//Something is wrong...
            {
                return -1*ones(1, xk.n_cols);
            }
            index += 2;//Two parameters: Z0 and length


        }
    }
    return xk;
}


mat GRABIM::GeneratingMatrix(int dim)
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

double GRABIM::CandidateEval(rowvec x)
{
    double fobj = -1e3;
    cx_mat S, ABCD;
    SparEngine S2PEngine;
    for (unsigned int i = 0; i < freq.n_elem; i++)
    {
        if (ObjFun == ObjectiveFunction::NINF_S11dB)
        {
            S = S2PEngine.getSparams(x, ZS.at(i), ZL.at(i), freq.at(i), topology);
            if (abs(S(0,0)) > fobj) fobj = abs(S(0,0));
        }
        if (ObjFun == ObjectiveFunction::NINF_POWERTRANS)
        {
            ABCD = S2PEngine.getABCDmatrix(x, freq.at(i), topology);
            fobj = CalcInvPowerTransfer(ABCD, ZS.at(i), ZL.at(i));
        }
    }
    if (ObjFun == ObjectiveFunction::NINF_S11dB)fobj = 20*log10(fobj);//|grad{log(x)}| > |grad{x}| when x < 1;
    return fobj;
}


typedef struct NLoptData {
    vec freq;
    cx_vec ZS, ZL;
    std::string topology;
    ObjectiveFunction objf;

} NLoptData;

double myfunc(const std::vector<double> &x, std::vector<double> &grad, void *n)
{
    GRABIM M;
    rowvec x_ = ones(1,x.size());
    NLoptData * N = (NLoptData *) n;
    M.SetTopology(N->topology);
    M.SetFrequency(N->freq);
    M.SetSourceImpedance(N->ZS);
    M.SetLoadImpedance(N->ZL);
    M.SetObjectiveFunction(N->objf);
    for (unsigned int i = 0; i < x.size(); i++) x_.at(i) = x[i];

    double eval = M.CandidateEval(x_);
    std::cout<< eval << " <= " << x_<< std::endl;
    return eval;
}

// Given the grid optimum, the local optimiser is supposed to refine the search
rowvec GRABIM::LocalOptimiser(rowvec x_grid)
{
    int dim = x_grid.n_cols;
    nlopt::opt opt(NLoptAlgo, dim);
    NLoptData n;
    n.topology = topology;
    n.ZL = ZL;
    n.ZS = ZS;
    n.objf = ObjFun;
    n.freq = freq;
    opt.set_min_objective(myfunc, &n);
    opt.set_stopval(this->GetThreshold()-10);
    opt.set_maxeval(1e6);
    // opt.set_xtol_rel(1e-4);

    //Bounds
    std::vector<double> lb(dim), ub(dim);
    unsigned int i, k;
    for (i = 0, k=0; i < topology.length(); i++, k++)//Automatically generates the boundaries
    {
        int element = atoi(topology.substr(i,1).c_str());
        if (element < 4)
        {
            lb[k] = 1e-20;
            ub[k] = 1e-6;
        }
        else
        {
            lb[k] = 5;   lb[k+1] = 1e-6;
            ub[k] = 1e4; ub[k+1] = 10;
            k++;
        }
    }


    // It seems that NLopt sometimes crashes because of the limits...
    // opt.set_lower_bounds(lb);
    // opt.set_upper_bounds(ub);

    std::vector<double> x(dim);
    for (int i = 0; i < dim; i++)x[i] = x_grid.at(i);
    double minf;

    cout << endl << "Local optimiser output:" << endl;
    nlopt::result result = opt.optimize(x, minf);
    std::cout << "NLopt status: " << std::endl;

    //NLopt exit status
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
    case -5: std::cout << "Forced termination" <<std::endl;
        break;
    }

    rowvec x_nlopt = ones(1, dim);
    for (int i = 0; i < dim; i++)x_nlopt.at(i) = x[i];

    return x_nlopt;

}

// Sets the local optimiser algorithm
int GRABIM::SetNLoptAlg(nlopt::algorithm NLA)
{
    NLoptAlgo = NLA;
    return 0;
}
//Gets the name of the local optimiser
nlopt::algorithm GRABIM::GetNLoptAlg()
{
    return NLoptAlgo;
}


// Sets the objective function. The problem can be thought in terms of either S11 or the inverse power transfer
int GRABIM::SetObjectiveFunction(ObjectiveFunction of)
{
    ObjFun = of;
    return 0;
}

//Returns the objective function
ObjectiveFunction GRABIM::GetObjectiveFunction()
{
    return ObjFun;
}

// Inverse power transfer from ABCD matrix
// [1] Eq (6.2.1), (6.2.2)
double GRABIM::CalcInvPowerTransfer(cx_mat ABCD, cx_double ZS, cx_double ZL)
{
    cx_double p = real(ZS)*real(ZL) - imag(ZS)*imag(ZL);
    cx_double q = imag(ZS)*real(ZL) + imag(ZL)*real(ZS);
    cx_double a = ABCD(0,0)*real(ZL) - ABCD(1,0)*q + ABCD(1,1)*real(ZS);
    cx_double b = ABCD(0,1) + ABCD(1,0)*p + ABCD(1,1)*imag(ZS) + ABCD(0,0)*imag(ZL);
    return abs((a*a + b*b)/(4*real(ZS)*real(ZL)));
}





void GRABIM::AutoSetInitialPivot()
{
    double meanf = .5*(freq.min()+freq.max());
    double meanw = 2*datum::pi*meanf;
    double lambda4 = c0/(4.*meanf);
    queue <double> XINI;
    for (unsigned int i = 0; i< topology.size();i++)
    {
        if (!topology.substr(i,1).compare("0")) XINI.push(20/meanw);//Series impedance at midband ~ 5 Ohm
        if (!topology.substr(i,1).compare("2")) XINI.push(300/meanw);//Parallel impedance at midband ~ 100 Ohm
        if (!topology.substr(i,1).compare("1")) XINI.push(1/(20*meanw));//Series impedance at midband ~ 5 Ohm
        if (!topology.substr(i,1).compare("3")) XINI.push(1/(300*meanw));//Parallel impedance at midband ~ 100 Ohm
        if((!topology.substr(i,1).compare("5"))||(!topology.substr(i,1).compare("6")))XINI.push(real(mean(ZS+ZL))),XINI.push(lambda4);

        if (!topology.substr(i,1).compare("4"))//Transmission line
        {
            // Here it is defined a linear impedance profile. In case the user selected 444, this tries
            // works real-to-real impedance transformer. Luckily, GRABIM will find a better result
            double Zi;
            double m = 2.*(mean(real(ZL))-mean(real(ZS)))/topology.length();

            Zi = m*i+mean(real(ZS));
            XINI.push(Zi);
            XINI.push(lambda4);
        }
    }
    x_ini = ones(1, XINI.size());
    for (unsigned int i = 0; i < x_ini.size();i++)
    {
        x_ini.at(i) = XINI.front();
        XINI.pop();
    }
}


void GRABIM::setTopoScript(std::string path)
{
    TopoScript_path = path;
}


//
void GRABIM::CheckDim(mat & C, mat & xkq, mat & best_pivot, unsigned int dim)
{
    if (xkq.n_cols != dim)//The current network has changed
    {
      C.resize(dim, pow(2, dim)+1);
      C = GeneratingMatrix(dim);
      xkq.resize(pow(2, dim)+1, dim);
      xkq.ones();
      best_pivot.resize(4, dim);
      best_pivot.ones();
    }
}
