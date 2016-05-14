#include "sparengine.h"

SparEngine::SparEngine()
{
    preTopoList.insert("444");
    preTopoList.insert("020313");
    preTopoList.insert("313020");
    preTopoList.insert("03030");
    preTopoList.insert("313202");
    preTopoList.insert("202313");
    preTopoList.insert("03030");
    preTopoList.insert("43434");
    preTopoList.insert("01230123");
    preTopoList.insert("23012301");
    preTopoList.insert("12031203");
    preTopoList.insert("21302130");
    preTopoList.insert("03120312");
    preTopoList.insert("30213021");
    preTopoList.insert("30303030");
    preTopoList.insert("03030303");
    preTopoList.insert("12121212");
    preTopoList.insert("21212121");
}


// Returns the S matrix at a given frequency
cx_mat SparEngine::getSparams(rowvec x, cx_double zs, cx_double zl, double f, std::string topology)
{
    cx_mat ABCD = getABCDmatrix(x, f, topology);
    cx_mat S;
    S << -1 << -1 << endr << -1 << -1 << endr;
    //Convert ABCD to S parameters
    S(0,0) = (ABCD(0,0)*zl+ABCD(0,1)-ABCD(1,0)*conj(zs)*zl-ABCD(1,1)*conj(zs))/(ABCD(0,0)*zl+ABCD(0,1)+ABCD(1,0)*zs*zl+ABCD(1,1)*zs);
    S(0,1) = (2.*(ABCD(0,0)*ABCD(1,1)-ABCD(0,1)*ABCD(1,0))*sqrt(real(zs)*real(zl)))/(ABCD(0,0)*zl+ABCD(0,1)+ABCD(1,0)*zs*zl+ABCD(1,1)*zs);
    S(1,0) = (2.*sqrt(real(zs)*real(zl)))/(ABCD(0,0)*zl+ABCD(0,1)+ABCD(1,0)*zs*zl+ABCD(1,1)*zs);
    S(1,1) = (-ABCD(0,0)*conj(zl)+ABCD(0,1)-ABCD(1,0)*conj(zl)*zs+ABCD(1,1)*zs)/(ABCD(0,0)*zl+ABCD(0,1)+ABCD(1,0)*zs*zl+ABCD(1,1)*zs);
    return S;
}


// Returns the ABCD matrix at a given frequency
cx_mat SparEngine::getABCDmatrix(rowvec x, double f, std::string topology)
{
    int element;
    double w = 2*datum::pi*f;
    double beta = w/c0;
    cx_double gamma = cx_double(0, beta);
    cx_mat ABCD, ABCD_t;
    ABCD << 1 << 0 << endr << 0 << 1 << endr;

    unsigned int i, k;



     ABCD = PreComputedABCD(x, w, topology);
     //cout << ABCD << endl << endl;
     if (preTopoList.find(topology) != preTopoList.end())
     {
         return ABCD;
     }

     ABCD << 1 << 0 << endr << 0 << 1 << endr;


    for (i = 0, k=0; i < topology.length(); i++, k++)
    {
        element = atoi(topology.substr(i,1).c_str());
        switch(element)
        {
        case 0: ABCD_t << 1. << cx_double(0,w*x.at(k)) << endr << 0 << 1. << endr;
            break;
        case 1: ABCD_t << 1. << cx_double(0,-1/(w*x.at(k))) << endr << 0 << 1. << endr;
            break;
        case 2: ABCD_t << 1. << 0 << endr << cx_double(0,-1./(w*x.at(k))) << 1. << endr;
            break;
        case 3: ABCD_t << 1. << 0 << endr << cx_double(0, w*x.at(k)) << 1. << endr;
            break;
        case 4: ABCD_t << cosh(gamma*x.at(k+1)) << x.at(k)*sinh(gamma*x.at(k+1)) << endr << sinh(gamma*x.at(k+1))/x(k) << cosh(gamma*x.at(k+1)) << endr;
            k++;//It involves two parameters, so we need to skip the next index
            break;
        case 5: ABCD_t << 1. << 0 << endr << (tanh(gamma*x.at(k+1)))/x.at(k) << 1. << endr;
            k++;
            break;
        case 6: ABCD_t << 1. << 0 << endr << 1./(x.at(k)*tanh(gamma*x.at(k+1))) << 1. << endr;
            k++;
            break;
        default: return -ABCD.eye();
        }

        ABCD = ABCD*ABCD_t;
    }
  //  cout << ABCD << endl;
        return ABCD;
}

// This function contains the ABCD matrices of the most common wideband matching networks. This way, it is not neccessay to
// compute the overall ABCD matrix as a matrix product of single ABCD matrices and consequently speeds up execution
cx_mat SparEngine::PreComputedABCD(rowvec x, double w, std::string topology)
{
    cx_mat ABCD(2,2);
    cx_double I = cx_double(0,1);

    ABCD << -1 << -1 << endr << -1 << -1 << endr;

    if (!topology.compare("444"))// Lowpass TL-C ladder
    {
        ABCD(0,0) = -(x.at(0)*x.at(2)*cos(x.at(3)*w/c0)*sin(x.at(1)*w/c0)*sin(x.at(5)*w/c0) + x.at(2)*x.at(2)*cos(x.at(1)*w/c0)*sin(x.at(3)*w/c0)*sin(x.at(5)*w/c0) - (x.at(2)*cos(x.at(1)*w/c0)*cos(x.at(3)*w/c0)*cos(x.at(5)*w/c0) - x.at(0)*cos(x.at(5)*w/c0)*sin(x.at(1)*w/c0)*sin(x.at(3)*w/c0))*x.at(4))/(x.at(2)*x.at(4));
        ABCD(0,1) = (-I*x.at(0)*x.at(2)*cos(x.at(3)*w/c0)*cos(x.at(5)*w/c0)*sin(x.at(1)*w/c0) - I*x.at(2)*x.at(2)*cos(x.at(1)*w/c0)*cos(x.at(5)*w/c0)*sin(x.at(3)*w/c0) + (-I*x.at(2)*cos(x.at(1)*w/c0)*cos(x.at(3)*w/c0)*sin(x.at(5)*w/c0) + I*x.at(0)*sin(x.at(1)*w/c0)*sin(x.at(3)*w/c0)*sin(x.at(5)*w/c0))*x.at(4))/x.at(2);
        ABCD(1,0) = (-I*x.at(0)*x.at(2)*cos(x.at(1)*w/c0)*cos(x.at(3)*w/c0)*sin(x.at(5)*w/c0) + I*x.at(2)*x.at(2)*sin(x.at(1)*w/c0)*sin(x.at(3)*w/c0)*sin(x.at(5)*w/c0) + (-I*x.at(2)*cos(x.at(3)*w/c0)*cos(x.at(5)*w/c0)*sin(x.at(1)*w/c0) - I*x.at(0)*cos(x.at(1)*w/c0)*cos(x.at(5)*w/c0)*sin(x.at(3)*w/c0))*x.at(4))/(x.at(0)*x.at(2)*x.at(4));
        ABCD(1,1) = (x.at(0)*x.at(2)*cos(x.at(1)*w/c0)*cos(x.at(3)*w/c0)*cos(x.at(5)*w/c0) - x.at(2)*x.at(2)*cos(x.at(5)*w/c0)*sin(x.at(1)*w/c0)*sin(x.at(3)*w/c0) - (x.at(2)*cos(x.at(3)*w/c0)*sin(x.at(1)*w/c0)*sin(x.at(5)*w/c0) + x.at(0)*cos(x.at(1)*w/c0)*sin(x.at(3)*w/c0)*sin(x.at(5)*w/c0))*x.at(4))/(x.at(0)*x.at(2));
        return ABCD;
    }

    if (!topology.compare("03030"))
    {
        ABCD(0,0) = -(x.at(1)*w*w + x.at(3)*w*w)*x.at(0) + (x.at(1)*x.at(3)*x.at(0)*w*w*w*w - x.at(3)*w*w)*x.at(2) + 1;
        ABCD(0,1) = (-I*x.at(1)*x.at(0)*w*w*w + I*w)*x.at(2) + ((-I*x.at(1)*w*w*w - I*x.at(3)*w*w*w)*x.at(0) + (I*x.at(1)*x.at(3)*x.at(0)*w*w*w*w*w - I*x.at(3)*w*w*w)*x.at(2) + I*w)*x.at(4) + I*x.at(0)*w;
        ABCD(1,0) = -I*x.at(1)*x.at(3)*x.at(2)*w*w*w + I*x.at(1)*w + I*x.at(3)*w;
        ABCD(1,1) =  -x.at(1)*x.at(2)*w*w + (x.at(1)*x.at(3)*x.at(2)*w*w*w*w - x.at(1)*w*w - x.at(3)*w*w)*x.at(4) + 1;
        return ABCD;
    }

    if (!topology.compare("313202"))
    {
        ABCD(0,0) = (((x.at(2)*w*w + x.at(1)*w*w)*x.at(3) - 1)*x.at(5) + ((x.at(2)*w*w + x.at(1)*w*w)*x.at(3) - 1)*x.at(4) - x.at(3))/(x.at(1)*x.at(3)*x.at(5)*w*w);
        ABCD(0,1) = (((I*x.at(2)*w*w + I*x.at(1)*w*w)*x.at(3) - I)*x.at(4) - I*x.at(3))/(x.at(1)*x.at(3)*w);
        ABCD(1,0) = ((-I*x.at(0) - I*x.at(1))*x.at(3) + ((I*x.at(0)*x.at(2)*w*w + (I*x.at(0)*w*w + I*x.at(2)*w*w)*x.at(1))*x.at(3) - I*x.at(0) - I*x.at(1))*x.at(5) + ((I*x.at(0)*x.at(2)*w*w + (I*x.at(0)*w*w + I*x.at(2)*w*w)*x.at(1))*x.at(3) - I*x.at(0) - I*x.at(1))*x.at(4))/(x.at(1)*x.at(3)*x.at(5)*w);
        ABCD(1,1) = ((x.at(0) + x.at(1))*x.at(3) - ((x.at(0)*x.at(2)*w*w + (x.at(0)*w*w + x.at(2)*w*w)*x.at(1))*x.at(3) - x.at(0) - x.at(1))*x.at(4))/(x.at(1)*x.at(3));
        return ABCD;
    }

    if (!topology.compare("202313"))
    {
        ABCD(0,0) = ((x.at(5) + x.at(4))*x.at(2) - ((x.at(3)*x.at(5)*w*w + (x.at(3)*w*w + x.at(5)*w*w)*x.at(4))*x.at(2) - x.at(5) - x.at(4))*x.at(1))/(x.at(4)*x.at(2));
        ABCD(0,1) = (((I*x.at(3)*w*w + I*x.at(4)*w*w)*x.at(2) - I)*x.at(1) - I*x.at(2))/(x.at(4)*x.at(2)*w);
        ABCD(1,0) = ((-I*x.at(5) - I*x.at(4))*x.at(0) + ((I*x.at(3)*x.at(5)*w*w + (I*x.at(3)*w*w + I*x.at(5)*w*w)*x.at(4))*x.at(0) - I*x.at(5) - I*x.at(4))*x.at(2) + ((I*x.at(3)*x.at(5)*w*w + (I*x.at(3)*w*w + I*x.at(5)*w*w)*x.at(4))*x.at(2) - I*x.at(5) - I*x.at(4))*x.at(1))/(x.at(4)*x.at(0)*x.at(2)*w);
        ABCD(1,1) = (((x.at(3)*w*w + x.at(4)*w*w)*x.at(0) - 1)*x.at(2) + ((x.at(3)*w*w + x.at(4)*w*w)*x.at(2) - 1)*x.at(1) - x.at(0))/(x.at(4)*x.at(0)*x.at(2)*w*w);
        return ABCD;
    }

    if (!topology.compare("020313"))
    {
        ABCD(0,0) = ((x.at(5) + x.at(4))*x.at(1) - ((x.at(3)*x.at(5)*w*w + (x.at(3)*w*w + x.at(5)*w*w)*x.at(4))*x.at(1) - x.at(5) - x.at(4))*x.at(0) - ((x.at(3)*x.at(5)*w*w + (x.at(3)*w*w + x.at(5)*w*w)*x.at(4))*x.at(1) + (x.at(3)*x.at(5)*w*w + (x.at(3)*w*w + x.at(5)*w*w)*x.at(4))*x.at(0))*x.at(2))/(x.at(4)*x.at(1));
        ABCD(0,1) = (((I*x.at(3)*w*w + I*x.at(4)*w*w)*x.at(1) - I)*x.at(0) + ((I*x.at(3)*w*w + I*x.at(4)*w*w)*x.at(1) + (I*x.at(3)*w*w + I*x.at(4)*w*w)*x.at(0))*x.at(2) - I*x.at(1))/(x.at(4)*x.at(1)*w);
        ABCD(1,0) = ((I*x.at(3)*x.at(5)*w*w + (I*x.at(3)*w*w + I*x.at(5)*w*w)*x.at(4))*x.at(1) + (I*x.at(3)*x.at(5)*w*w + (I*x.at(3)*w*w + I*x.at(5)*w*w)*x.at(4))*x.at(2) - I*x.at(5) - I*x.at(4))/(x.at(4)*x.at(1)*w);
        ABCD(1,1) = ((x.at(3)*w*w + x.at(4)*w*w)*x.at(1) + (x.at(3)*w*w + x.at(4)*w*w)*x.at(2) - 1)/(x.at(4)*x.at(1)*w*w);
        return ABCD;
    }

    if (!topology.compare("313020"))
    {
        ABCD(0,0) = ((x.at(2)*w*w + x.at(1)*w*w)*x.at(4) + (x.at(2)*w*w + x.at(1)*w*w)*x.at(3) - 1)/(x.at(1)*x.at(4)*w*w);
        ABCD(0,1) = ((I*x.at(2)*w*w + I*x.at(1)*w*w)*x.at(4)*x.at(3) + ((I*x.at(2)*w*w + I*x.at(1)*w*w)*x.at(4) + (I*x.at(2)*w*w + I*x.at(1)*w*w)*x.at(3) - I)*x.at(5) - I*x.at(4))/(x.at(1)*x.at(4)*w);
        ABCD(1,0) = ((I*x.at(0)*x.at(2)*w*w + (I*x.at(0)*w*w + I*x.at(2)*w*w)*x.at(1))*x.at(4) + (I*x.at(0)*x.at(2)*w*w + (I*x.at(0)*w*w + I*x.at(2)*w*w)*x.at(1))*x.at(3) - I*x.at(0) - I*x.at(1))/(x.at(1)*x.at(4)*w);
        ABCD(1,1) = -((x.at(0)*x.at(2)*w*w + (x.at(0)*w*w + x.at(2)*w*w)*x.at(1))*x.at(4)*x.at(3) - (x.at(0) + x.at(1))*x.at(4) + ((x.at(0)*x.at(2)*w*w + (x.at(0)*w*w + x.at(2)*w*w)*x.at(1))*x.at(4) + (x.at(0)*x.at(2)*w*w + (x.at(0)*w*w + x.at(2)*w*w)*x.at(1))*x.at(3) - x.at(0) - x.at(1))*x.at(5))/(x.at(1)*x.at(4));
        return ABCD;
    }
/*
    if (!topology.compare("43434"))// Lowpass TL-C ladder
    {
        ABCD(0,0) = -(x.at(0)*x.at(3)*cos(x.at(4)*w/c0)*sin(x.at(1)*w/c0)*sin(x.at(7)*w/c0) + (x.at(2)*x.at(0)*w*sin(x.at(1)*w/c0)*sin(x.at(4)*w/c0)*sin(x.at(7)*w/c0) + cos(x.at(1)*w/c0)*sin(x.at(4)*w/c0)*sin(x.at(7)*w/c0))*x.at(3)*x.at(3) + (x.at(0)*cos(x.at(7)*w/c0)*sin(x.at(1)*w/c0)*sin(x.at(4)*w/c0) - (x.at(2)*x.at(5)*x.at(0)*w*w*cos(x.at(7)*w/c0)*sin(x.at(1)*w/c0)*sin(x.at(4)*w/c0) + x.at(5)*w*cos(x.at(1)*w/c0)*cos(x.at(7)*w/c0)*sin(x.at(4)*w/c0))*x.at(3)*x.at(3) - ((x.at(2)*w + x.at(5)*w)*x.at(0)*cos(x.at(4)*w/c0)*cos(x.at(7)*w/c0)*sin(x.at(1)*w/c0) + cos(x.at(1)*w/c0)*cos(x.at(4)*w/c0)*cos(x.at(7)*w/c0))*x.at(3))*x.at(6))/(x.at(3)*x.at(6));
        ABCD(0,1) = (-I*x.at(0)*x.at(3)*cos(x.at(4)*w/c0)*cos(x.at(7)*w/c0)*sin(x.at(1)*w/c0) - (I*x.at(2)*x.at(0)*w*cos(x.at(7)*w/c0)*sin(x.at(1)*w/c0)*sin(x.at(4)*w/c0) + I*cos(x.at(1)*w/c0)*cos(x.at(7)*w/c0)*sin(x.at(4)*w/c0))*x.at(3)*x.at(3) + (I*x.at(0)*sin(x.at(1)*w/c0)*sin(x.at(4)*w/c0)*sin(x.at(7)*w/c0) - (I*x.at(2)*x.at(5)*x.at(0)*w*w*sin(x.at(1)*w/c0)*sin(x.at(4)*w/c0)*sin(x.at(7)*w/c0) + I*x.at(5)*w*cos(x.at(1)*w/c0)*sin(x.at(4)*w/c0)*sin(x.at(7)*w/c0))*x.at(3)*x.at(3) - ((I*x.at(2)*w + I*x.at(5)*w)*x.at(0)*cos(x.at(4)*w/c0)*sin(x.at(1)*w/c0)*sin(x.at(7)*w/c0) + I*cos(x.at(1)*w/c0)*cos(x.at(4)*w/c0)*sin(x.at(7)*w/c0))*x.at(3))*x.at(6))/x.at(3);
        ABCD(1,0) = (-I*x.at(0)*x.at(3)*cos(x.at(1)*w/c0)*cos(x.at(4)*w/c0)*sin(x.at(7)*w/c0) + (-I*x.at(2)*x.at(0)*w*cos(x.at(1)*w/c0)*sin(x.at(4)*w/c0)*sin(x.at(7)*w/c0) + I*sin(x.at(1)*w/c0)*sin(x.at(4)*w/c0)*sin(x.at(7)*w/c0))*x.at(3)*x.at(3) + (-I*x.at(0)*cos(x.at(1)*w/c0)*cos(x.at(7)*w/c0)*sin(x.at(4)*w/c0) + (I*x.at(2)*x.at(5)*x.at(0)*w*w*cos(x.at(1)*w/c0)*cos(x.at(7)*w/c0)*sin(x.at(4)*w/c0) - I*x.at(5)*w*cos(x.at(7)*w/c0)*sin(x.at(1)*w/c0)*sin(x.at(4)*w/c0))*x.at(3)*x.at(3) + ((I*x.at(2)*w + I*x.at(5)*w)*x.at(0)*cos(x.at(1)*w/c0)*cos(x.at(4)*w/c0)*cos(x.at(7)*w/c0) - I*cos(x.at(4)*w/c0)*cos(x.at(7)*w/c0)*sin(x.at(1)*w/c0))*x.at(3))*x.at(6))/(x.at(0)*x.at(3)*x.at(6));
        ABCD(1,1) = (x.at(0)*x.at(3)*cos(x.at(1)*w/c0)*cos(x.at(4)*w/c0)*cos(x.at(7)*w/c0) + (x.at(2)*x.at(0)*w*cos(x.at(1)*w/c0)*cos(x.at(7)*w/c0)*sin(x.at(4)*w/c0) - cos(x.at(7)*w/c0)*sin(x.at(1)*w/c0)*sin(x.at(4)*w/c0))*x.at(3)*x.at(3) - (x.at(0)*cos(x.at(1)*w/c0)*sin(x.at(4)*w/c0)*sin(x.at(7)*w/c0) - (x.at(2)*x.at(5)*x.at(0)*w*w*cos(x.at(1)*w/c0)*sin(x.at(4)*w/c0)*sin(x.at(7)*w/c0) - x.at(5)*w*sin(x.at(1)*w/c0)*sin(x.at(4)*w/c0)*sin(x.at(7)*w/c0))*x.at(3)*x.at(3) - ((x.at(2)*w + x.at(5)*w)*x.at(0)*cos(x.at(1)*w/c0)*cos(x.at(4)*w/c0)*sin(x.at(7)*w/c0) - cos(x.at(4)*w/c0)*sin(x.at(1)*w/c0)*sin(x.at(7)*w/c0))*x.at(3))*x.at(6))/(x.at(0)*x.at(3));
        return ABCD;
    }
*/
    if (!topology.compare("01230123")) //BPP3
    {
        ABCD(0,0) = -((x.at(3)*w*w + x.at(1)*w*w + x.at(5)*w*w)*x.at(2) + (x.at(7)*w*w + x.at(5)*w*w - (x.at(3)*x.at(7)*w*w*w*w + x.at(7)*x.at(1)*w*w*w*w + (x.at(3)*w*w*w*w + x.at(7)*w*w*w*w + x.at(1)*w*w*w*w)*x.at(5))*x.at(2))*x.at(6) + (x.at(1)*w*w - (x.at(3)*x.at(1)*w*w*w*w + x.at(1)*x.at(5)*w*w*w*w)*x.at(2) - (x.at(7)*x.at(1)*w*w*w*w + x.at(1)*x.at(5)*w*w*w*w - (x.at(3)*x.at(7)*x.at(1)*w*w*w*w*w*w + (x.at(3)*w*w*w*w*w*w + x.at(7)*w*w*w*w*w*w)*x.at(1)*x.at(5))*x.at(2))*x.at(6))*x.at(0) - ((x.at(3)*w*w*w*w + x.at(1)*w*w*w*w)*x.at(5)*x.at(2) - x.at(5)*w*w + (x.at(7)*x.at(5)*w*w*w*w - (x.at(3)*x.at(7)*w*w*w*w*w*w + x.at(7)*x.at(1)*w*w*w*w*w*w)*x.at(5)*x.at(2))*x.at(6) - (x.at(3)*x.at(1)*x.at(5)*x.at(2)*w*w*w*w*w*w - x.at(1)*x.at(5)*w*w*w*w - (x.at(3)*x.at(7)*x.at(1)*x.at(5)*x.at(2)*w*w*w*w*w*w*w*w - x.at(7)*x.at(1)*x.at(5)*w*w*w*w*w*w)*x.at(6))*x.at(0))*x.at(4) - 1)/(x.at(1)*x.at(5)*x.at(2)*x.at(6)*w*w*w*w);
        ABCD(0,1) = -((I*x.at(3)*w*w + I*x.at(1)*w*w + I*x.at(5)*w*w)*x.at(2) + (I*x.at(1)*w*w + (-I*x.at(3)*x.at(1)*w*w*w*w - I*x.at(1)*x.at(5)*w*w*w*w)*x.at(2))*x.at(0) + ((-I*x.at(3)*w*w*w*w - I*x.at(1)*w*w*w*w)*x.at(5)*x.at(2) + I*x.at(5)*w*w + (I*x.at(3)*x.at(1)*x.at(5)*x.at(2)*w*w*w*w*w*w - I*x.at(1)*x.at(5)*w*w*w*w)*x.at(0))*x.at(4) - I)/(x.at(1)*x.at(5)*x.at(2)*w*w*w);
        ABCD(1,0) = -((I*x.at(3)*w*w + I*x.at(5)*w*w)*x.at(2) + (I*x.at(7)*w*w + I*x.at(5)*w*w + (-I*x.at(3)*x.at(7)*w*w*w*w + (-I*x.at(3)*w*w*w*w - I*x.at(7)*w*w*w*w)*x.at(5))*x.at(2))*x.at(6) + (-I*x.at(3)*x.at(5)*x.at(2)*w*w*w*w + I*x.at(5)*w*w + (I*x.at(3)*x.at(7)*x.at(5)*x.at(2)*w*w*w*w*w*w - I*x.at(7)*x.at(5)*w*w*w*w)*x.at(6))*x.at(4) - I)/(x.at(5)*x.at(2)*x.at(6)*w*w*w);
        ABCD(1,1) = ((x.at(3)*w*w + x.at(5)*w*w)*x.at(2) - (x.at(3)*x.at(5)*x.at(2)*w*w*w*w - x.at(5)*w*w)*x.at(4) - 1)/(x.at(5)*x.at(2)*w*w);
        return ABCD;
    }

    if (!topology.compare("23012301")) //BPS3
    {
        ABCD(0,0) = ((x.at(5)*w*w + x.at(3)*w*w)*x.at(4) - (x.at(5)*x.at(3)*x.at(4)*w*w*w*w - x.at(3)*w*w)*x.at(2) - 1)/(x.at(3)*x.at(4)*w*w);
        ABCD(0,1) = -((I*x.at(5)*w*w + I*x.at(3)*w*w + I*x.at(7)*w*w)*x.at(4) + (I*x.at(3)*w*w + (-I*x.at(5)*x.at(3)*w*w*w*w - I*x.at(3)*x.at(7)*w*w*w*w)*x.at(4))*x.at(2) + ((-I*x.at(5)*w*w*w*w - I*x.at(3)*w*w*w*w)*x.at(7)*x.at(4) + I*x.at(7)*w*w + (I*x.at(5)*x.at(3)*x.at(7)*x.at(4)*w*w*w*w*w*w - I*x.at(3)*x.at(7)*w*w*w*w)*x.at(2))*x.at(6) - I)/(x.at(3)*x.at(7)*x.at(4)*w*w*w);
        ABCD(1,0) = -((I*x.at(1)*w*w + I*x.at(3)*w*w)*x.at(0) + (I*x.at(5)*w*w + I*x.at(3)*w*w + (-I*x.at(1)*x.at(5)*w*w*w*w + (-I*x.at(1)*w*w*w*w - I*x.at(5)*w*w*w*w)*x.at(3))*x.at(0))*x.at(4) + (-I*x.at(1)*x.at(3)*x.at(0)*w*w*w*w + I*x.at(3)*w*w + (I*x.at(1)*x.at(5)*x.at(3)*x.at(0)*w*w*w*w*w*w - I*x.at(5)*x.at(3)*w*w*w*w)*x.at(4))*x.at(2) - I)/(x.at(3)*x.at(0)*x.at(4)*w*w*w);
        ABCD(1,1) = -((x.at(1)*w*w + x.at(3)*w*w)*x.at(0) + (x.at(5)*w*w + x.at(3)*w*w + x.at(7)*w*w - (x.at(1)*x.at(5)*w*w*w*w + (x.at(1)*w*w*w*w + x.at(5)*w*w*w*w)*x.at(3) + (x.at(1)*w*w*w*w + x.at(3)*w*w*w*w)*x.at(7))*x.at(0))*x.at(4) - (x.at(1)*x.at(3)*x.at(0)*w*w*w*w - x.at(3)*w*w + (x.at(5)*x.at(3)*w*w*w*w + x.at(3)*x.at(7)*w*w*w*w - (x.at(1)*x.at(5)*x.at(3)*w*w*w*w*w*w + x.at(1)*x.at(3)*x.at(7)*w*w*w*w*w*w)*x.at(0))*x.at(4))*x.at(2) - ((x.at(1)*w*w*w*w + x.at(3)*w*w*w*w)*x.at(7)*x.at(0) - x.at(7)*w*w - ((x.at(1)*x.at(5)*w*w*w*w*w*w + (x.at(1)*w*w*w*w*w*w + x.at(5)*w*w*w*w*w*w)*x.at(3))*x.at(7)*x.at(0) - (x.at(5)*w*w*w*w + x.at(3)*w*w*w*w)*x.at(7))*x.at(4) - (x.at(1)*x.at(3)*x.at(7)*x.at(0)*w*w*w*w*w*w - x.at(3)*x.at(7)*w*w*w*w - (x.at(1)*x.at(5)*x.at(3)*x.at(7)*x.at(0)*w*w*w*w*w*w*w*w - x.at(5)*x.at(3)*x.at(7)*w*w*w*w*w*w)*x.at(4))*x.at(2))*x.at(6) - 1)/(x.at(3)*x.at(7)*x.at(0)*x.at(4)*w*w*w*w);
        return ABCD;
    }

    if (!topology.compare("12031203")) //LPHHP
    {
        ABCD(0,0) = -((x.at(3)*w*w + x.at(0)*w*w + x.at(4)*w*w)*x.at(1) + (x.at(7)*w*w + x.at(4)*w*w - (x.at(3)*x.at(7)*w*w*w*w + x.at(7)*x.at(0)*w*w*w*w + (x.at(3)*w*w*w*w + x.at(7)*w*w*w*w + x.at(0)*w*w*w*w)*x.at(4))*x.at(1))*x.at(5) + (x.at(3)*w*w + x.at(4)*w*w - (x.at(3)*x.at(0)*w*w*w*w + x.at(0)*x.at(4)*w*w*w*w)*x.at(1) - (x.at(3)*x.at(7)*w*w*w*w + (x.at(3)*w*w*w*w + x.at(7)*w*w*w*w)*x.at(4) - (x.at(3)*x.at(7)*x.at(0)*w*w*w*w*w*w + (x.at(3)*w*w*w*w*w*w + x.at(7)*w*w*w*w*w*w)*x.at(0)*x.at(4))*x.at(1))*x.at(5))*x.at(2) + (x.at(7)*w*w - (x.at(3)*x.at(7)*w*w*w*w + x.at(7)*x.at(0)*w*w*w*w + x.at(7)*x.at(4)*w*w*w*w)*x.at(1) - (x.at(7)*x.at(4)*w*w*w*w - (x.at(3)*x.at(7)*w*w*w*w*w*w + x.at(7)*x.at(0)*w*w*w*w*w*w)*x.at(4)*x.at(1))*x.at(5) - (x.at(3)*x.at(7)*w*w*w*w + x.at(7)*x.at(4)*w*w*w*w - (x.at(3)*x.at(7)*x.at(0)*w*w*w*w*w*w + x.at(7)*x.at(0)*x.at(4)*w*w*w*w*w*w)*x.at(1) + (x.at(3)*x.at(7)*x.at(0)*x.at(4)*x.at(1)*w*w*w*w*w*w*w*w - x.at(3)*x.at(7)*x.at(4)*w*w*w*w*w*w)*x.at(5))*x.at(2))*x.at(6) - 1)/(x.at(0)*x.at(4)*x.at(1)*x.at(5)*w*w*w*w);
        ABCD(0,1) = -((I*x.at(3)*w*w + I*x.at(4)*w*w + (-I*x.at(3)*x.at(0)*w*w*w*w - I*x.at(0)*x.at(4)*w*w*w*w)*x.at(1))*x.at(5)*x.at(2) + ((I*x.at(3)*w*w + I*x.at(0)*w*w + I*x.at(4)*w*w)*x.at(1) - I)*x.at(5) + ((I*x.at(3)*w*w + I*x.at(0)*w*w + I*x.at(4)*w*w)*x.at(1) + ((-I*x.at(3)*w*w*w*w - I*x.at(0)*w*w*w*w)*x.at(4)*x.at(1) + I*x.at(4)*w*w)*x.at(5) + (I*x.at(3)*w*w + I*x.at(4)*w*w + (-I*x.at(3)*x.at(0)*w*w*w*w - I*x.at(0)*x.at(4)*w*w*w*w)*x.at(1) + (I*x.at(3)*x.at(0)*x.at(4)*x.at(1)*w*w*w*w*w*w - I*x.at(3)*x.at(4)*w*w*w*w)*x.at(5))*x.at(2) - I)*x.at(6))/(x.at(0)*x.at(4)*x.at(1)*x.at(5)*w*w*w);
        ABCD(1,0) = -((I*x.at(3)*w*w + I*x.at(4)*w*w)*x.at(1) + (I*x.at(7)*w*w + I*x.at(4)*w*w + (-I*x.at(3)*x.at(7)*w*w*w*w + (-I*x.at(3)*w*w*w*w - I*x.at(7)*w*w*w*w)*x.at(4))*x.at(1))*x.at(5) + (I*x.at(3)*w*w + I*x.at(4)*w*w + (-I*x.at(3)*x.at(7)*w*w*w*w + (-I*x.at(3)*w*w*w*w - I*x.at(7)*w*w*w*w)*x.at(4))*x.at(5))*x.at(2) + (I*x.at(7)*w*w + (-I*x.at(3)*x.at(7)*w*w*w*w - I*x.at(7)*x.at(4)*w*w*w*w)*x.at(1) + (I*x.at(3)*x.at(7)*x.at(4)*x.at(1)*w*w*w*w*w*w - I*x.at(7)*x.at(4)*w*w*w*w)*x.at(5) + (I*x.at(3)*x.at(7)*x.at(4)*x.at(5)*w*w*w*w*w*w - I*x.at(3)*x.at(7)*w*w*w*w - I*x.at(7)*x.at(4)*w*w*w*w)*x.at(2))*x.at(6) - I)/(x.at(4)*x.at(1)*x.at(5)*w*w*w);
        ABCD(1,1) = ((x.at(3)*w*w + x.at(4)*w*w)*x.at(5)*x.at(2) + ((x.at(3)*w*w + x.at(4)*w*w)*x.at(1) - 1)*x.at(5) + ((x.at(3)*w*w + x.at(4)*w*w)*x.at(1) - (x.at(3)*x.at(4)*x.at(1)*w*w*w*w - x.at(4)*w*w)*x.at(5) - (x.at(3)*x.at(4)*x.at(5)*w*w*w*w - x.at(3)*w*w - x.at(4)*w*w)*x.at(2) - 1)*x.at(6))/(x.at(4)*x.at(1)*x.at(5)*w*w);
        return ABCD;
    }

    if (!topology.compare("21302130")) //LPHPS
    {
        ABCD(0,0) = ((x.at(2)*x.at(6)*w*w + x.at(6)*x.at(1)*w*w + (x.at(2)*w*w + x.at(6)*w*w + x.at(1)*w*w)*x.at(5))*x.at(4) + (x.at(2)*x.at(6)*w*w + x.at(6)*x.at(1)*w*w - (x.at(2)*x.at(6)*w*w*w*w + x.at(6)*x.at(1)*w*w*w*w)*x.at(5)*x.at(4) + (x.at(2)*w*w + x.at(1)*w*w)*x.at(5))*x.at(3) - x.at(6) - x.at(5))/(x.at(1)*x.at(5)*x.at(4)*w*w);
        ABCD(0,1) = -((I*x.at(2)*w*w + I*x.at(1)*w*w + I*x.at(5)*w*w)*x.at(4) + ((-I*x.at(2)*w*w*w*w - I*x.at(1)*w*w*w*w)*x.at(5)*x.at(4) + I*x.at(2)*w*w + I*x.at(1)*w*w)*x.at(3) + (I*x.at(6)*w*w + I*x.at(5)*w*w + (-I*x.at(2)*x.at(6)*w*w*w*w - I*x.at(6)*x.at(1)*w*w*w*w + (-I*x.at(2)*w*w*w*w - I*x.at(6)*w*w*w*w - I*x.at(1)*w*w*w*w)*x.at(5))*x.at(4) + (-I*x.at(2)*x.at(6)*w*w*w*w - I*x.at(6)*x.at(1)*w*w*w*w + (I*x.at(2)*x.at(6)*w*w*w*w*w*w + I*x.at(6)*x.at(1)*w*w*w*w*w*w)*x.at(5)*x.at(4) + (-I*x.at(2)*w*w*w*w - I*x.at(1)*w*w*w*w)*x.at(5))*x.at(3))*x.at(7) - I)/(x.at(1)*x.at(5)*x.at(4)*w*w*w);
        ABCD(1,0) = -((I*x.at(6)*x.at(1)*w*w + I*x.at(1)*x.at(5)*w*w)*x.at(0) + (I*x.at(2)*x.at(6)*w*w + I*x.at(6)*x.at(1)*w*w + (I*x.at(2)*w*w + I*x.at(6)*w*w + I*x.at(1)*w*w)*x.at(5) + (-I*x.at(2)*x.at(6)*x.at(1)*w*w*w*w + (-I*x.at(2)*w*w*w*w - I*x.at(6)*w*w*w*w)*x.at(1)*x.at(5))*x.at(0))*x.at(4) + (I*x.at(2)*x.at(6)*w*w + I*x.at(6)*x.at(1)*w*w + (I*x.at(2)*w*w + I*x.at(1)*w*w)*x.at(5) + (-I*x.at(2)*x.at(6)*x.at(1)*w*w*w*w - I*x.at(2)*x.at(1)*x.at(5)*w*w*w*w)*x.at(0) + (I*x.at(2)*x.at(6)*x.at(1)*x.at(5)*x.at(0)*w*w*w*w*w*w + (-I*x.at(2)*x.at(6)*w*w*w*w - I*x.at(6)*x.at(1)*w*w*w*w)*x.at(5))*x.at(4))*x.at(3) - I*x.at(6) - I*x.at(5))/(x.at(1)*x.at(5)*x.at(0)*x.at(4)*w*w*w);
        ABCD(1,1) = -(x.at(1)*x.at(0)*w*w + (x.at(2)*w*w + x.at(1)*w*w + x.at(5)*w*w - (x.at(2)*x.at(1)*w*w*w*w + x.at(1)*x.at(5)*w*w*w*w)*x.at(0))*x.at(4) - (x.at(2)*x.at(1)*x.at(0)*w*w*w*w - x.at(2)*w*w - x.at(1)*w*w - (x.at(2)*x.at(1)*x.at(5)*x.at(0)*w*w*w*w*w*w - (x.at(2)*w*w*w*w + x.at(1)*w*w*w*w)*x.at(5))*x.at(4))*x.at(3) + (x.at(6)*w*w + x.at(5)*w*w - (x.at(6)*x.at(1)*w*w*w*w + x.at(1)*x.at(5)*w*w*w*w)*x.at(0) - (x.at(2)*x.at(6)*w*w*w*w + x.at(6)*x.at(1)*w*w*w*w + (x.at(2)*w*w*w*w + x.at(6)*w*w*w*w + x.at(1)*w*w*w*w)*x.at(5) - (x.at(2)*x.at(6)*x.at(1)*w*w*w*w*w*w + (x.at(2)*w*w*w*w*w*w + x.at(6)*w*w*w*w*w*w)*x.at(1)*x.at(5))*x.at(0))*x.at(4) - (x.at(2)*x.at(6)*w*w*w*w + x.at(6)*x.at(1)*w*w*w*w + (x.at(2)*w*w*w*w + x.at(1)*w*w*w*w)*x.at(5) - (x.at(2)*x.at(6)*x.at(1)*w*w*w*w*w*w + x.at(2)*x.at(1)*x.at(5)*w*w*w*w*w*w)*x.at(0) + (x.at(2)*x.at(6)*x.at(1)*x.at(5)*x.at(0)*w*w*w*w*w*w*w*w - (x.at(2)*x.at(6)*w*w*w*w*w*w + x.at(6)*x.at(1)*w*w*w*w*w*w)*x.at(5))*x.at(4))*x.at(3))*x.at(7) - 1)/(x.at(1)*x.at(5)*x.at(0)*x.at(4)*w*w*w*w);
        return ABCD;
    }

    if (!topology.compare("03120312")) //HPLPP
    {
        ABCD(0,0) = -((x.at(5)*w*w + x.at(2)*w*w + x.at(6)*w*w)*x.at(3) - ((x.at(5)*w*w*w*w + x.at(2)*w*w*w*w)*x.at(6)*x.at(3) - x.at(6)*w*w)*x.at(7) + (x.at(1)*w*w + x.at(2)*w*w - (x.at(1)*x.at(5)*w*w*w*w + (x.at(1)*w*w*w*w + x.at(5)*w*w*w*w)*x.at(2) + (x.at(1)*w*w*w*w + x.at(2)*w*w*w*w)*x.at(6))*x.at(3) + ((x.at(1)*x.at(5)*w*w*w*w*w*w + (x.at(1)*w*w*w*w*w*w + x.at(5)*w*w*w*w*w*w)*x.at(2))*x.at(6)*x.at(3) - (x.at(1)*w*w*w*w + x.at(2)*w*w*w*w)*x.at(6))*x.at(7))*x.at(0) + (x.at(5)*w*w + x.at(6)*w*w - (x.at(5)*x.at(2)*w*w*w*w + x.at(2)*x.at(6)*w*w*w*w)*x.at(3) + (x.at(5)*x.at(2)*x.at(6)*x.at(3)*w*w*w*w*w*w - x.at(5)*x.at(6)*w*w*w*w)*x.at(7) - (x.at(1)*x.at(5)*w*w*w*w + x.at(5)*x.at(2)*w*w*w*w + (x.at(1)*w*w*w*w + x.at(2)*w*w*w*w)*x.at(6) - (x.at(1)*x.at(5)*x.at(2)*w*w*w*w*w*w + x.at(1)*x.at(2)*x.at(6)*w*w*w*w*w*w)*x.at(3) + (x.at(1)*x.at(5)*x.at(2)*x.at(6)*x.at(3)*w*w*w*w*w*w*w*w - (x.at(1)*x.at(5)*w*w*w*w*w*w + x.at(5)*x.at(2)*w*w*w*w*w*w)*x.at(6))*x.at(7))*x.at(0))*x.at(4) - 1)/(x.at(2)*x.at(6)*x.at(3)*x.at(7)*w*w*w*w);
        ABCD(0,1) = -((I*x.at(5)*w*w + I*x.at(2)*w*w + I*x.at(6)*w*w)*x.at(3) + (I*x.at(1)*w*w + I*x.at(2)*w*w + (-I*x.at(1)*x.at(5)*w*w*w*w + (-I*x.at(1)*w*w*w*w - I*x.at(5)*w*w*w*w)*x.at(2) + (-I*x.at(1)*w*w*w*w - I*x.at(2)*w*w*w*w)*x.at(6))*x.at(3))*x.at(0) + (I*x.at(5)*w*w + I*x.at(6)*w*w + (-I*x.at(5)*x.at(2)*w*w*w*w - I*x.at(2)*x.at(6)*w*w*w*w)*x.at(3) + (-I*x.at(1)*x.at(5)*w*w*w*w - I*x.at(5)*x.at(2)*w*w*w*w + (-I*x.at(1)*w*w*w*w - I*x.at(2)*w*w*w*w)*x.at(6) + (I*x.at(1)*x.at(5)*x.at(2)*w*w*w*w*w*w + I*x.at(1)*x.at(2)*x.at(6)*w*w*w*w*w*w)*x.at(3))*x.at(0))*x.at(4) - I)/(x.at(2)*x.at(6)*x.at(3)*w*w*w);
        ABCD(1,0) = -((I*x.at(1)*x.at(5)*w*w + (I*x.at(1)*w*w + I*x.at(5)*w*w)*x.at(2) + (I*x.at(1)*w*w + I*x.at(2)*w*w)*x.at(6))*x.at(3) + ((-I*x.at(1)*x.at(5)*w*w*w*w + (-I*x.at(1)*w*w*w*w - I*x.at(5)*w*w*w*w)*x.at(2))*x.at(6)*x.at(3) + (I*x.at(1)*w*w + I*x.at(2)*w*w)*x.at(6))*x.at(7) + (I*x.at(1)*x.at(5)*w*w + I*x.at(5)*x.at(2)*w*w + (I*x.at(1)*w*w + I*x.at(2)*w*w)*x.at(6) + (-I*x.at(1)*x.at(5)*x.at(2)*w*w*w*w - I*x.at(1)*x.at(2)*x.at(6)*w*w*w*w)*x.at(3) + (I*x.at(1)*x.at(5)*x.at(2)*x.at(6)*x.at(3)*w*w*w*w*w*w + (-I*x.at(1)*x.at(5)*w*w*w*w - I*x.at(5)*x.at(2)*w*w*w*w)*x.at(6))*x.at(7))*x.at(4) - I*x.at(1) - I*x.at(2))/(x.at(2)*x.at(6)*x.at(3)*x.at(7)*w*w*w);
        ABCD(1,1) = ((x.at(1)*x.at(5)*w*w + (x.at(1)*w*w + x.at(5)*w*w)*x.at(2) + (x.at(1)*w*w + x.at(2)*w*w)*x.at(6))*x.at(3) + (x.at(1)*x.at(5)*w*w + x.at(5)*x.at(2)*w*w + (x.at(1)*w*w + x.at(2)*w*w)*x.at(6) - (x.at(1)*x.at(5)*x.at(2)*w*w*w*w + x.at(1)*x.at(2)*x.at(6)*w*w*w*w)*x.at(3))*x.at(4) - x.at(1) - x.at(2))/(x.at(2)*x.at(6)*x.at(3)*w*w);
        return ABCD;
    }

    if (!topology.compare("30213021")) //HPLPS
    {
        ABCD(0,0) = ((x.at(4)*w*w + x.at(3)*w*w)*x.at(2)*x.at(6) + (x.at(3)*x.at(2)*w*w - (x.at(4)*x.at(3)*x.at(2)*w*w*w*w - x.at(4)*w*w - x.at(3)*w*w)*x.at(6) - 1)*x.at(1) + ((x.at(4)*w*w + x.at(3)*w*w)*x.at(2) - (x.at(4)*x.at(3)*x.at(2)*w*w*w*w - x.at(4)*w*w - x.at(3)*w*w)*x.at(1))*x.at(5) - x.at(2))/(x.at(3)*x.at(2)*x.at(6)*w*w);
        ABCD(0,1) = -((I*x.at(4)*w*w + I*x.at(3)*w*w + I*x.at(7)*w*w)*x.at(2)*x.at(6) + (I*x.at(3)*x.at(2)*w*w + (I*x.at(4)*w*w + I*x.at(3)*w*w + I*x.at(7)*w*w + (-I*x.at(4)*x.at(3)*w*w*w*w - I*x.at(3)*x.at(7)*w*w*w*w)*x.at(2))*x.at(6) - I)*x.at(1) + ((-I*x.at(4)*w*w*w*w - I*x.at(3)*w*w*w*w)*x.at(7)*x.at(2)*x.at(6) + (I*x.at(4)*w*w + I*x.at(3)*w*w)*x.at(2) + (-I*x.at(4)*x.at(3)*x.at(2)*w*w*w*w + I*x.at(4)*w*w + I*x.at(3)*w*w + (I*x.at(4)*x.at(3)*x.at(7)*x.at(2)*w*w*w*w*w*w + (-I*x.at(4)*w*w*w*w - I*x.at(3)*w*w*w*w)*x.at(7))*x.at(6))*x.at(1))*x.at(5) - I*x.at(2))/(x.at(3)*x.at(7)*x.at(2)*x.at(6)*w*w*w);
        ABCD(1,0) = -((I*x.at(0)*w*w + I*x.at(3)*w*w)*x.at(2) + (I*x.at(4)*w*w + I*x.at(3)*w*w + (-I*x.at(0)*x.at(4)*w*w*w*w + (-I*x.at(0)*w*w*w*w - I*x.at(4)*w*w*w*w)*x.at(3))*x.at(2))*x.at(6) + (-I*x.at(0)*x.at(3)*x.at(2)*w*w*w*w + I*x.at(0)*w*w + (I*x.at(0)*x.at(4)*x.at(3)*x.at(2)*w*w*w*w*w*w - I*x.at(0)*x.at(4)*w*w*w*w - I*x.at(0)*x.at(3)*w*w*w*w)*x.at(6))*x.at(1) + (I*x.at(4)*w*w + I*x.at(3)*w*w + (-I*x.at(0)*x.at(4)*w*w*w*w + (-I*x.at(0)*w*w*w*w - I*x.at(4)*w*w*w*w)*x.at(3))*x.at(2) + (I*x.at(0)*x.at(4)*x.at(3)*x.at(2)*w*w*w*w*w*w - I*x.at(0)*x.at(4)*w*w*w*w - I*x.at(0)*x.at(3)*w*w*w*w)*x.at(1))*x.at(5) - I)/(x.at(3)*x.at(2)*x.at(6)*w*w*w);
        ABCD(1,1) = -((x.at(0)*w*w + x.at(3)*w*w)*x.at(2) + (x.at(4)*w*w + x.at(3)*w*w + x.at(7)*w*w - (x.at(0)*x.at(4)*w*w*w*w + (x.at(0)*w*w*w*w + x.at(4)*w*w*w*w)*x.at(3) + (x.at(0)*w*w*w*w + x.at(3)*w*w*w*w)*x.at(7))*x.at(2))*x.at(6) - (x.at(0)*x.at(3)*x.at(2)*w*w*w*w - x.at(0)*w*w + (x.at(0)*x.at(4)*w*w*w*w + x.at(0)*x.at(3)*w*w*w*w + x.at(0)*x.at(7)*w*w*w*w - (x.at(0)*x.at(4)*x.at(3)*w*w*w*w*w*w + x.at(0)*x.at(3)*x.at(7)*w*w*w*w*w*w)*x.at(2))*x.at(6))*x.at(1) + (x.at(4)*w*w + x.at(3)*w*w - (x.at(0)*x.at(4)*w*w*w*w + (x.at(0)*w*w*w*w + x.at(4)*w*w*w*w)*x.at(3))*x.at(2) + ((x.at(0)*x.at(4)*w*w*w*w*w*w + (x.at(0)*w*w*w*w*w*w + x.at(4)*w*w*w*w*w*w)*x.at(3))*x.at(7)*x.at(2) - (x.at(4)*w*w*w*w + x.at(3)*w*w*w*w)*x.at(7))*x.at(6) + (x.at(0)*x.at(4)*x.at(3)*x.at(2)*w*w*w*w*w*w - x.at(0)*x.at(4)*w*w*w*w - x.at(0)*x.at(3)*w*w*w*w - (x.at(0)*x.at(4)*x.at(3)*x.at(7)*x.at(2)*w*w*w*w*w*w*w*w - (x.at(0)*x.at(4)*w*w*w*w*w*w + x.at(0)*x.at(3)*w*w*w*w*w*w)*x.at(7))*x.at(6))*x.at(1))*x.at(5) - 1)/(x.at(3)*x.at(7)*x.at(2)*x.at(6)*w*w*w*w);
        return ABCD;
    }


    if (!topology.compare("30303030")) //LPS4
    {
        ABCD(0,0) = -(x.at(2)*w*w + x.at(4)*w*w + x.at(6)*w*w)*x.at(1) - (x.at(4)*w*w + x.at(6)*w*w - (x.at(2)*x.at(4)*w*w*w*w + x.at(2)*x.at(6)*w*w*w*w)*x.at(1))*x.at(3) + ((x.at(2)*w*w*w*w + x.at(4)*w*w*w*w)*x.at(6)*x.at(1) - x.at(6)*w*w - (x.at(2)*x.at(4)*x.at(6)*x.at(1)*w*w*w*w*w*w - x.at(4)*x.at(6)*w*w*w*w)*x.at(3))*x.at(5) + 1;
        ABCD(0,1) = (-I*x.at(2)*x.at(1)*w*w*w + I*w)*x.at(3) + ((-I*x.at(2)*w*w*w - I*x.at(4)*w*w*w)*x.at(1) + (I*x.at(2)*x.at(4)*x.at(1)*w*w*w*w*w - I*x.at(4)*w*w*w)*x.at(3) + I*w)*x.at(5) + ((-I*x.at(2)*w*w*w - I*x.at(4)*w*w*w - I*x.at(6)*w*w*w)*x.at(1) + (-I*x.at(4)*w*w*w - I*x.at(6)*w*w*w + (I*x.at(2)*x.at(4)*w*w*w*w*w + I*x.at(2)*x.at(6)*w*w*w*w*w)*x.at(1))*x.at(3) + (-I*x.at(6)*w*w*w + (I*x.at(2)*w*w*w*w*w + I*x.at(4)*w*w*w*w*w)*x.at(6)*x.at(1) + (-I*x.at(2)*x.at(4)*x.at(6)*x.at(1)*w*w*w*w*w*w*w + I*x.at(4)*x.at(6)*w*w*w*w*w)*x.at(3))*x.at(5) + I*w)*x.at(7) + I*x.at(1)*w;
        ABCD(1,0) = (-I*x.at(0)*x.at(2)*w*w*w - I*x.at(0)*x.at(4)*w*w*w - I*x.at(0)*x.at(6)*w*w*w)*x.at(1) + ((-I*x.at(0)*w*w*w - I*x.at(2)*w*w*w)*x.at(4) + (-I*x.at(0)*w*w*w - I*x.at(2)*w*w*w)*x.at(6) + (I*x.at(0)*x.at(2)*x.at(4)*w*w*w*w*w + I*x.at(0)*x.at(2)*x.at(6)*w*w*w*w*w)*x.at(1))*x.at(3) + ((I*x.at(0)*x.at(2)*w*w*w*w*w + I*x.at(0)*x.at(4)*w*w*w*w*w)*x.at(6)*x.at(1) + (-I*x.at(0)*w*w*w - I*x.at(2)*w*w*w - I*x.at(4)*w*w*w)*x.at(6) + (-I*x.at(0)*x.at(2)*x.at(4)*x.at(6)*x.at(1)*w*w*w*w*w*w*w + (I*x.at(0)*w*w*w*w*w + I*x.at(2)*w*w*w*w*w)*x.at(4)*x.at(6))*x.at(3))*x.at(5) + I*x.at(0)*w + I*x.at(2)*w + I*x.at(4)*w + I*x.at(6)*w;
        ABCD(1,1) = -x.at(0)*x.at(1)*w*w + (x.at(0)*x.at(2)*x.at(1)*w*w*w*w - x.at(0)*w*w - x.at(2)*w*w)*x.at(3) - (x.at(0)*w*w + x.at(2)*w*w + x.at(4)*w*w - (x.at(0)*x.at(2)*w*w*w*w + x.at(0)*x.at(4)*w*w*w*w)*x.at(1) + (x.at(0)*x.at(2)*x.at(4)*x.at(1)*w*w*w*w*w*w - (x.at(0)*w*w*w*w + x.at(2)*w*w*w*w)*x.at(4))*x.at(3))*x.at(5) - (x.at(0)*w*w + x.at(2)*w*w + x.at(4)*w*w + x.at(6)*w*w - (x.at(0)*x.at(2)*w*w*w*w + x.at(0)*x.at(4)*w*w*w*w + x.at(0)*x.at(6)*w*w*w*w)*x.at(1) - ((x.at(0)*w*w*w*w + x.at(2)*w*w*w*w)*x.at(4) + (x.at(0)*w*w*w*w + x.at(2)*w*w*w*w)*x.at(6) - (x.at(0)*x.at(2)*x.at(4)*w*w*w*w*w*w + x.at(0)*x.at(2)*x.at(6)*w*w*w*w*w*w)*x.at(1))*x.at(3) + ((x.at(0)*x.at(2)*w*w*w*w*w*w + x.at(0)*x.at(4)*w*w*w*w*w*w)*x.at(6)*x.at(1) - (x.at(0)*w*w*w*w + x.at(2)*w*w*w*w + x.at(4)*w*w*w*w)*x.at(6) - (x.at(0)*x.at(2)*x.at(4)*x.at(6)*x.at(1)*w*w*w*w*w*w*w*w - (x.at(0)*w*w*w*w*w*w + x.at(2)*w*w*w*w*w*w)*x.at(4)*x.at(6))*x.at(3))*x.at(5))*x.at(7) + 1;
        return ABCD;
    }

    if (!topology.compare("03030303")) //LPP4
    {
        ABCD(0,0) = -(x.at(1)*w*w + x.at(3)*w*w + x.at(5)*w*w + x.at(7)*w*w)*x.at(0) - (x.at(3)*w*w + x.at(5)*w*w + x.at(7)*w*w - (x.at(1)*x.at(3)*w*w*w*w + x.at(1)*x.at(5)*w*w*w*w + x.at(1)*x.at(7)*w*w*w*w)*x.at(0))*x.at(2) - (x.at(5)*w*w + x.at(7)*w*w - ((x.at(1)*w*w*w*w + x.at(3)*w*w*w*w)*x.at(5) + (x.at(1)*w*w*w*w + x.at(3)*w*w*w*w)*x.at(7))*x.at(0) - (x.at(3)*x.at(5)*w*w*w*w + x.at(3)*x.at(7)*w*w*w*w - (x.at(1)*x.at(3)*x.at(5)*w*w*w*w*w*w + x.at(1)*x.at(3)*x.at(7)*w*w*w*w*w*w)*x.at(0))*x.at(2))*x.at(4) + ((x.at(1)*w*w*w*w + x.at(3)*w*w*w*w + x.at(5)*w*w*w*w)*x.at(7)*x.at(0) - x.at(7)*w*w - ((x.at(1)*x.at(3)*w*w*w*w*w*w + x.at(1)*x.at(5)*w*w*w*w*w*w)*x.at(7)*x.at(0) - (x.at(3)*w*w*w*w + x.at(5)*w*w*w*w)*x.at(7))*x.at(2) + (x.at(5)*x.at(7)*w*w*w*w - (x.at(1)*w*w*w*w*w*w + x.at(3)*w*w*w*w*w*w)*x.at(5)*x.at(7)*x.at(0) + (x.at(1)*x.at(3)*x.at(5)*x.at(7)*x.at(0)*w*w*w*w*w*w*w*w - x.at(3)*x.at(5)*x.at(7)*w*w*w*w*w*w)*x.at(2))*x.at(4))*x.at(6) + 1;
        ABCD(0,1) = (-I*x.at(1)*x.at(0)*w*w*w + I*w)*x.at(2) + ((-I*x.at(1)*w*w*w - I*x.at(3)*w*w*w)*x.at(0) + (I*x.at(1)*x.at(3)*x.at(0)*w*w*w*w*w - I*x.at(3)*w*w*w)*x.at(2) + I*w)*x.at(4) + ((-I*x.at(1)*w*w*w - I*x.at(3)*w*w*w - I*x.at(5)*w*w*w)*x.at(0) + (-I*x.at(3)*w*w*w - I*x.at(5)*w*w*w + (I*x.at(1)*x.at(3)*w*w*w*w*w + I*x.at(1)*x.at(5)*w*w*w*w*w)*x.at(0))*x.at(2) + (-I*x.at(5)*w*w*w + (I*x.at(1)*w*w*w*w*w + I*x.at(3)*w*w*w*w*w)*x.at(5)*x.at(0) + (-I*x.at(1)*x.at(3)*x.at(5)*x.at(0)*w*w*w*w*w*w*w + I*x.at(3)*x.at(5)*w*w*w*w*w)*x.at(2))*x.at(4) + I*w)*x.at(6) + I*x.at(0)*w;
        ABCD(1,0) = (-I*x.at(1)*x.at(3)*w*w*w - I*x.at(1)*x.at(5)*w*w*w - I*x.at(1)*x.at(7)*w*w*w)*x.at(2) + ((-I*x.at(1)*w*w*w - I*x.at(3)*w*w*w)*x.at(5) + (-I*x.at(1)*w*w*w - I*x.at(3)*w*w*w)*x.at(7) + (I*x.at(1)*x.at(3)*x.at(5)*w*w*w*w*w + I*x.at(1)*x.at(3)*x.at(7)*w*w*w*w*w)*x.at(2))*x.at(4) + ((I*x.at(1)*x.at(3)*w*w*w*w*w + I*x.at(1)*x.at(5)*w*w*w*w*w)*x.at(7)*x.at(2) + (-I*x.at(1)*w*w*w - I*x.at(3)*w*w*w - I*x.at(5)*w*w*w)*x.at(7) + (-I*x.at(1)*x.at(3)*x.at(5)*x.at(7)*x.at(2)*w*w*w*w*w*w*w + (I*x.at(1)*w*w*w*w*w + I*x.at(3)*w*w*w*w*w)*x.at(5)*x.at(7))*x.at(4))*x.at(6) + I*x.at(1)*w + I*x.at(3)*w + I*x.at(5)*w + I*x.at(7)*w;
        ABCD(1,1) = -x.at(1)*x.at(2)*w*w + (x.at(1)*x.at(3)*x.at(2)*w*w*w*w - x.at(1)*w*w - x.at(3)*w*w)*x.at(4) - (x.at(1)*w*w + x.at(3)*w*w + x.at(5)*w*w - (x.at(1)*x.at(3)*w*w*w*w + x.at(1)*x.at(5)*w*w*w*w)*x.at(2) + (x.at(1)*x.at(3)*x.at(5)*x.at(2)*w*w*w*w*w*w - (x.at(1)*w*w*w*w + x.at(3)*w*w*w*w)*x.at(5))*x.at(4))*x.at(6) + 1;
        return ABCD;
    }

    if (!topology.compare("12121212")) //HPP
    {
        ABCD(0,0) = -((x.at(0)*w*w + x.at(2)*w*w)*x.at(1) + (x.at(2)*w*w + x.at(4)*w*w - (x.at(0)*x.at(2)*w*w*w*w + (x.at(0)*w*w*w*w + x.at(2)*w*w*w*w)*x.at(4))*x.at(1))*x.at(3) + (x.at(4)*w*w + x.at(6)*w*w - ((x.at(0)*w*w*w*w + x.at(2)*w*w*w*w)*x.at(4) + (x.at(0)*w*w*w*w + x.at(2)*w*w*w*w)*x.at(6))*x.at(1) - (x.at(2)*x.at(4)*w*w*w*w + (x.at(2)*w*w*w*w + x.at(4)*w*w*w*w)*x.at(6) - (x.at(0)*x.at(2)*x.at(4)*w*w*w*w*w*w + (x.at(0)*x.at(2)*w*w*w*w*w*w + (x.at(0)*w*w*w*w*w*w + x.at(2)*w*w*w*w*w*w)*x.at(4))*x.at(6))*x.at(1))*x.at(3))*x.at(5) - ((x.at(0)*w*w*w*w + x.at(2)*w*w*w*w)*x.at(6)*x.at(1) - x.at(6)*w*w - ((x.at(0)*x.at(2)*w*w*w*w*w*w + (x.at(0)*w*w*w*w*w*w + x.at(2)*w*w*w*w*w*w)*x.at(4))*x.at(6)*x.at(1) - (x.at(2)*w*w*w*w + x.at(4)*w*w*w*w)*x.at(6))*x.at(3) + (x.at(4)*x.at(6)*w*w*w*w - (x.at(0)*w*w*w*w*w*w + x.at(2)*w*w*w*w*w*w)*x.at(4)*x.at(6)*x.at(1) + (x.at(0)*x.at(2)*x.at(4)*x.at(6)*x.at(1)*w*w*w*w*w*w*w*w - x.at(2)*x.at(4)*x.at(6)*w*w*w*w*w*w)*x.at(3))*x.at(5))*x.at(7) - 1)/(x.at(0)*x.at(2)*x.at(4)*x.at(6)*x.at(1)*x.at(3)*x.at(5)*x.at(7)*w*w*w*w*w*w*w*w);
        ABCD(0,1) = -((I*x.at(0)*w*w + I*x.at(2)*w*w)*x.at(1) + (I*x.at(2)*w*w + I*x.at(4)*w*w + (-I*x.at(0)*x.at(2)*w*w*w*w + (-I*x.at(0)*w*w*w*w - I*x.at(2)*w*w*w*w)*x.at(4))*x.at(1))*x.at(3) + (I*x.at(4)*w*w + I*x.at(6)*w*w + ((-I*x.at(0)*w*w*w*w - I*x.at(2)*w*w*w*w)*x.at(4) + (-I*x.at(0)*w*w*w*w - I*x.at(2)*w*w*w*w)*x.at(6))*x.at(1) + (-I*x.at(2)*x.at(4)*w*w*w*w + (-I*x.at(2)*w*w*w*w - I*x.at(4)*w*w*w*w)*x.at(6) + (I*x.at(0)*x.at(2)*x.at(4)*w*w*w*w*w*w + (I*x.at(0)*x.at(2)*w*w*w*w*w*w + (I*x.at(0)*w*w*w*w*w*w + I*x.at(2)*w*w*w*w*w*w)*x.at(4))*x.at(6))*x.at(1))*x.at(3))*x.at(5) - I)/(x.at(0)*x.at(2)*x.at(4)*x.at(6)*x.at(1)*x.at(3)*x.at(5)*w*w*w*w*w*w*w);
        ABCD(1,0) = (-I*x.at(2)*x.at(1)*w*w - (-I*x.at(2)*x.at(4)*x.at(1)*w*w*w*w + I*x.at(2)*w*w + I*x.at(4)*w*w)*x.at(3) - (I*x.at(4)*w*w + I*x.at(6)*w*w + (-I*x.at(2)*x.at(4)*w*w*w*w - I*x.at(2)*x.at(6)*w*w*w*w)*x.at(1) + (I*x.at(2)*x.at(4)*x.at(6)*x.at(1)*w*w*w*w*w*w - I*x.at(2)*x.at(4)*w*w*w*w + (-I*x.at(2)*w*w*w*w - I*x.at(4)*w*w*w*w)*x.at(6))*x.at(3))*x.at(5) + (I*x.at(2)*x.at(6)*x.at(1)*w*w*w*w - I*x.at(6)*w*w - (I*x.at(2)*x.at(4)*x.at(6)*x.at(1)*w*w*w*w*w*w + (-I*x.at(2)*w*w*w*w - I*x.at(4)*w*w*w*w)*x.at(6))*x.at(3) - (I*x.at(2)*x.at(4)*x.at(6)*x.at(1)*w*w*w*w*w*w + I*x.at(2)*x.at(4)*x.at(6)*x.at(3)*w*w*w*w*w*w - I*x.at(4)*x.at(6)*w*w*w*w)*x.at(5))*x.at(7) + I)/(x.at(2)*x.at(4)*x.at(6)*x.at(1)*x.at(3)*x.at(5)*x.at(7)*w*w*w*w*w*w*w);
        ABCD(1,1) = (x.at(2)*x.at(1)*w*w - (x.at(2)*x.at(4)*x.at(1)*w*w*w*w - x.at(2)*w*w - x.at(4)*w*w)*x.at(3) + (x.at(4)*w*w + x.at(6)*w*w - (x.at(2)*x.at(4)*w*w*w*w + x.at(2)*x.at(6)*w*w*w*w)*x.at(1) + (x.at(2)*x.at(4)*x.at(6)*x.at(1)*w*w*w*w*w*w - x.at(2)*x.at(4)*w*w*w*w - (x.at(2)*w*w*w*w + x.at(4)*w*w*w*w)*x.at(6))*x.at(3))*x.at(5) - 1)/(x.at(2)*x.at(4)*x.at(6)*x.at(1)*x.at(3)*x.at(5)*w*w*w*w*w*w);
        return ABCD;
    }

    if (!topology.compare("21212121")) //HPS
    {
        ABCD(0,0) = ((x.at(1)*w*w + x.at(3)*w*w)*x.at(2) + (x.at(3)*w*w + x.at(5)*w*w - (x.at(1)*x.at(3)*w*w*w*w + (x.at(1)*w*w*w*w + x.at(3)*w*w*w*w)*x.at(5))*x.at(2))*x.at(4) - ((x.at(1)*w*w*w*w + x.at(3)*w*w*w*w)*x.at(5)*x.at(2) - x.at(5)*w*w - (x.at(1)*x.at(3)*x.at(5)*x.at(2)*w*w*w*w*w*w - x.at(3)*x.at(5)*w*w*w*w)*x.at(4))*x.at(6) - 1)/(x.at(1)*x.at(3)*x.at(5)*x.at(2)*x.at(4)*x.at(6)*w*w*w*w*w*w);
        ABCD(0,1) = -((I*x.at(1)*w*w + I*x.at(3)*w*w)*x.at(2) + (I*x.at(3)*w*w + I*x.at(5)*w*w + (-I*x.at(1)*x.at(3)*w*w*w*w + (-I*x.at(1)*w*w*w*w - I*x.at(3)*w*w*w*w)*x.at(5))*x.at(2))*x.at(4) + (I*x.at(5)*w*w + I*x.at(7)*w*w + ((-I*x.at(1)*w*w*w*w - I*x.at(3)*w*w*w*w)*x.at(5) + (-I*x.at(1)*w*w*w*w - I*x.at(3)*w*w*w*w)*x.at(7))*x.at(2) + (-I*x.at(3)*x.at(5)*w*w*w*w + (-I*x.at(3)*w*w*w*w - I*x.at(5)*w*w*w*w)*x.at(7) + (I*x.at(1)*x.at(3)*x.at(5)*w*w*w*w*w*w + (I*x.at(1)*x.at(3)*w*w*w*w*w*w + (I*x.at(1)*w*w*w*w*w*w + I*x.at(3)*w*w*w*w*w*w)*x.at(5))*x.at(7))*x.at(2))*x.at(4))*x.at(6) - I)/(x.at(1)*x.at(3)*x.at(5)*x.at(7)*x.at(2)*x.at(4)*x.at(6)*w*w*w*w*w*w*w);
        ABCD(1,0) = (-I*x.at(1)*x.at(0)*w*w - (-I*x.at(1)*x.at(3)*x.at(0)*w*w*w*w + I*x.at(1)*w*w + I*x.at(3)*w*w)*x.at(2) - (I*x.at(3)*w*w + I*x.at(5)*w*w + (-I*x.at(1)*x.at(3)*w*w*w*w - I*x.at(1)*x.at(5)*w*w*w*w)*x.at(0) + (I*x.at(1)*x.at(3)*x.at(5)*x.at(0)*w*w*w*w*w*w - I*x.at(1)*x.at(3)*w*w*w*w + (-I*x.at(1)*w*w*w*w - I*x.at(3)*w*w*w*w)*x.at(5))*x.at(2))*x.at(4) + (I*x.at(1)*x.at(5)*x.at(0)*w*w*w*w - I*x.at(5)*w*w - (I*x.at(1)*x.at(3)*x.at(5)*x.at(0)*w*w*w*w*w*w + (-I*x.at(1)*w*w*w*w - I*x.at(3)*w*w*w*w)*x.at(5))*x.at(2) - (I*x.at(1)*x.at(3)*x.at(5)*x.at(0)*w*w*w*w*w*w + I*x.at(1)*x.at(3)*x.at(5)*x.at(2)*w*w*w*w*w*w - I*x.at(3)*x.at(5)*w*w*w*w)*x.at(4))*x.at(6) + I)/(x.at(1)*x.at(3)*x.at(5)*x.at(0)*x.at(2)*x.at(4)*x.at(6)*w*w*w*w*w*w*w);
        ABCD(1,1) = -(x.at(1)*x.at(0)*w*w - (x.at(1)*x.at(3)*x.at(0)*w*w*w*w - x.at(1)*w*w - x.at(3)*w*w)*x.at(2) + (x.at(3)*w*w + x.at(5)*w*w - (x.at(1)*x.at(3)*w*w*w*w + x.at(1)*x.at(5)*w*w*w*w)*x.at(0) + (x.at(1)*x.at(3)*x.at(5)*x.at(0)*w*w*w*w*w*w - x.at(1)*x.at(3)*w*w*w*w - (x.at(1)*w*w*w*w + x.at(3)*w*w*w*w)*x.at(5))*x.at(2))*x.at(4) + (x.at(5)*w*w + x.at(7)*w*w - (x.at(1)*x.at(5)*w*w*w*w + x.at(1)*x.at(7)*w*w*w*w)*x.at(0) - ((x.at(1)*w*w*w*w + x.at(3)*w*w*w*w)*x.at(5) + (x.at(1)*w*w*w*w + x.at(3)*w*w*w*w)*x.at(7) - (x.at(1)*x.at(3)*x.at(5)*w*w*w*w*w*w + x.at(1)*x.at(3)*x.at(7)*w*w*w*w*w*w)*x.at(0))*x.at(2) - (x.at(3)*x.at(5)*w*w*w*w + (x.at(3)*w*w*w*w + x.at(5)*w*w*w*w)*x.at(7) - (x.at(1)*x.at(3)*x.at(5)*w*w*w*w*w*w + (x.at(1)*x.at(3)*w*w*w*w*w*w + x.at(1)*x.at(5)*w*w*w*w*w*w)*x.at(7))*x.at(0) + (x.at(1)*x.at(3)*x.at(5)*x.at(7)*x.at(0)*w*w*w*w*w*w*w*w - x.at(1)*x.at(3)*x.at(5)*w*w*w*w*w*w - (x.at(1)*x.at(3)*w*w*w*w*w*w + (x.at(1)*w*w*w*w*w*w + x.at(3)*w*w*w*w*w*w)*x.at(5))*x.at(7))*x.at(2))*x.at(4))*x.at(6) - 1)/(x.at(1)*x.at(3)*x.at(5)*x.at(7)*x.at(0)*x.at(2)*x.at(4)*x.at(6)*w*w*w*w*w*w*w*w);
        return ABCD;
    }

    return ABCD;
}


