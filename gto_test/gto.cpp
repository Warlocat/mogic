#include<Eigen/Dense>
#include<string>
#include<iostream>
#include<iomanip>
#include<fstream>
#include"gto.h"
using namespace std;
using namespace Eigen;

double double_factorial(const int& n)
{
    switch (n)
    {
    case 0: return 1.0;
    case 2: return 2.0;
    case 4: return 8.0;
    case 6: return 48.0;
    case 8: return 384.0;
    case 10: return 3840.0;
    case 12: return 46080.0;
    case 14: return 645120.0;
    default:
        if(n < 0)
        {
            cout << "ERROR: double_factorial is called for a negative number!" << endl;
            exit(99);
        }
        double tmp = 645120.0;
        for(int ii = 16; ii <= n; ii = ii + 2)
        {
            tmp = tmp * ii;
        }
        return tmp;
    }
}
double factorial(const int& n)
{
    switch (n)
    {
    case 0: return 1.0;
    case 1: return 1.0;
    case 2: return 2.0;
    case 3: return 6.0;
    case 4: return 24.0;
    case 5: return 120.0;
    case 6: return 720.0;
    case 7: return 5040.0;
    case 8: return 40320.0;
    case 9: return 362880.0;
    case 10: return 3628800.0;
    
    default:
        if(n < 0)
        {
            cout << "ERROR: factorial is called for a negative number!" << endl;
            exit(99);
        }
        double tmp = 3628800.0;
        for(int ii = 11; ii <= n; ii = ii + 1)
        {
            tmp = tmp * ii;
        }
        return tmp;
    }
}

GTO::GTO()
{
}

GTO::~GTO()
{
}

void GTO::readBasis()
{
    string atmName;
    string orbAng;
    double tmp;
    ifstream ifs;
    ifs.open("ccpvdz");
        ifs >> atmName >> tmp;
        ifs >> orbAng >> size >> tmp;
        if(orbAng == "S")
        {
            angular = 0;
        }
        else
        {
            cout << "GOTO ERROR" << endl;
            exit(999);
        }

        gtos.resize(size);
        coeff.resize(size);
        for(int ii = 0; ii < size; ii++)
        {
            ifs >> gtos(ii).a >> coeff(ii);
            gtos(ii).l = 0;
            gtos(ii).m = 0;
        }

        
    ifs.close();
}

void GTO::normalization()
{
    double tmp = 0.0;
    for(int ii = 0; ii < size; ii++)
    for(int jj = 0; jj < size; jj++)
    {
        tmp += coeff(ii) * coeff(jj) * overlap(gtos(ii), gtos(jj));
    }
    cout << tmp << endl;
    coeff = coeff / sqrt(tmp);
    

    return;
}

double GTO::auxiliary(const int& l, const double& a)
{
    int n = l / 2;
    if(n*2 == l)    return double_factorial(2*n-1)/pow(a,n)/pow(2.0,n+1)*sqrt(PI/a);
    else    return factorial(n) /2.0/pow(a,n+1);
}

double GTO::overlap(const gto_single& gto1, const gto_single& gto2)
{
    if(gto1.l != gto2.l || gto1.m != gto2.m) return 0.0;
    else return auxiliary(2 + 2*gto1.l, gto1.a + gto2.a);
}

double GTO::nuc_attraction(const gto_single& gto1, const gto_single& gto2)
{
    if(gto1.l != gto2.l || gto1.m != gto2.m) return 0.0;
    else return -atomN * auxiliary(1 + 2*gto1.l, gto1.a + gto2.a);
}

double GTO::kinetic(const gto_single& gto1, const gto_single& gto2)
{
    if(gto1.l != gto2.l || gto1.m != gto2.m) return 0.0;
    else return gto2.a * (2*gto1.l + 3) * auxiliary(2 + 2*gto1.l, gto1.a + gto2.a) - 2 * gto2.a * gto2.a * auxiliary(4 + 2*gto1.l, gto1.a + gto2.a);
}

double GTO::h_1e(const gto_single& gto1, const gto_single& gto2)
{
    return kinetic(gto1, gto2) + nuc_attraction(gto1, gto2);
}