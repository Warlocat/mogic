#include<Eigen/Dense>
#include<string>
#include<iostream>
#include<iomanip>
#include<fstream>
#include"gto.h"
using namespace std;
using namespace Eigen;

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
            
        }

        gtos.resize(size);
        coeff.resize(size);
        for(int ii = 0; ii < size; ii++)
        {
            ifs >> gtos(ii).a >> coeff(ii);
            // cout << gtos(ii).a << endl;
            // cout << coeff(ii) << endl;
            gtos(ii).lx = 0;
            gtos(ii).ly = 0;
            gtos(ii).lz = 0;
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

double GTO::overlap(const gto_single& gto1, const gto_single& gto2)
{
    return sqrt(3.1415926535897932384626 / (gto1.a + gto2.a));
}