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
    case 1: return 1.0;
    case 2: return 2.0;
    case 3: return 3.0;
    case 4: return 8.0;
    case 5: return 15.0;
    case 6: return 48.0;
    case 7: return 105.0;
    case 8: return 384.0;
    case 9: return 945.0;
    case 10: return 3840.0;
    case 11: return 10395.0;
    case 12: return 46080.0;
    case 13: return 135135.0;
    case 14: return 645120.0;
    default:
        if(n < 0)
        {
            cout << "ERROR: double_factorial is called for a negative number!" << endl;
            exit(99);
        }
        return n * double_factorial(n - 2);
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
        return n * factorial(n - 1);
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
    int tmp_int;
    size = 5;
    gtos_c.resize(5);
    gtos_c(0).gto_list.resize(4);
    gtos_c(1).gto_list.resize(1);
    gtos_c(2).gto_list.resize(1);
    gtos_c(3).gto_list.resize(1);
    gtos_c(4).gto_list.resize(1);
    gtos_c(0).coeff.resize(4);
    gtos_c(1).coeff.resize(1);
    gtos_c(2).coeff.resize(1);
    gtos_c(3).coeff.resize(1);
    gtos_c(4).coeff.resize(1);
    ifs.open("ccpvdz");
        ifs >> atmName >> tmp;
        ifs >> orbAng >> tmp_int >> tmp;
        for(int ii = 0; ii < 4; ii++)
        {
            ifs >> gtos_c(0).gto_list(ii).a >> gtos_c(0).coeff(ii);
            gtos_c(0).gto_list(ii).l = 0;
            gtos_c(0).gto_list(ii).m = 0;
        }
        ifs >> orbAng >> tmp_int >> tmp;
        for(int ii = 0; ii < 1; ii++)
        {
            ifs >> gtos_c(1).gto_list(ii).a >> gtos_c(1).coeff(ii);
            gtos_c(1).gto_list(ii).l = 0;
            gtos_c(1).gto_list(ii).m = 0;
        }
        ifs >> orbAng >> tmp_int >> tmp;
        for(int ii = 0; ii < 1; ii++)
        {
            ifs >> gtos_c(2).gto_list(ii).a >> gtos_c(2).coeff(ii);
            gtos_c(3).gto_list(ii).a = gtos_c(4).gto_list(ii).a = gtos_c(2).gto_list(ii).a;
            gtos_c(3).coeff(ii) = gtos_c(4).coeff(ii) = gtos_c(2).coeff(ii);
            gtos_c(2).gto_list(ii).l = 1;    gtos_c(2).gto_list(ii).m = 0;
            gtos_c(3).gto_list(ii).l = 1;    gtos_c(3).gto_list(ii).m = -1;
            gtos_c(4).gto_list(ii).l = 1;    gtos_c(4).gto_list(ii).m = 1;
        }
    ifs.close();
}

void GTO::normalization()
{
    for(int ss = 0; ss < size; ss++)
    {
        double tmp = 0.0;
        int size_gtos = gtos_c(ss).coeff.rows();
        for(int ii = 0; ii < size_gtos; ii++)
        for(int jj = 0; jj < size_gtos; jj++)
        {
            tmp += gtos_c(ss).coeff(ii) * gtos_c(ss).coeff(jj) * overlap_single_gto(gtos_c(ss).gto_list(ii), gtos_c(ss).gto_list(jj));
        }
        // cout << tmp << endl;
        gtos_c(ss).coeff = gtos_c(ss).coeff / sqrt(tmp);
    }
    
    return;
}

MatrixXd GTO::get_overlap()
{
    MatrixXd s(size, size);
    s = s * 0.0;
    for(int ii = 0; ii < size; ii++)
    for(int jj = 0; jj < size; jj++)
    {
        int size_ii = gtos_c(ii).coeff.rows(), size_jj = gtos_c(jj).coeff.rows();
        for(int mm = 0; mm < size_ii; mm++)
        for(int nn = 0; nn < size_jj; nn++)
        {
            s(ii,jj) += gtos_c(ii).coeff(mm) * gtos_c(jj).coeff(nn) * overlap_single_gto(gtos_c(ii).gto_list(mm), gtos_c(jj).gto_list(nn));
        }
    }

    return s;
}


// double GTO::get_h1e()
// {
//     double tmp = 0.0;
//     for(int ii = 0; ii < size; ii++)
//     for(int jj = 0; jj < size; jj++)
//     {
//         tmp += coeff(ii) * coeff(jj) * h_1e(gtos(ii), gtos(jj));
//     }
    
//     return tmp;
// }





double GTO::auxiliary(const int& l, const double& a)
{
    int n = l / 2;
    if(n*2 == l)    return double_factorial(2*n-1)/pow(a,n)/pow(2.0,n+1)*sqrt(PI/a);
    else    return factorial(n) /2.0/pow(a,n+1);
}

double GTO::overlap_single_gto(const gto_single& gto1, const gto_single& gto2)
{
    if(gto1.l != gto2.l || gto1.m != gto2.m) return 0.0;
    else return auxiliary(2 + 2*gto1.l, gto1.a + gto2.a);
}

double GTO::nuc_attra_single_gto(const gto_single& gto1, const gto_single& gto2)
{
    if(gto1.l != gto2.l || gto1.m != gto2.m) return 0.0;
    else return -atomN * auxiliary(1 + 2*gto1.l, gto1.a + gto2.a);
}

double GTO::kinetic_single_gto(const gto_single& gto1, const gto_single& gto2)
{
    if(gto1.l != gto2.l || gto1.m != gto2.m) return 0.0;
    else return gto2.a * (2*gto1.l + 3) * auxiliary(2 + 2*gto1.l, gto1.a + gto2.a) - 2 * gto2.a * gto2.a * auxiliary(4 + 2*gto1.l, gto1.a + gto2.a);
}

double GTO::h_1e_single_gto(const gto_single& gto1, const gto_single& gto2)
{
    return kinetic_single_gto(gto1, gto2) + nuc_attra_single_gto(gto1, gto2);
}