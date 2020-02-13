#include<Eigen/Dense>
#include<string>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
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
        else return n * double_factorial(n - 2);
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
    gtos_c(0).gto_list.resize(3);
    gtos_c(1).gto_list.resize(1);
    gtos_c(2).gto_list.resize(1);
    gtos_c(3).gto_list.resize(1);
    gtos_c(4).gto_list.resize(1);
    gtos_c(0).coeff.resize(3);
    gtos_c(1).coeff.resize(1);
    gtos_c(2).coeff.resize(1);
    gtos_c(3).coeff.resize(1);
    gtos_c(4).coeff.resize(1);
    ifs.open("ccpvdz");
        ifs >> atmName >> tmp;
        ifs >> orbAng >> tmp_int >> tmp;
        for(int ii = 0; ii < 3; ii++)
        {
            ifs >> gtos_c(0).gto_list(ii).a >> gtos_c(0).coeff(ii);
            gtos_c(0).gto_list(ii).l = 0;
            gtos_c(0).gto_list(ii).m = 0;

            gtos_c(0).coeff(ii) = gtos_c(0).coeff(ii) / sqrt(auxiliary_1e(2*gtos_c(0).gto_list(ii).l + 2, 2 * gtos_c(0).gto_list(ii).a));
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
            tmp += gtos_c(ss).coeff(ii) * gtos_c(ss).coeff(jj) * auxiliary_1e(2+gtos_c(ss).gto_list(ii).l+gtos_c(ss).gto_list(jj).l, gtos_c(ss).gto_list(ii).a+gtos_c(ss).gto_list(jj).a);
        }
        gtos_c(ss).coeff = gtos_c(ss).coeff / sqrt(tmp);
    }
    
    return;
}

MatrixXd GTO::get_h1e(const string& intType)
{
    MatrixXd int_1e(size, size);
    int_1e = int_1e * 0.0;
    for(int ii = 0; ii < size; ii++)
    for(int jj = 0; jj < size; jj++)
    {
        int size_ii = gtos_c(ii).coeff.rows(), size_jj = gtos_c(jj).coeff.rows();
        for(int mm = 0; mm < size_ii; mm++)
        for(int nn = 0; nn < size_jj; nn++)
        {
            int_1e(ii,jj) += gtos_c(ii).coeff(mm) * gtos_c(jj).coeff(nn) * int1e_single_gto(gtos_c(ii).gto_list(mm), gtos_c(jj).gto_list(nn), intType);
        }
    }

    return int_1e;
}

Matrix<MatrixXd, -1, -1> GTO::get_h2e()
{
    Matrix<MatrixXd, -1, -1> int_2e(size, size);
    for(int ii = 0; ii < size; ii++)
    for(int jj = 0; jj < size; jj++)
    {
        int_2e(ii,jj).resize(size, size);
        for(int kk = 0; kk < size; kk++)
        for(int ll = 0; ll < size; ll++)
        {
            int size_i = gtos_c(ii).coeff.rows(), size_j = gtos_c(jj).coeff.rows(), size_k = gtos_c(kk).coeff.rows(), size_l = gtos_c(ll).coeff.rows();
            int_2e(ii,jj)(kk,ll) = 0.0;
            for(int ni = 0; ni < size_i; ni++)
            for(int nj = 0; nj < size_j; nj++)
            for(int nk = 0; nk < size_k; nk++)
            for(int nl = 0; nl < size_l; nl++)
            int_2e(ii,jj)(kk,ll) += gtos_c(ii).coeff(ni) * gtos_c(jj).coeff(nj) * gtos_c(kk).coeff(nk) * gtos_c(ll).coeff(nl) * int2e_single_gto(gtos_c(ii).gto_list(ni), gtos_c(jj).gto_list(nj),gtos_c(kk).gto_list(nk), gtos_c(ll).gto_list(nl));
        }
    }
    return int_2e;
}



/*
    auxiliary_1e is to evaluate \int_0^inf x^l exp(-ax^2) dx
*/
double GTO::auxiliary_1e(const int& l, const double& a)
{
    int n = l / 2;
    if(n*2 == l)    return double_factorial(2*n-1)/pow(a,n)/pow(2.0,n+1)*sqrt(M_PI/a);
    else    return factorial(n) /2.0/pow(a,n+1);
}


double GTO::int1e_single_gto(const gto_single& gto1, const gto_single& gto2, const string& integralTYPE)
{
    int l = gto1.l;
    double a1 = gto1.a, a2 = gto2.a;
    
    if(gto1.l != gto2.l || gto1.m != gto2.m) return 0.0;
    else if(integralTYPE == "overlap")  return auxiliary_1e(2 + 2*l, a1 + a2);
    else if(integralTYPE == "nuc_attra")  return -atomN * auxiliary_1e(1 + 2*l, a1 + a2);
    else if(integralTYPE == "kinetic")  return a2 * (2*l + 3) * auxiliary_1e(2 + 2*l, a1 + a2) - 2 * a2 * a2 * auxiliary_1e(4 + 2*l, a1 + a2);
    else if(integralTYPE == "h1e")  return int1e_single_gto(gto1, gto2, "kinetic") + int1e_single_gto(gto1, gto2, "nuc_attra");
    else if(integralTYPE == "p.Vp")  return ((2*l + 1) * l * auxiliary_1e(2*l-1, a1 + a2) - 2*l*(a1+a2)*auxiliary_1e(2*l + 1, a1 + a2) + 4*a1*a2 * auxiliary_1e(2*l + 3, a1 + a2)) * -atomN;
    else
    {
        cout << "ERROR: get_h1e is called for undefined type of integrals!" << endl;
        exit(99);
    }
}


/*
    auxiliary_2e is to evaluate \int_0^inf\int_0^inf x1^l1 x2^l2 |x1-x2|^-1 exp(-a1 * x1^2) exp(-a2 * x2^2) dx1dx2
*/
double GTO::auxiliary_2e(const int& l1, const int& l2, const double& a1, const double& a2)
{
    double res = 0.0, l1f = factorial(l1), l2f = factorial(l2);
    for(int ii = 0; ii <= l1; ii++) 
    {
        res += l1f * pow(a1, ii - l1 - 1) * pow(a1 + a2, -l2-ii-1.5) * double_factorial(2*l2 + 2*ii + 1)/ factorial(ii) / pow(2.0, l2 + ii + 1);
    }
    for(int ii = 0; ii <= l1; ii++) 
    {
        res += l2f * pow(a2, ii - l2 - 1) * pow(a1 + a2, -l1-ii-1.5) * double_factorial(2*l1 + 2*ii + 1)/ factorial(ii) / pow(2.0, l1 + ii + 1);
    }
    return res * 4.0 * pow(M_PI, 1.5);
}

double GTO::int2e_single_gto(const gto_single& gto1, const gto_single& gto2, const gto_single& gto3, const gto_single& gto4)
{
    if(gto1.l != gto2.l || gto1.m != gto2.m || gto3.l != gto4.l || gto3.m != gto4.m) return 0.0;
    else    return auxiliary_2e(2 + gto1.l, 2 + gto3.l, gto1.a + gto2.a, gto3.a + gto4.a);
}