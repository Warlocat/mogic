#include<Eigen/Dense>
#include<string>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<complex>
#include<omp.h>
#include<gsl/gsl_sf_coupling.h>
#include"gto.h"
using namespace std;
using namespace Eigen;

/*
    factorial and double_factorial
*/
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

/*
    Wigner 3j coefficients with J = l1 + l2 + l3 is even
*/
double wigner_3j(const int& l1, const int& l2, const int& l3, const int& m1, const int& m2, const int& m3)
{
    return gsl_sf_coupling_3j(2*l1,2*l2,2*l3,2*m1,2*m2,2*m3);


    if(l3 > l1 + l2 || l3 < abs(l1 - l2) || m1 + m2 + m3 != 0)
    {
        return 0.0;
    }
    else
    {
        Vector3i L(l1,l2,l3), M(m1,m2,m3);
        int tmp, Lmax = L.maxCoeff();
        for(int ii = 0; ii <= 1; ii++)
        {
            if(L(ii) == Lmax)
            {
                tmp = L(ii);
                L(ii) = L(2);
                L(2) = tmp;
                tmp = M(ii);
                M(ii) = M(2);
                M(2) = tmp;
                break;
            }
        }

        if(L(2) == L(0) + L(1))
        {
            return pow(-1, L(0) - L(1) - M(2)) * sqrt(factorial(2*L(0)) * factorial(2*L(1)) / factorial(2*L(2) + 1) * factorial(L(2) - M(2)) * factorial(L(2) + M(2)) / factorial(L(0)+M(0)) / factorial(L(0)-M(0)) / factorial(L(1)+M(1)) / factorial(L(1)-M(1)));
        }
        else
        {
            return gsl_sf_coupling_3j(L(0),L(1),L(2),M(0),M(1),M(2));
        }
        
    }
}
/*
    Wigner 3j coefficients with m1 = m2 = m3 = 0
*/
double wigner_3j_zeroM(const int& l1, const int& l2, const int& l3)
{
    int J = l1+l2+l3, g = J/2;
    if(J%2 || l3 > l1 + l2 || l3 < abs(l1 - l2))
    {
        return 0.0;
    }
    else
    {
        return pow(-1,g) * sqrt(factorial(J - 2*l1) * factorial(J - 2*l2) * factorial(J - 2*l3) / factorial(J + 1)) 
                * factorial(g) / factorial(g-l1) / factorial(g-l2) / factorial(g-l3);
    }
}




complex<double> U_SH_trans(const int& mu, const int& mm)
{
    complex<double> result;
    if(abs(mu) != abs(mm)) result = 0.0;
    else if (mu == 0)
    {
        result = 1.0;
    }
    else if (mu > 0)
    {
        if(mu == mm) result = pow(-1.0, mu) / sqrt(2.0);
        else result = 1.0 / sqrt(2.0);
    }
    else
    {
        if(mu == mm) result = 1.0 / sqrt(2.0) * complex<double>(0.0,1.0);
        else result = -pow(-1.0, mu) / sqrt(2.0) * complex<double>(0.0,1.0);
    }
    
    return result;
}


GTO::GTO()
{
}

GTO::~GTO()
{
}

/*
    Read basis file in gaussian format
*/
void GTO::readBasis(const string& atomName, const string& filename)
{
    if(atomName == "H") atomN = 1;
    else if(atomName == "C") atomN = 6;
    else if(atomName == "Cu") atomN = 29;

    ifstream ifs;
    int atom_position = 0, angularQN, nGTO, int_tmp;
    string flag_atom = atomName + "     0", flags, all_info = "", flags_tmp;
    size = 0;
    ifs.open(filename);
        while (!ifs.eof())
        {
            getline(ifs,flags);
            atom_position++;
            if(flags == flag_atom) break;
        }
        if(ifs.eof())
        {
            cout << "ERROR: can not find target atom (" + atomName + ") in the basis set file (" + filename + ")\n";
            exit(99);
        }
        else
        {
            while (true)
            {
                getline(ifs,flags);
                if(flags == "****") break;
                flags.resize(1);
                if(flags == "S") size += 1;
                else if(flags == "P") size += 3;
                else if(flags == "D") size += 5;
                else if(flags == "F") size += 7;
                else if(flags == "G") size += 9;
                else if(flags == "H") size += 11;
                else if(flags == "I") size += 13;
                else if(flags != " ")
                {
                    cout << "ERROR: " << flags << " orbital is not supported now." << endl;
                    exit(99);
                }
            }
            gtos_c.resize(size);
        }       
    ifs.close();
    // ifs.clear();
    ifs.open(filename);
        for(int ii = 0; ii < atom_position; ii++)   getline(ifs,flags);
        int_tmp = 0;
        while (true)
        {
            ifs >> flags;
            if(flags == "****") break;
            else ifs >> nGTO >> flags_tmp;

            if(flags == "S") angularQN = 0;
            else if(flags == "P") angularQN = 1;
            else if(flags == "D") angularQN = 2;
            else if(flags == "F") angularQN = 3;
            else if(flags == "G") angularQN = 4;
            else if(flags == "H") angularQN = 5;
            else if(flags == "I") angularQN = 6;
          
            gtos_c(int_tmp).exp_a.resize(nGTO);
            gtos_c(int_tmp).coeff.resize(nGTO);
            gtos_c(int_tmp).l = angularQN;
            gtos_c(int_tmp).m = -angularQN;

            for(int ii = 0 ; ii < nGTO; ii++)
            {
                ifs >> gtos_c(int_tmp).exp_a(ii) >> gtos_c(int_tmp).coeff(ii);
                gtos_c(int_tmp).coeff(ii) = gtos_c(int_tmp).coeff(ii) / sqrt(auxiliary_1e(2*gtos_c(int_tmp).l + 2, 2 * gtos_c(int_tmp).exp_a(ii)));
            }
            for(int ii = 1; ii < 2*angularQN + 1; ii++)
            {
                gtos_c(int_tmp + ii).exp_a.resize(nGTO);
                gtos_c(int_tmp + ii).coeff.resize(nGTO);
                gtos_c(int_tmp + ii).l = gtos_c(int_tmp).l;
                gtos_c(int_tmp + ii).m = gtos_c(int_tmp).m + ii;
                for(int jj = 0; jj < nGTO; jj++)
                {
                    gtos_c(int_tmp + ii).coeff(jj) = gtos_c(int_tmp).coeff(jj);
                    gtos_c(int_tmp + ii).exp_a(jj) = gtos_c(int_tmp).exp_a(jj);
                }
            }
            int_tmp += 2*angularQN + 1;
        }   
    ifs.close();
}


/*
    Normalization
*/
void GTO::normalization()
{
    for(int ss = 0; ss < size; ss++)
    {
        double tmp = 0.0;
        int size_gtos = gtos_c(ss).coeff.rows();
        for(int ii = 0; ii < size_gtos; ii++)
        for(int jj = 0; jj < size_gtos; jj++)
        {
            tmp += gtos_c(ss).coeff(ii) * gtos_c(ss).coeff(jj) * auxiliary_1e(2+2*gtos_c(ss).l, gtos_c(ss).exp_a(ii)+gtos_c(ss).exp_a(jj));
        }
        gtos_c(ss).coeff = gtos_c(ss).coeff / sqrt(tmp);
    }
    
    return;
}


/*
    Evaluate different one-electron integrals 
*/
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
            int_1e(ii,jj) += gtos_c(ii).coeff(mm) * gtos_c(jj).coeff(nn) * int1e_single_gto(gtos_c(ii).l, gtos_c(ii).m, gtos_c(ii).exp_a(mm), gtos_c(jj).l, gtos_c(jj).m, gtos_c(jj).exp_a(nn), intType);
        }
    }

    return int_1e;
}

Matrix<MatrixXd, -1, -1> GTO::get_h2e()
{
    Matrix<MatrixXd, -1, -1> int_2e(size, size);
    #pragma omp parallel for
    for(int ii = 0; ii < size; ii++)
    for(int jj = 0; jj < size; jj++)
    {
        int_2e(ii,jj).resize(size, size);
        for(int kk = 0; kk < size; kk++)
        for(int ll = 0; ll < size; ll++)
        {
            int_2e(ii,jj)(kk,ll) = 0.0;
            if((gtos_c(ii).l+gtos_c(jj).l+gtos_c(kk).l+gtos_c(ll).l) % 2 || (gtos_c(ii).m+gtos_c(jj).m+gtos_c(kk).m+gtos_c(ll).m) % 2 || (gtos_c(ii).l*gtos_c(jj).l*gtos_c(kk).l*gtos_c(ll).l) < 0)
            {
                continue;
            }
            int size_i = gtos_c(ii).coeff.rows(), size_j = gtos_c(jj).coeff.rows(), size_k = gtos_c(kk).coeff.rows(), size_l = gtos_c(ll).coeff.rows();
            for(int ni = 0; ni < size_i; ni++)
            for(int nj = 0; nj < size_j; nj++)
            for(int nk = 0; nk < size_k; nk++)
            for(int nl = 0; nl < size_l; nl++)
            int_2e(ii,jj)(kk,ll) += gtos_c(ii).coeff(ni) * gtos_c(jj).coeff(nj) * gtos_c(kk).coeff(nk) * gtos_c(ll).coeff(nl) * int2e_single_gto(gtos_c(ii).l, gtos_c(ii).m, gtos_c(ii).exp_a(ni), gtos_c(jj).l, gtos_c(jj).m, gtos_c(jj).exp_a(nj), gtos_c(kk).l, gtos_c(kk).m, gtos_c(kk).exp_a(nk), gtos_c(ll).l, gtos_c(ll).m, gtos_c(ll).exp_a(nl));
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
    else    return factorial(n)/2.0/pow(a,n+1);
}


double GTO::int1e_single_gto(const int& l1, const int& m1, const double& a1, const int& l2, const int& m2, const double& a2, const string& integralTYPE)
{
    if(l1 != l2 || m1 != m2) return 0.0;
    else if(integralTYPE == "overlap")  return auxiliary_1e(2 + 2*l1, a1 + a2);
    else if(integralTYPE == "nuc_attra")  return -atomN * auxiliary_1e(1 + 2*l1, a1 + a2);
    else if(integralTYPE == "kinetic")  return a2 * (2*l1 + 3) * auxiliary_1e(2 + 2*l1, a1 + a2) - 2 * a2 * a2 * auxiliary_1e(4 + 2*l1, a1 + a2);
    else if(integralTYPE == "h1e")  return a2 * (2*l1 + 3) * auxiliary_1e(2 + 2*l1, a1 + a2) - 2 * a2 * a2 * auxiliary_1e(4 + 2*l1, a1 + a2) + -atomN * auxiliary_1e(1 + 2*l1, a1 + a2);
    else if(integralTYPE == "p.Vp")  return ((2*l1 + 1) * l1 * auxiliary_1e(2*l1-1, a1 + a2) - 2*l1*(a1+a2)*auxiliary_1e(2*l1 + 1, a1 + a2) + 4*a1*a2 * auxiliary_1e(2*l1 + 3, a1 + a2)) * -atomN;
    else
    {
        cout << "ERROR: get_h1e is called for undefined type of integrals!" << endl;
        exit(99);
    }
}


/*
    auxiliary_2e_0_r is to evaluate \int_0^inf \int_0^r2 r1^l1 r2^l2 exp(-a1 * r1^2) exp(-a2 * r2^2) dr1dr2
*/
double GTO::auxiliary_2e_0_r(const int& l1, const int& l2, const double& a1, const double& a2)
{
    int n1 = l1 / 2;
    if(n1 * 2 == l1)
    {
        cout << "ERROR: When auxiliary_2e_0r is called, l1 must be set to an odd number!" << endl;
        exit(99);
    }
    else
    {
        double tmp = 0.5 / pow(a1, n1+1) * auxiliary_1e(l2, a2);
        for(int kk = 0; kk <= n1; kk++)
        {
            tmp -= 0.5 / factorial(kk) / pow(a1, n1 - kk + 1) * auxiliary_1e(l2 + 2*kk, a1 + a2);
        }
        return tmp * factorial(n1);
    }
    
}

/*
    auxiliary_2e_r_inf is to evaluate \int_0^inf \int_r2^inf r1^l1 r2^l2 exp(-a1 * r1^2) exp(-a2 * r2^2) dr1dr2
*/
double GTO::auxiliary_2e_r_inf(const int& l1, const int& l2, const double& a1, const double& a2)
{
    int n1 = l1 / 2;
    if(n1 * 2 == l1)
    {
        cout << "ERROR: When auxiliary_2e_0r is called, l1 must be set to an odd number!" << endl;
        exit(99);
    }
    else
    {
        double tmp = 0.0;
        for(int kk = 0; kk <= n1; kk++)
        {
            tmp += 0.5 / factorial(kk) / pow(a1, n1 - kk + 1) * auxiliary_1e(l2 + 2*kk, a1 + a2);
        }
        return tmp * factorial(n1);
    }
    
}

double GTO::int2e_single_gto(const int& l1, const int& m1, const double& a1, const int& l2, const int& m2, const double& a2, const int& l3, const int& m3, const double& a3, const int& l4, const int& m4, const double& a4)
{
    if((l1+l2+l3+l4) % 2 || (m1+m2+m3+m4) % 2 || m1 * m2 * m3 * m4 < 0) return 0.0;
    else
    {
        double result = 0.0, radial, angular;
        
        for(int ll = min(l1 + l2, l3 + l4); ll >= 0; ll = ll - 2)
        {
            if((l1 + l2 + 2 + ll) % 2)
            {
                radial = auxiliary_2e_0_r(l1 + l2 + 2 + ll, l3 + l4 + 1 - ll, a1 + a2, a3 + a4)
                        + auxiliary_2e_0_r(l3 + l4 + 2 + ll, l1 + l2 + 1 - ll, a3 + a4, a1 + a2);
            }
            else
            {
                radial = auxiliary_2e_r_inf(l3 + l4 + 1 - ll, l1 + l2 + 2 + ll, a3 + a4, a1 + a2)
                        + auxiliary_2e_r_inf(l1 + l2 + 1 - ll, l3 + l4 + 2 + ll, a1 + a2, a3 + a4);
            }

            angular = 0.0;
            for(int mm = -ll; mm <= ll; mm++)
            {
                double tmp = 0.0;
                for(int m_i = -abs(m1); m_i <= abs(m1); m_i+=2*abs(m1))
                {
                    for(int m_j = -abs(m2); m_j <= abs(m2); m_j+=2*abs(m2))
                    {
                        for(int m_k = -abs(m3); m_k <= abs(m3); m_k+=2*abs(m3))
                        {
                            for(int m_l = -abs(m4); m_l <= abs(m4); m_l+=2*abs(m4))
                            {
                                if(m_i + m_j - mm != 0 || m_k + m_l + mm != 0)
                                {
                                    tmp += 0.0;
                                }
                                else
                                {
                                    tmp += real(U_SH_trans(m1, m_i) * U_SH_trans(m2, m_j) * U_SH_trans(m3, m_k) * U_SH_trans(m4, m_l))
                                            * wigner_3j(l1, l2, ll, m_i, m_j, -mm) * wigner_3j(l3, l4, ll, m_k, m_l, mm);
                                }
                                if(m_l == 0) break;
                            }
                            if(m_k == 0) break;
                        }
                        if(m_j == 0) break;
                    }
                    if(m_i == 0) break;
                }
                angular += tmp * pow(-1, mm) * sqrt((2.0 * l1 + 1.0)*(2.0 * l2 + 1.0)*(2.0 * l3 + 1.0)*(2.0 * l4 + 1.0))
                            * wigner_3j_zeroM(l1, l2, ll) * wigner_3j_zeroM(l3, l4, ll);
            }

            
            result += radial * angular;
        }
        return result;
    }    
}

