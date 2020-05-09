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
          
            gtos_c(int_tmp).gto_list.resize(nGTO);
            gtos_c(int_tmp).coeff.resize(nGTO);

            for(int ii = 0 ; ii < nGTO; ii++)
            {
                ifs >> gtos_c(int_tmp).gto_list(ii).a >> gtos_c(int_tmp).coeff(ii);
                gtos_c(int_tmp).gto_list(ii).l = angularQN;
                gtos_c(int_tmp).gto_list(ii).m = -angularQN;
                gtos_c(int_tmp).coeff(ii) = gtos_c(int_tmp).coeff(ii) / sqrt(auxiliary_1e(2*gtos_c(int_tmp).gto_list(ii).l + 2, 2 * gtos_c(int_tmp).gto_list(ii).a));
            }
            for(int ii = 1; ii < 2*angularQN + 1; ii++)
            {
                gtos_c(int_tmp + ii).gto_list.resize(nGTO);
                gtos_c(int_tmp + ii).coeff.resize(nGTO);
                for(int jj = 0; jj < nGTO; jj++)
                {
                    gtos_c(int_tmp + ii).coeff(jj) = gtos_c(int_tmp).coeff(jj);
                    gtos_c(int_tmp + ii).gto_list(jj).a = gtos_c(int_tmp).gto_list(jj).a;
                    gtos_c(int_tmp + ii).gto_list(jj).l = gtos_c(int_tmp).gto_list(jj).l;
                    gtos_c(int_tmp + ii).gto_list(jj).m = gtos_c(int_tmp).gto_list(jj).m + ii;
                }
            }
            int_tmp += 2*angularQN + 1;
        }   
    ifs.close();
    // for(int ii = 0; ii < size; ii++)
    // for(int jj = 0; jj < gtos_c(ii).coeff.rows(); jj++)
    // {
    //     cout << gtos_c(ii).gto_list(jj).a << "\t" << gtos_c(ii).coeff(jj) << "\t" << gtos_c(ii).gto_list(jj).l << "\t" << gtos_c(ii).gto_list(jj).m << "\n";
    // }
    // exit(99);
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
    #pragma omp parallel for
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
    else    return factorial(n)/2.0/pow(a,n+1);
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

double GTO::int2e_single_gto(const gto_single& gto1, const gto_single& gto2, const gto_single& gto3, const gto_single& gto4)
{
    if((gto1.l + gto2.l + gto3.l + gto4.l) % 2 || (gto1.m + gto2.m + gto3.m + gto4.m) % 2 || gto1.m * gto2.m * gto3.m * gto4.m < 0) return 0.0;
    else
    {
        int l_i = gto1.l, l_j = gto2.l, l_k = gto3.l, l_l = gto4.l, m_i = gto1.m, m_j = gto2.m, m_k = gto3.m, m_l = gto4.m;
        double a_i = gto1.a, a_j = gto2.a, a_k = gto3.a, a_l = gto4.a;
        double result = 0.0, radial, angular;
        
        for(int ll = min(l_i + l_j, l_k + l_l); ll >= 0; ll = ll - 2)
        {
            if((l_i + l_j + 2 + ll) % 2)
            {
                radial = auxiliary_2e_0_r(l_i + l_j + 2 + ll, l_k + l_l + 1 - ll, a_i + a_j, a_k + a_l)
                        + auxiliary_2e_0_r(l_k + l_l + 2 + ll, l_i + l_j + 1 - ll, a_k + a_l, a_i + a_j);
            }
            else
            {
                radial = auxiliary_2e_r_inf(l_k + l_l + 1 - ll, l_i + l_j + 2 + ll, a_k + a_l, a_i + a_j)
                        + auxiliary_2e_r_inf(l_i + l_j + 1 - ll, l_k + l_l + 2 + ll, a_i + a_j, a_k + a_l);
            }

            angular = 0.0;
            for(int mm = -ll; mm <= ll; mm++)
            {
                double tmp = 0.0;
                for(int m1 = -abs(m_i); m1 <= abs(m_i); m1+=2*abs(m_i))
                {
                    for(int m2 = -abs(m_j); m2 <= abs(m_j); m2+=2*abs(m_j))
                    {
                        for(int m3 = -abs(m_k); m3 <= abs(m_k); m3+=2*abs(m_k))
                        {
                            for(int m4 = -abs(m_l); m4 <= abs(m_l); m4+=2*abs(m_l))
                            {
                                if(m1 + m2 - mm != 0 || m3 + m4 + mm != 0)
                                {
                                    tmp += 0.0;
                                }
                                else
                                {
                                    tmp += real(U_SH_trans(m_i, m1) * U_SH_trans(m_j, m2) * U_SH_trans(m_k, m3) * U_SH_trans(m_l, m4))
                                            * wigner_3j(l_i, l_j, ll, m1, m2, -mm) * wigner_3j(l_k, l_l, ll, m3, m4, mm);
                                }
                                if(m4 == 0) break;
                            }
                            if(m3 == 0) break;
                        }
                        if(m2 == 0) break;
                    }
                    if(m1 == 0) break;
                }
                angular += tmp * pow(-1, mm) * sqrt((2.0 * l_i + 1.0)*(2.0 * l_j + 1.0)*(2.0 * l_k + 1.0)*(2.0 * l_l + 1.0))
                            * wigner_3j_zeroM(l_i, l_j, ll) * wigner_3j_zeroM(l_k, l_l, ll);
            }

            
            result += radial * angular;
        }
        return result;
    }    
}

