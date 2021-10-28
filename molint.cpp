#include<Eigen/Dense>
#include<string>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<complex>
#include<omp.h>
#include"molint.h"
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
            cout << "n is " << n << endl;
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
            cout << "n is " << n << endl;
            exit(99);
        }
        return n * factorial(n - 1);
    }
}

//Calculate the complete Gamma function
double Gamma(double z)
{

    const int a = 12;
    static double c_space[12];
    static double *c = NULL;
    int k;
    double accm;

    if ( c == NULL ) {
        double k1_factrl = 1.0; /* (k - 1)!*(-1)^k with 0!==1*/
        c = c_space;
        c[0] = sqrt(2.0*M_PI);
        for(k=1; k < a; k++) {
            c[k] = exp(a-k) * pow(a-k, k-0.5) / k1_factrl;
            k1_factrl *= -k;
        }
    }
    accm = c[0];
    for(k=1; k < a; k++) {
        accm += c[k] / ( z + k );
    }
    accm *= exp(-(z+a)) * pow(z+a, z+0.5); /* Gamma(z+1) */
    return accm/z;
}

//Calculate the Boys function
double Boys(double x, int n)
{
    if(abs(x)<1e-8) return (1.0/(1.0+2.0 * (double) n));
    else return 0.5* pow(x,-0.5-n) * (gsl_sf_gamma(0.5 + n) - gsl_sf_gamma_inc(0.5+n,x));
}


MOLINT::MOLINT(const MOL& molecule, const bool& unc_):
Natom(molecule.Natom), uncontracted(unc_), atomList(molecule.atomList), coordList(molecule.coordCart), basisList(molecule.basisSet)
{
    readBasis();
    evaluate_l_xyz_coeff();
}

MOLINT::~MOLINT()
{
}


double MOLINT::recurrence_E(const int& t, const int& i, const int& j, const double& Ax, const double& Bx, const double& a, const double& b) const
{
    double p = a+b, mu = a*b/(a+b), Px = (a*Ax+b*Bx)/p;
    double XAB = Ax-Bx, XPA = Px-Ax, XPB = Px-Bx;
    if(i < 0 || j < 0 || t < 0 || t > (i+j)) return 0.0;
    else if(i == 0 && j==0 && t==0)
    {
        return exp(-mu*XAB*XAB);
    }
    else if(j == 0 && t != 0)
    {
        return pow(2.0*p,-t)*factorial(i)/factorial(t)/factorial(i-t)*recurrence_E(0,i-t,0,Ax,Bx,a,b);
    }
    else if(i == 0 && t != 0)
    {
        return pow(2.0*p,-t)*factorial(j)/factorial(t)/factorial(j-t)*recurrence_E(0,0,j-t,Ax,Bx,a,b);
    }
    else if(j == 0)
    {
        return 0.5/p*recurrence_E(t-1,i-1,j,Ax,Bx,a,b) + XPA*recurrence_E(t,i-1,j,Ax,Bx,a,b) + (t+1)*recurrence_E(t+1,i-1,j,Ax,Bx,a,b);
    }
    else
    {
        return 0.5/p*recurrence_E(t-1,i,j-1,Ax,Bx,a,b) + XPB*recurrence_E(t,i,j-1,Ax,Bx,a,b) + (t+1)*recurrence_E(t+1,i,j-1,Ax,Bx,a,b);
    }
}

double MOLINT::recurrence_R(const int& n, const int& t, const int& u, const int& v, const double& p, const Vector3d& Rpc) const
{
    if(t == 0 && u == 0 && v == 0)  return pow(-2.0*p,n) * Boys(p*pow(Rpc.norm(),2), n);
    else if(t == 0 && u == 0)
    {
        if(v > 1) return Rpc(2) * recurrence_R(n+1,t,u,v-1,p,Rpc) + (v-1) * recurrence_R(n+1,t,u,v-2,p,Rpc);
        else return Rpc(2) * recurrence_R(n+1,t,u,v-1,p,Rpc);
    }
    else if(t == 0)
    {
        if(u > 1) return Rpc(1) * recurrence_R(n+1,t,u-1,v,p,Rpc) + (u-1) * recurrence_R(n+1,t,u-2,v,p,Rpc);
        else return Rpc(1) * recurrence_R(n+1,t,u-1,v,p,Rpc);
    }
    else
    {
        if(t > 1) return Rpc(0) * recurrence_R(n+1,t-1,u,v,p,Rpc) + (t-1) * recurrence_R(n+1,t-2,u,v,p,Rpc);
        else return Rpc(0) * recurrence_R(n+1,t-1,u,v,p,Rpc);
    }
}



double MOLINT::overlap_xyz(const int& lx1, const int& ly1, const int& lz1, const double& a1, const Vector3d& X1, const int& lx2, const int& ly2, const int& lz2, const double& a2, const Vector3d& X2) const
{
    double AX = X1(0), AY = X1(1), AZ = X1(2), BX = X2(0), BY = X2(1), BZ = X2(2);
    return pow(M_PI/(a1+a2),1.5) * recurrence_E(0,lx1,lx2,AX,BX,a1,a2) * recurrence_E(0,ly1,ly2,AY,BY,a1,a2) * recurrence_E(0,lz1,lz2,AZ,BZ,a1,a2);
}

double MOLINT::kinetic_xyz(const int& lx1, const int& ly1, const int& lz1, const double& a1, const Vector3d& X1, const int& lx2, const int& ly2, const int& lz2, const double& a2, const Vector3d& X2) const
{
    double kinetic = 4.0*a2*a2*(overlap_xyz(lx1,ly1,lz1,a1,X1,lx2+2,ly2,lz2,a2,X2)
                            + overlap_xyz(lx1,ly1,lz1,a1,X1,lx2,ly2+2,lz2,a2,X2)
                            + overlap_xyz(lx1,ly1,lz1,a1,X1,lx2,ly2,lz2+2,a2,X2));
    kinetic += -2.0*a2*(2*(lx2+ly2+lz2)+3.0) * overlap_xyz(lx1,ly1,lz1,a1,X1,lx2,ly2,lz2,a2,X2);
    if(lx2 >=2) kinetic += lx2*(lx2-1) * overlap_xyz(lx1,ly1,lz1,a1,X1,lx2-2,ly2,lz2,a2,X2);
    if(ly2 >=2) kinetic += ly2*(ly2-1) * overlap_xyz(lx1,ly1,lz1,a1,X1,lx2,ly2-2,lz2,a2,X2);
    if(lz2 >=2) kinetic += lz2*(lz2-1) * overlap_xyz(lx1,ly1,lz1,a1,X1,lx2,ly2,lz2-2,a2,X2);

    return kinetic * -0.5;
}

double MOLINT::nucV_xyz(const int& lx1, const int& ly1, const int& lz1, const double& a1, const Vector3d& X1, const int& lx2, const int& ly2, const int& lz2, const double& a2, const Vector3d& X2) const
{
    double p = a1+a2;
    Vector3d XP = (a1 * X1 + a2 * X2) / p;
    double nucV = 0.0;
    for(int ia = 0; ia < Natom; ia++)
    {
        Vector3d Rpc = XP - coordList(ia);

        MatrixXd tmp1(lx1+lx2+1, ly1+ly2+1);
        for(int tt = 0; tt <= lx1 + lx2; tt++)
        for(int uu = 0; uu <= ly1 + ly2; uu++)
        {
            tmp1(tt,uu) = 0.0;
            for(int vv = 0; vv <= lz1 + lz2; vv++)
                tmp1(tt,uu) += recurrence_R(0, tt, uu, vv, p, Rpc) * recurrence_E(vv,lz1,lz2,X1(2),X2(2),a1,a2);
        }
        VectorXd tmp2(lx1+lx2+1);
        for(int tt = 0; tt <= lx1 + lx2; tt++)
        {
            tmp2(tt) = 0.0;
            for(int uu = 0; uu <= ly1 + ly2; uu++)
                tmp2(tt) += tmp1(tt, uu) * recurrence_E(uu,ly1,ly2,X1(1),X2(1),a1,a2);
        }
        double tmp = 0.0;
        for(int tt = 0; tt <= lx1 + lx2; tt++)
            tmp += tmp2(tt) * recurrence_E(tt,lx1,lx2,X1(0),X2(0),a1,a2);

        nucV += -atomList(ia) * tmp;
    }
    nucV = nucV * 2.0 * M_PI / p;
    return nucV;
}

double MOLINT::eri_xyz(const int& lx1, const int& ly1, const int& lz1, const double& a1, const Vector3d& X1, const int& lx2, const int& ly2, const int& lz2, const double& a2, const Vector3d& X2, const int& lx3, const int& ly3, const int& lz3, const double& a3, const Vector3d& X3, const int& lx4, const int& ly4, const int& lz4, const double& a4, const Vector3d& X4) const
{
    double p = a1+a2, q = a3+a4, alpha = p*q/(p+q);
    Vector3d Rpq = (a1 * X1 + a2 * X2) / p - (a3 * X3 + a4 * X4) / q;


    // double tmp1[lx1+lx2+1][ly1+ly2+1][lz1+lz2+1][lx3+lx4+1][ly3+ly4+1];
    // for(int tt = 0; tt <= lx1 + lx2; tt++)
    // for(int uu = 0; uu <= ly1 + ly2; uu++)
    // for(int vv = 0; vv <= lz1 + lz2; vv++)
    // for(int ll = 0; ll <= lx3 + lx4; ll++)
    // for(int mm = 0; mm <= ly3 + ly4; mm++)
    // {
    //     tmp1[tt][uu][vv][ll][mm] = 0.0;
    //     for(int nn = 0; nn <= lz3 + lz4; nn++)
    //     {
    //         tmp1[tt][uu][vv][ll][mm] += pow(-1,ll+mm+nn) * recurrence_R(0, tt+ll, uu+mm, vv+nn, alpha, Rpq) * recurrence_E(nn,lz3,lz4,X3(2),X4(2),a3,a4);
    //     }
    // }
    // double tmp2[lx1+lx2+1][ly1+ly2+1][lz1+lz2+1][lx3+lx4+1];
    // for(int tt = 0; tt <= lx1 + lx2; tt++)
    // for(int uu = 0; uu <= ly1 + ly2; uu++)
    // for(int vv = 0; vv <= lz1 + lz2; vv++)
    // for(int ll = 0; ll <= lx3 + lx4; ll++)
    // {
    //     tmp2[tt][uu][vv][ll] = 0.0;
    //     for(int mm = 0; mm <= ly3 + ly4; mm++)
    //     {
    //         tmp2[tt][uu][vv][ll] += tmp1[tt][uu][vv][ll][mm] * recurrence_E(mm,ly3,ly4,X3(1),X4(1),a3,a4);
    //     }
    // }
    // double tmp3[lx1+lx2+1][ly1+ly2+1][lz1+lz2+1];
    // for(int tt = 0; tt <= lx1 + lx2; tt++)
    // for(int uu = 0; uu <= ly1 + ly2; uu++)
    // for(int vv = 0; vv <= lz1 + lz2; vv++)
    // {
    //     tmp3[tt][uu][vv] = 0.0;
    //     for(int ll = 0; ll <= lx3 + lx4; ll++)
    //     {
    //         tmp3[tt][uu][vv] += tmp2[tt][uu][vv][ll] * recurrence_E(ll,lx3,lx4,X3(0),X4(0),a3,a4);
    //     }
    // }
    // double tmp4[lx1+lx2+1][ly1+ly2+1];
    // for(int tt = 0; tt <= lx1 + lx2; tt++)
    // for(int uu = 0; uu <= ly1 + ly2; uu++)
    // {
    //     tmp4[tt][uu] = 0.0;
    //     for(int vv = 0; vv <= lz1 + lz2; vv++)
    //     {
    //         tmp4[tt][uu] += tmp3[tt][uu][vv] * recurrence_E(vv,lz1,lz2,X1(2),X2(2),a1,a2);
    //     }
    // }
    // double tmp5[lx1+lx2+1];
    // for(int tt = 0; tt <= lx1 + lx2; tt++)
    // {
    //     tmp5[tt] = 0.0;
    //     for(int uu = 0; uu <= ly1 + ly2; uu++)
    //     {
    //         tmp5[tt] += tmp4[tt][uu] * recurrence_E(uu,ly1,ly2,X1(1),X2(1),a1,a2);
    //     }
    // }
    // double eri = 0.0;
    // for(int tt = 0; tt <= lx1 + lx2; tt++)
    // {
    //     eri += tmp5[tt] * recurrence_E(tt,lx1,lx2,X1(0),X2(0),a1,a2);
    // }

    double eri = 0.0;
    for(int tt = 0; tt <= lx1 + lx2; tt++)
    for(int uu = 0; uu <= ly1 + ly2; uu++)
    for(int vv = 0; vv <= lz1 + lz2; vv++)
    for(int ll = 0; ll <= lx3 + lx4; ll++)
    for(int mm = 0; mm <= ly3 + ly4; mm++)
    for(int nn = 0; nn <= lz3 + lz4; nn++)
    {
        eri += pow(-1,ll+mm+nn) * recurrence_R(0, tt+ll, uu+mm, vv+nn, alpha, Rpq) * recurrence_E(tt,lx1,lx2,X1(0),X2(0),a1,a2) * recurrence_E(uu,ly1,ly2,X1(1),X2(1),a1,a2) * recurrence_E(vv,lz1,lz2,X1(2),X2(2),a1,a2) * recurrence_E(ll,lx3,lx4,X3(0),X4(0),a3,a4) * recurrence_E(mm,ly3,ly4,X3(1),X4(1),a3,a4) * recurrence_E(nn,lz3,lz4,X3(2),X4(2),a3,a4);
    }
    
    eri = eri * 2.0 * pow(M_PI,2.5) / p / q / sqrt(p+q);
    return eri;
}


void MOLINT::readBasis()
{
    Nbasis_con = 0;
    Nbasis_unc = 0;
    string flags;
    Matrix<Matrix<gto_shell,-1,1>,-1,1> shell_list_tmp;
    shell_list_tmp.resize(Natom);
    NbasisList_con.resize(Natom);
    NbasisList_unc.resize(Natom);
    int Nshell = 0;
    for(int na = 0; na < Natom; na++)
    {
        string target = basisList(na);
        ifstream ifs;
        int int_tmp;
    
        ifs.open("GENBAS");
            while (!ifs.eof())
            {
                getline(ifs,flags);
                flags.resize(target.size());
                if(flags == target) 
                {
                    getline(ifs,flags);
                    break;
                }
            }
            if(ifs.eof())
            {
                cout << "ERROR: can not find target basis (" + target + ") in the basis set file (GENBAS)\n";
                exit(99);
            }
            else
            {
                int size_shell;
                ifs >> size_shell;
                Nshell += size_shell;
                MatrixXi orbitalInfo(3,size_shell);
                shell_list_tmp(na).resize(size_shell);

                for(int ii = 0; ii < 3; ii++)
                for(int jj = 0; jj < size_shell; jj++)
                {
                    ifs >> orbitalInfo(ii,jj);
                }
                NbasisList_con(na) = 0;
                NbasisList_unc(na) = 0;
                for(int ishell = 0; ishell < size_shell; ishell++)
                {
                    NbasisList_con(na) += (2 * orbitalInfo(0,ishell) + 1) * orbitalInfo(1,ishell);
                    NbasisList_unc(na) += (2 * orbitalInfo(0,ishell) + 1) * orbitalInfo(2,ishell);
                    
                    shell_list_tmp(na)(ishell).l = orbitalInfo(0,ishell);
                    shell_list_tmp(na)(ishell).coeff.resize(orbitalInfo(2,ishell),orbitalInfo(1,ishell));
                    shell_list_tmp(na)(ishell).exp_a.resize(orbitalInfo(2,ishell));
                    for(int ii = 0; ii < orbitalInfo(2,ishell); ii++)   
                    {    
                        ifs >> shell_list_tmp(na)(ishell).exp_a(ii);
                    }
                    for(int ii = 0; ii < orbitalInfo(2,ishell); ii++)
                    for(int jj = 0; jj < orbitalInfo(1,ishell); jj++)
                    {
                        ifs >> shell_list_tmp(na)(ishell).coeff(ii,jj);
                    }
                }
            }       
        ifs.close();
        Nbasis_con += NbasisList_con(na);
        Nbasis_unc += NbasisList_unc(na);
    }

    int tmp = 0;
    shell_list.resize(Nshell); 
    l_xyz_coeff.resize(Nshell);
    for(int ia = 0; ia < Natom; ia++)
    for(int ii = 0; ii < shell_list_tmp(ia).rows(); ii++)
    {
        shell_list(tmp) = shell_list_tmp(ia)(ii);
        shell_list(tmp).coord = coordList(ia);
        tmp++;
    }
}


inline double MOLINT::norm_factor(const int& lx, const int& ly, const int& lz, const double& a) const
{
    int ll = lx+ly+lz;
    int nn = ll + 1;
    return sqrt(sqrt(a/M_PI) * pow(a,nn) * pow(2.0,nn+1) / double_factorial(2*nn - 1));
}


MatrixXi MOLINT::l_xyz(const int& l) const
{
    MatrixXi xyz;
    switch (l)
    {
    case 0:
        xyz.resize(1,3);
        xyz << 0,0,0;
        return xyz;
        break;
    case 1:
        xyz.resize(3,3);
        xyz << 1,0,0,
               0,1,0,
               0,0,1;
        return xyz;
        break;
    case 2:
        xyz.resize(6,3);
        xyz << 2,0,0,
               1,1,0,
               1,0,1,
               0,2,0,
               0,1,1,
               0,0,2;
        return xyz;
        break;
    case 3:
        xyz.resize(10,3);
        xyz << 3,0,0,
               2,1,0,
               2,0,1,
               1,2,0,
               1,1,1,
               1,0,2,
               0,3,0,
               0,2,1,
               0,1,2,
               0,0,3;
        return xyz;
        break;
    default:
        cout << "Error: L above 4 is not supported." << endl;
        break;
    }
}

void MOLINT::evaluate_l_xyz_coeff()
{
    for(int ishell = 0; ishell < shell_list.rows(); ishell++)
    {
        int ll = shell_list(ishell).l;
        MatrixXd xyz_coeff;
        switch (ll)
        {
        case 0:
            xyz_coeff.resize(1,1);
            xyz_coeff << 1.0 / sqrt(4.0*M_PI);
            break;
        case 1:
            xyz_coeff.resize(3,3);
            xyz_coeff << sqrt(3.0/4.0*M_PI),0,0,
                         0,sqrt(3.0/4.0*M_PI),0,
                         0,0,sqrt(3.0/4.0*M_PI);
            break;
        case 2:
            xyz_coeff.resize(6,5);
            xyz_coeff << 0, 0,-0.25*sqrt(5.0/M_PI), 0, 0.25*sqrt(15.0/M_PI),
                         sqrt(15.0/4.0*M_PI), 0, 0, 0, 0,
                         0, 0, 0, sqrt(15.0/4.0*M_PI), 0,
                         0, 0,-0.25*sqrt(5.0/M_PI), 0,-0.25*sqrt(15.0/M_PI),
                         0, sqrt(15.0/4.0*M_PI), 0, 0, 0,
                         0, 0, 0.5*sqrt(5.0/M_PI), 0, 0;
            break;
        case 3:
            xyz_coeff.resize(10,7);
            xyz_coeff << 0, 0, 0, 0,-1*0.25*sqrt(21.0/2.0/M_PI), 0, 1*0.25*sqrt(35.0/2.0/M_PI),
                         3*0.25*sqrt(35.0/2.0/M_PI), 0,-1*0.25*sqrt(21.0/2.0/M_PI), 0, 0, 0, 0,
                         0, 0, 0,-3*0.25*sqrt(7.0/M_PI), 0, 0.25*sqrt(105.0/M_PI), 0,
                         0, 0, 0, 0,-1*0.25*sqrt(21.0/2.0/M_PI), 0,-3*0.25*sqrt(35.0/2.0/M_PI),
                         0, 0.5*sqrt(105.0/M_PI), 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 4*0.25*sqrt(21.0/2.0/M_PI), 0, 0,
                        -1*0.25*sqrt(35.0/2.0/M_PI), 0,-1*0.25*sqrt(21.0/2.0/M_PI), 0, 0, 0, 0,
                         0, 0, 0,-3*0.25*sqrt(7.0/M_PI), 0,-0.25*sqrt(105.0/M_PI), 0,
                         0, 0, 4*0.25*sqrt(21.0/2.0/M_PI), 0, 0, 0, 0,
                         0, 0, 0, 2*0.25*sqrt(7.0/M_PI), 0, 0, 0;
            break;
        default:
            cout << "Error: L above 3 is not supported." << endl;
            break;
        }
        MatrixXi xyz_i = l_xyz(ll);
        if(uncontracted)  
        {
            l_xyz_coeff(ishell).resize(shell_list(ishell).coeff.rows()*xyz_i.rows(), shell_list(ishell).coeff.rows()*(2*ll+1));
            l_xyz_coeff(ishell) = MatrixXd::Zero(shell_list(ishell).coeff.rows()*xyz_i.rows(), shell_list(ishell).coeff.rows()*(2*ll+1));
            for(int ii = 0; ii < shell_list(ishell).coeff.rows(); ii++)
            {
                for(int iii = 0; iii < xyz_i.rows(); iii++)
                for(int jjj = 0; jjj < 2*ll+1; jjj++)
                {
                    l_xyz_coeff(ishell)(ii*xyz_i.rows()+iii, ii*(2*ll+1)+jjj) = xyz_coeff(iii,jjj) * norm_factor(xyz_i(iii,0), xyz_i(iii,1), xyz_i(iii,2), shell_list(ishell).exp_a(ii));
                }
            }
        }
        else    
        {
            l_xyz_coeff(ishell).resize(shell_list(ishell).coeff.rows()*xyz_i.rows(), shell_list(ishell).coeff.cols()*(2*ll+1));
            l_xyz_coeff(ishell) = MatrixXd::Zero(shell_list(ishell).coeff.rows()*xyz_i.rows(), shell_list(ishell).coeff.cols()*(2*ll+1));
            for(int ii = 0; ii < shell_list(ishell).coeff.rows(); ii++)
            for(int jj = 0; jj < shell_list(ishell).coeff.cols(); jj++)
            {
                for(int iii = 0; iii < xyz_i.rows(); iii++)
                for(int jjj = 0; jjj < 2*ll+1; jjj++)
                {
                    l_xyz_coeff(ishell)(ii*xyz_i.rows()+iii, jj*(2*ll+1)+jjj) = shell_list(ishell).coeff(ii,jj) * xyz_coeff(iii,jjj) * norm_factor(xyz_i(iii,0), xyz_i(iii,1), xyz_i(iii,2), shell_list(ishell).exp_a(ii));
                }
            }
        }
        
        MatrixXd overlap_xyz_tmp(shell_list(ishell).coeff.rows()*xyz_i.rows(), shell_list(ishell).coeff.rows()*xyz_i.rows());
        Vector3d tmpCoord(0.0,0.0,0.0);
        int tmp_i = 0;
        for(int ii = 0; ii < shell_list(ishell).coeff.rows(); ii++)
        {
            for(int ixyz = 0; ixyz < xyz_i.rows(); ixyz++)
            {
                int tmp_j = 0;
                for(int jj = 0; jj < shell_list(ishell).coeff.rows(); jj++)
                {
                    for(int jxyz = 0; jxyz < xyz_i.rows(); jxyz++)
                    {
                        // double norm = norm_factor(xyz_i(ixyz,0), xyz_i(ixyz,1), xyz_i(ixyz,2), shell_list(ishell).exp_a(ii))*norm_factor(xyz_i(jxyz,0), xyz_i(jxyz,1), xyz_i(jxyz,2), shell_list(ishell).exp_a(jj));
                        overlap_xyz_tmp(tmp_i,tmp_j) = overlap_xyz(xyz_i(ixyz,0), xyz_i(ixyz,1), xyz_i(ixyz,2), shell_list(ishell).exp_a(ii), tmpCoord, xyz_i(jxyz,0), xyz_i(jxyz,1), xyz_i(jxyz,2), shell_list(ishell).exp_a(jj), tmpCoord);
                        tmp_j++;
                    }
                }
                tmp_i++;
            }
        }
        MatrixXd overlap_shell = l_xyz_coeff(ishell).adjoint() * overlap_xyz_tmp * l_xyz_coeff(ishell);
        for(int ii = 0; ii < l_xyz_coeff(ishell).rows(); ii++)
        for(int jj = 0; jj < l_xyz_coeff(ishell).cols(); jj++)
        {
            l_xyz_coeff(ishell)(ii,jj) = l_xyz_coeff(ishell)(ii,jj) / sqrt(overlap_shell(jj,jj));
        }
    }    
}

double MOLINT::get_V_RR() const
{
    double V_RR = 0.0;
    for(int ii = 0; ii < Natom; ii++)
    for(int jj = 0; jj < ii; jj++)
    {
        Vector3d R = coordList(ii) - coordList(jj);
        V_RR += atomList(ii) * atomList(jj) / R.norm();
    }
    return V_RR;
}

MatrixXd MOLINT::get_h1e(const string& intName) const
{
    MatrixXd h1e;
    if(uncontracted)    
    {
        h1e.resize(Nbasis_unc, Nbasis_unc);
        h1e = MatrixXd::Zero(Nbasis_unc, Nbasis_unc);
    }
    else
    {
        h1e.resize(Nbasis_con, Nbasis_con);
        h1e = MatrixXd::Zero(Nbasis_con, Nbasis_con);
    }

    int tmp_ii = 0, size_tmp_ii, size_tmp_jj;
    for(int ishell = 0; ishell < shell_list.rows(); ishell++)
    {
        int tmp_jj = 0;
        int l_i = shell_list(ishell).l;
        for(int jshell = 0; jshell < shell_list.rows(); jshell++)
        {
            int l_j = shell_list(jshell).l;
            MatrixXi xyz_i = l_xyz(l_i), xyz_j = l_xyz(l_j);
            int size_shell_i = shell_list(ishell).coeff.rows(), size_shell_j = shell_list(jshell).coeff.rows();
            MatrixXd h1e_unc_xyz_tmp(size_shell_i*xyz_i.rows(), size_shell_j*xyz_j.rows());
            int tmp1 = 0;
            for(int isubshell = 0; isubshell < size_shell_i; isubshell++)
            for(int ixyz = 0; ixyz < xyz_i.rows(); ixyz++)
            {
                int tmp2 = 0;
                for(int jsubshell = 0; jsubshell < size_shell_j; jsubshell++)
                for(int jxyz = 0; jxyz < xyz_j.rows(); jxyz++)
                {
                    if(intName == "overlap") h1e_unc_xyz_tmp(tmp1,tmp2) = overlap_xyz(xyz_i(ixyz,0), xyz_i(ixyz,1), xyz_i(ixyz,2), shell_list(ishell).exp_a(isubshell), shell_list(ishell).coord, xyz_j(jxyz,0), xyz_j(jxyz,1), xyz_j(jxyz,2), shell_list(jshell).exp_a(jsubshell), shell_list(jshell).coord);
                    else if(intName == "kinetic") h1e_unc_xyz_tmp(tmp1,tmp2) = kinetic_xyz(xyz_i(ixyz,0), xyz_i(ixyz,1), xyz_i(ixyz,2), shell_list(ishell).exp_a(isubshell), shell_list(ishell).coord, xyz_j(jxyz,0), xyz_j(jxyz,1), xyz_j(jxyz,2), shell_list(jshell).exp_a(jsubshell), shell_list(jshell).coord);
                    else if(intName == "nucV") h1e_unc_xyz_tmp(tmp1,tmp2) = nucV_xyz(xyz_i(ixyz,0), xyz_i(ixyz,1), xyz_i(ixyz,2), shell_list(ishell).exp_a(isubshell), shell_list(ishell).coord, xyz_j(jxyz,0), xyz_j(jxyz,1), xyz_j(jxyz,2), shell_list(jshell).exp_a(jsubshell), shell_list(jshell).coord);
                    else
                    {
                        cout << "Error: get_h1e called with an unknown intName." << endl;
                        exit(99);
                    }

                    tmp2++;
                }
                tmp1++;
            }

            if(uncontracted)
            {    
                size_tmp_ii = shell_list(ishell).coeff.rows();
                size_tmp_jj = shell_list(jshell).coeff.rows();
            }
            else
            {
                size_tmp_ii = shell_list(ishell).coeff.cols();
                size_tmp_jj = shell_list(jshell).coeff.cols();
            }
            MatrixXd h1e_tmp = l_xyz_coeff(ishell).adjoint() * h1e_unc_xyz_tmp * l_xyz_coeff(jshell);
            for(int ii = 0; ii < h1e_tmp.rows(); ii++)
            for(int jj = 0; jj < h1e_tmp.cols(); jj++)
            {
                h1e(tmp_ii + ii, tmp_jj + jj) = h1e_tmp(ii,jj);
            }
            tmp_jj += (2*l_j + 1) * size_tmp_jj;
        }
        tmp_ii += (2*l_i + 1) * size_tmp_ii;
    }

    return h1e;
}

VectorXd MOLINT::get_h2e(const string& intName) const
{
    int size_basis;
    VectorXd h2e;
    if(uncontracted)    size_basis = Nbasis_unc;
    else    size_basis = Nbasis_con;
    int tmp_int = size_basis * (size_basis + 1) / 2;
    h2e.resize(tmp_int*(tmp_int+1)/2);
    h2e = VectorXd::Zero(tmp_int*(tmp_int+1)/2);

    int tmp_ii = 0, size_tmp_ii, size_tmp_jj, size_tmp_kk, size_tmp_ll;
    for(int ishell = 0; ishell < shell_list.rows(); ishell++)
    {
        int tmp_jj = 0, l_i = shell_list(ishell).l, size_shell_i = shell_list(ishell).coeff.rows();
        for(int jshell = 0; jshell <= ishell; jshell++)
        {
            int tmp_kk = 0, l_j = shell_list(jshell).l, size_shell_j = shell_list(jshell).coeff.rows();
            for(int kshell = 0; kshell < shell_list.rows(); kshell++)
            {
                int tmp_ll = 0, l_k = shell_list(kshell).l, size_shell_k = shell_list(kshell).coeff.rows();
                for(int lshell = 0; lshell <= kshell; lshell++)
                {
                    if(uncontracted)
                    {    
                        size_tmp_ii = shell_list(ishell).coeff.rows();
                        size_tmp_jj = shell_list(jshell).coeff.rows();
                        size_tmp_kk = shell_list(kshell).coeff.rows();
                        size_tmp_ll = shell_list(lshell).coeff.rows();
                    }
                    else
                    {
                        size_tmp_ii = shell_list(ishell).coeff.cols();
                        size_tmp_jj = shell_list(jshell).coeff.cols();
                        size_tmp_kk = shell_list(kshell).coeff.cols();
                        size_tmp_ll = shell_list(lshell).coeff.cols();
                    }
                    
                    int l_l = shell_list(lshell).l, size_shell_l = shell_list(lshell).coeff.rows();
                    MatrixXi xyz_i = l_xyz(l_i), xyz_j = l_xyz(l_j), xyz_k = l_xyz(l_k), xyz_l = l_xyz(l_l);
                    double h2e_tmp_xyz[size_shell_i*xyz_i.rows()][size_shell_j*xyz_j.rows()][size_shell_k*xyz_k.rows()][size_shell_l*xyz_l.rows()];
                    int tmp1 = 0;
                    for(int isubshell = 0; isubshell < size_shell_i; isubshell++)
                    for(int ixyz = 0; ixyz < xyz_i.rows(); ixyz++)
                    {
                        int tmp2 = 0;
                        for(int jsubshell = 0; jsubshell < size_shell_j; jsubshell++)
                        for(int jxyz = 0; jxyz < xyz_j.rows(); jxyz++)
                        {
                            int tmp3 = 0;
                            for(int ksubshell = 0; ksubshell < size_shell_k; ksubshell++)
                            for(int kxyz = 0; kxyz < xyz_k.rows(); kxyz++)
                            {
                                int tmp4 = 0;
                                for(int lsubshell = 0; lsubshell < size_shell_l; lsubshell++)
                                for(int lxyz = 0; lxyz < xyz_l.rows(); lxyz++)
                                {
                                    if(intName == "eriLLLL") h2e_tmp_xyz[tmp1][tmp2][tmp3][tmp4] = eri_xyz(xyz_i(ixyz,0), xyz_i(ixyz,1), xyz_i(ixyz,2), shell_list(ishell).exp_a(isubshell), shell_list(ishell).coord, xyz_j(jxyz,0), xyz_j(jxyz,1), xyz_j(jxyz,2), shell_list(jshell).exp_a(jsubshell), shell_list(jshell).coord, xyz_k(kxyz,0), xyz_k(kxyz,1), xyz_k(kxyz,2), shell_list(kshell).exp_a(ksubshell), shell_list(kshell).coord, xyz_l(lxyz,0), xyz_l(lxyz,1), xyz_l(lxyz,2), shell_list(lshell).exp_a(lsubshell), shell_list(lshell).coord);
                                    else
                                    {
                                        cout << "Error: get_h2e called with an unknown intName." << endl;
                                    }
                                    tmp4++;
                                }
                                tmp3++;
                            }
                            tmp2++;
                        }
                        tmp1++;
                    }

    
                    for(int ii = 0; ii < l_xyz_coeff(ishell).cols(); ii++)
                    for(int jj = 0; jj < l_xyz_coeff(jshell).cols(); jj++)
                    for(int kk = 0; kk < l_xyz_coeff(kshell).cols(); kk++)
                    for(int ll = 0; ll < l_xyz_coeff(lshell).cols(); ll++)
                    {
                        int ij = max(tmp_ii + ii,tmp_jj + jj) * (max(tmp_ii + ii,tmp_jj + jj) + 1) / 2 + min(tmp_ii + ii,tmp_jj + jj);
                        int kl = max(tmp_kk + kk,tmp_ll + ll) * (max(tmp_kk + kk,tmp_ll + ll) + 1) / 2 + min(tmp_kk + kk,tmp_ll + ll);
                        if(tmp_ii+ii < tmp_jj+jj || tmp_kk+kk < tmp_ll+ll || ij < kl) continue;
                        int ijkl = ij*(ij+1)/2+kl;
                        h2e(ijkl) = 0.0;
                        for(int iii = 0; iii < size_shell_i*xyz_i.rows(); iii++)
                        for(int jjj = 0; jjj < size_shell_j*xyz_j.rows(); jjj++)
                        for(int kkk = 0; kkk < size_shell_k*xyz_k.rows(); kkk++)
                        for(int lll = 0; lll < size_shell_l*xyz_l.rows(); lll++)
                        {
                            h2e(ijkl) += l_xyz_coeff(ishell)(iii,ii)*l_xyz_coeff(jshell)(jjj,jj)*l_xyz_coeff(kshell)(kkk,kk)*l_xyz_coeff(lshell)(lll,ll)*h2e_tmp_xyz[iii][jjj][kkk][lll];
                        }
                    }
                    tmp_ll += (2*l_l+1) * size_tmp_ll;
                }
                tmp_kk += (2*l_k+1) * size_tmp_kk;
            }
            tmp_jj += (2*l_j + 1) * size_tmp_jj;
        }
        tmp_ii += (2*l_i + 1) * size_tmp_ii;
    }

    return h2e;
}




