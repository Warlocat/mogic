#include<Eigen/Dense>
#include<string>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<complex>
#include<omp.h>
#include<gsl/gsl_sf_coupling.h>
#include"gto_spinor.h"
using namespace std;
using namespace Eigen;

GTO_SPINOR::GTO_SPINOR(const string& atomName_, const string& basisSet_, const int& charge_, const int& spin_, const bool& uncontracted_):
GTO(atomName_, basisSet_, charge_, spin_, uncontracted_)
{
    size_gtoc_spinor = 2 * size_gtoc;
    size_gtou_spinor = 2 * size_gtou;
}

GTO_SPINOR::~GTO_SPINOR()
{
}


/*
    Evaluate different one-electron integrals 
*/
MatrixXd GTO_SPINOR::get_h1e(const string& intType, const bool& uncontracted_) const
{
    if(intType == "s_p_nuc_s_p" || intType == "i_s_pV_x_p" )
    {
        MatrixXd int_1e;
        int int_tmp = 0;
        if(!uncontracted_)
        {
            int_1e.resize(size_gtoc_spinor, size_gtoc_spinor);
            int_1e = MatrixXd::Zero(size_gtoc_spinor,size_gtoc_spinor);
        }
        else
        {
            cout << size_gtou_spinor << endl;
            int_1e.resize(size_gtou_spinor, size_gtou_spinor);
            int_1e = MatrixXd::Zero(size_gtou_spinor,size_gtou_spinor);
        }

        for(int ishell = 0; ishell < size_shell; ishell++)
        {
            int ll = shell_list(ishell).l;
            int size_gtos = shell_list(ishell).coeff.rows();
            for(int twojj = 2*ll+1; twojj >= 2*ll-1 && twojj > 0; twojj = twojj - 2)
            {
                double kappa = (twojj + 1.0) * (ll - twojj/2.0);
                MatrixXd h1e_single_shell(size_gtos, size_gtos);
                for(int ii = 0; ii < size_gtos; ii++)
                for(int jj = 0; jj < size_gtos; jj++)
                {
                    double a1 = shell_list(ishell).exp_a(ii), a2 = shell_list(ishell).exp_a(jj);
                    if(intType == "s_p_nuc_s_p")    h1e_single_shell(ii,jj) = (pow(ll + kappa + 1.0, 2) * auxiliary_1e(2*ll-1, a1 + a2) - 2.0*(ll + kappa + 1.0)*(a1 + a2)*auxiliary_1e(2*ll + 1, a1 + a2) + 4*a1*a2 * auxiliary_1e(2*ll + 3, a1 + a2)) * -atomNumber;
                    else if(intType == "i_s_pV_x_p" )   h1e_single_shell(ii,jj) = (kappa + 1.0) * auxiliary_1e(2*ll-1, a1 + a2) * -atomNumber;
                    
                    h1e_single_shell(ii,jj) = h1e_single_shell(ii,jj) / shell_list(ishell).norm(ii) / shell_list(ishell).norm(jj);
                }

                if(!uncontracted_)
                {
                    int size_subshell = shell_list(ishell).coeff.cols();
                    MatrixXd int_1e_shell(size_subshell,size_subshell);
                    for(int ii = 0; ii < size_subshell; ii++)
                    for(int jj = 0; jj < size_subshell; jj++)
                    {
                        int_1e_shell(ii,jj) = 0.0;
                        for(int mm = 0; mm < size_gtos; mm++)
                        for(int nn = 0; nn < size_gtos; nn++)
                        {
                            int_1e_shell(ii,jj) += shell_list(ishell).coeff(mm, ii) * shell_list(ishell).coeff(nn, jj) * h1e_single_shell(mm,nn);
                        }
                    }

                    for(int ii = 0; ii < size_subshell; ii++)
                    for(int jj = 0; jj < size_subshell; jj++)
                    for(int kk = 0; kk < twojj+1; kk++)
                    {
                        int_1e(int_tmp + kk + ii * (twojj+1), int_tmp + kk + jj * (twojj+1)) = int_1e_shell(ii,jj);
                    }
                    int_tmp += size_subshell * (twojj+1);
                }
                else
                {
                    for(int ii = 0; ii < size_gtos; ii++)
                    for(int jj = 0; jj < size_gtos; jj++)
                    for(int kk = 0; kk < twojj+1; kk++)
                    {
                        int_1e(int_tmp + kk + ii * (twojj+1), int_tmp + kk + jj * (twojj+1)) = h1e_single_shell(ii,jj);
                    }
                    int_tmp += size_gtos * (twojj+1);
                }
            }
        }

        return int_1e;
    }
    else
    {
        return get_h1e_spin_orbitals(intType, uncontracted_);
    }
}


/*
    Evaluate different one-electron integrals in spin orbital basis,
    i.e. A_so = A_nr \otimes I_2
*/
MatrixXd GTO_SPINOR::get_h1e_spin_orbitals(const string& intType, const bool& uncontracted_) const
{
    MatrixXd int_1e_nr = GTO::get_h1e(intType, uncontracted_);
    int size_tmp = int_1e_nr.rows();
    MatrixXd int_1e_so(2*size_tmp, 2*size_tmp);

    for(int ii = 0; ii < size_tmp; ii++)
    for(int jj = 0; jj < size_tmp; jj++)
    {
        int_1e_so(2*ii,2*jj) = int_1e_nr(ii,jj);
        int_1e_so(2*ii,2*jj+1) = 0.0;
        int_1e_so(2*ii+1,2*jj) = 0.0;
        int_1e_so(2*ii+1,2*jj+1) = int_1e_nr(ii,jj);
    }
    
    return int_1e_so;
}