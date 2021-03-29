#include"cis.h"
#include<omp.h>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>


/*
    Constructor for RHF
*/
CIS::CIS(const int& n_occ_, const int& n_vir_, const VectorXd& h2e_so_antiSym_, const VectorXd& ene_mo_):
n_occ(n_occ_), n_vir(n_vir_), h2e_so_antiSym(h2e_so_antiSym_)
{
    ene_so.resize(n_occ+n_vir);
    for(int ii = 0; ii < n_occ+n_vir; ii++)
        ene_so(ii) = ene_mo_(ii/2);
}

/*
    Constructor for UHF
*/
CIS::CIS(const int& n_occ_, const int& n_vir_, const VectorXd& h2e_so_antiSym_, const VectorXd& ene_mo_a, const VectorXd& ene_mo_b):
n_occ(n_occ_), n_vir(n_vir_), h2e_so_antiSym(h2e_so_antiSym_)
{
    cout << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << "! CAUTION: CIS from UHF is only correct when the system has the same number for alpha and beta electrons.  !" << endl;
    cout << "!          This module has NOT been finished.                                                              !" << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << endl;
    exit(99);
}

CIS::~CIS()
{
}


/*
    H_{ia,jb} = A_{ia,jb} = f_{ab}\delta_{ij} - f_{ij}\delta_{ab} - <aj||bi>
*/
double CIS::evaluateH_CIS(const int& i, const int& j, const int& a, const int& b)
{
    int aj = a*(a+1)/2+j, bi = b*(b+1)/2+i;
    int ajbi = max(aj,bi)*(max(aj,bi)+1)/2+min(aj,bi);
    double e = -h2e_so_antiSym(ajbi);
    if(i == j && a == b)    e += ene_so(a) - ene_so(i);
    return e;
}

/*
    B_{ia,jb} = <ab||ij>
*/
double CIS::evaluateH_TDHF_B(const int& i, const int& j, const int& a, const int& b)
{
    int sign = 1;
    if(a < b) sign *= -1;
    if(i < j) sign *= -1;
    int ab = max(a,b)*(max(a,b)+1)/2+min(a,b), ij = max(i,j)*(max(i,j)+1)/2+min(i,j);
    int abij = max(ab,ij)*(max(ab,ij)+1)/2+min(ab,ij);
    return sign * h2e_so_antiSym(abij);
}

void CIS::runCIS(const string& alg)
{
    MatrixXd H_CIS;
    if(alg == "explicit_D_L" || alg == "explicit")
    {
        H_CIS.resize(n_occ*n_vir, n_occ*n_vir);
        for(int ii = 0; ii < n_occ; ii++)
        for(int aa = 0; aa < n_vir; aa++)
        for(int jj = 0; jj < n_occ; jj++)
        for(int bb = 0; bb < n_vir; bb++)
        {
            int ia = ii*n_vir+aa, jb = jj*n_vir+bb;
            if(ia < jb) continue;
            else
            {
                H_CIS(ii*n_vir+aa, jj*n_vir+bb) = evaluateH_CIS(ii, jj, aa+n_occ, bb+n_occ);
                H_CIS(jj*n_vir+bb, ii*n_vir+aa) = H_CIS(ii*n_vir+aa, jj*n_vir+bb);
            }
        }
        if(alg == "explicit")
        {
            SelfAdjointEigenSolver<MatrixXd> solver(H_CIS); 
            ene_CIS = solver.eigenvalues();
            CIScoeff = solver.eigenvectors();
        }
        else
        {
            cout << "D-L for CIS has NOT been implemented." << endl;
            exit(99);
        }
    }
    else
    {
        cout << "Implicit construction of H_CIS has NOT been implemented." << endl;
        exit(99);
    }
}

void CIS::runTDHF(const string& alg)
{
    MatrixXd A_TDHF, B_TDHF;
    if(alg == "explicit_D_L" || alg == "explicit")
    {
        A_TDHF.resize(n_occ*n_vir, n_occ*n_vir);
        B_TDHF.resize(n_occ*n_vir, n_occ*n_vir);
        for(int ii = 0; ii < n_occ; ii++)
        for(int jj = 0; jj < n_occ; jj++)
        for(int aa = 0; aa < n_vir; aa++)
        for(int bb = 0; bb < n_vir; bb++)
        {
            int ia = ii*n_vir+aa, jb = jj*n_vir+bb;
            if(ia < jb) continue;
            else
            {
                A_TDHF(ii*n_vir+aa, jj*n_vir+bb) = evaluateH_CIS(ii, jj, aa+n_occ, bb+n_occ);
                B_TDHF(ii*n_vir+aa, jj*n_vir+bb) = evaluateH_TDHF_B(ii, jj, aa+n_occ, bb+n_occ);
                A_TDHF(jj*n_vir+bb, ii*n_vir+aa) = A_TDHF(ii*n_vir+aa, jj*n_vir+bb);
                B_TDHF(jj*n_vir+bb, ii*n_vir+aa) = B_TDHF(ii*n_vir+aa, jj*n_vir+bb);
            }
        }
        if(alg == "explicit")
        {
            SelfAdjointEigenSolver<MatrixXd> solver((A_TDHF + B_TDHF) * (A_TDHF - B_TDHF));
            ene_CIS = solver.eigenvalues();
            for(int ii = 0; ii < ene_CIS.rows(); ii++)  ene_CIS(ii) = sqrt(ene_CIS(ii));
            CIScoeff = solver.eigenvectors();
        }
        else
        {
            cout << "D-L for CIS has NOT been implemented." << endl;
            exit(99);
        }
    }
    else
    {
        cout << "Implicit construction of H_CIS has NOT been implemented." << endl;
        exit(99);
    }
}