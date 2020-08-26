#include"ccsd.h"
#include<omp.h>
#include<iostream>
#include<iomanip>
#include<fstream>

CCSD::CCSD(const int& n_occ_, const int& n_vir_, const VectorXd& h2e_mo_so_, const VectorXd& ene_mo_):
n_occ(n_occ_), n_vir(n_vir_), h2e_mo_so(h2e_mo_so_)
{
    ene_mo_so.resize(n_occ+n_vir);
    for(int ii = 0; ii < n_occ+n_vir; ii++)
        ene_mo_so(ii) = ene_mo_(ii/2);
    t1.resize(n_occ, n_vir);
    t2.resize(n_occ*n_occ, n_vir*n_vir);
    D1.resize(n_occ, n_vir);
    D2.resize(n_occ*n_occ, n_vir*n_vir);

    double ene_mp2 = 0.0;
    for(int ii = 0; ii < n_occ; ii++)
    for(int aa = n_occ; aa < n_occ + n_vir; aa++)
    {   
        t1(ii,aa-n_occ) = 0.0;
        D1(ii,aa-n_occ) = ene_mo_so(ii) - ene_mo_so(aa);
        for(int jj = 0; jj < n_occ; jj++)
        for(int bb = n_occ; bb < n_occ + n_vir; bb++)
        {
            D2(ii*n_occ+jj, (aa-n_occ)*n_vir+bb-n_occ) = ene_mo_so(ii)+ene_mo_so(jj)-ene_mo_so(aa)-ene_mo_so(bb);
            t2(ii*n_occ+jj, (aa-n_occ)*n_vir+bb-n_occ) = h2e_dirac_so(ii,jj,aa,bb)
                                                        / D2(ii*n_occ+jj, (aa-n_occ)*n_vir+bb-n_occ);
            ene_mp2 += 0.25*h2e_dirac_so(ii,jj,aa,bb)*t2(ii*n_occ+jj, (aa-n_occ)*n_vir+bb-n_occ);
        }
    }

    cout << "MP2 energy from CCSD initialization: " << ene_mp2 << endl;
}

CCSD::~CCSD()
{
}


inline double CCSD::h2e_dirac_so(const int& ii, const int& jj, const int& kk, const int& ll)
{
    int ik = max(ii,kk)*(max(ii,kk)+1)/2+min(ii,kk), il = max(ii,ll)*(max(ii,ll)+1)/2+min(ii,ll);
    int jk = max(jj,kk)*(max(jj,kk)+1)/2+min(jj,kk), jl = max(jj,ll)*(max(jj,ll)+1)/2+min(jj,ll);
    int ikjl = max(ik,jl)*(max(ik,jl)+1)/2+min(ik,jl), iljk=max(il,jk)*(max(il,jk)+1)/2+min(il,jk);
    return h2e_mo_so(ikjl) - h2e_mo_so(iljk);
}




void CCSD::evaluate_tau(double tau[][][][], double tau_tilde[][][][], const MatrixXd& t1, const MatrixXd& t2)
{
    for(int ii = 0; ii < n_occ; ii++)
    for(int aa = n_occ; aa < n_occ + n_vir; aa++)
    for(int jj = 0; jj < n_occ; jj++)
    for(int bb = n_occ; bb < n_occ + n_vir; bb++)
    {
        tau[ii][jj][aa-n_occ][bb-n_occ] = t2(ii*n_occ+jj, (aa-n_occ)*n_vir+bb-n_occ)
                                + t1(ii,aa-n_occ)*t1(jj,bb-n_occ) - t1(ii,bb-n_occ)*t1(jj,aa-n_occ);
        tau_tilde[ii][jj][aa-n_occ][bb-n_occ] = t2(ii*n_occ+jj, (aa-n_occ)*n_vir+bb-n_occ)
                                + 0.5*(t1(ii,aa-n_occ)*t1(jj,bb-n_occ) - t1(ii,bb-n_occ)*t1(jj,aa-n_occ));
    }

    return;
}



void CCSD::runCCSD()
{
    tau = new double[n_occ];
    // double tau[n_occ][n_occ][n_vir][n_vir], tau_tilde[n_occ][n_occ][n_vir][n_vir];
    for(int ii = 0; ii < n_occ; ii++)
    for(int aa = n_occ; aa < n_occ + n_vir; aa++)
    {
        denominatorD1[ii][aa - n_occ] = ene_mo_so(ii) - ene_mo_so(aa);
        for(int jj = 0; jj < n_occ; jj++)
        for(int bb = n_occ; bb < n_occ + n_vir; bb++)
        {
            denominatorD2[ii][jj][aa-n_occ][bb-n_occ] 
                = ene_mo_so(ii)+ene_mo_so(jj)-ene_mo_so(aa)-ene_mo_so(bb);
        }
    }
}