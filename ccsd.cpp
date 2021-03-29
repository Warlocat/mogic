#include"ccsd.h"
#include<omp.h>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>


/*
    Constructor for RHF
*/
CCSD::CCSD(const int& n_occ_, const int& n_vir_, const VectorXd& h2e_so_antiSym_, const VectorXd& ene_mo_):
n_occ(n_occ_), n_vir(n_vir_), h2e_so_antiSym(h2e_so_antiSym_)
{
    ene_so.resize(n_occ+n_vir);
    for(int ii = 0; ii < n_occ+n_vir; ii++)
        ene_so(ii) = ene_mo_(ii/2);
    t1.resize(n_occ, n_vir);
    t2.resize(n_occ*n_occ, n_vir*n_vir);
    t1_new.resize(n_occ, n_vir);
    t2_new.resize(n_occ*n_occ, n_vir*n_vir);
    D1.resize(n_occ, n_vir);
    D2.resize(n_occ*n_occ, n_vir*n_vir);

    double ene_mp2 = 0.0;
    for(int ii = 0; ii < n_occ; ii++)
    for(int aa = n_occ; aa < n_occ + n_vir; aa++)
    {   
        t1(ii,aa-n_occ) = 0.0;
        D1(ii,aa-n_occ) = ene_so(ii) - ene_so(aa);
        for(int jj = 0; jj < n_occ; jj++)
        for(int bb = n_occ; bb < n_occ + n_vir; bb++)
        {
            D2(ii*n_occ+jj, (aa-n_occ)*n_vir+bb-n_occ) = ene_so(ii)+ene_so(jj)-ene_so(aa)-ene_so(bb);
            t2(ii*n_occ+jj, (aa-n_occ)*n_vir+bb-n_occ) = h2e_dirac_so(ii,jj,aa,bb)
                                                        / D2(ii*n_occ+jj, (aa-n_occ)*n_vir+bb-n_occ);
            ene_mp2 += 0.25*h2e_dirac_so(ii,jj,aa,bb)*t2(ii*n_occ+jj, (aa-n_occ)*n_vir+bb-n_occ);
        }
    }

    cout << "MP2 energy from CCSD initialization: " << ene_mp2 << endl;
}

/*
    Constructor for UHF
*/
CCSD::CCSD(const int& n_occ_, const int& n_vir_, const VectorXd& h2e_so_antiSym_, const VectorXd& ene_mo_a, const VectorXd& ene_mo_b):
n_occ(n_occ_), n_vir(n_vir_), h2e_so_antiSym(h2e_so_antiSym_)
{
    cout << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << "! CAUTION: CCSD from UHF is only correct when the system has the same number for alpha and beta electrons. !" << endl;
    cout << "!          This module has NOT been finished.                                                              !" << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << endl;
    exit(99);

    ene_so.resize(n_occ+n_vir);
    for(int ii = 0; ii < n_occ+n_vir; ii++)
    {
        if(ii%2 == 0)    
            ene_so(ii) = ene_mo_a(ii/2);
        else
            ene_so(ii) = ene_mo_b(ii/2);
    }
    t1.resize(n_occ, n_vir);
    t2.resize(n_occ*n_occ, n_vir*n_vir);
    t1_new.resize(n_occ, n_vir);
    t2_new.resize(n_occ*n_occ, n_vir*n_vir);
    D1.resize(n_occ, n_vir);
    D2.resize(n_occ*n_occ, n_vir*n_vir);

    double ene_mp2 = 0.0;
    for(int ii = 0; ii < n_occ; ii++)
    for(int aa = n_occ; aa < n_occ + n_vir; aa++)
    {   
        t1(ii,aa-n_occ) = 0.0;
        D1(ii,aa-n_occ) = ene_so(ii) - ene_so(aa);
        for(int jj = 0; jj < n_occ; jj++)
        for(int bb = n_occ; bb < n_occ + n_vir; bb++)
        {
            D2(ii*n_occ+jj, (aa-n_occ)*n_vir+bb-n_occ) = ene_so(ii)+ene_so(jj)-ene_so(aa)-ene_so(bb);
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
    int sign = 1;
    if(ii < jj) sign *= -1;
    if(kk < ll) sign *= -1;
    int ij = max(ii,jj)*(max(ii,jj)+1)/2+min(ii,jj), kl = max(kk,ll)*(max(kk,ll)+1)/2+min(kk,ll);
    int ijkl = max(ij,kl)*(max(ij,kl)+1)/2+min(ij,kl);
    return sign * h2e_so_antiSym(ijkl);
}


void CCSD::memoryAllocation()
{
    double size = (pow(n_occ,2)*pow(n_vir,2)*5 + pow(n_occ,4) + pow(n_vir,4) + pow(n_occ,2)*pow(n_vir,2) + n_vir*n_vir + n_occ*n_occ + n_vir*n_occ + n_occ*n_vir*3) * sizeof(double) / 1024.0 / 1024.0 * (size_DIIS + 1);

    cout << "Memory requirement: " << size << " MB.\n";

    tau = new double***[n_occ];
    tau_tilde = new double***[n_occ];
    for(int ii = 0; ii < n_occ; ii++)
    {   
        tau[ii] = new double**[n_occ];
        tau_tilde[ii] = new double**[n_occ];
        for(int jj = 0; jj < n_occ; jj++)
        {
            tau[ii][jj] = new double*[n_vir];
            tau_tilde[ii][jj] = new double*[n_vir];
            for(int aa = 0; aa < n_vir; aa++)
            {
                tau[ii][jj][aa] = new double[n_vir];
                tau_tilde[ii][jj][aa] = new double[n_vir];
            }
        }
    }

    W_oooo = new double***[n_occ];
    for(int mm = 0; mm < n_occ; mm++)
    {
        W_oooo[mm] = new double**[n_occ];
        for(int nn = 0; nn < n_occ; nn++)
        {
            W_oooo[mm][nn] = new double*[n_occ];
            for(int ii = 0; ii < n_occ; ii++)
                W_oooo[mm][nn][ii] = new double[n_occ];
        }
    }
    W_vvvv = new double***[n_vir];
    for(int aa = 0; aa < n_vir; aa++)
    {
        W_vvvv[aa] = new double**[n_vir];
        for(int bb = 0; bb < n_vir; bb++)
        {
            W_vvvv[aa][bb] = new double*[n_vir];
            for(int ee = 0; ee < n_vir; ee++)
                W_vvvv[aa][bb][ee] = new double[n_vir];
        }
    }
    W_ovvo = new double***[n_occ];
    for(int mm = 0; mm < n_occ; mm++)
    {
        W_ovvo[mm] = new double**[n_vir];
        for(int bb = 0; bb < n_vir; bb++)
        {
            W_ovvo[mm][bb] = new double*[n_vir];
            for(int ee = 0; ee < n_vir; ee++)
                W_ovvo[mm][bb][ee] = new double[n_occ];
        }
    }

    F_vv = new double*[n_vir];
    for(int aa = 0; aa < n_vir; aa++)
        F_vv[aa] = new double[n_vir];
    F_oo = new double*[n_occ];
    for(int mm = 0; mm < n_occ; mm++)
        F_oo[mm] = new double[n_occ];
    F_ov = new double*[n_occ];
    for(int mm = 0; mm < n_occ; mm++)
        F_ov[mm] = new double[n_vir];
}
void CCSD::memoryDeallocation()
{
    for(int ii = 0; ii < n_occ; ii++)
    {   
        for(int jj = 0; jj < n_occ; jj++)
        {
            for(int aa = n_occ; aa < n_vir; aa++)
            {
                delete[] tau[ii][jj][aa];
                delete[] tau_tilde[ii][jj][aa];
            }
            delete[] tau[ii][jj];
            delete[] tau_tilde[ii][jj];
        }
        delete[] tau[ii];
        delete[] tau_tilde[ii];
    }
    delete[] tau;
    delete[] tau_tilde;
    for(int mm = 0; mm < n_occ; mm++)
    {
        for(int nn = 0; nn < n_occ; nn++)
        {    
            for(int ii = 0; ii < n_occ; ii++)
                delete[] W_oooo[mm][nn][ii];
            delete[] W_oooo[mm][nn];
        }
        delete[] W_oooo[mm];
    }
    delete[] W_oooo;    
    for(int aa = 0; aa < n_vir; aa++)
    {
        for(int bb = 0; bb < n_vir; bb++)
        {        
            for(int ee = 0; ee < n_vir; ee++)
                delete[] W_vvvv[aa][bb][ee];;
            delete[]W_vvvv[aa][bb];
        }
        delete[] W_vvvv[aa];
    }
    delete[] W_vvvv;
    for(int mm = 0; mm < n_occ; mm++)
    {
        for(int bb = 0; bb < n_vir; bb++)
        {
            for(int ee = 0; ee < n_vir; ee++)
                delete[] W_ovvo[mm][bb][ee];
            delete[] W_ovvo[mm][bb];
        }
        delete[] W_ovvo[mm];
    }
    delete[] W_ovvo;

    for(int aa = 0; aa < n_vir; aa++)
        delete[] F_vv[aa];
    delete[] F_vv;
    for(int mm = 0; mm < n_occ; mm++)
        delete[] F_oo[mm];
    delete[] F_oo;
    for(int mm = 0; mm < n_occ; mm++)
        delete[] F_ov[mm];
    delete[] F_ov;
}

MatrixXd CCSD::evaluateErrorDIIS(const MatrixXd& t1_old, const MatrixXd& t1_new, const MatrixXd& t2_old, const MatrixXd& t2_new)
{
    MatrixXd err(t1_old.rows()*t1_old.cols() + t2_old.rows()*t2_old.cols(),1);
    for(int ii = 0; ii < t1_old.rows(); ii++)
    for(int jj = 0; jj < t1_old.cols(); jj++)
    {
        err(ii*t1_old.cols() + jj,0) = t1_new(ii,jj) - t1_old(ii,jj);
    }
    for(int ii = 0; ii < t2_old.rows(); ii++)
    for(int jj = 0; jj < t2_old.cols(); jj++)
    {
        err(t1_old.rows()*t1_old.cols() + ii*t2_old.cols() + jj,0) = t2_new(ii,jj) - t2_old(ii,jj);
    }
    return err;
}

void CCSD::evaluate_tau()
{
    for(int ii = 0; ii < n_occ; ii++)
    for(int aa = 0; aa < n_vir; aa++)
    for(int jj = 0; jj < n_occ; jj++)
    for(int bb = 0; bb < n_vir; bb++)
    {
        tau[ii][jj][aa][bb] = t2(ii*n_occ+jj, aa*n_vir+bb)
                                + t1(ii,aa)*t1(jj,bb) - t1(ii,bb)*t1(jj,aa);
        tau_tilde[ii][jj][aa][bb] = t2(ii*n_occ+jj, aa*n_vir+bb)
                                + 0.5*(t1(ii,aa)*t1(jj,bb) - t1(ii,bb)*t1(jj,aa));
    }

    return;
}

void CCSD::evaluate_W_F()
{
    for(int mm = 0; mm < n_occ; mm++)
    for(int nn = 0; nn < n_occ; nn++)
    for(int ii = 0; ii < n_occ; ii++)
    for(int jj = 0; jj < n_occ; jj++)
    {
        W_oooo[mm][nn][ii][jj] = h2e_dirac_so(mm,nn,ii,jj);
        for(int ee = 0; ee < n_vir; ee++)        
        {
            W_oooo[mm][nn][ii][jj] += t1(jj,ee)*h2e_dirac_so(mm,nn,ii,ee+n_occ);
            W_oooo[mm][nn][ii][jj] -= t1(ii,ee)*h2e_dirac_so(mm,nn,jj,ee+n_occ);
            for(int ff = 0; ff < n_vir; ff++)
            {
                W_oooo[mm][nn][ii][jj] += 0.25*tau[ii][jj][ee][ff]*h2e_dirac_so(mm,nn,ee+n_occ,ff+n_occ);
            }
        }
    }
    for(int aa = 0; aa < n_vir; aa++)
    for(int bb = 0; bb < n_vir; bb++)
    for(int ee = 0; ee < n_vir; ee++)
    for(int ff = 0; ff < n_vir; ff++)
    {
        W_vvvv[aa][bb][ee][ff] = h2e_dirac_so(aa+n_occ,bb+n_occ,ee+n_occ,ff+n_occ);
        for(int mm = 0; mm < n_occ; mm++)
        {
            W_vvvv[aa][bb][ee][ff] -= t1(mm,bb)*h2e_dirac_so(n_occ+aa,mm,n_occ+ee,n_occ+ff); 
            W_vvvv[aa][bb][ee][ff] += t1(mm,aa)*h2e_dirac_so(n_occ+bb,mm,n_occ+ee,n_occ+ff);
            for(int nn = 0; nn < n_occ; nn++)
                W_vvvv[aa][bb][ee][ff] += 0.25*tau[mm][nn][aa][bb]*h2e_dirac_so(mm,nn,n_occ+ee,n_occ+ff);
        }
    }
    for(int mm = 0; mm < n_occ; mm++)
    for(int bb = 0; bb < n_vir; bb++)
    for(int ee = 0; ee < n_vir; ee++)
    for(int jj = 0; jj < n_occ; jj++)
    {
        W_ovvo[mm][bb][ee][jj] = h2e_dirac_so(mm,n_occ+bb,n_occ+ee,jj);
        for(int ff = 0; ff < n_vir; ff++)
        {
            W_ovvo[mm][bb][ee][jj] += t1(jj,ff)*h2e_dirac_so(mm,n_occ+bb,n_occ+ee,n_occ+ff);
        }
        for(int nn = 0; nn < n_occ; nn++)
        {
            W_ovvo[mm][bb][ee][jj] -= t1(nn,bb)*h2e_dirac_so(mm,nn,n_occ+ee,jj);
            for(int ff = 0; ff < n_vir; ff++)
                W_ovvo[mm][bb][ee][jj] -= (0.5*t2(jj*n_occ+nn,ff*n_vir+bb) + t1(jj,ff)*t1(nn,bb)) * h2e_dirac_so(mm,nn,n_occ+ee,n_occ+ff);
        }
    }

    for(int aa = 0; aa < n_vir; aa++)
    for(int ee = 0; ee < n_vir; ee++)
    {
        F_vv[aa][ee] = 0.0;
        for(int mm = 0; mm < n_occ; mm++)
        for(int ff = 0; ff < n_vir; ff++)
        {
            F_vv[aa][ee] += t1(mm,ff)*h2e_dirac_so(mm,n_occ+aa,n_occ+ff,n_occ+ee);
            for(int nn = 0; nn < n_occ; nn++)
                F_vv[aa][ee] -= 0.5*tau_tilde[mm][nn][aa][ff]*h2e_dirac_so(mm,nn,n_occ+ee,n_occ+ff);
        }
    }
    for(int mm = 0; mm < n_occ; mm++)
    for(int ii = 0; ii < n_occ; ii++)
    {
        F_oo[mm][ii] = 0.0;
        for(int nn = 0; nn < n_occ; nn++)
        for(int ee = 0; ee < n_vir; ee++)
        {
            F_oo[mm][ii] += t1(nn,ee)*h2e_dirac_so(mm,nn,ii,n_occ+ee);
            for(int ff = 0; ff < n_vir; ff++)
                F_oo[mm][ii] += 0.5*tau_tilde[ii][nn][ee][ff]*h2e_dirac_so(mm,nn,n_occ+ee,n_occ+ff);
        }
    }
    for(int mm = 0; mm < n_occ; mm++)
    for(int ee = 0; ee < n_vir; ee++)
    {
        F_ov[mm][ee] = 0.0;
        for(int nn = 0; nn < n_occ; nn++)
        for(int ff = 0; ff < n_vir; ff++)
        {
            F_ov[mm][ee] += t1(nn,ff)*h2e_dirac_so(mm,nn,n_occ+ee,n_occ+ff);
        } 
    }
}

void CCSD::evaluate_t1t2New()
{
    t1_new = t1;
    t2_new = t2;

    for(int ii = 0; ii < n_occ; ii++)
    for(int aa = 0; aa < n_vir; aa++)
    {
        // t1 amplitudes
        t1_new(ii,aa) = 0.0;
        for(int ee = 0; ee < n_vir; ee++)
            t1_new(ii,aa) += t1(ii,ee)*F_vv[aa][ee];
        for(int mm = 0; mm < n_occ; mm++)
            t1_new(ii,aa) -= t1(mm,aa)*F_oo[mm][ii];
        for(int mm = 0; mm < n_occ; mm++)
        for(int ee = 0; ee < n_vir; ee++)
        {
            t1_new(ii,aa) += t2(ii*n_occ+mm,aa*n_vir+ee)*F_ov[mm][ee];
            for(int ff = 0; ff < n_vir; ff++)
                t1_new(ii,aa) -= 0.5*t2(ii*n_occ+mm,ee*n_vir+ff)*h2e_dirac_so(mm,n_occ+aa,n_occ+ee,n_occ+ff);
            for(int nn = 0; nn < n_occ; nn++)
                t1_new(ii,aa) -= 0.5*t2(mm*n_occ+nn,aa*n_vir+ee)*h2e_dirac_so(nn,mm,n_occ+ee,ii);
        }
        for(int nn = 0; nn < n_occ; nn++)
        for(int ff = 0; ff < n_vir; ff++)
        {    
            t1_new(ii,aa) -= t1(nn,ff)*h2e_dirac_so(nn,n_occ+aa,ii,n_occ+ff);
        }
        t1_new(ii,aa) = t1_new(ii,aa) / D1(ii,aa);


        // t2 amplitudes
        for(int jj = 0; jj < n_occ; jj++)
        for(int bb = 0; bb < n_vir; bb++)
        {
            t2_new(ii*n_occ+jj,aa*n_vir+bb) = h2e_dirac_so(ii,jj,aa+n_occ,bb+n_occ);
            for(int ee = 0; ee < n_vir; ee++)
            {
                t2_new(ii*n_occ+jj,aa*n_vir+bb) += t2(ii*n_occ+jj,aa*n_vir+ee) * F_vv[bb][ee];
                t2_new(ii*n_occ+jj,aa*n_vir+bb) -= t2(ii*n_occ+jj,bb*n_vir+ee) * F_vv[aa][ee];
                t2_new(ii*n_occ+jj,aa*n_vir+bb) += t1(ii,ee)*h2e_dirac_so(aa+n_occ,bb+n_occ,ee+n_occ,jj);
                t2_new(ii*n_occ+jj,aa*n_vir+bb) -= t1(jj,ee)*h2e_dirac_so(aa+n_occ,bb+n_occ,ee+n_occ,ii);
                for(int mm = 0; mm < n_occ; mm++)
                {
                    t2_new(ii*n_occ+jj,aa*n_vir+bb) -= 0.5*t2(ii*n_occ+jj,aa*n_vir+ee)*t1(mm,bb)*F_ov[mm][ee];
                    t2_new(ii*n_occ+jj,aa*n_vir+bb) += 0.5*t2(ii*n_occ+jj,bb*n_vir+ee)*t1(mm,aa)*F_ov[mm][ee];

                    t2_new(ii*n_occ+jj,aa*n_vir+bb) += (t2(ii*n_occ+mm,aa*n_vir+ee)*W_ovvo[mm][bb][ee][jj] - t1(ii,ee)*t1(mm,aa)*h2e_dirac_so(mm,bb+n_occ,ee+n_occ,jj));
                    t2_new(ii*n_occ+jj,aa*n_vir+bb) -= (t2(ii*n_occ+mm,bb*n_vir+ee)*W_ovvo[mm][aa][ee][jj] - t1(ii,ee)*t1(mm,bb)*h2e_dirac_so(mm,aa+n_occ,ee+n_occ,jj));
                    t2_new(ii*n_occ+jj,aa*n_vir+bb) -= (t2(jj*n_occ+mm,aa*n_vir+ee)*W_ovvo[mm][bb][ee][ii] - t1(jj,ee)*t1(mm,aa)*h2e_dirac_so(mm,bb+n_occ,ee+n_occ,ii));
                    t2_new(ii*n_occ+jj,aa*n_vir+bb) += (t2(jj*n_occ+mm,bb*n_vir+ee)*W_ovvo[mm][aa][ee][ii] - t1(jj,ee)*t1(mm,bb)*h2e_dirac_so(mm,aa+n_occ,ee+n_occ,ii));
                }
                for(int ff = 0; ff < n_vir; ff++)
                    t2_new(ii*n_occ+jj,aa*n_vir+bb) += 0.5*tau[ii][jj][ee][ff]*W_vvvv[aa][bb][ee][ff];
            }
            for(int mm = 0; mm < n_occ; mm++)
            {
                t2_new(ii*n_occ+jj,aa*n_vir+bb) -= t2(ii*n_occ+mm,aa*n_vir+bb)*F_oo[mm][jj];
                t2_new(ii*n_occ+jj,aa*n_vir+bb) += t2(jj*n_occ+mm,aa*n_vir+bb)*F_oo[mm][ii];
                t2_new(ii*n_occ+jj,aa*n_vir+bb) -= t1(mm,aa)*h2e_dirac_so(mm,bb+n_occ,ii,jj);
                t2_new(ii*n_occ+jj,aa*n_vir+bb) += t1(mm,bb)*h2e_dirac_so(mm,aa+n_occ,ii,jj);
                for(int ee = 0; ee < n_vir; ee++)
                {
                    t2_new(ii*n_occ+jj,aa*n_vir+bb) -= 0.5 * t2(ii*n_occ+mm,aa*n_vir+bb)*t1(jj,ee)*F_ov[mm][ee];
                    t2_new(ii*n_occ+jj,aa*n_vir+bb) += 0.5 * t2(jj*n_occ+mm,aa*n_vir+bb)*t1(ii,ee)*F_ov[mm][ee];
                }
                for(int nn = 0; nn < n_occ; nn++)
                    t2_new(ii*n_occ+jj,aa*n_vir+bb) += 0.5*tau[mm][nn][aa][bb]*W_oooo[mm][nn][ii][jj];
            }                        
            t2_new(ii*n_occ+jj,aa*n_vir+bb) = t2_new(ii*n_occ+jj,aa*n_vir+bb) / D2(ii*n_occ+jj,aa*n_vir+bb);
        }
    }
}

double CCSD::evaluateChange(const MatrixXd& M1, const MatrixXd& M2)
{
    double tmp = 0.0;
    for(int ii = 0; ii < M1.rows(); ii++)
    for(int jj = 0; jj < M1.cols(); jj++)
        tmp += pow((M1(ii,jj) - M2(ii,jj)),2);

    return sqrt(tmp);
}

double CCSD::evaluate_ene_ccsd()
{
    double ene = 0.0;
    for(int ii = 0; ii < n_occ; ii++)
    for(int jj = 0; jj < n_occ; jj++)
    for(int aa = 0; aa < n_vir; aa++)
    for(int bb = 0; bb < n_vir; bb++)
    {
        ene += h2e_dirac_so(ii,jj,aa+n_occ,bb+n_occ)*(0.25*t2(ii*n_occ+jj,aa*n_vir+bb) + 0.5*t1(ii,aa)*t1(jj,bb));
    }

    return ene;
}

void CCSD::runCCSD()
{
    memoryAllocation();
    vector<MatrixXd> error_DIIS, t1_DIIS, t2_DIIS;
    double d_t1 = 10, d_t2 = 10, ene_new, de;
    cout << endl;
    cout << "###########################" << endl;
    cout << "#  Start CCSD iterations  #" << endl;
    cout << "###########################" << endl;
    cout << endl;
    for(int iteration = 1; iteration <= maxIter; iteration++)
    {
        if(iteration <= 2)
        {
            evaluate_tau();
            evaluate_W_F();
            evaluate_t1t2New();
        }
        else
        {
            if(iteration == 3)
            {
                cout << endl;
                cout << "###########################" << endl;
                cout << "#       Start DIIS        #" << endl;
                cout << "###########################" << endl;
                cout << endl;
            } 
            int tmp_size = error_DIIS.size();
            MatrixXd B4DIIS(tmp_size+1,tmp_size+1);
            VectorXd vec_b(tmp_size+1);
            for(int ii = 0; ii < tmp_size; ii++)
            {    
                for(int jj = 0; jj <= ii; jj++)
                {
                    B4DIIS(ii,jj) = (error_DIIS[ii].adjoint()*error_DIIS[jj])(0,0);
                    B4DIIS(jj,ii) = B4DIIS(ii,jj);
                }
                B4DIIS(tmp_size, ii) = -1.0;
                B4DIIS(ii, tmp_size) = -1.0;
                vec_b(ii) = 0.0;
            }
            B4DIIS(tmp_size, tmp_size) = 0.0;
            vec_b(tmp_size) = -1.0;
            VectorXd C = B4DIIS.partialPivLu().solve(vec_b);
            t1 = MatrixXd::Zero(t1.rows(),t1.cols());
            t2 = MatrixXd::Zero(t2.rows(),t2.cols());
            for(int ii = 0; ii < tmp_size; ii++)
            {
                t1 += C(ii) * t1_DIIS[ii];
                t2 += C(ii) * t2_DIIS[ii];
            }
            evaluate_tau();
            evaluate_W_F();
            evaluate_t1t2New();
        }
        d_t1 = evaluateChange(t1,t1_new);
        d_t2 = evaluateChange(t2,t2_new);
        t1 = t1_new;
        t2 = t2_new;
        ene_new = evaluate_ene_ccsd();
        de = ene_new - ene_ccsd;
        cout << "Iter " << iteration << " :\n";
        cout << "E(ccsd)         \td_t1            \td_t2            \tde              \n";
        cout << setprecision(16) << ene_new << "\t" << d_t1 << "\t" << d_t2 << "\t" << de << "\n";
        ene_ccsd = ene_new;

        if(max(d_t1,d_t2) < convControl)
        {
            converged = true;
            cout << endl << "CCSD converges after " << iteration << " iterations." << endl;
            cout << "Final CCSD energy is " << setprecision(16) << ene_ccsd << " hartree." << endl;
            memoryDeallocation();
            break;
        }

        evaluate_tau();
        evaluate_W_F();
        evaluate_t1t2New();
        if(error_DIIS.size() >= size_DIIS)
        {
            error_DIIS.erase(error_DIIS.begin());
            error_DIIS.push_back(evaluateErrorDIIS(t1,t1_new,t2,t2_new));
            t1_DIIS.erase(t1_DIIS.begin());
            t1_DIIS.push_back(t1_new);
            t2_DIIS.erase(t2_DIIS.begin());
            t2_DIIS.push_back(t2_new);
        }
        else
        {
            error_DIIS.push_back(evaluateErrorDIIS(t1,t1_new,t2,t2_new));
            t1_DIIS.push_back(t1_new);
            t2_DIIS.push_back(t2_new);
        }
    }
    if(!converged)
    {
        cout << "ERROR: CCSD did NOT converge!" << endl;
        exit(99);
    }
    
    return;
}

void CCSD::runCCSD_pT()
{
    if(!converged) runCCSD();
    double mem = pow(n_occ,3)*pow(n_vir,3) * 5 *sizeof(double) / 1024.0 / 1024.0;
    cout << endl << "Start CCSD(T) calculations......" << endl;
    cout << "CCSD(T) calculation is called. It requires " << mem << " MB memory." << endl;
    double D3[n_occ][n_occ][n_occ][n_vir][n_vir][n_vir];
    double tmp_c[n_occ][n_occ][n_occ][n_vir][n_vir][n_vir], tmp_d[n_occ][n_occ][n_occ][n_vir][n_vir][n_vir];
    double D3t3c[n_occ][n_occ][n_occ][n_vir][n_vir][n_vir], D3t3d[n_occ][n_occ][n_occ][n_vir][n_vir][n_vir];

    for(int ii = 0; ii < n_occ; ii++)
    for(int jj = 0; jj < n_occ; jj++)
    for(int kk = 0; kk < n_occ; kk++)
    for(int aa = 0; aa < n_vir; aa++)
    for(int bb = 0; bb < n_vir; bb++)
    for(int cc = 0; cc < n_vir; cc++)
    {
        D3[ii][jj][kk][aa][bb][cc] = ene_so(ii) + ene_so(jj) + ene_so(kk)
                                    - ene_so(aa + n_occ) - ene_so(bb + n_occ) - ene_so(cc + n_occ);
        tmp_d[ii][jj][kk][aa][bb][cc] = t1(ii,aa) * h2e_dirac_so(jj,kk,bb+n_occ,cc+n_occ);
        tmp_c[ii][jj][kk][aa][bb][cc] = 0.0;
        for(int ee = 0; ee < n_vir; ee++)
            tmp_c[ii][jj][kk][aa][bb][cc] += t2(jj*n_occ+kk,aa*n_vir+ee) * h2e_dirac_so(ee+n_occ,ii,bb+n_occ,cc+n_occ);
        for(int mm = 0; mm < n_occ; mm++)
            tmp_c[ii][jj][kk][aa][bb][cc] -= t2(ii*n_occ+mm,bb*n_vir+cc) * h2e_dirac_so(mm,aa+n_occ,jj,kk);
    }

    ene_pT = 0.0;
    for(int ii = 0; ii < n_occ; ii++)
    for(int jj = 0; jj < n_occ; jj++)
    for(int kk = 0; kk < n_occ; kk++)
    for(int aa = 0; aa < n_vir; aa++)
    for(int bb = 0; bb < n_vir; bb++)
    for(int cc = 0; cc < n_vir; cc++)
    {
        D3t3c[ii][jj][kk][aa][bb][cc] = tmp_c[ii][jj][kk][aa][bb][cc] - tmp_c[ii][jj][kk][bb][aa][cc] - tmp_c[ii][jj][kk][cc][bb][aa]
                                    - tmp_c[jj][ii][kk][aa][bb][cc] - tmp_c[kk][jj][ii][aa][bb][cc] + tmp_c[jj][ii][kk][bb][aa][cc]
                                    + tmp_c[jj][ii][kk][cc][bb][aa] + tmp_c[kk][jj][ii][bb][aa][cc] + tmp_c[kk][jj][ii][cc][bb][aa];
        D3t3d[ii][jj][kk][aa][bb][cc] = tmp_d[ii][jj][kk][aa][bb][cc] - tmp_d[ii][jj][kk][bb][aa][cc] - tmp_d[ii][jj][kk][cc][bb][aa]
                                    - tmp_d[jj][ii][kk][aa][bb][cc] - tmp_d[kk][jj][ii][aa][bb][cc] + tmp_d[jj][ii][kk][bb][aa][cc]
                                    + tmp_d[jj][ii][kk][cc][bb][aa] + tmp_d[kk][jj][ii][bb][aa][cc] + tmp_d[kk][jj][ii][cc][bb][aa];
        ene_pT += 1.0/36.0 * D3t3c[ii][jj][kk][aa][bb][cc] * (D3t3c[ii][jj][kk][aa][bb][cc] + D3t3d[ii][jj][kk][aa][bb][cc]) / D3[ii][jj][kk][aa][bb][cc];
    }

    cout << "E((T)) = " << ene_pT << " hartree." << endl;
    cout << "Final E(CCSD(T)) = " << ene_ccsd + ene_pT << " hartree." << endl; 

    return;
}