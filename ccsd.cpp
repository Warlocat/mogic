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
    t1.reset({n_occ, n_vir},0.0);
    t2.reset({n_occ,n_occ,n_vir,n_vir},0.0);
    t1_new.reset({n_occ, n_vir},0.0);
    t2_new.reset({n_occ,n_occ,n_vir,n_vir},0.0);
    D1.reset({n_occ, n_vir},0.0);
    D2.reset({n_occ,n_occ,n_vir,n_vir},0.0);

    double ene_mp2 = 0.0;
    for(int ii = 0; ii < n_occ; ii++)
    for(int aa = n_occ; aa < n_occ + n_vir; aa++)
    {   
        t1(ii,aa-n_occ) = 0.0;
        D1(ii,aa-n_occ) = ene_so(ii) - ene_so(aa);
        for(int jj = 0; jj < n_occ; jj++)
        for(int bb = n_occ; bb < n_occ + n_vir; bb++)
        {
            D2(ii,jj,aa-n_occ,bb-n_occ) = ene_so(ii)+ene_so(jj)-ene_so(aa)-ene_so(bb);
            t2(ii,jj,aa-n_occ,bb-n_occ) = h2e_dirac_so(ii,jj,aa,bb) / D2(ii,jj,aa-n_occ,bb-n_occ);
            ene_mp2 += 0.25*h2e_dirac_so(ii,jj,aa,bb)*t2(ii,jj,aa-n_occ,bb-n_occ);
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
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << "! CAUTION: CCSD from UHF is only correct when the system has !" << endl;
    cout << "! the same number for alpha and beta electrons.              !" << endl;
    cout << "!          This module has NOT been finished.                !" << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << endl;

    ene_so.resize(n_occ+n_vir);
    for(int ii = 0; ii < n_occ+n_vir; ii++)
    {
        if(ii%2 == 0)    
            ene_so(ii) = ene_mo_a(ii/2);
        else
            ene_so(ii) = ene_mo_b(ii/2);
    }
    t1.reset({n_occ, n_vir},0.0);
    t2.reset({n_occ,n_occ,n_vir,n_vir},0.0);
    t1_new.reset({n_occ, n_vir},0.0);
    t2_new.reset({n_occ,n_occ,n_vir,n_vir},0.0);
    D1.reset({n_occ, n_vir},0.0);
    D2.reset({n_occ,n_occ,n_vir,n_vir},0.0);

    double ene_mp2 = 0.0;
    for(int ii = 0; ii < n_occ; ii++)
    for(int aa = n_occ; aa < n_occ + n_vir; aa++)
    {   
        t1(ii,aa-n_occ) = 0.0;
        D1(ii,aa-n_occ) = ene_so(ii) - ene_so(aa);
        for(int jj = 0; jj < n_occ; jj++)
        for(int bb = n_occ; bb < n_occ + n_vir; bb++)
        {
            D2(ii,jj,aa-n_occ,bb-n_occ) = ene_so(ii)+ene_so(jj)-ene_so(aa)-ene_so(bb);
            t2(ii,jj,aa-n_occ,bb-n_occ) = h2e_dirac_so(ii,jj,aa,bb) / D2(ii,jj,aa-n_occ,bb-n_occ);
            ene_mp2 += 0.25*h2e_dirac_so(ii,jj,aa,bb)*t2(ii,jj,aa-n_occ,bb-n_occ);
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
    double size = (pow(n_occ,2)*pow(n_vir,2)*8 + pow(n_occ,4)*2 + pow(n_vir,4)*2 + pow(n_occ,1)*pow(n_vir,3) + pow(n_occ,3)*pow(n_vir,1) + n_vir*n_vir + n_occ*n_occ + n_vir*n_occ + n_occ*n_vir*3) * sizeof(double) / 1024.0 / 1024.0 * (size_DIIS + 1);

    cout << endl << endl << "CCSD memory requirement: " << size << " MB.\n";

    tau.reset({n_occ,n_occ,n_vir,n_vir},0.0);
    tau_tilde.reset({n_occ,n_occ,n_vir,n_vir},0.0);
    W_oooo.reset({n_occ,n_occ,n_occ,n_occ},0.0);
    W_vvvv.reset({n_vir,n_vir,n_vir,n_vir},0.0);
    W_ovvo.reset({n_occ,n_vir,n_vir,n_occ},0.0);

    F_vv.reset({n_vir,n_vir},0.0);
    F_oo.reset({n_occ,n_occ},0.0);
    F_ov.reset({n_occ,n_vir},0.0);

    h2e_oooo.reset({n_occ,n_occ,n_occ,n_occ},0.0);
    h2e_ooov.reset({n_occ,n_occ,n_occ,n_vir},0.0);
    h2e_oovv.reset({n_occ,n_occ,n_vir,n_vir},0.0);
    h2e_ovvo.reset({n_occ,n_vir,n_vir,n_occ},0.0);
    h2e_ovvv.reset({n_occ,n_vir,n_vir,n_vir},0.0);
    h2e_vvvv.reset({n_vir,n_vir,n_vir,n_vir},0.0);
    for(int ii = 0; ii < n_occ; ii++)
    for(int jj = 0; jj < n_occ; jj++)
    for(int kk = 0; kk < n_occ; kk++)
    for(int ll = 0; ll < n_occ; ll++)
        h2e_oooo(ii,jj,kk,ll) = h2e_dirac_so(ii,jj,kk,ll);
    for(int ii = 0; ii < n_occ; ii++)
    for(int jj = 0; jj < n_occ; jj++)
    for(int kk = 0; kk < n_occ; kk++)
    for(int ll = 0; ll < n_vir; ll++)
        h2e_ooov(ii,jj,kk,ll) = h2e_dirac_so(ii,jj,kk,ll+n_occ);
    for(int ii = 0; ii < n_occ; ii++)
    for(int jj = 0; jj < n_occ; jj++)
    for(int kk = 0; kk < n_vir; kk++)
    for(int ll = 0; ll < n_vir; ll++)
        h2e_oovv(ii,jj,kk,ll) = h2e_dirac_so(ii,jj,kk+n_occ,ll+n_occ);
    for(int ii = 0; ii < n_occ; ii++)
    for(int jj = 0; jj < n_vir; jj++)
    for(int kk = 0; kk < n_vir; kk++)
    for(int ll = 0; ll < n_occ; ll++)
        h2e_ovvo(ii,jj,kk,ll) = h2e_dirac_so(ii,jj+n_occ,kk+n_occ,ll);
    for(int ii = 0; ii < n_occ; ii++)
    for(int jj = 0; jj < n_vir; jj++)
    for(int kk = 0; kk < n_vir; kk++)
    for(int ll = 0; ll < n_vir; ll++)
        h2e_ovvv(ii,jj,kk,ll) = h2e_dirac_so(ii,jj+n_occ,kk+n_occ,ll+n_occ);
    for(int ii = 0; ii < n_vir; ii++)
    for(int jj = 0; jj < n_vir; jj++)
    for(int kk = 0; kk < n_vir; kk++)
    for(int ll = 0; ll < n_vir; ll++)
        h2e_vvvv(ii,jj,kk,ll) = h2e_dirac_so(ii+n_occ,jj+n_occ,kk+n_occ,ll+n_occ);
}

void CCSD::memoryDeallocation()
{
    tau.reset();
    tau_tilde.reset();
    W_oooo.reset();
    W_vvvv.reset();
    W_ovvo.reset();
    F_vv.reset();
    F_oo.reset();
    F_ov.reset();
    h2e_oooo.reset();
    h2e_ooov.reset();
    h2e_oovv.reset();
    h2e_ovvo.reset();
    h2e_ovvv.reset();
    h2e_vvvv.reset();
}

tensor<double> CCSD::evaluateErrorDIIS(const tensor<double>& t1_old, const tensor<double>& t1_new, const tensor<double>& t2_old, const tensor<double>& t2_new)
{
    tensor<double> err({t1_old.length(0)*t1_old.length(1) + t2_old.length(0)*t2_old.length(1)*t2_old.length(2)*t2_old.length(3)});
    for(int ii = 0; ii < t1_old.length(0); ii++)
    for(int jj = 0; jj < t1_old.length(1); jj++)
    {
        err(ii*t1_old.length(1) + jj) = t1_new(ii,jj) - t1_old(ii,jj);
    }
    for(int ii = 0; ii < t2_old.length(0); ii++)
    for(int jj = 0; jj < t2_old.length(1); jj++)
    for(int aa = 0; aa < t2_old.length(2); aa++)
    for(int bb = 0; bb < t2_old.length(3); bb++)
    {
        err(t1_old.length(0)*t1_old.length(1) + ii*t2_old.length(1)*t2_old.length(2)*t2_old.length(3) + jj*t2_old.length(2)*t2_old.length(3) + aa*t2_old.length(3) + bb) = t2_new(ii,jj,aa,bb) - t2_old(ii,jj,aa,bb);
    }
    return err;
}

void CCSD::evaluate_tau()
{
    tau = t2;
    tau_tilde = t2;
    mult<double>(1.0,t1,"ia",t1,"jb",1.0,tau,"ijab");
    mult<double>(-1.0,t1,"ib",t1,"ja",1.0,tau,"ijab");
    mult<double>(0.5,t1,"ia",t1,"jb",1.0,tau_tilde,"ijab");
    mult<double>(-0.5,t1,"ib",t1,"ja",1.0,tau_tilde,"ijab");
    return;
}

void CCSD::evaluate_W_F()
{
    W_oooo = h2e_oooo;
    mult<double>(1.0,t1,"je",h2e_ooov,"mnie",1.0,W_oooo,"mnij");
    mult<double>(-1.0,t1,"ie",h2e_ooov,"mnje",1.0,W_oooo,"mnij");
    mult<double>(0.25,tau,"ijef",h2e_oovv,"mnef",1.0,W_oooo,"mnij");

    W_vvvv = h2e_vvvv;
    mult<double>(1.0,t1,"mb",h2e_ovvv,"maef",1.0,W_vvvv,"abef");
    mult<double>(-1.0,t1,"ma",h2e_ovvv,"mbef",1.0,W_vvvv,"abef");
    mult<double>(0.25,tau,"mnab",h2e_oovv,"mnef",1.0,W_vvvv,"abef");

    W_ovvo = h2e_ovvo;
    mult<double>(1.0,t1,"jf",h2e_ovvv,"mbef",1.0,W_ovvo,"mbej");
    mult<double>(1.0,t1,"nb",h2e_ooov,"mnje",1.0,W_ovvo,"mbej");
    tensor<double> tensor_tmp({n_occ,n_occ,n_vir,n_vir});
    tensor_tmp = t2;
    mult<double>(1.0,t1,"jf",t1,"nb",0.5,tensor_tmp,"jnfb");
    mult<double>(-1.0,tensor_tmp,"jnfb",h2e_oovv,"mnef",1.0,W_ovvo,"mbej");


    mult<double>(1.0,t1,"mf",h2e_ovvv,"mafe",0.0,F_vv,"ae");
    mult<double>(-0.5,tau_tilde,"mnaf",h2e_oovv,"mnef",1.0,F_vv,"ae");
    mult<double>(1.0,t1,"ne",h2e_ooov,"mnie",0.0,F_oo,"mi");
    mult<double>(0.5,tau_tilde,"inef",h2e_oovv,"mnef",1.0,F_oo,"mi");
    mult<double>(1.0,t1,"nf",h2e_oovv,"mnef",0.0,F_ov,"me");    
}

void CCSD::evaluate_t1t2New()
{
    // t1 amplitudes
    mult<double>(1.0,t1,"ie",F_vv,"ae",0.0,t1_new,"ia");
    mult<double>(-1.0,t1,"ma",F_oo,"mi",1.0,t1_new,"ia");
    mult<double>(1.0,t2,"imae",F_ov,"me",1.0,t1_new,"ia");
    mult<double>(-0.5,t2,"imef",h2e_ovvv,"maef",1.0,t1_new,"ia");
    mult<double>(0.5,t2,"mnae",h2e_ooov,"nmie",1.0,t1_new,"ia");
    mult<double>(1.0,t1,"nf",h2e_ovvo,"nafi",1.0,t1_new,"ia");
    // t2 amplitudes
    t2_new = h2e_oovv;
    mult<double>(1.0,t2,"ijae",F_vv,"be",1.0,t2_new,"ijab");
    mult<double>(-1.0,t2,"ijbe",F_vv,"ae",1.0,t2_new,"ijab");
    mult<double>(-1.0,t1,"ie",h2e_ovvv,"jeab",1.0,t2_new,"ijab");
    mult<double>(1.0,t1,"je",h2e_ovvv,"ieab",1.0,t2_new,"ijab");
    tensor<double> tensor_tmp({n_vir,n_vir});
    mult<double>(1.0,t1,"mb",F_ov,"me",0.0,tensor_tmp,"be");
    mult<double>(-0.5,t2,"ijae",tensor_tmp,"be",1.0,t2_new,"ijab");
    mult<double>(0.5,t2,"ijbe",tensor_tmp,"ae",1.0,t2_new,"ijab");
    tensor_tmp.reset({n_occ,n_vir,n_occ,n_vir});
    mult<double>(1.0,t1,"ie",t1,"ma",0.0,tensor_tmp,"iema");
    mult<double>(1.0,t2,"imae",W_ovvo,"mbej",1.0,t2_new,"ijab");
    mult<double>(-1.0,t2,"imbe",W_ovvo,"maej",1.0,t2_new,"ijab");
    mult<double>(-1.0,t2,"jmae",W_ovvo,"mbei",1.0,t2_new,"ijab");
    mult<double>(1.0,t2,"jmbe",W_ovvo,"maei",1.0,t2_new,"ijab");
    mult<double>(-1.0,tensor_tmp,"iema",h2e_ovvo,"mbej",1.0,t2_new,"ijab");
    mult<double>(1.0,tensor_tmp,"iemb",h2e_ovvo,"maej",1.0,t2_new,"ijab");
    mult<double>(1.0,tensor_tmp,"jema",h2e_ovvo,"mbei",1.0,t2_new,"ijab");
    mult<double>(-1.0,tensor_tmp,"jemb",h2e_ovvo,"maei",1.0,t2_new,"ijab");
    mult<double>(0.5,tau,"ijef",W_vvvv,"abef",1.0,t2_new,"ijab");
    mult<double>(-1.0,t2,"imab",F_oo,"mj",1.0,t2_new,"ijab");
    mult<double>(1.0,t2,"jmab",F_oo,"mi",1.0,t2_new,"ijab");
    mult<double>(-1.0,t1,"ma",h2e_ooov,"ijmb",1.0,t2_new,"ijab");
    mult<double>(1.0,t1,"mb",h2e_ooov,"ijma",1.0,t2_new,"ijab");
    tensor_tmp.reset({n_occ,n_occ});
    mult<double>(1.0,t1,"je",F_ov,"me",0.0,tensor_tmp,"jm");
    mult<double>(-0.5,t2,"imab",tensor_tmp,"jm",1.0,t2_new,"ijab");
    mult<double>(0.5,t2,"jmab",tensor_tmp,"im",1.0,t2_new,"ijab");
    mult<double>(0.5,tau,"mnab",W_oooo,"mnij",1.0,t2_new,"ijab");
    for(int ii = 0; ii < n_occ; ii++)
    for(int aa = 0; aa < n_vir; aa++)
    {
        t1_new(ii,aa) = t1_new(ii,aa) / D1(ii,aa);
        for(int jj = 0; jj < n_occ; jj++)
        for(int bb = 0; bb < n_vir; bb++)
        {                     
            t2_new(ii,jj,aa,bb) = t2_new(ii,jj,aa,bb) / D2(ii,jj,aa,bb);
        }
    }
    return;
}

double CCSD::evaluateChange_t1(const tensor<double>& M1, const tensor<double>& M2)
{
    double tmp = 0.0;
    for(int ii = 0; ii < M1.length(0); ii++)
    for(int jj = 0; jj < M1.length(1); jj++)
        tmp += pow((M1(ii,jj) - M2(ii,jj)),2);
    return sqrt(tmp);
}
double CCSD::evaluateChange_t2(const tensor<double>& M1, const tensor<double>& M2)
{
    double tmp = 0.0;
    for(int ii = 0; ii < M1.length(0); ii++)
    for(int jj = 0; jj < M1.length(1); jj++)
    for(int aa = 0; aa < M1.length(2); aa++)
    for(int bb = 0; bb < M1.length(3); bb++)
        tmp += pow((M1(ii,jj,aa,bb) - M2(ii,jj,aa,bb)),2);
    return sqrt(tmp);
}

double CCSD::evaluate_ene_ccsd()
{
    double ene = 0.0;
    tensor<double> tmp = t2;
    mult<double>(0.5,t1,"ia",t1,"jb",0.25,tmp,"ijab");
    dot<double>(h2e_oovv,"ijab",tmp,"ijab",ene);

    return ene;
}

void CCSD::runCCSD()
{
    memoryAllocation();
    vector<tensor<double>> error_DIIS, t1_DIIS, t2_DIIS;
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
                    
                    dot<double>(error_DIIS[ii],"i",error_DIIS[jj],"i",B4DIIS(ii,jj));
                    B4DIIS(jj,ii) = B4DIIS(ii,jj);
                }
                B4DIIS(tmp_size, ii) = -1.0;
                B4DIIS(ii, tmp_size) = -1.0;
                vec_b(ii) = 0.0;
            }
            B4DIIS(tmp_size, tmp_size) = 0.0;
            vec_b(tmp_size) = -1.0;
            VectorXd C = B4DIIS.partialPivLu().solve(vec_b);
            t1 = 0.0;
            t2 = 0.0;
            for(int ii = 0; ii < tmp_size; ii++)
            {
                add<double>(C(ii),t1_DIIS[ii],"ij",1.0,t1,"ij");
                add<double>(C(ii),t2_DIIS[ii],"ijab",1.0,t2,"ijab");
            }
            evaluate_tau();
            evaluate_W_F();
            evaluate_t1t2New();
        }
        d_t1 = evaluateChange_t1(t1,t1_new);
        d_t2 = evaluateChange_t2(t2,t2_new);
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