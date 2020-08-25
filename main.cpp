#include<iostream>
#include<fstream>
#include<omp.h>
#include"mol.h"
#include"scf.h"
using namespace std;
using namespace Eigen;

VectorXd integralTransfermation(const VectorXd& h2e_ao, const MatrixXd& coeff, const string alg = "Noddy");
double get_energy_MP2(const VectorXd& h2e_mo, const VectorXd& ene_mo, const int& nelec_a, const int& nelec_b);

int main()
{
    int nelec_a = 5, nelec_b = 5, size_basis = 14;
    RHF rhf_test(nelec_a, nelec_b, size_basis);
    rhf_test.runSCF();
    VectorXd h2e_mo = integralTransfermation(rhf_test.get_h2e_vector(), rhf_test.coeff, "smart");
    double ene_mp2 = get_energy_MP2(h2e_mo, rhf_test.ene_orb, nelec_a, nelec_b);
    cout << ene_mp2 << "\t" << rhf_test.ene_scf + ene_mp2 << endl;

    cout << "Congratulation! This program finished successfully!"  << endl;
    return 0;
}


VectorXd integralTransfermation(const VectorXd& h2e_ao, const MatrixXd& coeff, const string alg)
{
    VectorXd h2e_mo;
    h2e_mo.resize(h2e_ao.rows());
    h2e_mo = VectorXd::Zero(h2e_ao.rows());

    if(alg == "Noddy")
    {
        #pragma omp parallel  for
        for(int ii = 0; ii < coeff.rows(); ii++)
        for(int jj = 0; jj <= ii; jj++)
        for(int kk = 0; kk < coeff.rows(); kk++)
        for(int ll = 0; ll <= kk; ll++)
        {
            int ij = ii*(ii+1)/2+jj, kl = kk*(kk+1)/2+ll;
            if(ij < kl) continue;
            else
            {
                int index = ij*(ij+1)/2+kl;
                h2e_mo(index) = 0.0;
                for(int aa = 0; aa < coeff.rows(); aa++)
                for(int bb = 0; bb < coeff.rows(); bb++)
                for(int cc = 0; cc < coeff.rows(); cc++)
                for(int dd = 0; dd < coeff.rows(); dd++)
                {
                    int ab = max(aa,bb)*(max(aa,bb)+1)/2+min(aa,bb), cd = max(cc,dd)*(max(cc,dd)+1)/2+min(cc,dd);
                    int abcd = max(ab,cd)*(max(ab,cd)+1)/2+min(ab,cd);
                    // h2e_mo += conj(coeff(aa,ii))*coeff(bb,jj)*conj(coeff(cc,kk))*coeff(dd,ll)*h2e_ao(abcd);
                    h2e_mo(index) += coeff(aa,ii)*coeff(bb,jj)*coeff(cc,kk)*coeff(dd,ll)*h2e_ao(abcd);
                }
            }
        }
    }
    else if(alg == "smart")
    {
        int size = coeff.rows();
        int size2 = coeff.rows()*coeff.rows(), size2p = coeff.rows()*(coeff.rows()+1)/2;
        MatrixXd h2e_tmp1, h2e_tmp2, h2e_tmp3;
        h2e_tmp1.resize(size2p, size2);   // abcl
        h2e_tmp2.resize(size2p, size2p);   // abkl
        h2e_tmp3.resize(size2, size2p);   // ajkl
        #pragma omp parallel  for
        for(int aa = 0; aa < coeff.rows(); aa++)
        for(int bb = 0; bb <= aa; bb++)
        for(int cc = 0; cc < coeff.rows(); cc++)
        for(int ll = 0; ll < coeff.rows(); ll++)
        {
            int ab = aa*(aa+1)/2+bb, cl = cc*size+ll;
            h2e_tmp1(ab,cl) = 0.0;
            for(int dd = 0; dd < coeff.rows(); dd++)
            {
                int cd = max(cc,dd)*(max(cc,dd)+1)/2 + min(cc,dd);
                int abcd = max(ab,cd)*(max(ab,cd)+1)/2 + min(ab,cd);
                h2e_tmp1(ab,cl) += coeff(dd,ll) * h2e_ao(abcd);
            }
        }

        for(int aa = 0; aa < coeff.rows(); aa++)
        for(int bb = 0; bb <= aa; bb++)
        for(int kk = 0; kk < coeff.rows(); kk++)
        for(int ll = 0; ll <= kk; ll++)
        {
            int ab = aa*(aa+1)/2+bb, kl = kk*(kk+1)/2+ll;
            h2e_tmp2(ab,kl) = 0.0;
            for(int cc = 0; cc < coeff.rows(); cc++)
            {
                int cl = cc*size+ll;
                h2e_tmp2(ab,kl) += coeff(cc,kk) * h2e_tmp1(ab,cl);
            }
        }

        for(int aa = 0; aa < coeff.rows(); aa++)
        for(int jj = 0; jj < coeff.rows(); jj++)
        for(int kk = 0; kk < coeff.rows(); kk++)
        for(int ll = 0; ll <= kk; ll++)
        {
            int aj = aa*size+jj, kl = kk*(kk+1)/2+ll;
            h2e_tmp3(aj,kl) = 0.0;
            for(int bb = 0; bb < coeff.rows(); bb++)
            {
                int ab = max(aa,bb)*(max(aa,bb)+1)/2 + min(aa,bb);
                h2e_tmp3(aj,kl) += coeff(bb,jj) * h2e_tmp2(ab,kl);
            }
        }

        for(int ii = 0; ii < coeff.rows(); ii++)
        for(int jj = 0; jj <= ii; jj++)
        for(int kk = 0; kk < coeff.rows(); kk++)
        for(int ll = 0; ll <= kk; ll++)
        {
            int ij = ii*(ii+1)/2+jj, kl = kk*(kk+1)/2+ll;
            if(ij < kl) continue;
            else
            {
                int ijkl = ij*(ij+1)/2+kl;
                h2e_mo(ijkl) = 0.0;
                for(int aa = 0; aa < coeff.rows(); aa++)
                {
                    int aj = aa*size + jj;
                    h2e_mo(ijkl) += coeff(aa,ii) * h2e_tmp3(aj,kl);
                }
            }
        }
    }
    else
    {
        cout << "ERROR: intTrans can be only use with alg = 'Noddy' or 'smart'." << endl;
        exit(99);
    }
    
    return h2e_mo;
}


double get_energy_MP2(const VectorXd& h2e_mo, const VectorXd& ene_mo, const int& nelec_a, const int& nelec_b)
{
    double ene_MP2 = 0.0;
    if(nelec_a == nelec_b)
    {
        for(int ii = 0; ii < nelec_a; ii++)
        for(int jj = 0; jj < nelec_a; jj++)
        for(int aa = nelec_a; aa < ene_mo.rows(); aa++)
        for(int bb = nelec_a; bb < ene_mo.rows(); bb++)
        {
            int ia = aa*(aa+1)/2+ii, jb = bb*(bb+1)/2+jj, ib = bb*(bb+1)/2+ii, ja = aa*(aa+1)/2+jj;
            int iajb = max(ia,jb)*(max(ia,jb)+1)/2+min(ia,jb), ibja = max(ib,ja)*(max(ib,ja)+1)/2+min(ib,ja);
            ene_MP2 += h2e_mo(iajb)*(2.0*h2e_mo(iajb) - h2e_mo(ibja)) / (ene_mo(ii) + ene_mo(jj) - ene_mo(aa) - ene_mo(bb));
        }
    }
    else
    {
        cout << "ERROR: MP2 for open-shell system has mot been implemented yet." << endl;
        exit(99);
    }
    
    return ene_MP2;
}