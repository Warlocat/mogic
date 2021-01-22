#include<iostream>
#include<fstream>
#include<omp.h>
#include<iomanip>
#include"mol.h"
#include"molint.h"
#include"scf.h"
#include"ccsd.h"
using namespace std;
using namespace Eigen;

VectorXd integralTransfermation(const VectorXd& h2e_ao, const MatrixXd& coeff, const string alg = "smart");
VectorXd integralTransfermation_spatial2spin(const VectorXd& h2e_mo, const int& size_basis);
double get_energy_MP2(const VectorXd& h2e_mo, const VectorXd& ene_mo, const int& nelec_a, const int& nelec_b);

int main()
{
    MOL mol_test;
    mol_test.readxyz("h2o.xyz");
    mol_test.readInput("inputFile");
    // cout << mol_test.atomList << endl;
    // cout << mol_test.basisSet << endl;
    // for(int ii = 0; ii < mol_test.coordCart.rows(); ii++)
    //     cout << mol_test.coordCart(ii).transpose() << endl;

    //MOLINT molint_test(mol_test,false);
    MOLINT molint_test(mol_test,true);

    // MatrixXd s = molint_test.get_h1e("nucV");
    // for(int ii = 0; ii < s.rows(); ii++)
    // for(int jj = 0; jj < ii+1 ; jj++)
    // {
    //    cout << s(ii,jj) << endl;
    // }
    // MatrixXd eri = molint_test.get_h2e("eriLLLL");
    // for(int ii = 0; ii < eri.rows(); ii++)
    //     cout << setprecision(15) << eri(ii) << endl;
// cout << molint_test.get_V_RR() << endl;
    MatrixXd s = molint_test.get_h1e("overlap");
    MatrixXd t = molint_test.get_h1e("kinetic");
    MatrixXd v = molint_test.get_h1e("nucV");
    VectorXd eri = molint_test.get_h2e("eriLLLL");
    double e_nuc = molint_test.get_V_RR();
    int size = s.rows();
    
    ofstream ofs;
    ofs.open("s.dat");
        for(int ii = 0; ii < size; ii++)
        for(int jj = 0; jj <= ii; jj++)
            ofs << ii+1 << "\t" << jj+1 << "\t" << setprecision(16) << s(ii,jj) << endl;
    ofs.close();
    ofs.open("t.dat");
        for(int ii = 0; ii < size; ii++)
        for(int jj = 0; jj <= ii; jj++)
            ofs << ii+1 << "\t" << jj+1 << "\t" << setprecision(16) << t(ii,jj) << endl;
    ofs.close();
    ofs.open("v.dat");
        for(int ii = 0; ii < size; ii++)
        for(int jj = 0; jj <= ii; jj++)
            ofs << ii+1 << "\t" << jj+1 << "\t" << setprecision(16) << v(ii,jj) << endl;
    ofs.close();
    ofs.open("enuc.dat");
        ofs << setprecision(16) << e_nuc << endl;
    ofs.close();
    ofs.open("eri.dat");
        for(int ii = 0; ii < size; ii++)
        for(int jj = 0; jj <= ii; jj++)
        for(int kk = 0; kk < size; kk++)
        for(int ll = 0; ll <= kk; ll++)
        {
            int ij = ii*(ii+1)/2+jj, kl = kk*(kk+1)/2+ll;
            if(ij < kl) continue;
            else
            {
                int ijkl = ij*(ij+1)/2+kl;
                ofs << ii+1 << "\t" << jj+1 << "\t" << kk+1 << "\t" << ll+1 << "\t" << eri(ijkl) << endl;
            }
            
        }
    ofs.close();


    // RHF rhf_test(mol_test, molint_test.get_h1e("overlap"), molint_test.get_h1e("kinetic"), molint_test.get_h1e("nucV"), molint_test.get_h2e("eriLLLL"), molint_test.get_V_RR());
    // rhf_test.runSCF();


    RHF rhf_test(5,5,s.rows());
    rhf_test.runSCF();


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

VectorXd integralTransfermation_spatial2spin(const VectorXd& h2e_mo, const int& size_basis)
{
    int tmp_i = 2*size_basis*(2*size_basis+1)/2 + 2*size_basis;
    int size = tmp_i*(tmp_i+1)/2+tmp_i;
    VectorXd h2e_mo_so;
    h2e_mo_so.resize(size);
    h2e_mo_so = VectorXd::Zero(size);
    for(int ii = 0; ii < 2*size_basis; ii++)
    for(int jj = 0; jj <= ii; jj++)
    for(int kk = 0; kk < 2*size_basis; kk++)
    for(int ll = 0; ll <= kk; ll++)
    {
        int ij = ii*(ii+1)/2+jj, kl = kk*(kk+1)/2+ll;
        if((ij < kl) || (ii%2 != jj%2) || (kk%2 != ll%2)) continue;
        else
        {
            int spa_ii = ii / 2, spa_jj = jj / 2, spa_kk = kk / 2, spa_ll = ll / 2;
            int spa_ij = spa_ii*(spa_ii+1)/2+spa_jj, spa_kl = spa_kk*(spa_kk+1)/2+spa_ll;
            int ijkl = ij*(ij+1)/2+kl, spa_ijkl = max(spa_ij,spa_kl)*(max(spa_ij,spa_kl)+1)/2+min(spa_kl,spa_ij);
            h2e_mo_so(ijkl) = h2e_mo(spa_ijkl);
        }
    }

    return h2e_mo_so;
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
