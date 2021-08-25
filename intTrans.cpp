#include<iostream>
#include<fstream>
#include<omp.h>
#include"intTrans.h"
using namespace std;
using namespace Eigen;

/*
    For RHF
*/

/* Integral transformation from AO to restricted MO (pq|rs) */
VectorXd integralTransfermation_MO(const VectorXd& h2e_ao, const MatrixXd& coeff, const string alg)
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
        h2e_tmp1.resize(0,0);

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
        h2e_tmp2.resize(0,0);

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
        h2e_tmp3.resize(0,0);
    }
    else
    {
        cout << "ERROR: intTrans can be only use with alg = 'Noddy' or 'smart'." << endl;
        exit(99);
    }
    
    return h2e_mo;
}


/* Integral transformation from AO to spin orbitals (pq|rs) */
VectorXd integralTransfermation_SO(const VectorXd& h2e_ao, const MatrixXd& coeff, const string alg)
{
    int size_basis = coeff.rows();
    VectorXd h2e_mo = integralTransfermation_MO(h2e_ao, coeff, alg);
    int tmp_i = 2*size_basis*(2*size_basis+1)/2 + 2*size_basis;
    int size = tmp_i*(tmp_i+1)/2+tmp_i;
    VectorXd h2e_so;
    h2e_so.resize(size);
    h2e_so = VectorXd::Zero(size);
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
            h2e_so(ijkl) = h2e_mo(spa_ijkl);
        }
    }

    return h2e_so;
}

/* Integral transformation from AO to anti-symmetrized spin orbitals <pq||rs> */
VectorXd integralTransfermation_SO_antiSym(const VectorXd& h2e_ao, const MatrixXd& coeff, const string alg)
{
    int size_basis = coeff.rows();
    VectorXd h2e_so = integralTransfermation_SO(h2e_ao, coeff, alg);
    VectorXd h2e_so_antiSym(h2e_so.rows());
    h2e_so_antiSym = VectorXd::Zero(h2e_so.rows());
    for(int ii = 0; ii < size_basis*2; ii++)
    for(int jj = 0; jj <= ii; jj++)
    for(int kk = 0; kk < size_basis*2; kk++)
    for(int ll = 0; ll <= kk; ll++)
    {
        int ij = ii*(ii+1)/2+jj, kl = kk*(kk+1)/2+ll;

        if(ij < kl) continue;
        else
        {
            int ijkl = ij*(ij+1)/2 + kl;
            int ik = max(ii,kk)*(max(ii,kk)+1)/2+min(ii,kk), il = max(ii,ll)*(max(ii,ll)+1)/2+min(ii,ll);
            int jk = max(jj,kk)*(max(jj,kk)+1)/2+min(jj,kk), jl = max(jj,ll)*(max(jj,ll)+1)/2+min(jj,ll);
            int ikjl = max(ik,jl)*(max(ik,jl)+1)/2+min(ik,jl), iljk=max(il,jk)*(max(il,jk)+1)/2+min(il,jk);
            h2e_so_antiSym(ijkl) = h2e_so(ikjl) - h2e_so(iljk);
        }
    }

    return h2e_so_antiSym;
}


/* Evaluate MP2 energy */
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



/*
    For UHF
*/

/* Integral transformation from AO to restricted MO (pq|rs) */
Matrix<VectorXd,4,1> integralTransfermation_MO(const VectorXd& h2e_ao, const MatrixXd& coeff_a, const MatrixXd& coeff_b, const string alg)
{
    /*
        (aa|bb) and (bb|aa) CAN be stored in the same array.
        However, my data structure assumes ij>=kl so that I have to treat them seperately.
    */
    Matrix<VectorXd,4,1> h2e_mo;
    h2e_mo(0).resize(h2e_ao.rows()); //aaaa
    h2e_mo(1).resize(h2e_ao.rows()); //bbbb
    h2e_mo(2).resize(h2e_ao.rows()); //aabb
    h2e_mo(3).resize(h2e_ao.rows()); //bbaa
    h2e_mo(0) = VectorXd::Zero(h2e_ao.rows());
    h2e_mo(1) = VectorXd::Zero(h2e_ao.rows());
    h2e_mo(2) = VectorXd::Zero(h2e_ao.rows());
    h2e_mo(3) = VectorXd::Zero(h2e_ao.rows());

    if(alg == "Noddy")
    {
        #pragma omp parallel  for
        for(int ii = 0; ii < coeff_a.rows(); ii++)
        for(int jj = 0; jj <= ii; jj++)
        for(int kk = 0; kk < coeff_a.rows(); kk++)
        for(int ll = 0; ll <= kk; ll++)
        {
            int ij = ii*(ii+1)/2+jj, kl = kk*(kk+1)/2+ll;
            if(ij < kl) continue;
            else
            {
                int index = ij*(ij+1)/2+kl;
                h2e_mo(0)(index) = 0.0;
                h2e_mo(1)(index) = 0.0;
                h2e_mo(2)(index) = 0.0;
                h2e_mo(3)(index) = 0.0;
                for(int aa = 0; aa < coeff_a.rows(); aa++)
                for(int bb = 0; bb < coeff_a.rows(); bb++)
                for(int cc = 0; cc < coeff_a.rows(); cc++)
                for(int dd = 0; dd < coeff_a.rows(); dd++)
                {
                    int ab = max(aa,bb)*(max(aa,bb)+1)/2+min(aa,bb), cd = max(cc,dd)*(max(cc,dd)+1)/2+min(cc,dd);
                    int abcd = max(ab,cd)*(max(ab,cd)+1)/2+min(ab,cd);
                    h2e_mo(0)(index) += coeff_a(aa,ii)*coeff_a(bb,jj)*coeff_a(cc,kk)*coeff_a(dd,ll)*h2e_ao(abcd);
                    h2e_mo(1)(index) += coeff_b(aa,ii)*coeff_b(bb,jj)*coeff_b(cc,kk)*coeff_b(dd,ll)*h2e_ao(abcd);
                    h2e_mo(2)(index) += coeff_a(aa,ii)*coeff_a(bb,jj)*coeff_b(cc,kk)*coeff_b(dd,ll)*h2e_ao(abcd);
                    h2e_mo(3)(index) += coeff_b(aa,ii)*coeff_b(bb,jj)*coeff_a(cc,kk)*coeff_a(dd,ll)*h2e_ao(abcd);
                }
            }
        }
    }
    else if(alg == "smart")
    {
        int size = coeff_a.rows();
        int size2 = coeff_a.rows()*coeff_a.rows(), size2p = coeff_a.rows()*(coeff_a.rows()+1)/2;
        MatrixXd h2e_tmp1_aa, h2e_tmp2_aa, h2e_tmp3_aa;
        MatrixXd h2e_tmp1_bb, h2e_tmp2_bb, h2e_tmp3_bb;
        MatrixXd h2e_tmp1_ab, h2e_tmp2_ab, h2e_tmp3_ab;
        MatrixXd h2e_tmp1_ba, h2e_tmp2_ba, h2e_tmp3_ba;
        h2e_tmp1_aa.resize(size2p, size2);   // abcl
        h2e_tmp2_aa.resize(size2p, size2p);   // abkl
        h2e_tmp3_aa.resize(size2, size2p);   // ajkl
        h2e_tmp1_bb.resize(size2p, size2);   // abcl
        h2e_tmp2_bb.resize(size2p, size2p);   // abkl
        h2e_tmp3_bb.resize(size2, size2p);   // ajkl
        h2e_tmp1_ab.resize(size2p, size2);   // abcl
        h2e_tmp2_ab.resize(size2p, size2p);   // abkl
        h2e_tmp3_ab.resize(size2, size2p);   // ajkl
        h2e_tmp1_ba.resize(size2p, size2);   // abcl
        h2e_tmp2_ba.resize(size2p, size2p);   // abkl
        h2e_tmp3_ba.resize(size2, size2p);   // ajkl
        #pragma omp parallel  for
        for(int aa = 0; aa < coeff_a.rows(); aa++)
        for(int bb = 0; bb <= aa; bb++)
        for(int cc = 0; cc < coeff_a.rows(); cc++)
        for(int ll = 0; ll < coeff_a.rows(); ll++)
        {
            int ab = aa*(aa+1)/2+bb, cl = cc*size+ll;
            h2e_tmp1_aa(ab,cl) = 0.0;
            h2e_tmp1_bb(ab,cl) = 0.0;
            h2e_tmp1_ab(ab,cl) = 0.0;
            h2e_tmp1_ba(ab,cl) = 0.0;
            for(int dd = 0; dd < coeff_a.rows(); dd++)
            {
                int cd = max(cc,dd)*(max(cc,dd)+1)/2 + min(cc,dd);
                int abcd = max(ab,cd)*(max(ab,cd)+1)/2 + min(ab,cd);
                h2e_tmp1_aa(ab,cl) += coeff_a(dd,ll) * h2e_ao(abcd);
                h2e_tmp1_bb(ab,cl) += coeff_b(dd,ll) * h2e_ao(abcd);
                h2e_tmp1_ab(ab,cl) += coeff_b(dd,ll) * h2e_ao(abcd);
                h2e_tmp1_ba(ab,cl) += coeff_a(dd,ll) * h2e_ao(abcd);
            }
        }

        for(int aa = 0; aa < coeff_a.rows(); aa++)
        for(int bb = 0; bb <= aa; bb++)
        for(int kk = 0; kk < coeff_a.rows(); kk++)
        for(int ll = 0; ll <= kk; ll++)
        {
            int ab = aa*(aa+1)/2+bb, kl = kk*(kk+1)/2+ll;
            h2e_tmp2_aa(ab,kl) = 0.0;
            h2e_tmp2_bb(ab,kl) = 0.0;
            h2e_tmp2_ab(ab,kl) = 0.0;
            h2e_tmp2_ba(ab,kl) = 0.0;
            for(int cc = 0; cc < coeff_a.rows(); cc++)
            {
                int cl = cc*size+ll;
                h2e_tmp2_aa(ab,kl) += coeff_a(cc,kk) * h2e_tmp1_aa(ab,cl);
                h2e_tmp2_bb(ab,kl) += coeff_b(cc,kk) * h2e_tmp1_bb(ab,cl);
                h2e_tmp2_ab(ab,kl) += coeff_b(cc,kk) * h2e_tmp1_ab(ab,cl);
                h2e_tmp2_ba(ab,kl) += coeff_a(cc,kk) * h2e_tmp1_ba(ab,cl);
            }
        }
        h2e_tmp1_aa.resize(0,0);
        h2e_tmp1_bb.resize(0,0);
        h2e_tmp1_ab.resize(0,0);
        h2e_tmp1_ba.resize(0,0);

        for(int aa = 0; aa < coeff_a.rows(); aa++)
        for(int jj = 0; jj < coeff_a.rows(); jj++)
        for(int kk = 0; kk < coeff_a.rows(); kk++)
        for(int ll = 0; ll <= kk; ll++)
        {
            int aj = aa*size+jj, kl = kk*(kk+1)/2+ll;
            h2e_tmp3_aa(aj,kl) = 0.0;
            h2e_tmp3_bb(aj,kl) = 0.0;
            h2e_tmp3_ab(aj,kl) = 0.0;
            h2e_tmp3_ba(aj,kl) = 0.0;
            for(int bb = 0; bb < coeff_a.rows(); bb++)
            {
                int ab = max(aa,bb)*(max(aa,bb)+1)/2 + min(aa,bb);
                h2e_tmp3_aa(aj,kl) += coeff_a(bb,jj) * h2e_tmp2_aa(ab,kl);
                h2e_tmp3_bb(aj,kl) += coeff_b(bb,jj) * h2e_tmp2_bb(ab,kl);
                h2e_tmp3_ab(aj,kl) += coeff_a(bb,jj) * h2e_tmp2_ab(ab,kl);
                h2e_tmp3_ba(aj,kl) += coeff_b(bb,jj) * h2e_tmp2_ba(ab,kl);
            }
        }
        h2e_tmp2_aa.resize(0,0);
        h2e_tmp2_bb.resize(0,0);
        h2e_tmp2_ab.resize(0,0);
        h2e_tmp2_ba.resize(0,0);

        for(int ii = 0; ii < coeff_a.rows(); ii++)
        for(int jj = 0; jj <= ii; jj++)
        for(int kk = 0; kk < coeff_a.rows(); kk++)
        for(int ll = 0; ll <= kk; ll++)
        {
            int ij = ii*(ii+1)/2+jj, kl = kk*(kk+1)/2+ll;
            if(ij < kl) continue;
            else
            {
                int ijkl = ij*(ij+1)/2+kl;
                h2e_mo(0)(ijkl) = 0.0;
                h2e_mo(1)(ijkl) = 0.0;
                h2e_mo(2)(ijkl) = 0.0;
                h2e_mo(3)(ijkl) = 0.0;
                for(int aa = 0; aa < coeff_a.rows(); aa++)
                {
                    int aj = aa*size + jj;
                    h2e_mo(0)(ijkl) += coeff_a(aa,ii) * h2e_tmp3_aa(aj,kl);
                    h2e_mo(1)(ijkl) += coeff_b(aa,ii) * h2e_tmp3_bb(aj,kl);
                    h2e_mo(2)(ijkl) += coeff_a(aa,ii) * h2e_tmp3_ab(aj,kl);
                    h2e_mo(3)(ijkl) += coeff_b(aa,ii) * h2e_tmp3_ba(aj,kl);
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

/* Integral transformation from AO to spin orbitals (pq|rs) */
VectorXd integralTransfermation_SO(const VectorXd& h2e_ao, const MatrixXd& coeff_a, const MatrixXd& coeff_b, const string alg)
{
    int size_basis = coeff_a.rows();
    Matrix<VectorXd,4,1> h2e_mo = integralTransfermation_MO(h2e_ao, coeff_a, coeff_b, alg);
    int tmp_i = 2*size_basis*(2*size_basis+1)/2 + 2*size_basis;
    int size = tmp_i*(tmp_i+1)/2+tmp_i;
    VectorXd h2e_so;
    h2e_so.resize(size);
    h2e_so = VectorXd::Zero(size);
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
            if(ii%2==0 && kk%2==0)
                h2e_so(ijkl) = h2e_mo(0)(spa_ijkl);
            else if(ii%2!=0 && kk%2!=0)
                h2e_so(ijkl) = h2e_mo(1)(spa_ijkl);
            else if(ii%2==0 && kk%2!=0)
                h2e_so(ijkl) = h2e_mo(2)(spa_ijkl);
            else
                h2e_so(ijkl) = h2e_mo(3)(spa_ijkl);
        }
    }

    return h2e_so;
}

/* Integral transformation from AO to anti-symmetrized spin orbitals <pq||rs> */
VectorXd integralTransfermation_SO_antiSym(const VectorXd& h2e_ao, const MatrixXd& coeff_a, const MatrixXd& coeff_b, const string alg)
{
    int size_basis = coeff_a.rows();
    VectorXd h2e_so = integralTransfermation_SO(h2e_ao, coeff_a, coeff_b, alg);
    VectorXd h2e_so_antiSym(h2e_so.rows());
    h2e_so_antiSym = VectorXd::Zero(h2e_so.rows());
    for(int ii = 0; ii < size_basis*2; ii++)
    for(int jj = 0; jj <= ii; jj++)
    for(int kk = 0; kk < size_basis*2; kk++)
    for(int ll = 0; ll <= kk; ll++)
    {
        int ij = ii*(ii+1)/2+jj, kl = kk*(kk+1)/2+ll;

        if(ij < kl) continue;
        else
        {
            int ijkl = ij*(ij+1)/2 + kl;
            int ik = max(ii,kk)*(max(ii,kk)+1)/2+min(ii,kk), il = max(ii,ll)*(max(ii,ll)+1)/2+min(ii,ll);
            int jk = max(jj,kk)*(max(jj,kk)+1)/2+min(jj,kk), jl = max(jj,ll)*(max(jj,ll)+1)/2+min(jj,ll);
            int ikjl = max(ik,jl)*(max(ik,jl)+1)/2+min(ik,jl), iljk=max(il,jk)*(max(il,jk)+1)/2+min(il,jk);
            h2e_so_antiSym(ijkl) = h2e_so(ikjl) - h2e_so(iljk);
        }
    }

    return h2e_so_antiSym;
}








