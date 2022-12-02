#include"cis.h"
#include<omp.h>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>


/*
    Constructor for RHF
*/
CIS::CIS(const int& n_occ_, const int& n_vir_, const int& nroot_, const VectorXd& h2e_so_antiSym_, const VectorXd& ene_mo_):
n_occ(n_occ_), n_vir(n_vir_), nroot(nroot_), h2e_so_antiSym(h2e_so_antiSym_)
{
    ene_so.resize(n_occ+n_vir);
    for(int ii = 0; ii < n_occ+n_vir; ii++)
        ene_so(ii) = ene_mo_(ii/2);
}

/*
    Constructor for UHF
*/
CIS::CIS(const int& n_occ_, const int& n_vir_, const int& nroot_, const VectorXd& h2e_so_antiSym_, const VectorXd& ene_mo_a, const VectorXd& ene_mo_b):
n_occ(n_occ_), n_vir(n_vir_), nroot(nroot_), h2e_so_antiSym(h2e_so_antiSym_)
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
VectorXd CIS::HdotC_CIS(const VectorXd& Cin)
{
    VectorXd Cout(n_occ*n_vir);
    for(int ii = 0; ii < n_occ; ii++)
    for(int aa = 0; aa < n_vir; aa++)
    {
        int ia = ii*n_vir+aa;
        Cout(ia) = 0.0;
        for(int jj = 0; jj < n_occ; jj++)
        for(int bb = 0; bb < n_vir; bb++)
        {
            int jb = jj*n_vir+bb;
            Cout(ia) += evaluateH_CIS(ii, jj, aa+n_occ, bb+n_occ)*Cin(jb);
        }
    }
    return Cout;
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

void CIS::runCIS()
{
    MatrixXd H_CIS;
    if(!davidsonLiu)
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
                H_CIS(ia, jb) = evaluateH_CIS(ii, jj, aa+n_occ, bb+n_occ);
                H_CIS(jb, ia) = H_CIS(ia, jb);
            }
        }
        SelfAdjointEigenSolver<MatrixXd> solver(H_CIS); 
        ene_CIS = solver.eigenvalues();
        CIScoeff = solver.eigenvectors();
    }
    else
    {
        VectorXd Hdiag(n_occ*n_vir);
        for(int ii = 0; ii < n_occ; ii++)
        for(int aa = 0; aa < n_vir; aa++)
            Hdiag(ii*n_vir+aa) = evaluateH_CIS(ii,ii,aa+n_occ,aa+n_occ);
        DLsolver_CIS(Hdiag, CIScoeff, ene_CIS, nroot, maxCycle);
    }
}

void CIS::runTDHF()
{
    MatrixXd A_TDHF, B_TDHF;
    if(!davidsonLiu)
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
        SelfAdjointEigenSolver<MatrixXd> solver((A_TDHF + B_TDHF) * (A_TDHF - B_TDHF));
        ene_CIS = solver.eigenvalues();
        for(int ii = 0; ii < ene_CIS.rows(); ii++)  ene_CIS(ii) = sqrt(ene_CIS(ii));
        CIScoeff = solver.eigenvectors();
    }
    else
    {
        cout << "Davidson-Liu algorithm has NOT been implemented for TDHF." << endl;
        exit(99);
    }
}

void CIS::DLsolver_CIS(const VectorXd& Hdiag, MatrixXd& eigenvectors, VectorXd& eigenvalues, const int& nroot, const int& maxCycle)
{
    cout << endl << "Start Davidson iterations for CIS." << endl;
    eigenvalues.resize(nroot);
    eigenvectors.resize(n_occ*n_vir,nroot);
    int DLsize = 1;
    vector<VectorXd> DLvectors;
    vector<VectorXd> HC;
    DLvectors.push_back(VectorXd::Zero(n_vir*n_occ));
    DLvectors[0]((n_occ-1)*n_vir) = 1.0;
    VectorXd HC_new = HdotC_CIS(DLvectors[0]);
    HC.push_back(HC_new);
    double norm12 = DLvectors[0].adjoint()*HC_new;
    // VectorXd HC_new_ortho = HC_new - norm12*DLvectors[0];
    // HC_new_ortho = HC_new_ortho/HC_new_ortho.norm();
    MatrixXd H_sub_old(1,1), H_sub_new;
    H_sub_old(0,0) = norm12;
    for(int ic = 1; ic < maxCycle; ic++)
    {
        cout << "Iter #" << ic << "\t\tResidual: ";
        SelfAdjointEigenSolver<MatrixXd> solver(H_sub_old);
        MatrixXd vectors = solver.eigenvectors();
        VectorXd values = solver.eigenvalues();

        int pick = 0;
        double ene = values(0);
        VectorXd target_sub = vectors.block(0,0,vectors.rows(),1);
        VectorXd target;
        for(int ii = 1; ii < H_sub_old.rows(); ii++)
        {
            if(values(ii) < ene)
            {
                ene = values(ii);
                target_sub = vectors.block(0,pick,vectors.rows(),1);
                pick = ii;
            }
        }

        HC_new = VectorXd::Zero(HC_new.rows());
        target = VectorXd::Zero(HC_new.rows());
        for(int ii = 0; ii < DLsize; ii++)
        {
            HC_new += target_sub(ii)*HC[ii];
            target += target_sub(ii)*DLvectors[ii];
        }
        
        VectorXd R = HC_new - ene*target;
        double Maxresidue = max(abs(R.maxCoeff()), abs(R.minCoeff()));
        cout << Maxresidue <<endl;
        converged = (Maxresidue < convControl);
        if(converged)
        {
            eigenvalues(0) = ene;
            // eigenvectors = target;
            break;
        }

        VectorXd newVec(R.rows());
        for(int ii = 0; ii < R.rows(); ii++)
        {
            if(abs(Hdiag(ii)-ene)<1e-8) newVec(ii) = target(ii);
            else newVec(ii) = target(ii) - R(ii)/(Hdiag(ii)-ene);
        }
        for(int ii = 0; ii < DLsize; ii++)
        {
            newVec = newVec - (DLvectors[ii].adjoint()*newVec)(0,0)*DLvectors[ii];
        }
        newVec = newVec/newVec.norm();
        DLvectors.push_back(newVec);
        HC.push_back(HdotC_CIS(newVec));
        DLsize++;

        H_sub_new = MatrixXd::Zero(DLsize,DLsize);
        for(int ii = 0; ii < DLsize-1; ii++)
        {
            for(int jj = 0; jj < DLsize-1; jj++)
                H_sub_new(ii,jj) = H_sub_old(ii,jj);
            H_sub_new(DLsize-1,ii) = DLvectors[DLsize-1].adjoint() * HC[ii];
            H_sub_new(ii,DLsize-1) = DLvectors[ii].adjoint() * HC[DLsize-1];
        }
        H_sub_new(DLsize-1,DLsize-1) = DLvectors[DLsize-1].adjoint() * HC[DLsize-1];
        H_sub_old = H_sub_new;
    }
}