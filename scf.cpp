#include"scf.h"
#include<omp.h>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>

SCF::SCF(const MOL& mol_, const MatrixXd& s_, const MatrixXd& t_, const MatrixXd& v_, const VectorXd& eri_, const double& V_RR_):
nelec_a(mol_.nelec_a), nelec_b(mol_.nelec_b), size_basis(s_.rows()), overlap(s_), h1e(t_+v_), h2e(eri_), V_RR(V_RR_)
{
    cout << "Special constructor of SCF class using MOL and integrals is called." << endl;
    overlap_half_i = matrix_half_inverse(overlap);
}
SCF::SCF(const int& nelec_a_, const int& nelec_b_, const int& size_basis_):
nelec_a(nelec_a_), nelec_b(nelec_b_), size_basis(size_basis_)
{
    cout << "Special constructor of SCF class used in CrawfordGroupProjects is called." << endl;
    MatrixXd tmpT(size_basis, size_basis);
    overlap.resize(size_basis, size_basis);
    h1e.resize(size_basis, size_basis);
    int tmp_i = size_basis*(size_basis+1)/2;
    h2e.resize(tmp_i*(tmp_i+1)/2);

    ifstream ifs;
    ifs.open("enuc.dat");
        ifs >> V_RR;
    ifs.close();
    
    readIntegrals_1e(overlap, "s.dat");
    overlap_half_i = matrix_half_inverse(overlap);
    
    readIntegrals_1e(tmpT, "t.dat");
    readIntegrals_1e(h1e, "v.dat");
    h1e = h1e + tmpT;

    readIntegrals_2e(h2e, "eri.dat");
}

SCF::~SCF()
{
}


MatrixXd SCF::evaluateErrorDIIS(const MatrixXd& fock_, const MatrixXd& density_)
{
    MatrixXd tmp = fock_*density_*overlap - overlap*density_*fock_;
    MatrixXd err(size_basis*size_basis,1);
    for(int ii = 0; ii < size_basis; ii++)
    for(int jj = 0; jj < size_basis; jj++)
    {
        err(ii*size_basis+jj,0) = tmp(ii,jj);
    }
    return err;
}




void SCF::readIntegrals_1e(MatrixXd& int_1e, const string& filename)
{
    int iii, jjj;
    double tmp_d;
    ifstream ifs;
    int_1e = MatrixXd::Zero(int_1e.rows(), int_1e.cols());
    ifs.open(filename);
        while (!ifs.eof())
        {
            ifs >> iii >> jjj >> tmp_d;
            int_1e(iii - 1,jjj - 1) = tmp_d;
            int_1e(jjj - 1,iii - 1) = tmp_d;
        }
    ifs.close();

    return;
}

void SCF::readIntegrals_2e(VectorXd& int_2e, const string& filename)
{
    int_2e = VectorXd::Zero(int_2e.rows());
    double tmp_d;
    int iii,jjj,kkk,lll;
    ifstream ifs;
    ifs.open(filename);
        while(!ifs.eof())
        {
            ifs >> iii >> jjj >> kkk >> lll >> tmp_d;
            int ij = (iii - 1) * iii / 2 + jjj - 1, kl = (kkk - 1) * kkk / 2 + lll - 1;
            int_2e(ij*(ij+1)/2 + kl) = tmp_d;
        }            
    ifs.close();

    return;
}

double SCF::evaluateChange(const MatrixXd& M1, const MatrixXd& M2)
{
    int size = M1.rows();
    double tmp = 0.0;
    for(int ii = 0; ii < size; ii++)
    for(int jj = 0; jj < size; jj++)
        tmp += pow((M1(ii,jj) - M2(ii,jj)),2);

    return sqrt(tmp);
}


MatrixXd SCF::matrix_half_inverse(const MatrixXd& inputM)
{
    int size = inputM.rows();
    SelfAdjointEigenSolver<MatrixXd> solver(inputM);
    VectorXd eigenvalues = solver.eigenvalues();
    MatrixXd eigenvectors = solver.eigenvectors();
    
    for(int ii = 0; ii < size; ii++)
    {
        if(eigenvalues(ii) < 0)
        {
            cout << "ERROR: Overlap Matrixhas negative eigenvalues!" << endl;
            exit(99);
        }
        else
        {
            eigenvalues(ii) = 1.0 / sqrt(eigenvalues(ii));
        }
    }
    MatrixXd tmp(size, size);
    for(int ii = 0; ii < size; ii++)
    for(int jj = 0; jj < size; jj++)
    {
        tmp(ii,jj) = 0.0;
        for(int kk = 0; kk < size; kk++)
            tmp(ii,jj) += eigenvectors(ii,kk) * eigenvalues(kk) * eigenvectors(jj,kk);
    }

    return tmp; 
}

MatrixXd SCF::matrix_half(const MatrixXd& inputM)
{
    int size = inputM.rows();
    SelfAdjointEigenSolver<MatrixXd> solver(inputM);
    VectorXd eigenvalues = solver.eigenvalues();
    MatrixXd eigenvectors = solver.eigenvectors();
    
    for(int ii = 0; ii < size; ii++)
    {
        if(eigenvalues(ii) < 0)
        {
            cout << "ERROR: Overlap Matrixhas negative eigenvalues!" << endl;
            exit(99);
        }
        else
        {
            eigenvalues(ii) = sqrt(eigenvalues(ii));
        }
    }
    MatrixXd tmp(size, size);
    for(int ii = 0; ii < size; ii++)
    for(int jj = 0; jj < size; jj++)
    {
        tmp(ii,jj) = 0.0;
        for(int kk = 0; kk < size; kk++)
            tmp(ii,jj) += eigenvectors(ii,kk) * eigenvalues(kk) * eigenvectors(jj,kk);
    }

    return tmp;
}

MatrixXd SCF::evaluateDensity(const MatrixXd& coeff_, const int& occ, const bool& spherical)
{
    int size = coeff_.rows();
    MatrixXd den(size, size);
    if(!spherical)
    {
        for(int aa = 0; aa < size; aa++)
        for(int bb = 0; bb < size; bb++)
        {
            den(aa,bb) = 0.0;
            for(int ii = 0; ii < occ; ii++)
                den(aa,bb) += coeff_(aa,ii) * coeff_(bb,ii);
        }
    }
    else
    {
        cout << "ERROR: Spherical occupation is NOT supported now!" << endl;
        exit(99);
    }
    
    
    return den;
}


void SCF::eigensolverG(const MatrixXd& inputM, const MatrixXd& s_h_i, VectorXd& values, MatrixXd& vectors)
{
    MatrixXd tmp = s_h_i * inputM * s_h_i;
    
    
    SelfAdjointEigenSolver<MatrixXd> solver(tmp);
    values = solver.eigenvalues();
    vectors = s_h_i * solver.eigenvectors();

    return;
}

VectorXd SCF::get_h2e_vector()
{
    return h2e;
}



/*********************************************************/
/**********    Member functions of class RHF    **********/
/*********************************************************/

RHF::RHF(const int& nelec_a_, const int& nelec_b_, const int& size_basis_):
SCF(nelec_a_, nelec_b_, size_basis_)
{
    if(nelec_a_ != nelec_b_)
    {
        cout << "ERROR: RHF is called when nelec_a != nelec_b!" << endl;
        exit(99);
    }
}

RHF::RHF(const MOL& mol_, const MatrixXd& s_, const MatrixXd& t_, const MatrixXd& v_, const VectorXd& eri_, const double& V_RR_):
SCF(mol_, s_, t_, v_, eri_, V_RR_)
{
}

RHF::~RHF()
{
}

void RHF::runSCF()
{
    vector<MatrixXd> error4DIIS, fock4DIIS;
    fock.resize(size_basis,size_basis);

    MatrixXd newDen;
    eigensolverG(h1e, overlap_half_i, ene_orb, coeff);
    density = evaluateDensity(coeff, nelec_a);

    for(int iter = 1; iter <= maxIter; iter++)
    {
        if(iter <= 2)
        {
            #pragma omp parallel  for
            for(int mm = 0; mm < size_basis; mm++)
            for(int nn = 0; nn <= mm; nn++)
            {
                fock(mm,nn) = h1e(mm,nn);
                for(int aa = 0; aa < size_basis; aa++)
                for(int bb = 0; bb < size_basis; bb++)
                {
                    int ab, mn, an, mb, abmn, anmb;
                    ab = max(aa,bb)*(max(aa,bb)+1)/2 + min(aa,bb);
                    mn = max(mm,nn)*(max(mm,nn)+1)/2 + min(mm,nn);
                    an = max(aa,nn)*(max(aa,nn)+1)/2 + min(aa,nn);
                    mb = max(mm,bb)*(max(mm,bb)+1)/2 + min(mm,bb);
                    abmn = max(ab,mn)*(max(ab,mn)+1)/2 + min(ab,mn);
                    anmb = max(an,mb)*(max(an,mb)+1)/2 + min(an,mb);
                    fock(mm,nn) += density(aa,bb) * (2.0 * h2e(abmn) - h2e(anmb));
                }
                fock(nn,mm) = fock(mm,nn);
            }
        }
        else
        {
            int tmp_size = fock4DIIS.size();
            MatrixXd B4DIIS(tmp_size+1,tmp_size+1);
            VectorXd vec_b(tmp_size+1);
            for(int ii = 0; ii < tmp_size; ii++)
            {    
                for(int jj = 0; jj <= ii; jj++)
                {
                    B4DIIS(ii,jj) = (error4DIIS[ii].adjoint()*error4DIIS[jj])(0,0);
                    B4DIIS(jj,ii) = B4DIIS(ii,jj);
                }
                B4DIIS(tmp_size, ii) = -1.0;
                B4DIIS(ii, tmp_size) = -1.0;
                vec_b(ii) = 0.0;
            }
            B4DIIS(tmp_size, tmp_size) = 0.0;
            vec_b(tmp_size) = -1.0;
            VectorXd C = B4DIIS.partialPivLu().solve(vec_b);
            fock = MatrixXd::Zero(size_basis,size_basis);
            for(int ii = 0; ii < tmp_size; ii++)
            {
                fock += C(ii) * fock4DIIS[ii];
            }
        }
        eigensolverG(fock, overlap_half_i, ene_orb, coeff);
        // newDen = 0.5*(evaluateDensity(coeff, nelec_a) + density);
        newDen = evaluateDensity(coeff, nelec_a);
        d_density = evaluateChange(density, newDen);
        cout << "Iter #" << iter << " maximum density difference: " << d_density << endl;
        density = newDen;
        if(d_density < convControl) 
        {
            converged = true;
            cout << endl << "RHF converges after " << iter << " iterations." << endl << endl;

            cout << "\tOrbital\t\tEnergy(in hartree)\n";
            cout << "\t*******\t\t******************\n";
            for(int ii = 1; ii <= size_basis; ii++)
            {
                cout << "\t" << ii << "\t\t" << setprecision(15) << ene_orb(ii - 1) << endl;
            }

            ene_scf = V_RR;
            for(int ii = 0; ii < size_basis; ii++)
            for(int jj = 0; jj < size_basis; jj++)
            {
                ene_scf += density(ii,jj) * (h1e(jj,ii) + fock(jj,ii));
            }
            cout << "Final RHF energy is " << setprecision(15) << ene_scf << " hartree." << endl;
            break;            
        }   
        #pragma omp parallel  for
        for(int mm = 0; mm < size_basis; mm++)
        for(int nn = 0; nn <= mm; nn++)
        {
            fock(mm,nn) = h1e(mm,nn);
            for(int aa = 0; aa < size_basis; aa++)
            for(int bb = 0; bb < size_basis; bb++)
            {
                int ab, mn, an, mb, abmn, anmb;
                ab = max(aa,bb)*(max(aa,bb)+1)/2 + min(aa,bb);
                mn = max(mm,nn)*(max(mm,nn)+1)/2 + min(mm,nn);
                an = max(aa,nn)*(max(aa,nn)+1)/2 + min(aa,nn);
                mb = max(mm,bb)*(max(mm,bb)+1)/2 + min(mm,bb);
                abmn = max(ab,mn)*(max(ab,mn)+1)/2 + min(ab,mn);
                anmb = max(an,mb)*(max(an,mb)+1)/2 + min(an,mb);
                fock(mm,nn) += density(aa,bb) * (2.0 * h2e(abmn) - h2e(anmb));
            }
            fock(nn,mm) = fock(mm,nn);
        }
        if(error4DIIS.size() >= size_DIIS)
        {
            error4DIIS.erase(error4DIIS.begin());
            error4DIIS.push_back(evaluateErrorDIIS(fock,density));
            fock4DIIS.erase(fock4DIIS.begin());
            fock4DIIS.push_back(fock);
        }
        else
        {
            error4DIIS.push_back(evaluateErrorDIIS(fock,density));
            fock4DIIS.push_back(fock);
        }        
    }
    if(!converged)
    {
        cout << "ERROR: SCF did NOT converge!" << endl;
        exit(99);
    }
}








/*********************************************************/
/**********    Member functions of class UHF    **********/
/*********************************************************/

UHF::UHF(const int& nelec_a_, const int& nelec_b_, const int& size_basis_):
SCF(nelec_a_, nelec_b_, size_basis_)
{
    if(nelec_a_ == nelec_b_)
    {
        cout << "Warning: UHF is called when nelec_a == nelec_b!" << endl;
    }
}

UHF::UHF(const MOL& mol_, const MatrixXd& s_, const MatrixXd& t_, const MatrixXd& v_, const VectorXd& eri_, const double& V_RR_):
SCF(mol_, s_, t_, v_, eri_, V_RR_)
{
}

UHF::~UHF()
{
}

void UHF::runSCF()
{
    vector<MatrixXd> error_a4DIIS, error_b4DIIS, fock_a4DIIS, fock_b4DIIS;
    fock_a.resize(size_basis, size_basis);
    fock_b.resize(size_basis, size_basis);
    MatrixXd newDen_a, newDen_b, density_total;
    
    eigensolverG(h1e, overlap_half_i, ene_orb_a, coeff_a);
    density_a = evaluateDensity(coeff_a, nelec_a);
    density_b = evaluateDensity(coeff_a, nelec_b);
    density_total = density_a + density_b;

    for(int iter = 1; iter <= maxIter; iter++)
    {
        if(iter <= 2)
        {
            #pragma omp parallel  for
            for(int mm = 0; mm < size_basis; mm++)
            for(int nn = 0; nn < size_basis; nn++)
            {
                fock_a(mm,nn) = h1e(mm,nn);
                fock_b(mm,nn) = h1e(mm,nn);
                for(int aa = 0; aa < size_basis; aa++)
                for(int bb = 0; bb < size_basis; bb++)
                {
                    int ab, mn, an, mb, abmn, anmb;
                    ab = max(aa,bb)*(max(aa,bb)+1)/2 + min(aa,bb);
                    mn = max(mm,nn)*(max(mm,nn)+1)/2 + min(mm,nn);
                    an = max(aa,nn)*(max(aa,nn)+1)/2 + min(aa,nn);
                    mb = max(mm,bb)*(max(mm,bb)+1)/2 + min(mm,bb);
                    abmn = max(ab,mn)*(max(ab,mn)+1)/2 + min(ab,mn);
                    anmb = max(an,mb)*(max(an,mb)+1)/2 + min(an,mb);
                    fock_a(mm,nn) += density_total(aa,bb) * h2e(abmn) - density_a(aa,bb) * h2e(anmb);
                    fock_b(mm,nn) += density_total(aa,bb) * h2e(abmn) - density_b(aa,bb) * h2e(anmb);
                }
            }
        }
        else
        {
            int tmp_size = fock_a4DIIS.size();
            MatrixXd B_a4DIIS(tmp_size+1,tmp_size+1), B_b4DIIS(tmp_size+1,tmp_size+1);
            VectorXd vec_b(tmp_size+1);
            for(int ii = 0; ii < tmp_size; ii++)
            {    
                for(int jj = 0; jj <= ii; jj++)
                {
                    B_a4DIIS(ii,jj) = (error_a4DIIS[ii].adjoint()*error_a4DIIS[jj])(0,0);
                    B_a4DIIS(jj,ii) = B_a4DIIS(ii,jj);
                    B_b4DIIS(ii,jj) = (error_b4DIIS[ii].adjoint()*error_b4DIIS[jj])(0,0);
                    B_b4DIIS(jj,ii) = B_b4DIIS(ii,jj);
                }
                B_a4DIIS(tmp_size, ii) = -1.0;
                B_a4DIIS(ii, tmp_size) = -1.0;
                B_b4DIIS(tmp_size, ii) = -1.0;
                B_b4DIIS(ii, tmp_size) = -1.0;
                vec_b(ii) = 0.0;
            }
            B_a4DIIS(tmp_size, tmp_size) = 0.0;
            B_b4DIIS(tmp_size, tmp_size) = 0.0;
            vec_b(tmp_size) = -1.0;
            VectorXd C_a = B_a4DIIS.partialPivLu().solve(vec_b);
            VectorXd C_b = B_b4DIIS.partialPivLu().solve(vec_b);
            fock_a = MatrixXd::Zero(size_basis,size_basis);
            fock_b = MatrixXd::Zero(size_basis,size_basis);
            for(int ii = 0; ii < tmp_size; ii++)
            {
                fock_a += C_a(ii) * fock_a4DIIS[ii];
                fock_b += C_b(ii) * fock_b4DIIS[ii];
            }
        }
        eigensolverG(fock_a, overlap_half_i, ene_orb_a, coeff_a);
        eigensolverG(fock_b, overlap_half_i, ene_orb_b, coeff_b);
        newDen_a = evaluateDensity(coeff_a, nelec_a);
        newDen_b = evaluateDensity(coeff_b, nelec_b);
        d_density_a = evaluateChange(density_a, newDen_a);
        d_density_b = evaluateChange(density_b, newDen_b);
        cout << "Iter #" << iter << " maximum density difference :\talpha " << d_density_a << "\t\tbeta " << d_density_b << endl;
        density_a = newDen_a;
        density_b = newDen_b;
        density_total = density_a + density_b;
        if(max(d_density_a,d_density_b) < convControl) 
        {
            converged = true;
            cout << endl << "UHF converges after " << iter << " iterations." << endl << endl;

            cout << "\tOrbital\t\tEnergy_a(in hartree)\t\tEnergy_b(in hartree)\n";
            cout << "\t*******\t\t********************\t\t********************\n";
            for(int ii = 1; ii <= size_basis; ii++)
            {
                cout << "\t" << ii << "\t\t" << setprecision(15) << ene_orb_a(ii - 1) << "\t\t" << ene_orb_b(ii - 1) << endl;
            }

            ene_scf = V_RR;
            for(int ii = 0; ii < size_basis; ii++)
            for(int jj = 0; jj < size_basis; jj++)
            {
                ene_scf += 0.5 * (density_total(ii,jj) * h1e(jj,ii) + density_a(ii,jj) * fock_a(jj,ii) + density_b(ii,jj) * fock_b(jj,ii));
            }
            cout << "Final UHF energy is " << setprecision(15) << ene_scf << " hartree." << endl;
            break;            
        }    
        #pragma omp parallel  for
        for(int mm = 0; mm < size_basis; mm++)
        for(int nn = 0; nn < size_basis; nn++)
        {
            fock_a(mm,nn) = h1e(mm,nn);
            fock_b(mm,nn) = h1e(mm,nn);
            for(int aa = 0; aa < size_basis; aa++)
            for(int bb = 0; bb < size_basis; bb++)
            {
                int ab, mn, an, mb, abmn, anmb;
                ab = max(aa,bb)*(max(aa,bb)+1)/2 + min(aa,bb);
                mn = max(mm,nn)*(max(mm,nn)+1)/2 + min(mm,nn);
                an = max(aa,nn)*(max(aa,nn)+1)/2 + min(aa,nn);
                mb = max(mm,bb)*(max(mm,bb)+1)/2 + min(mm,bb);
                abmn = max(ab,mn)*(max(ab,mn)+1)/2 + min(ab,mn);
                anmb = max(an,mb)*(max(an,mb)+1)/2 + min(an,mb);
                fock_a(mm,nn) += density_total(aa,bb) * h2e(abmn) - density_a(aa,bb) * h2e(anmb);
                fock_b(mm,nn) += density_total(aa,bb) * h2e(abmn) - density_b(aa,bb) * h2e(anmb);
            }
        }

        if(error_a4DIIS.size() >= size_DIIS)
        {
            error_a4DIIS.erase(error_a4DIIS.begin());
            error_a4DIIS.push_back(evaluateErrorDIIS(fock_a,density_a));
            fock_a4DIIS.erase(fock_a4DIIS.begin());
            fock_a4DIIS.push_back(fock_a);
            error_b4DIIS.erase(error_b4DIIS.begin());
            error_b4DIIS.push_back(evaluateErrorDIIS(fock_b,density_b));
            fock_b4DIIS.erase(fock_b4DIIS.begin());
            fock_b4DIIS.push_back(fock_b);
        }
        else
        {
            error_a4DIIS.push_back(evaluateErrorDIIS(fock_a,density_a));
            fock_a4DIIS.push_back(fock_a);
            error_b4DIIS.push_back(evaluateErrorDIIS(fock_b,density_b));
            fock_b4DIIS.push_back(fock_b);
        }       
    }
    if(!converged)
    {
        cout << "ERROR: SCF did NOT converge!" << endl;
        exit(99);
    }
}

