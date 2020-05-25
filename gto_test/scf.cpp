#include"scf.h"
#include<omp.h>
#include<iostream>
#include<iomanip>
#include<fstream>


SCF* scf_init(const GTO& gto_, const string& h2e_file)
{
    SCF* ptr_scf;
    if(gto_.nelec_a == gto_.nelec_b)
    {
        cout << "RHF program is used." << endl << endl;
        ptr_scf = new RHF(gto_, h2e_file);
    }
    else
    {
        cout << "UHF program is used." << endl << endl;
        ptr_scf = new UHF(gto_, h2e_file);
    }
    return ptr_scf;
}


/*********************************************************/
/**********    Member functions of class SCF    **********/
/*********************************************************/


SCF::SCF(const GTO& gto_, const string& h2e_file):
nelec_a(gto_.nelec_a), nelec_b(gto_.nelec_b), size_basis(gto_.size_gtoc)
{
    overlap = gto_.get_h1e("overlap");
    h1e = gto_.get_h1e("h1e");

    h2e.resize(size_basis * (size_basis + 1) / 2, size_basis * (size_basis + 1) / 2); h2e = h2e * 0.0;

    readIntegrals(h2e_file);
    overlap_half_i = inverse_half(overlap);
}

SCF::~SCF()
{
}


void SCF::readIntegrals(const string& filename)
{
    ifstream ifs;
    double tmp;
    VectorXi indices(4);
    ifs.open(filename);
        while(true)
        {
            ifs >> tmp >> indices(0) >> indices(1) >> indices(2) >> indices(3);
            if(indices(0) == 0) break;
            else
            {
                int ij = (indices(0) - 1) * indices(0) / 2 + indices(1) - 1, kl = (indices(2) - 1) * indices(2) / 2 + indices(3) - 1;
                h2e(ij,kl) = tmp;
                h2e(kl,ij) = tmp;
            }
        }            
    ifs.close();
}


MatrixXd SCF::inverse_half(const MatrixXd& inputM)
{
    SelfAdjointEigenSolver<MatrixXd> solver(inputM);
    MatrixXd eigenvalues(size_basis, size_basis);
    MatrixXd eigenvectors = solver.eigenvectors();
    
    for(int ii = 0; ii < size_basis; ii++)
    {
        eigenvalues(ii,ii) = solver.eigenvalues()(ii);
        if(eigenvalues(ii,ii) < 0)
        {
            cout << "ERROR: overlap matrix has negative eigenvalues!" << endl;
            exit(99);
        }
        else
        {
            eigenvalues(ii,ii) = 1.0 / sqrt(eigenvalues(ii,ii));
        }
    }
    return eigenvectors * eigenvalues * eigenvectors.transpose(); 
}

MatrixXd SCF::evaluateDensity(const MatrixXd& coeff_, const int& occ, const bool& spherical)
{
    MatrixXd den(size_basis, size_basis);
    if(!spherical)
    {
        for(int aa = 0; aa < size_basis; aa++)
        for(int bb = 0; bb < size_basis; bb++)
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
}




/*********************************************************/
/**********    Member functions of class RHF    **********/
/*********************************************************/

RHF::RHF(const GTO& gto_, const string& h2e_file):
SCF(gto_, h2e_file)
{
}

RHF::~RHF()
{
}

void RHF::runSCF()
{
    fock.resize(size_basis,size_basis);

    MatrixXd newDen;
    eigensolverG(h1e, overlap_half_i, ene_orb, coeff);
    density = evaluateDensity(coeff, nelec_a);

    for(int iter = 1; iter <= maxIter; iter++)
    {
        #pragma omp parallel  for
        for(int mm = 0; mm < size_basis; mm++)
        for(int nn = 0; nn < size_basis; nn++)
        {
            fock(mm,nn) = h1e(mm,nn);
            for(int aa = 0; aa < size_basis; aa++)
            for(int bb = 0; bb < size_basis; bb++)
            {
                int ab, mn, an, mb;
                ab = max(aa,bb)*(max(aa,bb)+1)/2 + min(aa,bb);
                mn = max(mm,nn)*(max(mm,nn)+1)/2 + min(mm,nn);
                an = max(aa,nn)*(max(aa,nn)+1)/2 + min(aa,nn);
                mb = max(mm,bb)*(max(mm,bb)+1)/2 + min(mm,bb);
                fock(mm,nn) += density(aa,bb) * (2.0 * h2e(ab,mn) - h2e(an, mb));
            }
        }
        eigensolverG(fock, overlap_half_i, ene_orb, coeff);
        newDen = evaluateDensity(coeff, nelec_a);
        d_density = abs((density - newDen).maxCoeff());
        cout << "Iter #" << iter << " maximum densitydifference: " << d_density << endl;
        
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

            ene_scf = 0.0;
            for(int ii = 0; ii < size_basis; ii++)
            for(int jj = 0; jj < size_basis; jj++)
            {
                ene_scf += density(ii,jj) * (h1e(jj,ii) + fock(jj,ii));
            }
            cout << "Final RHF energy is " << setprecision(15) << ene_scf << " hartree." << endl;
            break;            
        }           
    }
}








/*********************************************************/
/**********    Member functions of class UHF    **********/
/*********************************************************/

UHF::UHF(const GTO& gto_, const string& h2e_file):
SCF(gto_, h2e_file)
{
}

UHF::~UHF()
{
}

void UHF::runSCF()
{
    fock_a.resize(size_basis, size_basis);
    fock_b.resize(size_basis, size_basis);
    MatrixXd newDen_a, newDen_b, density_total;
    
    eigensolverG(h1e, overlap_half_i, ene_orb_a, coeff_a);
    density_a = evaluateDensity(coeff_a, nelec_a);
    density_b = evaluateDensity(coeff_a, nelec_b);
    density_total = density_a + density_b;

    for(int iter = 1; iter <= maxIter; iter++)
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
                int ab, mn, an, mb;
                ab = max(aa,bb)*(max(aa,bb)+1)/2 + min(aa,bb);
                mn = max(mm,nn)*(max(mm,nn)+1)/2 + min(mm,nn);
                an = max(aa,nn)*(max(aa,nn)+1)/2 + min(aa,nn);
                mb = max(mm,bb)*(max(mm,bb)+1)/2 + min(mm,bb);
                fock_a(mm,nn) += density_total(aa,bb) * h2e(ab,mn) - density_a(aa,bb) * h2e(an, mb);
                fock_b(mm,nn) += density_total(aa,bb) * h2e(ab,mn) - density_b(aa,bb) * h2e(an, mb);
            }
        }
        eigensolverG(fock_a, overlap_half_i, ene_orb_a, coeff_a);
        eigensolverG(fock_b, overlap_half_i, ene_orb_b, coeff_b);
        newDen_a = evaluateDensity(coeff_a, nelec_a);
        newDen_b = evaluateDensity(coeff_b, nelec_b);

        d_density_a = abs((density_a - newDen_a).maxCoeff());
        d_density_b = abs((density_b - newDen_b).maxCoeff());
        cout << "Iter #" << iter << " maximum density difference :\talpha " << d_density_a << "\t\tbeta " << d_density_b << endl;
            
        // density_a = 0.5 * (newDen_a + density_a);
        // density_b = 0.5 * (newDen_b + density_b);
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

            ene_scf = 0.0;
            for(int ii = 0; ii < size_basis; ii++)
            for(int jj = 0; jj < size_basis; jj++)
            {
                ene_scf += 0.5 * (density_total(ii,jj) * h1e(jj,ii) + density_a(ii,jj) * fock_a(jj,ii) + density_b(ii,jj) * fock_b(jj,ii));
            }
            cout << "Final UHF energy is " << setprecision(15) << ene_scf << " hartree." << endl;
            break;            
        }           
    }
}


