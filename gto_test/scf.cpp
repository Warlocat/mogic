#include"scf.h"
#include<iostream>
#include<iomanip>
#include<fstream>

SCF::SCF(const int& nelec_a_, const int& nelec_b_, const int& size_basis_):
nelec_a(nelec_a_), nelec_b(nelec_b_), size_basis(size_basis_)
{
    overlap.resize(size_basis, size_basis);
    fock_a.resize(size_basis, size_basis);
    h1e.resize(size_basis, size_basis); h1e = h1e * 0.0;
    h2e.resize(size_basis, size_basis);
    for(int ii = 0; ii < size_basis; ii++)
    for(int jj = 0; jj < size_basis; jj++)
    {
        h2e(ii,jj).resize(size_basis,size_basis);
        h2e(ii,jj) = h2e(ii,jj) * 0.0;
    }
    readIntegrals("integral.txt");
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
        for(int ii = 0; ii < size_basis; ii++)
        for(int jj = 0; jj <= ii; jj++)
        {
            ifs >> overlap(ii,jj);
            overlap(jj,ii) = overlap(ii,jj);
        }
        
        while(true)
        {
            ifs >> tmp >> indices(0) >> indices(1) >> indices(2) >> indices(3);
            if(indices(0) == 0) break;
            else if (indices(2) == 0)
            {
                h1e(indices(0) - 1, indices(1) - 1) = tmp;
                h1e(indices(1) - 1, indices(0) - 1) = tmp;
            }
            else
            {
                h2e(indices(0) - 1, indices(1) - 1)(indices(2) - 1, indices(3) - 1) = tmp;
                h2e(indices(1) - 1, indices(0) - 1)(indices(2) - 1, indices(3) - 1) = tmp;
                h2e(indices(0) - 1, indices(1) - 1)(indices(3) - 1, indices(2) - 1) = tmp;
                h2e(indices(1) - 1, indices(0) - 1)(indices(3) - 1, indices(2) - 1) = tmp;
                h2e(indices(2) - 1, indices(3) - 1)(indices(0) - 1, indices(1) - 1) = tmp;
                h2e(indices(2) - 1, indices(3) - 1)(indices(1) - 1, indices(0) - 1) = tmp;
                h2e(indices(3) - 1, indices(2) - 1)(indices(0) - 1, indices(1) - 1) = tmp;
                h2e(indices(3) - 1, indices(2) - 1)(indices(1) - 1, indices(0) - 1) = tmp;
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


void SCF::runSCF()
{
    if(nelec_a == nelec_b)
    {
        MatrixXd newDen;
        cout << "RHF program is used." << endl << endl;
        eigensolverG(h1e, overlap_half_i, ene_orb_a, coeff_a);
        density_a = evaluateDensity(coeff_a, nelec_a);

        for(int iter = 1; iter <= maxIter; iter++)
        {
            for(int mm = 0; mm < size_basis; mm++)
            for(int nn = 0; nn < size_basis; nn++)
            {
                fock_a(mm,nn) = h1e(mm,nn);
                for(int aa = 0; aa < size_basis; aa++)
                for(int bb = 0; bb < size_basis; bb++)
                    fock_a(mm,nn) += density_a(aa,bb) * (2.0 * h2e(aa,bb)(mm,nn) - h2e(aa,nn)(mm,bb));
            }
            eigensolverG(fock_a, overlap_half_i, ene_orb_a, coeff_a);
            newDen = evaluateDensity(coeff_a, nelec_a);
            d_density_a = abs((density_a - newDen).maxCoeff());
            cout << "Iter #" << iter << " maximum density difference: " << d_density_a << endl;
            
            density_a = newDen;
            if(d_density_a < convControl) 
            {
                converged = true;
                cout << endl << "RHF converges after " << iter << " iterations." << endl << endl;

                cout << "\tOrbital\t\tEnergy(in hartree)\n";
                cout << "\t*******\t\t******************\n";
                for(int ii = 1; ii <= size_basis; ii++)
                {
                    cout << "\t" << ii << "\t\t" << setprecision(15) << ene_orb_a(ii - 1) << endl;
                }

                ene_hf = 0.0;
                for(int ii = 0; ii < size_basis; ii++)
                for(int jj = 0; jj < size_basis; jj++)
                {
                    ene_hf += density_a(ii,jj) * (h1e(jj,ii) + fock_a(jj,ii));
                }
                cout << "Final RHF energy is " << setprecision(15) << ene_hf << " hartree." << endl;
                break;            
            }           
        }
    }
    else
    {
        fock_b.resize(size_basis, size_basis);
        MatrixXd newDen_a, newDen_b, density_total;
        cout << "UHF program is used (ROHF is not supported now)." << endl << endl;
        eigensolverG(h1e, overlap_half_i, ene_orb_a, coeff_a);
        density_a = evaluateDensity(coeff_a, nelec_a);
        density_b = evaluateDensity(coeff_a, nelec_b);
        density_total = density_a + density_b;

        for(int iter = 1; iter <= maxIter; iter++)
        {
            for(int mm = 0; mm < size_basis; mm++)
            for(int nn = 0; nn < size_basis; nn++)
            {
                fock_a(mm,nn) = h1e(mm,nn);
                fock_b(mm,nn) = h1e(mm,nn);
                for(int aa = 0; aa < size_basis; aa++)
                for(int bb = 0; bb < size_basis; bb++)
                {
                    fock_a(mm,nn) += density_total(aa,bb) * h2e(aa,bb)(mm,nn) - density_a(aa,bb) * h2e(aa,nn)(mm,bb);
                    fock_b(mm,nn) += density_total(aa,bb) * h2e(aa,bb)(mm,nn) - density_b(aa,bb) * h2e(aa,nn)(mm,bb);
                }    
            }
            eigensolverG(fock_a, overlap_half_i, ene_orb_a, coeff_a);
            eigensolverG(fock_b, overlap_half_i, ene_orb_b, coeff_b);
            newDen_a = evaluateDensity(coeff_a, nelec_a);
            newDen_b = evaluateDensity(coeff_b, nelec_b);
            d_density_a = abs((density_a - newDen_a).maxCoeff());
            d_density_b = abs((density_b - newDen_b).maxCoeff());
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

                ene_hf = 0.0;
                for(int ii = 0; ii < size_basis; ii++)
                for(int jj = 0; jj < size_basis; jj++)
                {
                    ene_hf += 0.5 * (density_total(ii,jj) * h1e(jj,ii) + density_a(ii,jj) * fock_a(jj,ii) + density_b(ii,jj) * fock_b(jj,ii));
                }
                cout << "Final UHF energy is " << setprecision(15) << ene_hf << " hartree." << endl;
                break;            
            }           
        }
    }

    if(!converged)
    {
        cout << "ERROR: HF scf reached max iteration without convergence!" << endl;
        exit(99);
    }
}