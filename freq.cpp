#include<iostream>
#include<fstream>
#include"mol.h"
using namespace std;

int main()
{
    MOL test;
    test.readxyz("test.xyz");
    MatrixXd hess(3*test.Natom,3*test.Natom);
    ifstream ifs;
    ifs.open("hess.txt");
    for(int ii = 0; ii < hess.rows(); ii++)
    for(int jj = 0; jj < hess.rows(); jj++)
    {
        ifs >> hess(ii, jj);
        hess(ii, jj) = hess(ii, jj) / sqrt(test.mass(ii/3) * test.mass(jj/3)) / ATOMIC_MASS_UNIT;
    }
    ifs.close();
    SelfAdjointEigenSolver<MatrixXd> solver(hess);
    VectorXd eigenv = solver.eigenvalues(), omega(3 * test.Natom);
    for(int ii = 0; ii < eigenv.rows(); ii++)
    {
        if(eigenv(ii) < 1e-6) omega(ii) = 0.0;
        else
        {
            omega(ii) = sqrt(eigenv(ii));
        }
    }

    cout << eigenv << endl;
    cout << omega * 5140.48714 << endl;

    cout << "Congratulation! This program finished successfully!"  << endl;
    return 0;
}