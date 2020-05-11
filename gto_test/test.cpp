#include<fstream>
#include<iostream>
#include<Eigen/Dense>
using namespace std;
using namespace Eigen;

int main()
{
	MatrixXd h1e;
	Matrix<MatrixXd, Dynamic, Dynamic> h2e;
	MatrixXd fock;
    
    int size = 14;
    double tmp;
    Vector4i indices;

    ifstream ifs;
    ifs.open("FCIDUMP_C_dz");
        h2e.resize(size, size);
        fock.resize(size, size);
        for(int ii = 0; ii < size; ii++)
        for(int jj = 0; jj < size; jj++)
        {
            h2e(ii,jj).resize(size, size);
        }

        for(int ii = 0; ii < 1447 ; ii++)
        {
            ifs >> tmp >> indices(0) >> indices(1) >> indices(2) >> indices(3);
            h2e(indices(0) - 1, indices(1) - 1)(indices(2) - 1, indices(3) - 1) = tmp;
            h2e(indices(1) - 1, indices(0) - 1)(indices(2) - 1, indices(3) - 1) = tmp;
            h2e(indices(0) - 1, indices(1) - 1)(indices(3) - 1, indices(2) - 1) = tmp;
            h2e(indices(1) - 1, indices(0) - 1)(indices(3) - 1, indices(2) - 1) = tmp;
            h2e(indices(2) - 1, indices(3) - 1)(indices(0) - 1, indices(1) - 1) = tmp;
            h2e(indices(2) - 1, indices(3) - 1)(indices(1) - 1, indices(0) - 1) = tmp;
            h2e(indices(3) - 1, indices(2) - 1)(indices(0) - 1, indices(1) - 1) = tmp;
            h2e(indices(3) - 1, indices(2) - 1)(indices(1) - 1, indices(0) - 1) = tmp;          
        }
    ifs.close();
	for(int ii = 0; ii < size; ii++)
	for(int jj = 0; jj < size; jj++)
	{
		fock(ii,jj) = 0.0;
		for(int kk = 0; kk < size; kk++)
		{
			fock(ii,jj) += h2e(ii,jj)(kk,kk) - h2e(ii,kk)(kk,jj);
		}
	}


	SelfAdjointEigenSolver<MatrixXd> solver(fock);
    cout << solver.eigenvalues() << endl;   


	return 0;
}

