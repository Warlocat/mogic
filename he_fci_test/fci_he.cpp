#include<iostream>
#include<fstream>
#include<iomanip>
#include<Eigen/Dense>

using namespace Eigen;
using namespace std;
inline int fact(const int& N);
double Hij(const int& iii, const int& jjj, const MatrixXi& occ);

MatrixXd h1e;
Matrix<MatrixXd, Dynamic, Dynamic> h2e;

int main()
{
    MatrixXd hamiltonian;
    Matrix<int, Dynamic, Dynamic> occ;
    
    int nelec, norb, size;
    double tmp;
    Vector4i indices;

    ifstream ifs;
    ifs.open("FCIDUMP_pvtz");
        ifs >> norb >> nelec;
        h1e.resize(norb, norb);
        h2e.resize(norb, norb);
        size = norb * norb;
        hamiltonian.resize(size, size);
        occ.resize(size, 2);

        for(int ii = 0; ii < norb; ii++)
        for(int jj = 0; jj < norb; jj++)
        {
            h2e(ii,jj).resize(norb, norb);
        }
        while(1)
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

    int tmp1 = 0;
    for(int ii = 0; ii < sqrt(size); ii++)
    for(int jj = 0; jj < sqrt(size); jj++)
    {
        occ(tmp1, 0) = ii;
        occ(tmp1, 1) = jj;
        tmp1++;
    }

    hamiltonian = 0.0 * hamiltonian;
    for(int ii = 0; ii < size; ii++)
    for(int jj = 0; jj < size; jj++)
    {
        hamiltonian(ii, jj) = Hij(ii,jj,occ);
    }

    SelfAdjointEigenSolver<MatrixXd> solver(hamiltonian);
    cout << fixed << setprecision(8) <<  solver.eigenvalues().block(0,0,10,1) << endl;   
    return 0;
}

double Hij(const int& iii, const int& jjj, const MatrixXi& occ)
{
    double res = 0.0;
    int oi1 = occ(iii, 0), oi2 = occ(iii, 1), oj1 = occ(jjj, 0), oj2 = occ(jjj, 1);
    if(oi1 == oj1 && oi2 == oj2)
    {
        res += h1e(oi1, oi1) + h1e(oi2, oi2);
        res += h2e(oi2, oi2)(oi1, oi1);
    }
    else if(oi1 != oj1 && oi2 != oj2)
    {
        res += h2e(oi1, oj1)(oi2, oj2);
    }
    else if(oi1 != oj1 && oi2 == oj2)
    {
        res += h1e(oi1, oj1);
        res += h2e(oi1, oj1)(oi2, oi2);
    }
    else if(oi1 == oj1 && oi2 != oj2)
    {
        res += h1e(oi2, oj2);
        res += h2e(oi2, oj2)(oi1, oi1);
    }
    
    
    return res;
}