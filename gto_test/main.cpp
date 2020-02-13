#include<iostream>
#include<fstream>
#include<Eigen/Dense>
#include"gto.h"
using namespace Eigen;
using namespace std;

int main()
{
    GTO test;
    test.readBasis();
    int size = test.size;   
    test.normalization();
    Matrix<MatrixXd, -1, -1> h2e;
    MatrixXd h1e = test.get_h1e("h1e");
    cout << h1e << endl;
    h2e = test.get_h2e();
    ofstream ofs;
    
    ofs.open("h2e.txt");
        for(int ii = 0; ii < size; ii++)
        for(int jj = 0; jj <= ii; jj++)
        for(int kk = 0; kk < size; kk++)
        for(int ll = 0; ll <= kk; ll++)
        {
            if(h2e(ii,jj)(kk,ll) > 1e-8)  ofs << h2e(ii,jj)(kk,ll) << "\t" << ii+1 << "\t" << jj+1 << "\t" << kk+1 << "\t" << ll+1 << endl;
        }
    ofs.close();
    
    
    

    return 0;
}