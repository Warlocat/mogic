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
    h2e = test.get_h2e();
    ofstream ofs;
    ofs.open("h2e.txt");
        for(int ii = 0; ii < size; ii++)
        for(int jj = 0; jj < size; jj++)
        for(int kk = 0; kk < size; kk++)
        for(int ll = 0; ll < size; ll++)
        {
            if(h2e(ii,jj)(kk,ll) > 1e-8)  ofs << ii << "\t" << jj << "\t" << kk << "\t" << ll << "\t" << h2e(ii,jj)(kk,ll) << endl;
        }
    ofs.close();
    
    
    

    return 0;
}