#include<iostream>
#include<fstream>
#include<string>
#include<Eigen/Dense>
#include<iomanip>
#include<cmath>
#include<ctime>
#include"gto.h"
using namespace Eigen;
using namespace std;

double f1(int l, double r2, double a);
double f2(int l, double r2, double a);

int main()
{
    GTO test;
    clock_t startTime, endTime;
    test.readBasis("Cu", "ccpvdz");


    int size = test.size;   
    test.normalization();
    Matrix<MatrixXd, -1, -1> h2e;
    
    startTime = clock();
    MatrixXd h1e = test.get_h1e("h1e");
    endTime = clock();
    cout << "1e integrals finished in " << (endTime - startTime) / (double)CLOCKS_PER_SEC << " seconds." << endl;

    startTime = clock();
    h2e = test.get_h2e();
    endTime = clock();
    cout << "2e integrals finished in " << (endTime - startTime) / (double)CLOCKS_PER_SEC << " seconds." << endl;

    ofstream ofs;
    ofs.open("h2e.txt");
        for(int ii = 0; ii < size; ii++)
        for(int jj = 0; jj <= ii; jj++)
        for(int kk = 0; kk < size; kk++)
        for(int ll = 0; ll <= kk; ll++)
        {
            if(abs(h2e(ii,jj)(kk,ll)) > 1e-12)  ofs << setprecision(16) << h2e(ii,jj)(kk,ll) << "\t" << ii+1 << "\t" << jj+1 << "\t" << kk+1 << "\t" << ll+1 << endl;
        }
    ofs.close();

    return 0;
}



