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
    test.normalization();
    // double h1e = test.get_h1e();
    MatrixXd s = test.get_overlap();
    cout << s << endl;
    

    return 0;
}