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
    test.normalization();
    

    return 0;
}