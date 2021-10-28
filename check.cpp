#include<iostream>
#include<fstream>
#include<omp.h>
#include<iomanip>
#include"mol.h"
#include"molint.h"
#include"scf.h"
#include"ccsd.h"
#include"cis.h"
#include"intTrans.h"
using namespace std;
using namespace Eigen;


int main()
{ 
    MOL mol_test;
    mol_test.readxyz("geom.dat");
    mol_test.readInput("inputFile");
    MOLINT molint_test(mol_test,true);
    MatrixXd V = molint_test.get_h1e("nucV");
    for(int ii = 0; ii < V.rows(); ii++)
    for(int jj = 0; jj <= ii; jj++)
    {
        cout << V(ii,jj) << endl;
    }

    return 0;
}
