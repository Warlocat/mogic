#include<iostream>
#include<fstream>
#include<string>
#include<Eigen/Dense>
#include<iomanip>
#include<cmath>
#include<ctime>
#include"gto.h"
#include"scf.h"
using namespace Eigen;
using namespace std;



int main()
{
    int charge, spin;
    string atomName, basisSet, flags, jobs;
    clock_t startTime, endTime;
    
    ifstream ifs;
    ifs.open("input");
        ifs >> atomName >> basisSet >> charge >> flags >> spin >> flags >> jobs;
    ifs.close();
    GTO gto_test(atomName, basisSet, charge, spin);

    int size = gto_test.size_gtoc;
    cout << "basis size: " << size << endl;

    MatrixXd overlap = gto_test.get_h1e("overlap");
    
    startTime = clock();
    MatrixXd h1e = gto_test.get_h1e("h1e");
    endTime = clock();
    cout << "1e integrals finished in " << (endTime - startTime) / (double)CLOCKS_PER_SEC << " seconds." << endl;

    startTime = clock();
    Matrix<MatrixXd, -1, -1> h2e = gto_test.get_h2e();
    endTime = clock();
    cout << "2e integrals finished in " << (endTime - startTime) / (double)CLOCKS_PER_SEC << " seconds." << endl;

    gto_test.writeIntegrals(overlap, h1e, h2e, "integral.txt");

    if(jobs == "SCF")
    {
        SCF scf_test(gto_test.nelec_a, gto_test.nelec_b, gto_test.size_gtoc);
        startTime = clock();
        scf_test.runSCF();
        endTime = clock();
        cout << "HF scf finished in " << (endTime - startTime) / (double)CLOCKS_PER_SEC << " seconds." << endl; 
    }
    

    return 0;
}



