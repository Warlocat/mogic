#include<iostream>
#include<fstream>
#include<string>
#include<Eigen/Dense>
#include<iomanip>
#include<cmath>
#include<ctime>
#include<memory>
#include"gto.h"
#include"scf.h"
#include"x2c.h"
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

    startTime = clock();
    const MatrixXd h2e = gto_test.get_h2e();
    endTime = clock();
    cout << "2e integrals finished in " << (endTime - startTime) / (double)CLOCKS_PER_SEC << " seconds." << endl;

    gto_test.writeIntegrals(h2e, "h2e_"+atomName+".txt");

    if(jobs == "SCF")
    {
        shared_ptr<SCF> ptr_scf(scf_init(gto_test,"h2e_"+atomName+".txt"));
        
        startTime = clock();
        ptr_scf->runSCF();
        endTime = clock();
        cout << "HF scf finished in " << (endTime - startTime) / (double)CLOCKS_PER_SEC << " seconds." << endl; 
    }
    else if(jobs == "SCF_SFX2C1E")
    {
        shared_ptr<SCF> ptr_scf(scf_init(gto_test,"h2e_"+atomName+".txt", "sfx2c1e"));
        
        startTime = clock();
        ptr_scf->runSCF();
        endTime = clock();
        cout << "HF (SFX2C 1E) scf finished in " << (endTime - startTime) / (double)CLOCKS_PER_SEC << " seconds." << endl; 
    }
    


    return 0;
}



