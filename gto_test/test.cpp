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
#include"gto_spinor.h"
using namespace Eigen;
using namespace std;

/* Global information */
int charge, spin;
string atomName, basisSet, flags, jobs, rel;
bool unc;

/* Read input file and set global variables */
void readInput(const string filename);

int main()
{
    readInput("input");
    clock_t startTime, endTime;
       
    GTO_SPINOR gto_spinor_test(atomName, basisSet, charge, spin);

    MatrixXd spnucsp = gto_spinor_test.get_h1e("s_p_nuc_s_p", unc);
    MatrixXd w_sf = gto_spinor_test.get_h1e("p.Vp", unc), w_sd = gto_spinor_test.get_h1e("i_s_pV_x_p", unc);
    // cout << (w_sf + w_sd - spnucsp).maxCoeff() << endl << endl;
    for(int ii = 0; ii < spnucsp.rows(); ii++)
    for(int jj = 0; jj < spnucsp.cols(); jj++)
        cout << setprecision(16) << spnucsp(ii,jj) << endl;

    

    return 0;
}



void readInput(const string filename)
{
    ifstream ifs;
    ifs.open(filename);
        ifs >> atomName >> flags;
        ifs >> basisSet >> flags;
        ifs >> charge >> flags;
        ifs >> spin >> flags;
        ifs >> jobs >> flags;
        ifs >> rel >> flags;
        ifs >> unc >> flags;
        cout << atomName << endl << basisSet <<endl << charge<< endl << spin << endl << jobs << endl << rel << endl << unc << endl;
    ifs.close();
}