#include<iostream>
#include<fstream>
#include<omp.h>
#include<iomanip>
#include"mol.h"
#include"molint.h"
#include"scf.h"
#include"ccsd.h"
#include"intTrans.h"
using namespace std;
using namespace Eigen;


int main()
{
    MOL mol_test;
    mol_test.readxyz("h2o.xyz");
    mol_test.readInput("inputFile");
    MOLINT molint_test(mol_test,false);
    MatrixXd int1e_s = molint_test.get_h1e("overlap");
    MatrixXd int1e_t = molint_test.get_h1e("kinetic");
    MatrixXd int1e_v = molint_test.get_h1e("nucV");
    VectorXd int2e_eri = molint_test.get_h2e("eriLLLL");
    double e_nuc = molint_test.get_V_RR();
    int size_basis = int1e_s.rows();
    
    RHF rhf_test(mol_test, int1e_s, int1e_t, int1e_v, int2e_eri, e_nuc);
    rhf_test.runSCF();
    // UHF uhf_test(mol_test, int1e_s, int1e_t, int1e_v, int2e_eri, e_nuc);
    // uhf_test.runSCF();

    VectorXd int2e_eri_mo = integralTransfermation(int2e_eri, rhf_test.coeff);
    VectorXd int2e_eri_so = integralTransfermation_spatial2spin(int2e_eri_mo, size_basis);
    CCSD ccsd_test(mol_test.nelec, size_basis*2 - mol_test.nelec, int2e_eri_so, rhf_test.ene_orb);
    ccsd_test.runCCSD();
    // ccsd_test.runCCSD_pT();

    return 0;
}
