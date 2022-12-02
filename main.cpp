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
    VectorXd int2e_eri_so_antiSym = integralTransfermation_SO_antiSym(int2e_eri, rhf_test.coeff);
    CIS cis_test(mol_test.nelec, size_basis*2 - mol_test.nelec, 1, int2e_eri_so_antiSym, rhf_test.ene_orb);
    // CCSD ccsd_test(mol_test.nelec, size_basis*2 - mol_test.nelec, int2e_eri_so_antiSym, rhf_test.ene_orb);

    // UHF uhf_test(mol_test, int1e_s, int1e_t, int1e_v, int2e_eri, e_nuc);
    // uhf_test.runSCF();
    // VectorXd int2e_eri_so_antiSym = integralTransfermation_SO_antiSym(int2e_eri, uhf_test.coeff_a, uhf_test.coeff_b);
    // CCSD ccsd_test(mol_test.nelec, size_basis*2 - mol_test.nelec, int2e_eri_so_antiSym, uhf_test.ene_orb_a, uhf_test.ene_orb_b);
    
    
    // ccsd_test.runCCSD();
    // ccsd_test.runCCSD_pT();

    cis_test.runCIS();
    cout << cis_test.ene_CIS(0) << endl;
    CIS cis_test2(mol_test.nelec, size_basis*2 - mol_test.nelec, 1, int2e_eri_so_antiSym, rhf_test.ene_orb);
    cis_test2.davidsonLiu = false;
    cis_test2.runCIS();
    cout << cis_test2.ene_CIS(0) << endl;
    // cis_test.runTDHF("explicit");
    // cout << 27.2114*cis_test.ene_CIS << endl;

    return 0;
}
