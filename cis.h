#ifndef CIS_H_
#define CIS_H_

#include<Eigen/Dense>
#include<complex>
#include<string>
using namespace std;
using namespace Eigen;


class CIS
{
protected:
    int n_occ, n_vir;
    double ene_scf;
    VectorXd h2e_so_antiSym, ene_so;
    double evaluateH_CIS(const int& i, const int& j, const int& a, const int& b);
    double evaluateH_TDHF_B(const int& i, const int& j, const int& a, const int& b);

public:
    bool converged = false;
    VectorXd ene_CIS;
    MatrixXd CIScoeff;

    /* Constructors for RHF and UHF */
    CIS(const int& n_occ_, const int& n_vir_, const VectorXd& h2e_so_antiSym_, const VectorXd& ene_mo_);
    CIS(const int& n_occ_, const int& n_vir_, const VectorXd& h2e_so_antiSym_, const VectorXd& ene_mo_a,  const VectorXd& ene_mo_b);
    ~CIS();

    void runCIS(const string& alg = "explicit_D_L");
    void runTDHF(const string& alg = "explicit_D_L");
};

#endif
