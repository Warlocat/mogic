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
    VectorXd HdotC_CIS(const VectorXd& Cin);

public:
    bool converged = false, davidsonLiu = true;
    int nDL = 20, nroot = 1, maxCycle = 100;
    double convControl = 1e-7;
    VectorXd ene_CIS;
    MatrixXd CIScoeff;

    /* Constructors for RHF and UHF */
    CIS(const int& n_occ_, const int& n_vir_, const int& nroot_, const VectorXd& h2e_so_antiSym_, const VectorXd& ene_mo_);
    CIS(const int& n_occ_, const int& n_vir_, const int& nroot_, const VectorXd& h2e_so_antiSym_, const VectorXd& ene_mo_a,  const VectorXd& ene_mo_b);
    ~CIS();

    void runCIS();
    void runTDHF();

    void DLsolver(const VectorXd& Hdiag, VectorXd (*HdotC)(const VectorXd&), MatrixXd& eigenvectors, VectorXd& eigenvalues, const int& nroot, const int& maxCycle);
    void DLsolver_CIS(const VectorXd& Hdiag, MatrixXd& eigenvectors, VectorXd& eigenvalues, const int& nroot, const int& maxCycle);
};

#endif
