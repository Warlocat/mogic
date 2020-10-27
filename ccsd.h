#ifndef CCSD_H_
#define CCSD_H_

#include<Eigen/Dense>
#include<complex>
#include<string>
using namespace std;
using namespace Eigen;


class CCSD
{
protected:
    int n_occ, n_vir;
    MatrixXd t1, t2, D1, D2, t1_new, t2_new;
    VectorXd h2e_mo_so, ene_mo_so;
    double ****tau, ****tau_tilde, **F_vv, **F_oo, **F_ov;
    double ****W_oooo, ****W_vvvv, ****W_ovvo;
    inline double h2e_dirac_so(const int& ii, const int& jj, const int& kk, const int& ll);
    void memoryAllocation();
    void memoryDeallocation();
    void evaluate_tau();
    void evaluate_W_F();
    void evaluate_t1t2New();
    static double evaluateChange(const MatrixXd& M1, const MatrixXd& M2);
    double evaluate_ene_ccsd();

public:
    bool converged = false;
    int maxIter = 200;
    double ene_ccsd, convControl = 1e-12;

    CCSD(const int& n_occ_, const int& n_vir_, const VectorXd& h2e_mo_so_, const VectorXd& ene_mo_);
    virtual ~CCSD();

    void runCCSD();
    void runCCSD_pT();
};

#endif
