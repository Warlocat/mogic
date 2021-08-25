#ifndef CCSD_H_
#define CCSD_H_

#include<Eigen/Dense>
#include<complex>
#include<string>
#include<tblis/tblis.h>
using namespace std;
using namespace Eigen;
using namespace tblis;


class CCSD
{
protected:
    int n_occ, n_vir;
    tensor<double> t1, t2, D1, D2, t1_new, t2_new, tau, tau_tilde, F_vv, F_oo, F_ov, W_oooo, W_vvvv, W_ovvo, h2e_oooo, h2e_ooov, h2e_oovv, h2e_ovvv, h2e_vvvv, h2e_ovvo;
    // MatrixXd t1, t2, D1, D2, t1_new, t2_new;
    VectorXd h2e_so_antiSym, ene_so;
    inline double h2e_dirac_so(const int& ii, const int& jj, const int& kk, const int& ll);
    void memoryAllocation();
    void memoryDeallocation();
    void evaluate_tau();
    void evaluate_W_F();
    void evaluate_t1t2New();
    static double evaluateChange_t1(const tensor<double>& M1, const tensor<double>& M2);
    static double evaluateChange_t2(const tensor<double>& M1, const tensor<double>& M2);
    double evaluate_ene_ccsd();

    static tensor<double> evaluateErrorDIIS(const tensor<double>& t1_old, const tensor<double>& t1_new, const tensor<double>& t2_old, const tensor<double>& t2_new);

public:
    bool converged = false;
    int maxIter = 150, size_DIIS = 8;
    double ene_ccsd, ene_pT, convControl = 1e-10;

    /* Constructors for RHF and UHF */
    CCSD(const int& n_occ_, const int& n_vir_, const VectorXd& h2e_so_antiSym_, const VectorXd& ene_mo_);
    CCSD(const int& n_occ_, const int& n_vir_, const VectorXd& h2e_so_antiSym_, const VectorXd& ene_mo_a,  const VectorXd& ene_mo_b);
    virtual ~CCSD();

    void runCCSD();
    void runCCSD_pT();
};

#endif
