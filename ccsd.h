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
    MatrixXd t1, t2, D1, D2;
    VectorXd h2e_mo_so, ene_mo_so;
    double tau[][][][], tau_tilde[][][][], F1, F2, F3;
    double W1[][][][], W2[][][][], W3[][][][];
    inline double h2e_dirac_so(const int& ii, const int& jj, const int& kk, const int& ll);
    static void evaluate_tau(double tau[][][][], double tau_tilde[][][][], const MatrixXd& t1, const MatrixXd& t2);

public:
    bool converged = false;
    int maxIter = 200;
    double ene_cc, convControl = 1e-10;

    CCSD(const int& n_occ_, const int& n_vir_, const VectorXd& h2e_mo_so_, const VectorXd& ene_mo_);
    virtual ~CCSD();

    void runCCSD();
};

#endif
