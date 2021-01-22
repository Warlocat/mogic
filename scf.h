#ifndef SCF_H_
#define SCF_H_

#include<Eigen/Dense>
#include<complex>
#include<string>
#include"mol.h"
using namespace std;
using namespace Eigen;


class SCF
{
protected:
    int nelec_a, nelec_b, size_basis;
    double V_RR;
    MatrixXd overlap, overlap_half_i, h1e;
    VectorXd h2e;

    static void readIntegrals_1e(MatrixXd& int_1e, const string& filename);
    static void readIntegrals_2e(VectorXd& int_2e, const string& filename);
    static double evaluateChange(const MatrixXd& M1, const MatrixXd& M2);
    static MatrixXd matrix_half_inverse(const MatrixXd& inputM);
    static MatrixXd matrix_half(const MatrixXd& inputM);
    static MatrixXd evaluateDensity(const MatrixXd& coeff_, const int& nocc, const bool& spherical = false);
    static void eigensolverG(const MatrixXd& inputM, const MatrixXd& s_h_i, VectorXd& values, MatrixXd& vectors);

public:
    bool converged = false;
    int maxIter = 200;
    double ene_scf, convControl = 1e-12;

    SCF(const MOL& mol_, const MatrixXd& s_, const MatrixXd& t_, const MatrixXd& v_, const VectorXd& eri_, const double& V_RR_);
    SCF(const int& nelec_a_, const int& nelec_b_, const int& size_basis_);
    virtual ~SCF();

    virtual void runSCF() = 0;
    VectorXd get_h2e_vector();
};

class RHF: public SCF
{
private:
    MatrixXd density, fock;
    double d_density;
public:
    MatrixXd coeff;
    VectorXd ene_orb;

    RHF(const MOL& mol_, const MatrixXd& s_, const MatrixXd& t_, const MatrixXd& v_, const VectorXd& eri_, const double& V_RR_);
    RHF(const int& nelec_a_, const int& nelec_b_, const int& size_basis_);
    virtual ~RHF();
    virtual void runSCF();
};


class UHF: public SCF
{
private:
    MatrixXd density_a, density_b, fock_a, fock_b;
    double d_density_a, d_density_b;
public:
    MatrixXd coeff_a, coeff_b;
    VectorXd ene_orb_a, ene_orb_b;

    UHF(const MOL& mol_, const MatrixXd& s_, const MatrixXd& t_, const MatrixXd& v_, const VectorXd& eri_, const double& V_RR_);
    UHF(const int& nelec_a_, const int& nelec_b_, const int& size_basis_);
    virtual ~UHF();
    virtual void runSCF();
};



#endif
