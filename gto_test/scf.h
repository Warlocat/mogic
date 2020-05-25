#ifndef SCF_H_
#define SCF_H_

#include<Eigen/Dense>
#include<complex>
#include<string>
#include"gto.h"
using namespace std;
using namespace Eigen;


class SCF
{
protected:
    int nelec_a, nelec_b, size_basis;
    MatrixXd overlap, overlap_half_i, h1e, h2e;

    void readIntegrals(const string& filename);
    MatrixXd inverse_half(const MatrixXd& inputM);
    MatrixXd evaluateDensity(const MatrixXd& coeff_, const int& nocc, const bool& spherical = false);
    void eigensolverG(const MatrixXd& inputM, const MatrixXd& s_h_i, VectorXd& values, MatrixXd& vectors);
public:
    bool converged = false;
    int maxIter = 200;
    double ene_scf, convControl = 1e-10;

    SCF(const GTO& gto_, const string& h2e_file);
    virtual ~SCF();

    virtual void runSCF() = 0;
};

class RHF: public SCF
{
private:
    MatrixXd density, fock;
    double d_density;
public:
    MatrixXd coeff;
    VectorXd ene_orb;

    RHF(const GTO& gto_, const string& h2e_file);
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

    UHF(const GTO& gto_, const string& h2e_file);
    virtual ~UHF();
    virtual void runSCF();
};


SCF* scf_init(const GTO& gto_, const string& h2e_file);



#endif