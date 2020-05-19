#ifndef SCF_H_
#define SCF_H_

#include<Eigen/Dense>
#include<complex>
#include<string>
using namespace std;
using namespace Eigen;

class SCF
{
private:
    int nelec_a, nelec_b, size_basis;
    double d_density_a, d_density_b;
    MatrixXd overlap, overlap_half_i, h1e, fock_a, fock_b, density_a, density_b;
    Matrix<MatrixXd,-1,-1> h2e;
public:
    MatrixXd coeff_a, coeff_b;
    bool converged = false;
    VectorXd ene_orb_a, ene_orb_b;
    int maxIter = 200;
    double ene_hf, convControl = 1e-10;

    SCF(const int& nelec_a_, const int& nelec_b_, const int& size_basis_);
    ~SCF();

    void readIntegrals(const string& filename);
    MatrixXd inverse_half(const MatrixXd& inputM);
    MatrixXd evaluateDensity(const MatrixXd& coeff_, const int& nocc, const bool& spherical = false);
    void eigensolverG(const MatrixXd& inputM, const MatrixXd& s_h_i, VectorXd& values, MatrixXd& vectors);

    void runSCF();
};


#endif