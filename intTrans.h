#include<Eigen/Dense>
#include<string>
using namespace std;
using namespace Eigen;


/*
    Integral transformation from AO to restricted MO or spin orbitals
*/
//  RHF 
VectorXd integralTransfermation_MO(const VectorXd& h2e_ao, const MatrixXd& coeff, const string alg = "smart");
VectorXd integralTransfermation_SO(const VectorXd& h2e_ao, const MatrixXd& coeff, const string alg = "smart");
VectorXd integralTransfermation_SO_antiSym(const VectorXd& h2e_ao, const MatrixXd& coeff, const string alg = "smart");
//  UHF 
Matrix<VectorXd,4,1> integralTransfermation_MO(const VectorXd& h2e_ao, const MatrixXd& coeff_a, const MatrixXd& coeff_b, const string alg = "smart");
VectorXd integralTransfermation_SO(const VectorXd& h2e_ao, const MatrixXd& coeff_a, const MatrixXd& coeff_b, const string alg = "smart");
VectorXd integralTransfermation_SO_antiSym(const VectorXd& h2e_ao, const MatrixXd& coeff_a, const MatrixXd& coeff_b, const string alg = "smart");

/*
    Evaluate MP2 energy
*/
// RHF
double get_energy_MP2(const VectorXd& h2e_mo, const VectorXd& ene_mo, const int& nelec_a, const int& nelec_b);
// RHF
double get_energy_MP2(const VectorXd& h2e_mo_a, const VectorXd& h2e_mo_b, const VectorXd& ene_mo_a, const VectorXd& ene_mo_b, const int& nelec_a, const int& nelec_b);