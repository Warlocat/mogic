#include<Eigen/Dense>
#include<string>
using namespace std;
using namespace Eigen;


/*
    Integral transformation from AO to restricted MO or spin orbitals
*/
VectorXd integralTransfermation(const VectorXd& h2e_ao, const MatrixXd& coeff, const string alg = "smart");
VectorXd integralTransfermation_spatial2spin(const VectorXd& h2e_mo, const int& size_basis);


/*
    Evaluate MP2 energy
*/
double get_energy_MP2(const VectorXd& h2e_mo, const VectorXd& ene_mo, const int& nelec_a, const int& nelec_b);