#ifndef MOL_H_
#define MOL_H_

#include<Eigen/Dense>
#include<string>
using namespace std;
using namespace Eigen;

const double ATOMIC_MASS_UNIT = 1.660539040/9.1093826*10000, au2cm_1 = 219474.63, ang2bohr = 1.8897161646320724;

class MOL
{
private:
    /* Get mass for atomic number an */
    double massInit(const int& an);
    double innerP(const VectorXd& v1, const VectorXd& v2);
    Vector3d crossP(const Vector3d& v1, const Vector3d& v2);

public:
    int Natom, charge, nelec, nelec_a, nelec_b, spin;
    string subgroup = "C1", unit = "BOHR";
    VectorXi atomList;
    Matrix<Vector3d,-1,1> coordCart;
    Matrix<string,-1,1> basisSet;
    VectorXd mass;

    MOL();
    ~MOL();

    /* Input and Output */
    void readInput(const string& filename);
    void readxyz(const string& filename);
    void writexyz(const string& filename);
    
    /* transfer coordCart to center of mass coordinates */
    void cart2COM();
    
    /* Calculate bond length r_ij between i and j atom */
    double length(const int& i, const int& j);
    /* Calculate bond angle phi_ijk between i, j and k atom */
    double angle(const int& i, const int& j, const int& k);
    /* Calculate sine of out-of-plane angle theta_ijkl between i, j, k and l atom */
    double angle_ofp(const int& i, const int& j, const int& k, const int& l);
    /* Calculate torsion angle tau_ijkl between i, j, k and l atom */
    double torsion(const int& i, const int& j, const int& k, const int& l);

    /* Calculate Moments of Inertia */
    Matrix3d rotationI();
};

#endif