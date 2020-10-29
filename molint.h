#ifndef MOLINT_H_
#define MOLINT_H_

#include<Eigen/Dense>
#include<complex>
#include<string>
#include<gsl/gsl_sf_gamma.h>
#include"mol.h"
using namespace std;
using namespace Eigen;

double double_factorial(const int& n);
double factorial(const int& n);
double Gamma(double z);
double Boys(double x, int n);


/*
    contracted gtos
*/
struct gto_shell
{
    Vector3d coord;
    VectorXd exp_a;
    MatrixXd coeff;
    int l;
};



class MOLINT
{
private:
    Matrix<gto_shell,-1,1> shell_list;

    MatrixXi l_xyz(const int& l) const;
    void evaluate_l_xyz_coeff();

    double recurrence_E(const int& t, const int& i, const int& j, const double& Ax, const double& Bx, const double& a, const double& b) const;
    double recurrence_R(const int& n, const int& t, const int& u, const int& v, const double& p, const Vector3d& Rpc) const;

    double overlap_xyz(const int& lx1, const int& ly1, const int& lz1, const double& a1, const Vector3d& X1, const int& lx2, const int& ly2, const int& lz2, const double& a2, const Vector3d& X2) const;
    double kinetic_xyz(const int& lx1, const int& ly1, const int& lz1, const double& a1, const Vector3d& X1, const int& lx2, const int& ly2, const int& lz2, const double& a2, const Vector3d& X2) const;
    double nucV_xyz(const int& lx1, const int& ly1, const int& lz1, const double& a1, const Vector3d& X1, const int& lx2, const int& ly2, const int& lz2, const double& a2, const Vector3d& X2) const;
    double eri_xyz(const int& lx1, const int& ly1, const int& lz1, const double& a1, const Vector3d& X1, const int& lx2, const int& ly2, const int& lz2, const double& a2, const Vector3d& X2, const int& lx3, const int& ly3, const int& lz3, const double& a3, const Vector3d& X3, const int& lx4, const int& ly4, const int& lz4, const double& a4, const Vector3d& X4) const;
public:
    int Nbasis_con, Nbasis_unc, Natom;
    MatrixXi atomList;
    VectorXi NbasisList_con, NbasisList_unc;
    Matrix<string,-1,1> basisList;
    Matrix<Vector3d,-1,1> coordList;
    Matrix<MatrixXd,-1,1> l_xyz_coeff;
    bool uncontracted;

    MOLINT(const MOL& molecule, const bool& unc_);
    ~MOLINT();

    void readBasis();
    double norm_factor(const int& lx, const int& ly, const int& lz, const double& a) const;
    double get_V_RR() const;
    MatrixXd get_h1e(const string& intName) const;
    VectorXd get_h2e(const string& intName) const;
};

#endif