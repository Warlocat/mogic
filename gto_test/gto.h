#ifndef GTO_H_
#define GTO_H_

#include<Eigen/Dense>
#include<complex>
#include<string>
using namespace std;
using namespace Eigen;

double factorial(const int& n);
double double_factorial(const int& n);
double wigner_3j(const int& l1, const int& l2, const int& l3, const int& m1, const int& m2, const int& m3);
double wigner_3j_zeroM(const int& l1, const int& l2, const int& l3);
complex<double> U_SH_trans(const int& mu, const int& mm);

struct gto_single
{
    int l, m;
    double a;
};
struct gto_contracted
{
    Matrix<gto_single, Dynamic, 1> gto_list;
    VectorXd coeff;
};


class GTO
{
private:
    int atomN = 1;
    Matrix<gto_contracted, Dynamic, 1> gtos_c; 
    Vector3d center;

public:
    int size;

    GTO();
    ~GTO();

    void readBasis(const string& atomName, const string& filename);

    void normalization();
    MatrixXd get_h1e(const string& integralTYPE);
    Matrix<MatrixXd, -1, -1> get_h2e();


    double auxiliary_1e(const int& l, const double& a);
    double auxiliary_2e_0_r(const int& l1, const int& l2, const double& a1, const double& a2);
    double auxiliary_2e_r_inf(const int& l1, const int& l2, const double& a1, const double& a2);
    double int1e_single_gto(const gto_single& gto1, const gto_single& gto2, const string& integralTYPE);
    double int2e_single_gto(const gto_single& gto1, const gto_single& gto2, const gto_single& gto3, const gto_single& gto4);
    
};

#endif