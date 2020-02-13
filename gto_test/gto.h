#ifndef GTO_H_
#define GTO_H_

#include<Eigen/Dense>
#include<string>
using namespace std;
using namespace Eigen;


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
    int angular, atomN = 1;
    Matrix<gto_contracted, Dynamic, 1> gtos_c; 
    Vector3d center;

public:
    int size;

    GTO();
    ~GTO();

    void readBasis();

    void normalization();
    MatrixXd get_h1e(const string& integralTYPE);
    Matrix<MatrixXd, -1, -1> get_h2e();


    double auxiliary_1e(const int& l, const double& a);
    double auxiliary_2e(const int& l1, const int& l2, const double& a1, const double& a2);
    double int1e_single_gto(const gto_single& gto1, const gto_single& gto2, const string& integralTYPE);
    double int2e_single_gto(const gto_single& gto1, const gto_single& gto2, const gto_single& gto3, const gto_single& gto4);
    
};

#endif