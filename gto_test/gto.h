#ifndef GTO_H_
#define GTO_H_

#include<Eigen/Dense>
#include<string>
using namespace std;
using namespace Eigen;

const double PI = asin(1.0) * 2.0;

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
    int size, angular, atomN = 1;
    Matrix<gto_contracted, Dynamic, 1> gtos_c; 
    Vector3d center;

public:
    GTO();
    ~GTO();

    void readBasis();

    void normalization();
    MatrixXd get_h1e(const string& integralTYPE);
    MatrixXd get_h1e();
    MatrixXd get_overlap();
    MatrixXd get_nuc_attra();
    MatrixXd get_kinetic();


    double auxiliary(const int& l, const double& a);
    double overlap_single_gto(const gto_single& gto1, const gto_single& gto2);
    double nuc_attra_single_gto(const gto_single& gto1, const gto_single& gto2);
    double kinetic_single_gto(const gto_single& gto1, const gto_single& gto2);
    double h1e_single_gto(const gto_single& gto1, const gto_single& gto2);
    
};

#endif