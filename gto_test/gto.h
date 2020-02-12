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
    Vector3d center;
};


class GTO
{
private:
    int size, angular, atomN = 1;
    Matrix<gto_single, Dynamic, 1> gtos; 
    Vector3d center;
    VectorXd coeff;

public:
    GTO();
    ~GTO();

    void readBasis();

    void normalization();

    double auxiliary(const int& l, const double& a);
    double overlap(const gto_single& gto1, const gto_single& gto2);
    double nuc_attraction(const gto_single& gto1, const gto_single& gto2);
    double kinetic(const gto_single& gto1, const gto_single& gto2);
    double h_1e(const gto_single& gto1, const gto_single& gto2);
    
};

#endif