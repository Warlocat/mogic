#ifndef GTO_H_
#define GTO_H_

#include<Eigen/Dense>
#include<string>
using namespace std;
using namespace Eigen;

struct gto_single
{
    int lx, ly, lz;
    double a;
    Vector3d center;
};


class GTO
{
private:
    int size, angular;
    Matrix<gto_single, Dynamic, 1> gtos; 
    Vector3d center;
    VectorXd coeff;

public:
    GTO();
    ~GTO();

    void readBasis();

    void normalization();
    double overlap(const gto_single& gto1, const gto_single& gto2);
    
};

#endif