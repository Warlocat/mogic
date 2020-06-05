#ifndef GTO_SPINOR_H_
#define GTO_SPINOR_H_

#include<Eigen/Dense>
#include<complex>
#include<string>
#include"gto.h"
using namespace std;
using namespace Eigen;


/*
    Class for single-atom (single-center) integrals in j-adopted spinor basis.
    Derived from class GTO.

    Variables:
        
*/
class GTO_SPINOR: public GTO
{
public:
    int size_gtoc_spinor, size_gtou_spinor;

    GTO_SPINOR(const string& atomName_, const string& basisSet_, const int& charge_ = 0, const int& spin_ = 1, const bool& uncontracted_ = false);
    ~GTO_SPINOR();

    /* return needed 1e and 2e integrals */
    MatrixXd get_h1e(const string& integralTYPE, const bool& uncontracted_ = false) const;
    MatrixXd get_h1e_spin_orbitals(const string& integralTYPE, const bool& uncontracted_ = false) const;
    MatrixXd get_h2e(const bool& uncontracted_ = false) const;
    
};

#endif