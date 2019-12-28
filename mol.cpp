#include<Eigen/Dense>
#include<string>
#include<iostream>
#include<iomanip>
#include<fstream>
#include"mol.h"
using namespace std;
using namespace Eigen;

/* Construction and destruction functions */
MOL::MOL()
{
}

MOL::~MOL()
{
}

/* Inner product and cross product */
double MOL::innerP(const VectorXd& v1, const VectorXd& v2)
{
    return v1.transpose()*v2;
}
Vector3d MOL::crossP(const Vector3d& v1, const Vector3d& v2)
{
    Vector3d cross;
    cross(0) = v1(1) * v2(2) - v2(1) * v1(2);
    cross(1) = v1(2) * v2(0) - v2(2) * v1(0);
    cross(2) = v1(0) * v2(1) - v2(0) * v1(1);
    return cross;
}

/* For mass initialization */
double MOL::massInit(const int& an)
{
    switch (an)
    {
        case 1: return 1.008;
        case 2: return 4.0026;
        case 3: return 6.94;
        case 4: return 9.0122;
        case 5: return 10.81;
        case 6: return 12.011;
        case 7: return 14.007;
        case 8: return 15.999;
        case 9: return 18.998;
        case 10: return 20.180;
        case 11: return 22.99;
        case 12: return 24.305;
        case 13: return 26.982;
        case 14: return 28.085;
        case 15: return 30.974;
        case 16: return 32.06;
        case 17: return 35.45;
        case 18: return 39.948;
        case 19: return 39.098;
        case 20: return 40.078;
        case 21: return 44.956;
        case 22: return 47.867;
        case 23: return 50.942;
        case 24: return 51.996;
        case 25: return 54.938;
        case 26: return 55.845;
        case 27: return 58.933;
        case 28: return 58.693;
        case 29: return 63.546;
        case 30: return 65.38;
        case 31: return 69.723;
        case 32: return 72.63;
        case 33: return 74.922;
        case 34: return 78.971;
        case 35: return 79.904;
        case 36: return 83,798;
        case 37: return 85.468;
        case 38: return 87.62;
        case 39: return 88.906;
        case 40: return 91.224;
        case 41: return 92.906;
        case 42: return 95.95;
        case 43: return 98;
        case 44: return 101.07;
        case 45: return 102.91;
        case 46: return 106.42;
        case 47: return 107.87;
        case 48: return 112.41;
        case 49: return 114.82;
        case 50: return 118.71;
        case 51: return 121.76;
        case 52: return 127.6;
        case 53: return 126.9;
        case 54: return 131.29;
        case 55: return 132.91;
        case 56: return 137.33;
        // case 57: return ;
        // case 58: return ;
        default:
            cout << "Input atomic number is larger than 56 or smaller than 1" << endl;
            exit(99);
    }
}

/* Read and Write XYZ */
void MOL::readxyz(const string& filename)
{
    ifstream ifs;
    ifs.open(filename);
        ifs >> Natom;
        Coord_cart.resize(Natom, 3);
        atomicNumber.resize(Natom);
        mass.resize(Natom);
        // ifs.getline(flagc, 99);
        for(int ii = 0; ii < Natom; ii++)
        {
            ifs >> atomicNumber(ii);
            mass(ii) = massInit(atomicNumber(ii)) / ATOMIC_MASS_UNIT;
            for(int jj = 0; jj < 3; jj++) ifs >> Coord_cart(ii, jj);
        }
    ifs.close();
}
void MOL::writexyz(const string& xyzfile)
{
    ofstream ofs;
    ofs.open(xyzfile);
        ofs << "\t" << Natom << endl << endl;
        for(int ii = 0; ii < Natom; ii++)
        {
            ofs << fixed << setprecision(10) << atomicNumber(ii) << "\t" << Coord_cart(ii, 0) << "\t" << Coord_cart(ii, 1) << "\t" << Coord_cart(ii, 2) << endl;
        }
    ofs.close();
}

/* transfer Coord_cart to center of mass coordinates */
void MOL::cart2COM()
{
    double totalmass = 0.0;
    Vector3d com;
    com = com * 0.0;
    for(int ii = 0; ii < Natom; ii++)
    {
        totalmass += mass(ii);
        for(int jj = 0; jj < 3; jj++)
        {
            com(jj) += mass(ii) * Coord_cart(ii, jj);
        }
    }
    com = com / totalmass;
    for(int ii = 0; ii < Natom; ii++)
    for(int jj = 0; jj < 3; jj++)
    {
        Coord_cart(ii, jj) = Coord_cart(ii, jj) - com(jj);
    }
}
    
/* Calculate bond length r_ij between i and j atom */
double MOL::length(const int& i, const int& j)
{
    Vector3d tmp;
    for(int ii = 0; ii < 3; ii++) tmp(ii) = Coord_cart(i, ii) - Coord_cart(j, ii);
    return tmp.norm();
}
/* Calculate bond angle phi_ijk between i, j and k atom */
double MOL::angle(const int& i, const int& j, const int& k)
{
    double tmp;
    Vector3d tmp1, tmp2;
    for(int ii = 0; ii < 3; ii++)
    {
        tmp1(ii) = Coord_cart(i, ii) - Coord_cart(j, ii);
        tmp2(ii) = Coord_cart(k, ii) - Coord_cart(j, ii);
    }
    tmp = tmp1.transpose() * tmp2;
    tmp = tmp / tmp1.norm() / tmp2.norm();
    if(tmp > 1.0) tmp = 1.0;
    else if(tmp < -1.0) tmp = -1.0;
    return acos(tmp);
}
/* Calculate sine of out-of-plane angle theta_ijkl between i, j, k and l atom */
double MOL::angle_ofp(const int& i, const int& j, const int& k, const int& l)
{
    Vector3d tmpi, tmpj, tmpl;
    for(int ii = 0; ii < 3; ii++)
    {
        tmpi(ii) = Coord_cart(i, ii) - Coord_cart(k, ii);
        tmpj(ii) = Coord_cart(j, ii) - Coord_cart(k, ii);
        tmpl(ii) = Coord_cart(l, ii) - Coord_cart(k, ii);
    }
    tmpi = tmpi / tmpi.norm();
    tmpj = tmpj / tmpj.norm();
    tmpl = tmpl / tmpl.norm();
    /* What returns here is sin theta */
    return innerP(crossP(tmpj, tmpl), tmpi) / sin(angle(j,k,l));
}
/* Calculate torsion angle tau_ijkl between i, j, k and l atom */
double MOL::torsion(const int& i, const int& j, const int& k, const int& l)
{
    double tmp;
    Vector3d tmp1, tmp2, tmp3;
    for(int ii = 0; ii < 3; ii++)
    {
        tmp1(ii) = Coord_cart(j, ii) - Coord_cart(i, ii);
        tmp2(ii) = Coord_cart(k, ii) - Coord_cart(j, ii);
        tmp3(ii) = Coord_cart(l, ii) - Coord_cart(k, ii);
    }
    tmp1 = tmp1 / tmp1.norm();
    tmp2 = tmp2 / tmp2.norm();
    tmp3 = tmp3 / tmp3.norm();
    tmp = innerP(crossP(tmp1, tmp2),crossP(tmp2, tmp3)) / sin(angle(i,j,k)) / sin(angle(j,k,l));
}

/* Calculate Moments of Inertia */
Matrix3d MOL::rotationI()
{
    Matrix3d I;
    I = I * 0.0;
    for(int ii = 0; ii < Natom; ii++)
    {
        I(0, 0) += mass(ii) * (pow(Coord_cart(ii, 1), 2) + pow(Coord_cart(ii, 2), 2));
        I(1, 1) += mass(ii) * (pow(Coord_cart(ii, 0), 2) + pow(Coord_cart(ii, 2), 2));
        I(2, 2) += mass(ii) * (pow(Coord_cart(ii, 0), 2) + pow(Coord_cart(ii, 1), 2));
        I(0, 1) += mass(ii) * Coord_cart(ii, 0) * Coord_cart(ii, 1);
        I(0, 2) += mass(ii) * Coord_cart(ii, 0) * Coord_cart(ii, 2);
        I(1, 2) += mass(ii) * Coord_cart(ii, 1) * Coord_cart(ii, 2);
    }
    I(1, 0) = I(0, 1);
    I(2, 0) = I(0, 2);
    I(2, 1) = I(1, 2);
    return I;
}