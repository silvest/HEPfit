/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "CKM.h"

CKM::CKM() : V(3, 3)
{
}


void CKM::setWolfenstein(double Lambda_v, double A_v, double Rho_v, double Eta_v)
{
    Rho = Rho_v;
    Eta = Eta_v;
    Lambda = Lambda_v;
    A = A_v;

    gslpp::complex num(Rho, Eta);
    num = num * sqrt(1. - pow(A, 2.) * pow(Lambda, 4.));
    gslpp::complex den = sqrt(1. - pow(Lambda, 2.)) * gslpp::complex(1. - pow(A, 2.) * pow(Lambda, 4.) * Rho, -pow(A, 2.) * pow(Lambda, 4.) * Eta);
    gslpp::complex ratio = num / den;

    double rho_nb = ratio.real();
    double eta_nb = ratio.imag();

    s12 = Lambda;
    s23 = A * pow(Lambda, 2.);
    s13 = (gslpp::complex(A * pow(Lambda, 3.) * rho_nb, -A * pow(Lambda, 3.) * eta_nb)).abs();
    delta = -(gslpp::complex(A * pow(Lambda, 3.) * rho_nb, -A * pow(Lambda, 3.) * eta_nb)).arg();

    c12 = sqrt(1. - s12 * s12);
    c13 = sqrt(1. - s13 * s13);
    c23 = sqrt(1. - s23 * s23);
    
    setCKMfromAngles();
}

void CKM::setCKMfromAngles()
{   
    V.assign(0, 0, c12*c13);
    V.assign(0, 1, s12*c13);
    V.assign(0, 2, gslpp::complex(s13, -delta, true));

    V.assign(1, 0 , -s12 * c23 - gslpp::complex(c12 * s23*s13, delta, true));
    V.assign(1, 1, c12 * c23 - gslpp::complex(s12 * s23*s13, delta, true));
    V.assign(1, 2, s23*c13);

    V.assign(2, 0, s12 * s23 - gslpp::complex(c12 * c23*s13, delta, true));
    V.assign(2, 1, -c12 * s23 - gslpp::complex(s12 * c23*s13, delta, true));
    V.assign(2, 2, c23*c13);  
}

void CKM::setCKM(double Vus_v, double Vcb_v, double Vub_v, double gamma_v)
{
    s13 = Vub_v;
    c13 = sqrt(1.-s13*s13);
    s12 = Vus_v / c13;
    s23 = Vcb_v / c13;

    c12 = sqrt(1. - s12 * s12);
    c23 = sqrt(1. - s23 * s23);
    
    double a = c12 * s13 * s23 / s12 / c23;
    if ( fabs(gamma_v) < 1.e-10 )
        delta = 0.;
    else
        delta = 2. * atan((1. + sqrt(1. - (a * a - 1.) * pow(tan(gamma_v), 2.))*(cos(gamma_v) < 0. ? 1. : (-1.))) / (a - 1.) / tan(gamma_v));

    setCKMfromAngles();
    
    // Wolfenstein to all orders
    Lambda = s12;
    A = s23 / Lambda / Lambda;
    gslpp::complex Rb = V(0, 0) * V(0, 2).conjugate() / (V(1, 0) * V(1, 2).conjugate());
    Rho = -Rb.real();
    Eta = -Rb.imag();

    //  std::cout << Lambda << " " << A << " " << Rho << " " << Eta << std::endl;
    return;
}

// Wolfenstein parameters

double CKM::getRho() const
{
    return Rho;
}

double CKM::getEta() const
{
    return Eta;
}

double CKM::getLambda() const
{
    return Lambda;
}

double CKM::getA() const
{
    return A;
}

double CKM::getRhoNB() const
{
    return (s13 * cos(delta) / s12 / s23);
}

double CKM::getEtaNB() const
{
    return (s13 * sin(delta) / s12 / s23);
}

// Gilman parameterization

double CKM::gets12() const
{
    return s12;
}

double CKM::gets13() const
{
    return s13;
}

double CKM::gets23() const
{
    return s23;
}

double CKM::getc12() const
{
    return c12;
}

double CKM::getc23() const
{
    return c23;
}

double CKM::getc13() const
{
    return c13;
}

double CKM::getdelta() const
{
    return delta;
}

// J_CP

double CKM::getJcp() const
{
    return Eta * pow(A * pow(Lambda, 3), 2);
}

//Absolute values of CKM elements

double CKM::getVud() const
{
    return V(0, 0).abs();
}

double CKM::getVus() const
{
    return V(0, 1).abs();
}

double CKM::getVub() const
{
    return V(0, 2).abs();
}

double CKM::getVcd() const
{
    return V(1, 0).abs();
}

double CKM::getVcs() const
{
    return V(1, 1).abs();
}

double CKM::getVcb() const
{
    return V(1, 2).abs();
}

double CKM::getVtd() const
{
    return V(2, 0).abs();
}

double CKM::getVts() const
{
    return V(2, 1).abs();
}

double CKM::getVtb() const
{
    return V(2, 2).abs();
}

// Phases

double CKM::getArgVud() const
{
    return V(0, 0).arg();
}

double CKM::getArgVus() const
{
    return V(0, 1).arg();
}

double CKM::getArgVub() const
{
    return V(0, 2).arg();
}

double CKM::getArgVcd() const
{
    return V(1, 0).arg();
}

double CKM::getArgVcs() const
{
    return V(1, 1).arg();
}

double CKM::getArgVcb() const
{
    return V(1, 2).arg();
}

double CKM::getArgVtd() const
{
    return V(2, 0).arg();
}

double CKM::getArgVts() const
{
    return V(2, 1).arg();
}

double CKM::getArgVtb() const
{
    return V(2, 2).arg();
}

// Angles

double CKM::computeBeta() const
{
    return (-V(1, 0)*V(1, 2).conjugate()/(V(2, 0)*V(2, 2).conjugate())).arg();
}

double CKM::computeGamma() const
{
    return (-V(0, 0)*V(0, 2).conjugate()/(V(1, 0)*V(1, 2).conjugate())).arg();
}

double CKM::computeAlpha() const 
{
    return (-V(2, 0)*V(2, 2).conjugate()/(V(0, 0)*V(0, 2).conjugate())).arg();
}

double CKM::computeBetas() const 
{
    return (-V(2, 1)*V(2, 2).conjugate()/(V(1, 1)*V(1, 2).conjugate())).arg();
}

// Lambda_q

gslpp::complex CKM::computelamt() const
{
    return V(2, 0)*V(2, 1).conjugate();
}

gslpp::complex CKM::computelamc() const 
{
    return V(1, 0)*V(1, 1).conjugate();
}

gslpp::complex CKM::computelamu() const
{
    return V(0, 0)*V(0, 1).conjugate();
}


gslpp::complex CKM::computelamt_d() const 
{
    return V(2, 0)*V(2, 2).conjugate();
}

gslpp::complex CKM::computelamc_d() const 
{
    return V(1, 0)*V(1, 2).conjugate();
}

gslpp::complex CKM::computelamu_d() const 
{
    return V(0, 0)*V(0, 2).conjugate();
}


gslpp::complex CKM::computelamt_s() const 
{
    return V(2, 1)*V(2, 2).conjugate();
}

gslpp::complex CKM::computelamc_s() const 
{
    return V(1, 1)*V(1, 2).conjugate();
}

gslpp::complex CKM::computelamu_s() const
{
    return V(0, 1)*V(0, 2).conjugate();
}


//Complex values of CKM elements

gslpp::complex CKM::V_ud() const
{
    return V(0, 0);
}

gslpp::complex CKM::V_us() const
{
    return V(0, 1);
}

gslpp::complex CKM::V_ub() const
{
    return V(0, 2);
}

gslpp::complex CKM::V_cd() const
{
    return V(1, 0);
}

gslpp::complex CKM::V_cs() const
{
    return V(1, 1);
}

gslpp::complex CKM::V_cb() const
{
    return V(1, 2);
}

gslpp::complex CKM::V_td() const
{
    return V(2, 0);
}

gslpp::complex CKM::V_ts() const
{
    return V(2, 1);
}

gslpp::complex CKM::V_tb() const
{
    return V(2, 2);
}

// Sides
double CKM::getRt() const
{
    return (V(2, 0)*V(2, 2).conjugate()/(V(1, 0)*V(1, 2).conjugate())).abs();
}
double CKM::getRts() const
{
    return (V(2, 1)*V(2, 2).conjugate()/(V(1, 1)*V(1, 2).conjugate())).abs();
}

double CKM::getRb() const
{
    return (V(0, 0)*V(0, 2).conjugate()/(V(1, 0)*V(1, 2).conjugate())).abs();
}
