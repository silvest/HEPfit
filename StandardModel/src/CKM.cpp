/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "CKM.h"

CKM::CKM() : V(3, 3)
{}

void CKM::computeCKMwithWolfenstein(double Lambda_v, double A_v, double Rho_v, double Eta_v)
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
    
    computeCKMfromAngles();
}

void CKM::computeCKMfromAngles()
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

void CKM::computeCKMfromAngles(double s12_in, double s23_in, double s13_in, double delta_in)
{   
    s12 = s12_in;
    s13 = s13_in;
    s23 = s23_in;
    delta = delta_in;
    
    c12 = sqrt(1.-s12*s12);
    c13 = sqrt(1.-s13*s13);
    c23 = sqrt(1.-s23*s23);
    
    V.assign(0, 0, c12*c13);
    V.assign(0, 1, s12*c13);
    V.assign(0, 2, gslpp::complex(s13, -delta, true));

    V.assign(1, 0 , -s12 * c23 - gslpp::complex(c12 * s23*s13, delta, true));
    V.assign(1, 1, c12 * c23 - gslpp::complex(s12 * s23*s13, delta, true));
    V.assign(1, 2, s23*c13);

    V.assign(2, 0, s12 * s23 - gslpp::complex(c12 * c23*s13, delta, true));
    V.assign(2, 1, -c12 * s23 - gslpp::complex(s12 * c23*s13, delta, true));
    V.assign(2, 2, c23*c13);  
    
    // Wolfenstein to all orders
    Lambda = s12;
    A = s23 / Lambda / Lambda;
    gslpp::complex Rb = V(0, 0) * V(0, 2).conjugate() / (V(1, 0) * V(1, 2).conjugate());
    Rho = -Rb.real();
    Eta = -Rb.imag();
}

void CKM::computeCKM(double Vus_v, double Vcb_v, double Vub_v, double gamma_v, bool useVud)
{
    s13 = Vub_v;
    c13 = sqrt(1.-s13*s13);
    if (useVud) {
        c12 = Vus_v / c13;
        s12 = sqrt(1. - c12 * c12);
    }
    else {
        s12 = Vus_v / c13;
        c12 = sqrt(1. - s12 * s12);
    }
    
    s23 = Vcb_v / c13;
    c23 = sqrt(1. - s23 * s23);
    
    double a = c12 * s13 * s23 / s12 / c23;
    if ( fabs(gamma_v) < 1.e-10 )
        delta = 0.;
    else
        delta = 2. * atan((1. + sqrt(1. - (a * a - 1.) * pow(tan(gamma_v), 2.))*(cos(gamma_v) < 0. ? 1. : (-1.))) / (a - 1.) / tan(gamma_v));

    computeCKMfromAngles();
    
    // Wolfenstein to all orders
    Lambda = s12;
    A = s23 / Lambda / Lambda;
    gslpp::complex Rb = V(0, 0) * V(0, 2).conjugate() / (V(1, 0) * V(1, 2).conjugate());
    Rho = -Rb.real();
    Eta = -Rb.imag();
    return;
}


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

// Sides
double CKM::computeRt() const
{
    return (V(2, 0)*V(2, 2).conjugate()/(V(1, 0)*V(1, 2).conjugate())).abs();
}
double CKM::computeRts() const
{
    return (V(2, 1)*V(2, 2).conjugate()/(V(1, 1)*V(1, 2).conjugate())).abs();
}

double CKM::computeRb() const
{
    return (V(0, 0)*V(0, 2).conjugate()/(V(1, 0)*V(1, 2).conjugate())).abs();
}
