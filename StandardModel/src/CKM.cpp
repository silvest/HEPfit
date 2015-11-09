/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "CKM.h"

CKM::CKM()
{
}

CKM::CKM(const CKM& orig)
{
    Rho = orig.Rho;
    Eta = orig.Eta;
    Lambda = orig.Lambda;
    A = orig.A;

    s12 = orig.s12;
    s13 = orig.s13;
    s23 = orig.s23;
    delta = orig.delta;
    c12 = orig.c12;
    c23 = orig.c23;
    c13 = orig.c13;

    Vud = orig.Vud;
    Vcd = orig.Vcd;
    Vtd = orig.Vtd;
    Vus = orig.Vus;
    Vcs = orig.Vcs;
    Vts = orig.Vts;
    Vub = orig.Vub;
    Vcb = orig.Vcb;
    Vtb = orig.Vtb;

}

CKM::~CKM()
{
}

void CKM::setCKM(gslpp::matrix<gslpp::complex> & x)
{
    Vud = x(0, 0);
    Vus = x(0, 1);
    Vub = x(0, 2);
    Vcd = x(1, 0);
    Vcs = x(1, 1);
    Vcb = x(1, 2);
    Vtd = x(2, 0);
    Vts = x(2, 1);
    Vtb = x(2, 2);
}

void CKM::getCKM(gslpp::matrix<gslpp::complex> & x) const
{
    x.assign(0, 0, Vud);
    x.assign(0, 1, Vus);
    x.assign(0, 2, Vub);
    x.assign(1, 0, Vcd);
    x.assign(1, 1, Vcs);
    x.assign(1, 2, Vcb);
    x.assign(2, 0, Vtd);
    x.assign(2, 1, Vts);
    x.assign(2, 2, Vtb);
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

    c12 = sqrt(1. - pow(s12, 2.));
    c23 = sqrt(1. - pow(s23, 2.));
    c13 = sqrt(1. - pow(s13, 2.));

    Vud = gslpp::complex(c12*c13, 0.);
    Vus = gslpp::complex(s12*c13, 0.);
    Vub = gslpp::complex(s13, -delta, true);

    Vcd = -s12 * c23 - gslpp::complex(c12 * s23*s13, delta, true);
    Vcs = c12 * c23 - gslpp::complex(s12 * s23*s13, delta, true);
    Vcb = gslpp::complex(s23*c13, 0.);

    Vtd = s12 * s23 - gslpp::complex(c12 * c23*s13, delta, true);
    Vts = -c12 * s23 - gslpp::complex(s12 * c23*s13, delta, true);
    Vtb = gslpp::complex(c23*c13, 0.);

    return;
}

void CKM::setCKM(double Vud_v, double Vcb_v, double Vub_v, double gamma_v)
{
    s13 = Vub_v;
    c13 = sqrt(1. - s13 * s13);
    c12 = Vud_v / c13;
    s12 = sqrt(1. - c12 * c12);
    s23 = Vcb_v / c13;
    c23 = sqrt(1. - s23 * s23);
    double a = c12 * s13 * s23 / s12 / c23;
    delta = 2. * atan((1. + sqrt(1. - (a * a - 1.) * pow(tan(gamma_v), 2.))*(cos(gamma_v) < 0. ? 1. : (-1.))) / (a - 1.) / tan(gamma_v));

    Vud = gslpp::complex(c12*c13, 0.);
    Vus = gslpp::complex(s12*c13, 0.);
    Vub = gslpp::complex(s13, -delta, true);

    Vcd = -s12 * c23 - gslpp::complex(c12 * s23*s13, delta, true);
    Vcs = c12 * c23 - gslpp::complex(s12 * s23*s13, delta, true);
    Vcb = gslpp::complex(s23*c13, 0.);

    Vtd = s12 * s23 - gslpp::complex(c12 * c23*s13, delta, true);
    Vts = -c12 * s23 - gslpp::complex(s12 * c23*s13, delta, true);
    Vtb = gslpp::complex(c23*c13, 0.);

    // Wolfenstein to all orders
    Lambda = s12;
    A = s23 / Lambda / Lambda;
    gslpp::complex Rb = Vud * Vub.conjugate() / (Vcd * Vcb.conjugate());
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

double CKM::getRhoNB()
{
    return (s13 * cos(delta) / s12 / s23);
}

double CKM::getEtaNB()
{
    return (s13 * sin(delta) / s12 / s23);
}

// Gilman parameterization

double CKM::gets12()
{
    return getLambda();
}

double CKM::gets13()
{
    return s13;
}

double CKM::gets23()
{
    return s23;
}

double CKM::getc12()
{
    return c12;
}

double CKM::getc23()
{
    return c23;
}

double CKM::getc13()
{
    return c13;
}

double CKM::getdelta()
{
    return delta;
}

// J_CP

double CKM::getJcp()
{
    return Eta * pow(A * pow(Lambda, 3), 2);
}

//Absolute values of CKM elements

double CKM::getVud()
{
    return Vud.abs();
}

double CKM::getVus()
{
    return Vus.abs();
}

double CKM::getVub()
{
    return Vub.abs();
}

double CKM::getVcd()
{
    return Vcd.abs();
}

double CKM::getVcs()
{
    return Vcs.abs();
}

double CKM::getVcb()
{
    return Vcb.abs();
}

double CKM::getVtd()
{
    return Vtd.abs();
}

double CKM::getVts()
{
    return Vts.abs();
}

double CKM::getVtb()
{
    return Vtb.abs();
}

// Phases

double CKM::getArgVud()
{
    return Vud.arg();
}

double CKM::getArgVus()
{
    return Vus.arg();
}

double CKM::getArgVub()
{
    return Vub.arg();
}

double CKM::getArgVcd()
{
    return Vcd.arg();
}

double CKM::getArgVcs()
{
    return Vcs.arg();
}

double CKM::getArgVcb()
{
    return Vcb.arg();
}

double CKM::getArgVtd()
{
    return Vtd.arg();
}

double CKM::getArgVts()
{
    return Vts.arg();
}

double CKM::getArgVtb()
{
    return Vtb.arg();
}

// Angles

double CKM::computeBeta()
{
    return (-Vcd*Vcb.conjugate()/(Vtd*Vtb.conjugate())).arg();
}

double CKM::computeGamma() 
{
    return (-Vud*Vub.conjugate()/(Vcd*Vcb.conjugate())).arg();
}

double CKM::computeAlpha() 
{
    return (-Vtd*Vtb.conjugate()/(Vud*Vub.conjugate())).arg();
}

double CKM::computeBetas() 
{
    return (-Vts*Vtb.conjugate()/(Vcs*Vcb.conjugate())).arg();
}

// Lambda_q

gslpp::complex CKM::computelamt()
{
    return Vtd*Vts.conjugate();
}

gslpp::complex CKM::computelamc() 
{
    return Vcd*Vcs.conjugate();
}

gslpp::complex CKM::computelamu()
{
    return Vud*Vus.conjugate();
}


gslpp::complex CKM::computelamt_d() 
{
    return Vtd*Vtb.conjugate();
}

gslpp::complex CKM::computelamc_d() 
{
    return Vcd*Vcb.conjugate();
}

gslpp::complex CKM::computelamu_d() 
{
    return Vud*Vub.conjugate();
}


gslpp::complex CKM::computelamt_s() 
{
    return Vts*Vtb.conjugate();
}

gslpp::complex CKM::computelamc_s() 
{
    return Vcs*Vcb.conjugate();
}

gslpp::complex CKM::computelamu_s() 
{
    return Vus*Vub.conjugate();
}


//Complex values of CKM elements

gslpp::complex CKM::V_ud()
{
    return Vud;
}

gslpp::complex CKM::V_us()
{
    return Vus;
}

gslpp::complex CKM::V_ub()
{
    return Vub;
}

gslpp::complex CKM::V_cd()
{
    return Vcd;
}

gslpp::complex CKM::V_cs()
{
    return Vcs;
}

gslpp::complex CKM::V_cb()
{
    return Vcb;
}

gslpp::complex CKM::V_td()
{
    return Vtd;
}

gslpp::complex CKM::V_ts()
{
    return Vts;
}

gslpp::complex CKM::V_tb()
{
    return Vtb;
}

// Sides
double CKM::getRt()
{
    return (Vtd*Vtb.conjugate()/(Vcd*Vcb.conjugate())).abs();
}
double CKM::getRts() 
{
    return (Vts*Vtb.conjugate()/(Vcs*Vcb.conjugate())).abs();
}

double CKM::getRb()
{
    return (Vud*Vub.conjugate()/(Vcd*Vcb.conjugate())).abs();
}
