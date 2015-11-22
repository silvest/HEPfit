/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include <gsl/gsl_sf.h>
#include "EWSMTwoLoopQCD.h"

EWSMTwoLoopQCD::EWSMTwoLoopQCD(const EWSMcache& cache_i)
: cache(cache_i)
{
}


////////////////////////////////////////////////////////////////////////

double EWSMTwoLoopQCD::DeltaAlpha_l(const double s) const
{
    return (0.0);
}

double EWSMTwoLoopQCD::DeltaAlpha_t(const double s) const
{
    double xt = s / cache.getSM().getMtpole() / cache.getSM().getMtpole();
    double als;
    if (s == cache.getSM().getMz() * cache.getSM().getMz())
        als = cache.getSM().getAlsMz();
    else
        als = cache.Als(sqrt(s), FULLNNLO);
    double tmp = (5.062 + xt * 0.8315) * als / M_PI;
    tmp *= -4.0 / 45.0 * cache.getSM().getAle() / M_PI*xt;
    return tmp;
}

double EWSMTwoLoopQCD::DeltaRho(const double Mw_i) const
{
    double Mw = Mw_i;
    return ( 3.0 * cache.Xt_alpha(Mw) * cache.alsMt() / M_PI * deltaQCD_2());
}

double EWSMTwoLoopQCD::DeltaR_rem(const double Mw_i) const
{
    double Mw = Mw_i;
    return ( (2.0 * DeltaR_ud(Mw) + DeltaR_tb(Mw))
            + cache.getSM().cW2(Mw) / cache.getSM().sW2(Mw) * DeltaRho(Mw));
}

gslpp::complex EWSMTwoLoopQCD::deltaRho_rem_f(const Particle f, const double Mw_i) const
{
    if (f.is("TOP")) return ( gslpp::complex(0.0, 0.0, false));
    double Mw = Mw_i;
    return ( (2.0 * DeltaRho_ud(Mw) + DeltaRho_tb(Mw)) - DeltaRho(Mw));
}

gslpp::complex EWSMTwoLoopQCD::deltaKappa_rem_f(const Particle f, const double Mw_i) const
{
    if (f.is("TOP")) return ( gslpp::complex(0.0, 0.0, false));
    double Mw = Mw_i;
    return ( (2.0 * DeltaKappa_ud(Mw) + DeltaKappa_tb(Mw))
            - cache.getSM().cW2(Mw) / cache.getSM().sW2(Mw) * DeltaRho(Mw));
}


////////////////////////////////////////////////////////////////////////

double EWSMTwoLoopQCD::deltaQCD_2() const
{
    return ( -2.0 / 3.0 * (1.0 + 2.0 * cache.getZeta2()));
}

double EWSMTwoLoopQCD::F1(const double x, const double Mw_i) const
{
    if (x < 0.0 || x >= 1.0) throw std::runtime_error("x is out of range in EWSMTwoLoopQCD::F1");

    /* Zeta functions */
    double zeta_2 = cache.getZeta2();
    double zeta_3 = cache.getZeta3();

    if (x == 0.0) return (23.0 / 16.0 - zeta_2 / 2.0 - 3.0 / 2.0 * zeta_3);

    /* Dilogarithm and Trilogarithm */
    double Li2_x, Li3_x, Li3_mx_1mx;
    double Mw = Mw_i;
    double Mt = cache.getSM().getMtpole();
    if (x == Mw * Mw / Mt / Mt) {
        Li2_x = cache.Li2_MW2toMTOP2(Mw);
        Li3_x = cache.Li3_MW2toMTOP2(Mw);
        Li3_mx_1mx = cache.Li3_for_F1(Mw);
    } else {
        Li2_x = cache.getPolyLog().Li2(x).real(); // x <= 1.0
        Li3_x = cache.getPolyLog().Li3(x);
        Li3_mx_1mx = cache.getPolyLog().Li3(-x / (1.0 - x));
    }

    double b = log(1.0 - x);

    double F1;
    F1 = (x - 3.0 / 2.0 + 1.0 / 2.0 / x / x)
            *(-Li3_x - Li3_mx_1mx + b / 3.0 * (2.0 * Li2_x - zeta_2) + b * b * b / 6.0)
            + 1.0 / 3.0 * (x + 1.0 / 2.0 - 1.0 / 2.0 / x) * Li2_x
            + b * b / 6.0 * (x - 3.0 / 4.0 - 3.0 / 2.0 / x + 5.0 / 4.0 / x / x)
            - b / 4.0 * (x - 5.0 / 2.0 + 2.0 / 3.0 / x + 5.0 / 6.0 / x / x)
            + zeta_3 * (x - 3.0 / 2.0) + zeta_2 / 3.0 * (x - 7.0 / 4.0 - 1.0 / 2.0 / x)
            + 13.0 / 12.0 - 5.0 / 24.0 / x;
    return F1;
}

double EWSMTwoLoopQCD::V1(const double r) const
{
    if (r < 0.0 || r >= 1.0) throw std::runtime_error("r is out of range in EWSMTwoLoopQCD::V1");

    /* Zeta functions */
    double zeta_3 = cache.getZeta3();

    if (r == 0.0) return (0.0);

    double Mz = cache.getSM().getMz();
    double Mt = cache.getSM().getMtpole();

    /* Logarithms etc */
    double Phi, gamma, h;
    if (r == Mz * Mz / 4.0 / Mt / Mt) {
        Phi = cache.Phi_QCD2();
        gamma = cache.gamma_QCD2();
        h = cache.h_QCD2();
    } else {
        Phi = asin(sqrt(r));
        gamma = log(2.0 * sqrt(r));
        h = log(2.0 * sqrt(1.0 - r));
    }

    /* Clausen functions */
    double Cl3_2Phi, Cl3_4Phi, Cl2_2Phi, Cl2_4Phi;
    if (r == Mz * Mz / 4.0 / Mt / Mt) {
        Cl3_2Phi = cache.Cl3_2Phi();
        Cl3_4Phi = cache.Cl3_4Phi();
        Cl2_2Phi = cache.Cl2_2Phi();
        Cl2_4Phi = cache.Cl2_4Phi();
    } else {
        ClausenFunctions* myClausen;
        myClausen = new ClausenFunctions();
        Cl3_2Phi = myClausen->Cl3(2.0 * Phi);
        Cl3_4Phi = myClausen->Cl3(4.0 * Phi);
        Cl2_2Phi = myClausen->Cl2(2.0 * Phi);
        Cl2_4Phi = myClausen->Cl2(4.0 * Phi);
        delete myClausen;
    }

    double V1;
    V1 = 4.0 * (r - 1.0 / 4.0 / r)
            *(2.0 * Cl3_2Phi - Cl3_4Phi + 8.0 / 3.0 * Phi * (Cl2_2Phi - Cl2_4Phi)
            - 4.0 / 3.0 * Phi * Phi * (gamma + 2.0 * h))
            + sqrt(1.0 / r - 1.0)*(8.0 / 3.0 * (r + 1.0 / 2.0)
            *(-Cl2_2Phi + Cl2_4Phi + 2.0 * Phi * (gamma + 2.0 * h))
            - 2.0 * Phi * (r + 3.0 / 2.0))
            + 8.0 * Phi * Phi * (r - 1.0 / 6.0 - 7.0 / 48.0 / r) + 13.0 / 6.0 + zeta_3 / r;
    return V1;
}

double EWSMTwoLoopQCD::A1(const double r) const
{
    if (r < 0.0 || r >= 1.0) throw std::runtime_error("r is out of range in EWSMTwoLoopQCD::A1");

    /* Zeta functions */
    double zeta_2 = cache.getZeta2();
    double zeta_3 = cache.getZeta3();

    if (r == 0.0) return (3.0 * (7.0 / 4.0 - zeta_2 - 2.0 * zeta_3));

    double Mz = cache.getSM().getMz();
    double Mt = cache.getSM().getMtpole();

    /* Logarithms etc */
    double Phi, gamma, h;
    if (r == Mz * Mz / 4.0 / Mt / Mt) {
        Phi = cache.Phi_QCD2();
        gamma = cache.gamma_QCD2();
        h = cache.h_QCD2();
    } else {
        Phi = asin(sqrt(r));
        gamma = log(2.0 * sqrt(r));
        h = log(2.0 * sqrt(1.0 - r));
    }

    /* Clausen functions */
    double Cl3_2Phi, Cl3_4Phi, Cl2_2Phi, Cl2_4Phi;
    if (r == Mz * Mz / 4.0 / Mt / Mt) {
        Cl3_2Phi = cache.Cl3_2Phi();
        Cl3_4Phi = cache.Cl3_4Phi();
        Cl2_2Phi = cache.Cl2_2Phi();
        Cl2_4Phi = cache.Cl2_4Phi();
    } else {
        ClausenFunctions* myClausen;
        myClausen = new ClausenFunctions();
        Cl3_2Phi = myClausen->Cl3(2.0 * Phi);
        Cl3_4Phi = myClausen->Cl3(4.0 * Phi);
        Cl2_2Phi = myClausen->Cl2(2.0 * Phi);
        Cl2_4Phi = myClausen->Cl2(4.0 * Phi);
        delete myClausen;
    }

    double A1;
    A1 = 4.0 * (r - 3.0 / 2.0 + 1.0 / 2.0 / r)*(2.0 * Cl3_2Phi - Cl3_4Phi
            + 8.0 / 3.0 * Phi * (Cl2_2Phi - Cl2_4Phi) - 4.0 / 3.0 * Phi * Phi * (gamma + 2.0 * h))
            + sqrt(1.0 / r - 1.0)*(8.0 / 3.0 * (r - 1.0)
            *(-Cl2_2Phi + Cl2_4Phi + 2.0 * Phi * (gamma + 2.0 * h))
            - 2.0 * Phi * (r - 3.0 + 1.0 / 4.0 / r))
            + 8.0 * Phi * Phi * (r - 11.0 / 12.0 + 5.0 / 48.0 / r + 1.0 / 32.0 / r / r)
            - 3.0 * zeta_2 + 13.0 / 6.0 + (-2.0 * zeta_3 + 1.0 / 4.0) / r;
    return A1;
}

double EWSMTwoLoopQCD::V1prime(const double r) const
{
    if (r < 0.0 || r >= 1.0) throw std::runtime_error("r is out of range in EWSMTwoLoopQCD::V1prime");

    /* Zeta functions */
    double zeta_3 = cache.getZeta3();

    if (r == 0.0) return (4.0 * zeta_3 - 5.0 / 6.0);

    double Mz = cache.getSM().getMz();
    double Mt = cache.getSM().getMtpole();

    /* Logarithms etc */
    double Phi, gamma, h;
    if (r == Mz * Mz / 4.0 / Mt / Mt) {
        Phi = cache.Phi_QCD2();
        gamma = cache.gamma_QCD2();
        h = cache.h_QCD2();
    } else {
        Phi = asin(sqrt(r));
        gamma = log(2.0 * sqrt(r));
        h = log(2.0 * sqrt(1.0 - r));
    }

    /* Clausen functions */
    double Cl3_2Phi, Cl3_4Phi, Cl2_2Phi, Cl2_4Phi;
    if (r == Mz * Mz / 4.0 / Mt / Mt) {
        Cl3_2Phi = cache.Cl3_2Phi();
        Cl3_4Phi = cache.Cl3_4Phi();
        Cl2_2Phi = cache.Cl2_2Phi();
        Cl2_4Phi = cache.Cl2_4Phi();
    } else {
        ClausenFunctions* myClausen;
        myClausen = new ClausenFunctions();
        Cl3_2Phi = myClausen->Cl3(2.0 * Phi);
        Cl3_4Phi = myClausen->Cl3(4.0 * Phi);
        Cl2_2Phi = myClausen->Cl2(2.0 * Phi);
        Cl2_4Phi = myClausen->Cl2(4.0 * Phi);
        delete myClausen;
    }

    /* for Bprime and Dprime below */
    gsl_complex OneMinusE2Iphi = gsl_complex_rect(1.0 - cos(2.0 * Phi), -sin(2.0 * Phi));
    gsl_complex OneMinusE4Iphi = gsl_complex_rect(1.0 - cos(4.0 * Phi), -sin(4.0 * Phi));
    double log_real;
    if (r == Mz * Mz / 4.0 / Mt / Mt) {
        log_real = cache.logV1primeAndA1prime();
    } else {
        log_real = GSL_REAL(gsl_complex_log(OneMinusE2Iphi))
                - 2.0 * GSL_REAL(gsl_complex_log(OneMinusE4Iphi));
    }

    double PhiPrime = 1.0 / 2.0 / sqrt(r * (1.0 - r));
    double Phi2Prime = Phi / sqrt(r * (1.0 - r));
    double gammaPrime = 1.0 / 2.0 / r;
    double hPrime = -1.0 / 2.0 / (1.0 - r);

    // V1(r) = 4.0*A*B + C*D + E
    double A = r - 1.0 / 4.0 / r;
    double Aprime = 1.0 + 1.0 / 4.0 / r / r;
    double B = 2.0 * Cl3_2Phi - Cl3_4Phi + 8.0 / 3.0 * Phi * (Cl2_2Phi - Cl2_4Phi)
            - 4.0 / 3.0 * Phi * Phi * (gamma + 2.0 * h);
    double Bprime = -2.0 / sqrt(r * (1.0 - r))*(Cl2_2Phi - Cl2_4Phi)
            + 8.0 / 3.0 * PhiPrime * (Cl2_2Phi - Cl2_4Phi)
            + 8.0 / 3.0 * Phi * (-log_real / sqrt(r * (1.0 - r)))
            - 4.0 / 3.0 * Phi2Prime * (gamma + 2.0 * h)
            - 4.0 / 3.0 * Phi * Phi * (gammaPrime + 2.0 * hPrime);
    double C = sqrt(1.0 / r - 1.0);
    double Cprime = -1.0 / 2.0 / r / sqrt(r * (1.0 - r));
    double D = 8.0 / 3.0 * (r + 1.0 / 2.0)
            *(-Cl2_2Phi + Cl2_4Phi + 2.0 * Phi * (gamma + 2.0 * h))
            - 2.0 * Phi * (r + 3.0 / 2.0);
    double Dprime = 8.0 / 3.0 * (-Cl2_2Phi + Cl2_4Phi + 2.0 * Phi * (gamma + 2.0 * h))
            + 8.0 / 3.0 * (r + 1.0 / 2.0)
            *(log_real / sqrt(r * (1.0 - r))
            + 2.0 * PhiPrime * (gamma + 2.0 * h)
            + 2.0 * Phi * (gammaPrime + 2.0 * hPrime))
            - 2.0 * PhiPrime * (r + 3.0 / 2.0) - 2.0 * Phi;
    double Eprime = 8.0 * Phi2Prime * (r - 1.0 / 6.0 - 7.0 / 48.0 / r)
            + 8.0 * Phi * Phi * (1.0 + 7.0 / 48.0 / r / r) - zeta_3 / r / r;
    return (4.0 * Aprime * B + 4.0 * A * Bprime + Cprime * D + C * Dprime + Eprime);


    /* TEST: Exact - Expansion */
    //std::cout << "V1: " 
    //          << (4.0*Aprime*B + 4.0*A*Bprime + Cprime*D + C*Dprime + Eprime)
    //              - (4.0*zeta_3 - 5.0/6.0 + 656.0/81.0*r)
    //          << std::endl;

    /* Expansion for r << 1 */
    //return (4.0*zeta_3 - 5.0/6.0 + 656.0/81.0*r);
}

double EWSMTwoLoopQCD::A1prime(const double r) const
{
    if (r < 0.0 || r >= 1.0) throw std::runtime_error("r is out of range in EWSMTwoLoopQCD::A1prime");

    /* Zeta functions */
    double zeta_2 = cache.getZeta2();
    double zeta_3 = cache.getZeta3();

    if (r == 0.0) return (3.0 * (7.0 / 4.0 - zeta_2 - 2.0 * zeta_3));

    double Mz = cache.getSM().getMz();
    double Mt = cache.getSM().getMtpole();

    /* Logarithms etc */
    double Phi, gamma, h;
    if (r == Mz * Mz / 4.0 / Mt / Mt) {
        Phi = cache.Phi_QCD2();
        gamma = cache.gamma_QCD2();
        h = cache.h_QCD2();
    } else {
        Phi = asin(sqrt(r));
        gamma = log(2.0 * sqrt(r));
        h = log(2.0 * sqrt(1.0 - r));
    }

    /* Clausen functions */
    double Cl3_2Phi, Cl3_4Phi, Cl2_2Phi, Cl2_4Phi;
    if (r == Mz * Mz / 4.0 / Mt / Mt) {
        Cl3_2Phi = cache.Cl3_2Phi();
        Cl3_4Phi = cache.Cl3_4Phi();
        Cl2_2Phi = cache.Cl2_2Phi();
        Cl2_4Phi = cache.Cl2_4Phi();
    } else {
        ClausenFunctions* myClausen;
        myClausen = new ClausenFunctions();
        Cl3_2Phi = myClausen->Cl3(2.0 * Phi);
        Cl3_4Phi = myClausen->Cl3(4.0 * Phi);
        Cl2_2Phi = myClausen->Cl2(2.0 * Phi);
        Cl2_4Phi = myClausen->Cl2(4.0 * Phi);
        delete myClausen;
    }

    /* for Bprime and Dprime below */
    gsl_complex OneMinusE2Iphi = gsl_complex_rect(1.0 - cos(2.0 * Phi), -sin(2.0 * Phi));
    gsl_complex OneMinusE4Iphi = gsl_complex_rect(1.0 - cos(4.0 * Phi), -sin(4.0 * Phi));
    double log_real;
    if (r == Mz * Mz / 4.0 / Mt / Mt) {
        log_real = cache.logV1primeAndA1prime();
    } else {
        log_real = GSL_REAL(gsl_complex_log(OneMinusE2Iphi))
                - 2.0 * GSL_REAL(gsl_complex_log(OneMinusE4Iphi));
    }

    double PhiPrime = 1.0 / 2.0 / sqrt(r * (1.0 - r));
    double Phi2Prime = Phi / sqrt(r * (1.0 - r));
    double gammaPrime = 1.0 / 2.0 / r;
    double hPrime = -1.0 / 2.0 / (1.0 - r);

    // A1(r) = 4.0*A*B + C*D + E
    double A = r - 3.0 / 2.0 + 1.0 / 2.0 / r;
    double Aprime = 1.0 - 1.0 / 2.0 / r / r;
    double B = 2.0 * Cl3_2Phi - Cl3_4Phi + 8.0 / 3.0 * Phi * (Cl2_2Phi - Cl2_4Phi)
            - 4.0 / 3.0 * Phi * Phi * (gamma + 2.0 * h);
    double Bprime = -2.0 / sqrt(r * (1.0 - r))*(Cl2_2Phi - Cl2_4Phi)
            + 8.0 / 3.0 * PhiPrime * (Cl2_2Phi - Cl2_4Phi)
            + 8.0 / 3.0 * Phi * (-log_real / sqrt(r * (1.0 - r)))
            - 4.0 / 3.0 * Phi2Prime * (gamma + 2.0 * h)
            - 4.0 / 3.0 * Phi * Phi * (gammaPrime + 2.0 * hPrime);
    double C = sqrt(1.0 / r - 1.0);
    double Cprime = -1.0 / 2.0 / r / sqrt(r * (1.0 - r));
    double D = 8.0 / 3.0 * (r - 1.0)*(-Cl2_2Phi + Cl2_4Phi + 2.0 * Phi * (gamma + 2.0 * h))
            - 2.0 * Phi * (r - 3.0 + 1.0 / 4.0 / r);
    double Dprime = 8.0 / 3.0 * (-Cl2_2Phi + Cl2_4Phi + 2.0 * Phi * (gamma + 2.0 * h))
            + 8.0 / 3.0 * (r - 1.0)*(log_real / sqrt(r * (1.0 - r))
            + 2.0 * PhiPrime * (gamma + 2.0 * h)
            + 2.0 * Phi * (gammaPrime + 2.0 * hPrime))
            - 2.0 * PhiPrime * (r - 3.0 + 1.0 / 4.0 / r)
            - 2.0 * Phi * (1.0 - 1.0 / 4.0 / r / r);
    double Eprime = 8.0 * Phi2Prime * (r - 11.0 / 12.0 + 5.0 / 48.0 / r + 1.0 / 32.0 / r / r)
            + 8.0 * Phi * Phi * (1.0 - 5.0 / 48.0 / r / r - 1.0 / 16.0 / r / r / r)
            - (-2.0 * zeta_3 + 1.0 / 4.0) / r / r;
    return (4.0 * Aprime * B + 4.0 * A * Bprime + Cprime * D + C * Dprime + Eprime);

    /* TEST: Exact - Expansion */
    //std::cout << "A1:" 
    //          <<  (4.0*Aprime*B + 4.0*A*Bprime + Cprime*D + C*Dprime + Eprime)
    //             - (4.0*zeta_3 - 49.0/18.0 + 1378.0/405.0*r )
    //          << std::endl;

    /* Expansion for r << 1 */
    //return (4.0*zeta_3 - 49.0/18.0 + 1378.0/405.0*r);    
}

double EWSMTwoLoopQCD::DeltaR_ud(const double Mw_i) const
{
    double Mw = Mw_i;
    double sW2 = cache.getSM().sW2(Mw);
    double cW2 = cache.getSM().cW2(Mw);

    /* Logarithm */
    double log_cW2 = cache.log_cW2(Mw);

    double DeltaR;
    DeltaR = -log_cW2;
    DeltaR *= (cW2 - sW2) / 4.0 / sW2 / sW2;
    DeltaR *= cache.getSM().getAle() * cache.getSM().getAlsMz() / M_PI / M_PI;
    return DeltaR;
}

double EWSMTwoLoopQCD::DeltaR_tb(const double Mw_i) const
{
    double Mw = Mw_i;
    double sW2 = cache.getSM().sW2(Mw);
    double cW2 = cache.getSM().cW2(Mw);
    double Mz = cache.getSM().getMz();
    double Mt = cache.getSM().getMtpole();
    double wt = Mt * Mt / Mw / Mw;
    double zt = Mt * Mt / Mz / Mz;
    double rZ4t = Mz * Mz / 4.0 / Mt / Mt;
    double xWt = Mw * Mw / Mt / Mt;

    double vt = cache.v_f(cache.getSM().getQuarks(QCD::TOP), Mw);
    double at = cache.a_f(cache.getSM().getQuarks(QCD::TOP));
    double vb = cache.v_f(cache.getSM().getQuarks(QCD::BOTTOM), Mw);
    double ab = cache.a_f(cache.getSM().getQuarks(QCD::BOTTOM));

    /* Zeta functions */
    double zeta_2 = cache.getZeta2();

    /* Logarithm */
    double log_zt = -2.0 * cache.logMZtoMTOP();

    double DeltaR;
    DeltaR = pow(cache.Q_f(cache.getSM().getQuarks(QCD::TOP)), 2.0) * V1prime(0.0)
            + cW2 / sW2 / sW2 * wt / 4.0 * (zeta_2 + 1.0 / 2.0)
            - zt / sW2 / sW2 * (vt * vt * V1(rZ4t) + at * at * (A1(rZ4t) - A1(0.0)))
            + (cW2 - sW2) / sW2 / sW2 * wt * (F1(xWt, Mw) - F1(0.0, Mw))
            - vb * ab / 2.0 / sW2 / sW2*log_zt;
    DeltaR *= cache.getSM().getAle() * cache.alsMt() / M_PI / M_PI;
    return DeltaR;
}

double EWSMTwoLoopQCD::DeltaRho_ud(const double Mw_i) const
{
    double Mw = Mw_i;
    double sW2 = cache.getSM().sW2(Mw);
    double cW2 = cache.getSM().cW2(Mw);

    double DeltaRho;
    DeltaRho = pow(cache.v_f(cache.getSM().getQuarks(QCD::UP), Mw), 2.0)
            + pow(cache.v_f(cache.getSM().getQuarks(QCD::DOWN), Mw), 2.0)
            + pow(cache.a_f(cache.getSM().getQuarks(QCD::UP)), 2.0)
            + pow(cache.a_f(cache.getSM().getQuarks(QCD::TOP)), 2.0);
    DeltaRho /= 4.0 * sW2*cW2;
    DeltaRho *= cache.getSM().getAle() * cache.getSM().getAlsMz() / M_PI / M_PI;
    return DeltaRho;
}

double EWSMTwoLoopQCD::DeltaRho_tb(const double Mw_i) const
{
    double Mw = Mw_i;
    double Mz = cache.getSM().getMz();
    double sW2 = cache.getSM().sW2(Mw);
    double cW2 = cache.getSM().cW2(Mw);
    double Mt = cache.getSM().getMtpole();
    double zt = Mt * Mt / Mz / Mz;
    double rZ4t = Mz * Mz / 4.0 / Mt / Mt;

    double vt = cache.v_f(cache.getSM().getQuarks(QCD::TOP), Mw);
    double at = cache.a_f(cache.getSM().getQuarks(QCD::TOP));
    double vb = cache.v_f(cache.getSM().getQuarks(QCD::BOTTOM), Mw);
    double ab = cache.a_f(cache.getSM().getQuarks(QCD::BOTTOM));

    double DeltaRho;
    DeltaRho = -(vt * vt * V1prime(rZ4t) + at * at * A1prime(rZ4t))
            + 4.0 * zt * (vt * vt * V1(rZ4t) + at * at * A1(rZ4t))
            + vb * vb + ab * ab - 4.0 * zt * F1(0.0, Mw);
    DeltaRho /= 4.0 * sW2*cW2;
    DeltaRho *= cache.getSM().getAle() * cache.alsMt() / M_PI / M_PI;
    return DeltaRho;
}

gslpp::complex EWSMTwoLoopQCD::DeltaKappa_ud(const double Mw_i) const
{
    double Mw = Mw_i;
    double sW2 = cache.getSM().sW2(Mw);
    double cW2 = cache.getSM().cW2(Mw);

    /* Logarithm */
    double log_cW2 = cache.log_cW2(Mw);

    gslpp::complex DeltaKappa(0.0, 0.0, false);
    DeltaKappa = cW2 / 4.0 / sW2 / sW2 * log_cW2
            + M_PI / 4.0 / sW2 * (1.0 - 20.0 / 9.0 * sW2)*(gslpp::complex::i());
    DeltaKappa *= cache.getSM().getAle() * cache.getSM().getAlsMz() / M_PI / M_PI;
    return DeltaKappa;
}

gslpp::complex EWSMTwoLoopQCD::DeltaKappa_tb(const double Mw_i) const
{
    double Mw = Mw_i;
    double Mz = cache.getSM().getMz();
    double sW2 = cache.getSM().sW2(Mw);
    double cW2 = cache.getSM().cW2(Mw);
    double Mt = cache.getSM().getMtpole();
    double wt = Mt * Mt / Mw / Mw;
    double zt = Mt * Mt / Mz / Mz;
    double rZ4t = Mz * Mz / 4.0 / Mt / Mt;
    double xWt = Mw * Mw / Mt / Mt;

    double vt = cache.v_f(cache.getSM().getQuarks(QCD::TOP), Mw);
    double at = cache.a_f(cache.getSM().getQuarks(QCD::TOP));
    double Qt = cache.Q_f(cache.getSM().getQuarks(QCD::TOP));
    double vb = cache.v_f(cache.getSM().getQuarks(QCD::BOTTOM), Mw);
    double ab = cache.a_f(cache.getSM().getQuarks(QCD::BOTTOM));
    double Qb = cache.Q_f(cache.getSM().getQuarks(QCD::BOTTOM));

    /* Logarithm */
    double log_zt = -2.0 * cache.logMZtoMTOP();

    gslpp::complex DeltaKappa(0.0, 0.0, false);
    DeltaKappa = 4.0 * cW2 * wt * (vt * vt * V1(rZ4t) + at * at * A1(rZ4t) - F1(xWt, Mw))
            + 4.0 * sW2 * (fabs(Qt) - 4.0 * sW2 * Qt * Qt) * zt * V1(rZ4t)
            + (vb * vb + ab * ab + sW2 * (fabs(Qb) - 4.0 * sW2 * Qb * Qb)) * log_zt;
    DeltaKappa += M_PI * sW2 * (1.0 / 3.0 - 4.0 / 9.0 * sW2)*(gslpp::complex::i());
    DeltaKappa /= 4.0 * sW2*sW2;
    DeltaKappa *= cache.getSM().getAle() * cache.alsMt() / M_PI / M_PI;
    return DeltaKappa;
}




