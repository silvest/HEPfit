/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include <TMath.h>
#include "EWSMTwoLoopEW.h"

/* include O(alpha^2 Mt2/Mz2) in addition to O(alpha^2 M_t^4/M_Z^4) */
#define EW_SUBLEADING_ALPHA2

EWSMTwoLoopEW::EWSMTwoLoopEW(const EWSMcache& cache_i)
: cache(cache_i), myOneLoopEW(cache_i)
{
}


////////////////////////////////////////////////////////////////////////

double EWSMTwoLoopEW::DeltaAlpha_l(const double s) const
{
    double xl[3] = {s / cache.mf(cache.getSM().getLeptons(StandardModel::ELECTRON)) / cache.mf(cache.getSM().getLeptons(StandardModel::ELECTRON)),
        s / cache.mf(cache.getSM().getLeptons(StandardModel::MU)) / cache.mf(cache.getSM().getLeptons(StandardModel::MU)),
        s / cache.mf(cache.getSM().getLeptons(StandardModel::TAU)) / cache.mf(cache.getSM().getLeptons(StandardModel::TAU))};
    double log_l[3];
    if (s == cache.getSM().getMz() * cache.getSM().getMz()) {
        log_l[0] = 2.0 * cache.logMZtoME();
        log_l[1] = 2.0 * cache.logMZtoMMU();
        log_l[2] = 2.0 * cache.logMZtoMTAU();
    } else {
        log_l[0] = log(xl[0]);
        log_l[1] = log(xl[1]);
        log_l[2] = log(xl[2]);
    }

    double twoLoop[3];
    for (int i = 0; i < 3; i++) {
        twoLoop[i] = -5.0 / 24.0 + cache.getZeta3() + log_l[i] / 4.0
                + 3.0 / xl[i] * log_l[i];
    }

    return ( pow(cache.getSM().getAle() / M_PI, 2.0)
            *(twoLoop[0] + twoLoop[1] + twoLoop[2]));
}

double EWSMTwoLoopEW::DeltaAlpha_t(const double s) const
{
    return (0.0);
}

double EWSMTwoLoopEW::DeltaRho(const double Mw_i) const
{
    double Mz = cache.getSM().getMz();
    double Mw = Mw_i;
    double sW2 = cache.getSM().sW2(Mw);
    double cW2 = cache.getSM().cW2(Mw);

    double DeltaRho = 0.0;

#ifndef EW_SUBLEADING_ALPHA2
    /* O(\alpha^2 Mt^4/Mz^4) */
    DeltaRho += 3.0 * rho_2();
    DeltaRho *= pow(cache.Xt_alpha(Mw), 2.0);
#else
    /* O(\alpha^2 Mt^4/Mz^4 + \alpha^2 Mt^2/Mz^2) */
    double zt = Mz * Mz / cache.getSM().getMtpole() / cache.getSM().getMtpole();
    DeltaRho += 3.0 * pow(cache.Xt_alpha(Mw), 2.0)
            *(DeltaRho2(Mw) + 4.0 * zt * cW2 * DeltaRho2Add(Mw));
#endif

    /* add O(alpha^2) contribution from the Z-gamma mixing */
    DeltaRho += -pow(cache.getSM().getAle() / 4.0 / M_PI, 2.0) * cW2 / sW2
            * pow(myOneLoopEW.PibarZgamma_fer(Mz, Mz*Mz, Mw).real(), 2.0);

    return DeltaRho;
}

double EWSMTwoLoopEW::DeltaR_rem(const double Mw_i) const
{
    double DeltaRrem = 0.0;

#ifdef EW_SUBLEADING_ALPHA2
    /* O(\alpha^2 Mt^4/Mz^4 + \alpha^2 Mt^2/Mz^2) */
    double Mw = Mw_i;
    double sW2 = cache.getSM().sW2(Mw);
    DeltaRrem += 3.0 * pow(cache.getSM().getAle() * cache.getSM().getMtpole() / 4.0 / M_PI / sW2 / Mw, 2.0)
            *(DeltaRw2(Mw) + sW2 * deltaEoverE2() + f2Add(Mw) / 4.0);
#endif    

    return DeltaRrem;
}

gslpp::complex EWSMTwoLoopEW::deltaRho_rem_f(const Particle f, const double Mw_i) const
{
    if (f.is("TOP")) return ( gslpp::complex(0.0, 0.0, false));

    gslpp::complex dRho = gslpp::complex(0.0, 0.0, false);

#ifdef EW_SUBLEADING_ALPHA2
    /* O(\alpha^2 Mt^4/Mz^4 + \alpha^2 Mt^2/Mz^2) */
    double Mz = cache.getSM().getMz();
    double Mw = Mw_i;
    double cW2 = cache.getSM().cW2(Mw);
    double zt = Mz * Mz / cache.getSM().getMtpole() / cache.getSM().getMtpole();
    dRho += 3.0 * pow(cache.Xt_alpha(Mw), 2.0)
            *(16.0 * zt * cW2 * DeltaEta2(Mw) + 4.0 * zt * cW2 * DeltaEta2Add_f(f, Mw));
#endif 

    return dRho;
}

gslpp::complex EWSMTwoLoopEW::deltaKappa_rem_f(const Particle f, const double Mw_i) const
{
    if (f.is("TOP")) return ( gslpp::complex(0.0, 0.0, false));

    gslpp::complex dKappa = gslpp::complex(0.0, 0.0, false);

#ifdef EW_SUBLEADING_ALPHA2
    /* O(\alpha^2 Mt^4/Mz^4 + \alpha^2 Mt^2/Mz^2) */
    double Mz = cache.getSM().getMz();
    double Mw = Mw_i;
    double cW2 = cache.getSM().cW2(Mw);
    double zt = Mz * Mz / cache.getSM().getMtpole() / cache.getSM().getMtpole();
    dKappa += 3.0 * pow(cache.Xt_alpha(Mw), 2.0)
            *(16.0 * zt * cW2 * DeltaKappa2(Mw) + 4.0 * zt * cW2 * DeltaKappa2Add_f(f, Mw));
#endif 

    return dKappa;
}


////////////////////////////////////////////////////////////////////////   `

double EWSMTwoLoopEW::rho_2() const
{
    double a = cache.getSM().getMHl() * cache.getSM().getMHl() / cache.getSM().getMtpole() / cache.getSM().getMtpole();
    if (a <= 0.0) throw std::runtime_error("a is out of range in EWSMTwoLoopEW::rho_2");
    double g_a = g(a);
    double f_a_0 = f0(a); // f(a,0)
    double f_a_1 = f1(a); // f(a,1)
    double log_a = -2.0 * cache.logMTOPtoMH();
    return ( 25.0 - 4.0 * a + 0.5 * (a * a - 12.0 * a - 12.0) * log_a
            + (a - 2.0) / 2.0 / a * M_PI * M_PI + 0.5 * (a - 4.0) * sqrt(a) * g_a
            - 3.0 / a * (a - 1.0)*(a - 1.0)*(a - 2.0) * f_a_0
            + 3.0 * (a * a - 6.0 * a + 10.0) * f_a_1);
}

double EWSMTwoLoopEW::tau_2() const
{
    double a = cache.getSM().getMHl() * cache.getSM().getMHl() / cache.getSM().getMtpole() / cache.getSM().getMtpole();
    if (a <= 0.0) throw std::runtime_error("a is out of range in EWSMTwoLoopEW::tau_2");
    double g_a = g(a);
    double f_a_0 = f0(a); // f(a,0)
    double f_a_1 = f1(a); // f(a,1)
    double log_a = -2.0 * cache.logMTOPtoMH();
    return ( 9.0 - 13.0 / 4.0 * a - 2.0 * a * a - a / 4.0 * (19.0 + 6.0 * a) * log_a
            - a * a / 4.0 * (7.0 - 6.0 * a) * log_a * log_a
            - (1.0 / 4.0 + 7.0 / 2.0 * a * a - 3.0 * a * a * a) * M_PI * M_PI / 6.0
            + (a / 2.0 - 2.0) * sqrt(a) * g_a
            + (a - 1.0)*(a - 1.0)*(4.0 * a - 7.0 / 4.0) * f_a_0
            - (a * a * a - 33.0 / 4.0 * a * a + 18.0 * a - 7.0) * f_a_1);
}

double EWSMTwoLoopEW::g(const double a) const
{
    if (a >= 0.0 && a <= 4.0) {
        double phi = 2.0 * asin(sqrt(a / 4.0));
        return ( sqrt(4.0 - a)*(M_PI - phi));
    } else if (a > 4.0) {
        double y = 4.0 / a;
        double xi = (sqrt(1.0 - y) - 1.0) / (sqrt(1.0 - y) + 1.0);
        return ( sqrt(a - 4.0) * log(-xi));
    } else
        throw std::runtime_error("Out of range in EWSMTwoLoopEW::g()");
}

double EWSMTwoLoopEW::f0(const double a) const
{
    if (a >= 0.0)
        return ( cache.getPolyLog().Li2(1.0 - a).real()); // 1-a<1
    else
        throw std::runtime_error("Out of range in EWSMTwoLoopEW::f0()");
}

double EWSMTwoLoopEW::f1(const double a) const
{
    if (a >= 0.0 && a <= 4.0) {
        double y = 4.0 / a;
        double phi = 2.0 * asin(sqrt(a / 4.0));
        return ( -2.0 / sqrt(y - 1.0) * cache.getClausen().Cl2(phi));
    } else if (a > 4.0) {
        double y = 4.0 / a; // 0<y<1
        double xi = (sqrt(1.0 - y) - 1.0) / (sqrt(1.0 - y) + 1.0); // -1<xi<0
        return ( -1.0 / sqrt(1.0 - y)*(cache.getPolyLog().Li2(xi).real()
                - cache.getPolyLog().Li2(1.0 / xi).real()));
    } else
        throw std::runtime_error("Out of range in EWSMTwoLoopEW::f1()");
}

double EWSMTwoLoopEW::DeltaRho2(const double Mw_i) const
{
    double Mt = cache.getSM().getMtpole();
    double mh = cache.getSM().getMHl();
    double ht = mh * mh / Mt / Mt;

    double rho2;
    if (mh < 3.8 * Mt) {
        rho2 = -15.642 + 0.036382 * Mt + pow(ht, 1.0 / 4.0)*(2.301 - 0.01343 * Mt)
                + sqrt(ht)*(0.01809 * Mt - 9.953) + ht * (5.687 - 0.01568 * Mt)
                + pow(ht, 3.0 / 2.0)*(0.005369 * Mt - 1.647)
                + ht * ht * (0.1852 - 0.000646 * Mt);
    } else {
        double Mz = cache.getSM().getMz();
        double Mw = Mw_i;
        double Mz2 = Mz*Mz, Mw2 = Mw*Mw, Mt2 = Mt*Mt;
        double cW2 = cache.getSM().cW2(Mw), sW2 = cache.getSM().sW2(Mw);
        double cW4 = cW2*cW2, cW6 = cW4*cW2;
        double zt = Mz * Mz / Mt / Mt;
        double ht2 = ht*ht, ht3 = ht2*ht, ht4 = ht3*ht, ht5 = ht4*ht, ht6 = ht5*ht;
        double mu = Mt; // renormalization scale

        double B0_Mt2_Mz2_Mw2_Mw2 = cache.getPV().B0(Mt2, Mz2, Mw2, Mw2).real();
        double B0_Mt2_Mw2_Mw2_Mz2 = cache.getPV().B0(Mt2, Mw2, Mw2, Mz2).real();
        double Li21mht = cache.getPolyLog().Li2(1.0 - ht).real();

        return ( 25.0 - 4.0 * ht + (1.0 / 2.0 - 1.0 / ht) * M_PI * M_PI
                + (ht - 4.0) * sqrt(ht) * g(ht) / 2.0
                + (-6.0 - 6.0 * ht + ht2 / 2.0) * log(ht)
                + (6.0 / ht - 15.0 + 12.0 * ht - 3.0 * ht2) * Li21mht
                + 3.0 / 2.0 * (-10.0 + 6.0 * ht - ht2) * phi(ht / 4.0)
                + zt * (1.0 / (54.0 * cW2 * (ht - 4.0) * ht)
                *(-1776.0 * cW4
                + (72.0 - 6250.0 * cW2 - 3056.0 * cW4 + 3696.0 * cW6) * ht
                + (-18.0 + 1283.0 * cW2 + 1371.0 * cW4 - 1436.0 * cW6) * ht2
                + (68.0 * cW2 - 124.0 * cW4 + 128.0 * cW6) * ht3)
                + (6.0 * cW2 * ht - 37.0 * cW2 - 119.0 * ht2 + 56.0 * cW2 * ht2)
                * M_PI * M_PI / 27.0 / ht2
                + (32.0 * cW4 / 3.0 - 2.0 / 3.0 - 12.0 * cW2) * B0_Mt2_Mz2_Mw2_Mw2
                + (20.0 / 3.0 + 1.0 / 3.0 / cW2 - 8.0 * cW2) * B0_Mt2_Mw2_Mw2_Mz2
                + (17.0 - 58.0 * cW2 + 32.0 * cW4)*(4.0 - ht) * sqrt(ht) * g(ht) / 27.0
                - 40.0 * sW2 * (4.0 - ht) * Lambda(ht) / 3.0 / ht
                + 2.0 * cW2 * (37.0 - 6.0 * ht - 12.0 * ht2 - 22.0 * ht3 + 9.0 * ht4)
                * Li21mht / 9.0 / ht2
                - (1.0 + 14.0 * cW2 + 16.0 * cW4) * log(cW2) / 3.0
                + (11520.0 - 15072.0 * cW2
                - (7170.0 - 8928.0 * cW2 - 768.0 * cW4) * ht
                + (3411.0 - 7062.0 * cW2 + 3264.0 * cW4) * ht2
                - (1259.0 - 3547.0 * cW2 + 2144.0 * cW4) * ht3
                + (238.0 - 758.0 * cW2 + 448.0 * cW4) * ht4
                - (17.0 - 58.0 * cW2 + 32.0 * cW4) * ht5)
                * log(ht) / 27.0 / (ht - 4.0) / (ht - 4.0) / ht
                + 8.0 / 9.0 * (4.0 - 26.0 * cW2 - 5.0 * cW4) * log(Mt * Mt / mu / mu)
                + (3.0 + 5.0 * cW2 - 26.0 * cW4 - 48.0 * cW6) * log(zt) / 9.0 / cW2
                + (3840.0 * sW2 - (4310.0 - 4224.0 * cW2 - 256.0 * cW4) * ht
                + (1706.0 - 1312.0 * cW2 - 320.0 * cW4) * ht2
                - (315.0 + 476.0 * cW2 - 64.0 * cW4) * ht3
                + (24.0 + 454.0 * cW2) * ht4 - 112.0 * cW2 * ht5
                + 9.0 * cW2 * ht6)
                * phi(ht / 4.0) / 9.0 / (ht - 4.0) / (ht - 4.0) / ht2));
    }
    return rho2;
}

double EWSMTwoLoopEW::DeltaRho2Add(const double Mw_i) const
{
    double Mz = cache.getSM().getMz(), Mz2 = Mz*Mz;
    double Mw = Mw_i, Mw2 = Mw*Mw;
    double cW2 = cache.getSM().cW2(Mw);
    double cW4 = cW2*cW2, cW6 = cW4*cW2;
    double Mt = cache.getSM().getMtpole(), Mt2 = Mt*Mt;
    double zt = Mz2 / Mt2;
    double mu = Mt; // renormalization scale

    double B0_Mt2_Mz2_Mw2_Mw2 = cache.getPV().B0(Mt2, Mz2, Mw2, Mw2).real();
    double B0_Mt2_Mw2_0_Mw2 = cache.getPV().B0(Mt2, Mw2, 0.0, Mw2).real();
    double B0_Mt2_Mw2_Mw2_Mz2 = cache.getPV().B0(Mt2, Mw2, Mw2, Mz2).real();

    double dRho2add = 542.0 / 27.0 - 2.0 / 3.0 / cW2 - 800.0 * cW2 / 27.0
            + 1.0 / 3.0 * (1.0 + 26.0 * cW2 + 24.0 * cW4) * B0_Mt2_Mz2_Mw2_Mw2
            + 4.0 * cW2 * B0_Mt2_Mw2_0_Mw2
            - (11.0 / 3.0 + 1.0 / 3.0 / cW2 + 4.0 * cW2) * B0_Mt2_Mw2_Mw2_Mz2
            - (2.0 / 3.0 + 4.0 * cW2 / 3.0 - 8.0 * cW4) * log(cW2)
            + (1.0 / cW2 - 38.0 / 3.0 + 34.0 * cW2 / 3.0) * log(Mt2 / mu / mu)
            + 2.0 * (3.0 - 62.0 * cW2 + 74.0 * cW4 + 36.0 * cW6) * log(zt) / 9.0 / cW2;
    return dRho2add;
}

double EWSMTwoLoopEW::DeltaRw2(const double Mw_i) const
{
    double Mz = cache.getSM().getMz(), Mz2 = Mz*Mz;
    double Mw = Mw_i, Mw2 = Mw*Mw;
    double cW2 = cache.getSM().cW2(Mw), sW2 = cache.getSM().sW2(Mw);
    double cW4 = cW2*cW2, cW6 = cW4*cW2;
    double Mt = cache.getSM().getMtpole(), Mt2 = Mt*Mt;
    double mh = cache.getSM().getMHl(), mh2 = mh*mh;
    double zt = Mz2 / Mt2;
    double zt2 = zt*zt;
    double ht = mh2 / Mt2;
    double ht2 = ht*ht, ht3 = ht2*ht, ht4 = ht3*ht, ht5 = ht4*ht;
    double mu = Mt; // renormalization scale

    double dRw2;
    if (mh < 0.3 * Mt) {
        double B0_Mt2_Mw2_mh2_Mw2 = cache.getPV().B0(Mt2, Mw2, mh2, Mw2).real();
        double B0_Mt2_Mw2_Mw2_Mz2 = cache.getPV().B0(Mt2, Mw2, Mw2, Mz2).real();
        dRw2 = -13.0 / 144.0 - 1.0 / 48.0 / cW4 - 41.0 / 96.0 / cW2 + 61.0 * cW2 / 72.0
                + (7.0 - 16.0 * cW2) / 27.0 * M_PI * sqrt(ht) - M_PI * M_PI / 36.0
                - 5.0 * ht2 / 144.0 / cW4 / zt2 + 35.0 * ht / 288.0 / cW2 / zt
                + 5.0 / 12.0 * (1.0 + ht2 / 12.0 / cW4 / zt2 - ht / 3.0 / cW2 / zt) * B0_Mt2_Mw2_mh2_Mw2
                + (1.0 + 20.0 * cW2 - 24.0 * cW4) / 48.0 / cW4 * B0_Mt2_Mw2_Mw2_Mz2
                - (5.0 * sW2 * ht2 + 3.0 * ht * zt + 48.0 * cW2 * ht * zt - 60.0 * cW4 * ht * zt
                - 3.0 * cW2 * zt2 - 8.0 * cW4 * zt2 + 20.0 * cW6 * zt2) * log(cW2)
                / 144.0 / cW2 / sW2 / zt / (ht - cW2 * zt)
                + 5.0 * ht * (ht2 - 4.0 * cW2 * ht * zt + 12.0 * cW4 * zt2) * log(ht)
                / 144.0 / cW4 / zt2 / (ht - cW2 * zt)
                + (17.0 / 36.0 - 13.0 * cW2 / 18.0) * log(Mt2 / mu / mu)
                - (5.0 * cW2 * ht2 - 3.0 * ht * zt - 60.0 * cW2 * ht * zt + 60.0 * cW4 * ht * zt
                + (3.0 * cW2 + 60.0 * cW4 - 20.0 * cW6) * zt2) * log(zt)
                / 144.0 / cW4 / zt / (ht - cW2 * zt);
    } else {
        double B0_Mt2_Mw2_Mw2_Mz2 = cache.getPV().B0(Mt2, Mw2, Mw2, Mz2).real();
        dRw2 = -121.0 / 288.0 - 1.0 / 48.0 / cW4 - 41.0 / 96.0 / cW2 + 77.0 * cW2 / 12.0
                + 19.0 / 72.0 / ht + (41.0 / 216.0 - 4.0 * cW2 / 27.0) * ht
                - (19.0 + 21.0 * ht) * M_PI * M_PI / 432.0 / ht2
                - (1.0 / 2.0 - 1.0 / 48.0 / cW4 - 5.0 / 12.0 / cW2) * B0_Mt2_Mw2_Mw2_Mz2
                + (16.0 * cW2 - 7.0) / 216.0 * (ht - 4.0) * sqrt(ht) * g(ht)
                - (1.0 / 12.0 - 1.0 / 3.0 / ht) * Lambda(ht)
                + (19.0 + 21.0 * ht - 12.0 * ht2 - 31.0 * ht3 + 9.0 * ht4) / 72.0 / ht2
                * cache.getPolyLog().Li2(1.0 - ht).real()
                - (1.0 + 21.0 * cW2 - 25.0 * cW4) * log(cW2) / 48.0 / cW2 / sW2
                + (17.0 / 36.0 - 13.0 * cW2 / 18.0) * log(Mt2 / mu / mu)
                + (1.0 + 20.0 * cW2 - 25.0 * cW4) * log(zt) / 48.0 / cW4
                + (372.0 + (96.0 * cW2 - 213.0) * ht + (432.0 * cW2 - 318.0) * ht2
                + (97.0 - 160.0 * cW2) * ht3 - (7.0 - 16.0 * cW2) * ht4)
                / 216.0 / (ht - 4.0) / ht * log(ht)
                + (96.0 - (384.0 - 64.0 * cW2) * ht - (2.0 + 64.0 * cW2) * ht2 + 231.0 * ht3
                - 85.0 * ht4 + 9.0 * ht5) / 144.0 / (ht - 4.0) / ht2 * phi(ht / 4.0);
    }
    return dRw2;
}

double EWSMTwoLoopEW::deltaEoverE2() const
{
    double Mt = cache.getSM().getMtpole(), Mt2 = Mt*Mt;
    double mh = cache.getSM().getMHl(), mh2 = mh*mh;
    double ht = mh2 / Mt2;
    double ht2 = ht*ht, ht3 = ht2*ht;
    double mu = Mt; // renormalization scale

    double dEoE2;
    if (mh < 0.3 * Mt) {
        dEoE2 = 61.0 / 72.0 - 16.0 * sqrt(ht) * M_PI / 27.0 - 13.0 / 18.0 * log(Mt2 / mu / mu);
    } else {
        dEoE2 = (231.0 - 32.0 * ht) / 216.0 - 2.0 / 27.0 * (4.0 - ht) * sqrt(ht) * g(ht)
                + 2.0 * (6.0 + 27.0 * ht - 10.0 * ht2 + ht3) / 27.0 / (ht - 4.0) * log(ht)
                - 13.0 / 18.0 * log(Mt2 / mu / mu)
                - 4.0 * (ht - 1.0) / 9.0 / (ht - 4.0) / ht * phi(ht / 4.0);
    }
    return dEoE2;
}

double EWSMTwoLoopEW::f2Add(const double Mw_i) const
{
    double Mz = cache.getSM().getMz(), Mz2 = Mz*Mz;
    double Mw = Mw_i, Mw2 = Mw*Mw;
    double cW2 = cache.getSM().cW2(Mw), sW2 = cache.getSM().sW2(Mw);
    double Mt = cache.getSM().getMtpole(), Mt2 = Mt*Mt;
    double zt = Mz2 / Mt2;

    double B0_Mt2_Mw2_0_Mw2 = cache.getPV().B0(Mt2, Mw2, 0.0, Mw2).real();
    double B0_Mt2_Mw2_Mw2_Mz2 = cache.getPV().B0(Mt2, Mw2, Mw2, Mz2).real();

    double f2a = 10.0 / 3.0 + 1.0 / 3.0 / cW2 + 4.0 * cW2 * B0_Mt2_Mw2_0_Mw2
            - (11.0 / 3.0 + 1.0 / 3.0 / cW2 + 4.0 * cW2) * B0_Mt2_Mw2_Mw2_Mz2
            + (11.0 - 8.0 * cW2) * log(cW2) / 6.0 / sW2
            - (11.0 / 3.0 + 1.0 / 3.0 / cW2) * log(zt);
    return f2a;
}

double EWSMTwoLoopEW::DeltaEta2(const double Mw_i) const
{
    double Mz = cache.getSM().getMz(), Mz2 = Mz*Mz;
    double Mw = Mw_i, Mw2 = Mw*Mw;
    double cW2 = cache.getSM().cW2(Mw);
    double cW4 = cW2*cW2, cW6 = cW4*cW2;
    double Mt = cache.getSM().getMtpole(), Mt2 = Mt*Mt;
    double mh = cache.getSM().getMHl(), mh2 = mh*mh;
    double zt = Mz2 / Mt2;
    double zt2 = zt*zt, zt3 = zt2*zt;
    double ht = mh2 / Mt2;
    double ht2 = ht*ht, ht3 = ht2*ht, ht4 = ht3*ht, ht5 = ht4*ht;
    double mu = Mt; // renormalization scale

    double dEta2;
    if (mh < 0.57 * Mt) {
        double B0_Mt2_Mz2_Mw2_Mw2 = cache.getPV().B0(Mt2, Mz2, Mw2, Mw2).real();
        double B0_Mt2_Mz2_mh2_Mz2 = cache.getPV().B0(Mt2, Mz2, mh2, Mz2).real();
        dEta2 = (ht3 - 6.0 * ht2 * zt + 11.0 * ht * zt2) / 9.0 / cW2 / (ht - 4.0 * zt) / zt2
                + (49.0 - 289.0 * cW2 - 349.0 * cW4 + 292.0 * cW6) / 216.0 / cW2 / (1.0 - 4.0 * cW2)
                + (1.0 + 18.0 * cW2 - 16.0 * cW4) / 12.0 / (1.0 - 4.0 * cW2) * log(cW2)
                - (17.0 - 40.0 * cW2 + 32.0 * cW4) / 54.0 / cW2 * (sqrt(ht) * M_PI - 2.0)
                + (11.0 * ht2 * zt - 2.0 * ht3 - 24.0 * ht * zt2 + 24.0 * zt3)
                / 18.0 / cW2 / (ht - 4.0 * zt) / zt2 * log(ht)
                + (1.0 - 4.0 * cW2 + 44.0 * cW4 - 32.0 * cW6) / 24.0 / cW2 / (1.0 - 4.0 * cW2)
                * B0_Mt2_Mz2_Mw2_Mw2
                + (13.0 * ht2 * zt - 2.0 * ht3 - 32.0 * ht * zt2 + 36.0 * zt3)
                / 18.0 / cW2 / (ht - 4.0 * zt) / zt2 * B0_Mt2_Mz2_mh2_Mz2
                - (17.0 - 34.0 * cW2 + 26.0 * cW4) / 36.0 / cW2 * log(Mt2 / mu / mu)
                + (ht * (2.0 * ht - 5.0 * zt) / 18.0 / cW2 / zt / (ht - 4.0 * zt)
                + (10.0 - 39.0 * cW2 - 70.0 * cW4 + 48.0 * cW6) / 36.0 / cW2 / (4.0 * cW2 - 1.0))
                * log(zt);
    } else {
        double B0_Mt2_Mz2_Mw2_Mw2 = cache.getPV().B0(Mt2, Mz2, Mw2, Mw2).real();
        dEta2 = (-17.0 + 40.0 * cW2 - 32.0 * cW4) * ht / 216.0 / cW2
                + 5.0 / 144.0 / cW2 / (ht - 4.0)
                + (707.0 - 4720.0 * cW2 + 5900.0 * cW4 - 3696.0 * cW6) / 864.0 / cW2 / (1.0 - 4.0 * cW2)
                + (10.0 / 27.0 - 17.0 / 108.0 / cW2 - 8.0 * cW2 / 27.0)*(1.0 - ht / 4.0) * sqrt(ht) * g(ht)
                + (1.0 + 18.0 * cW2 - 16.0 * cW4) / 12.0 / (1.0 - 4.0 * cW2) * log(cW2)
                + (4.0 - ht) / 12.0 / cW2 / ht * Lambda(ht)
                + (2.0 - 7.0 * cW2 - 70.0 * cW4 + 48.0 * cW6) / 36.0 / cW2 / (4.0 * cW2 - 1.0) * log(zt)
                + (1.0 - 4.0 * cW2 + 44.0 * cW4 - 32.0 * cW6) / 24.0 / cW2 / (1.0 - 4.0 * cW2) * B0_Mt2_Mz2_Mw2_Mw2
                - (17.0 - 34.0 * cW2 + 26.0 * cW4) / 36.0 / cW2 * log(Mt2 / mu / mu)
                + ((4.0 * cW2 - 5.0)*(6.0 + 27.0 * ht - 10.0 * ht2 + ht3) / 54.0 / (ht - 4.0)
                - (1152.0 + 606.0 * ht + 1467.0 * ht2 - 1097.0 * ht3 + 238.0 * ht4 - 17.0 * ht5)
                / 432.0 / cW2 / (ht - 4.0) / (ht - 4.0) / ht) * log(ht)
                + ((5.0 - 4.0 * cW2)*(ht - 1.0) / 9.0 / (ht - 4.0) / ht
                - (384.0 + 10.0 * ht - 238.0 * ht2 + 63.0 * ht3 - 3.0 * ht4)
                / 144.0 / cW2 / (ht - 4.0) / (ht - 4.0) / ht2) * phi(ht / 4.0);
    }
    return dEta2;
}

gslpp::complex EWSMTwoLoopEW::DeltaEta2Add_tmp(const double I3f, const double Qf,
        const double Mw_i) const
{
    double Mz = cache.getSM().getMz(), Mz2 = Mz*Mz;
    double Mw = Mw_i, Mw2 = Mw*Mw;
    double cW2 = cache.getSM().cW2(Mw);
    double cW4 = cW2*cW2, cW6 = cW4*cW2;
    double Mt = cache.getSM().getMtpole(), Mt2 = Mt*Mt;
    double zt = Mz2 / Mt2;
    double mu = Mt; // renormalization scale

    double B0_Mt2_Mz2_Mw2_Mw2 = cache.getPV().B0(Mt2, Mz2, Mw2, Mw2).real();

    gslpp::complex dEta2add = 16.0 * M_PI * M_PI * DeltaEtaf1(I3f, Qf, Mw_i) + Vadd(I3f, Qf, Mw_i)
            - (197.0 - 1378.0 * cW2 + 1064.0 * cW4) / 27.0 / (1.0 - 4.0 * cW2)
            - (1.0 + 16.0 * cW2 - 20.0 * cW4 + 48.0 * cW6) / 3.0 / (1.0 - 4.0 * cW2)
            * B0_Mt2_Mz2_Mw2_Mw2
            - 2.0 * cW2 * (1.0 + 26.0 * cW2 + 24.0 * cW4) / 3.0 / (1.0 - 4.0 * cW2)
            * log(cW2)
            + (41.0 / 3.0 - 46.0 * cW2 / 3.0) * log(Mt2 / mu / mu)
            + 2.0 * (50.0 - 283.0 * cW2 + 242.0 * cW4 - 72.0 * cW6)
            / 9.0 / (1.0 - 4.0 * cW2) * log(zt);
    return dEta2add;
}

gslpp::complex EWSMTwoLoopEW::DeltaEta2Add_f(const Particle f, const double Mw_i) const
{
    double I3f = cache.I3_f(f);
    double Qf = cache.Q_f(f);
    return DeltaEta2Add_tmp(I3f, Qf, Mw_i);
}

double EWSMTwoLoopEW::DeltaKappa2(const double Mw_i) const
{
    double Mz = cache.getSM().getMz(), Mz2 = Mz*Mz;
    double Mw = Mw_i, Mw2 = Mw*Mw;
    double cW2 = cache.getSM().cW2(Mw), sW2 = cache.getSM().sW2(Mw);
    double Mt = cache.getSM().getMtpole(), Mt2 = Mt*Mt;
    double mh = cache.getSM().getMHl(), mh2 = mh*mh;
    double zt = Mz2 / Mt2;
    double ht = mh2 / Mt2;
    double ht2 = ht*ht, ht3 = ht2*ht;
    double mu = Mt; // renormalization scale

    double B0_Mt2_Mz2_Mw2_Mw2 = cache.getPV().B0(Mt2, Mz2, Mw2, Mw2).real();

    double dKappa2;
    if (mh < 0.57 * Mt) {
        dKappa2 = (-175.0 + 366.0 * sW2) / 432.0
                + (3.0 / 8.0 - sW2 / 3.0) * B0_Mt2_Mz2_Mw2_Mw2 - cW2 / 6.0 * log(cW2)
                - 2.0 * M_PI / 27.0 * sqrt(ht)*(8.0 * sW2 - 3.0)
                - (1.0 / 4.0 + 2.0 / 9.0 * sW2) * log(Mt2 / mu / mu)
                + (3.0 * sW2 - 2.0) / 18.0 * log(zt);
    } else {
        dKappa2 = (-211.0 + 24.0 * ht + 462.0 * sW2 - 64.0 * ht * sW2) / 432.0
                + (3.0 / 8.0 - sW2 / 3.0) * B0_Mt2_Mz2_Mw2_Mw2
                - cW2 / 6.0 * log(cW2)
                + (ht - 4.0) * sqrt(ht)*(8.0 * sW2 - 3.0) * g(ht) / 108.0
                - (6.0 + 27.0 * ht - 10.0 * ht2 + ht3)*(3.0 - 8.0 * sW2)
                / 108.0 / (ht - 4.0) * log(ht)
                - (1.0 / 4.0 + 2.0 / 9.0 * sW2) * log(Mt2 / mu / mu)
                + (3.0 * sW2 - 2.0) / 18.0 * log(zt)
                + (ht - 1.0)*(8.0 * sW2 - 3.0) / 18.0 / (4.0 - ht) / ht * phi(ht / 4.0);
    }
    return dKappa2;
}

gslpp::complex EWSMTwoLoopEW::DeltaKappa2Add_tmp(const double I3f, const double Qf,
        const double Mw_i) const
{
    double Mz = cache.getSM().getMz(), Mz2 = Mz*Mz;
    double Mw = Mw_i;
    double cW2 = cache.getSM().cW2(Mw);
    double cW4 = cW2*cW2;
    double Mt = cache.getSM().getMtpole(), Mt2 = Mt*Mt;
    double zt = Mz2 / Mt2;
    double mu = Mt; // renormalization scale

    gslpp::complex i = gslpp::complex::i();

    gslpp::complex dKappa2add = -238.0 * cW2 / 27.0 + 8.0 * cW4
            - 2.0 * cW2 * sqrt(4.0 * cW2 - 1.0)*(3.0 + 4.0 * cW2)
            * atan(1.0 / sqrt(4.0 * cW2 - 1.0))
            - 16.0 / 9.0 * cW2 * log(zt)
            + (1.0 - 12.0 * I3f * Qf + 8.0 * Qf * Qf - 8.0 * cW4 * Qf * Qf)
            / 4.0 / cW2 * FV(1)
            + 4.0 * cW2 * GV(1.0 / cW2) - 7.0 * cW2 * log(cW2)
            - 17.0 / 3.0 * cW2 * log(mu * mu / Mz2)
            + cW2 * (1.0 - 2.0 * Qf * I3f) * FV(1.0 / cW2)
            - i * 80.0 / 9.0 * M_PI;
    return dKappa2add;
}

gslpp::complex EWSMTwoLoopEW::DeltaKappa2Add_f(const Particle f, const double Mw_i) const
{
    double I3f = cache.I3_f(f);
    double Qf = cache.Q_f(f);
    return DeltaKappa2Add_tmp(I3f, Qf, Mw_i);
}

gslpp::complex EWSMTwoLoopEW::Vadd(const double I3f, const double Qf,
        const double Mw_i) const
{
    double Mw = Mw_i, Mw2 = Mw*Mw;
    double cW2 = cache.getSM().cW2(Mw), sW2 = cache.getSM().sW2(Mw);
    double Mt = cache.getSM().getMtpole();
    double mu = Mt; // renormalization scale

    gslpp::complex V = 8.0 * cW2 * log(Mw2 / mu / mu)
            + 3.0 * (2.0 * I3f * Qf - 4.0 * sW2 * Qf * Qf) * FV(1.0)
            - 16.0 * cW2 * GV(1.0 / cW2)
            + (1.0 - 4.0 * cW2 - 4.0 * (1.0 - 2.0 * cW2) * I3f * Qf) * FV(1.0 / cW2);
    return V;
}

gslpp::complex EWSMTwoLoopEW::DeltaEtaf1(const double I3f, const double Qf,
        const double Mw_i) const
{
    double Mz = cache.getSM().getMz(), Mz2 = Mz*Mz;
    double Mw = Mw_i;
    double cW2 = cache.getSM().cW2(Mw);

    gslpp::complex SigmaPrime_ZZ = myOneLoopEW.SigmabarPrime_ZZ_bos_Mz2(Mz, Mw)
            + myOneLoopEW.SigmabarPrime_ZZ_fer_Mz2(Mz, Mw);

    gslpp::complex dEtaf1 = 1.0 / 16.0 / M_PI / M_PI
            * (-SigmaPrime_ZZ / cW2 - 4.0 * cW2 * log(cW2) + Vfi(I3f, Qf, Mz2, Mw));
    return dEtaf1;
}

gslpp::complex EWSMTwoLoopEW::Vfi(const double I3f, const double Qf,
        const double q2, const double Mw_i) const
{
    double I3i = cache.I3_f(cache.getSM().getLeptons(StandardModel::ELECTRON));
    double Qi = cache.Q_f(cache.getSM().getLeptons(StandardModel::ELECTRON));
    double I3aQaQa = I3i * Qi * Qi + I3f * Qf*Qf;
    double I3aQa = I3i * Qi + I3f*Qf;
    double QaQa = Qi * Qi + Qf*Qf;
    double Mz = cache.getSM().getMz(), Mz2 = Mz*Mz;
    double Mw = Mw_i, Mw2 = Mw*Mw;
    double cW2 = cache.getSM().cW2(Mw), sW2 = cache.getSM().sW2(Mw);

    gslpp::complex V = (1.0 - sW2 * (2.0 - 2.0 * I3aQaQa)) * FV(q2 / Mw2)
            + 8.0 * cW2 * GV(q2 / Mw2)
            - (1.0 - 6.0 * sW2 * I3aQa + 6.0 * sW2 * QaQa) / 2.0 / cW2 * FV(q2 / Mz2);
    return V;
}

double EWSMTwoLoopEW::Lambda(const double x) const
{
    if (x >= 0.0 && x <= 4.0) {
        return ( -1.0 / 2.0 / sqrt(x) * g(x) + M_PI / 2.0 * sqrt(4.0 / x - 1.0));
    } else if (x > 4.0) {
        return ( -1.0 / 2.0 / sqrt(x) * g(x));
    } else
        throw std::runtime_error("Out of range in EWSMTwoLoopEW::Lambda()");
}

double EWSMTwoLoopEW::phi(const double x) const
{
    if (x >= 0.0 && x <= 1.0) {
        return ( 4.0 * sqrt(x / (1.0 - x)) * cache.getClausen().Cl2(2.0 * asin(sqrt(x))));
    } else if (x > 1.0) {
        double lambda = sqrt(1.0 - 1.0 / x);
        return ( 1.0 / lambda * (-4.0 * cache.getPolyLog().Li2((1.0 - lambda) / 2.0).real()
                + 2.0 * pow(log((1.0 - lambda) / 2.0), 2.0)
                - pow(log(4.0 * x), 2.0) + M_PI * M_PI / 3.0));
    } else
        throw std::runtime_error("Out of range in EWSMTwoLoopEW::phi()");
}

gslpp::complex EWSMTwoLoopEW::FV(const double x) const
{
    if (x <= 0.0)
        throw std::runtime_error("Out of range in EWSMTwoLoopEW::FV()");

    gslpp::complex i = gslpp::complex::i();

    return ( i * M_PI * (2.0 / x + 3.0 - 2.0 * (1.0 + 1.0 / x)*(1.0 + 1.0 / x) * log(1.0 + x))
            + 2.0 / x + 7.0 / 2.0 - (3.0 + 2.0 / x) * log(x)
            + (1.0 + 1.0 / x)*(1.0 + 1.0 / x)
            *(2.0 * cache.getPolyLog().Li2(1.0 / (1.0 + x))
            - M_PI * M_PI / 3.0 + pow(log(1.0 + x), 2.0)));
}

gslpp::complex EWSMTwoLoopEW::GV(const double x) const
{
    if (x <= 0.0 || x >= 4.0)
        throw std::runtime_error("Out of range in EWSMTwoLoopEW::GV()");

    double atanX = atan(sqrt(x / (4.0 - x)));

    return ( (sqrt((4.0 - x) / x) * atanX - 1.0)*(1.0 / x + 1.0 / 2.0) + 9.0 / 8.0
            + 1.0 / 2.0 / x - (1.0 + 1.0 / 2.0 / x)*4.0 / x * atanX * atanX);
}




