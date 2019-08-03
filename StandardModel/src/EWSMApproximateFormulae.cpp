/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EWSMApproximateFormulae.h"
#include <stdexcept>
#include <sstream>

#define UpperBoundForApproximateFormulae 1000.0
//#define UpperBoundForApproximateFormulae 1500.0 // for test

EWSMApproximateFormulae::EWSMApproximateFormulae(const EWSMcache& cache_i)
: mycache(cache_i)
{
}


////////////////////////////////////////////////////////////////////////

double EWSMApproximateFormulae::Mw() const
{
    // Parametrization from arXiv:hep-ph/0311148v2 (updates from the journal version)
    double Mw0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11;
    if (mycache.getSM().getMHl() >= 100.0 && mycache.getSM().getMHl() <= UpperBoundForApproximateFormulae) {
        // applicable for 100 GeV <= mHl <= 1 TeV
        Mw0 = 80.3779;
        c1 = 0.05263;
        c2 = 0.010239;
        c3 = 0.000954;
        c4 = -0.000054;
        c5 = 1.077;
        c6 = 0.5252;
        c7 = 0.0700;
        c8 = 0.004102;
        c9 = 0.000111;
        c10 = 0.0774;
        c11 = 115.0;
    } else if (mycache.getSM().getMHl() >= 10.0 && mycache.getSM().getMHl() <= 1000.0) {
        // applicable for 10 GeV <= mHl <= 1 TeV
        Mw0 = 80.3799;
        c1 = 0.05427;
        c2 = 0.008931;
        c3 = 0.0000882;
        c4 = 0.000161;
        c5 = 1.070;
        c6 = 0.5237;
        c7 = 0.0679;
        c8 = 0.00179;
        c9 = 0.0000664;
        c10 = 0.0795;
        c11 = 114.9;
    } else {
        std::stringstream out;
        out << mycache.getSM().getMHl();
        throw std::runtime_error("ApproximateFormulae::Mw(): mh=" + out.str() + " is out of range");
    }

    double dH = log(mycache.getSM().getMHl() / 100.0);
    double dh = pow((mycache.getSM().getMHl() / 100.0), 2.0);
    double dt = pow((mycache.getSM().getMtpole() / 174.3), 2.0) - 1.0;
    double dZ = mycache.getSM().getMz() / 91.1875 - 1.0;
    double dalphae = mycache.getSM().DeltaAlphaL5q() / 0.05907 - 1.0;
    double dalphas = mycache.getSM().getAlsMz() / 0.119 - 1.0;

    return (Mw0 - c1 * dH - c2 * dH * dH + c3 * pow(dH, 4.0)
            + c4 * (dh - 1.0) - c5 * dalphae + c6 * dt - c7 * dt * dt
            - c8 * dH * dt + c9 * dh * dt - c10 * dalphas + c11 * dZ
            + mycache.getSM().getDelMw());
}

double EWSMApproximateFormulae::sin2thetaEff_l(const QCD::lepton l) const
{
    // applicable for 10 GeV <= mHl <= 1 TeV
    if (mycache.getSM().getMHl() < 10.0 || mycache.getSM().getMHl() > UpperBoundForApproximateFormulae) {
        std::stringstream out;
        out << mycache.getSM().getMHl();
        throw std::runtime_error("ApproximateFormulae::sin2thetaEff_l(): mh=" + out.str() + " is out of range");
    }

    double s0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10;
    switch (l) {
        case StandardModel::NEUTRINO_1:
        case StandardModel::NEUTRINO_2:
        case StandardModel::NEUTRINO_3:
            s0 = 0.2308772;
            d1 = 4.713 * 0.0001;
            d2 = 2.05 * 0.00001;
            d3 = 3.85 * 0.000001;
            d4 = -1.85 * 0.000001;
            d5 = 2.06 * 0.01;
            d6 = -2.850 * 0.001;
            d7 = 1.82 * 0.0001;
            d8 = -9.71 * 0.000001;
            d9 = 3.96 * 0.0001;
            d10 = -6.54 * 0.1;
            break;
        case StandardModel::ELECTRON:
        case StandardModel::MU:
        case StandardModel::TAU:
            s0 = 0.2312527;
            d1 = 4.729 * 0.0001;
            d2 = 2.07 * 0.00001;
            d3 = 3.85 * 0.000001;
            d4 = -1.85 * 0.000001;
            d5 = 2.07 * 0.01;
            d6 = -2.851 * 0.001;
            d7 = 1.82 * 0.0001;
            d8 = -9.74 * 0.000001;
            d9 = 3.98 * 0.0001;
            d10 = -6.55 * 0.1;
            break;
        default:
            throw std::runtime_error("Error in ApproximateFormulae::sin2thetaEff_l()");
    }

    double L_H = log(mycache.getSM().getMHl() / 100.0);
    double Delta_H = mycache.getSM().getMHl() / 100.0;
    double Delta_ale = mycache.getSM().DeltaAlphaL5q() / 0.05907 - 1.0;
    double Delta_t = pow((mycache.getSM().getMtpole() / 178.0), 2.0) - 1.0;
    double Delta_alphas = mycache.getSM().getAlsMz() / 0.117 - 1.0;
    double Delta_Z = mycache.getSM().getMz() / 91.1876 - 1.0;

    return (s0 + d1 * L_H + d2 * L_H * L_H + d3 * pow(L_H, 4.0)
            + d4 * (Delta_H * Delta_H - 1.0) + d5 * Delta_ale + d6 * Delta_t
            + d7 * Delta_t * Delta_t + d8 * Delta_t * (Delta_H - 1.0)
            + d9 * Delta_alphas + d10 * Delta_Z
            + mycache.getSM().getDelSin2th_l());
}

double EWSMApproximateFormulae::sin2thetaEff_q(const QCD::quark q) const
{
    // applicable for 10 GeV <= mHl <= 1 TeV
    if (mycache.getSM().getMHl() < 10.0 || mycache.getSM().getMHl() > UpperBoundForApproximateFormulae) {
        std::stringstream out;
        out << mycache.getSM().getMHl();
        throw std::runtime_error("ApproximateFormulae::sin2thetaEff_q(): mh=" + out.str() + " is out of range");
    }

    double s0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10;
    double ThError = 0.0; // Theoretical uncertainty
    
    switch (q) {
        case QCD::UP:
        case QCD::CHARM:
            s0 = 0.2311395;
            d1 = 4.726 * 0.0001;
            d2 = 2.07 * 0.00001;
            d3 = 3.85 * 0.000001;
            d4 = -1.85 * 0.000001;
            d5 = 2.07 * 0.01;
            d6 = -2.853 * 0.001;
            d7 = 1.83 * 0.0001;
            d8 = -9.73 * 0.000001;
            d9 = 3.98 * 0.0001;
            d10 = -6.55 * 0.1;
            ThError = mycache.getSM().getDelSin2th_q();
            break;
        case QCD::DOWN:
        case QCD::STRANGE:
            s0 = 0.2310286;
            d1 = 4.720 * 0.0001;
            d2 = 2.06 * 0.00001;
            d3 = 3.85 * 0.000001;
            d4 = -1.85 * 0.000001;
            d5 = 2.07 * 0.01;
            d6 = -2.848 * 0.001;
            d7 = 1.81 * 0.0001;
            d8 = -9.73 * 0.000001;
            d9 = 3.97 * 0.0001;
            d10 = -6.55 * 0.1;
            ThError = mycache.getSM().getDelSin2th_q();
            break;
        case QCD::BOTTOM:
            
            if (mycache.getSM().getMHl() < 120.1 || mycache.getSM().getMHl() > 130.1) {
                s0 = 0.2327580;
                d1 = 4.749 * 0.0001;
                d2 = 2.03 * 0.00001;
                d3 = 3.94 * 0.000001;
                d4 = -1.84 * 0.000001;
                d5 = 2.08 * 0.01;
                d6 = -0.993 * 0.001;
                d7 = 0.708 * 0.0001;
                d8 = -7.61 * 0.000001;
                d9 = 4.03 * 0.0001;
                d10 = 6.61 * 0.1;
                ThError = mycache.getSM().getDelSin2th_b();
                break;

            } else {
                return sin2thetaEff_b();
            }
        case QCD::TOP:
            return 0.0;
        default:
            throw std::runtime_error("Error in ApproximateFormulae::sin2thetaEff_q()");
    }

    double L_H = log(mycache.getSM().getMHl() / 100.0);
    double Delta_H = mycache.getSM().getMHl() / 100.0;
    double Delta_ale = mycache.getSM().DeltaAlphaL5q() / 0.05907 - 1.0;
    double Delta_t = pow((mycache.getSM().getMtpole() / 178.0), 2.0) - 1.0;
    double Delta_alphas = mycache.getSM().getAlsMz() / 0.117 - 1.0;
    double Delta_Z = mycache.getSM().getMz() / 91.1876 - 1.0;

    return (s0 + d1 * L_H + d2 * L_H * L_H + d3 * pow(L_H, 4.0)
            + d4 * (Delta_H * Delta_H - 1.0) + d5 * Delta_ale + d6 * Delta_t
            + d7 * Delta_t * Delta_t + d8 * Delta_t * (Delta_H - 1.0)
            + d9 * Delta_alphas + d10 * Delta_Z
            + ThError);
}

double EWSMApproximateFormulae::sin2thetaEff_b() const
{
    // applicable for 120.1 GeV <= mHl <= 130.1 GeV
    if (mycache.getSM().getMHl() < 120.1 || mycache.getSM().getMHl() > 130.1) {
        std::stringstream out;
        out << mycache.getSM().getMHl();
        throw std::runtime_error("ApproximateFormulae::sin2thetaEff_b(): mh=" + out.str() + " is out of range");
    }

    double s0, d1, d2, d3, d4, d5, d6, d7, d8, d9;

    s0 = 0.232704;
    d1 = 4.723 * 0.0001;
    d2 = 1.97 * 0.0001;
    d3 = 2.07 * 0.01;
    d4 = -9.733 * 0.0001;
    d5 = 3.93 * 0.0001;
    d6 = -1.38 * 0.0001;
    d7 = 2.42 * 0.0001;
    d8 = -8.10 * 0.0001;
    d9 = -0.664;

    double L_H = log(mycache.getSM().getMHl() / 125.7);
    double Delta_ale = mycache.getSM().DeltaAlphaL5q() / 0.059 - 1.0;
    double Delta_t = pow((mycache.getSM().getMtpole() / 173.2), 2.0) - 1.0;
    double Delta_alphas = mycache.getSM().getAlsMz() / 0.1184 - 1.0;
    double Delta_Z = mycache.getSM().getMz() / 91.1876 - 1.0;

    return (s0 + d1 * L_H + d2 * L_H * L_H 
            + d3 * Delta_ale
            + d4 * Delta_t + d5 * Delta_t * Delta_t + d6 * Delta_t * L_H
            + d7 * Delta_alphas + d8 * Delta_t * Delta_alphas
            + d9 * Delta_Z
            + mycache.getSM().getDelSin2th_b());
}

double EWSMApproximateFormulae::DeltaR_TwoLoopEW_rem(const double Mw_i) const
{
    // applicable for 10 GeV <= mHl <= 1 TeV
    if (mycache.getSM().getMHl() < 10.0 || mycache.getSM().getMHl() > UpperBoundForApproximateFormulae) {
        std::stringstream out;
        out << mycache.getSM().getMHl();
        throw std::runtime_error("ApproximateFormulae::DeltaR_TwoLoopEW_rem(): mh=" + out.str() + " is out of range");
    }

    double r0 = 0.003354;
    double r1 = -0.000209;
    double r2 = 0.0000254;
    double r3 = -0.00000785;
    double r4 = -0.00000233;
    double r5 = 0.00783;
    double r6 = 0.00338;
    double r7 = -0.00000989;
    double r8 = 0.0939;
    double r9 = 0.204;
    double r10 = -0.103;

    //double Mw = Mw(DeltaAlphaL5q_i); /* for test */
    double Mw = Mw_i;

    double L_H = log(mycache.getSM().getMHl() / 100.0);
    double Delta_H = mycache.getSM().getMHl() / 100.0;
    double Delta_t = pow((mycache.getSM().getMtpole() / 178.0), 2.0) - 1.0;
    double Delta_Z = mycache.getSM().getMz() / 91.1876 - 1.0;
    double Delta_W = Mw / 80.404 - 1.0;

    return ( r0 + r1 * L_H + r2 * L_H * L_H + r3 * pow(L_H, 4.0)
            + r4 * (Delta_H * Delta_H - 1.0) + r5 * Delta_t
            + r6 * Delta_t * Delta_t + r7 * Delta_t * L_H + r8 * Delta_W
            + r9 * Delta_W * Delta_t + r10 * Delta_Z);
}

double EWSMApproximateFormulae::DeltaKappa_l_TwoLoopEW_rem(const double Mw_i) const
{
    // applicable for 10 GeV <= mHl <= 1 TeV
    if (mycache.getSM().getMHl() < 10.0 || mycache.getSM().getMHl() > UpperBoundForApproximateFormulae) {
        std::stringstream out;
        out << mycache.getSM().getMHl();
        throw std::runtime_error("ApproximateFormulae::DeltaKappa_l_TwoLoopEW_rem(): mh=" + out.str() + " is out of range");
    }

    double k0 = -0.002711;
    double k1 = -0.0000312;
    double k2 = -0.0000412;
    double k3 = 0.00000528;
    double k4 = 0.00000375;
    double k5 = -0.00516;
    double k6 = -0.00206;
    double k7 = -0.000232;
    double k8 = -0.0647;
    double k9 = -0.129;
    double k10 = 0.0712;

    double L_H = log(mycache.getSM().getMHl() / 100.0);
    double Delta_H = mycache.getSM().getMHl() / 100.0;
    double Delta_t = pow((mycache.getSM().getMtpole() / 178.0), 2.0) - 1.0;
    double Delta_Z = mycache.getSM().getMz() / 91.1876 - 1.0;
    double Delta_W = Mw_i / 80.404 - 1.0;

    return ( k0 + k1 * L_H + k2 * L_H * L_H + k3 * pow(L_H, 4.0)
            + k4 * (Delta_H * Delta_H - 1.0) + k5 * Delta_t
            + k6 * Delta_t * Delta_t + k7 * Delta_t * L_H
            + k8 * Delta_W + k9 * Delta_W * Delta_t + k10 * Delta_Z);
}

double EWSMApproximateFormulae::DeltaKappa_b_TwoLoopEW_rem(const double Mw_i) const
{
    // applicable for 10 GeV <= mHl <= 1 TeV
    if (mycache.getSM().getMHl() < 10.0 || mycache.getSM().getMHl() > UpperBoundForApproximateFormulae) {
        std::stringstream out;
        out << mycache.getSM().getMHl();
        throw std::runtime_error("ApproximateFormulae::DeltaKappa_b_TwoLoopEW_rem(): mh=" + out.str() + " is out of range");
    }

    double k0 = -0.002666;
    double k1 = -0.0000592;
    double k2 = -0.00000329;
    double k3 = 0.00000349;
    double k4 = 0.00000283;
    double k5 = -0.00534;
    double k6 = -0.00210;
    double k7 = -0.000219;
    double k8 = -0.0631;
    double k9 = -0.126;
    double k10 = 0.0647;

    double L_H = log(mycache.getSM().getMHl() / 100.0);
    double Delta_H = mycache.getSM().getMHl() / 100.0;
    double Delta_t = pow((mycache.getSM().getMtpole() / 178.0), 2.0) - 1.0;
    double Delta_Z = mycache.getSM().getMz() / 91.1876 - 1.0;
    double Delta_W = Mw_i / 80.404 - 1.0;

    return ( k0 + k1 * L_H + k2 * L_H * L_H + k3 * pow(L_H, 4.0)
            + k4 * (Delta_H * Delta_H - 1.0) + k5 * Delta_t
            + k6 * Delta_t * Delta_t + k7 * Delta_t * L_H
            + k8 * Delta_W + k9 * Delta_W * Delta_t + k10 * Delta_Z);
}

double EWSMApproximateFormulae::R0_bottom_OLD() const
{
    // applicable for 10 GeV <= mHl <= 1 TeV
    if (mycache.getSM().getMHl() < 10.0 || mycache.getSM().getMHl() > UpperBoundForApproximateFormulae) {
        std::stringstream out;
        out << mycache.getSM().getMHl();
        throw std::runtime_error("ApproximateFormulae::R0_bottom(): mh=" + out.str() + " is out of range");
    }

    /*-----------------------------------------------*/
    /* arXiv:1205.0299v1 by Freitas and Huang */
    /*
    double Rb00 = 0.2147464;
    double c1 =  0.0000221;
    double c2 =  0.0000026;
    double c3 = -0.00000067;
    double c4 =  0.0000000911;
    double c5 =  0.000647;
    double c6 = -0.003239;
    double c7 =  0.0000673;
    double c8 = -0.000324;
    double c9 =  0.0610;
    
    double L_H = log(mycache.getSM().getMHl()/100.0);
    double Delta_H = mycache.getSM().getMHl()/100.0;
    double Delta_ale = mycache.getSM().mycache.getSM().DeltaAlphaL5q()/0.05900 - 1.0;
    double Delta_t = pow((mycache.getSM().getMtpole()/173.2), 2.0) - 1.0;
    double Delta_alphas = mycache.getSM().getAlsMz()/0.1184 - 1.0;
    double Delta_Z = mycache.getSM().getMz()/91.1876 - 1.0;

    return (Rb00 + c1*L_H + c2*L_H*L_H + c3*pow(L_H, 4.0)
            + c4*(Delta_H*Delta_H - 1.0) + c5*Delta_ale + c6*Delta_t
            + c7*Delta_t*L_H + c8*Delta_alphas + c9*Delta_Z );    
     */

    /*-----------------------------------------------*/
    /* arXiv:1205.0299v2 by Freitas and Huang */
    /*
    double Rb00 = 0.2149246;
    double c1 = 2.23 * pow(10.0, -5.); 
    double c2 = 2.6 * pow(10.0, -6.);
    double c3 = -6.8 * pow(10.0, -7.);
    double c4 = 9.19 *pow(10.0, -8.);
    double c5 = 6.58 * pow(10.0, -4.);
    double c6 = -3.363 * pow(10.0, -3.);
    double c7 = 6.74 * pow(10.0, -5.);
    double c8 = -1.688 * pow(10.0, -3.);
    double c9 = -9.26 * pow(10.0, -4.);
    double c10 = 5.93 * pow(10.0, -2.);

    double L_H = log(mycache.getSM().getMHl()/100.0);
    double Delta_H = mycache.getSM().getMHl()/100.0;
    double Delta_ale = mycache.getSM().mycache.getSM().DeltaAlphaL5q()/0.05900 - 1.0;
    double Delta_t = pow((mycache.getSM().getMtpole()/173.2), 2.0) - 1.0;
    double Delta_alphas = mycache.getSM().getAlsMz()/0.1184 - 1.0;
    double Delta_Z = mycache.getSM().getMz()/91.1876 - 1.0;    
    
    return (Rb00 + c1*L_H + c2*L_H*L_H + c3*pow(L_H, 4.0)
            + c4*(Delta_H*Delta_H - 1.0) + c5*Delta_ale + c6*Delta_t
            + c7*Delta_t*L_H + c8*Delta_alphas + c9*Delta_alphas*Delta_alphas 
            + c10*Delta_Z );
     */

    /*-----------------------------------------------*/
    /* arXiv:1205.0299v3 by Freitas and Huang */

    double Rb00 = 0.2154940;
    double c1 = 1.88 * pow(10.0, -5.);
    double c2 = 2.0 * pow(10.0, -6.);
    double c3 = -6.0 * pow(10.0, -7.);
    double c4 = 8.53 * pow(10.0, -8.);
    double c5 = 7.05 * pow(10.0, -4.);
    double c6 = -3.159 * pow(10.0, -3.);
    double c7 = 6.65 * pow(10.0, -5.);
    double c8 = -1.704 * pow(10.0, -3.);
    double c9 = -9.30 * pow(10.0, -4.);
    double c10 = 6.26 * pow(10.0, -2.);

    double L_H = log(mycache.getSM().getMHl() / 100.0);
    double Delta_H = mycache.getSM().getMHl() / 100.0;
    double Delta_ale = mycache.getSM().DeltaAlphaL5q() / 0.05900 - 1.0;
    double Delta_t = pow((mycache.getSM().getMtpole() / 173.2), 2.0) - 1.0;
    double Delta_alphas = mycache.getSM().getAlsMz() / 0.1184 - 1.0;
    double Delta_Z = mycache.getSM().getMz() / 91.1876 - 1.0;

    /* Debug (parameters in arXiv:1205.0299v3) */
    //double mHpaper = 100.0;
    //double mHpaper = 200.0;
    //double mHpaper = 400.0;
    //double mHpaper = 600.0;
    //double mHpaper = 1000.0;
    //L_H = log(mHpaper/100.0);
    //Delta_H = mHpaper/100.0;
    //Delta_ale = 0.0;//0.05900/0.05900 - 1.0;
    //Delta_t = 0.0;//pow((173.2/173.2), 2.0) - 1.0;
    //Delta_alphas = 0.0;//0.1184/0.1184 - 1.0;
    //Delta_Z = 0.0;//91.1876/91.1876 - 1.0;

    return (Rb00 + c1 * L_H + c2 * L_H * L_H + c3 * pow(L_H, 4.0)
            + c4 * (Delta_H * Delta_H - 1.0) + c5 * Delta_ale + c6 * Delta_t
            + c7 * Delta_t * L_H + c8 * Delta_alphas + c9 * Delta_alphas * Delta_alphas
            + c10 * Delta_Z);

}

double EWSMApproximateFormulae::Gu_over_Gb_OLD() const
{
    // applicable for 10 GeV <= mHl <= 1 TeV
    if (mycache.getSM().getMHl() < 10.0 || mycache.getSM().getMHl() > UpperBoundForApproximateFormulae) {
        std::stringstream out;
        out << mycache.getSM().getMHl();
        throw std::runtime_error("ApproximateFormulae::Gu_over_Gb(): mh=" + out.str() + " is out of range");
    }

    // obtained from Freitas on Apr. 23, 2013
    /*
    double R = 0.8024769; 
    double c1 = -1.9007e-4;
    double c2 = -2.112e-5;
    double c3 = 6.63e-6;
    double c4 = -1.0284e-6;
    double c5 = -8.081e-3;
    double c6 = 1.830e-4;
    double c7 = 1.7522e-2;
    double c8 = 4.440e-3;
    double c9 = -3.245e-4;
    double c10 = 1.8079e-2;
    double c11 = 1.0720e-2;
    double c12 = -0.129;
     */

    // obtained from Freitas on Sep. 21, 2013
    double R = 0.7997930;
    double c1 = -1.7991e-4;
    double c2 = -1.980e-5;
    double c3 = 6.24e-6;
    double c4 = -0.9829e-6;
    double c5 = -8.200e-3;
    double c6 = 1.657e-4;
    double c7 = 1.6476e-2;
    double c8 = 4.463e-3;
    double c9 = -3.187e-4;
    double c10 = 1.8113e-2;
    double c11 = 1.0720e-2;
    double c12 = -0.144;

    double LH = log(mycache.getSM().getMHl() / 100.0);
    double DH = pow((mycache.getSM().getMHl() / 100.0), 2.0) - 1.0;
    double Dal = mycache.getSM().DeltaAlphaL5q() / 0.059 - 1.0;
    double Dt = pow((mycache.getSM().getMtpole() / 173.2), 2.0) - 1.0;
    double Das = mycache.getSM().getAlsMz() / 0.1184 - 1.0;
    double Dmz = mycache.getSM().getMz() / 91.1876 - 1.0;

    return ( R + c1 * LH + c2 * LH * LH + c3 * pow(LH, 4.0) + c4 * DH + c5 * Dal + c6 * LH * Dal
            + c7 * Dt + c8 * Dt * Dt + c9 * LH * Dt + c10 * Das + c11 * Das * Das + c12 * Dmz);
}

double EWSMApproximateFormulae::Gd_over_Gb_OLD() const
{
    // applicable for 10 GeV <= mHl <= 1 TeV
    if (mycache.getSM().getMHl() < 10.0 || mycache.getSM().getMHl() > UpperBoundForApproximateFormulae) {
        std::stringstream out;
        out << mycache.getSM().getMHl();
        throw std::runtime_error("ApproximateFormulae::Gd_over_Gb(): mh=" + out.str() + " is out of range");
    }

    // obtained from Freitas on Apr. 23, 2013
    /*
    double R = 1.0239191;
    double c1 = -5.093e-5;
    double c2 = -7.08e-6;
    double c3 = 7.4e-7;
    double c4 = 3.27e-8;
    double c5 = 8.68e-4;
    double c6 = 1.064e-4;
    double c7 = 1.8875e-2;
    double c8 = 7.093e-3;
    double c9 = -4.128e-4;
    double c10 = 1.898e-4;
    double c11 = -8.0e-6;
    double c12 = -0.513;
     */

    // obtained from Freitas on Sep. 21, 2013
    double R = 1.0204024;
    double c1 = -2.242e-5;
    double c2 = -1.70e-6;
    double c3 = 2.1e-7;
    double c4 = 6.38e-8;
    double c5 = 5.28e-4;
    double c6 = 0.999e-4;
    double c7 = 1.7539e-2;
    double c8 = 7.138e-3;
    double c9 = -4.041e-4;
    double c10 = 2.290e-4;
    double c11 = -8.0e-6;
    double c12 = -0.530;

    double LH = log(mycache.getSM().getMHl() / 100.0);
    double DH = pow((mycache.getSM().getMHl() / 100.0), 2.0) - 1.0;
    double Dal = mycache.getSM().DeltaAlphaL5q() / 0.059 - 1.0;
    double Dt = pow((mycache.getSM().getMtpole() / 173.2), 2.0) - 1.0;
    double Das = mycache.getSM().getAlsMz() / 0.1184 - 1.0;
    double Dmz = mycache.getSM().getMz() / 91.1876 - 1.0;

    return ( R + c1 * LH + c2 * LH * LH + c3 * pow(LH, 4.0) + c4 * DH + c5 * Dal + c6 * LH * Dal
            + c7 * Dt + c8 * Dt * Dt + c9 * LH * Dt + c10 * Das + c11 * Das * Das + c12 * Dmz);
}

double EWSMApproximateFormulae::X(const std::string observable) const
{
    double LH = log(mycache.getSM().getMHl() / 125.7);
    double Dt = pow(mycache.getSM().getMtpole() / 173.2, 2.0) - 1.0;
    double Das = mycache.getSM().getAlsMz() / 0.1184 - 1.0;
    double Dal = mycache.getSM().DeltaAlphaL5q() / 0.059 - 1.0;
    double DZ = mycache.getSM().getMz() / 91.1876 - 1.0;

    double X0, c1, c2, c3, c4, c5, c6, c7;
    if (observable.compare("Gamma_nu") == 0) {
        X0 = 167.157;
        c1 = -0.055;
        c2 = 1.26;
        c3 = -0.19;
        c4 = -0.02;
        c5 = 0.36;
        c6 = -0.1;
        c7 = 503.0;
    } else if (observable.compare("Gamma_e_mu") == 0) {
        X0 = 83.966;
        c1 = -0.047;
        c2 = 0.807;
        c3 = -0.095;
        c4 = -0.01;
        c5 = 0.25;
        c6 = -1.1;
        c7 = 285.0;
    } else if (observable.compare("Gamma_tau") == 0) {
        X0 = 83.776;
        c1 = -0.047;
        c2 = 0.806;
        c3 = -0.095;
        c4 = -0.01;
        c5 = 0.25;
        c6 = -1.1;
        c7 = 285.0;
    } else if (observable.compare("Gamma_u") == 0) {
        X0 = 299.936;
        c1 = -0.34;
        c2 = 4.07;
        c3 = 14.27;
        c4 = 1.6;
        c5 = 1.8;
        c6 = -11.1;
        c7 = 1253.0;
    } else if (observable.compare("Gamma_c") == 0) {
        X0 = 299.860;
        c1 = -0.34;
        c2 = 4.07;
        c3 = 14.27;
        c4 = 1.6;
        c5 = 1.8;
        c6 = -11.1;
        c7 = 1253.0;
    } else if (observable.compare("Gamma_d_s") == 0) {
        X0 = 382.770;
        c1 = -0.34;
        c2 = 3.83;
        c3 = 10.20;
        c4 = -2.4;
        c5 = 0.67;
        c6 = -10.1;
        c7 = 1469.0;
    } else if (observable.compare("Gamma_b") == 0) {
        X0 = 375.724;
        c1 = -0.30;
        c2 = -2.28;
        c3 = 10.53;
        c4 = -2.4;
        c5 = 1.2;
        c6 = -10.0;
        c7 = 1458.0;
    } else if (observable.compare("GammaZ") == 0) {
        X0 = 2494.24;
        c1 = -2.0;
        c2 = 19.7;
        c3 = 58.60;
        c4 = -4.0;
        c5 = 8.0;
        c6 = -55.9;
        c7 = 9267.0;
    } else if (observable.compare("sigmaHadron") == 0) {
        X0 = 41488.4;
        c1 = 3.0;
        c2 = 60.9;
        c3 = -579.4;
        c4 = 38.0;
        c5 = 7.3;
        c6 = 85.0;
        c7 = -86027.0;
    } else if (observable.compare("R0_lepton") == 0) {
        X0 = 20750.9;
        c1 = -8.1;
        c2 = -39.0;
        c3 = 732.1;
        c4 = -44.0;
        c5 = 5.5;
        c6 = -358.0;
        c7 = 11702.0;
    } else if (observable.compare("R0_charm") == 0) {
        X0 = 172.23;
        c1 = -0.029;
        c2 = 1.0;
        c3 = 2.3;
        c4 = 1.3;
        c5 = 0.38;
        c6 = -1.2;
        c7 = 37.0;
    } else if (observable.compare("R0_bottom") == 0) {
        X0 = 215.80;
        c1 = 0.031;
        c2 = -2.98;
        c3 = -1.32;
        c4 = -0.84;
        c5 = 0.035;
        c6 = 0.73;
        c7 = -18.0;
    } else
        throw std::runtime_error("ApproximateFormulae::X(): " + observable + " is not defined");

    return ( 0.001
            * (X0 + c1 * LH + c2 * Dt + c3 * Das + c4 * Das * Das + c5 * Das * Dt + c6 * Dal + c7 * DZ));
}

double EWSMApproximateFormulae::X_extended(const std::string observable) const
{
    double LH = log(mycache.getSM().getMHl() / 125.7);
    double DH = mycache.getSM().getMHl() / 125.7 - 1.0;
    double Dt = pow(mycache.getSM().getMtpole() / 173.2, 2.0) - 1.0;
    double Das = mycache.getSM().getAlsMz() / 0.1184 - 1.0;
    double Dal = mycache.getSM().DeltaAlphaL5q() / 0.059 - 1.0;
    double DZ = mycache.getSM().getMz() / 91.1876 - 1.0;

    double X0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15;
    double ThError = 0.0; // Theoretical uncertainty
    if (observable.compare("Gamma_nu") == 0) {
        X0 = 167.157;
        a1 = -0.1567;
        a2 = -0.1194;
        a3 = 0.1031;
        a4 = -0.00269;
        a5 = 1.258;
        a6 = -0.13;
        a7 = -0.020;
        a8 = 0.0133;
        a9 = -0.19;
        a10 = -0.018;
        a11 = -0.021;
        a12 = 0.34;
        a13 = -0.084;
        a14 = 0.064;
        a15 = 503.0;
    } else if (observable.compare("Gamma_e_mu") == 0) {
        X0 = 83.966;
        a1 = -0.1017;
        a2 = -0.06352;
        a3 = 0.05500;
        a4 = -0.00145;
        a5 = 0.8051;
        a6 = -0.027;
        a7 = -0.017;
        a8 = 0.0066;
        a9 = -0.095;
        a10 = -0.010;
        a11 = -0.015;
        a12 = 0.23;
        a13 = -1.1;
        a14 = 0.064;
        a15 = 285.0;
    } else if (observable.compare("Gamma_tau") == 0) {
        X0 = 83.776;
        a1 = -0.1016;
        a2 = -0.06339;
        a3 = 0.05488;
        a4 = -0.00145;
        a5 = 0.8036;
        a6 = -0.026;
        a7 = -0.017;
        a8 = 0.0066;
        a9 = -0.095;
        a10 = -0.010;
        a11 = -0.015;
        a12 = 0.23;
        a13 = -1.1;
        a14 = 0.064;
        a15 = 285.0;
    } else if (observable.compare("Gamma_u") == 0) {
        X0 = 299.936;
        a1 = -0.5681;
        a2 = -0.2636;
        a3 = 0.2334;
        a4 = -0.00592;
        a5 = 4.057;
        a6 = -0.50;
        a7 = -0.058;
        a8 = 0.0352;
        a9 = 14.26;
        a10 = 1.6;
        a11 = -0.081;
        a12 = 1.7;
        a13 = -11.1;
        a14 = 0.19;
        a15 = 1251.0;
    } else if (observable.compare("Gamma_c") == 0) {
        X0 = 299.859;
        a1 = -0.5680;
        a2 = -0.2635;
        a3 = 0.2334;
        a4 = -0.00592;
        a5 = 4.056;
        a6 = -0.50;
        a7 = -0.058;
        a8 = 0.0352;
        a9 = 14.26;
        a10 = 1.6;
        a11 = -0.081;
        a12 = 1.7;
        a13 = -11.1;
        a14 = 0.19;
        a15 = 1251.0;
    } else if (observable.compare("Gamma_d_s") == 0) {
        X0 = 382.770;
        a1 = -0.6199;
        a2 = -0.3182;
        a3 = 0.2800;
        a4 = -0.00711;
        a5 = 3.810;
        a6 = -0.25;
        a7 = -0.060;
        a8 = 0.0420;
        a9 = 10.20;
        a10 = -2.4;
        a11 = -0.083;
        a12 = 0.65;
        a13 = -10.1;
        a14 = 0.19;
        a15 = 1468.0;
    } else if (observable.compare("Gamma_b") == 0) {
        X0 = 375.723;
        a1 = -0.5744;
        a2 = -0.3074;
        a3 = 0.2725;
        a4 = -0.00703;
        a5 = -2.292;
        a6 = -0.027;
        a7 = -0.013;
        a8 = 0.0428;
        a9 = 10.53;
        a10 = -2.4;
        a11 = -0.088;
        a12 = 1.2;
        a13 = -10.1;
        a14 = 0.19;
        a15 = 1456.0;
    } else if (observable.compare("Gamma_had") == 0) {
//  Removing leptonic contributions from GammaZ
//        X0 = 1741.06;
//        a1 = -2.9501;
//        a2 = -1.47064;
//        a3 = 1.29906;
//        a4 = -0.033105;
//        a5 = 13.4416;
//        a6 = -1.5285;
//        a7 = -0.249;
//        a8 = 0.19725;
//        a9 = 59.4525;
//        a10 = -4.008;
//        a11 = -0.419;
//        a12 = 5.895;
//        a13 = -52.474;
//        a14 = 0.933;
//        a15 = 6893.;

//  Summing all hadronic contributions        
        X0 = 1741.058;
        a1 = -2.9503;
        a2 = -1.4708999;
        a3 = 1.2993;
        a4 = -0.03308999;
        a5 = 13.440999;
        a6 = -1.527;
        a7 = -0.249;
        a8 = 0.1972;
        a9 = 59.44999;
        a10 = -3.9999;
        a11 = -0.416;
        a12 = 5.9;
        a13 = -52.5;
        a14 = 0.95;
        a15 = 6894.;        
        
    } else if (observable.compare("GammaZ") == 0) {
        X0 = 2494.24;
        a1 = -3.725;
        a2 = -2.019;
        a3 = 1.773;
        a4 = -0.04554;
        a5 = 19.63;
        a6 = -2.0;
        a7 = -0.36;
        a8 = 0.257;
        a9 = 58.60;
        a10 = -4.1;
        a11 = -0.53;
        a12 = 7.6;
        a13 = -56.0;
        a14 = 1.3;
        a15 = 9256.0;
        ThError = mycache.getSM().getDelGammaZ();
    } else if (observable.compare("sigmaHadron") == 0) {
        X0 = 41488.4;
        a1 = 3.88;
        a2 = 0.829;
        a3 = -0.911;
        a4 = 0.0076;
        a5 = 61.10;
        a6 = 16.0;
        a7 = -2.0;
        a8 = -0.59;
        a9 = -579.4;
        a10 = 38.0;
        a11 = -0.26;
        a12 = 6.5;
        a13 = 84.0;
        a14 = 9.5;
        a15 = -86152.0;
        ThError = mycache.getSM().getDelSigma0H();
    } else if (observable.compare("R0_lepton") == 0) {
        X0 = 20750.9;
        a1 = -10.00;
        a2 = -1.83;
        a3 = 1.878;
        a4 = -0.0343;
        a5 = -38.8;
        a6 = -11.0;
        a7 = 1.2;
        a8 = 0.72;
        a9 = 732.1;
        a10 = -44.0;
        a11 = -0.64;
        a12 = 5.6;
        a13 = -357.0;
        a14 = -4.7;
        a15 = 11771.0;
        ThError = mycache.getSM().getDelR0l();
    } else if (observable.compare("R0_electron") == 0) {
        ThError = mycache.getSM().getDelR0l();        
        return X_extended("Gamma_had")/X_extended("Gamma_e_mu") + ThError;

    } else if (observable.compare("R0_muon") == 0) {
        ThError = mycache.getSM().getDelR0l();        
        return X_extended("Gamma_had")/X_extended("Gamma_e_mu") + ThError;

    } else if (observable.compare("R0_tau") == 0) {
        ThError = mycache.getSM().getDelR0l();        
        return X_extended("Gamma_had")/X_extended("Gamma_tau") + ThError;

    } else if (observable.compare("R0_up") == 0) {
        ThError = 0.0; // Set to zero for the moment       
        return X_extended("Gamma_u")/X_extended("Gamma_had") + ThError;

    } else if (observable.compare("R0_strange") == 0) {
        ThError = 0.0; // Set to zero for the moment       
        return X_extended("Gamma_d_s")/X_extended("Gamma_had") + ThError;

    } else if (observable.compare("R0_charm") == 0) {
        X0 = 172.23;
        a1 = -0.034;
        a2 = -0.0058;
        a3 = 0.0054;
        a4 = -0.00012;
        a5 = 1.00;
        a6 = -0.15;
        a7 = -0.0074;
        a8 = 0.00091;
        a9 = 2.3;
        a10 = 1.3;
        a11 = -0.0013;
        a12 = 0.35;
        a13 = -1.2;
        a14 = 0.014;
        a15 = 37.0;
        ThError = mycache.getSM().getDelR0c();
    } else if (observable.compare("R0_bottom") == 0) {
        X0 = 215.80;
        a1 = 0.036;
        a2 = 0.0057;
        a3 = -0.0044;
        a4 = 0.000062;
        a5 = -2.98;
        a6 = 0.20;
        a7 = 0.020;
        a8 = -0.00036;
        a9 = -1.3;
        a10 = -0.84;
        a11 = -0.0019;
        a12 = 0.054;
        a13 = 0.73;
        a14 = -0.011;
        a15 = -18.0;
        ThError = mycache.getSM().getDelR0b();
    } else
        throw std::runtime_error("ApproximateFormulae::X_extended(): " + observable + " is not defined");

    return ( 0.001
            * (X0 + a1 * LH + a2 * LH * LH + a3 * DH + a4 * DH * DH + a5 * Dt + a6 * Dt * Dt
            + a7 * Dt * LH + a8 * Dt * LH * LH + a9 * Das + a10 * Das * Das + a11 * Das * LH
            + a12 * Das * Dt + a13 * Dal + a14 * Dal * LH + a15 * DZ) + ThError);
}


double EWSMApproximateFormulae::X_full_2_loop(const std::string observable) const
{
        
//  For MH not in [85,165] GeV there are significant differences with some predicions
//  of X_extended, which go well beyond the expected size of the bosonic corrections (>~2x).
//  Use EWSMApproximateFormulae::X_extended in that case
    if (mycache.getSM().getMHl() < 85.0 || mycache.getSM().getMHl() > 165.0) {
        return X_extended(observable);
    } 
    
//  Otherwise proceed with the full 2-loop code
    double LH = log(mycache.getSM().getMHl() / 125.7);
    double Dt = pow(mycache.getSM().getMtpole() / 173.2, 2.0) - 1.0;
    double Das = mycache.getSM().getAlsMz() / 0.1184 - 1.0;
    double Dal = mycache.getSM().DeltaAlphaL5q() / 0.059 - 1.0;
    double DZ = mycache.getSM().getMz() / 91.1876 - 1.0;

    double X0, c1, c2, c3, c4, c5, c6, c7;
    double ThError = 0.0; // Theoretical uncertainty
    if (observable.compare("Gamma_nu") == 0) {
        X0 = 167.176;
        c1 = -0.071;
        c2 = 1.26;
        c3 = -0.19;
        c4 = -0.02;
        c5 = 0.36;
        c6 = -0.1;
        c7 = 504.0;

    } else if (observable.compare("Gamma_e_mu") == 0) {
        X0 = 83.983;
        c1 = -0.061;
        c2 = 0.810;
        c3 = -0.096;
        c4 = -0.01;
        c5 = 0.25;
        c6 = -1.1;
        c7 = 286.0;

    } else if (observable.compare("Gamma_tau") == 0) {
        X0 = 83.793;
        c1 = -0.060;
        c2 = 0.810;
        c3 = -0.095;
        c4 = -0.01;
        c5 = 0.25;
        c6 = -1.1;
        c7 = 285.0;

    } else if (observable.compare("Gamma_u") == 0) {
        X0 = 299.993;
        c1 = -0.38;
        c2 = 4.08;
        c3 = 14.27;
        c4 = 1.6;
        c5 = 1.8;
        c6 = -11.1;
        c7 = 1253.0;

    } else if (observable.compare("Gamma_c") == 0) {
        X0 = 299.916;
        c1 = -0.38;
        c2 = 4.08;
        c3 = 14.27;
        c4 = 1.6;
        c5 = 1.8;
        c6 = -11.1;
        c7 = 1253.0;

    } else if (observable.compare("Gamma_d_s") == 0) {
        X0 = 382.828;
        c1 = -0.39;
        c2 = 3.83;
        c3 = 10.20;
        c4 = -2.4;
        c5 = 0.67;
        c6 = -10.1;
        c7 = 1470.0;

    } else if (observable.compare("Gamma_b") == 0) {
        X0 = 375.889;
        c1 = -0.36;
        c2 = -2.14;
        c3 = 10.53;
        c4 = -2.4;
        c5 = 1.2;
        c6 = -10.1;
        c7 = 1459.0;

    } else if (observable.compare("Gamma_had") == 0) {

//  Summing all hadronic contributions        
        X0 = 1741.454;
        c1 = -1.9;
        c2 = 13.68;
        c3 = 59.47;
        c4 = -4.0;
        c5 = 6.14;
        c6 = -52.5;
        c7 = 6905.0;       
        
    } else if (observable.compare("GammaZ") == 0) {
        X0 = 2494.74;
        c1 = -2.3;
        c2 = 19.9;
        c3 = 58.61;
        c4 = -4.0;
        c5 = 8.0;
        c6 = -56.0;
        c7 = 9273.0;

        ThError = mycache.getSM().getDelGammaZ();
    } else if (observable.compare("sigmaHadron") == 0) {
        X0 = 41489.6;
        c1 = 1.6;
        c2 = 60.0;
        c3 = -579.6;
        c4 = 38.0;
        c5 = 7.3;
        c6 = 85.0;
        c7 = -86011.0;

        ThError = mycache.getSM().getDelSigma0H();
    } else if (observable.compare("R0_lepton") == 0) {
        X0 = 20751.6;
        c1 = -7.8;
        c2 = -37.0;
        c3 = 732.3;
        c4 = -44.0;
        c5 = 5.5;
        c6 = -358.0;
        c7 = 11696.0;

        ThError = mycache.getSM().getDelR0l();
    } else if (observable.compare("R0_electron") == 0) {
        ThError = mycache.getSM().getDelR0l();        
        return X_full_2_loop("Gamma_had")/X_full_2_loop("Gamma_e_mu") + ThError;

    } else if (observable.compare("R0_muon") == 0) {
        ThError = mycache.getSM().getDelR0l();        
        return X_full_2_loop("Gamma_had")/X_full_2_loop("Gamma_e_mu") + ThError;

    } else if (observable.compare("R0_tau") == 0) {
        ThError = mycache.getSM().getDelR0l();        
        return X_full_2_loop("Gamma_had")/X_full_2_loop("Gamma_tau") + ThError;

    } else if (observable.compare("R0_neutrino") == 0) {
        ThError = 0.0;        
        return X_full_2_loop("Gamma_nu")/X_full_2_loop("Gamma_had") + ThError;

    } else if (observable.compare("R0_up") == 0) {
        ThError = 0.0; // Set to zero for the moment        
        return X_full_2_loop("Gamma_u")/X_full_2_loop("Gamma_had") + ThError;

    } else if (observable.compare("R0_strange") == 0) {
        ThError = 0.0; // Set to zero for the moment        
        return X_full_2_loop("Gamma_d_s")/X_full_2_loop("Gamma_had") + ThError;

    } else if (observable.compare("R0_charm") == 0) {
        X0 = 172.22;
        c1 = -0.031;
        c2 = 1.0;
        c3 = 2.3;
        c4 = 1.3;
        c5 = 0.38;
        c6 = -1.2;
        c7 = 37.0;

        ThError = mycache.getSM().getDelR0c();
    } else if (observable.compare("R0_bottom") == 0) {
        X0 = 215.85;
        c1 = 0.029;
        c2 = -2.92;
        c3 = -1.32;
        c4 = -0.84;
        c5 = 0.032;
        c6 = 0.72;
        c7 = -18.0;

        ThError = mycache.getSM().getDelR0b();
    } else
        throw std::runtime_error("ApproximateFormulae::X_full_2_loop(): " + observable + " is not defined");

    return ( 0.001
            * (X0 + c1 * LH + c2 * Dt 
            + c3 * Das + c4 * Das * Das
            + c5 * Das * Dt + c6 * Dal + c7 * DZ) + ThError);
}



double EWSMApproximateFormulae::X_full(const std::string observable) const
{
    
//  Full 2-loop implementation
    
    double LH = log(mycache.getSM().getMHl() / 125.7);
    double LH2 = LH * LH;
    
    double DH = pow(125.7 / (mycache.getSM().getMHl()), 4.0) - 1.0;
    
    double Dt = pow(mycache.getSM().getMtpole() / 173.2, 2.0) - 1.0;
    
    double Das = mycache.getSM().getAlsMz() / 0.1184 - 1.0;
    
    double Dal = mycache.getSM().DeltaAlphaL5q() / 0.059 - 1.0;
    
    double DZ = mycache.getSM().getMz() / 91.1876 - 1.0;

    double X0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16;
    
    double ThError = 0.0; // Theoretical uncertainty
    
    if (observable.compare("Gamma_e_mu") == 0) {
        X0 = 83.983;
        a1 = -0.1202;
        a2 = -0.06919;
        a3 = 0.00383;
        a4 = 0.0597;
        a5 = 0.8037;
        a6 = -0.015;
        a7 = -0.0195;
        a8 = 0.0032;
        a9 = -0.0956;
        a10 = -0.0078;
        a11 = -0.0095;
        a12 = 0.25;
        a13 = -1.08;
        a14 = 0.056;
        a15 = -0.37;
        a16 = 286.0;

    } else if (observable.compare("Gamma_tau") == 0) {
        X0 = 83.793;
        a1 = -0.1200;
        a2 = -0.06905;
        a3 = 0.00382;
        a4 = 0.0596;
        a5 = 0.8023;
        a6 = -0.015;
        a7 = -0.0195;
        a8 = 0.0032;
        a9 = -0.0954;
        a10 = -0.0078;
        a11 = -0.0094;
        a12 = 0.25;
        a13 = -1.08;
        a14 = 0.056;
        a15 = -0.37;
        a16 = 285.0;

    } else if (observable.compare("Gamma_nu") == 0) {
        X0 = 167.176;
        a1 = -0.1752;
        a2 = -0.1249;
        a3 = 0.00595;
        a4 = 0.1046;
        a5 = 1.253;
        a6 = -0.110;
        a7 = -0.0232;
        a8 = 0.0064;
        a9 = -0.187;
        a10 = -0.014;
        a11 = -0.014;
        a12 = 0.37;
        a13 = -0.085;
        a14 = 0.054;
        a15 = -0.30;
        a16 = 503.0;

    } else if (observable.compare("Gamma_u") == 0) {
        X0 = 299.994;
        a1 = -0.6152;
        a2 = -0.2771;
        a3 = 0.0174;
        a4 = 0.2341;
        a5 = 4.051;
        a6 = -0.467;
        a7 = -0.0676;
        a8 = 0.017;
        a9 = 14.26;
        a10 = 1.6;
        a11 = -0.046;
        a12 = 1.82;
        a13 = -11.1;
        a14 = 0.16;
        a15 = -1.0;
        a16 = 1253.0;

    } else if (observable.compare("Gamma_c") == 0) {
        X0 = 299.918;
        a1 = -0.6152;
        a2 = -0.2771;
        a3 = 0.0174;
        a4 = 0.2340;
        a5 = 4.051;
        a6 = -0.467;
        a7 = -0.0676;
        a8 = 0.017;
        a9 = 14.26;
        a10 = 1.6;
        a11 = -0.046;
        a12 = 1.82;
        a13 = -11.1;
        a14 = 0.16;
        a15 = -1.0;
        a16 = 1252.0;

    } else if (observable.compare("Gamma_d_s") == 0) {
        X0 = 382.829;
        a1 = -0.6685;
        a2 = -0.3322;
        a3 = 0.0193;
        a4 = 0.2792;
        a5 = 3.792;
        a6 = -0.18;
        a7 = -0.0706;
        a8 = 0.020;
        a9 = 10.20;
        a10 = -2.4;
        a11 = -0.052;
        a12 = 0.71;
        a13 = -10.1;
        a14 = 0.16;
        a15 = -0.92;
        a16 = 1469.0;

    } else if (observable.compare("Gamma_b") == 0) {
        X0 = 375.890;
        a1 = -0.6017;
        a2 = -0.3158;
        a3 = 0.0190;
        a4 = 0.227;
        a5 = -2.174;
        a6 = 0.042;
        a7 = -0.027;
        a8 = 0.021;
        a9 = 10.53;
        a10 = -2.4;
        a11 = -0.056;
        a12 = 1.2;
        a13 = -10.1;
        a14 = 0.15;
        a15 = -0.95;
        a16 = 1458.0;

    } else if (observable.compare("Gamma_had") == 0) {

//  Removing leptonic contributions from the total Z witdh        
        X0 = 1741.46;
        a1 = -3.169;
        a2 = -1.53487;
        a3 = 0.09267;
        a4 = 1.2532;
        a5 = 13.5113;
        a6 = -1.255;
        a7 = -0.3039;
        a8 = 0.0912;
        a9 = 59.4576;
        a10 = -3.9346;
        a11 = -0.2496;
        a12 = 6.24;
        a13 = -52.605;
        a14 = 0.77;
        a15 = -4.79;
        a16 = 6901.0;       
        
    } else if (observable.compare("GammaZ") == 0) {
        X0 = 2494.75;
        a1 = -4.055;
        a2 = -2.117;
        a3 = 0.122;
        a4 = 1.746;
        a5 = 19.68;
        a6 = -1.63;
        a7 = -0.432;
        a8 = 0.12;
        a9 = 58.61;
        a10 = -4.0;
        a11 = -0.32;
        a12 = 8.1;
        a13 = -56.1;
        a14 = 1.1;
        a15 = -6.8;
        a16 = 9267.0;

        ThError = mycache.getSM().getDelGammaZ();
    } else if (observable.compare("sigmaHadron") == 0) {
        X0 = 41489.6;
        a1 = 0.408;
        a2 = -0.320;
        a3 = 0.0424;
        a4 = 1.32;
        a5 = 60.17;
        a6 = 16.3;
        a7 = -2.31;
        a8 = -0.19;
        a9 = -579.58;
        a10 = 38.0;
        a11 = 0.010;
        a12 = 7.5;
        a13 = 85.2;
        a14 = 9.1;
        a15 = -68.0;
        a16 = -85957.0;

        ThError = mycache.getSM().getDelSigma0H();
    } else if (observable.compare("R0_lepton") == 0) {
        X0 = 20751.6;
        a1 = -8.112;
        a2 = -1.174;
        a3 = 0.155;
        a4 = 0.16;
        a5 = -37.59;
        a6 = -10.9;
        a7 = 1.27;
        a8 = 0.29;
        a9 = 732.30;
        a10 = -44.0;
        a11 = -0.61;
        a12 = 5.7;
        a13 = -358.0;
        a14 = -4.7;
        a15 = 37;
        a16 = 11649.0;

        ThError = mycache.getSM().getDelR0l();
    } else if (observable.compare("R0_electron") == 0) {
        ThError = mycache.getSM().getDelR0l();        
        return X_full("Gamma_had")/X_full("Gamma_e_mu") + ThError;

    } else if (observable.compare("R0_muon") == 0) {
        ThError = mycache.getSM().getDelR0l();        
        return X_full("Gamma_had")/X_full("Gamma_e_mu") + ThError;

    } else if (observable.compare("R0_tau") == 0) {
        ThError = mycache.getSM().getDelR0l();        
        return X_full("Gamma_had")/X_full("Gamma_tau") + ThError;

    } else if (observable.compare("R0_neutrino") == 0) {
        ThError = 0.0;        
        return X_full("Gamma_nu")/X_full("Gamma_had") + ThError;

    } else if (observable.compare("R0_up") == 0) {
        ThError = 0.0; // Set to zero for the moment        
        return X_full("Gamma_u")/X_full("Gamma_had") + ThError;

    } else if (observable.compare("R0_strange") == 0) {
        ThError = 0.0; // Set to zero for the moment        
        return X_full("Gamma_d_s")/X_full("Gamma_had") + ThError;

    } else if (observable.compare("R0_charm") == 0) {
        X0 = 172.222;
        a1 = -0.04049;
        a2 = -0.00749;
        a3 = 0.000832;
        a4 = 0.0108;
        a5 = 0.98956;
        a6 = -0.151;
        a7 = -0.00761;
        a8 = 0.00080;
        a9 = 2.309;
        a10 = 1.25;
        a11 = 0.00045;
        a12 = 0.369;
        a13 = -1.20;
        a14 = 0.012;
        a15 = -0.062;
        a16 = 36.67;

        ThError = mycache.getSM().getDelR0c();
    } else if (observable.compare("R0_bottom") == 0) {
        X0 = 215.850;
        a1 = 0.04904;
        a2 = 0.009149;
        a3 = -0.000535;
        a4 = -0.02676;
        a5 = -2.9221;
        a6 = 0.200;
        a7 = 0.0197;
        a8 = -0.0011;
        a9 = -1.319;
        a10 = -0.84;
        a11 = -0.0027;
        a12 = 0.044;
        a13 = 0.719;
        a14 = -0.0077;
        a15 = -0.044;
        a16 = -17.90;

        ThError = mycache.getSM().getDelR0b();
    } else
        throw std::runtime_error("ApproximateFormulae::X_full(): " + observable + " is not defined");

    return ( 0.001
            * ( X0 + a1 * LH + a2 * LH2 + a3 * LH2 * LH2 + a4 *DH  
            + a5 * Dt + a6 * Dt*Dt + a7 * Dt*LH + a8 * Dt * LH2  
            + a9 * Das + a10 * Das * Das + a11 * Das*DH + a12*Das*Dt
            + a13 * Dal + a14 * Dal * DH + a15 * Dal * Dt 
            + a16 * DZ ) + ThError);
}


double EWSMApproximateFormulae::sin2thetaEff_b_full() const
{
    // applicable for 25 GeV <= mHl <= 225 GeV. Remove boundaries for the moment
    //if (mycache.getSM().getMHl() < 25.0 || mycache.getSM().getMHl() > 225.0) {
    //    std::stringstream out;
    //    out << mycache.getSM().getMHl();
    //    throw std::runtime_error("ApproximateFormulae::sin2thetaEff_b_full(): mh=" + out.str() + " is out of range");
    //}

//  Full 2-loop implementation
    
    double LH = log(mycache.getSM().getMHl() / 125.7);
    double LH2 = LH * LH;
    
    double DH = pow(125.7 / (mycache.getSM().getMHl()), 4.0) - 1.0;
    
    double Dt = pow(mycache.getSM().getMtpole() / 173.2, 2.0) - 1.0;
    
    double Das = mycache.getSM().getAlsMz() / 0.1184 - 1.0;
    
    double Dal = mycache.getSM().DeltaAlphaL5q() / 0.059 - 1.0;
    
    double DZ = mycache.getSM().getMz() / 91.1876 - 1.0;

    double X0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10;
    
    double ThError = 0.0; // Theoretical uncertainty
    
    X0 = 2327.04;
    d1 = 4.638;
    d2 = 0.558;
    d3 = -0.0700;
    d4 = 207.0;
    d5 = -9.554;
    d6 = 3.83;
    d7 = 0.179;
    d8 = 2.41;
    d9 = -8.24;
    d10 = -6630.0;
    
    ThError = mycache.getSM().getDelSin2th_b();
    
    return ( 0.0001*( X0 + d1 * LH + d2 * LH2 + d3 * LH2 * LH2 
            + d4 * Dal
            + d5 * Dt + d6 * Dt * Dt 
            + d7 * Dt * LH
            + d8 * Das + d9 * Das * Dt
            + d10 * DZ ) + ThError);
}


double EWSMApproximateFormulae::sin2thetaEff_l_full() const
{
    // applicable for 25 GeV <= mHl <= 225 GeV. Remove boundaries for the moment
    //if (mycache.getSM().getMHl() < 25.0 || mycache.getSM().getMHl() > 225.0) {
    //    std::stringstream out;
    //    out << mycache.getSM().getMHl();
    //    throw std::runtime_error("ApproximateFormulae::sin2thetaEff_l_full(): mh=" + out.str() + " is out of range");
    //}

//  Full 2-loop implementation
    
    double LH = log(mycache.getSM().getMHl() / 125.7);
    double LH2 = LH * LH;
    
    double DH = pow(125.7 / (mycache.getSM().getMHl()), 4.0) - 1.0;
    
    double Dt = pow(mycache.getSM().getMtpole() / 173.2, 2.0) - 1.0;
    
    double Das = mycache.getSM().getAlsMz() / 0.1184 - 1.0;
    
    double Dal = mycache.getSM().DeltaAlphaL5q() / 0.059 - 1.0;
    
    double DZ = mycache.getSM().getMz() / 91.1876 - 1.0;

    double X0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10;
    
    double ThError = 0.0; // Theoretical uncertainty
    
    X0 = 2314.64;
    d1 = 4.616;
    d2 = 0.539;
    d3 = -0.0737;
    d4 = 206.0;
    d5 = -25.71;
    d6 = 4.00;
    d7 = 0.288;
    d8 = 3.88;
    d9 = -6.49;
    d10 = -6560.0;
    
    ThError = mycache.getSM().getDelSin2th_l();
    
    return ( 0.0001*( X0 + d1 * LH + d2 * LH2 + d3 * LH2 * LH2 
            + d4 * Dal
            + d5 * Dt + d6 * Dt * Dt 
            + d7 * Dt * LH
            + d8 * Das + d9 * Das * Dt
            + d10 * DZ ) + ThError);
}