/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EWSMApproximateFormulae.h"
#include <stdexcept>
#include <sstream>

// The latter value should be used in a test for the determination of m_h 
// from A_FB^b. 
#define UpperBoundForApproximateFormulae 1000.0
//#define UpperBoundForApproximateFormulae 1500.0


EWSMApproximateFormulae::EWSMApproximateFormulae(const StandardModel& SM_i) 
: SM(SM_i) 
{
    bDebug = SM_i.isBDebug();
}


////////////////////////////////////////////////////////////////////////

double EWSMApproximateFormulae::Mw(const double DeltaAlphaL5q_i) const 
{
    double Mw0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11;    
    if( SM.getMHl() >= 100.0 && SM.getMHl() <= UpperBoundForApproximateFormulae && !bDebug ) {
        // applicable for 100 GeV <= mHl <= 1 TeV
        Mw0 = 80.3800;
        c1 = 0.05253;
        c2 = 0.010345;
        c3 = 0.001021;
        c4 = -0.000070;
        c5 = 1.077;
        c6 = 0.5270;
        c7 = 0.0698;
        c8 = 0.004055;
        c9 = 0.000110;
        c10 = 0.0716;
        c11 = 115.0;
    } else if ( SM.getMHl() >= 10.0 && SM.getMHl() <= 1000.0 ) {
        // applicable for 10 GeV <= mHl <= 1 TeV
        Mw0 = 80.3799;
        c1 = 0.05429;
        c2 = 0.008939;
        c3 = 0.0000890;
        c4 = 0.000161;
        c5 = 1.070;
        c6 = 0.5256;
        c7 = 0.0678;
        c8 = 0.00179;
        c9 = 0.0000659;
        c10 = 0.0737;
        c11 = 114.9;      
    } else {
        std::stringstream out;
        out << SM.getMHl();
        throw std::runtime_error("ApproximateFormulae::Mw(): mh=" + out.str() + " is out of range");
    }
        
    // Inputs have to be varied within their combined 2 sigma region around 
    // their central values (year 2003) adopted below.
    double dH = log(SM.getMHl()/100.0);
    double dh = pow((SM.getMHl()/100.0), 2.0);
    double dt = pow((SM.getMtpole()/174.3), 2.0) - 1.0;
    double dZ = SM.getMz()/91.1875 - 1.0;
    double dalphae = DeltaAlphaL5q_i/0.05907 - 1.0;
    double dalphas = SM.getAlsMz()/0.119 - 1.0;

    return (Mw0 - c1*dH - c2*dH*dH + c3*pow(dH, 4.0)
            + c4*(dh - 1.0) - c5*dalphae + c6*dt - c7*dt*dt
            - c8*dH*dt + c9*dh*dt - c10*dalphas + c11*dZ
            + SM.getDelMw());
}
    
 
double EWSMApproximateFormulae::sin2thetaEff_l(const StandardModel::lepton l, 
                                               const double DeltaAlphaL5q_i) const 
{
    // applicable for 10 GeV <= mHl <= 1 TeV
    if( SM.getMHl() < 10.0 || SM.getMHl() > UpperBoundForApproximateFormulae ) {
        std::stringstream out;
        out << SM.getMHl();
        throw std::runtime_error("ApproximateFormulae::sin2thetaEff_l(): mh=" + out.str() + " is out of range");
    }
    
    double s0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10;
    switch (l) {
        case StandardModel::NEUTRINO_1: 
        case StandardModel::NEUTRINO_2: 
        case StandardModel::NEUTRINO_3: 
            s0 = 0.2308772;
            d1 = 4.713*0.0001;
            d2 = 2.05*0.00001;
            d3 = 3.85*0.000001;
            d4 = -1.85*0.000001;
            d5 = 2.06*0.01;
            d6 = -2.850*0.001;
            d7 = 1.82*0.0001;
            d8 = -9.71*0.000001;
            d9 = 3.96*0.0001;
            d10 = -6.54*0.1;
            break;
        case StandardModel::ELECTRON:
        case StandardModel::MU:
        case StandardModel::TAU:
            s0 = 0.2312527;
            d1 = 4.729*0.0001;
            d2 = 2.07*0.00001;
            d3 = 3.85*0.000001;
            d4 = -1.85*0.000001;
            d5 = 2.07*0.01;
            d6 = -2.851*0.001;
            d7 = 1.82*0.0001;
            d8 = -9.74*0.000001;
            d9 = 3.98*0.0001;
            d10 = -6.55*0.1;
            break;
        default:
            throw std::runtime_error("Error in ApproximateFormulae::sin2thetaEff_l()"); 
    }

    double L_H = log(SM.getMHl()/100.0);
    double Delta_H = SM.getMHl()/100.0;
    double Delta_ale = DeltaAlphaL5q_i/0.05907 - 1.0;
    double Delta_t = pow((SM.getMtpole()/178.0), 2.0) - 1.0;
    double Delta_alphas = SM.getAlsMz()/0.117 - 1.0;
    double Delta_Z = SM.getMz()/91.1876 - 1.0;
    
    return (s0 + d1*L_H + d2*L_H*L_H + d3*pow(L_H, 4.0)
            + d4*(Delta_H*Delta_H - 1.0) + d5*Delta_ale + d6*Delta_t
            + d7*Delta_t*Delta_t + d8*Delta_t*(Delta_H - 1.0)
            + d9*Delta_alphas + d10*Delta_Z
            + SM.getDelSin2th_l());
}

 
double EWSMApproximateFormulae::sin2thetaEff_q(const StandardModel::quark q, 
                                               const double DeltaAlphaL5q_i) const
{
    // applicable for 10 GeV <= mHl <= 1 TeV
    if( SM.getMHl() < 10.0 || SM.getMHl() > UpperBoundForApproximateFormulae ) {
        std::stringstream out;
        out << SM.getMHl();
        throw std::runtime_error("ApproximateFormulae::sin2thetaEff_q(): mh=" + out.str() + " is out of range");
    }
    
    double s0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10;
    switch (q) {
        case StandardModel::UP: 
        case StandardModel::CHARM:
            s0 = 0.2311395;
            d1 = 4.726*0.0001;
            d2 = 2.07*0.00001;
            d3 = 3.85*0.000001;
            d4 = -1.85*0.000001;
            d5 = 2.07*0.01;
            d6 = -2.853*0.001;
            d7 = 1.83*0.0001;
            d8 = -9.73*0.000001;
            d9 = 3.98*0.0001;
            d10 = -6.55*0.1;
            break;
        case StandardModel::DOWN: 
        case StandardModel::STRANGE:        
            s0 = 0.2310286;
            d1 = 4.720*0.0001;
            d2 = 2.06*0.00001;
            d3 = 3.85*0.000001;
            d4 = -1.85*0.000001;
            d5 = 2.07*0.01;
            d6 = -2.848*0.001;
            d7 = 1.81*0.0001;
            d8 = -9.73*0.000001;
            d9 = 3.97*0.0001;
            d10 = -6.55*0.1;        
            break;
        case StandardModel::BOTTOM:
            s0 = 0.2327580;
            d1 = 4.749*0.0001;
            d2 = 2.03*0.00001;
            d3 = 3.94*0.000001;
            d4 = -1.84*0.000001;
            d5 = 2.08*0.01;
            d6 = -0.993*0.001;
            d7 = 0.708*0.0001;
            d8 = -7.61*0.000001;
            d9 = 4.03*0.0001;
            d10 = 6.61*0.1;
            break;
        case StandardModel::TOP:
            return 0.0;
        default:
            throw std::runtime_error("Error in ApproximateFormulae::sin2thetaEff_q()"); 
    }

    double L_H = log(SM.getMHl()/100.0);
    double Delta_H = SM.getMHl()/100.0;
    double Delta_ale = DeltaAlphaL5q_i/0.05907 - 1.0;
    double Delta_t = pow((SM.getMtpole()/178.0), 2.0) - 1.0;
    double Delta_alphas = SM.getAlsMz()/0.117 - 1.0;
    double Delta_Z = SM.getMz()/91.1876 - 1.0;

    return (s0 + d1*L_H + d2*L_H*L_H + d3*pow(L_H, 4.0)
            + d4*(Delta_H*Delta_H - 1.0) + d5*Delta_ale + d6*Delta_t
            + d7*Delta_t*Delta_t + d8*Delta_t*(Delta_H - 1.0)
            + d9*Delta_alphas + d10*Delta_Z );
}


double EWSMApproximateFormulae::DeltaR_TwoLoopEW_rem(const double DeltaAlphaL5q_i,
                                                     const double Mw_i) const
{
    // applicable for 10 GeV <= mHl <= 1 TeV
    if( SM.getMHl() < 10.0 || SM.getMHl() > UpperBoundForApproximateFormulae ) {
        std::stringstream out;
        out << SM.getMHl();
        throw std::runtime_error("ApproximateFormulae::DeltaR_TwoLoopEW_rem(): mh=" + out.str() + " is out of range");
    }
    
    double r0 =  0.003354;
    double r1 = -0.000209;
    double r2 =  0.0000254;
    double r3 = -0.00000785;
    double r4 = -0.00000233;
    double r5 =  0.00783;
    double r6 =  0.00338;
    double r7 = -0.00000989;
    double r8 =  0.0939;
    double r9 =  0.204;
    double r10= -0.103;

    //double Mw = Mw(DeltaAlphaL5q_i); /* for test */
    double Mw = Mw_i;

    double L_H = log(SM.getMHl()/100.0);
    double Delta_H = SM.getMHl()/100.0;
    //double Delta_ale = DeltaAlphaL5q_i/0.05907 - 1.0;
    double Delta_t = pow((SM.getMtpole()/178.0), 2.0) - 1.0;
    //double Delta_alphas = SM.getAlsMz()/0.117 - 1.0;
    double Delta_Z = SM.getMz()/91.1876 - 1.0;
    double Delta_W = Mw/80.404 - 1.0;
    
    return ( r0 + r1*L_H + r2*L_H*L_H + r3*pow(L_H, 4.0)
             + r4*(Delta_H*Delta_H - 1.0) + r5*Delta_t 
             + r6*Delta_t*Delta_t + r7*Delta_t*L_H + r8*Delta_W
             + r9*Delta_W*Delta_t + r10*Delta_Z );
}


double EWSMApproximateFormulae::DeltaKappa_l_TwoLoopEW_rem(const double DeltaAlphaL5q_i,
                                                           const double Mw_i) const
{
    // applicable for 10 GeV <= mHl <= 1 TeV
    if( SM.getMHl() < 10.0 || SM.getMHl() > UpperBoundForApproximateFormulae ) {
        std::stringstream out;
        out << SM.getMHl();
        throw std::runtime_error("ApproximateFormulae::DeltaKappa_l_TwoLoopEW_rem(): mh=" + out.str() + " is out of range");
    }
    
    double k0 = -0.002711;
    double k1 = -0.0000312;
    double k2 = -0.0000412;
    double k3 =  0.00000528;
    double k4 =  0.00000375;
    double k5 = -0.00516;
    double k6 = -0.00206;
    double k7 = -0.000232;
    double k8 = -0.0647;
    double k9 = -0.129;
    double k10=  0.0712;
    
    double L_H = log(SM.getMHl()/100.0);
    double Delta_H = SM.getMHl()/100.0;
    //double Delta_ale = DeltaAlphaL5q_i/0.05907 - 1.0;
    double Delta_t = pow((SM.getMtpole()/178.0), 2.0) - 1.0;
    //double Delta_alphas = SM.getAlsMz()/0.117 - 1.0;
    double Delta_Z = SM.getMz()/91.1876 - 1.0;
    double Delta_W = Mw_i/80.404 - 1.0;
   
    return ( k0 + k1*L_H + k2*L_H*L_H + k3*pow(L_H, 4.0)
             + k4*(Delta_H*Delta_H - 1.0) + k5*Delta_t 
             + k6*Delta_t*Delta_t + k7*Delta_t*L_H 
             + k8*Delta_W + k9*Delta_W*Delta_t + k10*Delta_Z );
}


double EWSMApproximateFormulae::DeltaKappa_b_TwoLoopEW_rem(const double DeltaAlphaL5q_i,
                                                           const double Mw_i) const
{
    // applicable for 10 GeV <= mHl <= 1 TeV
    if( SM.getMHl() < 10.0 || SM.getMHl() > UpperBoundForApproximateFormulae ) {
        std::stringstream out;
        out << SM.getMHl();
        throw std::runtime_error("ApproximateFormulae::DeltaKappa_b_TwoLoopEW_rem(): mh=" + out.str() + " is out of range");
    }
    
    double k0 = -0.002666;
    double k1 = -0.0000592;
    double k2 = -0.00000329;
    double k3 =  0.00000349;
    double k4 =  0.00000283;
    double k5 = -0.00534;
    double k6 = -0.00210;
    double k7 = -0.000219;
    double k8 = -0.0631;
    double k9 = -0.126;
    double k10=  0.0647;
    
    double L_H = log(SM.getMHl()/100.0);
    double Delta_H = SM.getMHl()/100.0;
    //double Delta_ale = DeltaAlphaL5q_i/0.05907 - 1.0;
    double Delta_t = pow((SM.getMtpole()/178.0), 2.0) - 1.0;
    //double Delta_alphas = SM.getAlsMz()/0.117 - 1.0;
    double Delta_Z = SM.getMz()/91.1876 - 1.0;
    double Delta_W = Mw_i/80.404 - 1.0;    
   
    return ( k0 + k1*L_H + k2*L_H*L_H + k3*pow(L_H, 4.0)
             + k4*(Delta_H*Delta_H - 1.0) + k5*Delta_t 
             + k6*Delta_t*Delta_t + k7*Delta_t*L_H 
             + k8*Delta_W + k9*Delta_W*Delta_t + k10*Delta_Z );
}


double EWSMApproximateFormulae::R0_bottom(const double DeltaAlphaL5q_i) const 
{
    // applicable for 10 GeV <= mHl <= 1 TeV
    if( SM.getMHl() < 10.0 || SM.getMHl() > UpperBoundForApproximateFormulae ) {
        std::stringstream out;
        out << SM.getMHl();
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
    
    double L_H = log(SM.getMHl()/100.0);
    double Delta_H = SM.getMHl()/100.0;
    double Delta_ale = DeltaAlphaL5q_i/0.05900 - 1.0;
    double Delta_t = pow((SM.getMtpole()/173.2), 2.0) - 1.0;
    double Delta_alphas = SM.getAlsMz()/0.1184 - 1.0;
    double Delta_Z = SM.getMz()/91.1876 - 1.0;

    return (Rb00 + c1*L_H + c2*L_H*L_H + c3*pow(L_H, 4.0)
            + c4*(Delta_H*Delta_H - 1.0) + c5*Delta_ale + c6*Delta_t
            + c7*Delta_t*L_H + c8*Delta_alphas + c9*Delta_Z );    
     */

    /*-----------------------------------------------*/
    /* arXiv:1205.0299v2 by Freitas and Huang */

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

    double L_H = log(SM.getMHl()/100.0);
    double Delta_H = SM.getMHl()/100.0;
    double Delta_ale = DeltaAlphaL5q_i/0.05900 - 1.0;
    double Delta_t = pow((SM.getMtpole()/173.2), 2.0) - 1.0;
    double Delta_alphas = SM.getAlsMz()/0.1184 - 1.0;
    double Delta_Z = SM.getMz()/91.1876 - 1.0;    
    
    return (Rb00 + c1*L_H + c2*L_H*L_H + c3*pow(L_H, 4.0)
            + c4*(Delta_H*Delta_H - 1.0) + c5*Delta_ale + c6*Delta_t
            + c7*Delta_t*L_H + c8*Delta_alphas + c9*Delta_alphas*Delta_alphas 
            + c10*Delta_Z );
}


double EWSMApproximateFormulae::Gu_over_Gb(const double DeltaAlphaL5q_i) const 
{
    // applicable for 10 GeV <= mHl <= 1 TeV
    if( SM.getMHl() < 10.0 || SM.getMHl() > UpperBoundForApproximateFormulae ) {
        std::stringstream out;
        out << SM.getMHl();
        throw std::runtime_error("ApproximateFormulae::Gu_over_Gb(): mh=" + out.str() + " is out of range");
    }

    // obtained from Prof. Freitas on Apr. 23, 2013
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
    
    double LH = log(SM.getMHl()/100.0);
    double DH = pow((SM.getMHl()/100.0), 2.0) - 1.0;
    double Dal = DeltaAlphaL5q_i/0.059 - 1.0;
    double Dt = pow((SM.getMtpole()/173.2), 2.0) - 1.0;
    double Das = SM.getAlsMz()/0.1184 - 1.0;
    double Dmz = SM.getMz()/91.1876 - 1.0;
    
    return ( R + c1*LH + c2*LH*LH + c3*pow(LH, 4.0) + c4*DH + c5*Dal + c6*LH*Dal
             + c7*Dt + c8*Dt*Dt + c9*LH*Dt + c10*Das + c11*Das*Das + c12*Dmz );
}
    

double EWSMApproximateFormulae::Gd_over_Gb(const double DeltaAlphaL5q_i) const 
{
    // applicable for 10 GeV <= mHl <= 1 TeV
    if( SM.getMHl() < 10.0 || SM.getMHl() > UpperBoundForApproximateFormulae ) {
        std::stringstream out;
        out << SM.getMHl();
        throw std::runtime_error("ApproximateFormulae::Gd_over_Gb(): mh=" + out.str() + " is out of range");
    }

    // obtained from Prof. Freitas on Apr. 23, 2013
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
   
    double LH = log(SM.getMHl()/100.0);
    double DH = pow((SM.getMHl()/100.0), 2.0) - 1.0;
    double Dal = DeltaAlphaL5q_i/0.059 - 1.0;
    double Dt = pow((SM.getMtpole()/173.2), 2.0) - 1.0;
    double Das = SM.getAlsMz()/0.1184 - 1.0;
    double Dmz = SM.getMz()/91.1876 - 1.0;
    
    return ( R + c1*LH + c2*LH*LH + c3*pow(LH, 4.0) + c4*DH + c5*Dal + c6*LH*Dal
             + c7*Dt + c8*Dt*Dt + c9*LH*Dt + c10*Das + c11*Das*Das + c12*Dmz );
}



