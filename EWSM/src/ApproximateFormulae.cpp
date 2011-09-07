/* 
 * File:   ApproximateFormulae.cpp
 * Author: mishima
 */

#include "ApproximateFormulae.h"


ApproximateFormulae::ApproximateFormulae(const StandardModel& SM_i, 
                                         const double DeltaAlpha_i) : SM(SM_i) {
    myDeltaAlpha = DeltaAlpha_i;
}

//ApproximateFormulae::ApproximateFormulae(const ApproximateFormulae& orig) {
//}

ApproximateFormulae::~ApproximateFormulae() {
}


////////////////////////////////////////////////////////////////////////

double ApproximateFormulae::Mw() const {
    double Mw0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11;    
    if( SM.getMHl() >= 100.0 && SM.getMHl() <= 1000.0 ) {
#define TEST_DEBUG
#ifndef TEST_DEBUG
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
    } else if (SM.getMHl() >= 10.0 && SM.getMHl() < 100.0 ) {        
#endif
        
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
        throw "Higgs mass is out of range in ApproximateFormulae::Mw()";        
    }

    // Inputs have to be varied within their combined 2 sigma region around 
    // their central values (year 2003) adopted below.
    double dH = log(SM.getMHl()/100.0);
    double dh = pow((SM.getMHl()/100.0), 2.0);
    double dt = pow((SM.getQuarks(SM.TOP).getMass()/174.3), 2.0) - 1.0;
    double dZ = SM.getMz()/91.1875 - 1.0;
    double dalphae = myDeltaAlpha/0.05907 - 1.0;
    double dalphas = SM.getAlsMz()/0.119 - 1.0;

    return (Mw0 - c1*dH - c2*dH*dH + c3*pow(dH, 4.0)
            + c4*(dh - 1.0) - c5*dalphae + c6*dt - c7*dt*dt
            - c8*dH*dt + c9*dh*dt - c10*dalphas + c11*dZ);
}
    
 
double ApproximateFormulae::sin2thetaEff(const StandardModel::lepton l) const {
    // applicable for 10 GeV <= mHl <= 1 TeV
    if( SM.getMHl() < 100.0 || SM.getMHl() > 1000.0 ) {
        throw "Higgs mass is out of range in ApproximateFormulae::sin2thetaEff()";        
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
            throw "Error in ApproximateFormulae::sin2thetaEff()";
    }

    double L_H = log(SM.getMHl()/100.0);
    double Delta_H = SM.getMHl()/100.0;
    double Delta_ale = myDeltaAlpha/0.05907 - 1.0;
    double Delta_t = pow((SM.getQuarks(SM.TOP).getMass()/178.0), 2.0) - 1.0;
    double Delta_alphas = SM.getAlsMz()/0.117 - 1.0;
    double Delta_Z = SM.getMz()/91.1876 - 1.0;

    return (s0 + d1*L_H + d2*L_H*L_H + d3*pow(L_H, 4.0)
            + d4*(Delta_H*Delta_H - 1.0) + d5*Delta_ale + d6*Delta_t
            + d7*Delta_t*Delta_t + d8*Delta_t*(Delta_H - 1.0)
            + d9*Delta_alphas + d10*Delta_Z );
}
 
double ApproximateFormulae::sin2thetaEff(const StandardModel::quark q) const {
    // applicable for 10 GeV <= mHl <= 1 TeV
    if( SM.getMHl() < 100.0 || SM.getMHl() > 1000.0 ) {
        throw "Higgs mass is out of range in ApproximateFormulae::sin2thetaEff()";        
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
            throw "Error in ApproximateFormulae::sin2thetaEff()";
    }

    double L_H = log(SM.getMHl()/100.0);
    double Delta_H = SM.getMHl()/100.0;
    double Delta_ale = myDeltaAlpha/0.05907 - 1.0;
    double Delta_t = pow((SM.getQuarks(SM.TOP).getMass()/178.0), 2.0) - 1.0;
    double Delta_alphas = SM.getAlsMz()/0.117 - 1.0;
    double Delta_Z = SM.getMz()/91.1876 - 1.0;

    return (s0 + d1*L_H + d2*L_H*L_H + d3*pow(L_H, 4.0)
            + d4*(Delta_H*Delta_H - 1.0) + d5*Delta_ale + d6*Delta_t
            + d7*Delta_t*Delta_t + d8*Delta_t*(Delta_H - 1.0)
            + d9*Delta_alphas + d10*Delta_Z );
}

double ApproximateFormulae::DeltaR_TwoLoopEW() const {
    throw "Write codes!!";
}

double ApproximateFormulae::DeltaKappa_l_TwoLoopEW() const {
    throw "Write codes!!";    
}

double ApproximateFormulae::DeltaKappa_b_TwoLoopEW() const {
    throw "Write codes!!";
}


