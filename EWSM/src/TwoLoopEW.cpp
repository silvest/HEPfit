/* 
 * File:   TwoLoopEW.cpp
 * Author: mishima
 */

#include "TwoLoopEW.h"


TwoLoopEW::TwoLoopEW(const EWSMcommon& EWSMC_i) : EWSMcommon(EWSMC_i) {
}

//TwoLoopEW::TwoLoopEW(const TwoLoopEW& orig) {
//}

TwoLoopEW::~TwoLoopEW() {
}


////////////////////////////////////////////////////////////////////////

double TwoLoopEW::DeltaAlpha_l() const {
    double xl[3] = { pow(Mz/getLeptons(ELECTRON).getMass(), 2.0), 
                     pow(Mz/getLeptons(MU).getMass(), 2.0), 
                     pow(Mz/getLeptons(TAU).getMass(), 2.0) };
    double log_l[3] = { 2.0*logMZtoME, 2.0*logMZtoMMU, 2.0*logMZtoMTAU };    

    double twoLoop[3];
    for (int i = 0; i < 3; i++) {
        twoLoop[i] = - 5.0/24.0 + zeta3 + log_l[i]/4.0 
                     + 3.0/xl[i]*log_l[i];
    }
            
    return ( pow(ale/M_PI, 2.0)*(twoLoop[0] + twoLoop[1] + twoLoop[2]) );
}    

double TwoLoopEW::DeltaAlpha_t() const {   
    return (0.0);
}

double TwoLoopEW::DeltaRho() const {
    /* !! Write codes !!*/
    return (0.0);    
}

double TwoLoopEW::DeltaR_rem() const {
    /* !! Write codes !!*/
    return (0.0);     
}

complex TwoLoopEW::deltaRho_rem_l(const StandardModel::lepton l) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;      
}

complex TwoLoopEW::deltaRho_rem_q(const StandardModel::quark q) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;      
}

complex TwoLoopEW::deltaKappa_rem_l(const StandardModel::lepton l) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;      
}

complex TwoLoopEW::deltaKappa_rem_q(const StandardModel::quark q) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;      
}





