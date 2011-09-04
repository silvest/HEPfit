/* 
 * File:   OneLoopEW.cpp
 * Author: mishima
 */

#include "OneLoopEW.h"
#include <PVfunctions.h>


OneLoopEW::OneLoopEW(const EWSMcommon& EWSMC_i) : EWSMcommon(EWSMC_i) {
}

//OneLoopEW::OneLoopEW(const OneLoopEW& orig) {
//}

OneLoopEW::~OneLoopEW() {
}


////////////////////////////////////////////////////////////////////////

double OneLoopEW::DeltaAlpha_l() const {   
    double xl[3] = { pow(Mz/getLeptons(ELECTRON).getMass(), 2.0), 
                     pow(Mz/getLeptons(MU).getMass(), 2.0), 
                     pow(Mz/getLeptons(TAU).getMass(), 2.0) };
    double log_l[3] = { 2.0*logMZtoME, 2.0*logMZtoMMU, 2.0*logMZtoMTAU };    

    double oneLoop[3];
    for (int i = 0; i < 3; i++) {
        oneLoop[i] = - 5.0/9.0 + log_l[i]/3.0 - 2.0/xl[i];
    }
            
    return( ale/M_PI*(oneLoop[0] + oneLoop[1] + oneLoop[2]) );
}

double OneLoopEW::DeltaAlpha_t() const {   
    double xt = pow(Mz/getQuarks(TOP).getMass(), 2.0);
    double tmp = 1.0 + xt*0.1071;
    tmp *= -4.0/45.0*ale/M_PI*xt;
    return tmp;
}

double OneLoopEW::DeltaRho() const {
    /* !! Write codes !!*/
    return (0.0);
}

double OneLoopEW::DeltaR_rem() const {
    /* !! Write codes !!*/
    return (0.0);    
}

double OneLoopEW::DeltaRbar_rem() const {
    /* !! Write codes !!*/
    return (0.0);    
}

complex OneLoopEW::deltaRho_rem_l(const StandardModel::lepton l) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;    
}

complex OneLoopEW::deltaRho_rem_q(const StandardModel::quark q) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;      
}

complex OneLoopEW::deltaKappa_rem_l(const StandardModel::lepton l) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;      
}

complex OneLoopEW::deltaKappa_rem_q(const StandardModel::quark q) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;    
}


//////////////////////////////////////////////////////////////////////// 



