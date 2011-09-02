/* 
 * File:   ThreeLoopQCD.cpp
 * Author: mishima
 */

#include "ThreeLoopQCD.h"


ThreeLoopQCD::ThreeLoopQCD(const EWSMcommon& EWSMC_i) : EWSMcommon(EWSMC_i) {
}

//ThreeLoopQCD::ThreeLoopQCD(const ThreeLoopQCD& orig) {
//}

ThreeLoopQCD::~ThreeLoopQCD() {
}


////////////////////////////////////////////////////////////////////////

double ThreeLoopQCD::DeltaAlpha_l() const {
    return (0.0);
}    

double ThreeLoopQCD::DeltaAlpha_t() const {   
    double xt = pow(Mz/getQuarks(TOP).getMass(), 2.0);
    double log_t = 2.0*logMZtoMTOP;

    double tmp = ( (28.220 + 9.702*log_t) 
                   + xt*(6.924 + 1.594*log_t) )*AlsMz/M_PI*AlsMz/M_PI;
    tmp *= -4.0/45.0*ale/M_PI*xt;
    return tmp;
}

double ThreeLoopQCD::DeltaRho() const {
    /* !! Write codes !!*/
    return (0.0);     
}

double ThreeLoopQCD::DeltaR_rem() const {
    /* !! Write codes !!*/
    return (0.0);     
}

complex ThreeLoopQCD::deltaRho_rem_l(const StandardModel::lepton l) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;      
}

complex ThreeLoopQCD::deltaRho_rem_q(const StandardModel::quark q) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;      
}

complex ThreeLoopQCD::deltaKappa_rem_l(const StandardModel::lepton l) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;      
}

complex ThreeLoopQCD::deltaKappa_rem_q(const StandardModel::quark q) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;      
}










