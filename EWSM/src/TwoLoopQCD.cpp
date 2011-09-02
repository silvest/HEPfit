/* 
 * File:   TwoLoopQCD.cpp
 * Author: mishima
 */

#include "TwoLoopQCD.h"


TwoLoopQCD::TwoLoopQCD(const EWSMcommon& EWSMC_i) : EWSMcommon(EWSMC_i) {
}

//TwoLoopQCD::TwoLoopQCD(const TwoLoopQCD& orig) {
//}

TwoLoopQCD::~TwoLoopQCD() {
}


////////////////////////////////////////////////////////////////////////

double TwoLoopQCD::DeltaAlpha_l() const {
    return (0.0);
}    

double TwoLoopQCD::DeltaAlpha_t() const {   
    double xt = pow(Mz/getQuarks(TOP).getMass(), 2.0);
    double tmp = (5.062 + xt*0.8315)*AlsMz/M_PI;
    tmp *= -4.0/45.0*ale/M_PI*xt;
    return tmp;
}

double TwoLoopQCD::DeltaRho() const {
    /* !! Write codes !!*/
    return (0.0);     
}

double TwoLoopQCD::DeltaR_rem() const {
    /* !! Write codes !!*/
    return (0.0);     
}

complex TwoLoopQCD::deltaRho_rem_l(const StandardModel::lepton l) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;       
}

complex TwoLoopQCD::deltaRho_rem_q(const StandardModel::quark q) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;       
}

complex TwoLoopQCD::deltaKappa_rem_l(const StandardModel::lepton l) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;       
}

complex TwoLoopQCD::deltaKappa_rem_q(const StandardModel::quark q) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;       
}








