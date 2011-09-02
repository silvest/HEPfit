/* 
 * File:   ThreeLoopEW2QCD.cpp
 * Author: mishima
 */

#include "ThreeLoopEW2QCD.h"


ThreeLoopEW2QCD::ThreeLoopEW2QCD(const EWSMcommon& EWSMC_i) : EWSMcommon(EWSMC_i) {
}

//ThreeLoopEW2QCD::ThreeLoopEW2QCD(const ThreeLoopEW2QCD& orig) {
//}

ThreeLoopEW2QCD::~ThreeLoopEW2QCD() {
}


////////////////////////////////////////////////////////////////////////

double ThreeLoopEW2QCD::DeltaAlpha_l() const {
    return (0.0);
}    

double ThreeLoopEW2QCD::DeltaAlpha_t() const {   
    return (0.0);
}

double ThreeLoopEW2QCD::DeltaRho() const {
    /* !! Write codes !!*/
    return (0.0);     
}

double ThreeLoopEW2QCD::DeltaR_rem() const {
    /* !! Write codes !!*/
    return (0.0);     
}

complex ThreeLoopEW2QCD::deltaRho_rem_l(const StandardModel::lepton l) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;      
}

complex ThreeLoopEW2QCD::deltaRho_rem_q(const StandardModel::quark q) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;      
}

complex ThreeLoopEW2QCD::deltaKappa_rem_l(const StandardModel::lepton l) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;      
}

complex ThreeLoopEW2QCD::deltaKappa_rem_q(const StandardModel::quark q) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;      
}









