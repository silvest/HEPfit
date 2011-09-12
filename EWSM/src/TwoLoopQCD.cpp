/* 
 * File:   TwoLoopQCD.cpp
 * Author: mishima
 */

#include "TwoLoopQCD.h"


TwoLoopQCD::TwoLoopQCD(const EWSMcommon& EWSMC_i) : EWSMC(EWSMC_i) {
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
    double xt = pow(EWSMC.GetSM().getMz()
                /EWSMC.GetSM().getQuarks(EWSMC.GetSM().TOP).getMass(), 2.0);
    double tmp = (5.062 + xt*0.8315)*EWSMC.GetSM().getAlsMz()/M_PI;
    tmp *= -4.0/45.0*EWSMC.GetSM().getAle()/M_PI*xt;
    return tmp;
} 

double TwoLoopQCD::DeltaRho() const {
    return ( 3.0*EWSMC.GetXt_alpha()*EWSMC.GetAlsMt()/M_PI*deltaQCD_2() );     
}

double TwoLoopQCD::DeltaR_rem() const {
    return ( (2.0*DeltaR_ud() + DeltaR_tb())
             + EWSMC.GetCW2()/EWSMC.GetSW2()/EWSMC.GetF_AlphaToGF()*DeltaRho() );     
}

complex TwoLoopQCD::deltaRho_rem_l(const StandardModel::lepton l) const {
    return ( (2.0*DeltaRho_ud() + DeltaRho_tb()) - DeltaRho() );       
}

complex TwoLoopQCD::deltaRho_rem_q(const StandardModel::quark q) const {
    return ( (2.0*DeltaRho_ud() + DeltaRho_tb()) - DeltaRho() );    
}

complex TwoLoopQCD::deltaKappa_rem_l(const StandardModel::lepton l) const {
    return ( (2.0*DeltaKappa_ud() + DeltaKappa_tb())
             - EWSMC.GetCW2()/EWSMC.GetSW2()*DeltaRho() );  
}

complex TwoLoopQCD::deltaKappa_rem_q(const StandardModel::quark q) const {
    return ( (2.0*DeltaKappa_ud() + DeltaKappa_tb())
             - EWSMC.GetCW2()/EWSMC.GetSW2()*DeltaRho() );    
}


////////////////////////////////////////////////////////////////////////

double TwoLoopQCD::deltaQCD_2() const {
    return ( -2.0/3.0*(1.0+2.0*EWSMC.GetZeta2()) );
}

double TwoLoopQCD::DeltaR_ud() const {
    /* !! Write codes !!*/
    return (0.0);      
}

double TwoLoopQCD::DeltaR_tb() const {
    /* !! Write codes !!*/
    return (0.0);      
}

double TwoLoopQCD::DeltaRho_ud() const {
    /* !! Write codes !!*/
    return (0.0);      
}

complex TwoLoopQCD::DeltaRho_tb() const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;        
}

double TwoLoopQCD::DeltaKappa_ud() const {
    /* !! Write codes !!*/
    return (0.0);      
}

complex TwoLoopQCD::DeltaKappa_tb() const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;       
}




