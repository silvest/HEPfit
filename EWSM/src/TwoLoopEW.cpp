/* 
 * File:   TwoLoopEW.cpp
 * Author: mishima
 */

#include "TwoLoopEW.h"
#include "TwoLoopQCD.h"
#include "ApproximateFormulae.h"


TwoLoopEW::TwoLoopEW(const EWSMcommon& EWSMC_i) : EWSMC(EWSMC_i) {
}

//TwoLoopEW::TwoLoopEW(const TwoLoopEW& orig) {
//}

TwoLoopEW::~TwoLoopEW() {
}


////////////////////////////////////////////////////////////////////////

double TwoLoopEW::DeltaAlpha_l() const {
    double xl[3] = { pow(EWSMC.GetSM().getMz()
                     /EWSMC.GetSM().getLeptons(EWSMC.GetSM().ELECTRON).getMass(), 2.0), 
                     pow(EWSMC.GetSM().getMz()
                     /EWSMC.GetSM().getLeptons(EWSMC.GetSM().MU).getMass(), 2.0), 
                     pow(EWSMC.GetSM().getMz()
                     /EWSMC.GetSM().getLeptons(EWSMC.GetSM().TAU).getMass(), 2.0) };
    double log_l[3] = { 2.0*EWSMC.GetLogMZtoME(), 
                        2.0*EWSMC.GetLogMZtoMMU(), 
                        2.0*EWSMC.GetLogMZtoMTAU() };

    double twoLoop[3];
    for (int i = 0; i < 3; i++) {
        twoLoop[i] = - 5.0/24.0 + EWSMC.GetZeta3() + log_l[i]/4.0 
                     + 3.0/xl[i]*log_l[i];
    }
            
    return ( pow(EWSMC.GetSM().getAle()/M_PI, 2.0)
             *(twoLoop[0] + twoLoop[1] + twoLoop[2]) );
}    

double TwoLoopEW::DeltaAlpha_t() const {   
    return (0.0);
}

double TwoLoopEW::DeltaRho() const {

    /* added O(alpha^2) contribution from the Z-gamma mixing */
    //double Mw2 = pow(EWSMC.GetMw(), 2.0);
    //SigmaZZ += - Mw2*pow(PiZgamma_fer(mu,Mz2), 2.0)*EWSMC.GetSM().getAle()/4.0/M_PI;
    
    
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
    if(q==StandardModel::TOP) return ( complex(0.0,0.0,false) );
   
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
    if(q==StandardModel::TOP) return ( complex(0.0,0.0,false) );
    
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;      
}


////////////////////////////////////////////////////////////////////////   `

complex TwoLoopEW::rho_2() const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;    
}

complex TwoLoopEW::tau_2() const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;  
}







