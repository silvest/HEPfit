/* 
 * File:   EWSMTwoLoopEW.cpp
 * Author: mishima
 */

#include "EWSMTwoLoopEW.h"


EWSMTwoLoopEW::EWSMTwoLoopEW(const EWSMcache& cache_i) : cache(cache_i) {
}


////////////////////////////////////////////////////////////////////////

double EWSMTwoLoopEW::DeltaAlpha_l() const {
    double xl[3] = { pow(cache.Mz()/cache.mf(StandardModel::ELECTRON), 2.0), 
                     pow(cache.Mz()/cache.mf(StandardModel::MU), 2.0), 
                     pow(cache.Mz()/cache.mf(StandardModel::TAU), 2.0) };
    double log_l[3] = { 2.0*cache.logMZtoME(), 
                        2.0*cache.logMZtoMMU(), 
                        2.0*cache.logMZtoMTAU() };

    double twoLoop[3];
    for (int i = 0; i < 3; i++) {
        twoLoop[i] = - 5.0/24.0 + cache.GetZeta3() + log_l[i]/4.0 
                     + 3.0/xl[i]*log_l[i];
    }
            
    return ( pow(cache.ale()/M_PI, 2.0)
             *(twoLoop[0] + twoLoop[1] + twoLoop[2]) );
}    


double EWSMTwoLoopEW::DeltaAlpha_t() const {   
    return (0.0);
}


double EWSMTwoLoopEW::DeltaRho(const double Mw_i) const {

    /* added O(alpha^2) contribution from the Z-gamma mixing */
    //double Mw2 = pow(EWSMC.GetMw(), 2.0);
    //SigmaZZ += - Mw2*pow(PiZgamma_fer(mu,Mz2), 2.0)*EWSMC.GetSM().getAle()/4.0/M_PI;
    
    
    /* !! Write codes !!*/
    return (0.0);    
}


double EWSMTwoLoopEW::DeltaR_rem(const double Mw_i) const {
    /* !! Write codes !!*/
    return (0.0);     
}


template<typename T> 
complex EWSMTwoLoopEW::deltaRho_rem_f(const T f, const double Mw_i) const {
    cache.checkSMfermion(f, "EWSMTwoLoopEW::deltaRho_rem_f");
    if(f==StandardModel::TOP) return ( complex(0.0,0.0,false) );
   
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;      
}


template<typename T> 
complex EWSMTwoLoopEW::deltaKappa_rem_f(const T f, const double Mw_i) const {
    cache.checkSMfermion(f, "EWSMTwoLoopEW::deltaRho_rem_f");
    if(f==StandardModel::TOP) return ( complex(0.0,0.0,false) );
    
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;      
}


////////////////////////////////////////////////////////////////////////   `

complex EWSMTwoLoopEW::rho_2(const double Mw_i) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;    
}


complex EWSMTwoLoopEW::tau_2(const double Mw_i) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;  
}









