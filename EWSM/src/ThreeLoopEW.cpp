/* 
 * File:   ThreeLoopEW.cpp
 * Author: mishima
 */

#include "ThreeLoopEW.h"


ThreeLoopEW::ThreeLoopEW(const EWSMcommon& EWSMC_i) : EWSMC(EWSMC_i) {
}

//ThreeLoopEW::ThreeLoopEW(const ThreeLoopEW& orig) {
//}

ThreeLoopEW::~ThreeLoopEW() {
}


////////////////////////////////////////////////////////////////////////

double ThreeLoopEW::DeltaAlpha_l() const {
    double xl[3] = { pow(EWSMC.GetSM().getMz()
                     /EWSMC.GetSM().getLeptons(EWSMC.GetSM().ELECTRON).getMass(), 2.0), 
                     pow(EWSMC.GetSM().getMz()
                     /EWSMC.GetSM().getLeptons(EWSMC.GetSM().MU).getMass(), 2.0), 
                     pow(EWSMC.GetSM().getMz()
                     /EWSMC.GetSM().getLeptons(EWSMC.GetSM().TAU).getMass(), 2.0) };
    double log_l[3] = { 2.0*EWSMC.GetLogMZtoME(), 
                        2.0*EWSMC.GetLogMZtoMMU(), 
                        2.0*EWSMC.GetLogMZtoMTAU() };    

    double threeLoop[3];
    for (int i = 0; i < 3; i++) {
        threeLoop[i] = - 121.0/48.0 + (-5.0 + 8.0*EWSMC.GetLog2())*EWSMC.GetZeta2() 
                       - 99.0/16.0*EWSMC.GetZeta3() + 10.0*EWSMC.GetZeta5()
                       + log_l[i]/8.0;
        for (int j = 0; j < 3; j++) {
            if (i > j) { /* Pi^{(2)}_l */
                threeLoop[i] += - 116.0/27.0 + 4.0/3.0*EWSMC.GetZeta2() 
                                + 38.0/9.0*EWSMC.GetZeta3() + 14.0/9.0*log_l[i]
                                + (5.0/18.0 - 4.0/3.0*EWSMC.GetZeta3())*log_l[j]
                                + log_l[i]*log_l[i]/6.0
                                - log_l[i]*log_l[j]/3.0;
            } else if (i == j) { /* Pi^{(2)}_F */
                threeLoop[i] += - 307.0/216.0 - 8.0/3.0*EWSMC.GetZeta2() 
                                + 545.0/144.0*EWSMC.GetZeta3()
                                + (11.0/6.0 - 4.0/3.0*EWSMC.GetZeta3())*log_l[i]
                                - log_l[i]*log_l[i]/6.0;
            } else { /* Pi^{(2)}_h */
                threeLoop[i] += - 37.0/6.0 + 38.0/9.0*EWSMC.GetZeta3()
                                + (11.0/6.0 - 4.0/3.0*EWSMC.GetZeta3())*log_l[j]
                                - log_l[j]*log_l[j]/6.0;
            }
        }
        threeLoop[i] /= -4.0;
    }
            
    return ( pow(EWSMC.GetSM().getAle()/M_PI, 3.0)
             *(threeLoop[0] + threeLoop[1] + threeLoop[2]) );    
}    

double ThreeLoopEW::DeltaAlpha_t() const {   
    return (0.0);
}

double ThreeLoopEW::DeltaRho() const {
    /* !! Write codes !!*/
    return (0.0);    
}

double ThreeLoopEW::DeltaR_rem() const {
    /* !! Write codes !!*/
    return (0.0);    
}

complex ThreeLoopEW::deltaRho_rem_l(const StandardModel::lepton l) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;      
}

complex ThreeLoopEW::deltaRho_rem_q(const StandardModel::quark q) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;      
}

complex ThreeLoopEW::deltaKappa_rem_l(const StandardModel::lepton l) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;      
}

complex ThreeLoopEW::deltaKappa_rem_q(const StandardModel::quark q) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;      
}










