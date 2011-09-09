/* 
 * File:   EWSMcommon.cpp
 * Author: mishima
 * 
 * Created on August 31, 2011, 2:30 AM
 */

#include "EWSMcommon.h"
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <gsl/gsl_sf.h>


EWSMcommon::EWSMcommon(const StandardModel& SM_i) : SM(SM_i) {
    SetConstants();
}

//EWSMcommon::EWSMcommon(const EWSMcommon& orig) {
//}

EWSMcommon::~EWSMcommon() {
}


//////////////////////////////////////////////////////////////////////// 

void EWSMcommon::SetConstants() {
 
    /* zeta functions */
    zeta2 = gsl_sf_zeta_int(2);
    zeta3 = gsl_sf_zeta_int(3);
    zeta4 = gsl_sf_zeta_int(4);
    zeta5 = gsl_sf_zeta_int(5);
        
    /* logarithmic functions */
    log2 = log(2.0);
    logMZtoME = log( SM.getMz()/SM.getLeptons(SM.ELECTRON).getMass() );
    logMZtoMMU = log( SM.getMz()/SM.getLeptons(SM.MU).getMass() );    
    logMZtoMTAU = log( SM.getMz()/SM.getLeptons(SM.TAU).getMass() );    
    logMZtoMTOP = log( SM.getMz()/SM.getQuarks(SM.TOP).getMass() );

    Xt = SM.getGF()*pow(SM.getQuarks(SM.TOP).getMass(), 2.0)
         /8.0/sqrt(2.0)/M_PI/M_PI;
    
    /* TEST */
    AlsMt = 0.1074432788;

    double Cl2_Pi_3 = gsl_sf_clausen(M_PI/3.0); /* Clausen function */
    S2 = 4.0/9.0/sqrt(3)*Cl2_Pi_3;
    D3 = 6.0*zeta3 - 15.0/4.0*zeta4 - 6.0*Cl2_Pi_3*Cl2_Pi_3;
    //double Li4_1_2 = ;
    //B4 = 16.0*Li4_1_2 - 4.0*zeta2*log2*log2 + 2.0/3.0*pow(log2,4.0) - 13.0/2.0*zeta4;
    B4 = - 1.76280008707377;
    //int prec_def = std::cout.precision();
    //std::cout << std::setprecision(15) 
    //          << "TEST: S2=" << S2 << " D3=" << D3 << " B4=" << B4 << std::endl;
    //std::cout.precision(prec_def);

    
    
    
    
}

void EWSMcommon::Compute(const double Mw_i) {
    
    Mw = Mw_i;
    cW2 = Mw*Mw/SM.getMz()/SM.getMz();
    sW2 = 1.0 - cW2;    
    
    f_AlphaToGF = sqrt(2.0)*SM.getGF()*pow(SM.getMz(),2.0)*sW2*cW2/M_PI/SM.getAle();
      
    
}