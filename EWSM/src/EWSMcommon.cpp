/* 
 * File:   EWSMcommon.cpp
 * Author: mishima
 * 
 * Created on August 31, 2011, 2:30 AM
 */

#include "EWSMcommon.h"
#include <cmath>
#include <gsl/gsl_sf_zeta.h>


EWSMcommon::EWSMcommon(const StandardModel& SM_i) : StandardModel(SM_i) {
    SetConstants();
}

//EWSMcommon::EWSMcommon(const EWSMcommon& orig) {
//}

EWSMcommon::~EWSMcommon() {
}


//////////////////////////////////////////////////////////////////////// 

void EWSMcommon::SetConstants() {
 
    double tmp = 4.0*M_PI*ale/sqrt(2.0)/GF/Mz/Mz;
    Mw_tree = Mz/sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - tmp));
    
    /* zeta functions */
    zeta2 = gsl_sf_zeta_int(2);
    zeta3 = gsl_sf_zeta_int(3);
    zeta5 = gsl_sf_zeta_int(5);
    
    /* logarithmic functions */
    log2 = log(2.0);
    logMZtoME = log( Mz/getLeptons(ELECTRON).getMass() );
    logMZtoMMU = log( Mz/getLeptons(MU).getMass() );    
    logMZtoMTAU = log( Mz/getLeptons(TAU).getMass() );    
    logMZtoMTOP = log( Mz/getQuarks(TOP).getMass() );
}

void EWSMcommon::Compute(const double Mw_i) {
    
    Mw = Mw_i;
    cW2 = Mw*Mw/Mz/Mz;
    sW2 = 1.0 - cW2;    
    
    f_AlphaToGF = sqrt(2.0)*GF*Mz*Mz*sW2*cW2/M_PI/ale;
      
    
}