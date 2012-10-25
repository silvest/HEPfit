/* 
 * File:   epsilonb.cpp
 * Author: mishima
 */

#include <stdexcept>
#include "epsilonb.h"


double epsilonb::getThValue() {  
    double DeltaRhoPrime = 2.0*( sqrt(SM.rhoZ_l(SM.ELECTRON).abs()) - 1.0 );
    double eps1 = DeltaRhoPrime;
    double epsb = - 1.0 + sqrt(SM.rhoZ_q(SM.BOTTOM).abs())/(1.0 + eps1/2.0);
    
    if ( myEW.checkModelForSTU() ) {
        throw std::runtime_error("Write codes in epsilonb::getThValue()");   
    }
    
    return epsb;
}


