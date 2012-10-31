/* 
 * File:   epsilonb.cpp
 * Author: mishima
 */

#include <stdexcept>
#include "epsilonb.h"


double epsilonb::getThValue() {  
    double epsb = SM.epsilonb();
    
    if ( myEW.checkModelForSTU() ) {
        throw std::runtime_error("Write codes in epsilonb::getThValue()");   
    }
    
    return epsb;
}


