/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
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


