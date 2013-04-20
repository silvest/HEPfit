/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include "epsilonb.h"


double epsilonb::getThValue() 
{  
//    if (SM.IsFlagR0bApproximate() && !SM.IsFlagRhoZbFromR0b())
//        throw std::runtime_error("epsilonb::getThValue() cannot be used!");
//    else {
        double epsb = SM.epsilonb();
        
        //if ( myEW.checkModelForSTU() ) {
        //    throw std::runtime_error("Write codes in epsilonb::getThValue()");   
        //}
        
        return epsb;
  //  }
}


