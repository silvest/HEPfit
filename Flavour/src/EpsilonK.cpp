/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EpsilonK.h"

double EpsilonK::getThValue(){
    
    std::cout << "EpsK = " << 1./SM.getDeltaMK() * (AmpDK(NLO).imag()) * SM.getKbarEpsK() / sqrt(2.) << std::endl;
    
    return(1./SM.getDeltaMK() * (AmpDK(NLO).imag()) * SM.getKbarEpsK() / sqrt(2.));
}
