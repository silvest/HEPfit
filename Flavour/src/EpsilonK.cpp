/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EpsilonK.h"

double EpsilonK::getThValue(){
    return(1./SM.getDeltaMK() * (AmpDK(NLO).imag()) * SM.getKbarEpsK() / sqrt(2.));
}
