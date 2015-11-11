/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ZFDsigmaMuLEP2.h"


double ZFDsigmaMuLEP2::computeThValue() { 
    double myZFDsigmaMuLEP2;
    myZF.calcDXS(2, sqrt_s, cos_theta, &myZFDsigmaMuLEP2);
    myZFDsigmaMuLEP2 *= 1000.0;// nb --> pb
    return myZFDsigmaMuLEP2;
}


