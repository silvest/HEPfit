/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ZFsigmaMuLEP2.h"


double ZFsigmaMuLEP2::computeThValue() { 
    double myZFsigmaMuLEP2, dummy;
    myZF.calcXS_AFB(2, sqrt_s, &myZFsigmaMuLEP2, &dummy);
    myZFsigmaMuLEP2 *= 1000.0;// nb --> pb
    return myZFsigmaMuLEP2;
}


