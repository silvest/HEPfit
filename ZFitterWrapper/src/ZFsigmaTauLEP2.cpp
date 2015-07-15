/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ZFsigmaTauLEP2.h"


double ZFsigmaTauLEP2::computeThValue() { 
    double myZFsigmaTauLEP2, dummy;
    myZF.calcXS_AFB(3, sqrt_s, &myZFsigmaTauLEP2, &dummy);
    myZFsigmaTauLEP2 *= 1000.0;// nb --> pb
    return myZFsigmaTauLEP2;
}


