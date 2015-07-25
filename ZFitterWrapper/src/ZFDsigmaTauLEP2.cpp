/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ZFDsigmaTauLEP2.h"


double ZFDsigmaTauLEP2::computeThValue() { 
    double myZFDsigmaTauLEP2;
    myZF.calcDXS(3, sqrt_s, cos_theta, &myZFDsigmaTauLEP2);
    myZFDsigmaTauLEP2 *= 1000.0;// nb --> pb
    return myZFDsigmaTauLEP2;
}


