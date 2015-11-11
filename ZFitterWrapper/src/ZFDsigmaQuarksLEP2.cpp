/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ZFDsigmaQuarksLEP2.h"


double ZFDsigmaQuarksLEP2::computeThValue() { 
    double myZFDsigmaQuarksLEP2 = 0.0, tmp;
    for (int i=4; i<=9; i++) {
        if (i!=8) {
            myZF.calcDXS(i, sqrt_s, cos_theta, &tmp);
            myZFDsigmaQuarksLEP2 += tmp;
        }
    }
    myZFDsigmaQuarksLEP2 *= 1000.0;// nb --> pb
    return myZFDsigmaQuarksLEP2;
}
