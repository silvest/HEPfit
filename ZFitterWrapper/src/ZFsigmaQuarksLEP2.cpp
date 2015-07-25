/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ZFsigmaQuarksLEP2.h"

double ZFsigmaQuarksLEP2::computeThValue() { 
    double myZFsigmaQuarksLEP2 = 0.0, tmp, dummy;
    for (int i=4; i<=9; i++) {
        if (i!=8) {
            myZF.calcXS_AFB(i, sqrt_s, &tmp, &dummy);
            myZFsigmaQuarksLEP2 += tmp;
        }
    }
    myZFsigmaQuarksLEP2 *= 1000.0;// nb --> pb
    return myZFsigmaQuarksLEP2;
}


