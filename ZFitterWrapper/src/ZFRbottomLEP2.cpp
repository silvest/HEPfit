/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ZFRbottomLEP2.h"


double ZFRbottomLEP2::computeThValue() { 
    double myZFsigmaBottomLEP2;
    double myZFsigmaQuarksLEP2 = 0.0, tmp, dummy;
    for (int i=4; i<=9; i++) {
        if (i!=8) {
            myZF.calcXS_AFB(i, sqrt_s, &tmp, &dummy);
            myZFsigmaQuarksLEP2 += tmp;
            if (i==9) myZFsigmaBottomLEP2 = tmp;
        }
    }
    return (myZFsigmaBottomLEP2/myZFsigmaQuarksLEP2);
}


