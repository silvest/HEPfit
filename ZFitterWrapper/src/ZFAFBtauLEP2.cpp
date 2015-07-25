/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ZFAFBtauLEP2.h"


double ZFAFBtauLEP2::computeThValue() { 
    double myZFAFBtauLEP2, dummy;
    myZF.calcXS_AFB(3, sqrt_s, &dummy, &myZFAFBtauLEP2);
    return myZFAFBtauLEP2;
}


