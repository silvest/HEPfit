/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ZFAFBbottomLEP2.h"


double ZFAFBbottomLEP2::computeThValue() { 
    double myZFAFBbottomLEP2, dummy;
    myZF.calcXS_AFB(9, sqrt_s, &dummy, &myZFAFBbottomLEP2);
    return myZFAFBbottomLEP2;
}


