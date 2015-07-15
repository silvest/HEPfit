/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ZFsin2thetaEff.h"


double ZFsin2thetaEff::computeThValue() {
    return myZF.s2teff_f(1);
}

