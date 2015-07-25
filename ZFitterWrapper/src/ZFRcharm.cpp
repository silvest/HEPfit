/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ZFRcharm.h"


double ZFRcharm::computeThValue() {
    return myZF.Gamma_f(6)/myZF.Gamma_had();
}


