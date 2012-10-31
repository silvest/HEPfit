/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ZFRcharm.h"


double ZFRcharm::getThValue() {
    return myZF.Gamma_f(6)/myZF.Gamma_had();
}


