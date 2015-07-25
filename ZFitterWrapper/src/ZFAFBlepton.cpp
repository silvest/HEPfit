/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ZFAFBlepton.h"


double ZFAFBlepton::computeThValue() {
    return ( 3.0/4.0*myZF.Af(1)*myZF.Af(1) );
}

