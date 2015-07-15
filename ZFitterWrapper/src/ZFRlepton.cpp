/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ZFRlepton.h"


double ZFRlepton::computeThValue() {
    return myZF.Gamma_had()/myZF.Gamma_f(1);
}

