/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ZFRlepton.h"


double ZFRlepton::getThValue() {
    return myZF.Gamma_had()/myZF.Gamma_f(1);
}

