/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ZFRbottom.h"


double ZFRbottom::getThValue() {
    return myZF.Gamma_f(9)/myZF.Gamma_had();
}

