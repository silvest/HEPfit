/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ZFMh0.h"


double ZFMh0::computeThValue() {   
    return myZF.getModel().getMHl();
}

