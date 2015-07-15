/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "DmK.h"

double DmK::computeThValue() 
{
    return(2.*AmpMK(FULLNLO).real() + SM.getDmk());
}
