/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "DmBs.h"

double  DmBs::computeThValue() 
{
    return(2. * AmpBs(FULLNLO).abs());
}