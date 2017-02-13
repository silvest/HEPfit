/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "DmBs.h"
#include "StandardModel.h"

double  DmBs::computeThValue() 
{
    return(2. * SM.getCBs() * AmpBs(FULLNLO).abs());
}