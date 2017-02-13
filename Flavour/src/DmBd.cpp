/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "DmBd.h"
#include "StandardModel.h"
 
double  DmBd::computeThValue() 
{
    return(2. * SM.getCBd() * AmpBd(FULLNLO).abs());
}