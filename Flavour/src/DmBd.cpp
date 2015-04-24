/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "DmBd.h"
 
double  DmBd::computeThValue() 
{
    return(2. * AmpBd(FULLNLO).abs());
}