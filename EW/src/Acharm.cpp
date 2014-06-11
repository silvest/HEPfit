/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Acharm.h"

double Acharm::computeThValue()
{
    return SM.A_f(SM.getQuarks(SM.CHARM));
}



