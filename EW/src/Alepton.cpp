/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Alepton.h"

double Alepton::computeThValue()
{
    return SM.A_f(SM.getLeptons(SM.ELECTRON));
}


