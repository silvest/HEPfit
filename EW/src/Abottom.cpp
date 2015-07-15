/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Abottom.h"

double Abottom::computeThValue()
{
    return SM.A_f(SM.getQuarks(SM.BOTTOM));
}


