/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Rbottom.h"


double Rbottom::computeThValue() 
{ 
    double R0_b = SM.R0_f(SM.getQuarks(SM.BOTTOM));

    return R0_b;
}
        

