/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "LEP2Rbottom.h"


double LEP2Rbottom::computeThValue() 
{ 
    
    double R_b = SM.LEP2Rbottom(s);

    #ifdef LEP2TEST
    R_b = myTEST.RbottomTEST(sqrt_s);
    #endif

    return R_b;
}





