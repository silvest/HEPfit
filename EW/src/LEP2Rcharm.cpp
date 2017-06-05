/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "LEP2Rcharm.h"


double LEP2Rcharm::computeThValue() 
{ 
    
    double R_c = SM.LEP2Rcharm(s);

    #ifdef LEP2TEST
    R_c = myTEST.RcharmTEST(sqrt_s);
    #endif

    return R_c;
}


