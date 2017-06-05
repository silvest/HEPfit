/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "LEP2AFBcharm.h"


double LEP2AFBcharm::computeThValue() 
{ 
    double AFB_c = SM.LEP2AFBcharm(s);

    #ifdef LEP2TEST
    AFB_c = myTEST.AFBcharmTEST(sqrt_s);
    #endif


    return AFB_c;
}
        


