/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "LEP2AFBbottom.h"


double LEP2AFBbottom::computeThValue()
{ 
    double AFB_b = SM.LEP2AFBbottom(s);

    #ifdef LEP2TEST
    AFB_b = myTEST.AFBbottomTEST(sqrt_s);
    #endif
       
    return AFB_b;
}
        
