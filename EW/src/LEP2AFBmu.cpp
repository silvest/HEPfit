/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "LEP2AFBmu.h"


double LEP2AFBmu::computeThValue() 
{ 
    double AFB_mu = SM.LEP2AFBmu(s);
  
    #ifdef LEP2TEST
    AFB_mu = myTEST.AFBmuTEST(sqrt_s);
    #endif

    return AFB_mu;
}
        
double eeffAFBmu::computeThValue()
{ 
    double AFB_mu = SM.eeffAFBmu(pol_e, pol_p, s);
       
    return AFB_mu;
}