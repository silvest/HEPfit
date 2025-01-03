/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "LEP2AFBtau.h"


double LEP2AFBtau::computeThValue() 
{ 

    double AFB_tau = SM.LEP2AFBtau(s);

    #ifdef LEP2TEST
    AFB_tau = myTEST.AFBtauTEST(sqrt_s);
    #endif


    return AFB_tau;
}
        
double eeffAFBtau::computeThValue()
{ 
    double AFB_tau = SM.eeffAFBtau(pol_e, pol_p, s);
       
    return AFB_tau;
}