/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "eeffAFBelectron.h"
        
double eeffAFBelectron::computeThValue()
{ 
    double AFB_e = SM.eeffAFBe(pol_e, pol_p, s);
       
    return AFB_e;
}