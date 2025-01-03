/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "eeffAFBstrange.h"


double eeffAFBstrange::computeThValue()
{ 
    double AFB_s = SM.eeffAFBstrange(pol_e, pol_p, s);
       
    return AFB_s;
}
        
