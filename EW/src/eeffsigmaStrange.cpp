/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */


#include "eeffsigmaStrange.h"


double eeffsigmaStrange::computeThValue()
{ 
    double sigma_strange = SM.eeffsigmaStrange(pol_e, pol_p, s);
    
    return ( sigma_strange );
}
        
