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
    
    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_AFBeejj");
       
    return ( AFB_s * ( 1.0 + SM_intr_err ) );
}
        
