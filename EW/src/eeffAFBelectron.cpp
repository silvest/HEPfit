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
    
    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_AFBeeee");
       
    return ( AFB_e * ( 1.0 + SM_intr_err ) );
}


double eeffAFBelectrontsub::computeThValue()
{ 
    double AFB_e = SM.eeffAFBetsub(pol_e, pol_p, s);
    
    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_AFBeeee");
       
    return ( AFB_e * ( 1.0 + SM_intr_err ) );
}