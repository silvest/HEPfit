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
    
    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_AFBeellLEP2");
  
    #ifdef LEP2TEST
    AFB_mu = myTEST.AFBmuTEST(sqrt_s);
    #endif

    return ( AFB_mu * ( 1.0 + SM_intr_err ) );
}
        
double eeffAFBmu::computeThValue()
{ 
    double AFB_mu = SM.eeffAFBmu(pol_e, pol_p, s);
    
    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_AFBeell");
       
    return ( AFB_mu * ( 1.0 + SM_intr_err ) );
}