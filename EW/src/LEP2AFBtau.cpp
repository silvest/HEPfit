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
    
    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_AFBeellLEP2");

    #ifdef LEP2TEST
    AFB_tau = myTEST.AFBtauTEST(sqrt_s);
    #endif


    return ( AFB_tau * ( 1.0 + SM_intr_err ) );
}
        
double eeffAFBtau::computeThValue()
{ 
    double AFB_tau = SM.eeffAFBtau(pol_e, pol_p, s);
    
    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_AFBeell");
       
    return ( AFB_tau * ( 1.0 + SM_intr_err ) );
}