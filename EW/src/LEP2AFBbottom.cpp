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


double eeffAFBbottom::computeThValue()
{ 
    double AFB_b = SM.eeffAFBbottom(pol_e, pol_p, s);
    
    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_AFBeebb");
       
    return ( AFB_b * ( 1.0 + SM_intr_err ) );
}
        
