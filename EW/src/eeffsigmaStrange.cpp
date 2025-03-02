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
    
    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_xseejj");
    
    return ( sigma_strange * ( 1.0 + SM_intr_err ) );
}


double eeffRstrange::computeThValue()
{ 
    double R_strange = SM.eeffRstrange(pol_e, pol_p, s);

    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_Reejj");
    
    return ( R_strange * ( 1.0 + SM_intr_err ) );
}
        
