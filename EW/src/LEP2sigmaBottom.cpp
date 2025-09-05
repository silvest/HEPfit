/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */


#include "LEP2sigmaBottom.h"


double LEP2sigmaBottom::computeThValue()
{ 
    double sigma_bottom = SM.LEP2sigmaBottom(s);
    
    return ( sigma_bottom*SM.GeVminus2_to_nb*1000.0 );
}


double eeffsigmaBottom::computeThValue()
{ 
    //double sigma_bottom = SM.eeffsigmaBottom(pol_e, pol_p, s);
    double coscut = 0.9;
    double sigma_bottom = SM.eeffsigma(SM.getQuarks(SM.BOTTOM), pol_e, pol_p, s, -coscut, coscut);
    
    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_xseebb");
    
    return ( sigma_bottom * ( 1.0 + SM_intr_err ) );
}


double eeffRbottom::computeThValue()
{ 
    double R_bottom = SM.eeffRbottom(pol_e, pol_p, s);
    
    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_Reebb");
    
    return ( R_bottom * ( 1.0 + SM_intr_err ) );
}
        
