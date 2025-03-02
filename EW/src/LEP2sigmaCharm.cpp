/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "LEP2sigmaCharm.h"


double LEP2sigmaCharm::computeThValue() 
{ 
    double sigma_charm = SM.LEP2sigmaCharm(s);
    
    return ( sigma_charm*SM.GeVminus2_to_nb*1000.0 );
}

double eeffsigmaCharm::computeThValue() 
{ 
    double sigma_charm = SM.eeffsigmaCharm(pol_e, pol_p, s);
    
    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_xseejj");
    
    return ( sigma_charm * ( 1.0 + SM_intr_err ) );
}
        

double eeffRcharm::computeThValue()
{ 
    double R_charm = SM.eeffRcharm(pol_e, pol_p, s);
    
    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_Reejj");
    
    return ( R_charm * ( 1.0 + SM_intr_err ) );
}
     