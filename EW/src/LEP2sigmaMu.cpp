/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "LEP2sigmaMu.h"


double LEP2sigmaMu::computeThValue() 
{ 
    double sigma_mu = SM.LEP2sigmaMu(s);
    
    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_xseellLEP2");
    
    #ifdef LEP2TEST
    sigma_mu = myTEST.sigmaMuTEST(sqrt_s)/SM.GeVminus2_to_nb/1000.0;
    #endif 
        
//    return ( sigma_mu*SM.GeVminus2_to_nb*1000.0 );
    //  Approximate formulae returns result in pb
    return ( sigma_mu * ( 1.0 + SM_intr_err ) );
}



double LEP2dsigmadcosMu::computeThValue() 
{ 
    double dsigmadcos_mu = SM.LEP2dsigmadcosMu(s,cos);
    
    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_deelldcosLEP2");
    
    #ifdef LEP2TEST
    sigma_mu = myTEST.sigmaMuTEST(sqrt_s)/SM.GeVminus2_to_nb/1000.0;
    #endif 

    //  Approximate formulae returns result in pb
    return ( dsigmadcos_mu * ( 1.0 + SM_intr_err ) );
}


double eeffsigmaMu::computeThValue() 
{ 
    double sigma_mu = SM.eeffsigmaMu(pol_e, pol_p, s);
    
    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_xseell");
        
//    return ( sigma_mu*SM.GeVminus2_to_nb*1000.0 );
    //  Approximate formulae returns result in pb
    return ( sigma_mu * ( 1.0 + SM_intr_err ) );
}
        

double eeffRmuon::computeThValue()
{ 
    double R_muon = SM.eeffRmuon(pol_e, pol_p, s);
    
    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_Reell");
    
    return ( R_muon * ( 1.0 + SM_intr_err ) );
}
     