/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */


#include "LEP2sigmaTau.h"


double LEP2sigmaTau::computeThValue() 
{ 
    double sigma_tau = SM.LEP2sigmaTau(s);
    
    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_xseellLEP2");
    
    #ifdef LEP2TEST
    sigma_tau = myTEST.sigmaTauTEST(sqrt_s)/SM.GeVminus2_to_nb/1000.0;
    #endif 
    
//    return ( sigma_tau*SM.GeVminus2_to_nb*1000.0 );
    
    //  Approximate formulae returns result in pb
    return ( sigma_tau * ( 1.0 + SM_intr_err ) );
}
        


double LEP2dsigmadcosTau::computeThValue() 
{ 
    double dsigmadcos_tau = SM.LEP2dsigmadcosBinTau(s, cos, cosmin, cosmax);
    
    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_deelldcosLEP2");
    
    #ifdef LEP2TEST
    sigma_tau = myTEST.sigmaTauTEST(sqrt_s)/SM.GeVminus2_to_nb/1000.0;
    #endif 
    
    //  Approximate formulae returns result in pb
    return ( dsigmadcos_tau * ( 1.0 + SM_intr_err ) );
}


double eeffsigmaTau::computeThValue() 
{ 
    //double sigma_tau = SM.eeffsigmaTau(pol_e, pol_p, s);
    double coscut = 0.9;
    double sigma_tau = SM.eeffsigma(SM.getLeptons(SM.TAU), pol_e, pol_p, s, -coscut, coscut);
    
    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_xseell");
        
//    return ( sigma_mu*SM.GeVminus2_to_nb*1000.0 );
    //  Approximate formulae returns result in pb
    return ( sigma_tau * ( 1.0 + SM_intr_err ) );
}
  

double eeffRtau::computeThValue()
{ 
    double R_tau = SM.eeffRtau(pol_e, pol_p, s);
    
    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_Reell");
    
    return ( R_tau * ( 1.0 + SM_intr_err ) );
}
     