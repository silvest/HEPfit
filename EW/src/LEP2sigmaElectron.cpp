/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "LEP2sigmaElectron.h"


double LEP2dsigmadcosElectron::computeThValue() 
{ 
    double dsigmadcos_e = SM.LEP2dsigmadcosE(s,cos);
    
    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_deeeedcosLEP2");
    
    #ifdef LEP2TEST
    sigma_mu = myTEST.sigmaMuTEST(sqrt_s)/SM.GeVminus2_to_nb/1000.0;
    #endif 

    //  Approximate formulae returns result in pb
    return ( dsigmadcos_e * ( 1.0 + SM_intr_err ) );
}
        

double eeffsigmaElectron::computeThValue()
{ 
    double sigma_e = SM.eeffsigmaE(pol_e, pol_p, s);
    
    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_xseeee");
    
    return ( sigma_e * ( 1.0 + SM_intr_err ) );
}


double eeffRelectron::computeThValue()
{ 
    double R_electron = SM.eeffRelectron(pol_e, pol_p, s);

    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_Reeee");
    
    return ( R_electron * ( 1.0 + SM_intr_err ) );
}


double eeffsigmatsubElectron::computeThValue()
{ 
    double sigma_e = SM.eeffsigmaEtsub(pol_e, pol_p, s);
    
    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_xseeee");
    
    return ( sigma_e * ( 1.0 + SM_intr_err ) );
}


double eeffRelectrontsub::computeThValue()
{ 
    double R_electron = SM.eeffRelectrontsub(pol_e, pol_p, s);

    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_Reeee");
    
    return ( R_electron * ( 1.0 + SM_intr_err ) );
}
     