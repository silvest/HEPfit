/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */


#include "LEP2sigmaHadron.h"


double LEP2sigmaHadron::computeThValue()
{ 
    double sigmaH = SM.LEP2sigmaHadron(s);
    
    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_xseejjLEP2");
    
    #ifdef LEP2TEST
    sigmaH = myTEST.sigmaHadronTEST(sqrt_s)/SM.GeVminus2_to_nb/1000.0;
    #endif
    
//    return ( sigmaH*SM.GeVminus2_to_nb*1000.0 );
    
    //  Approximate formulae returns result in pb    
    return ( sigmaH * ( 1.0 + SM_intr_err ) );
}
        

double eeffsigmaHadron::computeThValue()
{ 
    double sigmaH = SM.eeffsigmaHadron(pol_e, pol_p, s);
    
    // Intrinsic SM theory uncertainty
    double SM_intr_err = SM.getOptionalParameter("errSMint_xseejj");
    
    return ( sigmaH * ( 1.0 + SM_intr_err ) );
}

