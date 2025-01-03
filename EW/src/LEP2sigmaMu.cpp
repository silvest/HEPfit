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
    
    #ifdef LEP2TEST
    sigma_mu = myTEST.sigmaMuTEST(sqrt_s)/SM.GeVminus2_to_nb/1000.0;
    #endif 
        
//    return ( sigma_mu*SM.GeVminus2_to_nb*1000.0 );
    //  Approximate formulae returns result in pb
    return ( sigma_mu );
}



double LEP2dsigmadcosMu::computeThValue() 
{ 
    double dsigmadcos_mu = SM.LEP2dsigmadcosMu(s,cos);
    
    #ifdef LEP2TEST
    sigma_mu = myTEST.sigmaMuTEST(sqrt_s)/SM.GeVminus2_to_nb/1000.0;
    #endif 

    //  Approximate formulae returns result in pb
    return ( dsigmadcos_mu );
}


double eeffsigmaMu::computeThValue() 
{ 
    double sigma_mu = SM.eeffsigmaMu(pol_e, pol_p, s);
        
//    return ( sigma_mu*SM.GeVminus2_to_nb*1000.0 );
    //  Approximate formulae returns result in pb
    return ( sigma_mu );
}
        

