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
    
    #ifdef LEP2TEST
    sigma_tau = myTEST.sigmaTauTEST(sqrt_s)/SM.GeVminus2_to_nb/1000.0;
    #endif 
    
//    return ( sigma_tau*SM.GeVminus2_to_nb*1000.0 );
    
    //  Approximate formulae returns result in pb
    return ( sigma_tau );
}
        


double LEP2dsigmadcosTau::computeThValue() 
{ 
    double dsigmadcos_tau = SM.LEP2dsigmadcosTau(s,cos);
    
    #ifdef LEP2TEST
    sigma_tau = myTEST.sigmaTauTEST(sqrt_s)/SM.GeVminus2_to_nb/1000.0;
    #endif 
    
    //  Approximate formulae returns result in pb
    return ( dsigmadcos_tau );
}


double eeffsigmaTau::computeThValue() 
{ 
    double sigma_tau = SM.eeffsigmaTau(pol_e, pol_p, s);
        
//    return ( sigma_mu*SM.GeVminus2_to_nb*1000.0 );
    //  Approximate formulae returns result in pb
    return ( sigma_tau );
}
  

double eeffRtau::computeThValue()
{ 
    double R_tau = SM.eeffRtau(pol_e, pol_p, s);
    
    return ( R_tau );
}
     