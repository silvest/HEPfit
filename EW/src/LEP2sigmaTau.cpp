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
    
    return ( sigma_tau*SM.GeVminus2_to_nb*1000.0 );
}
        

