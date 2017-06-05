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
        
    return ( sigma_mu*SM.GeVminus2_to_nb*1000.0 );
}
        

