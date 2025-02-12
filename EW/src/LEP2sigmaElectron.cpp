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
    
    #ifdef LEP2TEST
    sigma_mu = myTEST.sigmaMuTEST(sqrt_s)/SM.GeVminus2_to_nb/1000.0;
    #endif 

    //  Approximate formulae returns result in pb
    return ( dsigmadcos_e );
}
        

double eeffsigmaElectron::computeThValue()
{ 
    double sigma_e = SM.eeffsigmaE(pol_e, pol_p, s);
    
    return ( sigma_e );
}


double eeffRelectron::computeThValue()
{ 
    double R_electron = SM.eeffRelectron(pol_e, pol_p, s);
    
    return ( R_electron );
}


double eeffsigmatsubElectron::computeThValue()
{ 
    double sigma_e = SM.eeffsigmaEtsub(pol_e, pol_p, s);
    
    return ( sigma_e );
}


double eeffRelectrontsub::computeThValue()
{ 
    double R_electron = SM.eeffRelectrontsub(pol_e, pol_p, s);
    
    return ( R_electron );
}
     