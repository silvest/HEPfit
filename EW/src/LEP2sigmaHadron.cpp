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
    
    #ifdef LEP2TEST
    sigmaH = myTEST.sigmaHadronTEST(sqrt_s)/SM.GeVminus2_to_nb/1000.0;
    #endif
    
    return ( sigmaH*SM.GeVminus2_to_nb*1000.0 );
}
        

