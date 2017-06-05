/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "LEP2sigmaCharm.h"


double LEP2sigmaCharm::computeThValue() 
{ 
    double sigma_charm = SM.LEP2sigmaCharm(s);
    
    return ( sigma_charm*SM.GeVminus2_to_nb*1000.0 );
}
        
