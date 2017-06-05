/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */


#include "LEP2sigmaBottom.h"


double LEP2sigmaBottom::computeThValue()
{ 
    double sigma_bottom = SM.LEP2sigmaBottom(s);
    
    return ( sigma_bottom*SM.GeVminus2_to_nb*1000.0 );
}
        
