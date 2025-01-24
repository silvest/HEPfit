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


double eeffsigmaBottom::computeThValue()
{ 
    double sigma_bottom = SM.eeffsigmaBottom(pol_e, pol_p, s);
    
    return ( sigma_bottom );
}


double eeffRbottom::computeThValue()
{ 
    double R_bottom = SM.eeffRbottom(pol_e, pol_p, s);
    
    return ( R_bottom );
}
        
