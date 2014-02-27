/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "sigmaHadron.h"


double sigmaHadron::computeThValue() 
{ 
    double sigma_had = myEW.sigma0_had();

    /* NP contribution to the Zff vertex */
    if (myEW.checkNPZff_linearized() && SM.ModelName().compare("StandardModel") != 0)
        sigma_had = myEW.getMyEW_NPZff().sigmaHadron(sigma_had);
    
    return ( sigma_had*GeVminus2_to_nb );
}
        


