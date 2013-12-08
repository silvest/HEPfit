/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <EWSM.h>
#include <NPZbbbar.h>
#include "sigmaHadron.h"


double sigmaHadron::computeThValue() 
{ 
    double sigma_had;
    if (SM.IsFlagNoApproximateSigmaH())
        sigma_had = myEW.sigma0_had();
    else
        sigma_had = SM.getEWSM()->sigmaHadron_SM()/GeVminus2_to_nb;

    /* NP contribution to the Zff vertex */
    if (myEW.checkNPZff())
        sigma_had = myEW.getMyEW_NPZff().sigmaHadron(sigma_had);

    /* Debug: extract pure NP contribution */
    //sigma_had -= myEW.sigma0_had();
    
    return ( sigma_had*GeVminus2_to_nb );
}
        


