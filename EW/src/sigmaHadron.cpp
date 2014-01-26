/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
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
    if (!SM.IsFlagNoApproximateSigmaH()  && SM.ModelName() != "NPEpsilons")
        sigma_had = SM.getEWSM()->getMyApproximateFormulae()->X_extended("sigmaHadron", SM.DeltaAlphaL5q())/GeVminus2_to_nb;
    else
        sigma_had = myEW.sigma0_had();

    /* NP contribution to the Zff vertex */
    if (myEW.checkNPZff())
        sigma_had = myEW.getMyEW_NPZff().sigmaHadron(sigma_had);

    /* Debug: extract pure NP contribution */
    //sigma_had -= myEW.sigma0_had();
    
    return ( sigma_had*GeVminus2_to_nb );
}
        


