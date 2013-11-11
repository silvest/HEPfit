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
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        sigma_had = myEW.getMyEW_CHMN().sigma0_had();
    else if (myEWTYPE==EW::EWABC) 
        sigma_had = myEW.getMyEW_ABC().sigma0_had(SM.epsilon1(),SM.epsilon3(),SM.epsilonb(),false);
    else if (myEWTYPE==EW::EWABC2)
        sigma_had = myEW.getMyEW_ABC().sigma0_had(SM.epsilon1(),SM.epsilon3(),SM.epsilonb(),true)/GeVminus2_to_nb;
    else {   
        if (SM.IsFlagApproximateSigmaH())
            sigma_had = SM.getEWSM()->sigmaHadron_SM()/GeVminus2_to_nb;
        else
            sigma_had = myEW.sigma0_had();
        
        if (myEWTYPE==EW::EWBURGESS) {
            sigma_had = myEW.getMyEW_BURGESS().sigmaHadron(sigma_had,
                    myEW.Gamma_Z(), myEW.Gamma_had(), myEW.Gamma_l(SM.ELECTRON));
            return ( sigma_had*GeVminus2_to_nb );
        }

        /* NP contribution to the Zff vertex */
        if (myEW.checkLEP1NP())
            sigma_had = myEW.getMyEW_NPZff().sigmaHadron(sigma_had);
        
        /* Debug: extract pure NP contribution */
        //sigma_had -= myEW.sigma0_had();
    }
    
    return ( sigma_had*GeVminus2_to_nb );
}
        


