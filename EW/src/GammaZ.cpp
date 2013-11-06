/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <EWSM.h>
#include <NPZbbbar.h>
#include "GammaZ.h"


double GammaZ::computeThValue() 
{ 
    double Gamma_Z;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        Gamma_Z = myEW.getMyEW_CHMN().GammaZ();
    else if (myEWTYPE==EW::EWABC) 
        Gamma_Z = myEW.getMyEW_ABC().GammaZ(SM.epsilon1(),SM.epsilon3(),SM.epsilonb(),false);
    else if (myEWTYPE==EW::EWABC2)
        Gamma_Z = myEW.getMyEW_ABC().GammaZ(SM.epsilon1(),SM.epsilon3(),SM.epsilonb(),true);
    else {
        if (SM.IsFlagApproximateGammaZ())
            Gamma_Z = SM.getEWSM()->GammaZ_SM();
        else
            Gamma_Z = myEW.Gamma_Z();

        /* Theoretical uncertainty */
        Gamma_Z += SM.getDelGammaZ();

        if (myEWTYPE==EW::EWBURGESS)
            return myEW.getMyEW_BURGESS().GammaZ(Gamma_Z);

        /* NP contribution to the Zff vertex */
        if (myEW.checkLEP1NP())
            Gamma_Z = myEW.getMyEW_NPZff().GammaZ(Gamma_Z);

        /* Debug: extract pure NP contribution */
        //Gamma_Z -= myEW.Gamma_Z();
    }
      
    return Gamma_Z;
}
        
