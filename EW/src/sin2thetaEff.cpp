/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <NPZbbbar.h>
#include "sin2thetaEff.h"


double sin2thetaEff::computeThValue() 
{ 
    double sin2_theta_eff;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        sin2_theta_eff = myEW.getMyEW_CHMN().sin2thetaEff();
    else if (myEWTYPE==EW::EWABC) 
        sin2_theta_eff = myEW.getMyEW_ABC().sin2thetaEff(false);
    else if (myEWTYPE==EW::EWABC2)
        sin2_theta_eff = myEW.getMyEW_ABC().sin2thetaEff(true);
    else { 
        sin2_theta_eff = myEW.sin2thetaEff(SM.ELECTRON);
    
        if (myEWTYPE==EW::EWBURGESS)
            return myEW.getMyEW_BURGESS().sin2thetaEff(sin2_theta_eff);

        /* NP contribution to the Zff vertex */
        if (myEW.checkLEP1NP())
            sin2_theta_eff = myEW.getMyEW_NPZff().sin2thetaEff(sin2_theta_eff);
        
        /* Debug: extract pure NP contribution */
        //sin2_theta_eff -= myEW.sin2thetaEff(SM.ELECTRON);
    }

    return sin2_theta_eff;
}

