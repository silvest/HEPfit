/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "sin2thetaEff.h"


double sin2thetaEff::computeThValue() 
{ 
    double sin2_theta_eff = myEW.sin2thetaEff(SM.ELECTRON);
    
    /* NP contribution to the Zff vertex */
    if (myEW.checkNPZff_linearized() && SM.ModelName().compare("StandardModel") != 0)
        sin2_theta_eff = myEW.getMyEW_NPZff().sin2thetaEff(sin2_theta_eff);

    return sin2_theta_eff;
}

