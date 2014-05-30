/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "sin2thetaEff.h"


double sin2thetaEff::computeThValue() 
{ 
    double sin2_theta_eff = SM.sin2thetaEff(SM.getLeptons(SM.ELECTRON));
    
    return sin2_theta_eff;
}

