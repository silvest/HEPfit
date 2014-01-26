/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <EWSM.h>
#include <NPZbbbar.h>
#include "Rcharm.h"


double Rcharm::computeThValue() 
{   
    double R0_c;
    if (!SM.IsFlagNoApproximateGqOverGb() && SM.ModelName() != "NPEpsilons")
        R0_c = SM.getEWSM()->getMyApproximateFormulae()->X_extended("R0_charm", SM.DeltaAlphaL5q());
    else
        R0_c = myEW.Gamma_q(SM.CHARM)/myEW.Gamma_had();

    /* NP contribution to the Zff vertex */
    if (myEW.checkNPZff())
        R0_c = myEW.getMyEW_NPZff().Rcharm(R0_c);

    /* Debug: extract pure NP contribution */
    //R0_c -= myEW.Gamma_q(SM.CHARM)/myEW.Gamma_had();

    return R0_c;
}
        

