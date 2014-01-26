/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <EWSM.h>
#include <NPZbbbar.h>
#include "Rbottom.h"
#include "sigmaHadron.h"


double Rbottom::computeThValue() 
{ 
    double R0_b;
    if (!SM.IsFlagNoApproximateRb() && SM.ModelName() != "NPEpsilons")
        R0_b = SM.getEWSM()->getMyApproximateFormulae()->X_extended("R0_bottom");
    else
        R0_b = myEW.Gamma_q(SM.BOTTOM)/myEW.Gamma_had();

    /* NP contribution to the Zff vertex */
    if (myEW.checkNPZff())
        R0_b = myEW.getMyEW_NPZff().Rbottom(R0_b);

    /* Debug: extract pure NP contribution */
    //R0_b -= myEW.Gamma_q(SM.BOTTOM)/myEW.Gamma_had();
    
    return R0_b;
}
        

