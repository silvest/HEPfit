/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
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
    if (SM.IsFlagApproximateGqOverGb()
            //&& !SM.IsFlagRhoZbFromGuOverGb()
            //&& !SM.IsFlagRhoZbFromGdOverGb()
            //&& !SM.IsFlagTestSubleadingTwoLoopEW()
            && SM.ModelName() != "NPEpsilons"
            ) {
        double Gu_over_Gb = SM.getEWSM()->Gu_over_Gb_SM();
        double Gd_over_Gb = SM.getEWSM()->Gd_over_Gb_SM();
        R0_b = 1.0/(1.0 + 2.0*(Gd_over_Gb + Gu_over_Gb));
    } else
        R0_b = myEW.Gamma_q(SM.BOTTOM)/myEW.Gamma_had();

    /* NP contribution to the Zff vertex */
    if (myEW.checkLEP1NP())
        R0_b = myEW.getMyEW_NPZff().Rbottom(R0_b);

    /* Debug: extract pure NP contribution */
    //R0_b -= myEW.Gamma_q(SM.BOTTOM)/myEW.Gamma_had();
    
    return R0_b;
}
        

