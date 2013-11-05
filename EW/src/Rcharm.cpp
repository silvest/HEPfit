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
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        R0_c = myEW.getMyEW_CHMN().R_c();
    else if (myEWTYPE==EW::EWABC || myEWTYPE==EW::EWABC2) 
        R0_c = myEW.getMyEW_ABC().R_c(SM.epsilon1(),SM.epsilon3(),SM.epsilonb());
    else {    
        if (SM.IsFlagApproximateGqOverGb() 
                //&& !SM.IsFlagRhoZbFromGuOverGb()
                //&& !SM.IsFlagRhoZbFromGdOverGb()
                //&& !SM.IsFlagTestSubleadingTwoLoopEW()
                && SM.ModelName() != "NPEpsilons"
                ) {
            double Gu_over_Gb = SM.getEWSM()->Gu_over_Gb_SM();
            double Gd_over_Gb = SM.getEWSM()->Gd_over_Gb_SM();
            R0_c = Gu_over_Gb/(1.0 + 2.0*(Gd_over_Gb + Gu_over_Gb));
        } else
            R0_c = myEW.Gamma_q(SM.CHARM)/myEW.Gamma_had();
        
        if(myEWTYPE==EW::EWBURGESS)
            return myEW.getMyEW_BURGESS().Rcharm(R0_c, myEW.Gamma_had());

        /* NP contribution to the Zff vertex */
         if (SM.ModelName().compare("NPZbbbar") == 0) {
            if (!(static_cast<const NPZbbbar*> (&SM))->IsFlagNotLinearizedNP())
                R0_c = myEW.getMyEW_NPZff().Rcharm(R0_c);
        } else
            R0_c = myEW.getMyEW_NPZff().Rcharm(R0_c);
        
        /* Debug: extract pure NP contribution */
        //R0_c -= myEW.Gamma_q(SM.CHARM)/myEW.Gamma_had();
    }

    return R0_c;
}
        

