/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <NPZbbbar.h>
#include "Rlepton.h"


double Rlepton::computeThValue() 
{
    double R0_l;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        R0_l = myEW.getMyEW_CHMN().R_l(SM.ELECTRON);
    else if (myEWTYPE==EW::EWABC) 
        R0_l = myEW.getMyEW_ABC().R_l(SM.epsilon1(),SM.epsilon3(),SM.epsilonb(),false);
    else if (myEWTYPE==EW::EWABC2)
        R0_l = myEW.getMyEW_ABC().R_l(SM.epsilon1(),SM.epsilon3(),SM.epsilonb(),true);
    else {       
        R0_l = myEW.Gamma_had()/myEW.Gamma_l(SM.ELECTRON);
        
        if (myEWTYPE==EW::EWBURGESS)
            return myEW.getMyEW_BURGESS().Rlepton(R0_l, myEW.Gamma_had(), myEW.Gamma_l(SM.ELECTRON));

        /* NP contribution to the Zff vertex */
        if (myEW.checkLEP1NP())
            R0_l = myEW.getMyEW_NPZff().Rlepton(R0_l);
        
        /* Debug: extract pure NP contribution */
        //R0_l -= myEW.Gamma_had()/myEW.Gamma_l(SM.ELECTRON);
    }
 
    return R0_l;
}
        

