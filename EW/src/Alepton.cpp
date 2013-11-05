/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <NPZbbbar.h>
#include "Alepton.h"


double Alepton::computeThValue()
{
    double A_l;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)
        A_l = myEW.getMyEW_CHMN().A_l(SM.ELECTRON);
    else if (myEWTYPE==EW::EWABC)
        A_l = myEW.getMyEW_ABC().A_l(SM.ELECTRON,SM.epsilon1(),SM.epsilon3(),false);
    else if (myEWTYPE==EW::EWABC2)
        A_l = myEW.getMyEW_ABC().A_l(SM.ELECTRON,SM.epsilon1(),SM.epsilon3(),true);
    else {
        A_l = myEW.A_l(SM.ELECTRON);

        if(myEWTYPE==EW::EWBURGESS)
            return myEW.getMyEW_BURGESS().Alepton(A_l);

        /* NP contribution to the Zff vertex */
        if (SM.ModelName().compare("NPZbbbar") == 0) {
            if (!(static_cast<const NPZbbbar*> (&SM))->IsFlagNotLinearizedNP())
                A_l = myEW.getMyEW_NPZff().Alepton(A_l);
        } else
            A_l = myEW.getMyEW_NPZff().Alepton(A_l);

        /* Debug: extract pure NP contribution */
        //A_l -= myEW.A_l(SM.ELECTRON);
    }

    return A_l;
}


