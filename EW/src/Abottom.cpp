/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <NPZbbbar.h>
#include "Abottom.h"


double Abottom::computeThValue() 
{
    double A_b;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        A_b = myEW.getMyEW_CHMN().A_q(SM.BOTTOM);
    else if (myEWTYPE==EW::EWABC || myEWTYPE==EW::EWABC2) 
        A_b = myEW.getMyEW_ABC().A_b(SM.epsilon1(),SM.epsilon3(),SM.epsilonb());
    else {
        A_b = myEW.A_q(SM.BOTTOM);

        if(myEWTYPE==EW::EWBURGESS)
            return myEW.getMyEW_BURGESS().Abottom(A_b, myEW.A_l(SM.ELECTRON));

        /* NP contribution to the Zff vertex */
        if (SM.ModelName().compare("NPZbbbar") == 0) {
            if (!(static_cast<const NPZbbbar*> (&SM))->IsFlagNotLinearizedNP())
                A_b = myEW.getMyEW_NPZff().Abottom(A_b);
        } else
            A_b = myEW.getMyEW_NPZff().Abottom(A_b);

        /* Debug: extract pure NP contribution */
        //A_b -= myEW.A_q(SM.BOTTOM);
    }
    
    return A_b;
}
        

