/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <NPZbbbar.h>
#include "Acharm.h"


double Acharm::computeThValue() 
{ 
    double A_c;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        A_c = myEW.getMyEW_CHMN().A_q(SM.CHARM);
    else if (myEWTYPE==EW::EWABC || myEWTYPE==EW::EWABC2) 
        A_c = myEW.getMyEW_ABC().A_q(SM.CHARM,SM.epsilon1(),SM.epsilon3());
    else {
        A_c = myEW.A_q(SM.CHARM);

        if(myEWTYPE==EW::EWBURGESS)
            return myEW.getMyEW_BURGESS().Acharm(A_c, myEW.A_l(SM.ELECTRON));
        
        /* NP contribution to the Zff vertex */
        if (SM.ModelName().compare("NPZbbbar") == 0) {
            if (!(static_cast<const NPZbbbar*> (&SM))->IsFlagNotLinearizedNP())
                A_c = myEW.getMyEW_NPZff().Acharm(A_c);
        } else
            A_c = myEW.getMyEW_NPZff().Acharm(A_c);

        /* Debug: extract pure NP contribution */
        //A_c -= myEW.A_q(SM.CHARM);
    }
    
    return A_c;
}
        


