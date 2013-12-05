/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <NPZbbbar.h>
#include "PtauPol.h"
#include "EW_BURGESS.h"


double PtauPol::computeThValue() 
{  
    double P_tau_pol = myEW.A_l(SM.TAU);

    /* NP contribution to the Zff vertex */
    if (myEW.checkLEP1NP())
        P_tau_pol = myEW.getMyEW_NPZff().PtauPol(P_tau_pol);

    /* Debug: extract pure NP contribution */
    //P_tau_pol -= myEW.A_l(SM.TAU);
 
    return P_tau_pol;
}
        
