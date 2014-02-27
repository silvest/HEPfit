/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "PtauPol.h"


double PtauPol::computeThValue() 
{  
    double P_tau_pol = myEW.A_l(SM.TAU);

    /* NP contribution to the Zff vertex */
    if (myEW.checkNPZff_linearized() && SM.ModelName().compare("StandardModel") != 0)
        P_tau_pol = myEW.getMyEW_NPZff().PtauPol(P_tau_pol);
 
    return P_tau_pol;
}
        
