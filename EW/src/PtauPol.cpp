/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "PtauPol.h"


double PtauPol::computeThValue() 
{
    double P_tau_pol = SM.A_f(SM.getLeptons(SM.TAU));

    return P_tau_pol;
}
        
