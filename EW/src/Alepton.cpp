/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Alepton.h"


double Alepton::computeThValue()
{
//    double A_l = SM.A_l(SM.ELECTRON);
//
//    /* NP contribution to the Zff vertex */
//    if (SM.checkNPZff_linearized() && SM.ModelName().compare("StandardModel") != 0)
//        A_l = SM.getMyEW_NPZff().Alepton(A_l);

    return SM.A_f(SM.getLeptons(SM.ELECTRON));
}


