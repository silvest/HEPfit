/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Abottom.h"


double Abottom::computeThValue() 
{
//    double A_b = SM.A_q(SM.BOTTOM);
//
//    /* NP contribution to the Zff vertex */
//    if (SM.checkNPZff_linearized() && SM.ModelName().compare("StandardModel") != 0)
//        A_b = SM.getMyEW_NPZff().Abottom(A_b);
    
    return SM.A_f(SM.getQuarks(SM.BOTTOM));
}
        

