/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Acharm.h"


double Acharm::computeThValue() 
{ 
//    double A_c = SM.A_q(SM.CHARM);
//
//    /* NP contribution to the Zff vertex */
//    if (SM.checkNPZff_linearized() && SM.ModelName().compare("StandardModel") != 0)
//        A_c = SM.getMyEW_NPZff().Acharm(A_c);
    
    return SM.A_f(SM.getQuarks(SM.CHARM));
}
        


