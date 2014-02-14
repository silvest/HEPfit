/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Abottom.h"


double Abottom::computeThValue() 
{
    double A_b = myEW.A_q(SM.BOTTOM);

    /* NP contribution to the Zff vertex */
    if (myEW.checkNPZff_linearized() && SM.ModelName().compare("StandardModel") != 0)
        A_b = myEW.getMyEW_NPZff().Abottom(A_b);
    
    return A_b;
}
        

