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
    double A_b = myEW.A_q(SM.BOTTOM);

    /* NP contribution to the Zff vertex */
    if (myEW.checkNPZff())
        A_b = myEW.getMyEW_NPZff().Abottom(A_b);

    /* Debug: extract pure NP contribution */
    //A_b -= myEW.A_q(SM.BOTTOM);
    
    return A_b;
}
        

