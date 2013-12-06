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
    double A_c = myEW.A_q(SM.CHARM);

    /* NP contribution to the Zff vertex */
    if (myEW.checkNPZff())
        A_c = myEW.getMyEW_NPZff().Acharm(A_c);

    /* Debug: extract pure NP contribution */
    //A_c -= myEW.A_q(SM.CHARM);
    
    return A_c;
}
        


