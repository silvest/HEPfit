/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <NPZbbbar.h>
#include "Rlepton.h"


double Rlepton::computeThValue() 
{
    double R0_l = myEW.Gamma_had()/myEW.Gamma_l(SM.ELECTRON);

    /* NP contribution to the Zff vertex */
    if (myEW.checkLEP1NP())
        R0_l = myEW.getMyEW_NPZff().Rlepton(R0_l);

    /* Debug: extract pure NP contribution */
    //R0_l -= myEW.Gamma_had()/myEW.Gamma_l(SM.ELECTRON);
 
    return R0_l;
}
        

