/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <EWSM.h>
#include <NPZbbbar.h>
#include "GammaZ.h"


double GammaZ::computeThValue() 
{ 
    double Gamma_Z;
    if (!SM.IsFlagNoApproximateGammaZ() && SM.ModelName() != "NPEpsilons")
        Gamma_Z = SM.getEWSM()->getMyApproximateFormulae()->X_extended("GammaZ");
    else
        Gamma_Z = myEW.Gamma_Z();

    /* Theoretical uncertainty */
    Gamma_Z += SM.getDelGammaZ();

    /* NP contribution to the Zff vertex */
    if (myEW.checkNPZff())
        Gamma_Z = myEW.getMyEW_NPZff().GammaZ(Gamma_Z);

    /* Debug: extract pure NP contribution */
    //Gamma_Z -= myEW.Gamma_Z();
      
    return Gamma_Z;
}
        
