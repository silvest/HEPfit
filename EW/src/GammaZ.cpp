/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GammaZ.h"


double GammaZ::computeThValue() 
{ 
    double Gamma_Z = myEW.Gamma_Z();

    /* Theoretical uncertainty */
    Gamma_Z += SM.getDelGammaZ();

    /* NP contribution to the Zff vertex */
    if (myEW.checkNPZff_linearized() && SM.ModelName().compare("StandardModel") != 0)
        Gamma_Z = myEW.getMyEW_NPZff().GammaZ(Gamma_Z);
      
    return Gamma_Z;
}
        
