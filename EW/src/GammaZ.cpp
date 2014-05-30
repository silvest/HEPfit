/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GammaZ.h"


double GammaZ::computeThValue() 
{ 
//    double Gamma_Z = SM.Gamma_Z();
//
//    /* NP contribution to the Zff vertex */
//    if (SM.checkNPZff_linearized() && SM.ModelName().compare("StandardModel") != 0)
//        Gamma_Z = SM.getMyEW_NPZff().GammaZ(Gamma_Z);
      
    return SM.Gamma_Z();
}
        
