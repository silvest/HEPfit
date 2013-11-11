/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GammaW.h"


double GammaW::computeThValue() 
{  
    double Gamma_W;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        Gamma_W = myEW.getMyEW_CHMN().GammaW();
    else if (myEWTYPE==EW::EWABC || myEWTYPE==EW::EWABC2) 
        throw std::runtime_error("GammaW::computeThValue() is not implemented for EW::EWABC");  
    else
        Gamma_W = SM.GammaW();
 
    return Gamma_W;
}
