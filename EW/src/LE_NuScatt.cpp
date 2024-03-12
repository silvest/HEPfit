/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "LE_NuScatt.h"
#include "StandardModel.h"

double gLnuN2::computeThValue()
{    
    return SM.gLnuN2();
}


double gRnuN2::computeThValue()
{    
    return SM.gRnuN2();
}


double ThetaLnuN::computeThValue()
{    
    return SM.ThetaLnuN();
}


double ThetaRnuN::computeThValue()
{    
    return SM.ThetaRnuN();
}


double gVnue::computeThValue()
{    
    return SM.gVnue();
}


double gAnue::computeThValue()
{    
    return SM.gAnue();
}