/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Rlepton.h"


double Rlepton::computeThValue() 
{
    double R0_l = SM.R0_f(SM.getLeptons(SM.ELECTRON));

    return R0_l;
}
        

