/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Rlepton.h"

double Rlepton::computeThValue()
{
    return SM.R0_f(SM.getLeptons(SM.ELECTRON));
}


