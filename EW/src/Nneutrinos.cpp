/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Nneutrinos.h"
#include "StandardModel.h"

double Nneutrinos::computeThValue()
{
    return SM.N_nu();
}

