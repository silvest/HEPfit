/*
 * Copyright (C) 2015 HEPfit Collaboration
 *
 * For the licensing terms see doc/COPYING.
 */

#include "MtMSbar.h"
#include "StandardModel.h"

double MtMSbar::computeThValue() 
{
    // three-loop pole <-> MSbar relation (see QCD::Mbar2Mp)
    return SM.Mp2Mbar(SM.getMtpole(), QCD::TOP, FULLNNNLO);
}
