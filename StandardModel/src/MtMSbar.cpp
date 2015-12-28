/*
 * Copyright (C) 2015 HEPfit Collaboration
 *
 * For the licensing terms see doc/COPYING.
 */

#include "MtMSbar.h"

double MtMSbar::computeThValue() 
{
    return SM.Mp2Mbar(SM.getMtpole());
}
