/*
 * Copyright (C) 2019 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "C2beta.h"
#include "StandardModel.h"

double C2beta::computeThValue() 
{
    return cos(-M21_Bd(FULLNLO).arg() + 2.*SM.getPhiBd());
}
