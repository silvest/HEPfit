/*
 * Copyright (C) 2014 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ParamObs.h"

ParamObs::ParamObs(const StandardModel& SM_i, std::string name) 
: ThObservable(SM_i), param(SM_i.getModelParam(name))
{
}

double ParamObs::computeThValue()
{
    return param;
}
