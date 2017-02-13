/*
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ParamObs.h"
#include "StandardModel.h"

ParamObs::ParamObs(const StandardModel& SM_i, std::string name) 
: ThObservable(SM_i), param(SM_i.getModelParam(name))
{
}

double ParamObs::computeThValue()
{
    return param;
}
