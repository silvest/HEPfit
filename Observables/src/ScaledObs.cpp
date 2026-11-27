/*
 * Copyright (C) 2026 HEPfit Collaboration
 * For the licensing terms see doc/COPYING.
 */

#include "ScaledObs.h"
#include "StandardModel.h"

ScaledObs::ScaledObs(const StandardModel& SM_i, const std::string parameter, ThObservable* obs_i)
: ThObservable(SM_i), scale(SM_i.getModelParam(parameter)), obs(obs_i)
{
}

double ScaledObs::computeThValue()
{
    obs->setBinMin(getBinMin());
    obs->setBinMax(getBinMax());
    return scale * obs->computeThValue();
}
