/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ModelParaVsObs.h"

ModelParaVsObs::ModelParaVsObs(const std::string name_i,
        const std::string ParaName_i, const std::string ParaLabel_i,
        const double ParaMin_i, const double ParaMax_i,
        const std::string ObsName_i, const std::string ObsLabel_i,
        const double ObsMin_i, const double ObsMax_i,
        ThObservable * Obs_i)
: Observable(name_i, ObsName_i, ObsLabel_i, false, ObsMin_i, ObsMax_i, Obs_i)
{
    ParaName = ParaName_i;
    ParaLabel = ParaLabel_i;
    ParaMin = ParaMin_i;
    ParaMax = ParaMax_i;
}

