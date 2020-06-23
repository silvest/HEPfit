/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "RVllp.h"
#include "StandardModel.h"

double RWmue::computeThValue()
{    
    return SM.RWlilj(SM.getLeptons(StandardModel::MU),SM.getLeptons(StandardModel::ELECTRON));
}

double RWtaue::computeThValue()
{    
    return SM.RWlilj(SM.getLeptons(StandardModel::TAU),SM.getLeptons(StandardModel::ELECTRON));
}

double RWtaumu::computeThValue()
{    
    return SM.RWlilj(SM.getLeptons(StandardModel::TAU),SM.getLeptons(StandardModel::MU));
}


double RZmue::computeThValue()
{    
    return SM.RZlilj(SM.getLeptons(StandardModel::MU),SM.getLeptons(StandardModel::ELECTRON));
}


double RZtaue::computeThValue()
{    
    return SM.RZlilj(SM.getLeptons(StandardModel::TAU),SM.getLeptons(StandardModel::ELECTRON));
}

double RZtaumu::computeThValue()
{    
    return SM.RZlilj(SM.getLeptons(StandardModel::TAU),SM.getLeptons(StandardModel::MU));
}


