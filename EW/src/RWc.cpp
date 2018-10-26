/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "RWc.h"
#include "StandardModel.h"

double RWc::computeThValue()
{    
    double GammWcs = SM.GammaW(SM.getQuarks(StandardModel::CHARM), SM.getQuarks(StandardModel::STRANGE));
    double Gammhad = GammWcs;
    
//  Add the ud decays into the hadronic part
    Gammhad += SM.GammaW(SM.getQuarks(StandardModel::UP), SM.getQuarks(StandardModel::DOWN));

    return GammWcs/Gammhad;
}


