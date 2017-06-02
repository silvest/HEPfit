/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EffectiveEDMcouplings.h"

EffectiveEDMcouplings::EffectiveEDMcouplings(const StandardModel& SM_i) 
/*: SM(SM_i)*/
{}

EffectiveEDMcouplings::~EffectiveEDMcouplings() 
{}

double EffectiveEDMcouplings::kappa_t() 
{
    return 0.;
}

double EffectiveEDMcouplings::kappa_b() 
{
    return 0.;
}