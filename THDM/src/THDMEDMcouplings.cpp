/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "THDMEDMcouplings.h"
#include "THDM.h"

THDMEDMcouplings::THDMEDMcouplings(const THDM& THDM_i) 
: EffectiveEDMcouplings(THDM_i)/*, myTHDM(THDM_i)*/
{}

THDMEDMcouplings::~THDMEDMcouplings() 
{}

double THDMEDMcouplings::kappa_t() 
{
    return 0.;
}

double THDMEDMcouplings::kappa_b() 
{
    return 0.;
}

