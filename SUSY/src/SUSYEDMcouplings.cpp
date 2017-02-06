/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "SUSYEDMcouplings.h"
#include "SUSY.h"

SUSYEDMcouplings::SUSYEDMcouplings(const SUSY& SUSY_i) 
: EffectiveEDMcouplings(SUSY_i), mySUSY(SUSY_i)
{}

SUSYEDMcouplings::~SUSYEDMcouplings() 
{}

double SUSYEDMcouplings::kappa_t() 
{
    return 0.;
}

double SUSYEDMcouplings::kappa_b() 
{
    return 0.;
}

