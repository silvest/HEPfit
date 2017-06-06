/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMEDMcouplings.h"
#include "GeneralTHDM.h"

GeneralTHDMEDMcouplings::GeneralTHDMEDMcouplings(const GeneralTHDM& GeneralTHDM_i) 
: EffectiveEDMcouplings(GeneralTHDM_i)/*, myGeneralTHDM(GeneralTHDM_i)*/
{}

GeneralTHDMEDMcouplings::~GeneralTHDMEDMcouplings() 
{}

double GeneralTHDMEDMcouplings::kappa_t() 
{
    return 0.;
}

double GeneralTHDMEDMcouplings::kappa_b() 
{
    return 0.;
}
