/* 
 * Copyright (C) 2023 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMLowMass.h"
#include "GeneralTHDM.h"
#include "GeneralTHDMcache.h"

Hobs_pp_h_phi3phi3_mumutautau_CMS13::Hobs_pp_h_phi3phi3_mumutautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_h_phi3phi3_mumutautau_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_h_phi3phi3_mumutautau_CMS13;
}