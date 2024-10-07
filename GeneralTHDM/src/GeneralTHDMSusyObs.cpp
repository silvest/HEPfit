/* 
 * Copyright (C) 2024 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMSusyObs.h"
#include "GeneralTHDM.h"
#include "GeneralTHDMcache.h"

Hobs_pp_HpHm_taunutaunu_ATLAS13::Hobs_pp_HpHm_taunutaunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_HpHm_taunutaunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_HpHm_taunutaunu_ATLAS13;
}

Hobs_pp_HpHm_taunutaunu_CMS13::Hobs_pp_HpHm_taunutaunu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_HpHm_taunutaunu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_HpHm_taunutaunu_CMS13;
}

Hobs_pp_HpHm_munumunu_ATLAS13::Hobs_pp_HpHm_munumunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_HpHm_munumunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_HpHm_munumunu_ATLAS13;
}

Hobs_pp_HpHm_munumunu_CMS13::Hobs_pp_HpHm_munumunu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_HpHm_munumunu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_HpHm_munumunu_CMS13;
}

Hobs_HpHm_munumunu_LEP208::Hobs_HpHm_munumunu_LEP208(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_HpHm_munumunu_LEP208::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_HpHm_munumunu_LEP208;
}
