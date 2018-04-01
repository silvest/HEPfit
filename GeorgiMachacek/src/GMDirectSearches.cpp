/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GMDirectSearches.h"
#include "StandardModel.h"



BR_H1_hh_GM::BR_H1_hh_GM(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double BR_H1_hh_GM::computeThValue()
{
    return myGM.getMyGMCache()->Br_Htohh;
}



Hobs_ggF_H1_tautau_ATLAS8::Hobs_ggF_H1_tautau_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_ggF_H1_tautau_ATLAS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_ggF_H_tautau_ATLAS8;
}

Robs_ggF_H1_tautau_ATLAS8::Robs_ggF_H1_tautau_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Robs_ggF_H1_tautau_ATLAS8::computeThValue()
{
    return myGM.getMyGMCache()->R_ggF_H_tautau_ATLAS8;
}




Hobs_pp_H1_hh_bbbb_CMS13::Hobs_pp_H1_hh_bbbb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H1_hh_bbbb_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H_hh_bbbb_CMS13;
}

Robs_pp_H1_hh_bbbb_CMS13::Robs_pp_H1_hh_bbbb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Robs_pp_H1_hh_bbbb_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->R_pp_H_hh_bbbb_CMS13;
}




log10_ggF_H1_tautau_TH8::log10_ggF_H1_tautau_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_ggF_H1_tautau_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->ggF_H_tautau_TH8);
}





log10_pp_H1_hh_bbbb_TH13::log10_pp_H1_hh_bbbb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H1_hh_bbbb_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H_hh_bbbb_TH13);
}
