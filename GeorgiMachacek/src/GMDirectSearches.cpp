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



Hobs_tt_H1_tt_ATLAS13::Hobs_tt_H1_tt_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_tt_H1_tt_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_tt_H_tt_ATLAS13;
}



Hobs_bb_H1_tt_ATLAS13::Hobs_bb_H1_tt_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_bb_H1_tt_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_bb_H_tt_ATLAS13;
}



Hobs_tt_H3_tt_ATLAS13::Hobs_tt_H3_tt_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_tt_H3_tt_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_tt_A_tt_ATLAS13;
}



Hobs_bb_H3_tt_ATLAS13::Hobs_bb_H3_tt_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_bb_H3_tt_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_bb_A_tt_ATLAS13;
}



Hobs_bb_H1_bb_CMS8::Hobs_bb_H1_bb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_bb_H1_bb_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_bb_H_bb_CMS8;
}



Hobs_gg_H1_bb_CMS8::Hobs_gg_H1_bb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H1_bb_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_H_bb_CMS8;
}



Hobs_pp_H1_bb_CMS13::Hobs_pp_H1_bb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H1_bb_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H_bb_CMS13;
}



Hobs_bb_H1_bb_CMS13::Hobs_bb_H1_bb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_bb_H1_bb_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_bb_H_bb_CMS13;
}



Hobs_bb_H3_bb_CMS8::Hobs_bb_H3_bb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_bb_H3_bb_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_bb_A_bb_CMS8;
}



Hobs_gg_H3_bb_CMS8::Hobs_gg_H3_bb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H3_bb_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_A_bb_CMS8;
}



Hobs_pp_H3_bb_CMS13::Hobs_pp_H3_bb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H3_bb_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_A_bb_CMS13;
}



Hobs_bb_H3_bb_CMS13::Hobs_bb_H3_bb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_bb_H3_bb_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_bb_A_bb_CMS13;
}



Hobs_gg_H1_tautau_CMS8::Hobs_gg_H1_tautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H1_tautau_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_H_tautau_CMS8;
}



Hobs_bb_H1_tautau_CMS8::Hobs_bb_H1_tautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_bb_H1_tautau_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_bb_H_tautau_CMS8;
}



Hobs_gg_H1_tautau_ATLAS13::Hobs_gg_H1_tautau_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H1_tautau_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_H_tautau_ATLAS13;
}



Hobs_gg_H1_tautau_CMS13::Hobs_gg_H1_tautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H1_tautau_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_H_tautau_CMS13;
}



Hobs_bb_H1_tautau_ATLAS13::Hobs_bb_H1_tautau_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_bb_H1_tautau_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_bb_H_tautau_ATLAS13;
}



Hobs_bb_H1_tautau_CMS13::Hobs_bb_H1_tautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_bb_H1_tautau_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_bb_H_tautau_CMS13;
}



Hobs_gg_H1_tautau_ATLAS8::Hobs_gg_H1_tautau_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H1_tautau_ATLAS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_H_tautau_ATLAS8;
}



Hobs_bb_H1_tautau_ATLAS8::Hobs_bb_H1_tautau_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_bb_H1_tautau_ATLAS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_bb_H_tautau_ATLAS8;
}



Hobs_gg_H3_tautau_ATLAS8::Hobs_gg_H3_tautau_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H3_tautau_ATLAS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_A_tautau_ATLAS8;
}



Hobs_gg_H3_tautau_CMS8::Hobs_gg_H3_tautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H3_tautau_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_A_tautau_CMS8;
}



Hobs_bb_H3_tautau_ATLAS8::Hobs_bb_H3_tautau_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_bb_H3_tautau_ATLAS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_bb_A_tautau_ATLAS8;
}



Hobs_bb_H3_tautau_CMS8::Hobs_bb_H3_tautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_bb_H3_tautau_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_bb_A_tautau_CMS8;
}



Hobs_gg_H3_tautau_ATLAS13::Hobs_gg_H3_tautau_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H3_tautau_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_A_tautau_ATLAS13;
}



Hobs_gg_H3_tautau_CMS13::Hobs_gg_H3_tautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H3_tautau_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_A_tautau_CMS13;
}



Hobs_bb_H3_tautau_ATLAS13::Hobs_bb_H3_tautau_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_bb_H3_tautau_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_bb_A_tautau_ATLAS13;
}



Hobs_bb_H3_tautau_CMS13::Hobs_bb_H3_tautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_bb_H3_tautau_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_bb_A_tautau_CMS13;
}



Hobs_gg_H1_gaga_ATLAS8::Hobs_gg_H1_gaga_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H1_gaga_ATLAS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_H_gaga_ATLAS8;
}



Hobs_pp_H1_gaga_ATLAS13::Hobs_pp_H1_gaga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H1_gaga_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H_gaga_ATLAS13;
}



Hobs_gg_H1_gaga_CMS13::Hobs_gg_H1_gaga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H1_gaga_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_H_gaga_CMS13;
}



Hobs_gg_H3_gaga_ATLAS8::Hobs_gg_H3_gaga_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H3_gaga_ATLAS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_A_gaga_ATLAS8;
}



Hobs_pp_H3_gaga_ATLAS13::Hobs_pp_H3_gaga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H3_gaga_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_A_gaga_ATLAS13;
}



Hobs_gg_H3_gaga_CMS13::Hobs_gg_H3_gaga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H3_gaga_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_A_gaga_CMS13;
}



Hobs_pp_H5_gaga_ATLAS13::Hobs_pp_H5_gaga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H5_gaga_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H5_gaga_ATLAS13;
}



Hobs_pp_H1_Zga_llga_ATLAS8::Hobs_pp_H1_Zga_llga_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H1_Zga_llga_ATLAS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H_Zga_llga_ATLAS8;
}



Hobs_pp_H1_Zga_llga_CMS8::Hobs_pp_H1_Zga_llga_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H1_Zga_llga_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H_Zga_llga_CMS8;
}



Hobs_gg_H1_Zga_llga_ATLAS13::Hobs_gg_H1_Zga_llga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H1_Zga_llga_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_H_Zga_llga_ATLAS13;
}



Hobs_gg_H1_Zga_CMS13::Hobs_gg_H1_Zga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H1_Zga_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_H_Zga_CMS13;
}



Hobs_pp_H3_Zga_llga_ATLAS8::Hobs_pp_H3_Zga_llga_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H3_Zga_llga_ATLAS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_A_Zga_llga_ATLAS8;
}



Hobs_pp_H3_Zga_llga_CMS8::Hobs_pp_H3_Zga_llga_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H3_Zga_llga_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_A_Zga_llga_CMS8;
}



Hobs_gg_H3_Zga_llga_ATLAS13::Hobs_gg_H3_Zga_llga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H3_Zga_llga_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_A_Zga_llga_ATLAS13;
}



Hobs_gg_H3_Zga_CMS13::Hobs_gg_H3_Zga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H3_Zga_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_A_Zga_CMS13;
}



Hobs_pp_H5_Zga_llga_ATLAS8::Hobs_pp_H5_Zga_llga_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H5_Zga_llga_ATLAS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H5_Zga_llga_ATLAS8;
}



Hobs_pp_H5_Zga_llga_CMS8::Hobs_pp_H5_Zga_llga_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H5_Zga_llga_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H5_Zga_llga_CMS8;
}



Hobs_gg_H1_ZZ_ATLAS8::Hobs_gg_H1_ZZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H1_ZZ_ATLAS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_H_ZZ_ATLAS8;
}



Hobs_VV_H1_ZZ_ATLAS8::Hobs_VV_H1_ZZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_VV_H1_ZZ_ATLAS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_VV_H_ZZ_ATLAS8;
}



Hobs_gg_H1_ZZ_llllnunu_ATLAS13::Hobs_gg_H1_ZZ_llllnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H1_ZZ_llllnunu_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_H_ZZ_llllnunu_ATLAS13;
}



Hobs_VV_H1_ZZ_llllnunu_ATLAS13::Hobs_VV_H1_ZZ_llllnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_VV_H1_ZZ_llllnunu_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_VV_H_ZZ_llllnunu_ATLAS13;
}



Hobs_gg_H1_ZZ_qqllnunu_ATLAS13::Hobs_gg_H1_ZZ_qqllnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H1_ZZ_qqllnunu_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_H_ZZ_qqllnunu_ATLAS13;
}



Hobs_VV_H1_ZZ_qqllnunu_ATLAS13::Hobs_VV_H1_ZZ_qqllnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_VV_H1_ZZ_qqllnunu_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_VV_H_ZZ_qqllnunu_ATLAS13;
}



Hobs_pp_H1_ZZ_llqqnunull_CMS13::Hobs_pp_H1_ZZ_llqqnunull_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H1_ZZ_llqqnunull_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H_ZZ_llqqnunull_CMS13;
}



Hobs_VV_H1_ZZ_llqqnunull_CMS13::Hobs_VV_H1_ZZ_llqqnunull_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_VV_H1_ZZ_llqqnunull_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_VV_H_ZZ_llqqnunull_CMS13;
}



Hobs_pp_H1_ZZ_qqnunu_CMS13::Hobs_pp_H1_ZZ_qqnunu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H1_ZZ_qqnunu_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H_ZZ_qqnunu_CMS13;
}



Hobs_VV_H5_ZZ_ATLAS8::Hobs_VV_H5_ZZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_VV_H5_ZZ_ATLAS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_VV_H5_ZZ_ATLAS8;
}



Hobs_VV_H5_ZZ_llllnunu_ATLAS13::Hobs_VV_H5_ZZ_llllnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_VV_H5_ZZ_llllnunu_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_VV_H5_ZZ_llllnunu_ATLAS13;
}



Hobs_VV_H5_ZZ_qqllnunu_ATLAS13::Hobs_VV_H5_ZZ_qqllnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_VV_H5_ZZ_qqllnunu_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_VV_H5_ZZ_qqllnunu_ATLAS13;
}



Hobs_pp_H5_ZZ_llqqnunull_CMS13::Hobs_pp_H5_ZZ_llqqnunull_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H5_ZZ_llqqnunull_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H5_ZZ_llqqnunull_CMS13;
}



Hobs_VV_H5_ZZ_llqqnunull_CMS13::Hobs_VV_H5_ZZ_llqqnunull_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_VV_H5_ZZ_llqqnunull_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_VV_H5_ZZ_llqqnunull_CMS13;
}



Hobs_pp_H5_ZZ_qqnunu_CMS13::Hobs_pp_H5_ZZ_qqnunu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H5_ZZ_qqnunu_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H5_ZZ_qqnunu_CMS13;
}



Hobs_gg_H1_WW_ATLAS8::Hobs_gg_H1_WW_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H1_WW_ATLAS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_H_WW_ATLAS8;
}



Hobs_VV_H1_WW_ATLAS8::Hobs_VV_H1_WW_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_VV_H1_WW_ATLAS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_VV_H_WW_ATLAS8;
}



Hobs_gg_H1_WW_enumunu_ATLAS13::Hobs_gg_H1_WW_enumunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H1_WW_enumunu_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_H_WW_enumunu_ATLAS13;
}



Hobs_VV_H1_WW_enumunu_ATLAS13::Hobs_VV_H1_WW_enumunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_VV_H1_WW_enumunu_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_VV_H_WW_enumunu_ATLAS13;
}



Hobs_gg_H1_WW_lnuqq_ATLAS13::Hobs_gg_H1_WW_lnuqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H1_WW_lnuqq_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_H_WW_lnuqq_ATLAS13;
}



Hobs_VV_H1_WW_lnuqq_ATLAS13::Hobs_VV_H1_WW_lnuqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_VV_H1_WW_lnuqq_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_VV_H_WW_lnuqq_ATLAS13;
}



Hobs_ggVV_H1_WW_lnulnu_CMS13::Hobs_ggVV_H1_WW_lnulnu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_ggVV_H1_WW_lnulnu_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_ggVV_H_WW_lnulnu_CMS13;
}



Hobs_pp_H1_WW_lnuqq_CMS13::Hobs_pp_H1_WW_lnuqq_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H1_WW_lnuqq_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H_WW_lnuqq_CMS13;
}



Hobs_VV_H5_WW_ATLAS8::Hobs_VV_H5_WW_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_VV_H5_WW_ATLAS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_VV_H5_WW_ATLAS8;
}



Hobs_VV_H5_WW_enumunu_ATLAS13::Hobs_VV_H5_WW_enumunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_VV_H5_WW_enumunu_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_VV_H5_WW_enumunu_ATLAS13;
}



Hobs_VV_H5_WW_lnuqq_ATLAS13::Hobs_VV_H5_WW_lnuqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_VV_H5_WW_lnuqq_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_VV_H5_WW_lnuqq_ATLAS13;
}



Hobs_ggVV_H5_WW_lnulnu_CMS13::Hobs_ggVV_H5_WW_lnulnu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_ggVV_H5_WW_lnulnu_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_ggVV_H5_WW_lnulnu_CMS13;
}



Hobs_pp_H5_WW_lnuqq_CMS13::Hobs_pp_H5_WW_lnuqq_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H5_WW_lnuqq_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H5_WW_lnuqq_CMS13;
}



Hobs_mu_pp_H1_VV_CMS8::Hobs_mu_pp_H1_VV_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_mu_pp_H1_VV_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_mu_pp_H_VV_CMS8;
}



Hobs_pp_H1_VV_qqqq_ATLAS13::Hobs_pp_H1_VV_qqqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H1_VV_qqqq_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H_VV_qqqq_ATLAS13;
}



Hobs_mu_pp_H5_VV_CMS8::Hobs_mu_pp_H5_VV_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_mu_pp_H5_VV_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_mu_pp_H5_VV_CMS8;
}



Hobs_pp_H5_VV_qqqq_ATLAS13::Hobs_pp_H5_VV_qqqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H5_VV_qqqq_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H5_VV_qqqq_ATLAS13;
}



Hobs_gg_H1_hh_ATLAS8::Hobs_gg_H1_hh_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H1_hh_ATLAS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_H_hh_ATLAS8;
}



Hobs_pp_H1_hh_bbbb_CMS8::Hobs_pp_H1_hh_bbbb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H1_hh_bbbb_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H_hh_bbbb_CMS8;
}



Hobs_pp_H1_hh_gagabb_CMS8::Hobs_pp_H1_hh_gagabb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H1_hh_gagabb_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H_hh_gagabb_CMS8;
}



Hobs_gg_H1_hh_bbtautau_CMS8::Hobs_gg_H1_hh_bbtautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H1_hh_bbtautau_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_H_hh_bbtautau_CMS8;
}



Hobs_pp_H1_hh_bbtautau_CMS8::Hobs_pp_H1_hh_bbtautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H1_hh_bbtautau_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H_hh_bbtautau_CMS8;
}



Hobs_pp_H1_hh_bbbb_ATLAS13::Hobs_pp_H1_hh_bbbb_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H1_hh_bbbb_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H_hh_bbbb_ATLAS13;
}



Hobs_pp_H1_hh_bbbb_CMS13::Hobs_pp_H1_hh_bbbb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H1_hh_bbbb_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H_hh_bbbb_CMS13;
}



Hobs_gg_H1_hh_bbbb_CMS13::Hobs_gg_H1_hh_bbbb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H1_hh_bbbb_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_H_hh_bbbb_CMS13;
}



Hobs_pp_H1_hh_gagabb_ATLAS13::Hobs_pp_H1_hh_gagabb_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H1_hh_gagabb_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H_hh_gagabb_ATLAS13;
}



Hobs_pp_H1_hh_gagabb_CMS13::Hobs_pp_H1_hh_gagabb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H1_hh_gagabb_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H_hh_gagabb_CMS13;
}



Hobs_pp_H1_hh_bbtautau_CMS13::Hobs_pp_H1_hh_bbtautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H1_hh_bbtautau_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H_hh_bbtautau_CMS13;
}



Hobs_pp_H1_hh_bblnulnu_CMS13::Hobs_pp_H1_hh_bblnulnu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H1_hh_bblnulnu_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H_hh_bblnulnu_CMS13;
}



Hobs_gg_H1_hh_gagaWW_ATLAS13::Hobs_gg_H1_hh_gagaWW_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H1_hh_gagaWW_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_H_hh_gagaWW_ATLAS13;
}



Hobs_pp_H5_hh_bbbb_CMS8::Hobs_pp_H5_hh_bbbb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H5_hh_bbbb_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H5_hh_bbbb_CMS8;
}



Hobs_pp_H5_hh_gagabb_CMS8::Hobs_pp_H5_hh_gagabb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H5_hh_gagabb_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H5_hh_gagabb_CMS8;
}



Hobs_pp_H5_hh_bbtautau_CMS8::Hobs_pp_H5_hh_bbtautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H5_hh_bbtautau_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H5_hh_bbtautau_CMS8;
}



Hobs_pp_H5_hh_bbbb_ATLAS13::Hobs_pp_H5_hh_bbbb_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H5_hh_bbbb_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H5_hh_bbbb_ATLAS13;
}



Hobs_pp_H5_hh_bbbb_CMS13::Hobs_pp_H5_hh_bbbb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H5_hh_bbbb_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H5_hh_bbbb_CMS13;
}



Hobs_pp_H5_hh_gagabb_ATLAS13::Hobs_pp_H5_hh_gagabb_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H5_hh_gagabb_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H5_hh_gagabb_ATLAS13;
}



Hobs_pp_H5_hh_gagabb_CMS13::Hobs_pp_H5_hh_gagabb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H5_hh_gagabb_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H5_hh_gagabb_CMS13;
}



Hobs_pp_H5_hh_bbtautau_CMS13::Hobs_pp_H5_hh_bbtautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H5_hh_bbtautau_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H5_hh_bbtautau_CMS13;
}



Hobs_pp_H5_hh_bblnulnu_CMS13::Hobs_pp_H5_hh_bblnulnu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H5_hh_bblnulnu_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H5_hh_bblnulnu_CMS13;
}



Hobs_gg_H3_hZ_bbZ_ATLAS8::Hobs_gg_H3_hZ_bbZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H3_hZ_bbZ_ATLAS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_A_hZ_bbZ_ATLAS8;
}



Hobs_gg_H3_hZ_bbll_CMS8::Hobs_gg_H3_hZ_bbll_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H3_hZ_bbll_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_A_hZ_bbll_CMS8;
}



Hobs_gg_H3_hZ_tautauZ_ATLAS8::Hobs_gg_H3_hZ_tautauZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H3_hZ_tautauZ_ATLAS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_A_hZ_tautauZ_ATLAS8;
}



Hobs_gg_H3_hZ_tautaull_CMS8::Hobs_gg_H3_hZ_tautaull_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H3_hZ_tautaull_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_A_hZ_tautaull_CMS8;
}



Hobs_gg_H3_hZ_bbZ_ATLAS13::Hobs_gg_H3_hZ_bbZ_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_gg_H3_hZ_bbZ_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_gg_A_hZ_bbZ_ATLAS13;
}



Hobs_bb_H3_hZ_bbZ_ATLAS13::Hobs_bb_H3_hZ_bbZ_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_bb_H3_hZ_bbZ_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_bb_A_hZ_bbZ_ATLAS13;
}



Hobs_pp_H3_H1Z_bbll_CMS8::Hobs_pp_H3_H1Z_bbll_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H3_H1Z_bbll_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_A_HZ_bbll_CMS8;
}



Hobs_pp_H3_H5Z_bbll_CMS8::Hobs_pp_H3_H5Z_bbll_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H3_H5Z_bbll_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_A_H5Z_bbll_CMS8;
}



Hobs_pp_H1_H3Z_bbll_CMS8::Hobs_pp_H1_H3Z_bbll_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H1_H3Z_bbll_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H_AZ_bbll_CMS8;
}



Hobs_pp_H5_H3Z_bbll_CMS8::Hobs_pp_H5_H3Z_bbll_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H5_H3Z_bbll_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H5_AZ_bbll_CMS8;
}



Hobs_pp_H3pm_taunu_ATLAS8::Hobs_pp_H3pm_taunu_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H3pm_taunu_ATLAS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_Hpm_taunu_ATLAS8;
}



Hobs_pp_H3p_taunu_CMS8::Hobs_pp_H3p_taunu_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H3p_taunu_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_Hp_taunu_CMS8;
}



Hobs_pp_H3pm_taunu_ATLAS13::Hobs_pp_H3pm_taunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H3pm_taunu_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_Hpm_taunu_ATLAS13;
}



Hobs_pp_H3pm_taunu_CMS13::Hobs_pp_H3pm_taunu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H3pm_taunu_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_Hpm_taunu_CMS13;
}



Hobs_pp_H3pm_tb_ATLAS8::Hobs_pp_H3pm_tb_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H3pm_tb_ATLAS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_Hpm_tb_ATLAS8;
}



Hobs_pp_H3p_tb_CMS8::Hobs_pp_H3p_tb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H3p_tb_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_Hp_tb_CMS8;
}



Hobs_pp_H3p_tb1_ATLAS13::Hobs_pp_H3p_tb1_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H3p_tb1_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_Hp_tb1_ATLAS13;
}



Hobs_pp_H3p_tb2_ATLAS13::Hobs_pp_H3p_tb2_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H3p_tb2_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_Hp_tb2_ATLAS13;
}



Hobs_WZ_H5pm_WZ_qqll_ATLAS8::Hobs_WZ_H5pm_WZ_qqll_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_WZ_H5pm_WZ_qqll_ATLAS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_WZ_H5pm_WZ_qqll_ATLAS8;
}



Hobs_WZ_H5pm_WZ_lnull_CMS13::Hobs_WZ_H5pm_WZ_lnull_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_WZ_H5pm_WZ_lnull_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_WZ_H5pm_WZ_lnull_CMS13;
}



Hobs_pp_H5ppmmH5mmpp_eeee_ATLAS8::Hobs_pp_H5ppmmH5mmpp_eeee_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H5ppmmH5mmpp_eeee_ATLAS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H5ppmmH5mmpp_eeee_ATLAS8;
}



Hobs_pp_H5ppmmH5mmpp_emuemu_ATLAS8::Hobs_pp_H5ppmmH5mmpp_emuemu_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H5ppmmH5mmpp_emuemu_ATLAS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H5ppmmH5mmpp_emuemu_ATLAS8;
}



Hobs_pp_H5ppmmH5mmpp_mumumumu_ATLAS8::Hobs_pp_H5ppmmH5mmpp_mumumumu_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H5ppmmH5mmpp_mumumumu_ATLAS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H5ppmmH5mmpp_mumumumu_ATLAS8;
}



Hobs_pp_H5ppmmH5mmpp_llll_ATLAS13::Hobs_pp_H5ppmmH5mmpp_llll_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H5ppmmH5mmpp_llll_ATLAS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H5ppmmH5mmpp_llll_ATLAS13;
}



Hobs_pp_H5ppmm_WW_jjll_CMS8::Hobs_pp_H5ppmm_WW_jjll_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H5ppmm_WW_jjll_CMS8::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H5ppmm_WW_jjll_CMS8;
}



Hobs_pp_H5ppmm_WW_jjll_CMS13::Hobs_pp_H5ppmm_WW_jjll_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double Hobs_pp_H5ppmm_WW_jjll_CMS13::computeThValue()
{
    return myGM.getMyGMCache()->THoEX_pp_H5ppmm_WW_jjll_CMS13;
}



log10_tt_H1_tt_TH13::log10_tt_H1_tt_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_tt_H1_tt_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->tt_H_tt_TH13);
}



log10_bb_H1_tt_TH13::log10_bb_H1_tt_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_bb_H1_tt_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->bb_H_tt_TH13);
}



log10_tt_H3_tt_TH13::log10_tt_H3_tt_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_tt_H3_tt_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->tt_A_tt_TH13);
}



log10_bb_H3_tt_TH13::log10_bb_H3_tt_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_bb_H3_tt_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->bb_A_tt_TH13);
}



log10_bb_H1_bb_TH8::log10_bb_H1_bb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_bb_H1_bb_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->bb_H_bb_TH8);
}



log10_gg_H1_bb_TH8::log10_gg_H1_bb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H1_bb_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->gg_H_bb_TH8);
}



log10_pp_H1_bb_TH13::log10_pp_H1_bb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H1_bb_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H_bb_TH13);
}



log10_bb_H1_bb_TH13::log10_bb_H1_bb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_bb_H1_bb_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->bb_H_bb_TH13);
}



log10_bb_H3_bb_TH8::log10_bb_H3_bb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_bb_H3_bb_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->bb_A_bb_TH8);
}



log10_gg_H3_bb_TH8::log10_gg_H3_bb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H3_bb_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->gg_A_bb_TH8);
}



log10_pp_H3_bb_TH13::log10_pp_H3_bb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H3_bb_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_A_bb_TH13);
}



log10_bb_H3_bb_TH13::log10_bb_H3_bb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_bb_H3_bb_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->bb_A_bb_TH13);
}



log10_gg_H1_tautau_TH8::log10_gg_H1_tautau_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H1_tautau_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->gg_H_tautau_TH8);
}



log10_bb_H1_tautau_TH8::log10_bb_H1_tautau_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_bb_H1_tautau_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->bb_H_tautau_TH8);
}



log10_gg_H1_tautau_TH13::log10_gg_H1_tautau_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H1_tautau_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->gg_H_tautau_TH13);
}



log10_bb_H1_tautau_TH13::log10_bb_H1_tautau_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_bb_H1_tautau_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->bb_H_tautau_TH13);
}



log10_gg_H3_tautau_TH8::log10_gg_H3_tautau_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H3_tautau_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->gg_A_tautau_TH8);
}



log10_bb_H3_tautau_TH8::log10_bb_H3_tautau_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_bb_H3_tautau_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->bb_A_tautau_TH8);
}



log10_gg_H3_tautau_TH13::log10_gg_H3_tautau_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H3_tautau_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->gg_A_tautau_TH13);
}



log10_bb_H3_tautau_TH13::log10_bb_H3_tautau_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_bb_H3_tautau_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->bb_A_tautau_TH13);
}



log10_gg_H1_gaga_TH8::log10_gg_H1_gaga_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H1_gaga_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->gg_H_gaga_TH8);
}



log10_pp_H1_gaga_TH13::log10_pp_H1_gaga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H1_gaga_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H_gaga_TH13);
}



log10_gg_H1_gaga_TH13::log10_gg_H1_gaga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H1_gaga_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->gg_H_gaga_TH13);
}



log10_gg_H3_gaga_TH8::log10_gg_H3_gaga_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H3_gaga_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->gg_A_gaga_TH8);
}



log10_pp_H3_gaga_TH13::log10_pp_H3_gaga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H3_gaga_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_A_gaga_TH13);
}



log10_gg_H3_gaga_TH13::log10_gg_H3_gaga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H3_gaga_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->gg_A_gaga_TH13);
}



log10_pp_H5_gaga_TH13::log10_pp_H5_gaga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H5_gaga_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H5_gaga_TH13);
}



log10_pp_H1_Zga_llga_TH8::log10_pp_H1_Zga_llga_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H1_Zga_llga_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H_Zga_llga_TH8);
}



log10_gg_H1_Zga_llga_TH13::log10_gg_H1_Zga_llga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H1_Zga_llga_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->gg_H_Zga_llga_TH13);
}



log10_gg_H1_Zga_TH13::log10_gg_H1_Zga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H1_Zga_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->gg_H_Zga_TH13);
}



log10_pp_H3_Zga_llga_TH8::log10_pp_H3_Zga_llga_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H3_Zga_llga_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_A_Zga_llga_TH8);
}



log10_gg_H3_Zga_llga_TH13::log10_gg_H3_Zga_llga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H3_Zga_llga_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->gg_A_Zga_llga_TH13);
}



log10_gg_H3_Zga_TH13::log10_gg_H3_Zga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H3_Zga_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->gg_A_Zga_TH13);
}



log10_pp_H5_Zga_llga_TH8::log10_pp_H5_Zga_llga_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H5_Zga_llga_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H5_Zga_llga_TH8);
}



log10_gg_H1_ZZ_TH8::log10_gg_H1_ZZ_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H1_ZZ_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->gg_H_ZZ_TH8);
}



log10_VV_H1_ZZ_TH8::log10_VV_H1_ZZ_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_VV_H1_ZZ_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->VV_H_ZZ_TH8);
}



log10_gg_H1_ZZ_llllnunu_TH13::log10_gg_H1_ZZ_llllnunu_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H1_ZZ_llllnunu_TH13::computeThValue()
{
    return log10(1.);
}



log10_VV_H1_ZZ_llllnunu_TH13::log10_VV_H1_ZZ_llllnunu_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_VV_H1_ZZ_llllnunu_TH13::computeThValue()
{
    return log10(1.);
}



log10_gg_H1_ZZ_qqllnunu_TH13::log10_gg_H1_ZZ_qqllnunu_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H1_ZZ_qqllnunu_TH13::computeThValue()
{
    return log10(1.);
}



log10_VV_H1_ZZ_qqllnunu_TH13::log10_VV_H1_ZZ_qqllnunu_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_VV_H1_ZZ_qqllnunu_TH13::computeThValue()
{
    return log10(1.);
}



log10_pp_H1_ZZ_llqqnunull_TH13::log10_pp_H1_ZZ_llqqnunull_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H1_ZZ_llqqnunull_TH13::computeThValue()
{
    return log10(1.);
}



log10_VV_H1_ZZ_llqqnunull_TH13::log10_VV_H1_ZZ_llqqnunull_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_VV_H1_ZZ_llqqnunull_TH13::computeThValue()
{
    return log10(1.);
}



log10_pp_H1_ZZ_qqnunu_TH13::log10_pp_H1_ZZ_qqnunu_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H1_ZZ_qqnunu_TH13::computeThValue()
{
    return log10(1.);
}



log10_VV_H5_ZZ_TH8::log10_VV_H5_ZZ_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_VV_H5_ZZ_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->VV_H5_ZZ_TH8);
}



log10_VV_H5_ZZ_llllnunu_TH13::log10_VV_H5_ZZ_llllnunu_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_VV_H5_ZZ_llllnunu_TH13::computeThValue()
{
    return log10(1.);
}



log10_VV_H5_ZZ_qqllnunu_TH13::log10_VV_H5_ZZ_qqllnunu_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_VV_H5_ZZ_qqllnunu_TH13::computeThValue()
{
    return log10(1.);
}



log10_pp_H5_ZZ_llqqnunull_TH13::log10_pp_H5_ZZ_llqqnunull_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H5_ZZ_llqqnunull_TH13::computeThValue()
{
    return log10(1.);
}



log10_VV_H5_ZZ_llqqnunull_TH13::log10_VV_H5_ZZ_llqqnunull_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_VV_H5_ZZ_llqqnunull_TH13::computeThValue()
{
    return log10(1.);
}



log10_pp_H5_ZZ_qqnunu_TH13::log10_pp_H5_ZZ_qqnunu_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H5_ZZ_qqnunu_TH13::computeThValue()
{
    return log10(1.);
}



log10_gg_H1_WW_TH8::log10_gg_H1_WW_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H1_WW_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->gg_H_WW_TH8);
}



log10_VV_H1_WW_TH8::log10_VV_H1_WW_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_VV_H1_WW_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->VV_H_WW_TH8);
}



log10_gg_H1_WW_enumunu_TH13::log10_gg_H1_WW_enumunu_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H1_WW_enumunu_TH13::computeThValue()
{
    return log10(1.);
}



log10_VV_H1_WW_enumunu_TH13::log10_VV_H1_WW_enumunu_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_VV_H1_WW_enumunu_TH13::computeThValue()
{
    return log10(1.);
}



log10_gg_H1_WW_lnuqq_TH13::log10_gg_H1_WW_lnuqq_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H1_WW_lnuqq_TH13::computeThValue()
{
    return log10(1.);
}



log10_VV_H1_WW_lnuqq_TH13::log10_VV_H1_WW_lnuqq_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_VV_H1_WW_lnuqq_TH13::computeThValue()
{
    return log10(1.);
}



log10_ggVV_H1_WW_lnulnu_TH13::log10_ggVV_H1_WW_lnulnu_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_ggVV_H1_WW_lnulnu_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->ggVV_H_WW_lnulnu_TH13);
}



log10_pp_H1_WW_lnuqq_TH13::log10_pp_H1_WW_lnuqq_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H1_WW_lnuqq_TH13::computeThValue()
{
    return log10(1.);
}



log10_VV_H5_WW_TH8::log10_VV_H5_WW_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_VV_H5_WW_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->VV_H5_WW_TH8);
}



log10_VV_H5_WW_enumunu_TH13::log10_VV_H5_WW_enumunu_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_VV_H5_WW_enumunu_TH13::computeThValue()
{
    return log10(1.);
}



log10_VV_H5_WW_lnuqq_TH13::log10_VV_H5_WW_lnuqq_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_VV_H5_WW_lnuqq_TH13::computeThValue()
{
    return log10(1.);
}



log10_ggVV_H5_WW_lnulnu_TH13::log10_ggVV_H5_WW_lnulnu_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_ggVV_H5_WW_lnulnu_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->ggVV_H5_WW_lnulnu_TH13);
}



log10_pp_H5_WW_lnuqq_TH13::log10_pp_H5_WW_lnuqq_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H5_WW_lnuqq_TH13::computeThValue()
{
    return log10(1.);
}



log10_pp_H1_VV_TH8::log10_pp_H1_VV_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H1_VV_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H_VV_TH8);
}



log10_mu_pp_H1_VV_TH8::log10_mu_pp_H1_VV_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_mu_pp_H1_VV_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->mu_pp_H_VV_TH8);
}



log10_pp_H1_VV_TH13::log10_pp_H1_VV_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H1_VV_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H_VV_TH13);
}



log10_pp_H5_VV_TH8::log10_pp_H5_VV_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H5_VV_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H5_VV_TH8);
}



log10_mu_pp_H5_VV_TH8::log10_mu_pp_H5_VV_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_mu_pp_H5_VV_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->mu_pp_H5_VV_TH8);
}



log10_pp_H5_VV_TH13::log10_pp_H5_VV_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H5_VV_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H5_VV_TH13);
}



log10_gg_H1_hh_TH8::log10_gg_H1_hh_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H1_hh_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->gg_H_hh_TH8);
}



log10_pp_H1_hh_bbbb_TH8::log10_pp_H1_hh_bbbb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H1_hh_bbbb_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H_hh_bbbb_TH8);
}



log10_pp_H1_hh_gagabb_TH8::log10_pp_H1_hh_gagabb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H1_hh_gagabb_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H_hh_gagabb_TH8);
}



log10_gg_H1_hh_bbtautau_TH8::log10_gg_H1_hh_bbtautau_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H1_hh_bbtautau_TH8::computeThValue()
{
    return log10(1.);
}



log10_pp_H1_hh_bbtautau_TH8::log10_pp_H1_hh_bbtautau_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H1_hh_bbtautau_TH8::computeThValue()
{
    return log10(1.);
}



log10_pp_H1_hh_bbbb_TH13::log10_pp_H1_hh_bbbb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H1_hh_bbbb_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H_hh_bbbb_TH13);
}



log10_gg_H1_hh_bbbb_TH13::log10_gg_H1_hh_bbbb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H1_hh_bbbb_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->gg_H_hh_bbbb_TH13);
}



log10_pp_H1_hh_gagabb_TH13::log10_pp_H1_hh_gagabb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H1_hh_gagabb_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H_hh_gagabb_TH13);
}



log10_pp_H1_hh_bbtautau_TH13::log10_pp_H1_hh_bbtautau_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H1_hh_bbtautau_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H_hh_bbtautau_TH13);
}



log10_pp_H1_hh_bblnulnu_TH13::log10_pp_H1_hh_bblnulnu_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H1_hh_bblnulnu_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H_hh_bblnulnu_TH13);
}



log10_gg_H1_hh_gagaWW_TH13::log10_gg_H1_hh_gagaWW_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H1_hh_gagaWW_TH13::computeThValue()
{
    return log10(1.);
}



log10_pp_H5_hh_bbbb_TH8::log10_pp_H5_hh_bbbb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H5_hh_bbbb_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H5_hh_bbbb_TH8);
}



log10_pp_H5_hh_gagabb_TH8::log10_pp_H5_hh_gagabb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H5_hh_gagabb_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H5_hh_gagabb_TH8);
}



log10_pp_H5_hh_bbtautau_TH8::log10_pp_H5_hh_bbtautau_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H5_hh_bbtautau_TH8::computeThValue()
{
    return log10(1.);
}



log10_pp_H5_hh_bbbb_TH13::log10_pp_H5_hh_bbbb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H5_hh_bbbb_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H5_hh_bbbb_TH13);
}



log10_pp_H5_hh_gagabb_TH13::log10_pp_H5_hh_gagabb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H5_hh_gagabb_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H5_hh_gagabb_TH13);
}



log10_pp_H5_hh_bbtautau_TH13::log10_pp_H5_hh_bbtautau_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H5_hh_bbtautau_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H5_hh_bbtautau_TH13);
}



log10_pp_H5_hh_bblnulnu_TH13::log10_pp_H5_hh_bblnulnu_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H5_hh_bblnulnu_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H5_hh_bblnulnu_TH13);
}



log10_gg_H3_hZ_bbZ_TH8::log10_gg_H3_hZ_bbZ_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H3_hZ_bbZ_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->gg_A_hZ_bbZ_TH8);
}



log10_gg_H3_hZ_bbll_TH8::log10_gg_H3_hZ_bbll_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H3_hZ_bbll_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->gg_A_hZ_bbll_TH8);
}



log10_gg_H3_hZ_tautauZ_TH8::log10_gg_H3_hZ_tautauZ_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H3_hZ_tautauZ_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->gg_A_hZ_tautauZ_TH8);
}



log10_gg_H3_hZ_tautaull_TH8::log10_gg_H3_hZ_tautaull_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H3_hZ_tautaull_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->gg_A_hZ_tautaull_TH8);
}



log10_gg_H3_hZ_bbZ_TH13::log10_gg_H3_hZ_bbZ_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_gg_H3_hZ_bbZ_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->gg_A_hZ_bbZ_TH13);
}



log10_bb_H3_hZ_bbZ_TH13::log10_bb_H3_hZ_bbZ_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_bb_H3_hZ_bbZ_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->bb_A_hZ_bbZ_TH13);
}



log10_pp_H3_H1Z_bbll_TH8::log10_pp_H3_H1Z_bbll_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H3_H1Z_bbll_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_A_HZ_bbll_TH8);
}



log10_pp_H3_H5Z_bbll_TH8::log10_pp_H3_H5Z_bbll_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H3_H5Z_bbll_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_A_H5Z_bbll_TH8);
}



log10_pp_H1_H3Z_bbll_TH8::log10_pp_H1_H3Z_bbll_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H1_H3Z_bbll_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H_AZ_bbll_TH8);
}



log10_pp_H5_H3Z_bbll_TH8::log10_pp_H5_H3Z_bbll_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H5_H3Z_bbll_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H5_AZ_bbll_TH8);
}



log10_pp_H3pm_taunu_TH8::log10_pp_H3pm_taunu_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H3pm_taunu_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_Hpm_taunu_TH8);
}



log10_pp_H3p_taunu_TH8::log10_pp_H3p_taunu_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H3p_taunu_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_Hp_taunu_TH8);
}



log10_pp_H3pm_taunu_TH13::log10_pp_H3pm_taunu_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H3pm_taunu_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_Hpm_taunu_TH13);
}



log10_pp_H3pm_tb_TH8::log10_pp_H3pm_tb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H3pm_tb_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_Hpm_tb_TH8);
}



log10_pp_H3p_tb_TH8::log10_pp_H3p_tb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H3p_tb_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_Hp_tb_TH8);
}



log10_pp_H3p_tb_TH13::log10_pp_H3p_tb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H3p_tb_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_Hp_tb_TH13);
}



log10_WZ_H5pm_WZ_TH8::log10_WZ_H5pm_WZ_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_WZ_H5pm_WZ_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->WZ_H5pm_WZ_TH8);
}



log10_WZ_H5pm_WZ_TH13::log10_WZ_H5pm_WZ_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_WZ_H5pm_WZ_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->WZ_H5pm_WZ_TH13);
}



log10_pp_H5ppmmH5mmpp_TH8::log10_pp_H5ppmmH5mmpp_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H5ppmmH5mmpp_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H5ppmmH5mmpp_TH8);
}



log10_pp_H5ppmmH5mmpp_TH13::log10_pp_H5ppmmH5mmpp_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H5ppmmH5mmpp_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H5ppmmH5mmpp_TH13);
}



log10_pp_H5ppmmH5mmpp_WWWW_TH13::log10_pp_H5ppmmH5mmpp_WWWW_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H5ppmmH5mmpp_WWWW_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H5ppmmH5mmpp_WWWW_TH13);
}



log10_pp_H5ppmm_WW_TH8::log10_pp_H5ppmm_WW_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H5ppmm_WW_TH8::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H5ppmm_WW_TH8);
}



log10_pp_H5ppmm_WW_TH13::log10_pp_H5ppmm_WW_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double log10_pp_H5ppmm_WW_TH13::computeThValue()
{
    return log10(myGM.getMyGMCache()->pp_H5ppmm_WW_TH13);
}
