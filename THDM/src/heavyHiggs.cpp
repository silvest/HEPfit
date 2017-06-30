/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "heavyHiggs.h"
#include "StandardModel.h"



Hobs_ggF_H_tautau_ATLAS8::Hobs_ggF_H_tautau_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_H_tautau_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_H_tautau_ATLAS8;
}

Robs_ggF_H_tautau_ATLAS8::Robs_ggF_H_tautau_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_H_tautau_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_H_tautau_ATLAS8;
}



Hobs_ggF_H_tautau_CMS8::Hobs_ggF_H_tautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_H_tautau_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_H_tautau_CMS8;
}

Robs_ggF_H_tautau_CMS8::Robs_ggF_H_tautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_H_tautau_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_H_tautau_CMS8;
}



Hobs_bbF_H_tautau_ATLAS8::Hobs_bbF_H_tautau_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_bbF_H_tautau_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_bbF_H_tautau_ATLAS8;
}

Robs_bbF_H_tautau_ATLAS8::Robs_bbF_H_tautau_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_bbF_H_tautau_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_bbF_H_tautau_ATLAS8;
}



Hobs_bbF_H_tautau_CMS8::Hobs_bbF_H_tautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_bbF_H_tautau_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_bbF_H_tautau_CMS8;
}

Robs_bbF_H_tautau_CMS8::Robs_bbF_H_tautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_bbF_H_tautau_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_bbF_H_tautau_CMS8;
}



Hobs_pp_H_gaga_ATLAS8::Hobs_pp_H_gaga_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_H_gaga_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_H_gaga_ATLAS8;
}

Robs_pp_H_gaga_ATLAS8::Robs_pp_H_gaga_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_H_gaga_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_H_gaga_ATLAS8;
}



Hobs_ggF_H_gaga_CMS8::Hobs_ggF_H_gaga_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_H_gaga_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_H_gaga_CMS8;
}

Robs_ggF_H_gaga_CMS8::Robs_ggF_H_gaga_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_H_gaga_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_H_gaga_CMS8;
}



Hobs_pp_H_Zga_llga_ATLAS8::Hobs_pp_H_Zga_llga_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_H_Zga_llga_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_H_Zga_llga_ATLAS8;
}

Robs_pp_H_Zga_llga_ATLAS8::Robs_pp_H_Zga_llga_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_H_Zga_llga_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_H_Zga_llga_ATLAS8;
}



Hobs_pp_H_Zga_llga_CMS8::Hobs_pp_H_Zga_llga_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_H_Zga_llga_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_H_Zga_llga_CMS8;
}

Robs_pp_H_Zga_llga_CMS8::Robs_pp_H_Zga_llga_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_H_Zga_llga_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_H_Zga_llga_CMS8;
}



Hobs_mu_pp_H_VV_CMS8::Hobs_mu_pp_H_VV_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_mu_pp_H_VV_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_mu_pp_H_VV_CMS8;
}

Robs_mu_pp_H_VV_CMS8::Robs_mu_pp_H_VV_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_mu_pp_H_VV_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_mu_pp_H_VV_CMS8;
}



Hobs_ggF_H_ZZ_ATLAS8::Hobs_ggF_H_ZZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_H_ZZ_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_H_ZZ_ATLAS8;
}

Robs_ggF_H_ZZ_ATLAS8::Robs_ggF_H_ZZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_H_ZZ_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_H_ZZ_ATLAS8;
}



Hobs_VBF_H_ZZ_ATLAS8::Hobs_VBF_H_ZZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_VBF_H_ZZ_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_VBF_H_ZZ_ATLAS8;
}

Robs_VBF_H_ZZ_ATLAS8::Robs_VBF_H_ZZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_VBF_H_ZZ_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_VBF_H_ZZ_ATLAS8;
}



Hobs_ggF_H_WW_ATLAS8::Hobs_ggF_H_WW_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_H_WW_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_H_WW_ATLAS8;
}

Robs_ggF_H_WW_ATLAS8::Robs_ggF_H_WW_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_H_WW_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_H_WW_ATLAS8;
}



Hobs_VBF_H_WW_ATLAS8::Hobs_VBF_H_WW_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_VBF_H_WW_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_VBF_H_WW_ATLAS8;
}

Robs_VBF_H_WW_ATLAS8::Robs_VBF_H_WW_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_VBF_H_WW_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_VBF_H_WW_ATLAS8;
}



Hobs_ggF_H_hh_ATLAS8::Hobs_ggF_H_hh_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_H_hh_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_H_hh_ATLAS8;
}

Robs_ggF_H_hh_ATLAS8::Robs_ggF_H_hh_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_H_hh_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_H_hh_ATLAS8;
}



Hobs_pp_H_hh_CMS8::Hobs_pp_H_hh_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_H_hh_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_H_hh_CMS8;
}

Robs_pp_H_hh_CMS8::Robs_pp_H_hh_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_H_hh_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_H_hh_CMS8;
}



Hobs_ggF_H_hh_bbtautau_CMS8::Hobs_ggF_H_hh_bbtautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_H_hh_bbtautau_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_H_hh_bbtautau_CMS8;
}

Robs_ggF_H_hh_bbtautau_CMS8::Robs_ggF_H_hh_bbtautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_H_hh_bbtautau_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_H_hh_bbtautau_CMS8;
}



Hobs_pp_H_hh_bbbb_CMS8::Hobs_pp_H_hh_bbbb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_H_hh_bbbb_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_H_hh_bbbb_CMS8;
}

Robs_pp_H_hh_bbbb_CMS8::Robs_pp_H_hh_bbbb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_H_hh_bbbb_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_H_hh_bbbb_CMS8;
}



Hobs_pp_H_hh_gagabb_CMS8::Hobs_pp_H_hh_gagabb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_H_hh_gagabb_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_H_hh_gagabb_CMS8;
}

Robs_pp_H_hh_gagabb_CMS8::Robs_pp_H_hh_gagabb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_H_hh_gagabb_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_H_hh_gagabb_CMS8;
}



Hobs_ggF_H_tt_ATLAS8::Hobs_ggF_H_tt_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_H_tt_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_H_tt_ATLAS8;
}

Robs_ggF_H_tt_ATLAS8::Robs_ggF_H_tt_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_H_tt_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_H_tt_ATLAS8;
}



Hobs_bbF_H_bb_CMS8::Hobs_bbF_H_bb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_bbF_H_bb_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_bbF_H_bb_CMS8;
}

Robs_bbF_H_bb_CMS8::Robs_bbF_H_bb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_bbF_H_bb_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_bbF_H_bb_CMS8;
}



Hobs_ttF_H_tt_ATLAS13::Hobs_ttF_H_tt_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ttF_H_tt_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ttF_H_tt_ATLAS13;
}

Robs_ttF_H_tt_ATLAS13::Robs_ttF_H_tt_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ttF_H_tt_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ttF_H_tt_ATLAS13;
}



Hobs_bbF_H_tt_ATLAS13::Hobs_bbF_H_tt_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_bbF_H_tt_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_bbF_H_tt_ATLAS13;
}

Robs_bbF_H_tt_ATLAS13::Robs_bbF_H_tt_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_bbF_H_tt_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_bbF_H_tt_ATLAS13;
}



Hobs_ggF_H_tautau_ATLAS13::Hobs_ggF_H_tautau_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_H_tautau_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_H_tautau_ATLAS13;
}

Robs_ggF_H_tautau_ATLAS13::Robs_ggF_H_tautau_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_H_tautau_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_H_tautau_ATLAS13;
}



Hobs_bbF_H_tautau_ATLAS13::Hobs_bbF_H_tautau_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_bbF_H_tautau_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_bbF_H_tautau_ATLAS13;
}

Robs_bbF_H_tautau_ATLAS13::Robs_bbF_H_tautau_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_bbF_H_tautau_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_bbF_H_tautau_ATLAS13;
}



Hobs_ggF_H_tautau_CMS13::Hobs_ggF_H_tautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_H_tautau_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_H_tautau_CMS13;
}

Robs_ggF_H_tautau_CMS13::Robs_ggF_H_tautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_H_tautau_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_H_tautau_CMS13;
}



Hobs_bbF_H_tautau_CMS13::Hobs_bbF_H_tautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_bbF_H_tautau_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_bbF_H_tautau_CMS13;
}

Robs_bbF_H_tautau_CMS13::Robs_bbF_H_tautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_bbF_H_tautau_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_bbF_H_tautau_CMS13;
}



Hobs_pp_H_gaga_ATLAS13::Hobs_pp_H_gaga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_H_gaga_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_H_gaga_ATLAS13;
}

Robs_pp_H_gaga_ATLAS13::Robs_pp_H_gaga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_H_gaga_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_H_gaga_ATLAS13;
}



Hobs_ggF_H_gaga_CMS13::Hobs_ggF_H_gaga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_H_gaga_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_H_gaga_CMS13;
}

Robs_ggF_H_gaga_CMS13::Robs_ggF_H_gaga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_H_gaga_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_H_gaga_CMS13;
}



Hobs_pp_H_Zga_llga_ATLAS13::Hobs_pp_H_Zga_llga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_H_Zga_llga_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_H_Zga_ATLAS13;
}

Robs_pp_H_Zga_llga_ATLAS13::Robs_pp_H_Zga_llga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_H_Zga_llga_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_H_Zga_ATLAS13;
}



Hobs_pp_H_Zga_llga_CMS13::Hobs_pp_H_Zga_llga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_H_Zga_llga_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_H_Zga_llga_CMS13;
}

Robs_pp_H_Zga_llga_CMS13::Robs_pp_H_Zga_llga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_H_Zga_llga_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_H_Zga_llga_CMS13;
}



Hobs_pp_H_Zga_qqga_CMS13::Hobs_pp_H_Zga_qqga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_H_Zga_qqga_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_H_Zga_qqga_CMS13;
}

Robs_pp_H_Zga_qqga_CMS13::Robs_pp_H_Zga_qqga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_H_Zga_qqga_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_H_Zga_qqga_CMS13;
}



Hobs_ggF_H_Zga_CMS13::Hobs_ggF_H_Zga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_H_Zga_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_H_Zga_CMS13;
}

Robs_ggF_H_Zga_CMS13::Robs_ggF_H_Zga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_H_Zga_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_H_Zga_CMS13;
}



Hobs_ggF_H_ZZ_llnunu_ATLAS13::Hobs_ggF_H_ZZ_llnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_H_ZZ_llnunu_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_H_ZZ_llnunu_ATLAS13;
}

Robs_ggF_H_ZZ_llnunu_ATLAS13::Robs_ggF_H_ZZ_llnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_H_ZZ_llnunu_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_H_ZZ_llnunu_ATLAS13;
}



Hobs_ggF_H_ZZ_llll_ATLAS13::Hobs_ggF_H_ZZ_llll_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_H_ZZ_llll_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_H_ZZ_llll_ATLAS13;
}

Robs_ggF_H_ZZ_llll_ATLAS13::Robs_ggF_H_ZZ_llll_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_H_ZZ_llll_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_H_ZZ_llll_ATLAS13;
}



Hobs_VBF_H_ZZ_llll_ATLAS13::Hobs_VBF_H_ZZ_llll_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_VBF_H_ZZ_llll_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_VBF_H_ZZ_llll_ATLAS13;
}

Robs_VBF_H_ZZ_llll_ATLAS13::Robs_VBF_H_ZZ_llll_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_VBF_H_ZZ_llll_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_VBF_H_ZZ_llll_ATLAS13;
}



Hobs_pp_H_ZZ_llll_CMS13::Hobs_pp_H_ZZ_llll_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_H_ZZ_llll_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_H_ZZ_llll_CMS13;
}

Robs_pp_H_ZZ_llll_CMS13::Robs_pp_H_ZZ_llll_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_H_ZZ_llll_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_H_ZZ_llll_CMS13;
}



Hobs_VBF_VH_H_ZZ_llll_CMS13::Hobs_VBF_VH_H_ZZ_llll_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_VBF_VH_H_ZZ_llll_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_VBF_VH_H_ZZ_llll_CMS13;
}

Robs_VBF_VH_H_ZZ_llll_CMS13::Robs_VBF_VH_H_ZZ_llll_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_VBF_VH_H_ZZ_llll_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_VBF_VH_H_ZZ_llll_CMS13;
}



Hobs_ggF_H_ZZ_llqq_ATLAS13::Hobs_ggF_H_ZZ_llqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_H_ZZ_llqq_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_H_ZZ_llqq_ATLAS13;
}

Robs_ggF_H_ZZ_llqq_ATLAS13::Robs_ggF_H_ZZ_llqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_H_ZZ_llqq_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_H_ZZ_llqq_ATLAS13;
}



Hobs_VBF_H_ZZ_llqq_ATLAS13::Hobs_VBF_H_ZZ_llqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_VBF_H_ZZ_llqq_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_VBF_H_ZZ_llqq_ATLAS13;
}

Robs_VBF_H_ZZ_llqq_ATLAS13::Robs_VBF_H_ZZ_llqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_VBF_H_ZZ_llqq_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_VBF_H_ZZ_llqq_ATLAS13;
}



Hobs_ggF_H_ZZ_nunuqq_ATLAS13::Hobs_ggF_H_ZZ_nunuqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_H_ZZ_nunuqq_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_H_ZZ_nunuqq_ATLAS13;
}

Robs_ggF_H_ZZ_nunuqq_ATLAS13::Robs_ggF_H_ZZ_nunuqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_H_ZZ_nunuqq_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_H_ZZ_nunuqq_ATLAS13;
}



Hobs_pp_H_ZZ_llqq_CMS13::Hobs_pp_H_ZZ_llqq_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_H_ZZ_llqq_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_H_ZZ_llqq_CMS13;
}

Robs_pp_H_ZZ_llqq_CMS13::Robs_pp_H_ZZ_llqq_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_H_ZZ_llqq_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_H_ZZ_llqq_CMS13;
}



Hobs_ggF_H_WW_lnuqq_ATLAS13::Hobs_ggF_H_WW_lnuqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_H_WW_lnuqq_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_H_WW_lnuqq_ATLAS13;
}

Robs_ggF_H_WW_lnuqq_ATLAS13::Robs_ggF_H_WW_lnuqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_H_WW_lnuqq_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_H_WW_lnuqq_ATLAS13;
}



Hobs_ggF_H_WW_enumunu_ATLAS13::Hobs_ggF_H_WW_enumunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_H_WW_enumunu_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_H_WW_enumunu_ATLAS13;
}

Robs_ggF_H_WW_enumunu_ATLAS13::Robs_ggF_H_WW_enumunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_H_WW_enumunu_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_H_WW_enumunu_ATLAS13;
}



Hobs_VBF_H_WW_enumunu_ATLAS13::Hobs_VBF_H_WW_enumunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_VBF_H_WW_enumunu_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_VBF_H_WW_enumunu_ATLAS13;
}

Robs_VBF_H_WW_enumunu_ATLAS13::Robs_VBF_H_WW_enumunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_VBF_H_WW_enumunu_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_VBF_H_WW_enumunu_ATLAS13;
}



Hobs_ggF_VBF_H_WW_lnulnu_CMS13::Hobs_ggF_VBF_H_WW_lnulnu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_VBF_H_WW_lnulnu_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_VBF_H_WW_lnulnu_CMS13;
}

Robs_ggF_VBF_H_WW_lnulnu_CMS13::Robs_ggF_VBF_H_WW_lnulnu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_VBF_H_WW_lnulnu_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_VBF_H_WW_lnulnu_CMS13;
}



Hobs_pp_H_hh_bbgaga_ATLAS13::Hobs_pp_H_hh_bbgaga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_H_hh_bbgaga_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_H_hh_bbgaga_ATLAS13;
}

Robs_pp_H_hh_bbgaga_ATLAS13::Robs_pp_H_hh_bbgaga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_H_hh_bbgaga_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_H_hh_bbgaga_ATLAS13;
}



Hobs_pp_H_hh_bbgaga_CMS13::Hobs_pp_H_hh_bbgaga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_H_hh_bbgaga_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_H_hh_bbgaga_CMS13;
}

Robs_pp_H_hh_bbgaga_CMS13::Robs_pp_H_hh_bbgaga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_H_hh_bbgaga_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_H_hh_bbgaga_CMS13;
}



Hobs_pp_H_hh_bbbb_ATLAS13::Hobs_pp_H_hh_bbbb_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_H_hh_bbbb_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_H_hh_bbbb_ATLAS13;
}

Robs_pp_H_hh_bbbb_ATLAS13::Robs_pp_H_hh_bbbb_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_H_hh_bbbb_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_H_hh_bbbb_ATLAS13;
}



Hobs_pp_H_hh_bbbb_CMS13::Hobs_pp_H_hh_bbbb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_H_hh_bbbb_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_H_hh_bbbb_CMS13;
}

Robs_pp_H_hh_bbbb_CMS13::Robs_pp_H_hh_bbbb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_H_hh_bbbb_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_H_hh_bbbb_CMS13;
}



Hobs_ggF_H_hh_bbbb_CMS13::Hobs_ggF_H_hh_bbbb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_H_hh_bbbb_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_H_hh_bbbb_CMS13;
}

Robs_ggF_H_hh_bbbb_CMS13::Robs_ggF_H_hh_bbbb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_H_hh_bbbb_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_H_hh_bbbb_CMS13;
}



Hobs_ggF_H_hh_gagaWW_ATLAS13::Hobs_ggF_H_hh_gagaWW_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_H_hh_gagaWW_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_H_hh_gagaWW_ATLAS13;
}

Robs_ggF_H_hh_gagaWW_ATLAS13::Robs_ggF_H_hh_gagaWW_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_H_hh_gagaWW_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_H_hh_gagaWW_ATLAS13;
}



Hobs_pp_H_hh_bbtautau_CMS13::Hobs_pp_H_hh_bbtautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_H_hh_bbtautau_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_H_hh_bbtautau_CMS13;
}

Robs_pp_H_hh_bbtautau_CMS13::Robs_pp_H_hh_bbtautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_H_hh_bbtautau_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_H_hh_bbtautau_CMS13;
}



Hobs_pp_H_hh_bbtautau1_CMS13::Hobs_pp_H_hh_bbtautau1_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_H_hh_bbtautau1_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_H_hh_bbtautau1_CMS13;
}

Robs_pp_H_hh_bbtautau1_CMS13::Robs_pp_H_hh_bbtautau1_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_H_hh_bbtautau1_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_H_hh_bbtautau1_CMS13;
}



Hobs_pp_H_hh_bblnulnu_CMS13::Hobs_pp_H_hh_bblnulnu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_H_hh_bblnulnu_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_H_hh_bblnulnu_CMS13;
}

Robs_pp_H_hh_bblnulnu_CMS13::Robs_pp_H_hh_bblnulnu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_H_hh_bblnulnu_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_H_hh_bblnulnu_CMS13;
}



Hobs_pp_H_hh_bbVV_CMS13::Hobs_pp_H_hh_bbVV_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_H_hh_bbVV_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_H_hh_bbVV_CMS13;
}

Robs_pp_H_hh_bbVV_CMS13::Robs_pp_H_hh_bbVV_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_H_hh_bbVV_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_H_hh_bbVV_CMS13;
}



Hobs_pp_H_bb_CMS13::Hobs_pp_H_bb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_H_bb_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_H_bb_CMS13;
}

Robs_pp_H_bb_CMS13::Robs_pp_H_bb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_H_bb_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_H_bb_CMS13;
}



log10_ggF_H_tautau_TH8::log10_ggF_H_tautau_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_ggF_H_tautau_TH8::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->ggF_H_tautau_TH8);
}



log10_bbF_H_tautau_TH8::log10_bbF_H_tautau_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_bbF_H_tautau_TH8::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->bbF_H_tautau_TH8);
}



log10_pp_H_gaga_TH8::log10_pp_H_gaga_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_pp_H_gaga_TH8::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->pp_H_gaga_TH8);
}



log10_ggF_H_gaga_TH8::log10_ggF_H_gaga_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_ggF_H_gaga_TH8::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->ggF_H_gaga_TH8);
}



log10_pp_H_Zga_llga_TH8::log10_pp_H_Zga_llga_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_pp_H_Zga_llga_TH8::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->pp_H_Zga_llga_TH8);
}



log10_mu_pp_H_VV_TH8::log10_mu_pp_H_VV_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_mu_pp_H_VV_TH8::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->mu_pp_H_VV_TH8);
}



log10_ggF_H_ZZ_TH8::log10_ggF_H_ZZ_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_ggF_H_ZZ_TH8::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->ggF_H_ZZ_TH8);
}



log10_VBF_H_ZZ_TH8::log10_VBF_H_ZZ_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_VBF_H_ZZ_TH8::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->VBF_H_ZZ_TH8);
}



log10_ggF_H_WW_TH8::log10_ggF_H_WW_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_ggF_H_WW_TH8::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->ggF_H_WW_TH8);
}



log10_VBF_H_WW_TH8::log10_VBF_H_WW_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_VBF_H_WW_TH8::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->VBF_H_WW_TH8);
}



log10_ggF_H_hh_TH8::log10_ggF_H_hh_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_ggF_H_hh_TH8::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->ggF_H_hh_TH8);
}



log10_pp_H_hh_TH8::log10_pp_H_hh_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_pp_H_hh_TH8::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->pp_H_hh_TH8);
}



log10_ggF_H_hh_bbtautau_TH8::log10_ggF_H_hh_bbtautau_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_ggF_H_hh_bbtautau_TH8::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->ggF_H_hh_bbtautau_TH8);
}



log10_pp_H_hh_bbbb_TH8::log10_pp_H_hh_bbbb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_pp_H_hh_bbbb_TH8::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->pp_H_hh_bbbb_TH8);
}



log10_pp_H_hh_gagabb_TH8::log10_pp_H_hh_gagabb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_pp_H_hh_gagabb_TH8::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->pp_H_hh_gagabb_TH8);
}



log10_ggF_H_tt_TH8::log10_ggF_H_tt_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_ggF_H_tt_TH8::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->ggF_H_tt_TH8);
}



log10_bbF_H_bb_TH8::log10_bbF_H_bb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_bbF_H_bb_TH8::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->bbF_H_bb_TH8);
}



log10_ggF_H_tautau_TH13::log10_ggF_H_tautau_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_ggF_H_tautau_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->ggF_H_tautau_TH13);
}



log10_bbF_H_tautau_TH13::log10_bbF_H_tautau_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_bbF_H_tautau_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->bbF_H_tautau_TH13);
}



log10_pp_H_gaga_TH13::log10_pp_H_gaga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_pp_H_gaga_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->pp_H_gaga_TH13);
}



log10_ggF_H_gaga_TH13::log10_ggF_H_gaga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_ggF_H_gaga_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->ggF_H_gaga_TH13);
}



log10_pp_H_Zga_TH13::log10_pp_H_Zga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_pp_H_Zga_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->pp_H_Zga_TH13);
}



log10_ggF_H_Zga_TH13::log10_ggF_H_Zga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_ggF_H_Zga_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->ggF_H_Zga_TH13);
}



log10_pp_H_ZZ_TH13::log10_pp_H_ZZ_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_pp_H_ZZ_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->pp_H_ZZ_TH13);
}



log10_ggF_H_ZZ_TH13::log10_ggF_H_ZZ_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_ggF_H_ZZ_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->ggF_H_ZZ_TH13);
}



log10_VBF_H_ZZ_TH13::log10_VBF_H_ZZ_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_VBF_H_ZZ_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->VBF_H_ZZ_TH13);
}



log10_ggF_H_ZZ_llll_TH13::log10_ggF_H_ZZ_llll_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_ggF_H_ZZ_llll_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->ggF_H_ZZ_llll_TH13);
}



log10_VBF_H_ZZ_llll_TH13::log10_VBF_H_ZZ_llll_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_VBF_H_ZZ_llll_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->VBF_H_ZZ_llll_TH13);
}



log10_pp_H_ZZ_llll_TH13::log10_pp_H_ZZ_llll_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_pp_H_ZZ_llll_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->pp_H_ZZ_llll_TH13);
}



log10_VBF_VH_H_ZZ_llll_TH13::log10_VBF_VH_H_ZZ_llll_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_VBF_VH_H_ZZ_llll_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->VBF_VH_H_ZZ_llll_TH13);
}



log10_ggF_H_WW_TH13::log10_ggF_H_WW_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_ggF_H_WW_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->ggF_H_WW_TH13);
}



log10_VBF_H_WW_TH13::log10_VBF_H_WW_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_VBF_H_WW_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->VBF_H_WW_TH13);
}



log10_ggF_VBF_H_WW_lnulnu_TH13::log10_ggF_VBF_H_WW_lnulnu_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_ggF_VBF_H_WW_lnulnu_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->ggF_VBF_H_WW_lnulnu_TH13);
}



log10_ggF_H_hh_TH13::log10_ggF_H_hh_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_ggF_H_hh_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->ggF_H_hh_TH13);
}



log10_pp_H_hh_TH13::log10_pp_H_hh_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_pp_H_hh_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->pp_H_hh_TH13);
}



log10_pp_H_hh_bbbb_TH13::log10_pp_H_hh_bbbb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_pp_H_hh_bbbb_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->pp_H_hh_bbbb_TH13);
}



log10_ggF_H_hh_bbbb_TH13::log10_ggF_H_hh_bbbb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_ggF_H_hh_bbbb_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->ggF_H_hh_bbbb_TH13);
}



log10_pp_H_hh_gagabb_TH13::log10_pp_H_hh_gagabb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_pp_H_hh_gagabb_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->pp_H_hh_gagabb_TH13);
}



log10_pp_H_hh_bbtautau_TH13::log10_pp_H_hh_bbtautau_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_pp_H_hh_bbtautau_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->pp_H_hh_bbtautau_TH13);
}



log10_pp_H_hh_bblnulnu_TH13::log10_pp_H_hh_bblnulnu_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_pp_H_hh_bblnulnu_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->pp_H_hh_bblnulnu_TH13);
}



log10_pp_H_hh_bbVV_TH13::log10_pp_H_hh_bbVV_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_pp_H_hh_bbVV_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->pp_H_hh_bbVV_TH13);
}



log10_tt_H_tt_TH13::log10_tt_H_tt_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_tt_H_tt_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->ttF_H_tt_TH13);
}



log10_bb_H_tt_TH13::log10_bb_H_tt_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_bb_H_tt_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->bbF_H_tt_TH13);
}



log10_pp_H_bb_TH13::log10_pp_H_bb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_pp_H_bb_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->pp_H_bb_TH13);
}



Gamma_HH_THDM::Gamma_HH_THDM(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Gamma_HH_THDM::computeThValue()
{
    return myTHDM.getMyTHDMCache()->GammaHtot;
}



//rHH_gaga_THDM::rHH_gaga_THDM(const StandardModel& SM_i)
//: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
//{}
//
//double rHH_gaga_THDM::computeThValue()
//{
//    return myTHDM.getMyTHDMCache()->rHH_gaga;
//}



rHH_gg_THDM::rHH_gg_THDM(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double rHH_gg_THDM::computeThValue()
{
    return myTHDM.getMyTHDMCache()->rHH_gg;
}



BR_HH_hh_THDM::BR_HH_hh_THDM(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double BR_HH_hh_THDM::computeThValue()
{
    return myTHDM.getMyTHDMCache()->Br_Htohh;
}



BR_HH_AA_THDM::BR_HH_AA_THDM(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double BR_HH_AA_THDM::computeThValue()
{
    return myTHDM.getMyTHDMCache()->Br_HtoAA;
}



BR_HH_HpHm_THDM::BR_HH_HpHm_THDM(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double BR_HH_HpHm_THDM::computeThValue()
{
    return myTHDM.getMyTHDMCache()->Br_HtoHpHm;
}



BR_HH_AZ_THDM::BR_HH_AZ_THDM(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double BR_HH_AZ_THDM::computeThValue()
{
    return myTHDM.getMyTHDMCache()->Br_HtoAZ;
}



BR_HH_HpW_THDM::BR_HH_HpW_THDM(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double BR_HH_HpW_THDM::computeThValue()
{
    return myTHDM.getMyTHDMCache()->Br_HtoHpW;
}









//LIMITmLIMEST::LIMITmLIMEST(const StandardModel& SM_i)
//: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
//{}
//
//double LIMITmLIMEST::computeThValue()
//{
//    return (myTHDM.getMyTHDMCache()->LIMIT_ggF_H_gaga_CMS8) - (myTHDM.getMyTHDMCache()->LIMEST_ggF_H_gaga_CMS8);
//}
//
//DEVIATIONoBANDSIZE::DEVIATIONoBANDSIZE(const StandardModel& SM_i)
//: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
//{}
//
//double DEVIATIONoBANDSIZE::computeThValue()
//{
//    return (myTHDM.getMyTHDMCache()->DEVIATION_ggF_H_gaga_CMS8) / (myTHDM.getMyTHDMCache()->BANDSIZE_ggF_H_gaga_CMS8);
//}
