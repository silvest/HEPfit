/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMHeavyHiggs.h"
#include "StandardModel.h"


Hobs_tt_phi2_tt_ATLAS13::Hobs_tt_phi2_tt_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_tt_phi2_tt_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_tt_phi2_tt_ATLAS13;
}

Hobs_bb_phi2_tt_ATLAS13::Hobs_bb_phi2_tt_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_bb_phi2_tt_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_bb_phi2_tt_ATLAS13;
}

Hobs_bb_phi2_bb_CMS8::Hobs_bb_phi2_bb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_bb_phi2_bb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_bb_phi2_bb_CMS8;
}

Hobs_gg_phi2_bb_CMS8::Hobs_gg_phi2_bb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi2_bb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_bb_CMS8;
}

Hobs_pp_phi2_bb_CMS13::Hobs_pp_phi2_bb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_bb_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_bb_CMS13;
}


Hobs_bb_phi2_bb_CMS13::Hobs_bb_phi2_bb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_bb_phi2_bb_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_bb_phi2_bb_CMS13;
}

Hobs_gg_phi2_tautau_ATLAS8::Hobs_gg_phi2_tautau_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Hobs_gg_phi2_tautau_ATLAS8::computeThValue() 
{
    return myGTHDM.getMyGTHDMCache() -> THoEX_gg_phi2_tautau_ATLAS8;
}

Hobs_bb_phi2_tautau_ATLAS8::Hobs_bb_phi2_tautau_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Hobs_bb_phi2_tautau_ATLAS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->THoEX_bb_phi2_tautau_ATLAS8;
}


Hobs_gg_phi2_tautau_CMS8::Hobs_gg_phi2_tautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Hobs_gg_phi2_tautau_CMS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_tautau_CMS8;
}

Hobs_bb_phi2_tautau_CMS8::Hobs_bb_phi2_tautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Hobs_bb_phi2_tautau_CMS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->THoEX_bb_phi2_tautau_CMS8;
}


Hobs_gg_phi2_tautau_ATLAS13::Hobs_gg_phi2_tautau_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi2_tautau_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_tautau_ATLAS13;
}

Hobs_bb_phi2_tautau_ATLAS13::Hobs_bb_phi2_tautau_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_bb_phi2_tautau_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_bb_phi2_tautau_ATLAS13;
}

Hobs_gg_phi2_tautau_CMS13::Hobs_gg_phi2_tautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi2_tautau_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_tautau_CMS13;
}



Hobs_bb_phi2_tautau_CMS13::Hobs_bb_phi2_tautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_bb_phi2_tautau_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_bb_phi2_tautau_CMS13;
}

Hobs_gg_phi2_gaga_ATLAS8::Hobs_gg_phi2_gaga_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Hobs_gg_phi2_gaga_ATLAS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_gaga_ATLAS8;
}

Hobs_pp_phi2_gaga_ATLAS13::Hobs_pp_phi2_gaga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_gaga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_gaga_ATLAS13;
}


Hobs_gg_phi2_gaga_CMS13::Hobs_gg_phi2_gaga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi2_gaga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_gaga_CMS13;
}

Hobs_pp_phi2_Zga_llga_ATLAS8::Hobs_pp_phi2_Zga_llga_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_Zga_llga_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_Zga_llga_ATLAS8;
}

Hobs_pp_phi2_Zga_llga_CMS8::Hobs_pp_phi2_Zga_llga_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_Zga_llga_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_Zga_llga_CMS8;
}

Hobs_gg_phi2_Zga_llga_ATLAS13::Hobs_gg_phi2_Zga_llga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi2_Zga_llga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_Zga_llga_ATLAS13;
}

Hobs_gg_phi2_Zga_qqga_ATLAS13::Hobs_gg_phi2_Zga_qqga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi2_Zga_qqga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_Zga_qqga_ATLAS13;
}

Hobs_gg_phi2_Zga_CMS13::Hobs_gg_phi2_Zga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi2_Zga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_Zga_CMS13;
}

Hobs_gg_phi2_ZZ_ATLAS8::Hobs_gg_phi2_ZZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi2_ZZ_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_ZZ_ATLAS8;
}


Hobs_VV_phi2_ZZ_ATLAS8::Hobs_VV_phi2_ZZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VV_phi2_ZZ_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VV_phi2_ZZ_ATLAS8;
}

Hobs_gg_phi2_ZZ_llllnunu_ATLAS13::Hobs_gg_phi2_ZZ_llllnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi2_ZZ_llllnunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_ZZ_llllnunu_ATLAS13;
}

Hobs_VV_phi2_ZZ_llllnunu_ATLAS13::Hobs_VV_phi2_ZZ_llllnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VV_phi2_ZZ_llllnunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VV_phi2_ZZ_llllnunu_ATLAS13;
}


Hobs_gg_phi2_ZZ_llnunuqq_ATLAS13::Hobs_gg_phi2_ZZ_llnunuqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi2_ZZ_llnunuqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_ZZ_qqllnunu_ATLAS13;
}

Hobs_VV_phi2_ZZ_llnunuqq_ATLAS13::Hobs_VV_phi2_ZZ_llnunuqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VV_phi2_ZZ_llnunuqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VV_phi2_ZZ_qqllnunu_ATLAS13;
}

Hobs_pp_phi2_ZZ_llqqnunull_CMS13::Hobs_pp_phi2_ZZ_llqqnunull_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_ZZ_llqqnunull_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_ZZ_llqqnunull_CMS13;
}

Hobs_pp_phi2_ZZ_qqnunu_CMS13::Hobs_pp_phi2_ZZ_qqnunu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_ZZ_qqnunu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_ZZ_qqnunu_CMS13;
}

Hobs_gg_phi2_WW_ATLAS8::Hobs_gg_phi2_WW_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi2_WW_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_WW_ATLAS8;
}



Hobs_VV_phi2_WW_ATLAS8::Hobs_VV_phi2_WW_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VV_phi2_WW_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VV_phi2_WW_ATLAS8;
}



Hobs_gg_phi2_WW_enumunu_ATLAS13::Hobs_gg_phi2_WW_enumunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi2_WW_enumunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_WW_enumunu_ATLAS13;
}

Hobs_VV_phi2_WW_enumunu_ATLAS13::Hobs_VV_phi2_WW_enumunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VV_phi2_WW_enumunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VV_phi2_WW_enumunu_ATLAS13;
}


Hobs_ggVV_phi2_WW_lnulnu_CMS13::Hobs_ggVV_phi2_WW_lnulnu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggVV_phi2_WW_lnulnu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggVV_phi2_WW_lnulnu_CMS13;
}


Hobs_gg_phi2_WW_lnuqq_ATLAS13::Hobs_gg_phi2_WW_lnuqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi2_WW_lnuqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_WW_lnuqq_ATLAS13;
}


Hobs_VV_phi2_WW_lnuqq_ATLAS13::Hobs_VV_phi2_WW_lnuqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VV_phi2_WW_lnuqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VV_phi2_WW_lnuqq_ATLAS13;
}

Hobs_pp_phi2_WW_lnuqq_CMS13::Hobs_pp_phi2_WW_lnuqq_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_WW_lnuqq_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_WW_lnuqq_CMS13;
}

Hobs_mu_pp_phi2_VV_CMS8::Hobs_mu_pp_phi2_VV_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_mu_pp_phi2_VV_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_mu_pp_phi2_VV_CMS8;
}

Hobs_pp_phi2_VV_qqqq_ATLAS13::Hobs_pp_phi2_VV_qqqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_VV_qqqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_VV_qqqq_ATLAS13;
}

Hobs_gg_phi2_phi1phi1_ATLAS8::Hobs_gg_phi2_phi1phi1_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi2_phi1phi1_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_phi1phi1_ATLAS8;
}

Hobs_pp_phi2_phi1phi1_bbbb_CMS8::Hobs_pp_phi2_phi1phi1_bbbb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi1phi1_bbbb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi1phi1_bbbb_CMS8;
}

Hobs_pp_phi2_phi1phi1_bbgaga_CMS8::Hobs_pp_phi2_phi1phi1_bbgaga_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi1phi1_bbgaga_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi1phi1_bbgaga_CMS8;
}

Hobs_gg_phi2_phi1phi1_bbtautau_CMS8::Hobs_gg_phi2_phi1phi1_bbtautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi2_phi1phi1_bbtautau_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_phi1phi1_bbtautau_CMS8;
}

Hobs_pp_phi2_phi1phi1_bbtautau_CMS8::Hobs_pp_phi2_phi1phi1_bbtautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi1phi1_bbtautau_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi1phi1_bbtautau_CMS8;
}

Hobs_pp_phi2_phi1phi1_bbbb_ATLAS13::Hobs_pp_phi2_phi1phi1_bbbb_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi1phi1_bbbb_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi1phi1_bbbb_ATLAS13;
}

Hobs_pp_phi2_phi1phi1_bbbb_1_CMS13::Hobs_pp_phi2_phi1phi1_bbbb_1_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi1phi1_bbbb_1_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi1phi1_bbbb_1_CMS13;
}

Hobs_pp_phi2_phi1phi1_bbbb_2_CMS13::Hobs_pp_phi2_phi1phi1_bbbb_2_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi1phi1_bbbb_2_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi1phi1_bbbb_2_CMS13;
}

Hobs_pp_phi2_phi1phi1_bbgaga_ATLAS13::Hobs_pp_phi2_phi1phi1_bbgaga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi1phi1_bbgaga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi1phi1_bbgaga_ATLAS13;
}

Hobs_pp_phi2_phi1phi1_bbgaga_CMS13::Hobs_pp_phi2_phi1phi1_bbgaga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi1phi1_bbgaga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi1phi1_bbgaga_CMS13;
}

Hobs_pp_phi2_phi1phi1_bbtautau_ATLAS13::Hobs_pp_phi2_phi1phi1_bbtautau_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi1phi1_bbtautau_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi1phi1_bbtautau_ATLAS13;
}

Hobs_pp_phi2_phi1phi1_bbtautau_1_CMS13::Hobs_pp_phi2_phi1phi1_bbtautau_1_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi1phi1_bbtautau_1_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi1phi1_bbtautau_1_CMS13;
}

Hobs_pp_phi2_phi1phi1_bbtautau_2_CMS13::Hobs_pp_phi2_phi1phi1_bbtautau_2_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi1phi1_bbtautau_2_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi1phi1_bbtautau_2_CMS13;
}

Hobs_pp_phi2_phi1phi1_bbVV_CMS13::Hobs_pp_phi2_phi1phi1_bbVV_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi1phi1_bbVV_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi1phi1_bbVV_CMS13;
}


Hobs_pp_phi2_phi1phi1_bbWW_ATLAS13::Hobs_pp_phi2_phi1phi1_bbWW_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi1phi1_bbWW_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi1phi1_bbWW_ATLAS13;
}

Hobs_gg_phi2_phi1phi1_gagaWW_ATLAS13::Hobs_gg_phi2_phi1phi1_gagaWW_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi2_phi1phi1_gagaWW_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_phi1phi1_gagaWW_ATLAS13;
}

Hobs_gg_phi2_phi1Z_bbZ_ATLAS8::Hobs_gg_phi2_phi1Z_bbZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi2_phi1Z_bbZ_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_phi1Z_bbZ_ATLAS8;
}

Hobs_gg_phi2_phi1Z_bbll_CMS8::Hobs_gg_phi2_phi1Z_bbll_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi2_phi1Z_bbll_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_phi1Z_bbll_CMS8;
}

Hobs_gg_phi2_phi1Z_tautauZ_ATLAS8::Hobs_gg_phi2_phi1Z_tautauZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi2_phi1Z_tautauZ_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_phi1Z_tautauZ_ATLAS8;
}

Hobs_gg_phi2_phi1Z_tautaull_CMS8::Hobs_gg_phi2_phi1Z_tautaull_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi2_phi1Z_tautaull_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_phi1Z_tautaull_CMS8;
}

Hobs_gg_phi2_phi1Z_bbZ_ATLAS13::Hobs_gg_phi2_phi1Z_bbZ_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi2_phi1Z_bbZ_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_phi1Z_bbZ_ATLAS13;
}

Hobs_gg_phi2_phi1Z_bbZ_1_CMS13::Hobs_gg_phi2_phi1Z_bbZ_1_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi2_phi1Z_bbZ_1_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_phi1Z_bbZ_1_CMS13;
}

Hobs_gg_phi2_phi1Z_bbZ_2_CMS13::Hobs_gg_phi2_phi1Z_bbZ_2_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi2_phi1Z_bbZ_2_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_phi1Z_bbZ_2_CMS13;
}

Hobs_bb_phi2_phi1Z_bbZ_ATLAS13::Hobs_bb_phi2_phi1Z_bbZ_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_bb_phi2_phi1Z_bbZ_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_phi1Z_bbZ_ATLAS13;
}

Hobs_bb_phi2_phi1Z_bbZ_1_CMS13::Hobs_bb_phi2_phi1Z_bbZ_1_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_bb_phi2_phi1Z_bbZ_1_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_phi1Z_bbZ_1_CMS13;
}

Hobs_bb_phi2_phi1Z_bbZ_2_CMS13::Hobs_bb_phi2_phi1Z_bbZ_2_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_bb_phi2_phi1Z_bbZ_2_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_phi1Z_bbZ_2_CMS13;
}



Hobs_tt_phi3_tt_ATLAS13::Hobs_tt_phi3_tt_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_tt_phi3_tt_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_tt_phi3_tt_ATLAS13;
}

Hobs_bb_phi3_tt_ATLAS13::Hobs_bb_phi3_tt_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_bb_phi3_tt_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_bb_phi3_tt_ATLAS13;
}

Hobs_bb_phi3_bb_CMS8::Hobs_bb_phi3_bb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_bb_phi3_bb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_bb_phi3_bb_CMS8;
}

Hobs_gg_phi3_bb_CMS8::Hobs_gg_phi3_bb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi3_bb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_bb_CMS8;
}

Hobs_pp_phi3_bb_CMS13::Hobs_pp_phi3_bb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_bb_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_bb_CMS13;
}


Hobs_bb_phi3_bb_CMS13::Hobs_bb_phi3_bb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_bb_phi3_bb_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_bb_phi3_bb_CMS13;
}

Hobs_gg_phi3_tautau_ATLAS8::Hobs_gg_phi3_tautau_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Hobs_gg_phi3_tautau_ATLAS8::computeThValue() 
{
    return myGTHDM.getMyGTHDMCache() -> THoEX_gg_phi3_tautau_ATLAS8;
}

Hobs_bb_phi3_tautau_ATLAS8::Hobs_bb_phi3_tautau_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Hobs_bb_phi3_tautau_ATLAS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->THoEX_bb_phi3_tautau_ATLAS8;
}


Hobs_gg_phi3_tautau_CMS8::Hobs_gg_phi3_tautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Hobs_gg_phi3_tautau_CMS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_tautau_CMS8;
}

Hobs_bb_phi3_tautau_CMS8::Hobs_bb_phi3_tautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Hobs_bb_phi3_tautau_CMS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->THoEX_bb_phi3_tautau_CMS8;
}


Hobs_gg_phi3_tautau_ATLAS13::Hobs_gg_phi3_tautau_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi3_tautau_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_tautau_ATLAS13;
}

Hobs_bb_phi3_tautau_ATLAS13::Hobs_bb_phi3_tautau_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_bb_phi3_tautau_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_bb_phi3_tautau_ATLAS13;
}

Hobs_gg_phi3_tautau_CMS13::Hobs_gg_phi3_tautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi3_tautau_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_tautau_CMS13;
}

Hobs_bb_phi3_tautau_CMS13::Hobs_bb_phi3_tautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_bb_phi3_tautau_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_bb_phi3_tautau_CMS13;
}

Hobs_gg_phi3_gaga_ATLAS8::Hobs_gg_phi3_gaga_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Hobs_gg_phi3_gaga_ATLAS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_gaga_ATLAS8;
}

Hobs_pp_phi3_gaga_ATLAS13::Hobs_pp_phi3_gaga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_gaga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_gaga_ATLAS13;
}

Hobs_gg_phi3_gaga_CMS13::Hobs_gg_phi3_gaga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi3_gaga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_gaga_CMS13;
}

Hobs_pp_phi3_Zga_llga_ATLAS8::Hobs_pp_phi3_Zga_llga_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_Zga_llga_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_Zga_llga_ATLAS8;
}

Hobs_pp_phi3_Zga_llga_CMS8::Hobs_pp_phi3_Zga_llga_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_Zga_llga_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_Zga_llga_CMS8;
}

Hobs_gg_phi3_Zga_llga_ATLAS13::Hobs_gg_phi3_Zga_llga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi3_Zga_llga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_Zga_llga_ATLAS13;
}

Hobs_gg_phi3_Zga_qqga_ATLAS13::Hobs_gg_phi3_Zga_qqga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi3_Zga_qqga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_Zga_qqga_ATLAS13;
}

Hobs_gg_phi3_Zga_CMS13::Hobs_gg_phi3_Zga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi3_Zga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_Zga_CMS13;
}

Hobs_gg_phi3_ZZ_ATLAS8::Hobs_gg_phi3_ZZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi3_ZZ_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_ZZ_ATLAS8;
}


Hobs_VV_phi3_ZZ_ATLAS8::Hobs_VV_phi3_ZZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VV_phi3_ZZ_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VV_phi3_ZZ_ATLAS8;
}

Hobs_gg_phi3_ZZ_llllnunu_ATLAS13::Hobs_gg_phi3_ZZ_llllnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi3_ZZ_llllnunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_ZZ_llllnunu_ATLAS13;
}

Hobs_VV_phi3_ZZ_llllnunu_ATLAS13::Hobs_VV_phi3_ZZ_llllnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VV_phi3_ZZ_llllnunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VV_phi3_ZZ_llllnunu_ATLAS13;
}


Hobs_gg_phi3_ZZ_llnunuqq_ATLAS13::Hobs_gg_phi3_ZZ_llnunuqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi3_ZZ_llnunuqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_ZZ_qqllnunu_ATLAS13;
}

Hobs_VV_phi3_ZZ_llnunuqq_ATLAS13::Hobs_VV_phi3_ZZ_llnunuqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VV_phi3_ZZ_llnunuqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VV_phi3_ZZ_qqllnunu_ATLAS13;
}

Hobs_pp_phi3_ZZ_llqqnunull_CMS13::Hobs_pp_phi3_ZZ_llqqnunull_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_ZZ_llqqnunull_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_ZZ_llqqnunull_CMS13;
}

Hobs_pp_phi3_ZZ_qqnunu_CMS13::Hobs_pp_phi3_ZZ_qqnunu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_ZZ_qqnunu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_ZZ_qqnunu_CMS13;
}

Hobs_gg_phi3_WW_ATLAS8::Hobs_gg_phi3_WW_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi3_WW_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_WW_ATLAS8;
}



Hobs_VV_phi3_WW_ATLAS8::Hobs_VV_phi3_WW_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VV_phi3_WW_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VV_phi3_WW_ATLAS8;
}



Hobs_gg_phi3_WW_enumunu_ATLAS13::Hobs_gg_phi3_WW_enumunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi3_WW_enumunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_WW_enumunu_ATLAS13;
}

Hobs_VV_phi3_WW_enumunu_ATLAS13::Hobs_VV_phi3_WW_enumunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VV_phi3_WW_enumunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VV_phi3_WW_enumunu_ATLAS13;
}


Hobs_ggVV_phi3_WW_lnulnu_CMS13::Hobs_ggVV_phi3_WW_lnulnu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggVV_phi3_WW_lnulnu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggVV_phi3_WW_lnulnu_CMS13;
}


Hobs_gg_phi3_WW_lnuqq_ATLAS13::Hobs_gg_phi3_WW_lnuqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi3_WW_lnuqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_WW_lnuqq_ATLAS13;
}


Hobs_VV_phi3_WW_lnuqq_ATLAS13::Hobs_VV_phi3_WW_lnuqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VV_phi3_WW_lnuqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VV_phi3_WW_lnuqq_ATLAS13;
}

Hobs_pp_phi3_WW_lnuqq_CMS13::Hobs_pp_phi3_WW_lnuqq_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_WW_lnuqq_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_WW_lnuqq_CMS13;
}

Hobs_mu_pp_phi3_VV_CMS8::Hobs_mu_pp_phi3_VV_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_mu_pp_phi3_VV_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_mu_pp_phi3_VV_CMS8;
}

Hobs_pp_phi3_VV_qqqq_ATLAS13::Hobs_pp_phi3_VV_qqqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_VV_qqqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_VV_qqqq_ATLAS13;
}

Hobs_gg_phi3_phi1phi1_ATLAS8::Hobs_gg_phi3_phi1phi1_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi3_phi1phi1_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_phi1phi1_ATLAS8;
}

Hobs_pp_phi3_phi1phi1_bbbb_CMS8::Hobs_pp_phi3_phi1phi1_bbbb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi1_bbbb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi1_bbbb_CMS8;
}

Hobs_pp_phi3_phi1phi1_bbgaga_CMS8::Hobs_pp_phi3_phi1phi1_bbgaga_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi1_bbgaga_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi1_bbgaga_CMS8;
}

Hobs_gg_phi3_phi1phi1_bbtautau_CMS8::Hobs_gg_phi3_phi1phi1_bbtautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi3_phi1phi1_bbtautau_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_phi1phi1_bbtautau_CMS8;
}

Hobs_pp_phi3_phi1phi1_bbtautau_CMS8::Hobs_pp_phi3_phi1phi1_bbtautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi1_bbtautau_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi1_bbtautau_CMS8;
}

Hobs_pp_phi3_phi1phi1_bbbb_ATLAS13::Hobs_pp_phi3_phi1phi1_bbbb_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi1_bbbb_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi1_bbbb_ATLAS13;
}

Hobs_pp_phi3_phi1phi1_bbbb_1_CMS13::Hobs_pp_phi3_phi1phi1_bbbb_1_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi1_bbbb_1_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi1_bbbb_1_CMS13;
}

Hobs_pp_phi3_phi1phi1_bbbb_2_CMS13::Hobs_pp_phi3_phi1phi1_bbbb_2_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi1_bbbb_2_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi1_bbbb_2_CMS13;
}

Hobs_pp_phi3_phi1phi1_bbgaga_ATLAS13::Hobs_pp_phi3_phi1phi1_bbgaga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi1_bbgaga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi1_bbgaga_ATLAS13;
}

Hobs_pp_phi3_phi1phi1_bbgaga_CMS13::Hobs_pp_phi3_phi1phi1_bbgaga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi1_bbgaga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi1_bbgaga_CMS13;
}

Hobs_pp_phi3_phi1phi1_bbtautau_ATLAS13::Hobs_pp_phi3_phi1phi1_bbtautau_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi1_bbtautau_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi1_bbtautau_ATLAS13;
}

Hobs_pp_phi3_phi1phi1_bbtautau_1_CMS13::Hobs_pp_phi3_phi1phi1_bbtautau_1_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi1_bbtautau_1_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi1_bbtautau_1_CMS13;
}

Hobs_pp_phi3_phi1phi1_bbtautau_2_CMS13::Hobs_pp_phi3_phi1phi1_bbtautau_2_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi1_bbtautau_2_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi1_bbtautau_2_CMS13;
}

Hobs_pp_phi3_phi1phi1_bbVV_CMS13::Hobs_pp_phi3_phi1phi1_bbVV_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi1_bbVV_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi1_bbVV_CMS13;
}


Hobs_pp_phi3_phi1phi1_bbWW_ATLAS13::Hobs_pp_phi3_phi1phi1_bbWW_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi1_bbWW_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi1_bbWW_ATLAS13;
}

Hobs_gg_phi3_phi1phi1_gagaWW_ATLAS13::Hobs_gg_phi3_phi1phi1_gagaWW_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi3_phi1phi1_gagaWW_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_phi1phi1_gagaWW_ATLAS13;
}

Hobs_gg_phi3_phi1Z_bbZ_ATLAS8::Hobs_gg_phi3_phi1Z_bbZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi3_phi1Z_bbZ_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_phi1Z_bbZ_ATLAS8;
}

Hobs_gg_phi3_phi1Z_bbll_CMS8::Hobs_gg_phi3_phi1Z_bbll_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi3_phi1Z_bbll_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_phi1Z_bbll_CMS8;
}

Hobs_gg_phi3_phi1Z_tautauZ_ATLAS8::Hobs_gg_phi3_phi1Z_tautauZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi3_phi1Z_tautauZ_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_phi1Z_tautauZ_ATLAS8;
}

Hobs_gg_phi3_phi1Z_tautaull_CMS8::Hobs_gg_phi3_phi1Z_tautaull_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi3_phi1Z_tautaull_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_phi1Z_tautaull_CMS8;
}

Hobs_gg_phi3_phi1Z_bbZ_ATLAS13::Hobs_gg_phi3_phi1Z_bbZ_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi3_phi1Z_bbZ_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_phi1Z_bbZ_ATLAS13;
}

Hobs_gg_phi3_phi1Z_bbZ_1_CMS13::Hobs_gg_phi3_phi1Z_bbZ_1_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi3_phi1Z_bbZ_1_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_phi1Z_bbZ_1_CMS13;
}

Hobs_gg_phi3_phi1Z_bbZ_2_CMS13::Hobs_gg_phi3_phi1Z_bbZ_2_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi3_phi1Z_bbZ_2_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_phi1Z_bbZ_2_CMS13;
}

Hobs_bb_phi3_phi1Z_bbZ_ATLAS13::Hobs_bb_phi3_phi1Z_bbZ_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_bb_phi3_phi1Z_bbZ_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_phi1Z_bbZ_ATLAS13;
}

Hobs_bb_phi3_phi1Z_bbZ_1_CMS13::Hobs_bb_phi3_phi1Z_bbZ_1_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_bb_phi3_phi1Z_bbZ_1_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_phi1Z_bbZ_1_CMS13;
}

Hobs_bb_phi3_phi1Z_bbZ_2_CMS13::Hobs_bb_phi3_phi1Z_bbZ_2_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_bb_phi3_phi1Z_bbZ_2_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_phi1Z_bbZ_2_CMS13;
}

Hobs_pp_phi3_phi2Z_bbll_1_CMS8::Hobs_pp_phi3_phi2Z_bbll_1_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi2Z_bbll_1_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi2Z_bbll_1_CMS8;
}

Hobs_pp_phi2_phi3Z_bbll_1_CMS8::Hobs_pp_phi2_phi3Z_bbll_1_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi3Z_bbll_1_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi3Z_bbll_1_CMS8;
}

Hobs_pp_phi3_phi2Z_bbll_2_CMS8::Hobs_pp_phi3_phi2Z_bbll_2_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi2Z_bbll_2_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi2Z_bbll_2_CMS8;
}
Hobs_pp_phi2_phi3Z_bbll_2_CMS8::Hobs_pp_phi2_phi3Z_bbll_2_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi3Z_bbll_2_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi3Z_bbll_2_CMS8;
}


Hobs_pp_phi3_phi2Z_tautaull_1_CMS8::Hobs_pp_phi3_phi2Z_tautaull_1_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi2Z_tautaull_1_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi2Z_tautaull_1_CMS8;
}

Hobs_pp_phi2_phi3Z_tautaull_1_CMS8::Hobs_pp_phi2_phi3Z_tautaull_1_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi3Z_tautaull_1_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi3Z_tautaull_1_CMS8;
}
Hobs_pp_phi3_phi2Z_tautaull_2_CMS8::Hobs_pp_phi3_phi2Z_tautaull_2_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi2Z_tautaull_2_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi2Z_tautaull_2_CMS8;
}
Hobs_pp_phi2_phi3Z_tautaull_2_CMS8::Hobs_pp_phi2_phi3Z_tautaull_2_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi3Z_tautaull_2_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi3Z_tautaull_2_CMS8;
}


Hobs_gg_phi3_phi2Z_bbZ_ATLAS13::Hobs_gg_phi3_phi2Z_bbZ_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi3_phi2Z_bbZ_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_phi2Z_bbZ_ATLAS13;
}
Hobs_gg_phi2_phi3Z_bbZ_ATLAS13::Hobs_gg_phi2_phi3Z_bbZ_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi2_phi3Z_bbZ_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi2_phi3Z_bbZ_ATLAS13;
}

Hobs_bb_phi3_phi2Z_bbZ_ATLAS13::Hobs_bb_phi3_phi2Z_bbZ_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_bb_phi3_phi2Z_bbZ_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_bb_phi3_phi2Z_bbZ_ATLAS13;
}
Hobs_bb_phi2_phi3Z_bbZ_ATLAS13::Hobs_bb_phi2_phi3Z_bbZ_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_bb_phi2_phi3Z_bbZ_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_bb_phi3_phi2Z_bbZ_ATLAS13;
}


Hobs_pp_Hpm_taunu_ATLAS8_GTHDM::Hobs_pp_Hpm_taunu_ATLAS8_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_Hpm_taunu_ATLAS8_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_Hpm_taunu_ATLAS8;
}

Hobs_pp_Hp_taunu_CMS8_GTHDM::Hobs_pp_Hp_taunu_CMS8_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_Hp_taunu_CMS8_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_Hp_taunu_CMS8;
}

Hobs_pp_Hpm_taunu_ATLAS13_GTHDM::Hobs_pp_Hpm_taunu_ATLAS13_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_Hpm_taunu_ATLAS13_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_Hpm_taunu_ATLAS13;
}

Hobs_pp_Hpm_taunu_CMS13_GTHDM::Hobs_pp_Hpm_taunu_CMS13_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_Hpm_taunu_CMS13_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_Hpm_taunu_CMS13;
}

Hobs_pp_Hpm_tb_ATLAS8_GTHDM::Hobs_pp_Hpm_tb_ATLAS8_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_Hpm_tb_ATLAS8_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_Hpm_tb_ATLAS8;
}


Hobs_pp_Hp_tb_CMS8_GTHDM::Hobs_pp_Hp_tb_CMS8_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_Hp_tb_CMS8_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_Hp_tb_CMS8;
}

Hobs_pp_Hpm_tb_ATLAS13::Hobs_pp_Hpm_tb_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_Hpm_tb_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_Hpm_tb_ATLAS13;
}


log10_tt_phi2_tt_TH13::log10_tt_phi2_tt_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_tt_phi2_tt_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->tt_phi2_tt_TH13);
}

log10_tt_phi3_tt_TH13::log10_tt_phi3_tt_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_tt_phi3_tt_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->tt_phi3_tt_TH13);
}

log10_bb_phi2_tt_TH13::log10_bb_phi2_tt_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_bb_phi2_tt_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->bb_phi2_tt_TH13);
}

log10_bb_phi3_tt_TH13::log10_bb_phi3_tt_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_bb_phi3_tt_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->bb_phi3_tt_TH13);
}


log10_bb_phi2_bb_TH8::log10_bb_phi2_bb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_bb_phi2_bb_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->bb_phi2_bb_TH8);
}

log10_bb_phi3_bb_TH8::log10_bb_phi3_bb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_bb_phi3_bb_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->bb_phi3_bb_TH8);
}

log10_gg_phi2_bb_TH8::log10_gg_phi2_bb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi2_bb_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi2_bb_TH8);
}


log10_gg_phi3_bb_TH8::log10_gg_phi3_bb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi3_bb_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi3_bb_TH8);
}

log10_pp_phi2_bb_TH13::log10_pp_phi2_bb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_bb_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_bb_TH13);
}

log10_pp_phi3_bb_TH13::log10_pp_phi3_bb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_bb_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_bb_TH13);
}


log10_bb_phi2_bb_TH13::log10_bb_phi2_bb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_bb_phi2_bb_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->bb_phi2_bb_TH13);
}

log10_bb_phi3_bb_TH13::log10_bb_phi3_bb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_bb_phi3_bb_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->bb_phi3_bb_TH13);
}


log10_gg_phi2_tautau_TH8::log10_gg_phi2_tautau_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi2_tautau_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi2_tautau_TH8);
}

log10_gg_phi3_tautau_TH8::log10_gg_phi3_tautau_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi3_tautau_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi3_tautau_TH8);
}

log10_bb_phi2_tautau_TH8::log10_bb_phi2_tautau_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_bb_phi2_tautau_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->bb_phi2_tautau_TH8);
}

log10_bb_phi3_tautau_TH8::log10_bb_phi3_tautau_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_bb_phi3_tautau_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->bb_phi3_tautau_TH8);
}

log10_gg_phi2_tautau_TH13::log10_gg_phi2_tautau_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi2_tautau_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi2_tautau_TH13);
}

log10_gg_phi3_tautau_TH13::log10_gg_phi3_tautau_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi3_tautau_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi3_tautau_TH13);
}

log10_bb_phi2_tautau_TH13::log10_bb_phi2_tautau_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_bb_phi2_tautau_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->bb_phi2_tautau_TH13);
}

log10_bb_phi3_tautau_TH13::log10_bb_phi3_tautau_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_bb_phi3_tautau_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->bb_phi3_tautau_TH13);
}


log10_gg_phi2_gaga_TH8::log10_gg_phi2_gaga_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi2_gaga_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi2_gaga_TH8);
}

log10_gg_phi3_gaga_TH8::log10_gg_phi3_gaga_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi3_gaga_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi3_gaga_TH8);
}

log10_pp_phi2_gaga_TH13::log10_pp_phi2_gaga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_gaga_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_gaga_TH13);
}

log10_pp_phi3_gaga_TH13::log10_pp_phi3_gaga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_gaga_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_gaga_TH13);
}

log10_gg_phi2_gaga_TH13::log10_gg_phi2_gaga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi2_gaga_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi2_gaga_TH13);
}

log10_gg_phi3_gaga_TH13::log10_gg_phi3_gaga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi3_gaga_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi3_gaga_TH13);
}


log10_pp_phi2_Zga_llga_TH8::log10_pp_phi2_Zga_llga_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_Zga_llga_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_Zga_llga_TH8);
}


log10_pp_phi3_Zga_llga_TH8::log10_pp_phi3_Zga_llga_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_Zga_llga_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_Zga_llga_TH8);
}

log10_gg_phi2_Zga_TH13::log10_gg_phi2_Zga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi2_Zga_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi2_Zga_TH13);
}

log10_gg_phi3_Zga_TH13::log10_gg_phi3_Zga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi3_Zga_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi3_Zga_TH13);
}


log10_gg_phi2_ZZ_TH8::log10_gg_phi2_ZZ_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi2_ZZ_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi2_ZZ_TH8);
}

log10_gg_phi3_ZZ_TH8::log10_gg_phi3_ZZ_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi3_ZZ_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi3_ZZ_TH8);
}

log10_VV_phi2_ZZ_TH8::log10_VV_phi2_ZZ_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_VV_phi2_ZZ_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->VV_phi2_ZZ_TH8);
}


log10_VV_phi3_ZZ_TH8::log10_VV_phi3_ZZ_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_VV_phi3_ZZ_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->VV_phi3_ZZ_TH8);
}

log10_gg_phi2_ZZ_TH13::log10_gg_phi2_ZZ_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi2_ZZ_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi2_ZZ_TH13);
}

log10_gg_phi3_ZZ_TH13::log10_gg_phi3_ZZ_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi3_ZZ_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi3_ZZ_TH13);
}

log10_VV_phi2_ZZ_TH13::log10_VV_phi2_ZZ_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_VV_phi2_ZZ_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->VV_phi2_ZZ_TH13);
}

log10_VV_phi3_ZZ_TH13::log10_VV_phi3_ZZ_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_VV_phi3_ZZ_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->VV_phi3_ZZ_TH13);
}

log10_pp_phi2_ZZ_TH13::log10_pp_phi2_ZZ_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_ZZ_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_ZZ_TH13);
}


log10_pp_phi3_ZZ_TH13::log10_pp_phi3_ZZ_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_ZZ_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_ZZ_TH13);
}


log10_gg_phi2_WW_TH8::log10_gg_phi2_WW_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi2_WW_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi2_WW_TH8);
}


log10_gg_phi3_WW_TH8::log10_gg_phi3_WW_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi3_WW_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi3_WW_TH8);
}

log10_VV_phi2_WW_TH8::log10_VV_phi2_WW_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_VV_phi2_WW_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->VV_phi2_WW_TH8);
}

log10_VV_phi3_WW_TH8::log10_VV_phi3_WW_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_VV_phi3_WW_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->VV_phi3_WW_TH8);
}

log10_gg_phi2_WW_TH13::log10_gg_phi2_WW_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi2_WW_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi2_WW_TH13);
}

log10_gg_phi3_WW_TH13::log10_gg_phi3_WW_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi3_WW_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi3_WW_TH13);
}

log10_VV_phi2_WW_TH13::log10_VV_phi2_WW_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_VV_phi2_WW_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->VV_phi2_WW_TH13);
}

log10_VV_phi3_WW_TH13::log10_VV_phi3_WW_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_VV_phi3_WW_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->VV_phi3_WW_TH13);
}


log10_ggVV_phi2_WW_lnulnu_TH13::log10_ggVV_phi2_WW_lnulnu_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggVV_phi2_WW_lnulnu_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggVV_phi2_WW_lnulnu_TH13);
}

log10_ggVV_phi3_WW_lnulnu_TH13::log10_ggVV_phi3_WW_lnulnu_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggVV_phi3_WW_lnulnu_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggVV_phi3_WW_lnulnu_TH13);
}

log10_pp_phi2_WW_TH13::log10_pp_phi2_WW_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_WW_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_WW_TH13);
}


log10_pp_phi3_WW_TH13::log10_pp_phi3_WW_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_WW_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_WW_TH13);
}


log10_mu_pp_phi2_VV_TH8::log10_mu_pp_phi2_VV_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_mu_pp_phi2_VV_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->mu_pp_phi2_VV_TH8);
}

log10_mu_pp_phi3_VV_TH8::log10_mu_pp_phi3_VV_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_mu_pp_phi3_VV_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->mu_pp_phi3_VV_TH8);
}

log10_pp_phi2_VV_TH13::log10_pp_phi2_VV_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_VV_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_VV_TH13);
}

log10_pp_phi3_VV_TH13::log10_pp_phi3_VV_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_VV_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_VV_TH13);
}


log10_gg_phi2_phi1phi1_TH8::log10_gg_phi2_phi1phi1_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi2_phi1phi1_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi2_phi1phi1_TH8);
}

log10_gg_phi3_phi1phi1_TH8::log10_gg_phi3_phi1phi1_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi3_phi1phi1_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi3_phi1phi1_TH8);
}

log10_pp_phi2_phi1phi1_bbbb_TH8::log10_pp_phi2_phi1phi1_bbbb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_phi1phi1_bbbb_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_phi1phi1_bbbb_TH8);
}

log10_pp_phi3_phi1phi1_bbbb_TH8::log10_pp_phi3_phi1phi1_bbbb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1phi1_bbbb_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1phi1_bbbb_TH8);
}

log10_pp_phi2_phi1phi1_bbgaga_TH8::log10_pp_phi2_phi1phi1_bbgaga_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_phi1phi1_bbgaga_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_phi1phi1_bbgaga_TH8);
}

log10_pp_phi3_phi1phi1_bbgaga_TH8::log10_pp_phi3_phi1phi1_bbgaga_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1phi1_bbgaga_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1phi1_bbgaga_TH8);
}

log10_gg_phi2_phi1phi1_bbtautau_TH8::log10_gg_phi2_phi1phi1_bbtautau_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi2_phi1phi1_bbtautau_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi2_phi1phi1_bbtautau_TH8);
}

log10_gg_phi3_phi1phi1_bbtautau_TH8::log10_gg_phi3_phi1phi1_bbtautau_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi3_phi1phi1_bbtautau_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi3_phi1phi1_bbtautau_TH8);
}


log10_pp_phi2_phi1phi1_TH8::log10_pp_phi2_phi1phi1_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_phi1phi1_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_phi1phi1_TH8);
}

log10_pp_phi3_phi1phi1_TH8::log10_pp_phi3_phi1phi1_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1phi1_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1phi1_TH8);
}

log10_pp_phi2_phi1phi1_bbbb_TH13::log10_pp_phi2_phi1phi1_bbbb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_phi1phi1_bbbb_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_phi1phi1_bbbb_TH13);
}


log10_pp_phi3_phi1phi1_bbbb_TH13::log10_pp_phi3_phi1phi1_bbbb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1phi1_bbbb_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1phi1_bbbb_TH13);
}

log10_pp_phi2_phi1phi1_TH13::log10_pp_phi2_phi1phi1_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_phi1phi1_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_phi1phi1_TH13);
}


log10_pp_phi3_phi1phi1_TH13::log10_pp_phi3_phi1phi1_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1phi1_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1phi1_TH13);
}


log10_pp_phi2_phi1phi1_bbgaga_TH13::log10_pp_phi2_phi1phi1_bbgaga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_phi1phi1_bbgaga_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_phi1phi1_bbgaga_TH13);
}

log10_pp_phi3_phi1phi1_bbgaga_TH13::log10_pp_phi3_phi1phi1_bbgaga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1phi1_bbgaga_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1phi1_bbgaga_TH13);
}

log10_pp_phi2_phi1phi1_bbtautau_TH13::log10_pp_phi2_phi1phi1_bbtautau_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_phi1phi1_bbtautau_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_phi1phi1_bbtautau_TH13);
}

log10_pp_phi3_phi1phi1_bbtautau_TH13::log10_pp_phi3_phi1phi1_bbtautau_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1phi1_bbtautau_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1phi1_bbtautau_TH13);
}

log10_pp_phi2_phi1phi1_bbVV_TH13::log10_pp_phi2_phi1phi1_bbVV_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_phi1phi1_bbVV_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_phi1phi1_bbVV_TH13);
}

log10_pp_phi3_phi1phi1_bbVV_TH13::log10_pp_phi3_phi1phi1_bbVV_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1phi1_bbVV_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1phi1_bbVV_TH13);
}

log10_gg_phi2_phi1phi1_gagaWW_TH13::log10_gg_phi2_phi1phi1_gagaWW_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi2_phi1phi1_gagaWW_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi2_phi1phi1_gagaWW_TH13);
}

log10_gg_phi3_phi1phi1_gagaWW_TH13::log10_gg_phi3_phi1phi1_gagaWW_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi3_phi1phi1_gagaWW_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi3_phi1phi1_gagaWW_TH13);
}



log10_gg_phi2_phi1Z_bbZ_TH8::log10_gg_phi2_phi1Z_bbZ_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi2_phi1Z_bbZ_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi2_phi1Z_bbZ_TH8);
}

log10_gg_phi3_phi1Z_bbZ_TH8::log10_gg_phi3_phi1Z_bbZ_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi3_phi1Z_bbZ_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi3_phi1Z_bbZ_TH8);
}

log10_gg_phi2_phi1Z_bbll_TH8::log10_gg_phi2_phi1Z_bbll_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi2_phi1Z_bbll_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi2_phi1Z_bbll_TH8);
}

log10_gg_phi3_phi1Z_bbll_TH8::log10_gg_phi3_phi1Z_bbll_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi3_phi1Z_bbll_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi3_phi1Z_bbll_TH8);
}

log10_gg_phi2_phi1Z_tautauZ_TH8::log10_gg_phi2_phi1Z_tautauZ_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi2_phi1Z_tautauZ_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi2_phi1Z_tautauZ_TH8);
}


log10_gg_phi3_phi1Z_tautauZ_TH8::log10_gg_phi3_phi1Z_tautauZ_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi3_phi1Z_tautauZ_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi3_phi1Z_tautauZ_TH8);
}

log10_gg_phi2_phi1Z_tautaull_TH8::log10_gg_phi2_phi1Z_tautaull_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi2_phi1Z_tautaull_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi2_phi1Z_tautaull_TH8);
}

log10_gg_phi3_phi1Z_tautaull_TH8::log10_gg_phi3_phi1Z_tautaull_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi3_phi1Z_tautaull_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi3_phi1Z_tautaull_TH8);
}

log10_gg_phi2_phi1Z_bbZ_TH13::log10_gg_phi2_phi1Z_bbZ_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi2_phi1Z_bbZ_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi2_phi1Z_bbZ_TH13);
}

log10_gg_phi3_phi1Z_bbZ_TH13::log10_gg_phi3_phi1Z_bbZ_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi3_phi1Z_bbZ_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi3_phi1Z_bbZ_TH13);
}

log10_bb_phi2_phi1Z_bbZ_TH13::log10_bb_phi2_phi1Z_bbZ_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_bb_phi2_phi1Z_bbZ_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->bb_phi2_phi1Z_bbZ_TH13);
}

log10_bb_phi3_phi1Z_bbZ_TH13::log10_bb_phi3_phi1Z_bbZ_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_bb_phi3_phi1Z_bbZ_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->bb_phi3_phi1Z_bbZ_TH13);
}


log10_pp_phi3_phi2Z_bbll_TH8::log10_pp_phi3_phi2Z_bbll_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi2Z_bbll_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi2Z_bbll_TH8);
}

log10_pp_phi3_phi2Z_tautaull_TH8::log10_pp_phi3_phi2Z_tautaull_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi2Z_tautaull_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi2Z_tautaull_TH8);
}

log10_gg_phi3_phi2Z_bbZ_TH13::log10_gg_phi3_phi2Z_bbZ_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_gg_phi3_phi2Z_bbZ_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->gg_phi3_phi2Z_bbZ_TH13);
}


log10_bb_phi3_phi2Z_bbZ_TH13::log10_bb_phi3_phi2Z_bbZ_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_bb_phi3_phi2Z_bbZ_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->bb_phi3_phi2Z_bbZ_TH13);
}


log10_pp_Hpm_taunu_TH8_GTHDM::log10_pp_Hpm_taunu_TH8_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_Hpm_taunu_TH8_GTHDM::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_Hpm_taunu_TH8);
}

log10_pp_Hp_taunu_TH8_GTHDM::log10_pp_Hp_taunu_TH8_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_Hp_taunu_TH8_GTHDM::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_Hp_taunu_TH8);
}

log10_pp_Hpm_taunu_TH13_GTHDM::log10_pp_Hpm_taunu_TH13_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_Hpm_taunu_TH13_GTHDM::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_Hpm_taunu_TH13);
}


log10_pp_Hpm_tb_TH8_GTHDM::log10_pp_Hpm_tb_TH8_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_Hpm_tb_TH8_GTHDM::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_Hpm_tb_TH8);
}


log10_pp_Hp_tb_TH8_GTHDM::log10_pp_Hp_tb_TH8_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_Hp_tb_TH8_GTHDM::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_Hp_tb_TH8);
}

log10_pp_Hpm_tb_TH13_GTHDM::log10_pp_Hpm_tb_TH13_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_Hpm_tb_TH13_GTHDM::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_Hpm_tb_TH13);
}
/*NOT CLEANED YET*/

Gamma_phi3_GTHDM::Gamma_phi3_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Gamma_phi3_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Gammaphi3tot;
}



//rHH_gaga_THDM::rHH_gaga_THDM(const StandardModel& SM_i)
//: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
//{}
//
//double rHH_gaga_THDM::computeThValue()
//{
//    return myGTHDM.getMyGTHDMCache()->rHH_gaga;
//}



rphi3_ggE_GTHDM::rphi3_ggE_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double rphi3_ggE_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->rphi3_ggE;
}

rphi3_ggO_GTHDM::rphi3_ggO_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double rphi3_ggO_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->rphi3_ggO;
}



// phi3 -> phi1phi1 and phi3 -> phi2phi2

BR_phi3_phi1phi1_GTHDM::BR_phi3_phi1phi1_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double BR_phi3_phi1phi1_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Br_phi3tophi1phi1;
}


BR_phi3_phi2phi2_GTHDM::BR_phi3_phi2phi2_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double BR_phi3_phi2phi2_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Br_phi3tophi2phi2;
}



//phi3 -> H+H-


BR_phi3_HpHm_GTHDM::BR_phi3_HpHm_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double BR_phi3_HpHm_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Br_phi3toHpHm;
}



BR_phi3_phi1Z_GTHDM::BR_phi3_phi1Z_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double BR_phi3_phi1Z_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Br_phi3tophi1Z;
}

BR_phi3_phi2Z_GTHDM::BR_phi3_phi2Z_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double BR_phi3_phi2Z_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Br_phi3tophi2Z;
}



BR_phi3_HpW_GTHDM::BR_phi3_HpW_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double BR_phi3_HpW_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Br_phi3toHpW;
}



Gamma_phi2_GTHDM::Gamma_phi2_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Gamma_phi2_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Gammaphi2tot;
}



//rHH_gaga_THDM::rHH_gaga_THDM(const StandardModel& SM_i)
//: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
//{}
//
//double rHH_gaga_THDM::computeThValue()
//{
//    return myGTHDM.getMyGTHDMCache()->rHH_gaga;
//}



rphi2_ggE_GTHDM::rphi2_ggE_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double rphi2_ggE_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->rphi2_ggE;
}

rphi2_ggO_GTHDM::rphi2_ggO_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double rphi2_ggO_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->rphi2_ggO;
}

// phi2 -> phi1phi1 and phi2 -> phi2phi2

BR_phi2_phi1phi1_GTHDM::BR_phi2_phi1phi1_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double BR_phi2_phi1phi1_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Br_phi2tophi1phi1;
}

//phi2 -> H+H-


BR_phi2_HpHm_GTHDM::BR_phi2_HpHm_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double BR_phi2_HpHm_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Br_phi2toHpHm;
}


BR_phi2_phi1Z_GTHDM::BR_phi2_phi1Z_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double BR_phi2_phi1Z_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Br_phi2tophi1Z;
}

BR_phi2_HpW_GTHDM::BR_phi2_HpW_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double BR_phi2_HpW_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Br_phi2toHpW;
}
