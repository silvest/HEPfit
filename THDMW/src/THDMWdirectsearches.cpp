/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "THDMWdirectsearches.h"
#include "StandardModel.h"



Hobs_pp_Sr_tt_ATLAS13::Hobs_pp_Sr_tt_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double Hobs_pp_Sr_tt_ATLAS13::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->THoEX_pp_Sr_tt;
}


log10_pp_Sr_tt_TH13::log10_pp_Sr_tt_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double log10_pp_Sr_tt_TH13::computeThValue()
{
    return log10(myTHDMW.getMyTHDMWCache()->pp_Sr_tt_TH13);
}

Hobs_pp_Srtt_tttt_ATLAS13::Hobs_pp_Srtt_tttt_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double Hobs_pp_Srtt_tttt_ATLAS13::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->THoEX_pp_Srtt_tttt;
}


log10_pp_Srtt_tttt_TH13::log10_pp_Srtt_tttt_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double log10_pp_Srtt_tttt_TH13::computeThValue()
{
    return log10(myTHDMW.getMyTHDMWCache()->pp_Srtt_tttt_TH13);
}


Hobs_pp_Sr_jj_CMS13::Hobs_pp_Sr_jj_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double Hobs_pp_Sr_jj_CMS13::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->THoEX_pp_Sr_jj;
}


log10_pp_Sr_jj_TH13::log10_pp_Sr_jj_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double log10_pp_Sr_jj_TH13::computeThValue()
{
    return log10(myTHDMW.getMyTHDMWCache()->pp_Sr_jj_TH13);
}







Hobs_pp_SrSr_jjjj_ATLAS13::Hobs_pp_SrSr_jjjj_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double Hobs_pp_SrSr_jjjj_ATLAS13::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->THoEX_pp_SrSr_jjjj;
}


log10_pp_SrSr_jjjj_TH13::log10_pp_SrSr_jjjj_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double log10_pp_SrSr_jjjj_TH13::computeThValue()
{
    return log10(myTHDMW.getMyTHDMWCache()->pp_SrSr_jjjj_TH13);
    //return myTHDMW.getMyTHDMWCache()->pp_SrSr_jjjj_TH13;
}










Hobs_pp_Stb_tbtb_ATLAS13::Hobs_pp_Stb_tbtb_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double Hobs_pp_Stb_tbtb_ATLAS13::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->THoEX_pp_Stb_tbtb;
}


log10_pp_Stb_tbtb_TH13::log10_pp_Stb_tbtb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double log10_pp_Stb_tbtb_TH13::computeThValue()
{
    return log10(myTHDMW.getMyTHDMWCache()->pp_Stb_tbtb_TH13);
}




Hobs_pp_Sitt_tttt_ATLAS13::Hobs_pp_Sitt_tttt_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double Hobs_pp_Sitt_tttt_ATLAS13::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->THoEX_pp_Sitt_tttt;
}


log10_pp_Sitt_tttt_TH13::log10_pp_Sitt_tttt_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double log10_pp_Sitt_tttt_TH13::computeThValue()
{
    return log10(myTHDMW.getMyTHDMWCache()->pp_Sitt_tttt_TH13);
}






Hobs_pp_Srbb_bbbb_CMS13::Hobs_pp_Srbb_bbbb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double Hobs_pp_Srbb_bbbb_CMS13::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->THoEX_pp_Srbb_bbbb;
}


log10_pp_Srbb_bbbb_TH13::log10_pp_Srbb_bbbb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double log10_pp_Srbb_bbbb_TH13::computeThValue()
{
    return log10(myTHDMW.getMyTHDMWCache()->pp_Srbb_bbbb_TH13);
}




Hobs_pp_Sibb_bbbb_CMS13::Hobs_pp_Sibb_bbbb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double Hobs_pp_Sibb_bbbb_CMS13::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->THoEX_pp_Sibb_bbbb;
}


log10_pp_Sibb_bbbb_TH13::log10_pp_Sibb_bbbb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double log10_pp_Sibb_bbbb_TH13::computeThValue()
{
    return log10(myTHDMW.getMyTHDMWCache()->pp_Sibb_bbbb_TH13);
}







Hobs_pp_Sr_bb_CMS13::Hobs_pp_Sr_bb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double Hobs_pp_Sr_bb_CMS13::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->THoEX_pp_Sr_bb;
}


log10_pp_Sr_bb_TH13::log10_pp_Sr_bb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double log10_pp_Sr_bb_TH13::computeThValue()
{
    return log10(myTHDMW.getMyTHDMWCache()->pp_Sr_bb_TH13);
}




Hobs_pp_Sr_bb_CMS8::Hobs_pp_Sr_bb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double Hobs_pp_Sr_bb_CMS8::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->THoEX_pp_Sr_bb_8TeV;
}


log10_pp_Sr_bb_TH8::log10_pp_Sr_bb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double log10_pp_Sr_bb_TH8::computeThValue()
{
    return log10(myTHDMW.getMyTHDMWCache()->pp_Sr_bb_TH8);
}






Hobs_pp_Si_bb_CMS13::Hobs_pp_Si_bb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double Hobs_pp_Si_bb_CMS13::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->THoEX_pp_Si_bb;
}


log10_pp_Si_bb_TH13::log10_pp_Si_bb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double log10_pp_Si_bb_TH13::computeThValue()
{
    return log10(myTHDMW.getMyTHDMWCache()->pp_Si_bb_TH13);
}


//logpp_SrSr_jjjj_TH13::logpp_SrSr_jjjj_TH13(const StandardModel& SM_i)
//: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
//{}

/*
double logpp_SrSr_jjjj_TH13::computeThValue()
{
    //return log10(myTHDMW.getMyTHDMWCache()->pp_SrSr_jjjj_TH13);
    return myTHDMW.getMyTHDMWCache()->logpp_SrSr_jjjj_TH13;
}
*/