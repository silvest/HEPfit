/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "CPoddHiggs.h"
#include "StandardModel.h"



Hobs_ggF_A_tautau_ATLAS8::Hobs_ggF_A_tautau_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_A_tautau_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_A_tautau_ATLAS8;
}

Robs_ggF_A_tautau_ATLAS8::Robs_ggF_A_tautau_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_A_tautau_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_A_tautau_ATLAS8;
}



Hobs_ggF_A_tautau_CMS8::Hobs_ggF_A_tautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_A_tautau_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_A_tautau_CMS8;
}

Robs_ggF_A_tautau_CMS8::Robs_ggF_A_tautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_A_tautau_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_A_tautau_CMS8;
}



Hobs_bbF_A_tautau_ATLAS8::Hobs_bbF_A_tautau_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_bbF_A_tautau_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_bbF_A_tautau_ATLAS8;
}

Robs_bbF_A_tautau_ATLAS8::Robs_bbF_A_tautau_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_bbF_A_tautau_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_bbF_A_tautau_ATLAS8;
}



Hobs_bbF_A_tautau_CMS8::Hobs_bbF_A_tautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_bbF_A_tautau_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_bbF_A_tautau_CMS8;
}

Robs_bbF_A_tautau_CMS8::Robs_bbF_A_tautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_bbF_A_tautau_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_bbF_A_tautau_CMS8;
}



Hobs_pp_A_gaga_ATLAS8::Hobs_pp_A_gaga_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_A_gaga_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_A_gaga_ATLAS8;
}

Robs_pp_A_gaga_ATLAS8::Robs_pp_A_gaga_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_A_gaga_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_A_gaga_ATLAS8;
}



Hobs_ggF_A_gaga_CMS8::Hobs_ggF_A_gaga_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_A_gaga_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_A_gaga_CMS8;
}

Robs_ggF_A_gaga_CMS8::Robs_ggF_A_gaga_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_A_gaga_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_A_gaga_CMS8;
}



Hobs_pp_A_Zga_llga_CMS8::Hobs_pp_A_Zga_llga_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_A_Zga_llga_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_A_Zga_llga_CMS8;
}

Robs_pp_A_Zga_llga_CMS8::Robs_pp_A_Zga_llga_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_A_Zga_llga_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_A_Zga_llga_CMS8;
}



Hobs_ggF_A_hZ_bbll_CMS8::Hobs_ggF_A_hZ_bbll_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_A_hZ_bbll_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_A_hZ_bbll_CMS8;
}

Robs_ggF_A_hZ_bbll_CMS8::Robs_ggF_A_hZ_bbll_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_A_hZ_bbll_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_A_hZ_bbll_CMS8;
}



Hobs_ggF_A_hZ_bbZ_ATLAS8::Hobs_ggF_A_hZ_bbZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_A_hZ_bbZ_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_A_hZ_bbZ_ATLAS8;
}

Robs_ggF_A_hZ_bbZ_ATLAS8::Robs_ggF_A_hZ_bbZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_A_hZ_bbZ_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_A_hZ_bbZ_ATLAS8;
}



Hobs_ggF_A_hZ_tautaull_CMS8::Hobs_ggF_A_hZ_tautaull_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_A_hZ_tautaull_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_A_hZ_tautaull_CMS8;
}

Robs_ggF_A_hZ_tautaull_CMS8::Robs_ggF_A_hZ_tautaull_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_A_hZ_tautaull_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_A_hZ_tautaull_CMS8;
}



Hobs_ggF_A_hZ_tautauZ_ATLAS8::Hobs_ggF_A_hZ_tautauZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_A_hZ_tautauZ_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_A_hZ_tautauZ_ATLAS8;
}

Robs_ggF_A_hZ_tautauZ_ATLAS8::Robs_ggF_A_hZ_tautauZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_A_hZ_tautauZ_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_A_hZ_tautauZ_ATLAS8;
}



Hobs_ggF_A_tt_ATLAS8::Hobs_ggF_A_tt_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_A_tt_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_A_tt_ATLAS8;
}

Robs_ggF_A_tt_ATLAS8::Robs_ggF_A_tt_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_A_tt_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_A_tt_ATLAS8;
}



Hobs_bbF_A_bb_CMS8::Hobs_bbF_A_bb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_bbF_A_bb_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_bbF_A_bb_CMS8;
}

Robs_bbF_A_bb_CMS8::Robs_bbF_A_bb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_bbF_A_bb_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_bbF_A_bb_CMS8;
}



log10_ggF_A_tautau_TH::log10_ggF_A_tautau_TH(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_ggF_A_tautau_TH::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->ggF_A_tautau_TH);
}



log10_bbF_A_tautau_TH::log10_bbF_A_tautau_TH(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_bbF_A_tautau_TH::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->bbF_A_tautau_TH);
}



log10_pp_A_gaga_TH::log10_pp_A_gaga_TH(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_pp_A_gaga_TH::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->pp_A_gaga_TH);
}



log10_ggF_A_gaga_TH::log10_ggF_A_gaga_TH(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_ggF_A_gaga_TH::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->ggF_A_gaga_TH);
}



log10_pp_A_Zga_llga_TH::log10_pp_A_Zga_llga_TH(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_pp_A_Zga_llga_TH::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->pp_A_Zga_llga_TH);
}



log10_ggF_A_hZ_bbll_TH::log10_ggF_A_hZ_bbll_TH(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_ggF_A_hZ_bbll_TH::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->ggF_A_hZ_bbll_TH);
}



log10_ggF_A_hZ_bbZ_TH::log10_ggF_A_hZ_bbZ_TH(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_ggF_A_hZ_bbZ_TH::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->ggF_A_hZ_bbZ_TH);
}



log10_ggF_A_hZ_tautaull_TH::log10_ggF_A_hZ_tautaull_TH(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_ggF_A_hZ_tautaull_TH::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->ggF_A_hZ_tautaull_TH);
}



log10_ggF_A_hZ_tautauZ_TH::log10_ggF_A_hZ_tautauZ_TH(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_ggF_A_hZ_tautauZ_TH::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->ggF_A_hZ_tautauZ_TH);
}



log10_ggF_A_tt_TH::log10_ggF_A_tt_TH(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_ggF_A_tt_TH::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->ggF_A_tt_TH);
}



log10_bbF_A_bb_TH::log10_bbF_A_bb_TH(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_bbF_A_bb_TH::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->bbF_A_bb_TH);
}



Gamma_A_THDM::Gamma_A_THDM(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Gamma_A_THDM::computeThValue()
{
    return myTHDM.getMyTHDMCache()->GammaAtot;
}



//rA_gaga_THDM::rA_gaga_THDM(const StandardModel& SM_i)
//: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
//{}
//
//double rA_gaga_THDM::computeThValue()
//{
//    return myTHDM.getMyTHDMCache()->rA_gaga;
//}



rA_gg_THDM::rA_gg_THDM(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double rA_gg_THDM::computeThValue()
{
    return myTHDM.getMyTHDMCache()->rA_gg;
}



BR_A_HZ_THDM::BR_A_HZ_THDM(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double BR_A_HZ_THDM::computeThValue()
{
    return myTHDM.getMyTHDMCache()->Br_AtoHZ;
}



BR_A_hZ_THDM::BR_A_hZ_THDM(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double BR_A_hZ_THDM::computeThValue()
{
    return myTHDM.getMyTHDMCache()->Br_AtohZ;
}



BR_A_HpW_THDM::BR_A_HpW_THDM(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double BR_A_HpW_THDM::computeThValue()
{
    return myTHDM.getMyTHDMCache()->Br_AtoHpW;
}
