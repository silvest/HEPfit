/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "CPoddHiggs.h"
#include "StandardModel.h"



Hobs_ggF_A_tautau_ATLAS::Hobs_ggF_A_tautau_ATLAS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_A_tautau_ATLAS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_A_tautau_ATLAS;
}

Robs_ggF_A_tautau_ATLAS::Robs_ggF_A_tautau_ATLAS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_A_tautau_ATLAS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_A_tautau_ATLAS;
}



Hobs_ggF_A_tautau_CMS::Hobs_ggF_A_tautau_CMS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_A_tautau_CMS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_A_tautau_CMS;
}

Robs_ggF_A_tautau_CMS::Robs_ggF_A_tautau_CMS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_A_tautau_CMS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_A_tautau_CMS;
}



Hobs_bbF_A_tautau_ATLAS::Hobs_bbF_A_tautau_ATLAS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_bbF_A_tautau_ATLAS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_bbF_A_tautau_ATLAS;
}

Robs_bbF_A_tautau_ATLAS::Robs_bbF_A_tautau_ATLAS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_bbF_A_tautau_ATLAS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_bbF_A_tautau_ATLAS;
}



Hobs_bbF_A_tautau_CMS::Hobs_bbF_A_tautau_CMS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_bbF_A_tautau_CMS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_bbF_A_tautau_CMS;
}

Robs_bbF_A_tautau_CMS::Robs_bbF_A_tautau_CMS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_bbF_A_tautau_CMS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_bbF_A_tautau_CMS;
}



Hobs_pp_A_gaga_ATLAS::Hobs_pp_A_gaga_ATLAS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_A_gaga_ATLAS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_A_gaga_ATLAS;
}

Robs_pp_A_gaga_ATLAS::Robs_pp_A_gaga_ATLAS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_A_gaga_ATLAS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_A_gaga_ATLAS;
}



Hobs_ggF_A_gaga_CMS::Hobs_ggF_A_gaga_CMS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_A_gaga_CMS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_A_gaga_CMS;
}

Robs_ggF_A_gaga_CMS::Robs_ggF_A_gaga_CMS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_A_gaga_CMS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_A_gaga_CMS;
}



Hobs_pp_A_Zga_llga_CMS::Hobs_pp_A_Zga_llga_CMS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_A_Zga_llga_CMS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_A_Zga_llga_CMS;
}

Robs_pp_A_Zga_llga_CMS::Robs_pp_A_Zga_llga_CMS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_A_Zga_llga_CMS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_A_Zga_llga_CMS;
}



Hobs_ggF_A_hZ_bbll_CMS::Hobs_ggF_A_hZ_bbll_CMS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_A_hZ_bbll_CMS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_A_hZ_bbll_CMS;
}

Robs_ggF_A_hZ_bbll_CMS::Robs_ggF_A_hZ_bbll_CMS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_A_hZ_bbll_CMS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_A_hZ_bbll_CMS;
}



Hobs_ggF_A_hZ_bbZ_ATLAS::Hobs_ggF_A_hZ_bbZ_ATLAS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_A_hZ_bbZ_ATLAS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_A_hZ_bbZ_ATLAS;
}

Robs_ggF_A_hZ_bbZ_ATLAS::Robs_ggF_A_hZ_bbZ_ATLAS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_A_hZ_bbZ_ATLAS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_A_hZ_bbZ_ATLAS;
}



Hobs_ggF_A_hZ_tautaull_CMS::Hobs_ggF_A_hZ_tautaull_CMS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_A_hZ_tautaull_CMS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_A_hZ_tautaull_CMS;
}

Robs_ggF_A_hZ_tautaull_CMS::Robs_ggF_A_hZ_tautaull_CMS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_A_hZ_tautaull_CMS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_A_hZ_tautaull_CMS;
}



Hobs_ggF_A_hZ_tautauZ_ATLAS::Hobs_ggF_A_hZ_tautauZ_ATLAS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_A_hZ_tautauZ_ATLAS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_A_hZ_tautauZ_ATLAS;
}

Robs_ggF_A_hZ_tautauZ_ATLAS::Robs_ggF_A_hZ_tautauZ_ATLAS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_A_hZ_tautauZ_ATLAS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_A_hZ_tautauZ_ATLAS;
}



Hobs_ggF_A_tt_ATLAS::Hobs_ggF_A_tt_ATLAS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_A_tt_ATLAS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_A_tt_ATLAS;
}

Robs_ggF_A_tt_ATLAS::Robs_ggF_A_tt_ATLAS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_ggF_A_tt_ATLAS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_ggF_A_tt_ATLAS;
}



Hobs_bbF_A_bb_CMS::Hobs_bbF_A_bb_CMS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_bbF_A_bb_CMS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_bbF_A_bb_CMS;
}

Robs_bbF_A_bb_CMS::Robs_bbF_A_bb_CMS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_bbF_A_bb_CMS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_bbF_A_bb_CMS;
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
