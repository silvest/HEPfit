/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "CPoddHiggs.h"
//#include "lightHiggs.h"
#include "StandardModel.h"

//CPoddHiggs::CPoddHiggs(const StandardModel& SM_i):
//
//        ThObservable(SM_i), 
//        myTHDM(static_cast<const THDM&> (SM_i))
//{
//    mycache = new THDMcache(SM_i);
////    mylightHiggs = new lightHiggs(SM_i);
//}
//
//CPoddHiggs::~CPoddHiggs()
//{}
//
//void CPoddHiggs::updateCPoddHiggsParameters() {
//    bma=myTHDM->getbma();
//    tanb=myTHDM->gettanb();
//    mHh2=myTHDM->getmHh2();
//    mA2=myTHDM->getmA2();
//    mHp2=myTHDM->getmHp2();
//    MW=mycache->MWTHDM(myTHDM->Mw_tree());
//    cW2=mycache->cW2THDM(myTHDM->c02());
//    mHl=myTHDM->getMHl();
//    vev=myTHDM->v();
//    Ale=myTHDM->getAle();
//    Als=myTHDM->getAlsMz();
//    Mt=myTHDM->getQuarks(QCD::TOP).getMass();
//    Mb=myTHDM->getQuarks(QCD::BOTTOM).getMass();   
//    Mtau=myTHDM->getLeptons(StandardModel::TAU).getMass();
//    Mc=myTHDM->getQuarks(QCD::CHARM).getMass();
//    Ms=myTHDM->getQuarks(QCD::STRANGE).getMass();
//    Mmu=myTHDM->getLeptons(StandardModel::MU).getMass();
//    MZ=myTHDM->getMz();
//}
//
//void CPoddHiggs::updateCPoddHiggsQuantities(
//    const double bma, const double tanb, const double mHh2, const double mA2,
//    const double mHp2, const double MW, const double cW2,
//    const double mHl, const double vev, const double Ale, const double Als, const double Mt,
//    const double Mb, const double Mtau, const double Mc,
//    const double Ms, const double Mmu, const double MZ)
//{
//    int NumPar = 18;
//    double params[] = {bma,tanb,mHh2,mA2,mHp2,MW,cW2,mHl,vev,Ale,Als,Mt,Mb,Mtau,Mc,Ms,Mmu,MZ};
//
//    bool bCache = true;
//    for(int j=0; j<NumPar; j++) {
//        bCache &= (params[j] == CPoddHiggscache[j]);
//    }
//    if (!bCache)
//    {
//        computeAquantities();
//        for(int j=0; j<NumPar; j++) {
//            CPoddHiggscache[j] = params[j];
//        }
////        std::cout<<"new Acalculation"<<std::endl;
////        std::cout<<bma<<","<<tanb<<","<<mHh2<<","<<mA2<<","<<mHp2<<","<<MW<<","<<cW2<<","<<mHl<<","<<vev<<","<<Ale<<","<<Als<<","<<Mt<<","<<Mb<<","<<Mtau<<std::endl;
//    }
//    else
//    {
////        std::cout<<"Acache"<<std::endl;
//    }
//}
//
//double CPoddHiggs::computeThValue()
//{
//    return 0;
//}



/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/



Hobs_ggF_A_tautau_ATLAS::Hobs_ggF_A_tautau_ATLAS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_A_tautau_ATLAS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_A_tautau_ATLAS;
}



Hobs_ggF_A_tautau_CMS::Hobs_ggF_A_tautau_CMS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_A_tautau_CMS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_A_tautau_CMS;
}



Hobs_bbF_A_tautau_ATLAS::Hobs_bbF_A_tautau_ATLAS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_bbF_A_tautau_ATLAS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_bbF_A_tautau_ATLAS;
}



Hobs_bbF_A_tautau_CMS::Hobs_bbF_A_tautau_CMS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_bbF_A_tautau_CMS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_bbF_A_tautau_CMS;
}



Hobs_ggF_A_gaga_ATLAS::Hobs_ggF_A_gaga_ATLAS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_A_gaga_ATLAS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_A_gaga_ATLAS;
}



Hobs_ggF_A_gaga_CMS::Hobs_ggF_A_gaga_CMS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_A_gaga_CMS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_A_gaga_CMS;
}



Hobs_ggF_A_hZ_bbll_CMS::Hobs_ggF_A_hZ_bbll_CMS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_A_hZ_bbll_CMS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_A_hZ_bbll_CMS;
}



Hobs_ggF_A_hZ_bbZ_ATLAS::Hobs_ggF_A_hZ_bbZ_ATLAS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_A_hZ_bbZ_ATLAS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_A_hZ_bbZ_ATLAS;
}



Hobs_ggF_A_hZ_tautaull_CMS::Hobs_ggF_A_hZ_tautaull_CMS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_A_hZ_tautaull_CMS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_A_hZ_tautaull_CMS;
}



Hobs_ggF_A_hZ_tautauZ_ATLAS::Hobs_ggF_A_hZ_tautauZ_ATLAS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_ggF_A_hZ_tautauZ_ATLAS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_ggF_A_hZ_tautauZ_ATLAS;
}



Hobs_pp_A_tt_ATLAS::Hobs_pp_A_tt_ATLAS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_A_tt_ATLAS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_A_tt_ATLAS;
}



Hobs_bbF_A_bb_CMS::Hobs_bbF_A_bb_CMS(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_bbF_A_bb_CMS::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_bbF_A_bb_CMS;
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



log10_ggF_A_gaga_TH::log10_ggF_A_gaga_TH(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_ggF_A_gaga_TH::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->ggF_A_gaga_TH);
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



log10_pp_A_tt_TH::log10_pp_A_tt_TH(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_pp_A_tt_TH::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->pp_A_tt_TH);
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
