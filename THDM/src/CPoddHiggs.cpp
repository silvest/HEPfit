/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "CPoddHiggs.h"
#include "lightHiggs.h"
#include "StandardModel.h"

CPoddHiggsCache::CPoddHiggsCache(const StandardModel& SM_i):

        ThObservable(SM_i),
        myTHDM(static_cast<const THDM*> (&SM_i)),
        mySM (SM_i)

{
    mycache = new THDMcache();
    mylightHiggs = new lightHiggs(SM_i);
}


CPoddHiggsCache::~CPoddHiggsCache()
{}


void CPoddHiggsCache::computeParameters()
{

    THDMfunctions myfunctions(mySM);

    double Mc=myTHDM->getQuarks(QCD::CHARM).getMass();
    double Ms=myTHDM->getQuarks(QCD::STRANGE).getMass();
    double Mt=myTHDM->getQuarks(QCD::TOP).getMass();
    double Mb=myTHDM->getQuarks(QCD::BOTTOM).getMass();
    double Mmu=myTHDM->getLeptons(StandardModel::MU).getMass();
    double Mtau=myTHDM->getLeptons(StandardModel::TAU).getMass();

    double GF=myTHDM->getGF();
    double Ale=myTHDM->getAle();
    double Als=myTHDM->getAlsMz();
    double MZ=myTHDM->getMz();
    double MW=myTHDM->Mw();
    double vev=myTHDM->v();
    double s02=myTHDM->s02();
    double c02=myTHDM->c02();

    double mHh=myTHDM->getmHh();
    double mHl=myTHDM->getMHl();
    double mHp=myTHDM->getmHp();
    double mA=myTHDM->getmA();

    double m12_2=myTHDM->getm12_2();

    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double sin_ba=myTHDM->getsin_ba();
    double cos_ab=cosa*cosb+sina*sinb;

    //These cross sections are necessary for rA_gg
    //SM gg -> h production cross section at 8 TeV
    SigmaggF = mycache->cs_ggFtoHP(mA);
    //SM gg -> h productmycaion cross section at 8 TeV, top loop only
    Sigmaggh_tt = mycache->cs_ggHP_tt(mA);
    //SM gg -> h production cross section at 8 TeV, bottom loop only
    Sigmaggh_bb = mycache->cs_ggHP_bb(mA);

    /* r_ii is the ratio between the squared 2HDM vertex coupling of the CP-odd
     * Higgs to the particle i and the corresponding coupling of the SM Higgs boson.*/
    rA_QuQu=cosb/sinb*cosb/sinb;
    rA_QdQd=0.0;//It depends on the modelType
    rA_ll=0.0;//It depends on the modelType
    rA_gg=0.0;//It depends on the modelType

    //The ga-ga decay width of the CP-odd Higgs boson
    Gamma_Agaga=0.0;//It depends on the modelType
    //The Z-ga decay width of the CP-odd Higgs boson
    Gamma_AZga=0.0;//It depends on the modelType
    //The g-g decay width of the CP-odd Higgs boson
    Gamma_Agg=0.0;//It depends on the modelType

    /*Calulation of rA_QdQd, rA_ll, rA_gg, Gamma_Agaga, Gamma_AZga, Gamma_Agg
     * (they depend on the model type): START*/

    /*Gamma_Agaga and Gamma_AZga expressions can be found in
     "The Higgs Hunter's Guide", Appendix C and in arXiv:0902.4665v3, Appendix A;
     *Gamma_Agg expression can be found in arXiv:0902.4665v3, Appendix A*/

    /*TAUi is necessary for Gamma_Agaga, Gamma_AZga and Gamma_Agg*/
    TAUc=4.0*Mc*Mc/(mA*mA);
    TAUt=4.0*Mt*Mt/(mA*mA);
    TAUs=4.0*Ms*Ms/(mA*mA);
    TAUb=4.0*Mb*Mb/(mA*mA);
    TAUmu=4.0*Mmu*Mmu/(mA*mA);
    TAUtau=4.0*Mtau*Mtau/(mA*mA);
    TAUw=4.0*MW*MW/(mA*mA);
    TAUhp=4.0*mHp*mHp/(mA*mA);

    /*I_A_f, I_A_W and I_A_Hp are necessary for Gamma_Agaga, Gamma_AZga and Gamma_Agg;
     * their expressions can be found in "The Higgs Hunter's Guide", Appendix C, C.4*/
    I_A_f=0.0;//It depends on the modelType
    I_A_fU=-(8./3.)*(//TAUu*myfunctions.f_func(TAUu)+
                          TAUc*myfunctions.f_func(TAUc)+TAUt*myfunctions.f_func(TAUt));
    I_A_fD=-(2./3.)*(//TAUd*myfunctions.f_func(TAUd)+
                          TAUs*myfunctions.f_func(TAUs)+TAUb*myfunctions.f_func(TAUb));
    I_A_fL=-2.*(//TAUe*myfunctions.f_func(TAUe)+
                          TAUmu*myfunctions.f_func(TAUmu)+TAUtau*myfunctions.f_func(TAUtau));

    //LAMi is necessary for Gamma_HZga
    LAMc=4.0*Mc*Mc/(MZ*MZ);
    LAMt=4.0*Mt*Mt/(MZ*MZ);
    LAMs=4.0*Ms*Ms/(MZ*MZ);
    LAMb=4.0*Mb*Mb/(MZ*MZ);
    LAMmu=4.0*Mmu*Mmu/(MZ*MZ);
    LAMtau=4.0*Mtau*Mtau/(MZ*MZ);
    LAMw=4.0*MW*MW/(MZ*MZ);
    LAMhp=4.0*mHp*mHp/(MZ*MZ);

    /*A_A_F, A_A_W and A_A_Hp are necessary for Gamma_AZga*/
    /*A_A_F expression can be found in "The Higgs Hunter's Guide", Appendix C, C.12*/
    A_A_F = 0.0;//It depends on the modelType
    A_A_U = -4.0*(1.0/2.0-4.0/3.0*s02)
                           *(//-myfunctions.Int2(TAUu,LAMu)
                           -myfunctions.Int2(TAUc,LAMc)-myfunctions.Int2(TAUt,LAMt))/sqrt(s02*c02);
    A_A_D = +2.0*(-1.0/2.0+2.0/3.0*s02)
                           *(//-myfunctions.Int2(TAUd,LAMd)
                           -myfunctions.Int2(TAUs,LAMs)-myfunctions.Int2(TAUb,LAMb))/sqrt(s02*c02);
    A_A_L  = +2.0*(-1.0/2.0+2.0*s02)
                            *(//-myfunctions.Int2(TAUe,LAMe)
                            -myfunctions.Int2(TAUmu,LAMmu)-myfunctions.Int2(TAUtau,LAMtau))/sqrt(s02*c02);

    std::string modelflag=myTHDM->getModelTypeflag();

    if( modelflag == "type1" ) {
        rA_gg=-cosb/sinb*cosb/sinb+2.0*cosb/sinb*cosb/sinb*(Sigmaggh_tt+Sigmaggh_bb)/SigmaggF;
        rA_QdQd=cosb/sinb*cosb/sinb;
        rA_ll=cosb/sinb*cosb/sinb;
        I_A_f=cosb/sinb*(I_A_fU-I_A_fD-I_A_fL);
        A_A_F = cosb/sinb*(A_A_U-A_A_D-A_A_L);
    }
    else if( modelflag == "type2" ) {
        rA_gg= 1.0+(cosb/sinb-sinb/cosb)*(Sigmaggh_tt*cosb/sinb-Sigmaggh_bb*sinb/cosb)/SigmaggF;
        rA_QdQd=sinb/cosb*sinb/cosb;
        rA_ll=sinb/cosb*sinb/cosb;
        I_A_f=cosb/sinb*I_A_fU+sinb/cosb*(I_A_fD+I_A_fL);
        A_A_F = (cosb/sinb*A_A_U+sinb/cosb*(A_A_D+A_A_L));
    }
    else if( modelflag == "typeX" ) {
        rA_gg=-cosb/sinb*cosb/sinb+2.0*cosb/sinb*cosb/sinb*(Sigmaggh_tt+Sigmaggh_bb)/SigmaggF;
        rA_QdQd=cosb/sinb*cosb/sinb;
        rA_ll=sinb/cosb*sinb/cosb;
        I_A_f = cosb/sinb*(I_A_fU-I_A_fD)+sinb/cosb*I_A_fL;
        A_A_F = (cosb/sinb*(A_A_U-A_A_D)+sinb/cosb*A_A_L);
    }
    else if( modelflag == "typeY" ) {
        rA_gg=1.0+(cosb/sinb-sinb/cosb)*(Sigmaggh_tt*cosb/sinb-Sigmaggh_bb*sinb/cosb)/SigmaggF;
        rA_QdQd=sinb/cosb*sinb/cosb;
        rA_ll=cosb/sinb*cosb/sinb;
        I_A_f = cosb/sinb*(I_A_fU-I_A_fL)+sinb/cosb*I_A_fD;
        A_A_F = (cosb/sinb*(A_A_U-A_A_L)+sinb/cosb*A_A_D);
    }
    else {
        throw std::runtime_error("modelflag can be only any of \"type1\", \"type2\", \"typeX\" or \"typeY\"");
    }

//    std::cout<<"rA_gg: "<<rA_gg<<std::endl;

    /*Gamma_Agaga expression can be found in in arXiv:0902.4665v3, Appendix A, A.8*/
    Gamma_Agaga=GF*Ale*Ale*mA*mA*mA/(sqrt(2)*128.0*M_PI*M_PI*M_PI)
                *(I_A_f).abs()*(I_A_f).abs();
    /*Gamma_AZga expression can be found in in arXiv:0902.4665v3, Appendix A, A.9*/
    Gamma_AZga=GF*Ale*Ale*mA*mA*mA/(sqrt(2)*64.0*M_PI*M_PI*M_PI)
               *(1.0-MZ*MZ/(mA*mA))*(1.0-MZ*MZ/(mA*mA))*(1.0-MZ*MZ/(mA*mA))
               *(A_A_F).abs()*(A_A_F).abs();
    /*Gamma_Agg expression can be found in in arXiv:0902.4665v3, Appendix A, A.10*/
    Gamma_Agg=GF*Als*Als*mA*mA*mA/(sqrt(2)*64.0*M_PI*M_PI*M_PI)*rA_gg;

    /*Calulation of rA_QdQd, rA_ll, rA_gg, Gamma_Hgaga, Gamma_HZga, Gamma_Hgg
     * (they depend on the model type): FINISH*/

    SigmaggF_A=mycache->cs_ggFtoA(mA)*rA_gg;
    SigmabbF_A=mycache->cs_bbFtoHP(mA)*rA_QdQd;

    SigmaSum = SigmaggF_A + SigmabbF_A; //+ SigmattF_A ;

    BrSM_Atocc=mycache->Br_HPtocc(mA);
    BrSM_Atobb=mycache->Br_HPtobb(mA);
    BrSM_Atott=mycache->Br_HPtott(mA);
    BrSM_Atomumu=mycache->Br_HPtomumu(mA);
    BrSM_Atotautau=mycache->Br_HPtotautau(mA);

    GammaAtotSM=mycache->GammaHPtotSM(mA);

    GammaAHZ=myfunctions.HSTheta(mA-MZ-mHh)*pow(myfunctions.KaellenFunction(mA,MZ,mHh),3)
                    *sin_ba*sin_ba/(2.0*M_PI*vev*vev);

    GammaAhZ=myfunctions.HSTheta(mA-MZ-mHl)*pow(myfunctions.KaellenFunction(mA,MZ,mHl),3)
                    *cos_ab*cos_ab/(2.0*M_PI*vev*vev);

    GammaAHpW=2.*myfunctions.HSTheta(mA-MW-mHp)*pow(myfunctions.KaellenFunction(mA,MW,mHp),3)
                     /(2.0*M_PI*vev*vev);

    GammaAtot= ((BrSM_Atott+BrSM_Atocc)*rA_QuQu
                    +BrSM_Atobb*rA_QdQd
                    +(BrSM_Atotautau+BrSM_Atomumu)*rA_ll)*GammaAtotSM
               +GammaAHZ+GammaAhZ+GammaAHpW+Gamma_Agaga+Gamma_AZga+Gamma_Agg;

    Br_Atott=BrSM_Atott*rA_QuQu*GammaAtotSM/GammaAtot;
    Br_Atobb=BrSM_Atobb*rA_QdQd*GammaAtotSM/GammaAtot;
    Br_Atotautau=BrSM_Atotautau*rA_ll*GammaAtotSM/GammaAtot;
    Br_Atogaga=Gamma_Agaga/GammaAtot;
    Br_AtohZ=GammaAhZ/GammaAtot;

    Br_htobb=mylightHiggs->THDM_BR_h_bb();
    Br_htotautau=mylightHiggs->THDM_BR_h_tautau();

    Br_Ztoee=3.363*0.01; //K.A. Olive et al. (Particle Data Group), Chin. Phys. C38, 090001 (2014)
    Br_Ztomumu=3.366*0.01;
    Br_Ztotautau=3.370*0.01;

    //CP-odd Higgs Signals, theoretical expressions

    ggF_A_tautau_TH=SigmaggF_A*Br_Atotautau;
    bbF_A_tautau_TH=SigmabbF_A*Br_Atotautau;
    ggF_A_gaga_TH=SigmaggF_A*Br_Atogaga;
    ggF_A_hZ_bbll_TH=SigmaggF_A*Br_AtohZ*Br_htobb*(Br_Ztoee+Br_Ztomumu);
    ggF_A_hZ_bbZ_TH=SigmaggF_A*Br_AtohZ*Br_htobb;
    ggF_A_hZ_tautaull_TH=SigmaggF_A*Br_AtohZ*Br_htotautau*(Br_Ztoee+Br_Ztomumu+Br_Ztotautau);
    ggF_A_hZ_tautauZ_TH=SigmaggF_A*Br_AtohZ*Br_htotautau;
    pp_A_tt_TH=SigmaSum*Br_Atott;
    bbF_A_bb_TH=SigmabbF_A*Br_Atobb;   
    
    //CP-odd Higgs Signals, experimental expressions

    ggF_A_tautau_EX_ATLAS=mycache->ex_ggF_phi_tautau_ATLAS(mA);
    ggF_A_tautau_EX_CMS=mycache->ex_ggF_phi_tautau_CMS(mA);
    bbF_A_tautau_EX_ATLAS=mycache->ex_bbF_phi_tautau_ATLAS(mA);
    bbF_A_tautau_EX_CMS=mycache->ex_bbF_phi_tautau_CMS(mA);
    ggF_A_gaga_EX_ATLAS=mycache->ex_ggF_phi_gaga_ATLAS(mA);
    ggF_A_gaga_EX_CMS=mycache->ex_ggF_phi_gaga_CMS(mA);
    ggF_A_hZ_bbll_EX_CMS=mycache->ex_ggF_A_hZ_bbll_CMS(mA);
    ggF_A_hZ_bbZ_EX_ATLAS=mycache->ex_ggF_A_hZ_bbZ_ATLAS(mA);
    ggF_A_hZ_tautaull_EX_CMS=mycache->ex_ggF_A_hZ_tautaull_CMS(mA);
    ggF_A_hZ_tautauZ_EX_ATLAS=mycache->ex_ggF_A_hZ_tautauZ_ATLAS(mA);
    pp_A_tt_EX_ATLAS=mycache->ex_pp_phi_tt_ATLAS(mA);
    bbF_A_bb_EX_CMS=mycache->ex_bbF_phi_bb_CMS(mA);
}


double CPoddHiggsCache::computeThValue()
{
    return 0;
}



/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/



Hobs_ggF_A_tautau_ATLAS::Hobs_ggF_A_tautau_ATLAS(const StandardModel& SM_i)
: CPoddHiggsCache(SM_i)
{}

double Hobs_ggF_A_tautau_ATLAS::computeThValue()
{
    computeParameters();
    return ggF_A_tautau_TH/ggF_A_tautau_EX_ATLAS;
}



Hobs_ggF_A_tautau_CMS::Hobs_ggF_A_tautau_CMS(const StandardModel& SM_i)
: CPoddHiggsCache(SM_i)
{}

double Hobs_ggF_A_tautau_CMS::computeThValue()
{
    computeParameters();
    return ggF_A_tautau_TH/ggF_A_tautau_EX_CMS;
}



Hobs_bbF_A_tautau_ATLAS::Hobs_bbF_A_tautau_ATLAS(const StandardModel& SM_i)
: CPoddHiggsCache(SM_i)
{}

double Hobs_bbF_A_tautau_ATLAS::computeThValue()
{
    computeParameters();
    return bbF_A_tautau_TH/bbF_A_tautau_EX_ATLAS;
}



Hobs_bbF_A_tautau_CMS::Hobs_bbF_A_tautau_CMS(const StandardModel& SM_i)
: CPoddHiggsCache(SM_i)
{}

double Hobs_bbF_A_tautau_CMS::computeThValue()
{
    computeParameters();
    return bbF_A_tautau_TH/bbF_A_tautau_EX_CMS;
}



Hobs_ggF_A_gaga_ATLAS::Hobs_ggF_A_gaga_ATLAS(const StandardModel& SM_i)
: CPoddHiggsCache(SM_i)
{}

double Hobs_ggF_A_gaga_ATLAS::computeThValue()
{
    computeParameters();
    return ggF_A_gaga_TH/ggF_A_gaga_EX_ATLAS;
}



Hobs_ggF_A_gaga_CMS::Hobs_ggF_A_gaga_CMS(const StandardModel& SM_i)
: CPoddHiggsCache(SM_i)
{}

double Hobs_ggF_A_gaga_CMS::computeThValue()
{
    computeParameters();
    return ggF_A_gaga_TH/ggF_A_gaga_EX_CMS;
}



Hobs_ggF_A_hZ_bbll_CMS::Hobs_ggF_A_hZ_bbll_CMS(const StandardModel& SM_i)
: CPoddHiggsCache(SM_i)
{}

double Hobs_ggF_A_hZ_bbll_CMS::computeThValue()
{
    computeParameters();
    return ggF_A_hZ_bbll_TH/ggF_A_hZ_bbll_EX_CMS;
}



Hobs_ggF_A_hZ_bbZ_ATLAS::Hobs_ggF_A_hZ_bbZ_ATLAS(const StandardModel& SM_i)
: CPoddHiggsCache(SM_i)
{}

double Hobs_ggF_A_hZ_bbZ_ATLAS::computeThValue()
{
    computeParameters();
    return ggF_A_hZ_bbZ_TH/ggF_A_hZ_bbZ_EX_ATLAS;
}



Hobs_ggF_A_hZ_tautaull_CMS::Hobs_ggF_A_hZ_tautaull_CMS(const StandardModel& SM_i)
: CPoddHiggsCache(SM_i)
{}

double Hobs_ggF_A_hZ_tautaull_CMS::computeThValue()
{
    computeParameters();
    return ggF_A_hZ_tautaull_TH/ggF_A_hZ_tautaull_EX_CMS;
}



Hobs_ggF_A_hZ_tautauZ_ATLAS::Hobs_ggF_A_hZ_tautauZ_ATLAS(const StandardModel& SM_i)
: CPoddHiggsCache(SM_i)
{}

double Hobs_ggF_A_hZ_tautauZ_ATLAS::computeThValue()
{
    computeParameters();
    return ggF_A_hZ_tautauZ_TH/ggF_A_hZ_tautauZ_EX_ATLAS;
}



Hobs_pp_A_tt_ATLAS::Hobs_pp_A_tt_ATLAS(const StandardModel& SM_i)
: CPoddHiggsCache(SM_i)
{}

double Hobs_pp_A_tt_ATLAS::computeThValue()
{
    computeParameters();
    return pp_A_tt_TH/pp_A_tt_EX_ATLAS;
}



Hobs_bbF_A_bb_CMS::Hobs_bbF_A_bb_CMS(const StandardModel& SM_i)
: CPoddHiggsCache(SM_i)
{}

double Hobs_bbF_A_bb_CMS::computeThValue()
{
    computeParameters();
    return bbF_A_bb_TH/bbF_A_bb_EX_CMS;
}



log10_ggF_A_tautau_TH::log10_ggF_A_tautau_TH(const StandardModel& SM_i)
: CPoddHiggsCache(SM_i)
{}

double log10_ggF_A_tautau_TH::computeThValue()
{
    computeParameters();
    return log10(ggF_A_tautau_TH);
}



log10_bbF_A_tautau_TH::log10_bbF_A_tautau_TH(const StandardModel& SM_i)
: CPoddHiggsCache(SM_i)
{}

double log10_bbF_A_tautau_TH::computeThValue()
{
    computeParameters();
    return log10(bbF_A_tautau_TH);
}



log10_ggF_A_gaga_TH::log10_ggF_A_gaga_TH(const StandardModel& SM_i)
: CPoddHiggsCache(SM_i)
{}

double log10_ggF_A_gaga_TH::computeThValue()
{
    computeParameters();
    return log10(ggF_A_gaga_TH);
}



log10_ggF_A_hZ_bbll_TH::log10_ggF_A_hZ_bbll_TH(const StandardModel& SM_i)
: CPoddHiggsCache(SM_i)
{}

double log10_ggF_A_hZ_bbll_TH::computeThValue()
{
    computeParameters();
    return log10(ggF_A_hZ_bbll_TH);
}



log10_ggF_A_hZ_bbZ_TH::log10_ggF_A_hZ_bbZ_TH(const StandardModel& SM_i)
: CPoddHiggsCache(SM_i)
{}

double log10_ggF_A_hZ_bbZ_TH::computeThValue()
{
    computeParameters();
    return log10(ggF_A_hZ_bbZ_TH);
}



log10_ggF_A_hZ_tautaull_TH::log10_ggF_A_hZ_tautaull_TH(const StandardModel& SM_i)
  : CPoddHiggsCache(SM_i)
{}

double log10_ggF_A_hZ_tautaull_TH::computeThValue()
{
  computeParameters();
  return log10(ggF_A_hZ_tautaull_TH);
}



log10_ggF_A_hZ_tautauZ_TH::log10_ggF_A_hZ_tautauZ_TH(const StandardModel& SM_i)
: CPoddHiggsCache(SM_i)
{}

double log10_ggF_A_hZ_tautauZ_TH::computeThValue()
{
    computeParameters();
    return log10(ggF_A_hZ_tautauZ_TH);
}



log10_pp_A_tt_TH::log10_pp_A_tt_TH(const StandardModel& SM_i)
: CPoddHiggsCache(SM_i)
{}

double log10_pp_A_tt_TH::computeThValue()
{
    computeParameters();
    return log10(pp_A_tt_TH);
}



log10_bbF_A_bb_TH::log10_bbF_A_bb_TH(const StandardModel& SM_i)
: CPoddHiggsCache(SM_i)
{}

double log10_bbF_A_bb_TH::computeThValue()
{
    computeParameters();
    return log10(bbF_A_bb_TH);
}
