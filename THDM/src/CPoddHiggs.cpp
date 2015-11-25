/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "CPoddHiggs.h"
#include "lightHiggs.h"
#include "StandardModel.h"

CPoddHiggs::CPoddHiggs(const StandardModel& SM_i):

        ThObservable(SM_i),
        myTHDM(static_cast<const THDM*> (&SM_i)),
        mySM (SM_i)
{
    mycache = new THDMcache();
    mylightHiggs = new lightHiggs(SM_i);
}

CPoddHiggs::~CPoddHiggs()
{}


void CPoddHiggs::computeParameters()
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
    double sW2=myTHDM->sW2();
    double cW2=myTHDM->cW2();
    double mHh=myTHDM->getmHh();
    double mHl=myTHDM->getMHl();
    double mHp=myTHDM->getmHp();
    mA=myTHDM->getmA();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double sin_ba=myTHDM->getsin_ba();
    double cos_ab=cosa*cosb+sina*sinb;

    //These cross sections are necessary for rA_gg
    //SM gg -> A production cross section at 8 TeV
    double SigmaggA = mycache->ip_cs_ggA(mA);
    //SM gg -> A productmycaion cross section at 8 TeV, top loop only
    double SigmaggA_tt = mycache->ip_cs_ggA_tt(mA);
    //SM gg -> A production cross section at 8 TeV, bottom loop only
    double SigmaggA_bb = mycache->ip_cs_ggA_bb(mA);

    /* r_ii is the ratio between the squared 2HDM vertex coupling of the CP-odd
     * Higgs to the particle i and the corresponding coupling of the SM Higgs boson.*/
    double rA_QuQu=(cosb*cosb)/(sinb*sinb);
    double rA_QdQd=0.0;//It depends on the modelType
    double rA_ll=0.0;//It depends on the modelType
    double rA_gg=0.0;//It depends on the modelType

    /*Calulation of rA_QdQd, rA_ll, rA_gg, Gamma_Agaga, Gamma_AZga, Gamma_Agg
     * (they depend on the model type): START*/

    /*Gamma_Agaga and Gamma_AZga expressions can be found in
     "The Higgs Hunter's Guide", Appendix C and in arXiv:0902.4665v3, Appendix A;
     *Gamma_Agg expression can be found in arXiv:0902.4665v3, Appendix A*/

    double TAUc=4.0*Mc*Mc/(mA*mA);
    double TAUt=4.0*Mt*Mt/(mA*mA);
    double TAUs=4.0*Ms*Ms/(mA*mA);
    double TAUb=4.0*Mb*Mb/(mA*mA);
    double TAUmu=4.0*Mmu*Mmu/(mA*mA);
    double TAUtau=4.0*Mtau*Mtau/(mA*mA);

    /*I_A_F is needed for Gamma_Agaga;
     * The expression can be found in "The Higgs Hunter's Guide", Appendix C, C.4*/
    gslpp::complex I_A_F=0.0;//It depends on the modelType
    gslpp::complex I_A_U=-(8./3.)*(TAUc*myfunctions.f_func(TAUc)+TAUt*myfunctions.f_func(TAUt));
    gslpp::complex I_A_D=-(2./3.)*(TAUs*myfunctions.f_func(TAUs)+TAUb*myfunctions.f_func(TAUb));
    gslpp::complex I_A_L=-2.*(TAUmu*myfunctions.f_func(TAUmu)+TAUtau*myfunctions.f_func(TAUtau));

    double LAMc=4.0*Mc*Mc/(MZ*MZ);
    double LAMt=4.0*Mt*Mt/(MZ*MZ);
    double LAMs=4.0*Ms*Ms/(MZ*MZ);
    double LAMb=4.0*Mb*Mb/(MZ*MZ);
    double LAMmu=4.0*Mmu*Mmu/(MZ*MZ);
    double LAMtau=4.0*Mtau*Mtau/(MZ*MZ);

    /*A_A_F is needed for Gamma_AZga*/
    /*The expression can be found in "The Higgs Hunter's Guide", Appendix C, C.12*/
    gslpp::complex A_A_F = 0.0;//It depends on the modelType
    gslpp::complex A_A_U = -4.0*(1.0/2.0-4.0/3.0*sW2)*(-myfunctions.Int2(TAUc,LAMc)-myfunctions.Int2(TAUt,LAMt))/sqrt(sW2*cW2);
    gslpp::complex A_A_D = 2.0*(-1.0/2.0+2.0/3.0*sW2)*(-myfunctions.Int2(TAUs,LAMs)-myfunctions.Int2(TAUb,LAMb))/sqrt(sW2*cW2);
    gslpp::complex A_A_L = 2.0*(-1.0/2.0+2.0*sW2)*(-myfunctions.Int2(TAUmu,LAMmu)-myfunctions.Int2(TAUtau,LAMtau))/sqrt(sW2*cW2);

    std::string modelflag=myTHDM->getModelTypeflag();

    if( modelflag == "type1" ) {
        rA_gg=-cosb/sinb*cosb/sinb+2.0*cosb/sinb*cosb/sinb*(SigmaggA_tt+SigmaggA_bb)/SigmaggA;
        rA_QdQd=cosb/sinb*cosb/sinb;
        rA_ll=cosb/sinb*cosb/sinb;
        I_A_F=cosb/sinb*(I_A_U-I_A_D-I_A_L);
        A_A_F=cosb/sinb*(A_A_U-A_A_D-A_A_L);
    }
    else if( modelflag == "type2" ) {
        rA_gg= 1.0+(cosb/sinb-sinb/cosb)*(SigmaggA_tt*cosb/sinb-SigmaggA_bb*sinb/cosb)/SigmaggA;
        rA_QdQd=sinb/cosb*sinb/cosb;
        rA_ll=sinb/cosb*sinb/cosb;
        I_A_F=cosb/sinb*I_A_U+sinb/cosb*(I_A_D+I_A_L);
        A_A_F=cosb/sinb*A_A_U+sinb/cosb*(A_A_D+A_A_L);
    }
    else if( modelflag == "typeX" ) {
        rA_gg=-cosb/sinb*cosb/sinb+2.0*cosb/sinb*cosb/sinb*(SigmaggA_tt+SigmaggA_bb)/SigmaggA;
        rA_QdQd=cosb/sinb*cosb/sinb;
        rA_ll=sinb/cosb*sinb/cosb;
        I_A_F=cosb/sinb*(I_A_U-I_A_D)+sinb/cosb*I_A_L;
        A_A_F=cosb/sinb*(A_A_U-A_A_D)+sinb/cosb*A_A_L;
    }
    else if( modelflag == "typeY" ) {
        rA_gg=1.0+(cosb/sinb-sinb/cosb)*(SigmaggA_tt*cosb/sinb-SigmaggA_bb*sinb/cosb)/SigmaggA;
        rA_QdQd=sinb/cosb*sinb/cosb;
        rA_ll=cosb/sinb*cosb/sinb;
        I_A_F=cosb/sinb*(I_A_U-I_A_L)+sinb/cosb*I_A_D;
        A_A_F=cosb/sinb*(A_A_U-A_A_L)+sinb/cosb*A_A_D;
    }
    else {
        throw std::runtime_error("modelflag can be only any of \"type1\", \"type2\", \"typeX\" or \"typeY\"");
    }

    /*Gamma_Agaga expression can be found in in arXiv:0902.4665v3, Appendix A, A.8*/
    double Gamma_Agaga=GF*Ale*Ale*mA*mA*mA/(sqrt(2)*128.0*M_PI*M_PI*M_PI)
                *(I_A_F).abs2();
    /*Gamma_AZga expression can be found in in arXiv:0902.4665v3, Appendix A, A.9*/
    double Gamma_AZga=GF*Ale*Ale*mA*mA*mA/(sqrt(2)*64.0*M_PI*M_PI*M_PI)
               *(1.0-MZ*MZ/(mA*mA))*(1.0-MZ*MZ/(mA*mA))*(1.0-MZ*MZ/(mA*mA))
               *(A_A_F).abs2();
    /*Gamma_Agg expression can be found in in arXiv:0902.4665v3, Appendix A, A.10*/
    double Gamma_Agg=GF*Als*Als*mA*mA*mA/(sqrt(2)*64.0*M_PI*M_PI*M_PI)*rA_gg;

    /*Calulation of rA_QdQd, rA_ll, rA_gg, Gamma_Agaga, Gamma_AZga, Gamma_Agg: END*/

    double SigmaggF_A=SigmaggA*rA_gg;
    double SigmabbF_A=mycache->ip_cs_bbFtoHP(mA)*rA_QdQd;
    double SigmaSum = SigmaggF_A + SigmabbF_A; //+ SigmattF_A ;

    double BrSM_Atocc=mycache->ip_Br_HPtocc(mA);
    double BrSM_Atobb=mycache->ip_Br_HPtobb(mA);
    double BrSM_Atott=mycache->ip_Br_HPtott(mA);
    double BrSM_Atomumu=mycache->ip_Br_HPtomumu(mA);
    double BrSM_Atotautau=mycache->ip_Br_HPtotautau(mA);

    double GammaAtotSM=mycache->ip_GammaHPtotSM(mA);

    double GammaAHZ=myfunctions.HSTheta(mA-MZ-mHh)*pow(myfunctions.KaellenFunction(mA,MZ,mHh),3)
                    *sin_ba*sin_ba/(2.0*M_PI*vev*vev);

    double GammaAhZ=myfunctions.HSTheta(mA-MZ-mHl)*pow(myfunctions.KaellenFunction(mA,MZ,mHl),3)
                    *cos_ab*cos_ab/(2.0*M_PI*vev*vev);

    double GammaAHpW=2.*myfunctions.HSTheta(mA-MW-mHp)*pow(myfunctions.KaellenFunction(mA,MW,mHp),3)
                     /(2.0*M_PI*vev*vev);

    GammaAtot= ((BrSM_Atott+BrSM_Atocc)*rA_QuQu
                    +BrSM_Atobb*rA_QdQd
                    +(BrSM_Atotautau+BrSM_Atomumu)*rA_ll)*GammaAtotSM
               +Gamma_Agaga+Gamma_AZga+Gamma_Agg+GammaAHZ+GammaAhZ+GammaAHpW;

    double Br_Atott=BrSM_Atott*rA_QuQu*GammaAtotSM/GammaAtot;
    double Br_Atobb=BrSM_Atobb*rA_QdQd*GammaAtotSM/GammaAtot;
    double Br_Atotautau=BrSM_Atotautau*rA_ll*GammaAtotSM/GammaAtot;
    double Br_Atogaga=Gamma_Agaga/GammaAtot;
    double Br_AtohZ=GammaAhZ/GammaAtot;

    double Br_htobb=mylightHiggs->THDM_BR_h_bb();
    double Br_htotautau=mylightHiggs->THDM_BR_h_tautau();

    double Br_Ztoee=3.363*0.01; //K.A. Olive et al. (Particle Data Group), Chin. Phys. C38, 090001 (2014)
    double Br_Ztomumu=3.366*0.01; //K.A. Olive et al. (Particle Data Group), Chin. Phys. C38, 090001 (2014)
    double Br_Ztotautau=3.370*0.01; //K.A. Olive et al. (Particle Data Group), Chin. Phys. C38, 090001 (2014)

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
}


double CPoddHiggs::computeThValue()
{
    return 0;
}



/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/



Hobs_ggF_A_tautau_ATLAS::Hobs_ggF_A_tautau_ATLAS(const StandardModel& SM_i)
: CPoddHiggs(SM_i)
{}

double Hobs_ggF_A_tautau_ATLAS::computeThValue()
{
    if(mA>=100.0 && mA<=1000.0)
    {
        computeParameters();
        return ggF_A_tautau_TH/mycache->ip_ex_ggF_phi_tautau_ATLAS(mA);
    }
    else
    {
        return 0.0;
    }
}



Hobs_ggF_A_tautau_CMS::Hobs_ggF_A_tautau_CMS(const StandardModel& SM_i)
: CPoddHiggs(SM_i)
{}

double Hobs_ggF_A_tautau_CMS::computeThValue()
{
    if(mA>=90.0 && mA<=1000.0)
    {
        computeParameters();
        return ggF_A_tautau_TH/mycache->ip_ex_ggF_phi_tautau_CMS(mA);
    }
    else
    {
        return 0.0;
    }
}



Hobs_bbF_A_tautau_ATLAS::Hobs_bbF_A_tautau_ATLAS(const StandardModel& SM_i)
: CPoddHiggs(SM_i)
{}

double Hobs_bbF_A_tautau_ATLAS::computeThValue()
{
    if(mA>=100.0 && mA<=1000.0)
    {
        computeParameters();
        return bbF_A_tautau_TH/mycache->ip_ex_bbF_phi_tautau_ATLAS(mA);
    }
    else
    {
        return 0.0;
    }
}



Hobs_bbF_A_tautau_CMS::Hobs_bbF_A_tautau_CMS(const StandardModel& SM_i)
: CPoddHiggs(SM_i)
{}

double Hobs_bbF_A_tautau_CMS::computeThValue()
{
    if(mA>=90.0 && mA<=1000.0)
    {
        computeParameters();
        return bbF_A_tautau_TH/mycache->ip_ex_bbF_phi_tautau_CMS(mA);
    }
    else
    {
        return 0.0;
    }
}



Hobs_ggF_A_gaga_ATLAS::Hobs_ggF_A_gaga_ATLAS(const StandardModel& SM_i)
: CPoddHiggs(SM_i)
{}

double Hobs_ggF_A_gaga_ATLAS::computeThValue()
{
    if(mA>=100.0 && mA<=600.0)
    {
        computeParameters();
        return ggF_A_gaga_TH/mycache->ip_ex_ggF_phi_gaga_ATLAS(mA);
    }
    else
    {
        return 0.0;
    }
}



Hobs_ggF_A_gaga_CMS::Hobs_ggF_A_gaga_CMS(const StandardModel& SM_i)
: CPoddHiggs(SM_i)
{}

double Hobs_ggF_A_gaga_CMS::computeThValue()
{
    if(mA>=150.0 && mA<=850.0)
    {
        computeParameters();
        return ggF_A_gaga_TH/mycache->ip_ex_ggF_phi_gaga_CMS(mA);
    }
    else
    {
        return 0.0;
    }
}



Hobs_ggF_A_hZ_bbll_CMS::Hobs_ggF_A_hZ_bbll_CMS(const StandardModel& SM_i)
: CPoddHiggs(SM_i)
{}

double Hobs_ggF_A_hZ_bbll_CMS::computeThValue()
{
    if(mA>=170.0 && mA<=600.0)
    {
        computeParameters();
        return ggF_A_hZ_bbll_TH/mycache->ip_ex_ggF_A_hZ_bbll_CMS(mA);
    }
    else
    {
        return 0.0;
    }
}



Hobs_ggF_A_hZ_bbZ_ATLAS::Hobs_ggF_A_hZ_bbZ_ATLAS(const StandardModel& SM_i)
: CPoddHiggs(SM_i)
{}

double Hobs_ggF_A_hZ_bbZ_ATLAS::computeThValue()
{
    if(mA>=220.0 && mA<=1000.0)
    {
        computeParameters();
        return ggF_A_hZ_bbZ_TH/mycache->ip_ex_ggF_A_hZ_bbZ_ATLAS(mA);
    }
    else
    {
        return 0.0;
    }
}



Hobs_ggF_A_hZ_tautaull_CMS::Hobs_ggF_A_hZ_tautaull_CMS(const StandardModel& SM_i)
: CPoddHiggs(SM_i)
{}

double Hobs_ggF_A_hZ_tautaull_CMS::computeThValue()
{
    if(mA>=220.0 && mA<=350.0)
    {
        computeParameters();
        return ggF_A_hZ_tautaull_TH/mycache->ip_ex_ggF_A_hZ_tautaull_CMS(mA);
    }
    else
    {
        return 0.0;
    }
}



Hobs_ggF_A_hZ_tautauZ_ATLAS::Hobs_ggF_A_hZ_tautauZ_ATLAS(const StandardModel& SM_i)
: CPoddHiggs(SM_i)
{}

double Hobs_ggF_A_hZ_tautauZ_ATLAS::computeThValue()
{
    if(mA>=220.0 && mA<=1000.0)
    {
        computeParameters();
        return ggF_A_hZ_tautauZ_TH/mycache->ip_ex_ggF_A_hZ_tautauZ_ATLAS(mA);
    }
    else
    {
        return 0.0;
    }
}



Hobs_pp_A_tt_ATLAS::Hobs_pp_A_tt_ATLAS(const StandardModel& SM_i)
: CPoddHiggs(SM_i)
{}

double Hobs_pp_A_tt_ATLAS::computeThValue()
{
    if(mA>=400.0 && mA<=3000.0)
    {
        computeParameters();
        return pp_A_tt_TH/mycache->ip_ex_pp_phi_tt_ATLAS(mA);
    }
    else
    {
        return 0.0;
    }
}



Hobs_bbF_A_bb_CMS::Hobs_bbF_A_bb_CMS(const StandardModel& SM_i)
: CPoddHiggs(SM_i)
{}

double Hobs_bbF_A_bb_CMS::computeThValue()
{
    if(mA>=100.0 && mA<=900.0)
    {
        computeParameters();
        return bbF_A_bb_TH/mycache->ip_ex_bbF_phi_bb_CMS(mA);
    }
    else
    {
        return 0.0;
    }
}



log10_ggF_A_tautau_TH::log10_ggF_A_tautau_TH(const StandardModel& SM_i)
: CPoddHiggs(SM_i)
{}

double log10_ggF_A_tautau_TH::computeThValue()
{
    computeParameters();
    return log10(ggF_A_tautau_TH);
}



log10_bbF_A_tautau_TH::log10_bbF_A_tautau_TH(const StandardModel& SM_i)
: CPoddHiggs(SM_i)
{}

double log10_bbF_A_tautau_TH::computeThValue()
{
    computeParameters();
    return log10(bbF_A_tautau_TH);
}



log10_ggF_A_gaga_TH::log10_ggF_A_gaga_TH(const StandardModel& SM_i)
: CPoddHiggs(SM_i)
{}

double log10_ggF_A_gaga_TH::computeThValue()
{
    computeParameters();
    return log10(ggF_A_gaga_TH);
}



log10_ggF_A_hZ_bbll_TH::log10_ggF_A_hZ_bbll_TH(const StandardModel& SM_i)
: CPoddHiggs(SM_i)
{}

double log10_ggF_A_hZ_bbll_TH::computeThValue()
{
    computeParameters();
    return log10(ggF_A_hZ_bbll_TH);
}



log10_ggF_A_hZ_bbZ_TH::log10_ggF_A_hZ_bbZ_TH(const StandardModel& SM_i)
: CPoddHiggs(SM_i)
{}

double log10_ggF_A_hZ_bbZ_TH::computeThValue()
{
    computeParameters();
    return log10(ggF_A_hZ_bbZ_TH);
}



log10_ggF_A_hZ_tautaull_TH::log10_ggF_A_hZ_tautaull_TH(const StandardModel& SM_i)
  : CPoddHiggs(SM_i)
{}

double log10_ggF_A_hZ_tautaull_TH::computeThValue()
{
  computeParameters();
  return log10(ggF_A_hZ_tautaull_TH);
}



log10_ggF_A_hZ_tautauZ_TH::log10_ggF_A_hZ_tautauZ_TH(const StandardModel& SM_i)
: CPoddHiggs(SM_i)
{}

double log10_ggF_A_hZ_tautauZ_TH::computeThValue()
{
    computeParameters();
    return log10(ggF_A_hZ_tautauZ_TH);
}



log10_pp_A_tt_TH::log10_pp_A_tt_TH(const StandardModel& SM_i)
: CPoddHiggs(SM_i)
{}

double log10_pp_A_tt_TH::computeThValue()
{
    computeParameters();
    return log10(pp_A_tt_TH);
}



log10_bbF_A_bb_TH::log10_bbF_A_bb_TH(const StandardModel& SM_i)
: CPoddHiggs(SM_i)
{}

double log10_bbF_A_bb_TH::computeThValue()
{
    computeParameters();
    return log10(bbF_A_bb_TH);
}
