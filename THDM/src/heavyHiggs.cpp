/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "heavyHiggs.h"
#include "lightHiggs.h"
#include "StandardModel.h"
//#include "gslpp.h"

heavyHiggsCache::heavyHiggsCache(const StandardModel& SM_i): 

        ThObservable(SM_i), 
        myTHDM(static_cast<const THDM*> (&SM_i)),
        mySM (SM_i)
        
{
    mycache = new THDMcache();
    mylightHiggs = new lightHiggs(SM_i);
}


heavyHiggsCache::~heavyHiggsCache()
{}


void heavyHiggsCache::computeParameters()
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
    double cos_2b=cosb*cosb-sinb*sinb;
    double cos_ab=cosa*cosb+sina*sinb;

    //These cross sections are necessary for rHH_gg
    //SM gg -> h production cross section at 8 TeV
    SigmaggF = mycache->cs_ggFtoHP(mHh);
    //SM gg -> h production cross section at 8 TeV, top loop only
    Sigmaggh_tt = mycache->cs_ggHP_tt(mHh);
    //SM gg -> h production cross section at 8 TeV, bottom loop only
    Sigmaggh_bb = mycache->cs_ggHP_bb(mHh);

    /* r_ii is the ratio between the squared 2HDM vertex coupling of the heavy Higgs to
     * the particle i and the corresponding coupling of the SM Higgs boson.*/
    rHH_QuQu=sina/sinb*sina/sinb;
    rHH_QdQd=0.0;//It depends on the modelType
    rHH_ll=0.0;//It depends on the modelType
    rHH_gg=0.0;//It depends on the modelType
    rHH_VV=cos_ab*cos_ab;

    //The ga-ga decay width of the heavy Higgs boson
    Gamma_Hgaga=0.0;//It depends on the modelType
    //The Z-ga decay width of the heavy Higgs boson
    Gamma_HZga=0.0;//It depends on the modelType
    //The g-g decay width of the heavy Higgs boson
    Gamma_Hgg=0.0;//It depends on the modelType

    /*Calulation of rHH_QdQd, rHH_ll, rHH_gg, Gamma_Hgaga, Gamma_HZga, Gamma_Hgg
     * (they depend on the model type): START*/

    /*Gamma_Hgaga and Gamma_HZga expressions can be found in
     "The Higgs Hunter's Guide", Appendix C and in arXiv:0902.4665v3, Appendix A;
     *Gamma_Hgg expression can be found in arXiv:0902.4665v3, Appendix A*/

    /*TAUi is necessary for Gamma_Hgaga, Gamma_HZga and Gamma_Hgg*/
    TAUc=4.0*Mc*Mc/(mHh*mHh);
    TAUt=4.0*Mt*Mt/(mHh*mHh);
    TAUs=4.0*Ms*Ms/(mHh*mHh);
    TAUb=4.0*Mb*Mb/(mHh*mHh);
    TAUmu=4.0*Mmu*Mmu/(mHh*mHh);
    TAUtau=4.0*Mtau*Mtau/(mHh*mHh);
    TAUw=4.0*MW*MW/(mHh*mHh);
    TAUhp=4.0*mHp*mHp/(mHh*mHh);

    /*I_HH_f, I_HH_W and I_HH_Hp are necessary for Gamma_Hgaga, Gamma_HZga and Gamma_Hgg;
     * their expressions can be found in "The Higgs Hunter's Guide", Appendix C, C.4*/
    I_HH_f=0.0;//It depends on the modelType
    I_HH_fU=-(8./3.)*(TAUc*(1+(1-TAUc)*myfunctions.f_func(TAUc))
                      +TAUt*(1+(1-TAUt)*myfunctions.f_func(TAUt)));
    I_HH_fD=-(2./3.)*(TAUs*(1+(1-TAUs)*myfunctions.f_func(TAUs))
                      +TAUb*(1+(1-TAUb)*myfunctions.f_func(TAUb)));
    I_HH_fL=-2.*(TAUmu*(1+(1-TAUmu)*myfunctions.f_func(TAUmu))
                 +TAUtau*(1+(1-TAUtau)*myfunctions.f_func(TAUtau)));
    I_HH_W=cos_ab*(2.0 + 3.0*TAUw + 3.0*TAUw*(2.0-TAUw)
                          *myfunctions.f_func(TAUw));
    /* g_HH_HpHm is the coupling of the heavy Higgs boson to Hp and Hm; its
     * expression can be found in arXiv:1403.1264v2, formula 5*/
    g_HH_HpHm = (cos_ab*(mHh*mHh-2.0*mHp*mHp)
                        -(cosa/cosb+sina/sinb)*(mHh*mHh
                        -m12_2/(cosb*sinb)))/vev;
    I_HH_Hp=-TAUhp*(1.0-TAUhp*myfunctions.f_func(TAUhp))*vev
                            /(2.*mHp*mHp)*g_HH_HpHm;

    //LAMi is necessary for Gamma_HZga.
    LAMc=4.0*Mc*Mc/(MZ*MZ);
    LAMt=4.0*Mt*Mt/(MZ*MZ);
    LAMs=4.0*Ms*Ms/(MZ*MZ);
    LAMb=4.0*Mb*Mb/(MZ*MZ);
    LAMmu=4.0*Mmu*Mmu/(MZ*MZ);
    LAMtau=4.0*Mtau*Mtau/(MZ*MZ);
    LAMw=4.0*MW*MW/(MZ*MZ);
    LAMhp=4.0*mHp*mHp/(MZ*MZ);

    /*A_HH_F, A_HH_W and A_HH_Hp are necessary for Gamma_HZga*/
    /*A_HH_F expression can be found in "The Higgs Hunter's Guide", Appendix C, C.12*/
    A_HH_F = 0.0;//It depends on the modelType
    A_HH_U = -4.0*(1.0/2.0-4.0/3.0*s02)*(myfunctions.Int1(TAUc,LAMc)-myfunctions.Int2(TAUc,LAMc)
                                         +myfunctions.Int1(TAUt,LAMt)-myfunctions.Int2(TAUt,LAMt));
    A_HH_D = +2.0*(-1.0/2.0+2.0/3.0*s02)*(myfunctions.Int1(TAUs,LAMs)-myfunctions.Int2(TAUs,LAMs)
                                          +myfunctions.Int1(TAUb,LAMb)-myfunctions.Int2(TAUb,LAMb));
    A_HH_L = +2.0*(-1.0/2.0+2.0*s02)*(myfunctions.Int1(TAUmu,LAMmu)-myfunctions.Int2(TAUmu,LAMmu)
                                      +myfunctions.Int1(TAUtau,LAMtau)-myfunctions.Int2(TAUtau,LAMtau));
    /*A_HH_W expression can be found in "The Higgs Hunter's Guide", Appendix C, C.13*/
    A_HH_W = -cos_ab*sqrt(c02/s02)*(4*(3-s02/c02)*myfunctions.Int2(TAUw,LAMw)
                            +((1.0+2.0/TAUw)*s02/c02-(5.0+2.0/TAUw))*myfunctions.Int1(TAUw,LAMw));
    /*A_HH_Hp expression can be found in "The Higgs Hunter's Guide", Appendix C, C.14*/
    A_HH_Hp= g_HH_HpHm*(1-2.0*s02)/sqrt(c02*s02)*myfunctions.Int1(TAUhp,LAMhp)
                            *vev/(2.*mHp*mHp);
    std::string modelflag=myTHDM->getModelTypeflag();

    if( modelflag == "type1" ) {
        rHH_gg=sina/sinb*sina/sinb;
        rHH_QdQd=sina/sinb*sina/sinb;
        rHH_ll=sina/sinb*sina/sinb;
        I_HH_f=sina/sinb*(I_HH_fU+I_HH_fD+I_HH_fL);
        A_HH_F = sina/sinb*(A_HH_U+A_HH_D+A_HH_L)/sqrt(s02*c02);
    }
    else if( modelflag == "type2" ) {
        rHH_gg=sina/sinb*cosa/cosb+(Sigmaggh_tt*sina/sinb*(sina/sinb-cosa/cosb)
             +Sigmaggh_bb*cosa/cosb*(cosa/cosb-sina/sinb))/SigmaggF;
        rHH_QdQd=cosa/cosb*cosa/cosb;
        rHH_ll=cosa/cosb*cosa/cosb;
        I_HH_f=sina/sinb*I_HH_fU+cosa/cosb*(I_HH_fD+I_HH_fL);
        A_HH_F = (sina/sinb*A_HH_U+cosa/cosb*(A_HH_D+A_HH_L))/sqrt(s02*c02);
    }
    else if( modelflag == "typeX" ) {
        rHH_gg=sina/sinb*sina/sinb;
        rHH_QdQd=sina/sinb*sina/sinb;
        rHH_ll=cosa/cosb*cosa/cosb;
        I_HH_f=sina/sinb*(I_HH_fU+I_HH_fD)+cosa/cosb*I_HH_fL;
        A_HH_F = (sina/sinb*(A_HH_U+A_HH_D)+cosa/cosb*A_HH_L)/sqrt(s02*c02);
    }
    else if( modelflag == "typeY" ) {
        rHH_gg=sina/sinb*cosa/cosb+(Sigmaggh_tt*sina/sinb*(sina/sinb-cosa/cosb)
             +Sigmaggh_bb*cosa/cosb*(cosa/cosb-sina/sinb))/SigmaggF;
        rHH_QdQd=cosa/cosb*cosa/cosb;
        rHH_ll=sina/sinb*sina/sinb;
        I_HH_f=sina/sinb*(I_HH_fU+I_HH_fL)+cosa/cosb*I_HH_fD;
        A_HH_F = (sina/sinb*(A_HH_U+A_HH_L)+cosa/cosb*A_HH_D)/sqrt(s02*c02);
    }
    else {
        throw std::runtime_error("modelflag can be only any of \"type1\", \"type2\", \"typeX\" or \"typeY\"");
    }

    /*Gamma_Hgaga expression can be found in in arXiv:0902.4665v3, Appendix A, A.8*/
    Gamma_Hgaga=GF*Ale*Ale*mHh*mHh*mHh/(sqrt(2)*128.0*M_PI*M_PI*M_PI)
                *(I_HH_f+I_HH_W+I_HH_Hp).abs()*(I_HH_f+I_HH_W+I_HH_Hp).abs();

//    std::cout<<"rHH_gaga: "<<(I_HH_f+I_HH_W+I_HH_Hp).abs()/(I_HH_f+I_HH_W).abs()<<std::endl;

    /*Gamma_HZga expression can be found in in arXiv:0902.4665v3, Appendix A, A.9*/
    Gamma_HZga=GF*Ale*Ale*mHh*mHh*mHh/(sqrt(2)*64.0*M_PI*M_PI*M_PI)
               *(1.0-MZ*MZ/(mHh*mHh))*(1.0-MZ*MZ/(mHh*mHh))*(1.0-MZ*MZ/(mHh*mHh))
               *(A_HH_F+A_HH_W+A_HH_Hp).abs()*(A_HH_F+A_HH_W+A_HH_Hp).abs();
    /*Gamma_Hgg expression can be found in in arXiv:0902.4665v3, Appendix A, A.10; relative coupling see above*/
    Gamma_Hgg=GF*Als*Als*mHh*mHh*mHh/(sqrt(2)*64.0*M_PI*M_PI*M_PI)*rHH_gg;

    /*Calulation of rHH_QdQd, rHH_ll, rHH_gg, Gamma_Hgaga, Gamma_HZga, Gamma_Hgg
     * (they depend on the model type): FINISH*/

    SigmaTotSM_H=mycache->cs_ggFtoHP(mHh)/mycache->pc_ggFtoHP(mHh);
    SigmaggF_H=mycache->cs_ggFtoHP(mHh)*rHH_gg;
    SigmabbF_H=mycache->cs_bbFtoHP(mHh)*rHH_QdQd;
    SigmaVBF_H=mycache->pc_VBFtoHP(mHh)*SigmaTotSM_H*rHH_VV;
    SigmattF_H=mycache->pc_ttFtoHP(mHh)*SigmaTotSM_H*rHH_QuQu;
    SigmaVH_H=(mycache->pc_WHP_HP(mHh)+mycache->pc_ZHP_HP(mHh))*SigmaTotSM_H*rHH_VV;

    SigmaSum = SigmaggF_H + SigmaVBF_H + SigmaVH_H + SigmattF_H + SigmabbF_H;

    BrSM_Htott=mycache->Br_HPtott(mHh);
    BrSM_Htocc=mycache->Br_HPtocc(mHh);
    BrSM_Htobb=mycache->Br_HPtobb(mHh);
    BrSM_Htotautau=mycache->Br_HPtotautau(mHh);
    BrSM_Htomumu=mycache->Br_HPtomumu(mHh);
    BrSM_HtoWW =mycache->Br_HPtoWW(mHh);
    BrSM_HtoZZ =mycache->Br_HPtoZZ(mHh);

    GammaHtotSM=mycache->GammaHPtotSM(mHh);

    GammaHhh=myfunctions.HSTheta(mHh - 2.0*mHl)*sqrt(std::abs(1.0 - (4.0*mHl*mHl)/(mHh*mHh)))
                    *std::abs((cos_ab*cos_ab/(4.0*sinb*cosb*sinb*cosb)
                    *pow(m12_2 + mHh*mHh*cosa*sina + (2.0* //Otto added a factor of 2 - check this!
            mHl*mHl - 3.0*m12_2/(sinb*cosb))
                    *sina*cosa,2))/(vev*vev))/(8.0*mHh*M_PI);

    GammaHHpHm=myfunctions.HSTheta(mHh - 2.0*mHp)*sqrt(std::abs(1.0 - (4.0*mHp*mHp)/(mHh*mHh)))
                      *std::abs(pow(cos_ab*(mHh*mHh + 2.0*mHp*mHp - 2.0*m12_2/sinb/cosb)
                      - cos_2b/(sinb*cosb)*(mHh*mHh - m12_2/sinb/cosb)*sin_ba,2)/(vev*vev))
                      /(16.0*mHh*M_PI);

    GammaHAA=myfunctions.HSTheta(-2.0*mA + mHh)*sqrt(std::abs(1.0 - (4.0*mA*mA)/(mHh*mHh)))
                    *std::abs(pow(cos_ab*(2*mA*mA + mHh*mHh - 2.0*m12_2/sinb/cosb)
                    - cos_2b/(sinb*cosb)*(mHh*mHh - m12_2/sinb/cosb)*sin_ba,2)/(vev*vev))
                    /(32.0*mHh*M_PI);

    GammaHAZ=myfunctions.HSTheta(mHh-mA-MZ)*pow(myfunctions.KaellenFunction(mHh,MZ,mA),3)
                    *sin_ba*sin_ba/(2.0*M_PI*vev*vev);

    GammaHHpW=myfunctions.HSTheta(mHh-mHp-MW)*pow(myfunctions.KaellenFunction(mHh,MW,mHp),3)*sin_ba*sin_ba/(M_PI*vev*vev);

    GammaHtot= ((BrSM_Htott+BrSM_Htocc)*rHH_QuQu
                    +BrSM_Htobb*rHH_QdQd
                    +(BrSM_Htotautau+BrSM_Htomumu)*rHH_ll
                    +(BrSM_HtoWW+BrSM_HtoZZ)*rHH_VV)*GammaHtotSM
               +GammaHhh+GammaHHpHm+GammaHAA+GammaHAZ+GammaHHpW+Gamma_Hgaga
               +Gamma_HZga+Gamma_Hgg;

    Br_Htott=BrSM_Htott*rHH_QuQu*GammaHtotSM/GammaHtot;
    Br_Htobb=BrSM_Htobb*rHH_QdQd*GammaHtotSM/GammaHtot;
    Br_Htotautau=BrSM_Htotautau*rHH_ll*GammaHtotSM/GammaHtot;
    Br_HtoWW=BrSM_HtoWW*rHH_VV*GammaHtotSM/GammaHtot;
    Br_HtoZZ=BrSM_HtoZZ*rHH_VV*GammaHtotSM/GammaHtot;
    Br_Htohh=GammaHhh/GammaHtot;
    Br_Htogaga=Gamma_Hgaga/GammaHtot;

    Br_htobb=mylightHiggs->THDM_BR_h_bb();
    Br_htogaga=mylightHiggs->THDM_BR_h_gaga();
    Br_htotautau=mylightHiggs->THDM_BR_h_tautau();

    //Heavy Higgs Signals, theoretical expressions

    ggF_H_tautau_TH=SigmaggF_H*Br_Htotautau;
    bbF_H_tautau_TH=SigmabbF_H*Br_Htotautau;
    ggF_H_gaga_TH=SigmaggF_H*Br_Htogaga;
    pp_H_ZZ_TH=SigmaSum/SigmaTotSM_H*rHH_VV*GammaHtotSM/GammaHtot;
    ggF_H_WW_TH=SigmaggF_H*Br_HtoWW;
    VBF_H_WW_TH=SigmaVBF_H*Br_HtoWW;
    ggF_H_hh_TH=SigmaggF_H*Br_Htohh;
    ggF_H_hh_bbtautau_TH=SigmaggF_H*Br_Htohh*Br_htobb*Br_htotautau;    
    pp_H_hh_bbbb_TH=SigmaSum*Br_Htohh*Br_htobb*Br_htobb;
    pp_H_hh_gagabb_TH=SigmaSum*Br_Htohh*Br_htobb*Br_htogaga;
    pp_H_tt_TH=SigmaSum*Br_Htott;
    bbF_H_bb_TH=SigmabbF_H*Br_Htobb; 

    //Heavy Higgs Signals, experimental expressions

    ggF_H_tautau_EX_ATLAS=mycache->ex_ggF_phi_tautau_ATLAS(mHh);
    ggF_H_tautau_EX_CMS=mycache->ex_ggF_phi_tautau_CMS(mHh);    
    bbF_H_tautau_EX_ATLAS=mycache->ex_bbF_phi_tautau_ATLAS(mHh);  
    bbF_H_tautau_EX_CMS=mycache->ex_bbF_phi_tautau_CMS(mHh);       
    ggF_H_gaga_EX_ATLAS=mycache->ex_ggF_phi_gaga_ATLAS(mHh); 
    ggF_H_gaga_EX_CMS=mycache->ex_ggF_phi_gaga_CMS(mHh);
    pp_H_ZZ_EX_CMS=mycache->ex_pp_H_ZZ_CMS(mHh);    
    ggF_H_WW_EX_ATLAS=mycache->ex_ggF_H_WW_ATLAS(mHh);   
    VBF_H_WW_EX_ATLAS=mycache->ex_VBF_H_WW_ATLAS(mHh);           
    ggF_H_hh_EX_ATLAS=mycache->ex_ggF_H_hh_ATLAS(mHh);  
    ggF_H_hh_bbtautau_EX_CMS=mycache->ex_ggF_H_hh_bbtautau_CMS(mHh);     
    pp_H_hh_bbbb_EX_CMS=mycache->ex_pp_phi_hh_bbbb_CMS(mHh);   
    pp_H_hh_gagabb_EX_CMS=mycache->ex_pp_phi_hh_gagabb_CMS(mHh);    
    pp_H_tt_EX_ATLAS=mycache->ex_pp_phi_tt_ATLAS(mHh);    
    bbF_H_bb_EX_CMS=mycache->ex_bbF_phi_bb_CMS(mHh);    
    
}

double heavyHiggsCache::computeThValue()
{
    return 0;
}

/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/

Hobs_ggF_H_tautau_ATLAS::Hobs_ggF_H_tautau_ATLAS(const StandardModel& SM_i)
: heavyHiggsCache(SM_i)
{}

double Hobs_ggF_H_tautau_ATLAS::computeThValue()
{
    computeParameters();
    return ggF_H_tautau_TH/ggF_H_tautau_EX_ATLAS;
}



Hobs_ggF_H_tautau_CMS::Hobs_ggF_H_tautau_CMS(const StandardModel& SM_i)
: heavyHiggsCache(SM_i)
{}

double Hobs_ggF_H_tautau_CMS::computeThValue()
{
    computeParameters();
    return ggF_H_tautau_TH/ggF_H_tautau_EX_CMS;
}



Hobs_bbF_H_tautau_ATLAS::Hobs_bbF_H_tautau_ATLAS(const StandardModel& SM_i)
: heavyHiggsCache(SM_i)
{}

double Hobs_bbF_H_tautau_ATLAS::computeThValue()
{
    computeParameters();
    return bbF_H_tautau_TH/bbF_H_tautau_EX_ATLAS;
}



Hobs_bbF_H_tautau_CMS::Hobs_bbF_H_tautau_CMS(const StandardModel& SM_i)
: heavyHiggsCache(SM_i)
{}

double Hobs_bbF_H_tautau_CMS::computeThValue()
{
    computeParameters();
    return bbF_H_tautau_TH/bbF_H_tautau_EX_CMS;
}



Hobs_ggF_H_gaga_ATLAS::Hobs_ggF_H_gaga_ATLAS(const StandardModel& SM_i)
: heavyHiggsCache(SM_i)
{}

double Hobs_ggF_H_gaga_ATLAS::computeThValue()
{
    computeParameters();
    return ggF_H_gaga_TH/ggF_H_gaga_EX_ATLAS;
}



Hobs_ggF_H_gaga_CMS::Hobs_ggF_H_gaga_CMS(const StandardModel& SM_i)
: heavyHiggsCache(SM_i)
{}

double Hobs_ggF_H_gaga_CMS::computeThValue()
{
    computeParameters();
    return ggF_H_gaga_TH/ggF_H_gaga_EX_CMS;
}



Hobs_pp_H_ZZ_CMS::Hobs_pp_H_ZZ_CMS(const StandardModel& SM_i)
: heavyHiggsCache(SM_i)
{}

double Hobs_pp_H_ZZ_CMS::computeThValue()
{
    computeParameters();
    return pp_H_ZZ_TH/pp_H_ZZ_EX_CMS;
}



Hobs_ggF_H_WW_ATLAS::Hobs_ggF_H_WW_ATLAS(const StandardModel& SM_i)
: heavyHiggsCache(SM_i)
{}

double Hobs_ggF_H_WW_ATLAS::computeThValue()
{
    computeParameters();
    return ggF_H_WW_TH/ggF_H_WW_EX_ATLAS;
}



Hobs_VBF_H_WW_ATLAS::Hobs_VBF_H_WW_ATLAS(const StandardModel& SM_i)
: heavyHiggsCache(SM_i)
{}

double Hobs_VBF_H_WW_ATLAS::computeThValue()
{
    computeParameters();
    return VBF_H_WW_TH/VBF_H_WW_EX_ATLAS;
}



Hobs_ggF_H_hh_ATLAS::Hobs_ggF_H_hh_ATLAS(const StandardModel& SM_i)
: heavyHiggsCache(SM_i)
{}

double Hobs_ggF_H_hh_ATLAS::computeThValue()
{
    computeParameters();
    return ggF_H_hh_TH/ggF_H_hh_EX_ATLAS;
}



Hobs_ggF_H_hh_bbtautau_CMS::Hobs_ggF_H_hh_bbtautau_CMS(const StandardModel& SM_i)
: heavyHiggsCache(SM_i)
{}

double Hobs_ggF_H_hh_bbtautau_CMS::computeThValue()
{
    computeParameters();
    return ggF_H_hh_bbtautau_TH/ggF_H_hh_bbtautau_EX_CMS;
}



Hobs_pp_H_hh_bbbb_CMS::Hobs_pp_H_hh_bbbb_CMS(const StandardModel& SM_i)
: heavyHiggsCache(SM_i)
{}

double Hobs_pp_H_hh_bbbb_CMS::computeThValue()
{
    computeParameters();
    return pp_H_hh_bbbb_TH/pp_H_hh_bbbb_EX_CMS;
}



Hobs_pp_H_hh_gagabb_CMS::Hobs_pp_H_hh_gagabb_CMS(const StandardModel& SM_i)
: heavyHiggsCache(SM_i)
{}

double Hobs_pp_H_hh_gagabb_CMS::computeThValue()
{
    computeParameters();
    return pp_H_hh_gagabb_TH/pp_H_hh_gagabb_EX_CMS;
}



Hobs_pp_H_tt_ATLAS::Hobs_pp_H_tt_ATLAS(const StandardModel& SM_i)
: heavyHiggsCache(SM_i)
{}

double Hobs_pp_H_tt_ATLAS::computeThValue()
{
    computeParameters();
    return pp_H_tt_TH/pp_H_tt_EX_ATLAS;
}



Hobs_bbF_H_bb_CMS::Hobs_bbF_H_bb_CMS(const StandardModel& SM_i)
: heavyHiggsCache(SM_i)
{}

double Hobs_bbF_H_bb_CMS::computeThValue()
{
    computeParameters();
    return bbF_H_bb_TH/bbF_H_bb_EX_CMS;
}



log10_ggF_H_tautau_TH::log10_ggF_H_tautau_TH(const StandardModel& SM_i)
: heavyHiggsCache(SM_i)
{}

double log10_ggF_H_tautau_TH::computeThValue()
{
    computeParameters();
    return log10(ggF_H_tautau_TH);
}



log10_bbF_H_tautau_TH::log10_bbF_H_tautau_TH(const StandardModel& SM_i)
: heavyHiggsCache(SM_i)
{}

double log10_bbF_H_tautau_TH::computeThValue()
{
    computeParameters();
    return log10(bbF_H_tautau_TH);
}



log10_ggF_H_gaga_TH::log10_ggF_H_gaga_TH(const StandardModel& SM_i)
: heavyHiggsCache(SM_i)
{}

double log10_ggF_H_gaga_TH::computeThValue()
{
    computeParameters();
    return log10(ggF_H_gaga_TH);
}



log10_pp_H_ZZ_TH::log10_pp_H_ZZ_TH(const StandardModel& SM_i)
: heavyHiggsCache(SM_i)
{}

double log10_pp_H_ZZ_TH::computeThValue()
{
    computeParameters();
    return log10(pp_H_ZZ_TH);
}



log10_ggF_H_WW_TH::log10_ggF_H_WW_TH(const StandardModel& SM_i)
: heavyHiggsCache(SM_i)
{}

double log10_ggF_H_WW_TH::computeThValue()
{
    computeParameters();
    return log10(ggF_H_WW_TH);
}



log10_VBF_H_WW_TH::log10_VBF_H_WW_TH(const StandardModel& SM_i)
: heavyHiggsCache(SM_i)
{}

double log10_VBF_H_WW_TH::computeThValue()
{
    computeParameters();
    return log10(VBF_H_WW_TH);
}



log10_ggF_H_hh_TH::log10_ggF_H_hh_TH(const StandardModel& SM_i)
: heavyHiggsCache(SM_i)
{}

double log10_ggF_H_hh_TH::computeThValue()
{
    computeParameters();
    return log10(ggF_H_hh_TH);
}



log10_ggF_H_hh_bbtautau_TH::log10_ggF_H_hh_bbtautau_TH(const StandardModel& SM_i)
  : heavyHiggsCache(SM_i)
{}

double log10_ggF_H_hh_bbtautau_TH::computeThValue()
{
  computeParameters();
  return log10(ggF_H_hh_bbtautau_TH);
}



log10_pp_H_hh_bbbb_TH::log10_pp_H_hh_bbbb_TH(const StandardModel& SM_i)
: heavyHiggsCache(SM_i)
{}

double log10_pp_H_hh_bbbb_TH::computeThValue()
{
    computeParameters();
    return log10(pp_H_hh_bbbb_TH);
}



log10_pp_H_hh_gagabb_TH::log10_pp_H_hh_gagabb_TH(const StandardModel& SM_i)
: heavyHiggsCache(SM_i)
{}

double log10_pp_H_hh_gagabb_TH::computeThValue()
{
    computeParameters();
    return log10(pp_H_hh_gagabb_TH);
}



log10_pp_H_tt_TH::log10_pp_H_tt_TH(const StandardModel& SM_i)
: heavyHiggsCache(SM_i)
{}

double log10_pp_H_tt_TH::computeThValue()
{
    computeParameters();
    return log10(pp_H_tt_TH);
}



log10_bbF_H_bb_TH::log10_bbF_H_bb_TH(const StandardModel& SM_i)
: heavyHiggsCache(SM_i)
{}

double log10_bbF_H_bb_TH::computeThValue()
{
    computeParameters();
    return log10(bbF_H_bb_TH);
}
