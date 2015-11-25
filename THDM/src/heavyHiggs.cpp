/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "heavyHiggs.h"
#include "lightHiggs.h"
#include "StandardModel.h"

heavyHiggs::heavyHiggs(const StandardModel& SM_i): 

        ThObservable(SM_i), 
        myTHDM(static_cast<const THDM*> (&SM_i)),
        mySM (SM_i)
        
{
    mycache = new THDMcache();
    mylightHiggs = new lightHiggs(SM_i);
}


heavyHiggs::~heavyHiggs()
{}


void heavyHiggs::computeParameters()
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
    double SigmaggF = mycache->ip_cs_ggFtoHP(mHh);
    //SM gg -> h production cross section at 8 TeV, top loop only
    double Sigmaggh_tt = mycache->ip_cs_ggHP_tt(mHh);
    //SM gg -> h production cross section at 8 TeV, bottom loop only
    double Sigmaggh_bb = mycache->ip_cs_ggHP_bb(mHh);

    /* r_ii is the ratio between the squared 2HDM vertex coupling of the heavy Higgs to
     * the particle i and the corresponding coupling of the SM Higgs boson.*/
    double rHH_QuQu=sina/sinb*sina/sinb;
    double rHH_QdQd=0.0;//It depends on the modelType
    double rHH_ll=0.0;//It depends on the modelType
    double rHH_gg=0.0;//It depends on the modelType
    double rHH_VV=cos_ab*cos_ab;

    /*Calulation of rHH_QdQd, rHH_ll, rHH_gg, Gamma_Hgaga, Gamma_HZga, Gamma_Hgg
     * (they depend on the model type): START*/

    /*Gamma_Hgaga and Gamma_HZga expressions can be found in
     "The Higgs Hunter's Guide", Appendix C and in arXiv:0902.4665v3, Appendix A;
     *Gamma_Hgg expression can be found in arXiv:0902.4665v3, Appendix A*/

    double TAUc=4.0*Mc*Mc/(mHh*mHh);
    double TAUt=4.0*Mt*Mt/(mHh*mHh);
    double TAUs=4.0*Ms*Ms/(mHh*mHh);
    double TAUb=4.0*Mb*Mb/(mHh*mHh);
    double TAUmu=4.0*Mmu*Mmu/(mHh*mHh);
    double TAUtau=4.0*Mtau*Mtau/(mHh*mHh);
    double TAUw=4.0*MW*MW/(mHh*mHh);
    double TAUhp=4.0*mHp*mHp/(mHh*mHh);

    /*I_HH_F, I_HH_W and I_HH_Hp are needed for Gamma_Hgaga;
     * their expressions can be found in "The Higgs Hunter's Guide", Appendix C, C.4*/
    gslpp::complex I_HH_F=0.0;//It depends on the modelType
    gslpp::complex I_HH_U=-(8./3.)*(TAUc*(1+(1-TAUc)*myfunctions.f_func(TAUc))
                      +TAUt*(1+(1-TAUt)*myfunctions.f_func(TAUt)));
    gslpp::complex I_HH_D=-(2./3.)*(TAUs*(1+(1-TAUs)*myfunctions.f_func(TAUs))
                      +TAUb*(1+(1-TAUb)*myfunctions.f_func(TAUb)));
    gslpp::complex I_HH_L=-2.*(TAUmu*(1+(1-TAUmu)*myfunctions.f_func(TAUmu))
                 +TAUtau*(1+(1-TAUtau)*myfunctions.f_func(TAUtau)));
    gslpp::complex I_HH_W=cos_ab*(2.0 + 3.0*TAUw + 3.0*TAUw*(2.0-TAUw)
                          *myfunctions.f_func(TAUw));
    /* g_HH_HpHm is the coupling of the heavy Higgs boson to Hp and Hm; its
     * expression can be found in arXiv:1403.1264v2, formula 5*/
    double g_HH_HpHm = (cos_ab*(mHh*mHh-2.0*mHp*mHp)
                        -(cosa/cosb+sina/sinb)*(mHh*mHh
                        -m12_2/(cosb*sinb)))/vev;
    gslpp::complex I_HH_Hp=-TAUhp*(1.0-TAUhp*myfunctions.f_func(TAUhp))*vev
                            /(2.*mHp*mHp)*g_HH_HpHm;

    double LAMc=4.0*Mc*Mc/(MZ*MZ);
    double LAMt=4.0*Mt*Mt/(MZ*MZ);
    double LAMs=4.0*Ms*Ms/(MZ*MZ);
    double LAMb=4.0*Mb*Mb/(MZ*MZ);
    double LAMmu=4.0*Mmu*Mmu/(MZ*MZ);
    double LAMtau=4.0*Mtau*Mtau/(MZ*MZ);
    double LAMw=4.0*MW*MW/(MZ*MZ);
    double LAMhp=4.0*mHp*mHp/(MZ*MZ);

    /*A_HH_F, A_HH_W and A_HH_Hp are needed for Gamma_HZga;
     * their expressions can be found in "The Higgs Hunter's Guide", Appendix C, C.12*/
    gslpp::complex A_HH_F = 0.0;//It depends on the modelType
    gslpp::complex A_HH_U = -4.0*(1.0/2.0-4.0/3.0*sW2)*(myfunctions.Int1(TAUc,LAMc)-myfunctions.Int2(TAUc,LAMc)
                                         +myfunctions.Int1(TAUt,LAMt)-myfunctions.Int2(TAUt,LAMt));
    gslpp::complex A_HH_D = +2.0*(-1.0/2.0+2.0/3.0*sW2)*(myfunctions.Int1(TAUs,LAMs)-myfunctions.Int2(TAUs,LAMs)
                                          +myfunctions.Int1(TAUb,LAMb)-myfunctions.Int2(TAUb,LAMb));
    gslpp::complex A_HH_L = +2.0*(-1.0/2.0+2.0*sW2)*(myfunctions.Int1(TAUmu,LAMmu)-myfunctions.Int2(TAUmu,LAMmu)
                                      +myfunctions.Int1(TAUtau,LAMtau)-myfunctions.Int2(TAUtau,LAMtau));
    /*A_HH_W expression can be found in "The Higgs Hunter's Guide", Appendix C, C.13*/
    gslpp::complex A_HH_W = -cos_ab*sqrt(cW2/sW2)*(4*(3-sW2/cW2)*myfunctions.Int2(TAUw,LAMw)
                            +((1.0+2.0/TAUw)*sW2/cW2-(5.0+2.0/TAUw))*myfunctions.Int1(TAUw,LAMw));
    /*A_HH_Hp expression can be found in "The Higgs Hunter's Guide", Appendix C, C.14*/
    gslpp::complex A_HH_Hp= g_HH_HpHm*(1-2.0*sW2)/sqrt(cW2*sW2)*myfunctions.Int1(TAUhp,LAMhp)
                            *vev/(2.*mHp*mHp);

    std::string modelflag=myTHDM->getModelTypeflag();

    if( modelflag == "type1" ) {
        rHH_gg=sina/sinb*sina/sinb;
        rHH_QdQd=sina/sinb*sina/sinb;
        rHH_ll=sina/sinb*sina/sinb;
        I_HH_F=sina/sinb*(I_HH_U+I_HH_D+I_HH_L);
        A_HH_F = sina/sinb*(A_HH_U+A_HH_D+A_HH_L)/sqrt(sW2*cW2);
    }
    else if( modelflag == "type2" ) {
        rHH_gg=sina/sinb*cosa/cosb+(Sigmaggh_tt*sina/sinb*(sina/sinb-cosa/cosb)
             +Sigmaggh_bb*cosa/cosb*(cosa/cosb-sina/sinb))/SigmaggF;
        rHH_QdQd=cosa/cosb*cosa/cosb;
        rHH_ll=cosa/cosb*cosa/cosb;
        I_HH_F=sina/sinb*I_HH_U+cosa/cosb*(I_HH_D+I_HH_L);
        A_HH_F = (sina/sinb*A_HH_U+cosa/cosb*(A_HH_D+A_HH_L))/sqrt(sW2*cW2);
    }
    else if( modelflag == "typeX" ) {
        rHH_gg=sina/sinb*sina/sinb;
        rHH_QdQd=sina/sinb*sina/sinb;
        rHH_ll=cosa/cosb*cosa/cosb;
        I_HH_F=sina/sinb*(I_HH_U+I_HH_D)+cosa/cosb*I_HH_L;
        A_HH_F = (sina/sinb*(A_HH_U+A_HH_D)+cosa/cosb*A_HH_L)/sqrt(sW2*cW2);
    }
    else if( modelflag == "typeY" ) {
        rHH_gg=sina/sinb*cosa/cosb+(Sigmaggh_tt*sina/sinb*(sina/sinb-cosa/cosb)
             +Sigmaggh_bb*cosa/cosb*(cosa/cosb-sina/sinb))/SigmaggF;
        rHH_QdQd=cosa/cosb*cosa/cosb;
        rHH_ll=sina/sinb*sina/sinb;
        I_HH_F=sina/sinb*(I_HH_U+I_HH_L)+cosa/cosb*I_HH_D;
        A_HH_F = (sina/sinb*(A_HH_U+A_HH_L)+cosa/cosb*A_HH_D)/sqrt(sW2*cW2);
    }
    else {
        throw std::runtime_error("modelflag can be only any of \"type1\", \"type2\", \"typeX\" or \"typeY\"");
    }

    /*Gamma_Hgaga expression can be found in arXiv:0902.4665v3, Appendix A, A.8*/
    double Gamma_Hgaga=GF*Ale*Ale*mHh*mHh*mHh/(sqrt(2)*128.0*M_PI*M_PI*M_PI)
                *(I_HH_F+I_HH_W+I_HH_Hp).abs()*(I_HH_F+I_HH_W+I_HH_Hp).abs();

    /*Gamma_HZga expression can be found in arXiv:0902.4665v3, Appendix A, A.9*/
    double Gamma_HZga=GF*Ale*Ale*mHh*mHh*mHh/(sqrt(2)*64.0*M_PI*M_PI*M_PI)
               *(1.0-MZ*MZ/(mHh*mHh))*(1.0-MZ*MZ/(mHh*mHh))*(1.0-MZ*MZ/(mHh*mHh))
               *(A_HH_F+A_HH_W+A_HH_Hp).abs()*(A_HH_F+A_HH_W+A_HH_Hp).abs();

    /*Gamma_Hgg expression can be found in arXiv:0902.4665v3, Appendix A, A.10; relative coupling see above*/
    double Gamma_Hgg=GF*Als*Als*mHh*mHh*mHh/(sqrt(2)*64.0*M_PI*M_PI*M_PI)*rHH_gg;

    /*Calulation of rHH_QdQd, rHH_ll, rHH_gg, Gamma_Hgaga, Gamma_HZga, Gamma_Hgg: END*/

    double SigmaTotSM_H=mycache->ip_cs_ggFtoHP(mHh)/mycache->ip_pc_ggFtoHP(mHh);
    double SigmaggF_H=mycache->ip_cs_ggFtoHP(mHh)*rHH_gg;
    double SigmabbF_H=mycache->ip_cs_bbFtoHP(mHh)*rHH_QdQd;
    double SigmaVBF_H=mycache->ip_pc_VBFtoHP(mHh)*SigmaTotSM_H*rHH_VV;
    double SigmattF_H=mycache->ip_pc_ttFtoHP(mHh)*SigmaTotSM_H*rHH_QuQu;
    double SigmaVH_H=(mycache->ip_pc_WHP_HP(mHh)+mycache->ip_pc_ZHP_HP(mHh))*SigmaTotSM_H*rHH_VV;

    double SigmaSum = SigmaggF_H + SigmaVBF_H + SigmaVH_H + SigmattF_H + SigmabbF_H;

    double BrSM_Htott=mycache->ip_Br_HPtott(mHh);
    double BrSM_Htocc=mycache->ip_Br_HPtocc(mHh);
    double BrSM_Htobb=mycache->ip_Br_HPtobb(mHh);
    double BrSM_Htotautau=mycache->ip_Br_HPtotautau(mHh);
    double BrSM_Htomumu=mycache->ip_Br_HPtomumu(mHh);
    double BrSM_HtoWW =mycache->ip_Br_HPtoWW(mHh);
    double BrSM_HtoZZ =mycache->ip_Br_HPtoZZ(mHh);

    double GammaHtotSM=mycache->ip_GammaHPtotSM(mHh);

    double GammaHhh=myfunctions.HSTheta(mHh - 2.0*mHl)*sqrt(std::abs(1.0 - (4.0*mHl*mHl)/(mHh*mHh)))
                    *std::abs((cos_ab*cos_ab/(4.0*sinb*cosb*sinb*cosb)
                    *pow(m12_2 + mHh*mHh*cosa*sina + (2.0* //Otto added a factor of 2 - check this!
            mHl*mHl - 3.0*m12_2/(sinb*cosb))
                    *sina*cosa,2))/(vev*vev))/(8.0*mHh*M_PI);

    double GammaHHpHm=myfunctions.HSTheta(mHh - 2.0*mHp)*sqrt(std::abs(1.0 - (4.0*mHp*mHp)/(mHh*mHh)))
                      *std::abs(pow(cos_ab*(mHh*mHh + 2.0*mHp*mHp - 2.0*m12_2/sinb/cosb)
                                -cos_2b/(sinb*cosb)*(mHh*mHh - m12_2/sinb/cosb)*sin_ba,2)/(vev*vev))
                      /(16.0*mHh*M_PI);

    double GammaHAA=myfunctions.HSTheta(-2.0*mA + mHh)*sqrt(std::abs(1.0 - (4.0*mA*mA)/(mHh*mHh)))
                    *std::abs(pow(cos_ab*(2*mA*mA + mHh*mHh - 2.0*m12_2/sinb/cosb)
                    - cos_2b/(sinb*cosb)*(mHh*mHh - m12_2/sinb/cosb)*sin_ba,2)/(vev*vev))
                    /(32.0*mHh*M_PI);

    double GammaHAZ=myfunctions.HSTheta(mHh-mA-MZ)*pow(myfunctions.KaellenFunction(mHh,MZ,mA),3)
                    *sin_ba*sin_ba/(2.0*M_PI*vev*vev);

    double GammaHHpW=myfunctions.HSTheta(mHh-mHp-MW)*pow(myfunctions.KaellenFunction(mHh,MW,mHp),3)*sin_ba*sin_ba/(M_PI*vev*vev);

    GammaHtot= ((BrSM_Htott+BrSM_Htocc)*rHH_QuQu
                    +BrSM_Htobb*rHH_QdQd
                    +(BrSM_Htotautau+BrSM_Htomumu)*rHH_ll
                    +(BrSM_HtoWW+BrSM_HtoZZ)*rHH_VV)*GammaHtotSM
               +Gamma_Hgaga+Gamma_HZga+Gamma_Hgg
               +GammaHhh+GammaHHpHm+GammaHAA+GammaHAZ+GammaHHpW;

    double Br_Htott=BrSM_Htott*rHH_QuQu*GammaHtotSM/GammaHtot;
    double Br_Htobb=BrSM_Htobb*rHH_QdQd*GammaHtotSM/GammaHtot;
    double Br_Htotautau=BrSM_Htotautau*rHH_ll*GammaHtotSM/GammaHtot;
    double Br_HtoWW=BrSM_HtoWW*rHH_VV*GammaHtotSM/GammaHtot;
//    double Br_HtoZZ=BrSM_HtoZZ*rHH_VV*GammaHtotSM/GammaHtot;
    double Br_Htohh=GammaHhh/GammaHtot;
    double Br_Htogaga=Gamma_Hgaga/GammaHtot;

    double Br_htobb=mylightHiggs->THDM_BR_h_bb();
    double Br_htogaga=mylightHiggs->THDM_BR_h_gaga();
    double Br_htotautau=mylightHiggs->THDM_BR_h_tautau();

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
    pp_H_hh_gagabb_TH=SigmaSum*Br_Htohh*Br_htogaga*Br_htobb;
    pp_H_tt_TH=SigmaSum*Br_Htott;
    bbF_H_bb_TH=SigmabbF_H*Br_Htobb;
}

double heavyHiggs::computeThValue()
{
    return 0;
}

/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/

Hobs_ggF_H_tautau_ATLAS::Hobs_ggF_H_tautau_ATLAS(const StandardModel& SM_i)
: heavyHiggs(SM_i)
{}

double Hobs_ggF_H_tautau_ATLAS::computeThValue()
{
    if(mHh>=100.0 && mHh<=1000.0)
    {
        computeParameters();
        return ggF_H_tautau_TH/mycache->ip_ex_ggF_phi_tautau_ATLAS(mHh);
    }
    else
    {
        return 0.0;
    }
}



Hobs_ggF_H_tautau_CMS::Hobs_ggF_H_tautau_CMS(const StandardModel& SM_i)
: heavyHiggs(SM_i)
{}

double Hobs_ggF_H_tautau_CMS::computeThValue()
{
    if(mHh>=90.0 && mHh<=1000.0)
    {
        computeParameters();
        return ggF_H_tautau_TH/mycache->ip_ex_ggF_phi_tautau_CMS(mHh);
    }
    else
    {
        return 0.0;
    }
}



Hobs_bbF_H_tautau_ATLAS::Hobs_bbF_H_tautau_ATLAS(const StandardModel& SM_i)
: heavyHiggs(SM_i)
{}

double Hobs_bbF_H_tautau_ATLAS::computeThValue()
{
    if(mHh>=100.0 && mHh<=1000.0)
    {
        computeParameters();
        return bbF_H_tautau_TH/mycache->ip_ex_bbF_phi_tautau_ATLAS(mHh);
    }
    else
    {
        return 0.0;
    }
}



Hobs_bbF_H_tautau_CMS::Hobs_bbF_H_tautau_CMS(const StandardModel& SM_i)
: heavyHiggs(SM_i)
{}

double Hobs_bbF_H_tautau_CMS::computeThValue()
{
    if(mHh>=90.0 && mHh<=1000.0)
    {
        computeParameters();
        return bbF_H_tautau_TH/mycache->ip_ex_bbF_phi_tautau_CMS(mHh);
    }
    else
    {
        return 0.0;
    }
}



Hobs_ggF_H_gaga_ATLAS::Hobs_ggF_H_gaga_ATLAS(const StandardModel& SM_i)
: heavyHiggs(SM_i)
{}

double Hobs_ggF_H_gaga_ATLAS::computeThValue()
{
    if(mHh>=100.0 && mHh<=600.0)
    {
        computeParameters();
        return ggF_H_gaga_TH/mycache->ip_ex_ggF_phi_gaga_ATLAS(mHh);
    }
    else
    {
        return 0.0;
    }
}



Hobs_ggF_H_gaga_CMS::Hobs_ggF_H_gaga_CMS(const StandardModel& SM_i)
: heavyHiggs(SM_i)
{}

double Hobs_ggF_H_gaga_CMS::computeThValue()
{
    if(mHh>=150.0 && mHh<=850.0)
    {
        computeParameters();
        return ggF_H_gaga_TH/mycache->ip_ex_ggF_phi_gaga_CMS(mHh);
    }
    else
    {
        return 0.0;
    }
}



Hobs_pp_H_ZZ_CMS::Hobs_pp_H_ZZ_CMS(const StandardModel& SM_i)
: heavyHiggs(SM_i)
{}

double Hobs_pp_H_ZZ_CMS::computeThValue()
{
    if(mHh>=150.0 && mHh<=991.0)
    {
        computeParameters();
        return pp_H_ZZ_TH/mycache->ip_ex_pp_H_ZZ_CMS(mHh);
    }
    else
    {
        return 0.0;
    }
}



Hobs_ggF_H_WW_ATLAS::Hobs_ggF_H_WW_ATLAS(const StandardModel& SM_i)
: heavyHiggs(SM_i)
{}

double Hobs_ggF_H_WW_ATLAS::computeThValue()
{
    if(mHh>=300.0 && mHh<=1500.0)
    {
        computeParameters();
        return ggF_H_WW_TH/mycache->ip_ex_ggF_H_WW_ATLAS(mHh);
    }
    else
    {
        return 0.0;
    }
}



Hobs_VBF_H_WW_ATLAS::Hobs_VBF_H_WW_ATLAS(const StandardModel& SM_i)
: heavyHiggs(SM_i)
{}

double Hobs_VBF_H_WW_ATLAS::computeThValue()
{
    if(mHh>=300.0 && mHh<=1500.0)
    {
        computeParameters();
        return VBF_H_WW_TH/mycache->ip_ex_VBF_H_WW_ATLAS(mHh);
    }
    else
    {
        return 0.0;
    }
}



Hobs_ggF_H_hh_ATLAS::Hobs_ggF_H_hh_ATLAS(const StandardModel& SM_i)
: heavyHiggs(SM_i)
{}

double Hobs_ggF_H_hh_ATLAS::computeThValue()
{
    if(mHh>=260.0 && mHh<=1000.0)
    {
        computeParameters();
        return ggF_H_hh_TH/mycache->ip_ex_ggF_H_hh_ATLAS(mHh);
    }
    else
    {
        return 0.0;
    }
}



Hobs_ggF_H_hh_bbtautau_CMS::Hobs_ggF_H_hh_bbtautau_CMS(const StandardModel& SM_i)
: heavyHiggs(SM_i)
{}

double Hobs_ggF_H_hh_bbtautau_CMS::computeThValue()
{
    if(mHh>=260.0 && mHh<=350.0)
    {
        computeParameters();
        return ggF_H_hh_bbtautau_TH/mycache->ip_ex_ggF_H_hh_bbtautau_CMS(mHh);
    }
    else
    {
        return 0.0;
    }
}



Hobs_pp_H_hh_bbbb_CMS::Hobs_pp_H_hh_bbbb_CMS(const StandardModel& SM_i)
: heavyHiggs(SM_i)
{}

double Hobs_pp_H_hh_bbbb_CMS::computeThValue()
{
    if(mHh>=260.0 && mHh<=1100.0)
    {
        computeParameters();
        return pp_H_hh_bbbb_TH/mycache->ip_ex_pp_phi_hh_bbbb_CMS(mHh);
    }
    else
    {
        return 0.0;
    }
}



Hobs_pp_H_hh_gagabb_CMS::Hobs_pp_H_hh_gagabb_CMS(const StandardModel& SM_i)
: heavyHiggs(SM_i)
{}

double Hobs_pp_H_hh_gagabb_CMS::computeThValue()
{
    if(mHh>=260.0 && mHh<=1100.0)
    {
        computeParameters();
        return pp_H_hh_gagabb_TH/mycache->ip_ex_pp_phi_hh_gagabb_CMS(mHh);
    }
    else
    {
        return 0.0;
    }
}



Hobs_pp_H_tt_ATLAS::Hobs_pp_H_tt_ATLAS(const StandardModel& SM_i)
: heavyHiggs(SM_i)
{}

double Hobs_pp_H_tt_ATLAS::computeThValue()
{
    if(mHh>=400.0 && mHh<=3000.0)
    {
        computeParameters();
        return pp_H_tt_TH/mycache->ip_ex_pp_phi_tt_ATLAS(mHh);
    }
    else
    {
        return 0.0;
    }
}



Hobs_bbF_H_bb_CMS::Hobs_bbF_H_bb_CMS(const StandardModel& SM_i)
: heavyHiggs(SM_i)
{}

double Hobs_bbF_H_bb_CMS::computeThValue()
{
    if(mHh>=100.0 && mHh<=900.0)
    {
        computeParameters();
        return bbF_H_bb_TH/mycache->ip_ex_bbF_phi_bb_CMS(mHh);
    }
    else
    {
        return 0.0;
    }
}



log10_ggF_H_tautau_TH::log10_ggF_H_tautau_TH(const StandardModel& SM_i)
: heavyHiggs(SM_i)
{}

double log10_ggF_H_tautau_TH::computeThValue()
{
    computeParameters();
    return log10(ggF_H_tautau_TH);
}



log10_bbF_H_tautau_TH::log10_bbF_H_tautau_TH(const StandardModel& SM_i)
: heavyHiggs(SM_i)
{}

double log10_bbF_H_tautau_TH::computeThValue()
{
    computeParameters();
    return log10(bbF_H_tautau_TH);
}



log10_ggF_H_gaga_TH::log10_ggF_H_gaga_TH(const StandardModel& SM_i)
: heavyHiggs(SM_i)
{}

double log10_ggF_H_gaga_TH::computeThValue()
{
    computeParameters();
    return log10(ggF_H_gaga_TH);
}



log10_pp_H_ZZ_TH::log10_pp_H_ZZ_TH(const StandardModel& SM_i)
: heavyHiggs(SM_i)
{}

double log10_pp_H_ZZ_TH::computeThValue()
{
    computeParameters();
    return log10(pp_H_ZZ_TH);
}



log10_ggF_H_WW_TH::log10_ggF_H_WW_TH(const StandardModel& SM_i)
: heavyHiggs(SM_i)
{}

double log10_ggF_H_WW_TH::computeThValue()
{
    computeParameters();
    return log10(ggF_H_WW_TH);
}



log10_VBF_H_WW_TH::log10_VBF_H_WW_TH(const StandardModel& SM_i)
: heavyHiggs(SM_i)
{}

double log10_VBF_H_WW_TH::computeThValue()
{
    computeParameters();
    return log10(VBF_H_WW_TH);
}



log10_ggF_H_hh_TH::log10_ggF_H_hh_TH(const StandardModel& SM_i)
: heavyHiggs(SM_i)
{}

double log10_ggF_H_hh_TH::computeThValue()
{
    computeParameters();
    return log10(ggF_H_hh_TH);
}



log10_ggF_H_hh_bbtautau_TH::log10_ggF_H_hh_bbtautau_TH(const StandardModel& SM_i)
  : heavyHiggs(SM_i)
{}

double log10_ggF_H_hh_bbtautau_TH::computeThValue()
{
  computeParameters();
  return log10(ggF_H_hh_bbtautau_TH);
}



log10_pp_H_hh_bbbb_TH::log10_pp_H_hh_bbbb_TH(const StandardModel& SM_i)
: heavyHiggs(SM_i)
{}

double log10_pp_H_hh_bbbb_TH::computeThValue()
{
    computeParameters();
    return log10(pp_H_hh_bbbb_TH);
}



log10_pp_H_hh_gagabb_TH::log10_pp_H_hh_gagabb_TH(const StandardModel& SM_i)
: heavyHiggs(SM_i)
{}

double log10_pp_H_hh_gagabb_TH::computeThValue()
{
    computeParameters();
    return log10(pp_H_hh_gagabb_TH);
}



log10_pp_H_tt_TH::log10_pp_H_tt_TH(const StandardModel& SM_i)
: heavyHiggs(SM_i)
{}

double log10_pp_H_tt_TH::computeThValue()
{
    computeParameters();
    return log10(pp_H_tt_TH);
}



log10_bbF_H_bb_TH::log10_bbF_H_bb_TH(const StandardModel& SM_i)
: heavyHiggs(SM_i)
{}

double log10_bbF_H_bb_TH::computeThValue()
{
    computeParameters();
    return log10(bbF_H_bb_TH);
}
