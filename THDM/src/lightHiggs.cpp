/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "lightHiggs.h"
#include "StandardModel.h"


lightHiggs::lightHiggs(const StandardModel& SM_i): 

        ThObservable(SM_i), 
        myTHDM(static_cast<const THDM*> (&SM_i)),
        mySM (SM_i)
{}

lightHiggs::~lightHiggs()
{}

void lightHiggs::computeSignalStrengthQuantities()
{

    THDMfunctions myfunctions(mySM);

    double Mu=myTHDM->getQuarks(QCD::UP).getMass();
    double Md=myTHDM->getQuarks(QCD::DOWN).getMass();
    double Mc=myTHDM->getQuarks(QCD::CHARM).getMass();
    double Ms=myTHDM->getQuarks(QCD::STRANGE).getMass();
    double Mt=myTHDM->getQuarks(QCD::TOP).getMass();
    double Mb=myTHDM->getQuarks(QCD::BOTTOM).getMass();   
    double Me=myTHDM->getLeptons(StandardModel::ELECTRON).getMass();
    double Mmu=myTHDM->getLeptons(StandardModel::MU).getMass();
    double Mtau=myTHDM->getLeptons(StandardModel::TAU).getMass();
    double MZ=myTHDM->getMz();
    double MW=myTHDM->Mw();
    double vev=myTHDM->v();
    double sW2=myTHDM->sW2();
    double cW2=myTHDM->cW2();
    double mHl=myTHDM->getMHl();
    double mHp=myTHDM->getmHp();
    double m12_2=myTHDM->getm12_2();    
    double sin_ba=myTHDM->getsin_ba();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double cos_bpa=cosb*cosa-sinb*sina;

    //The Standard Model h branching ratios

    BrSM_htobb = myTHDM->computeBrHtobb();
    BrSM_htotautau = myTHDM->computeBrHtotautau();
    BrSM_htogaga = myTHDM->computeBrHtogaga();
    double BrSM_htoWW = myTHDM->computeBrHtoWW();
    double BrSM_htoZZ = myTHDM->computeBrHtoZZ();
    double BrSM_htogg = myTHDM->computeBrHtogg();
    double BrSM_htoZga = myTHDM->computeBrHtoZga();
    double BrSM_htocc = myTHDM->computeBrHtocc();

    //The ggH cross section in the SM.
    double SigmaggF = myTHDM->computeSigmaggH(8.0);
    //The square of the top-quark contribution to the ggH cross section in the SM
    double Sigmaggh_tt = myTHDM->computeSigmaggH_tt(8.0);
    //The square of the bottom-quark contribution to the ggH cross section in the SM
    double Sigmaggh_bb = myTHDM->computeSigmaggH_bb(8.0);
    //The ttH production cross section in the SM
    double Sigmatth = myTHDM->computeSigmattH(8.0);
    //The VBF cross section in the SM
    double SigmaVBF = myTHDM->computeSigmaVBF(8.0);
    //The WH production cross section in the SM
    double SigmaWh = myTHDM->computeSigmaWH(8.0);
    //The ZH production cross section in the SM
    double SigmaZh = myTHDM->computeSigmaZH(8.0);
    double SigmaVh = SigmaWh + SigmaZh;

    /* r_ii is the ratio of the squared 2HDM vertex coupling of h to
     * the particle i and the respective squared SM coupling.*/
    rh_QuQu=cosa*cosa/(sinb*sinb);
    rh_VV=sin_ba*sin_ba;
    rh_QdQd=0.0;//It depends on the modelType
    rh_ll=0.0;//It depends on the modelType
    rh_gg=0.0;//It depends on the modelType 

    //Calulation of rh_gg, rh_QdQd, rh_ll, rh_gaga, rh_Zga (depending on the model type): START

    //rh_gaga formula = abs(I_h_F+I_h_W+I_h_Hp)^2 / abs(I_hSM_F+I_hSM_W)^2

    double TAUu=4.0*Mu*Mu/(mHl*mHl);
    double TAUc=4.0*Mc*Mc/(mHl*mHl);
    double TAUt=4.0*Mt*Mt/(mHl*mHl);
    double TAUd=4.0*Md*Md/(mHl*mHl);
    double TAUs=4.0*Ms*Ms/(mHl*mHl);
    double TAUb=4.0*Mb*Mb/(mHl*mHl);
    double TAUe=4.0*Me*Me/(mHl*mHl);
    double TAUmu=4.0*Mmu*Mmu/(mHl*mHl);
    double TAUtau=4.0*Mtau*Mtau/(mHl*mHl);
    double TAUw=4.0*MW*MW/(mHl*mHl);
    double TAUhp=4.0*mHp*mHp/(mHl*mHl);
//
    gslpp::complex I_h_F=0.0;//It depends on the modelType
    gslpp::complex fermU=-(8./3.)*(TAUu*(1+(1-TAUu)*myfunctions.f_func(TAUu))
                         +TAUc*(1+(1-TAUc)*myfunctions.f_func(TAUc))+TAUt*(1+(1-TAUt)*myfunctions.f_func(TAUt)));
    gslpp::complex fermD=-(2./3.)*(TAUd*(1+(1-TAUd)*myfunctions.f_func(TAUd))
                         +TAUs*(1+(1-TAUs)*myfunctions.f_func(TAUs))+TAUb*(1+(1-TAUb)*myfunctions.f_func(TAUb)));
    gslpp::complex fermL=-2.*(TAUe*(1+(1-TAUe)*myfunctions.f_func(TAUe))
                         +TAUmu*(1+(1-TAUmu)*myfunctions.f_func(TAUmu))
                         +TAUtau*(1+(1-TAUtau)*myfunctions.f_func(TAUtau)));
    gslpp::complex I_hSM_W=2.0 + 3.0*TAUw + 3.0*TAUw*(2.0-TAUw)*myfunctions.f_func(TAUw);
    gslpp::complex I_h_W=sin_ba*I_hSM_W;
    double ghHpHm = ((mHl*mHl -2.0*mHp*mHp)*sin_ba
                    -(mHl*mHl -m12_2/(cosb*sinb))/(cosb*sinb)*cos_bpa)/vev;
    gslpp::complex I_h_Hp=-TAUhp*(1.0-TAUhp*myfunctions.f_func(TAUhp))*vev/(2.*mHp*mHp)*ghHpHm;

    double ABSgagaTHDM=0.0;
    double ABSgagaSM=0.0;

    //rh_Zga formula = abs(A_h_F+A_h_W+A_h_Hp)^2 / abs(A_hSM_F+A_hSM_W)^2

    double LAMu=4.0*Mu*Mu/(MZ*MZ);
    double LAMc=4.0*Mc*Mc/(MZ*MZ);
    double LAMt=4.0*Mt*Mt/(MZ*MZ);
    double LAMd=4.0*Md*Md/(MZ*MZ);
    double LAMs=4.0*Ms*Ms/(MZ*MZ);
    double LAMb=4.0*Mb*Mb/(MZ*MZ);
    double LAMe=4.0*Me*Me/(MZ*MZ);
    double LAMmu=4.0*Mmu*Mmu/(MZ*MZ);
    double LAMtau=4.0*Mtau*Mtau/(MZ*MZ);
    double LAMw=4.0*MW*MW/(MZ*MZ);
    double LAMhp=4.0*mHp*mHp/(MZ*MZ);

    gslpp::complex A_h_F = 0.0;//It depends on the modelType
    gslpp::complex A_h_U = -4.0*(1.0/2.0-4.0/3.0*sW2)*(myfunctions.Int1(TAUu,LAMu)+myfunctions.Int1(TAUc,LAMc)
                           +myfunctions.Int1(TAUt,LAMt)-myfunctions.Int2(TAUu,LAMu)-myfunctions.Int2(TAUc,LAMc)-myfunctions.Int2(TAUt,LAMt));
    gslpp::complex A_h_D = +2.0*(-1.0/2.0+2.0/3.0*sW2)*(myfunctions.Int1(TAUd,LAMd)+myfunctions.Int1(TAUs,LAMs)
                           +myfunctions.Int1(TAUb,LAMb)-myfunctions.Int2(TAUd,LAMd)-myfunctions.Int2(TAUs,LAMs)-myfunctions.Int2(TAUb,LAMb));
    gslpp::complex A_h_L  = +2.0*(-1.0/2.0+2.0*sW2)*(myfunctions.Int1(TAUe,LAMe)+myfunctions.Int1(TAUmu,LAMmu)
                            +myfunctions.Int1(TAUtau,LAMtau)-myfunctions.Int2(TAUe,LAMe)-myfunctions.Int2(TAUmu,LAMmu)
                            -myfunctions.Int2(TAUtau,LAMtau));
    gslpp::complex A_hSM_W = -sqrt(cW2/sW2)*(4*(3-sW2/cW2)*myfunctions.Int2(TAUw,LAMw)
                           +((1.0+2.0/TAUw)*sW2/cW2-(5.0+2.0/TAUw))*myfunctions.Int1(TAUw,LAMw));
    gslpp::complex A_h_W = sin_ba*A_hSM_W;
    gslpp::complex A_h_Hp = ghHpHm*(1-2.0*sW2)/sqrt(cW2*sW2)*myfunctions.Int1(TAUhp,LAMhp)
                            *vev/(2.*mHp*mHp);

    double ABSZgaTHDM=0.0;
    double ABSZgaSM=0.0;

    std::string modelflag=myTHDM->getModelTypeflag();

    if( modelflag == "type1" ) {
        rh_gg=cosa/sinb*cosa/sinb;
        rh_QdQd=cosa/sinb*cosa/sinb;
        rh_ll=cosa/sinb*cosa/sinb;
        I_h_F=cosa/sinb*(fermU+fermD+fermL);
        A_h_F = cosa/sinb*(A_h_U+A_h_D+A_h_L)/sqrt(sW2*cW2);
    }
    else if( modelflag == "type2" ) {
        rh_gg=-cosa/sinb*sina/cosb+(cosa/sinb+sina/cosb)
             *(Sigmaggh_tt*cosa/sinb+Sigmaggh_bb*sina/cosb)/SigmaggF;
        rh_QdQd=sina/cosb*sina/cosb;
        rh_ll=sina/cosb*sina/cosb;
        I_h_F=cosa/sinb*fermU -sina/cosb*(fermD+fermL);
        A_h_F = (cosa/sinb*A_h_U-sina/cosb*(A_h_D+A_h_L))/sqrt(sW2*cW2);
    }
    else if( modelflag == "typeX" ) {
        rh_gg=cosa/sinb*cosa/sinb;
        rh_QdQd=cosa/sinb*cosa/sinb;
        rh_ll=sina/cosb*sina/cosb;
        I_h_F = cosa/sinb*(fermU+fermD) -sina/cosb*fermL;
        A_h_F = (cosa/sinb*(A_h_U+A_h_D)-sina/cosb*A_h_L)/sqrt(sW2*cW2);
    }
    else if( modelflag == "typeY" ) {
        rh_gg=-cosa/sinb*sina/cosb+(cosa/sinb+sina/cosb)
             *(Sigmaggh_tt*cosa/sinb+Sigmaggh_bb*sina/cosb)/SigmaggF;
        rh_QdQd=sina/cosb*sina/cosb;
        rh_ll=cosa/sinb*cosa/sinb;
        I_h_F = cosa/sinb*(fermU+fermL) -sina/cosb*fermD;
        A_h_F = (cosa/sinb*(A_h_U+A_h_L)-sina/cosb*A_h_D)/sqrt(sW2*cW2);
    }
    else {
        throw std::runtime_error("modelflag can be only any of \"type1\", \"type2\", \"typeX\" or \"typeY\"");
    }

    ABSgagaTHDM=(I_h_F+I_h_W+I_h_Hp).abs2();
    ABSgagaSM=(fermU+fermL+fermD+I_h_W).abs2();
    rh_gaga=ABSgagaTHDM/ABSgagaSM;

    ABSZgaTHDM=(A_h_F+A_h_W+A_h_Hp).abs2();
    ABSZgaSM=(A_h_U+A_h_L+A_h_D+A_hSM_W).abs2();
    rh_Zga=ABSZgaTHDM/ABSZgaSM;

    //Calulation of rh_gg, rh_QdQd, rh_ll, rh_gaga, rh_Zga (they depend on the model type): END

    /* ggF_tth is the ratio of the THDM and SM cross sections for ggF or tth production */
    ggF_tth = (SigmaggF*rh_gg + Sigmatth*rh_QuQu)/(SigmaggF + Sigmatth);
    /* VBF_Vh is the ratio of the THDM and SM cross sections for VBF or Vh production */
    VBF_Vh = rh_VV;

    sumModBRs = rh_QdQd*BrSM_htobb + rh_VV*(BrSM_htoWW+BrSM_htoZZ) + rh_ll*BrSM_htotautau +
          rh_gaga*BrSM_htogaga + rh_gg*BrSM_htogg + rh_Zga*BrSM_htoZga + rh_QuQu*BrSM_htocc;
}

double lightHiggs::computeThValue()
{
    return 0;
}

double lightHiggs::THDM_BR_h_bb()
{
   computeSignalStrengthQuantities();
   return rh_QdQd*BrSM_htobb/sumModBRs;
}

double lightHiggs::THDM_BR_h_gaga()
{
   computeSignalStrengthQuantities();
   return rh_gaga*BrSM_htogaga/sumModBRs;
}

double lightHiggs::THDM_BR_h_tautau()
{
   computeSignalStrengthQuantities();
   return rh_ll*BrSM_htotautau/sumModBRs;
}


/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/



ggF_tth_htobb::ggF_tth_htobb(const StandardModel& SM_i)
: lightHiggs(SM_i)
{}

double ggF_tth_htobb::computeThValue()
{
    computeSignalStrengthQuantities();
    
    return ggF_tth*rh_QdQd/sumModBRs;
}



ggF_tth_htoWW::ggF_tth_htoWW(const StandardModel& SM_i)
: lightHiggs(SM_i)
{}

double ggF_tth_htoWW::computeThValue()
{
    computeSignalStrengthQuantities();
    
    return ggF_tth*rh_VV/sumModBRs;
}



ggF_tth_htotautau::ggF_tth_htotautau(const StandardModel& SM_i)
: lightHiggs(SM_i)
{}

double ggF_tth_htotautau::computeThValue()
{
    computeSignalStrengthQuantities();
    
    return ggF_tth*rh_ll/sumModBRs;
}



ggF_tth_htoZZ::ggF_tth_htoZZ(const StandardModel& SM_i)
: lightHiggs(SM_i)
{}

double ggF_tth_htoZZ::computeThValue()
{
    computeSignalStrengthQuantities();
    
    return ggF_tth*rh_VV/sumModBRs;
}



ggF_tth_htogaga::ggF_tth_htogaga(const StandardModel& SM_i)
: lightHiggs(SM_i)
{}

double ggF_tth_htogaga::computeThValue()
{
    computeSignalStrengthQuantities();
    
    return ggF_tth*rh_gaga/sumModBRs;
}



VBF_Vh_htobb::VBF_Vh_htobb(const StandardModel& SM_i)
: lightHiggs(SM_i)
{}

double VBF_Vh_htobb::computeThValue()
{
    computeSignalStrengthQuantities();
    
    return VBF_Vh*rh_QdQd/sumModBRs;
}



VBF_Vh_htoWW::VBF_Vh_htoWW(const StandardModel& SM_i)
: lightHiggs(SM_i)
{}

double VBF_Vh_htoWW::computeThValue()
{
    computeSignalStrengthQuantities();
    
    return VBF_Vh*rh_VV/sumModBRs;
}



VBF_Vh_htotautau::VBF_Vh_htotautau(const StandardModel& SM_i)
: lightHiggs(SM_i)
{}

double VBF_Vh_htotautau::computeThValue()
{
    computeSignalStrengthQuantities();
    
    return VBF_Vh*rh_ll/sumModBRs;
}



VBF_Vh_htoZZ::VBF_Vh_htoZZ(const StandardModel& SM_i)
: lightHiggs(SM_i)
{}

double VBF_Vh_htoZZ::computeThValue()
{
    computeSignalStrengthQuantities();
    
    return VBF_Vh*rh_VV/sumModBRs;
}

VBF_Vh_htogaga::VBF_Vh_htogaga(const StandardModel& SM_i)
: lightHiggs(SM_i)
{}

double VBF_Vh_htogaga::computeThValue()
{
    computeSignalStrengthQuantities();
    
    return VBF_Vh*rh_gaga/sumModBRs;
}
