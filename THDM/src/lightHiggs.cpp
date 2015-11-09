/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "lightHiggs.h"
#include "StandardModel.h"
//#include "gslpp.h"


lightHiggs::lightHiggs(const StandardModel& SM_i): 

        ThObservable(SM_i), 
        myTHDM(static_cast<const THDM*> (&SM_i)),
        mySM (SM_i)
{

}

lightHiggs::~lightHiggs()
{}

void lightHiggs::computeParameters()
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
    double s02=myTHDM->s02();
    double c02=myTHDM->c02();

    double mHl=myTHDM->getMHl();
    double mHp=myTHDM->getMHp();
    
    double m12_2=myTHDM->getM12_2();
    
    double sin_ba=myTHDM->getsin_ba();
    double sina=myTHDM->computeSina();
    double cosa=myTHDM->computeCosa();
    double sinb=myTHDM->getSinb();
    double cosb=myTHDM->getCosb();
    cos_bpa=cosb*cosa-sinb*sina;
    
    Br_htobb = myTHDM->computeBrHtobb();
    Br_htoWW = myTHDM->computeBrHtoWW();
    Br_htotautau = myTHDM->computeBrHtotautau();
    Br_htoZZ = myTHDM->computeBrHtoZZ();
    Br_htogaga = myTHDM->computeBrHtogaga();
    Br_htogg = myTHDM->computeBrHtogg();
    Br_htoZga = myTHDM->computeBrHtoZga();
    Br_htocc = myTHDM->computeBrHtocc();
    //The ggH cross section in the SM.
    SigmaggF = myTHDM->computeSigmaggH(8.0);
    //The square of the top-quark contribution to the ggH cross section in the SM
    Sigmaggh_tt = myTHDM->computeSigmaggH_tt(8.0);
    //The square of the bottom-quark contribution to the ggH cross section in the SM
    Sigmaggh_bb = myTHDM->computeSigmaggH_bb(8.0);
    //The ttH production cross section in the SM
    Sigmatth = myTHDM->computeSigmattH(8.0);
    //The VBF cross section in the SM
    SigmaVBF = myTHDM->computeSigmaVBF(8.0);
    //The WH production cross section in the SM
    SigmaWh = myTHDM->computeSigmaWH(8.0);
    //The ZH production cross section in the SM
    SigmaZh = myTHDM->computeSigmaZH(8.0);
    SigmaVh = SigmaWh + SigmaZh;
    /* pc_X is the percentage contribution of the X Higgs production mode in the SM 
     * assuming the same efficiencies in the 2HDM */
    pc_ggF = SigmaggF/(SigmaggF + Sigmatth);
    pc_tth = Sigmatth/(SigmaggF + Sigmatth);
    pc_VBF = SigmaVBF/(SigmaVBF + SigmaVh);
    pc_Vh = SigmaVh/(SigmaVBF + SigmaVh);
    /* r_ii is the ratio of the squared 2HDM vertex coupling of the SM Higgs to 
     * the particle i and the respective squared SM coupling.*/
    rh_QdQd=0.0;//It depends on the modelType
    rh_VV=sin_ba*sin_ba;
    rh_ll=0.0;//It depends on the modelType
    rh_gaga=0.0;//It depends on the modelType
    rh_gg=0.0;//It depends on the modelType 
    rh_Zga=0.0;//It depends on the modelType 
    rh_QuQu=cosa*cosa/(sinb*sinb);

    //Calulation of rh_gg, rh_QdQd, rh_ll, rh_gaga, rh_Zga (they depend on the model type): START

    //rh_gaga formula = pow(abs(I_h_ferm+I_h_W+I_h_Hp),2)/pow(abs(I_h_ferm+I_h_W),2)

    TAUu=4.0*Mu*Mu/(mHl*mHl);
    TAUc=4.0*Mc*Mc/(mHl*mHl);
    TAUt=4.0*Mt*Mt/(mHl*mHl);
    TAUd=4.0*Md*Md/(mHl*mHl);
    TAUs=4.0*Ms*Ms/(mHl*mHl);
    TAUb=4.0*Mb*Mb/(mHl*mHl);
    TAUe=4.0*Me*Me/(mHl*mHl);
    TAUmu=4.0*Mmu*Mmu/(mHl*mHl);
    TAUtau=4.0*Mtau*Mtau/(mHl*mHl);
    TAUw=4.0*MW*MW/(mHl*mHl);
    TAUhp=4.0*mHp*mHp/(mHl*mHl);

    gslpp::complex I_h_ferm=0.0;//It depends on the modelType
    gslpp::complex fermU=-(8./3.)*(TAUu*(1+(1-TAUu)*myfunctions.f_func(TAUu))
                         +TAUc*(1+(1-TAUc)*myfunctions.f_func(TAUc))+TAUt*(1+(1-TAUt)*myfunctions.f_func(TAUt)));
    gslpp::complex fermD=-(2./3.)*(TAUd*(1+(1-TAUd)*myfunctions.f_func(TAUd))
                         +TAUs*(1+(1-TAUs)*myfunctions.f_func(TAUs))+TAUb*(1+(1-TAUb)*myfunctions.f_func(TAUb)));
    gslpp::complex fermL=-2.*(TAUe*(1+(1-TAUe)*myfunctions.f_func(TAUe))
                         +TAUmu*(1+(1-TAUmu)*myfunctions.f_func(TAUmu))
                         +TAUtau*(1+(1-TAUtau)*myfunctions.f_func(TAUtau)));
    gslpp::complex I_h_W=sin_ba*(2.0 + 3.0*TAUw + 3.0*TAUw*(2.0-TAUw)*myfunctions.f_func(TAUw));
    double ghHpHm = ((mHl*mHl -2.0*mHp*mHp)*sin_ba
                    -(mHl*mHl -m12_2/(cosb*sinb))/(cosb*sinb)*cos_bpa)/vev;
    gslpp::complex I_h_Hp=-TAUhp*(1.0-TAUhp*myfunctions.f_func(TAUhp))*vev/(2.*mHp*mHp)*ghHpHm;

    double ABS1=0.0;
    double ABS2=0.0;

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
    
    A_h_F = 0.0;//It depends on the modelType
    A_h_U = -4.0*(1.0/2.0-4.0/3.0*s02)*(myfunctions.Int1(TAUu,LAMu)+myfunctions.Int1(TAUc,LAMc)
                           +myfunctions.Int1(TAUt,LAMt)-myfunctions.Int2(TAUu,LAMu)-myfunctions.Int2(TAUc,LAMc)-myfunctions.Int2(TAUt,LAMt));
    A_h_D = +2.0*(-1.0/2.0+2.0/3.0*s02)*(myfunctions.Int1(TAUd,LAMd)+myfunctions.Int1(TAUs,LAMs)
                           +myfunctions.Int1(TAUb,LAMb)-myfunctions.Int2(TAUd,LAMd)-myfunctions.Int2(TAUs,LAMs)-myfunctions.Int2(TAUb,LAMb));
    A_h_L  = +2.0*(-1.0/2.0+2.0*s02)*(myfunctions.Int1(TAUe,LAMe)+myfunctions.Int1(TAUmu,LAMmu)
                            +myfunctions.Int1(TAUtau,LAMtau)-myfunctions.Int2(TAUe,LAMe)-myfunctions.Int2(TAUmu,LAMmu)
                            -myfunctions.Int2(TAUtau,LAMtau));
    A_h_W = -sin_ba*sqrt(c02/s02)*(4*(3-s02/c02)*myfunctions.Int2(TAUw,LAMw)
                           +((1.0+2.0/TAUw)*s02/c02-(5.0+2.0/TAUw))*myfunctions.Int1(TAUw,LAMw));
    A_h_Hp = ghHpHm*(1-2.0*s02)/sqrt(c02*s02)*myfunctions.Int1(TAUhp,LAMhp)
                            *vev/(2.*mHp*mHp);

    ABS3=0.0;
    ABS4=0.0;

    modelType=myTHDM->getModelType();

    switch(modelType){
        case 1 ://type 1
        rh_gg=cosa/sinb*cosa/sinb;
        rh_QdQd=cosa/sinb*cosa/sinb;
        rh_ll=cosa/sinb*cosa/sinb;
        I_h_ferm=cosa/sinb*(fermU+fermD+fermL);
        A_h_F = cosa/sinb*(A_h_U+A_h_D+A_h_L)/sqrt(s02*c02); 
        break;
        case 2 ://type 2
        rh_gg=-cosa/sinb*sina/cosb+(cosa/sinb+sina/cosb)
             *(Sigmaggh_tt*cosa/sinb+Sigmaggh_bb*sina/cosb)/SigmaggF;
        rh_QdQd=sina/cosb*sina/cosb;
        rh_ll=sina/cosb*sina/cosb;
        I_h_ferm=cosa/sinb*fermU -sina/cosb*(fermD+fermL);
        A_h_F = (cosa/sinb*A_h_U-sina/cosb*(A_h_D+A_h_L))/sqrt(s02*c02);
        break;
        case 3 ://Lepton-specific
        rh_gg=cosa/sinb*cosa/sinb;
        rh_QdQd=cosa/sinb*cosa/sinb;
        rh_ll=sina/cosb*sina/cosb;
        I_h_ferm = cosa/sinb*(fermU+fermD) -sina/cosb*fermL;
        A_h_F = (cosa/sinb*(A_h_U+A_h_D)-sina/cosb*A_h_L)/sqrt(s02*c02);
        break;
        case 4 ://Flipped
        rh_gg=-cosa/sinb*sina/cosb+(cosa/sinb+sina/cosb)
             *(Sigmaggh_tt*cosa/sinb+Sigmaggh_bb*sina/cosb)/SigmaggF;
        rh_QdQd=sina/cosb*sina/cosb;
        rh_ll=cosa/sinb*cosa/sinb;
        I_h_ferm = cosa/sinb*(fermU+fermL) -sina/cosb*fermD;
        A_h_F = (cosa/sinb*(A_h_U+A_h_L)-sina/cosb*A_h_D)/sqrt(s02*c02);
        break;
        default :
        std::stringstream out;
        out << modelType;
        //throw std::runtime_error(“lightHiggs::modelType " + out.str() + " not implemented: it can only be 1,2,3 or 4.”);
        throw std::runtime_error("modelType can be only any of 1,2,3,4");
    }
    
//    std::cout<<"rh_gg: "<<rh_gg<<std::endl;
    
    ABS1=(I_h_ferm+I_h_W+I_h_Hp).abs();
    ABS2=(I_h_ferm+I_h_W).abs();
    rh_gaga=ABS1*ABS1/(ABS2*ABS2);
    
//    std::cout<<"rh_gaga: "<<rh_gaga<<std::endl;
    
    ABS3=(A_h_F+A_h_W+A_h_Hp).abs();
    ABS4=(A_h_F+A_h_W).abs();
    rh_Zga=ABS3*ABS3/(ABS4*ABS4);

    //Calulation of rh_gg, rh_QdQd, rh_ll, rh_gaga, rh_Zga (they depend on the model type): FINISH

    /* ggF_ttH is in the formulas of the Higgs signal strengths with ggF or ttH
     * as Higgs production mode*/
    ggF_tth = pc_ggF*rh_gg + pc_tth*rh_QuQu;
    /* VBF_VH is in the formulas of the Higgs signal strengths with VBF or VH
     * as Higgs production mode*/
    VBF_Vh = (pc_VBF + pc_Vh ) * rh_VV;

    sum = rh_QdQd*Br_htobb + rh_VV*(Br_htoWW+Br_htoZZ) + rh_ll*Br_htotautau +
          rh_gaga*Br_htogaga + rh_gg*Br_htogg + rh_Zga*Br_htoZga + rh_QuQu*Br_htocc;
  
}



double lightHiggs::computeThValue()
{
    return 0;
}

double lightHiggs::THDM_BR_h_bb()
{
   computeParameters();
   return rh_QdQd*Br_htobb/sum;
}

double lightHiggs::THDM_BR_h_gaga()
{
   computeParameters();
   return rh_gaga*Br_htogaga/sum;
}

double lightHiggs::THDM_BR_h_tautau()
{
   computeParameters();
   return rh_ll*Br_htotautau/sum;
}


/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/



ggF_tth_htobb::ggF_tth_htobb(const StandardModel& SM_i)
: lightHiggs(SM_i)
{}

double ggF_tth_htobb::computeThValue()
{
    computeParameters();
    
    return ggF_tth*rh_QdQd/sum;
}



ggF_tth_htoWW::ggF_tth_htoWW(const StandardModel& SM_i)
: lightHiggs(SM_i)
{}

double ggF_tth_htoWW::computeThValue()
{
    computeParameters();
    
    return ggF_tth*rh_VV/sum;
}



ggF_tth_htotautau::ggF_tth_htotautau(const StandardModel& SM_i)
: lightHiggs(SM_i)
{}

double ggF_tth_htotautau::computeThValue()
{
    computeParameters();
    
    return ggF_tth*rh_ll/sum;
}



ggF_tth_htoZZ::ggF_tth_htoZZ(const StandardModel& SM_i)
: lightHiggs(SM_i)
{}

double ggF_tth_htoZZ::computeThValue()
{
    computeParameters();
    
    return ggF_tth*rh_VV/sum;
}



ggF_tth_htogaga::ggF_tth_htogaga(const StandardModel& SM_i)
: lightHiggs(SM_i)
{}

double ggF_tth_htogaga::computeThValue()
{
    computeParameters();
    
    return ggF_tth*rh_gaga/sum;
}



VBF_Vh_htobb::VBF_Vh_htobb(const StandardModel& SM_i)
: lightHiggs(SM_i)
{}

double VBF_Vh_htobb::computeThValue()
{
    computeParameters();
    
    return VBF_Vh*rh_QdQd/sum;
}



VBF_Vh_htoWW::VBF_Vh_htoWW(const StandardModel& SM_i)
: lightHiggs(SM_i)
{}

double VBF_Vh_htoWW::computeThValue()
{
    computeParameters();
    
    return VBF_Vh*rh_VV/sum;
}



VBF_Vh_htotautau::VBF_Vh_htotautau(const StandardModel& SM_i)
: lightHiggs(SM_i)
{}

double VBF_Vh_htotautau::computeThValue()
{
    computeParameters();
    
    return VBF_Vh*rh_ll/sum;
}



VBF_Vh_htoZZ::VBF_Vh_htoZZ(const StandardModel& SM_i)
: lightHiggs(SM_i)
{}

double VBF_Vh_htoZZ::computeThValue()
{
    computeParameters();
    
    return VBF_Vh*rh_VV/sum;
}

VBF_Vh_htogaga::VBF_Vh_htogaga(const StandardModel& SM_i)
: lightHiggs(SM_i)
{}

double VBF_Vh_htogaga::computeThValue()
{
    computeParameters();
    
    return VBF_Vh*rh_gaga/sum;
}
