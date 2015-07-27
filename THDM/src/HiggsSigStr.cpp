/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "HiggsSigStr.h"
#include "StandardModel.h"
#include "gslpp.h"

HiggsSigStr::HiggsSigStr(const StandardModel& SM_i, int obsFlag) 
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{
    if (obsFlag > 0 and obsFlag < 11) obs = obsFlag;
    else throw std::runtime_error("obsFlag in HiggsSigStr() called from "
            "ThFactory::ThFactory() can only be 1 (ggF_tth_htobb) "
            "or 2 (ggF_tth_htoWW) or 3 (ggF_tth_htotautau) "
            "or 4 (ggF_tth_htoZZ) or 5 (ggF_tth_htogaga) "
            "or 6 (VBF_Vh_htobb) or 7 (VBF_Vh_htoWW) "
            "or 8 (VBF_Vh_htotautau) or 9 (VBF_Vh_htoZZ) "
            "or 10 (VBF_Vh_htogaga) ");
};

double HiggsSigStr::computeThValue()
{

    double Mu=myTHDM->getQuarks(QCD::UP).getMass();
    double Md=myTHDM->getQuarks(QCD::DOWN).getMass();
    double Mc=myTHDM->getQuarks(QCD::CHARM).getMass();
    double Ms=myTHDM->getQuarks(QCD::STRANGE).getMass();
    double Mt=myTHDM->getQuarks(QCD::TOP).getMass();
    double Mb=myTHDM->getQuarks(QCD::BOTTOM).getMass();
    double Me=myTHDM->getLeptons(StandardModel::ELECTRON).getMass();
    double Mmu=myTHDM->getLeptons(StandardModel::MU).getMass();
    double Mtau=myTHDM->getLeptons(StandardModel::TAU).getMass();
    double MW=myTHDM->Mw();
    double MZ=myTHDM->getMz();
    double s02=myTHDM->s02();
    double c02=myTHDM->c02();
    double mHl=myTHDM->getMHl();
    double mHp=myTHDM->getMHp();
    double sin_ba=myTHDM->getsin_ba();
    double sina=myTHDM->computeSina();
    double cosa=myTHDM->computeCosa();
    double sinb=myTHDM->getSinb();
    double cosb=myTHDM->getCosb();
    double cos_bpa=cosb*cosa-sinb*sina;
    double m12_2=myTHDM->getM12_2();
    double vev=myTHDM->v();
    double Br_htobb = myTHDM->computeBrHtobb();
    double Br_htoWW = myTHDM->computeBrHtoWW();
    double Br_htotautau = myTHDM->computeBrHtotautau();
    double Br_htoZZ = myTHDM->computeBrHtoZZ();
    double Br_htogaga = myTHDM->computeBrHtogaga();
    double Br_htogg = myTHDM->computeBrHtogg();
    double Br_htoZga = myTHDM->computeBrHtoZga();
    double Br_htocc = myTHDM->computeBrHtocc();
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
    /* pc_X is the percentage contribution of the X Higgs production mode in the SM 
     * assuming the same efficiencies in the 2HDM */
    double pc_ggF = SigmaggF/(SigmaggF + Sigmatth);
    double pc_tth = Sigmatth/(SigmaggF + Sigmatth);
    double pc_VBF = SigmaVBF/(SigmaVBF + SigmaVh);
    double pc_Vh = SigmaVh/(SigmaVBF + SigmaVh);
    /* r_ii is the ratio of the squared 2HDM vertex coupling of the SM Higgs to 
     * the particle i and the respective squared SM coupling.*/
    double rh_QdQd=0.0;//It depends on the modelType
    double rh_VV=sin_ba*sin_ba;
    double rh_ll=0.0;//It depends on the modelType
    double rh_gaga=0.0;//It depends on the modelType
    double rh_gg=0.0;//It depends on the modelType 
    double rh_Zga=0.0;//It depends on the modelType 
    double rh_QuQu=cosa*cosa/(sinb*sinb);

    //Calulation of rh_gg, rh_QdQd, rh_ll, rh_gaga, rh_Zga (they depend on the model type): START

    //rh_gaga formula = pow(abs(I_h_ferm+I_h_W+I_h_Hp),2)/pow(abs(I_h_ferm+I_h_W),2)

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

    gslpp::complex I_h_ferm=0.0;//It depends on the modelType
    gslpp::complex fermU=-(8./3.)*(TAUu*(1+(1-TAUu)*f_func(TAUu))
                         +TAUc*(1+(1-TAUc)*f_func(TAUc))+TAUt*(1+(1-TAUt)*f_func(TAUt)));
    gslpp::complex fermD=-(2./3.)*(TAUd*(1+(1-TAUd)*f_func(TAUd))
                         +TAUs*(1+(1-TAUs)*f_func(TAUs))+TAUb*(1+(1-TAUb)*f_func(TAUb)));
    gslpp::complex fermL=-2.*(TAUe*(1+(1-TAUe)*f_func(TAUe))
                         +TAUmu*(1+(1-TAUmu)*f_func(TAUmu))
                         +TAUtau*(1+(1-TAUtau)*f_func(TAUtau)));
    gslpp::complex I_h_W=sin_ba*(2.0 + 3.0*TAUw + 3.0*TAUw*(2.0-TAUw)*f_func(TAUw));
    double ghHpHm = ((mHl*mHl -2.0*mHp*mHp)*sin_ba
                    -(mHl*mHl -m12_2/(cosb*sinb))/(cosb*sinb)*cos_bpa)/vev;
    gslpp::complex I_h_Hp=-TAUhp*(1.0-TAUhp*f_func(TAUhp))*vev/(2.*mHp*mHp)*ghHpHm;

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
    
    gslpp::complex A_h_F = 0.0;//It depends on the modelType
    gslpp::complex A_h_U = -4.0*(1.0/2.0-4.0/3.0*s02)*(Int1(TAUu,LAMu)+Int1(TAUc,LAMc)
                           +Int1(TAUt,LAMt)-Int2(TAUu,LAMu)-Int2(TAUc,LAMc)-Int2(TAUt,LAMt));
    gslpp::complex A_h_D = +2.0*(-1.0/2.0+2.0/3.0*s02)*(Int1(TAUd,LAMd)+Int1(TAUs,LAMs)
                           +Int1(TAUb,LAMb)-Int2(TAUd,LAMd)-Int2(TAUs,LAMs)-Int2(TAUb,LAMb));
    gslpp::complex A_h_L  = +2.0*(-1.0/2.0+2.0*s02)*(Int1(TAUe,LAMe)+Int1(TAUmu,LAMmu)
                            +Int1(TAUtau,LAMtau)-Int2(TAUe,LAMe)-Int2(TAUmu,LAMmu)
                            -Int2(TAUtau,LAMtau));
    gslpp::complex A_h_W = -sin_ba*sqrt(c02/s02)*(4*(3-s02/c02)*Int2(TAUw,LAMw)
                           +((1.0+2.0/TAUw)*s02/c02-(5.0+2.0/TAUw))*Int1(TAUw,LAMw));
    gslpp::complex A_h_Hp = ghHpHm*(1-2.0*s02)/sqrt(c02*s02)*Int1(TAUhp,LAMhp)
                            *vev/(2.*mHp*mHp);

    double ABS3=0.0;
    double ABS4=0.0;

    int modelType=myTHDM->getModelType();

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
        //throw std::runtime_error(“HiggsSigStr::modelType " + out.str() + " not implemented: it can only be 1,2,3 or 4.”);
        throw std::runtime_error("modelType can be only any of 1,2,3,4");
    }
    
    ABS1=(I_h_ferm+I_h_W+I_h_Hp).abs();
    ABS2=(I_h_ferm+I_h_W).abs();
    rh_gaga=ABS1*ABS1/(ABS2*ABS2);
    
    ABS3=(A_h_F+A_h_W+A_h_Hp).abs();
    ABS4=(A_h_F+A_h_W).abs();
    rh_Zga=ABS3*ABS3/(ABS4*ABS4);

//    std::cout << "I_h_ferm = " << I_h_ferm << std::endl;
//    std::cout << "I_h_W = " << I_h_W << std::endl;
//    std::cout << "I_h_Hp = " << I_h_Hp << std::endl;
//
//    std::cout << "A_h_F = " << A_h_F << std::endl;
//    std::cout << "A_h_W = " << A_h_W << std::endl;
//    std::cout << "A_h_Hp = " << A_h_Hp << std::endl;

    //Calulation of rh_gg, rh_QdQd, rh_ll, rh_gaga, rh_Zga (they depend on the model type): FINISH

    /* ggF_ttH is in the formulas of the Higgs signal strengths with ggF or ttH
     * as Higgs production mode*/
    double ggF_tth = pc_ggF*rh_gg + pc_tth*rh_QuQu;
    /* VBF_VH is in the formulas of the Higgs signal strengths with VBF or VH
     * as Higgs production mode*/
    double VBF_Vh = (pc_VBF + pc_Vh ) * rh_VV;

    double sum = rh_QdQd*Br_htobb + rh_VV*(Br_htoWW+Br_htoZZ) + rh_ll*Br_htotautau +
                 rh_gaga*Br_htogaga + rh_gg*Br_htogg + rh_Zga*Br_htoZga + rh_QuQu*Br_htocc;
    
    //HIGGS SIGNAL STRENGTHS
    
    double ggF_tth_htobb= ggF_tth*rh_QdQd/sum;//ggF_ttH means ggF+ttH
    double ggF_tth_htoWW= ggF_tth*rh_VV/sum;
    double ggF_tth_htotautau= ggF_tth*rh_ll/sum;
    double ggF_tth_htoZZ= ggF_tth*rh_VV/sum;
    double ggF_tth_htogaga= ggF_tth*rh_gaga/sum;
    double VBF_Vh_htobb= VBF_Vh*rh_QdQd/sum;//VBF_VH means VBF+VH
    double VBF_Vh_htoWW= VBF_Vh*rh_VV/sum;
    double VBF_Vh_htotautau= VBF_Vh*rh_ll/sum;
    double VBF_Vh_htoZZ= VBF_Vh*rh_VV/sum;
    double VBF_Vh_htogaga= VBF_Vh*rh_gaga/sum;

    if (obs ==1) return(ggF_tth_htobb);
    if (obs ==2) return(ggF_tth_htoWW);
    if (obs ==3) return(ggF_tth_htotautau);
    if (obs ==4) return(ggF_tth_htoZZ);
    if (obs ==5) return(ggF_tth_htogaga);
    if (obs ==6) return(VBF_Vh_htobb);
    if (obs ==7) return(VBF_Vh_htoWW);
    if (obs ==8) return(VBF_Vh_htotautau);
    if (obs ==9) return(VBF_Vh_htoZZ);
    if (obs ==10) return(VBF_Vh_htogaga);

    throw std::runtime_error("HiggsSigStr::computeThValue(): Observable type not "
            "defined. Can be only any of (1,2,3,4,5,6,7,8,9,10)");
    return (EXIT_FAILURE);

}

gslpp::complex HiggsSigStr::f_func(const double x) const {
    if(x<1) {
    gslpp::complex z = -gslpp::complex::i()*M_PI;
    return -pow(log((1+sqrt(1-x))/(1-sqrt(1-x)))+z,2)/4.0;
    }
    else {
        return pow(asin(sqrt(1.0/x)),2);
    }
}

gslpp::complex HiggsSigStr::g_func(const double x) const {
    if(x<1) {
    gslpp::complex z = -gslpp::complex::i()*M_PI;
    gslpp::complex gs1 = sqrt(1-x)*(log((1+sqrt(1-x))/(1-sqrt(1-x)))+z)/2.0;
    return gs1;
    }
    else {
        gslpp::complex gg1 = sqrt(x-1)*asin(sqrt(1.0/x));
        return gg1;
    }
}

gslpp::complex HiggsSigStr::Int1(const double tau, const double lambda) const {
    return tau*lambda/(tau-lambda)/2.0+tau*tau*lambda*lambda/((tau-lambda)
           *(tau-lambda))/2.0*(f_func(tau)-f_func(lambda))+tau*tau*lambda/((tau-lambda)
           *(tau-lambda))*(g_func(tau)-g_func(lambda));
    }

gslpp::complex HiggsSigStr::Int2(const double tau, const double lambda) const {
    return -tau*lambda/(tau-lambda)/2.0*(f_func(tau)-f_func(lambda));
    }
