/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "THDMMatching.h"
#include "THDM.h"
#include <math.h>
#include <stdexcept>

THDMMatching::THDMMatching(const THDM & THDM_i) :

    StandardModelMatching(THDM_i),
    myCKM(3, 3, 0.),
    myTHDM(THDM_i),
    mcdbs2(5, NDR, NLO),
    mcbtaunu(3, NDR, LO)
{}

std::vector<WilsonCoefficient>& THDMMatching::CMdbs2() {

    double Mut = myTHDM.getMut();
    double xt = x_t(Mut);
    double GF=myTHDM.getGF();
    double MW=myTHDM.Mw();
    gslpp::complex co = GF / 4. / M_PI * MW * myTHDM.computelamt_s();
    double tanb = myTHDM.gettanb();
    double mHp2=myTHDM.getmHp2();
    double xHW=mHp2/(MW*MW);
    double xtH=xt/xHW;
    double SWH=xtH*((2.0*xHW-8.0)*log(xtH)/((1.0-xHW)*(1.0-xtH)*(1.0-xtH))+6.0*xHW*log(xt)/((1.0-xHW)*(1.0-xt)*(1.0-xt))-(8.0-2.0*xt)/((1.0-xt)*(1.0-xtH)))/(tanb*tanb);
    double SHH=xtH*((1.0+xtH)/((1.0-xtH)*(1.0-xtH))+2.0*xtH*log(xtH)/((1.0-xtH)*(1.0-xtH)*(1.0-xtH)))/(tanb*tanb*tanb*tanb);

    vmcds = StandardModelMatching::CMdbs2();
    mcdbs2.setMu(Mut);

    switch (mcdbs2.getOrder()) {
        case NNLO:
        case NLO:
        case LO:
            mcdbs2.setCoeff(0, co * co * xt * (SWH+SHH), LO);
            break;
        default:
            std::stringstream out;
            out << mcdbs2.getOrder();
            throw std::runtime_error("THDMMatching::CMdbs2(): order " + out.str() + "not implemented");
    }

    vmcds.push_back(mcdbs2);
    return(vmcds);
}

std::vector<WilsonCoefficient>& THDMMatching::CMbtaunu() {

    double Muw = myTHDM.getMuw();
    double GF = myTHDM.getGF();
    myCKM = myTHDM.getVCKM();
    double mB = myTHDM.getMesons(QCD::B_P).getMass();
    double tanb = myTHDM.gettanb();
    double mHp2=myTHDM.getmHp2();

    vmcbtaunu = StandardModelMatching::CMbtaunu();
    mcbtaunu.setMu(Muw);
 
    switch (mcbtaunu.getOrder()) {
        case NNLO:
        case NLO:
        case LO:
            mcbtaunu.setCoeff(0, -4.*GF * myCKM(0,2) / sqrt(2.) * mB*mB*tanb*tanb/mHp2, LO);
            break;
        default:
            std::stringstream out;
            out << mcbtaunu.getOrder();
            throw std::runtime_error("THDMMatching::CMbtaunu(): order " + out.str() + "not implemented");
    }

    vmcbtaunu.push_back(mcbtaunu);
    return(vmcbtaunu);

}
