/*
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Btaunu.h"
#include "StandardModel.h"

Btaunu::Btaunu(const StandardModel& SM_i, QCD::meson meson_i)
: ThObservable(SM_i)
{
    meson = meson_i;
    SM.initializeMeson(meson);
};


double Btaunu::computeThValue()
{
    gslpp::vector<gslpp::complex> ** allcoeff;
    double mlight;
    
    if(meson == QCD::B_P) {
        mlight=SM.getQuarks(QCD::UP).getMass();
        allcoeff = SM.getFlavour().ComputeCoeffdiujlknu(2,0,2,SM.getQuarks(QCD::BOTTOM).getMass());
    }
    else if (meson == QCD::B_C) { 
        mlight=SM.getQuarks(QCD::CHARM).getMass();
        allcoeff = SM.getFlavour().ComputeCoeffdiujlknu(2,1,2,SM.getQuarks(QCD::BOTTOM).getMass());
    }
    else 
        throw std::runtime_error("The observable Btaunu is only implemented for QCD::B_P and QCD::B_C " );
    //The WC are written in the LEFT basis of arxiv:1709.04486 the expressions can be found in arxiv:1706.00410 and arxiv:1605.07114 in a similar basis

    double mtau = SM.getLeptons(StandardModel::TAU).getMass();
    double mB = SM.getMesons(meson).getMass();
    double mb = SM.getQuarks(QCD::BOTTOM).getMass();

    double fact = 1.; /*factor introduced to scale the decay constant from that of the neutral B to the charged B.*/
    //double fact = 0.989;
    
    return 1. / (64. * M_PI) * mtau * mtau * pow(fact * SM.getMesons(meson).getDecayconst(), 2.) * mB * (1. - mtau * mtau / mB / mB) *  (1. - mtau * mtau / mB / mB) / SM.getMesons(meson).computeWidth() * ((*(allcoeff[LO]))(0) 
            -(*(allcoeff[LO]))(1) + mB * mB / (mb+mlight) / mtau * ((*(allcoeff[LO]))(2) - (*(allcoeff[LO]))(3))).abs2(); // PLEASE NOTE THE DECAY CONST
    

}