/*
 * Copyright (C) 2023 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Dleptonnu.h"
#include "StandardModel.h"


Dleptonnu::Dleptonnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::lepton lepton_i)
: ThObservable(SM_i)
{
    lepton = lepton_i;
    meson = meson_i;
    SM.initializeMeson(meson);
    
};


double Dleptonnu::computeThValue()
{
    //The WC are written in the LEFT basis of arxiv:1709.04486 the expressions can be found in arxiv:1706.00410 and arxiv:1605.07114 in a similar basis
    gslpp::vector<gslpp::complex> ** allcoeff = SM.getFlavour().ComputeCoeffcleptonnu(meson, lepton);


    double mD = SM.getMesons(meson).getMass();
    double mc = SM.getQuarks(QCD::CHARM).getMass();


    double mlight=0;
    double mlepton=0;
    double fact = 1.; /*factor introduced to scale the decay constant from that of the neutral B to the charged B.*/
    //double fact = 0.989;
    
    switch (lepton){
        case QCD::TAU:      
            mlepton = SM.getLeptons(StandardModel::TAU).getMass();
        break;
        case QCD::MU:      
            mlepton = SM.getLeptons(StandardModel::MU).getMass();
        break;
        default:
            std::runtime_error("The observable Dleptonnu is not added for that lepton " );
    }
    
    switch(meson) {
        case QCD::D_P:
            mlight=SM.getQuarks(QCD::DOWN).getMass();
        break;
        case QCD::D_S:
            mlight=SM.getQuarks(QCD::STRANGE).getMass();
        break;
        default:
            throw std::runtime_error("The observable Dleptonnu is only implemented for QCD::D_S " );
}
    
    //std::cout<< "SM value for D_S = "<< 1. / (64. * M_PI) * mlepton * mlepton * pow(fact * SM.getMesons(meson).getDecayconst(), 2.) * mD * pow(1. - mlepton * mlepton / mD / mD, 2.) / SM.getMesons(meson).computeWidth() * (4.*SM.getGF() * SM.getCKM().getV_cs() / sqrt(2.)).abs2() <<std::endl;
    
    return 1. / (64. * M_PI) * mlepton * mlepton * pow(fact * SM.getMesons(meson).getDecayconst(), 2.) * mD * pow(1. - mlepton * mlepton / mD / mD, 2.) / SM.getMesons(meson).computeWidth() * ((*(allcoeff[LO]))(0) 
            -(*(allcoeff[LO]))(1) + mD * mD / (mc+mlight) / mlepton * ((*(allcoeff[LO]))(2) - (*(allcoeff[LO]))(3))).abs2(); // PLEASE NOTE THE DECAY CONST
    

}