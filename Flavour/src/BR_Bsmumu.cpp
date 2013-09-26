/*
 * Copyright (C) 2012 SsuyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "BR_Bsmumu.h"

double BR_Bsmumu::getThValue(){
    double mmu = myFlavour.getModel().getLeptons(StandardModel::MU).getMass();
    double mBs = myFlavour.getModel().getMesons(QCD::B_S).getMass();
    double FBs = myFlavour.getModel().getMesons(QCD::B_S).getDecayconst();
    double beta = sqrt(1. - pow(2. * mmu / mBs, 2.));
    double coupling = myFlavour.getModel().getGF() * myFlavour.getModel().getAle() / 4. / M_PI / myFlavour.getModel().sW2();
    double PRF = pow(coupling, 2.) / M_PI / myFlavour.getModel().getMesons(QCD::B_S).getWidth() * pow(FBs, 2.) * pow(mmu, 2.) * mBs * beta;
    /*    double theta = asin(sqrt( (M_PI * myFlavour.getModel().getAle() )/( sqrt(2) * myFlavour.getModel().getGF() *
     myFlavour.getModel().Mw_tree() * myFlavour.getModel().Mw_tree()) ));
     
     return( myFlavour.getModel().getMesons(QCD::B_S).Lifetime() * myFlavour.getModel().getGF()*myFlavour.getModel().getGF()/M_PI
     * myFlavour.getModel().getAle()*myFlavour.getModel().getAle()/(16.*M_PI*M_PI*pow(sin(theta),4.))
     * mBs * myFlavour.getModel().getMesons(QCD::B_S).getDecayconst()
     * mmu * mmu * sqrt(1.-4.*mmu*mmu/mBs/mBs) * BRBsmumu(NLO).real());*/
    double Bsmumu = PRF * BRBsmumu(FULLNLO).real();
    
    return(Bsmumu);
}

complex BR_Bsmumu::BRBsmumu(orders order){
    if (myFlavour.getHDB1().getCoeffsmumu().getOrder() < order % 3){
        std::stringstream out;
        out << order;
        throw std::runtime_error("BRBsmumu::getThValue(): required cofficient of "
                                 "order " + out.str() + " not computed");
    }
    
    vector<complex> ** allcoeff = myFlavour.ComputeCoeffsmumu();
    
    switch(order) {
        case FULLNLO:
            return((*(allcoeff[LO]) + *(allcoeff[NLO])).conjugate() *
                   (*(allcoeff[LO]) + *(allcoeff[NLO])));
        case LO:
            return((*(allcoeff[LO])).conjugate() *
                   (*(allcoeff[LO])));
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("BRBsmumu::BRBsmumu(): order " + out.str() + "not implemented");;
    }
}