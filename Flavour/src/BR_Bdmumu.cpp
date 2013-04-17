/* 
 * Copyright (C) 2012 SsuyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "BR_Bdmumu.h"

double BR_Bdmumu::getThValue(){
    double mmu = myFlavour.getModel().getLeptons(StandardModel::MU).getMass();
    double mbd = myFlavour.getModel().getMesons(QCD::B_D).getMass();
    double theta = asin(sqrt( (M_PI * myFlavour.getModel().getAle() )/( sqrt(2) * myFlavour.getModel().getGF() * 
                   myFlavour.getModel().Mw_tree() * myFlavour.getModel().Mw_tree()) ));
    
    return( myFlavour.getModel().getMesons(QCD::B_D).Lifetime() * myFlavour.getModel().getGF()*myFlavour.getModel().getGF()/M_PI
            * myFlavour.getModel().getAle()*myFlavour.getModel().getAle()/(16.*M_PI*M_PI*pow(sin(theta),4.)) 
            * mbd * myFlavour.getModel().getMesons(QCD::B_D).getDecayconst() 
            * mmu * mmu * sqrt(1.-4.*mmu*mmu/mbd/mbd) * BRBdmumu(NLO).real());
}

complex BR_Bdmumu::BRBdmumu(orders order){
    if (myFlavour.getHDB1().getCoeffsmumu().getOrder() < order){
        std::stringstream out;
        out << order;
        throw std::runtime_error("BRBdmumu::getThValue(): requires cofficient of "
                                 "order" + out.str() + "not computed");
    }
    
    vector<complex> ** allcoeff = myFlavour.ComputeCoeffdmumu();
    
    switch(order) {
        case NLO:
                return((*(allcoeff[LO]) + *(allcoeff[NLO])) *
                       (*(allcoeff[LO]) + *(allcoeff[NLO])));
        case LO:
                return((*(allcoeff[LO])) *
                       (*(allcoeff[LO])));
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("BRBdmumu::BRBdmumu(): order " + out.str() + "not implemented");;
    }
}