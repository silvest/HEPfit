/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "BR_Kmumu.h"

double BR_Kmumu::getThValue(){
    double theta= asin(sqrt( (M_PI * myFlavour.getModel().getAle() )/( sqrt(2) * myFlavour.getModel().getGF() * 
                   myFlavour.getModel().Mw_tree() * myFlavour.getModel().Mw_tree()) ));
    
    return((myFlavour.getModel().getMesons(QCD::K_0).getLifetime() / HCUT / myFlavour.getModel().getMesons(QCD::K_P).getLifetime() / HCUT)
           * myFlavour.getModel().getAle()*myFlavour.getModel().getAle()/(2.*M_PI*M_PI*pow(sin(theta),4.)) 
           * myFlavour.getModel().getBr_Kp_munu() * BRKmumu(NLO).real());
}

complex BR_Kmumu::BRKmumu(orders order){
    if (myFlavour.getHDS1().getCoeffDS1mumu().getOrder() < order){
        std::stringstream out;
        out << order;
        throw std::runtime_error("BRKmumu::getThValue(): requires cofficient of "
                                 "order" + out.str() + "not computed");
    }
    
    vector<complex> ** allcoeff = myFlavour.ComputeCoeffDS1mumu();
    
    switch(order) {
        case NLO:
                return((*(allcoeff[LO]) + *(allcoeff[NLO])) *
                       (*(allcoeff[LO]) + *(allcoeff[NLO]))
                       + CPB.X_ch() );
        case LO:
                return((*(allcoeff[LO])) *
                       (*(allcoeff[LO]))
                       + CPB.X_ch() );
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("BRKmumu::BRKmumu(): order " + out.str() + "not implemented");;
    }
}