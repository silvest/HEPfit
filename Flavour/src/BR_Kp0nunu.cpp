/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "BR_Kp0nunu.h"

double BR_Kp0nunu::getThValue(){
    double theta = asin(sqrt( (M_PI * SM.getAle() )/( sqrt(2) * SM.getGF() * 
                   SM.Mw_tree() * SM.Mw_tree()) ));
    
    return(SM.getIB_Kl() * (SM.getMesons(QCD::K_0).Lifetime()/SM.getMesons(QCD::K_P).Lifetime())
           * 3. * SM.getAle() * SM.getAle() / (2.*M_PI*M_PI*pow(sin(theta),4.)) * SM.getBr_Kp_ppenu() *
           BRKp0nunu(NLO, NLO_ew).real());
}

complex BR_Kp0nunu::BRKp0nunu(orders order, orders_ew order_ew){
    if (myFlavour.getHDS1().getCoeffDS1p0nunu().getOrder() < order){
        std::stringstream out;
        out << order;
        throw std::runtime_error("BRKp0nunu::getThValue(): requires cofficient of "
                                 "order" + out.str() + "not computed");
    }
    
    vector<complex> ** allcoeff = myFlavour.ComputeCoeffDS1p0nunu();
    
    switch(order_ew) {
        case NLO_ew:
            return((*(allcoeff[LO]) + *(allcoeff[NLO]) + *(allcoeff[NLO_ew])) *
                   (*(allcoeff[LO]) + *(allcoeff[NLO]) + *(allcoeff[NLO_ew])));
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("BRKp0nunu::BRKp0nunu(): order " + out.str() + "not implemented");;
    }
    
    switch(order) {
        case NLO:
            return((*(allcoeff[LO]) + *(allcoeff[NLO])) * 
                   (*(allcoeff[LO]) + *(allcoeff[NLO])));
        case LO:
            return((*(allcoeff[LO])) * (*(allcoeff[LO])));
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("BRKp0nunu::BRKp0nunu(): order " + out.str() + "not implemented");;
    }
}


