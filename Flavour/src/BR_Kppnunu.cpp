/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "BR_Kppnunu.h"

double BR_Kppnunu::getThValue(){
    double theta= asin(sqrt( (M_PI * SM.getAle() )/( sqrt(2) * SM.getGF() * 
                   SM.Mw_tree() * SM.Mw_tree()) ));
    
    return( SM.getIB_Kp() * 3.*SM.getAle()*SM.getAle()/(2.*M_PI*M_PI*pow(sin(theta),4.))
           * SM.getBr_Kp_P0enu() * BRKppnunu(NLO, NLO_ew).real());
}

complex BR_Kppnunu::BRKppnunu(orders order, orders_ew order_ew){
    if (myFlavour.getHDS1().getCoeffDS1ppnunu().getOrder() < order){
        std::stringstream out;
        out << order;
        throw std::runtime_error("BRKppnunu::getThValue(): requires cofficient of "
                                 "order" + out.str() + "not computed");
    }
    
    vector<complex> ** allcoeff = myFlavour.ComputeCoeffDS1ppnunu();
    
    switch(order) {
        case NLO:
            switch(order_ew) {
                case NLO_ew:
                    return((*(allcoeff[LO]) + *(allcoeff[NLO]) + *(allcoeff[NLO_ew])) *
                           (*(allcoeff[LO]) + *(allcoeff[NLO]) + *(allcoeff[NLO_ew]))
                           + CPB.CT_tot(NNLO, NLO_ew) * CPB.CT_tot(NNLO, NLO_ew) );
                case LO_ew:
                    return((*(allcoeff[LO]) + *(allcoeff[NLO])) *
                           (*(allcoeff[LO]) + *(allcoeff[NLO]))
                           + CPB.CT_tot(NNLO, LO_ew) * CPB.CT_tot(NNLO, LO_ew));
                default:
                    std::stringstream out;
                    out << order_ew;
                    throw std::runtime_error("BRKppnunu::BRKppnunu(): order_ew " + out.str() + "not implemented");;
            }
        case LO:
            switch(order_ew) {
                case NLO_ew:
                    return((*(allcoeff[LO]) + *(allcoeff[NLO_ew])) *
                           (*(allcoeff[LO]) + *(allcoeff[NLO_ew]))
                           + CPB.CT_tot(NNLO, NLO_ew) *CPB.CT_tot(NNLO, NLO_ew));
                case LO_ew:
                    return((*(allcoeff[LO]) * (*(allcoeff[LO])) 
                           + CPB.CT_tot(NNLO, LO_ew)) * CPB.CT_tot(NNLO, LO_ew));
                default:
                    std::stringstream out;
                    out << order_ew;
                    throw std::runtime_error("BRKppnunu::BRKppnunu(): order_ew " + out.str() + "not implemented");;
            }
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("BRKppnunu::BRKppnunu(): order " + out.str() + "not implemented");;
    }
}

