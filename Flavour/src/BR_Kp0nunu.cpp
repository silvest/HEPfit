/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "BR_Kp0nunu.h"

double BR_Kp0nunu::computeThValue()
{
    double theta = asin(sqrt( (M_PI * SM.getAle() )/( sqrt(2) * SM.getGF() * 
                   SM.Mw_tree() * SM.Mw_tree()) ));
    
    return(SM.getIB_Kl() * (SM.getMesons(QCD::K_0).getLifetime() / HCUT / SM.getMesons(QCD::K_P).getLifetime() / HCUT)
           * 3. * SM.getAle() * SM.getAle() / (2.*M_PI*M_PI*pow(sin(theta),4.)) * SM.getBr_Kp_P0enu() *
           BRKp0nunu(NLO, NLO_ew).real());
}

gslpp::complex BR_Kp0nunu::BRKp0nunu(orders order, orders_ew order_ew)
{
    if (mySM.getMyFlavour()->getHDS1().getCoeffDS1pnunu().getOrder() < order){
        std::stringstream out;
        out << order;
        throw std::runtime_error("BRKp0nunu::computeThValue(): requires cofficient of "
                                 "order" + out.str() + "not computed");
    }
    
    gslpp::vector<gslpp::complex> ** allcoeff = mySM.getMyFlavour()->ComputeCoeffDS1pnunu();
    
    switch(order_ew) {
        case NLO_ew:
            return((*(allcoeff[LO]) + *(allcoeff[NLO]) + *(allcoeff[NLO_ew])) *
                   (*(allcoeff[LO]) + *(allcoeff[NLO]) + *(allcoeff[NLO_ew])));
        case LO_ew:
            switch(order) {
                case NLO:
                    return((*(allcoeff[LO]) + *(allcoeff[NLO])) * 
                           (*(allcoeff[LO]) + *(allcoeff[NLO])));
                case LO:
                    return((*(allcoeff[LO])) * (*(allcoeff[LO])));
                default:
                    std::stringstream out;
                    out << order;
                    throw std::runtime_error("BRKp0nunu::BRKp0nunu(): order " + out.str() + "not implemented");
            }
         default:
            std::stringstream out;
            out << order_ew;
            throw std::runtime_error("BRKp0nunu::BRKp0nunu(): order_ew " + out.str() + "not implemented");
    }
}


