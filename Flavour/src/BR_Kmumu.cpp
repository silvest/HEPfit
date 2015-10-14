/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "BR_Kmumu.h"

double BR_Kmumu::computeThValue()
{
    double theta= asin(sqrt( (M_PI * mySM.getAle() )/( sqrt(2) * mySM.getGF() * 
                   mySM.Mw_tree() * mySM.Mw_tree()) ));
    
    return((mySM.getMesons(QCD::K_0).getLifetime() / HCUT / mySM.getMesons(QCD::K_P).getLifetime() / HCUT)
           * mySM.getAle()*mySM.getAle()/(2.*M_PI*M_PI*pow(sin(theta),4.)) 
           * mySM.getBr_Kp_munu() * BRKmumu(NLO).real());
}

gslpp::complex BR_Kmumu::BRKmumu(orders order)
{
    if (mySM.getMyFlavour()->getHDS1().getCoeffDS1mumu().getOrder() < order){
        std::stringstream out;
        out << order;
        throw std::runtime_error("BRKmumu::computeThValue(): requires cofficient of "
                                 "order" + out.str() + "not computed");
    }
    
    gslpp::vector<gslpp::complex> ** allcoeff = mySM.getMyFlavour()->ComputeCoeffDS1mumu();
    
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