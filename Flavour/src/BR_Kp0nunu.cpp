/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "BR_Kp0nunu.h"
#include "StandardModel.h"
#include "std_make_vector.h"
#include "HeffDS1.h"

BR_Kp0nunu::BR_Kp0nunu(StandardModel& SM_i)
: ThObservable(SM_i), mySM(SM_i)
{
    setParametersForObservable(make_vector<std::string>() << "Br_Kp_P0enu" << "IB_Kl");
    mySM.initializeMeson(QCD::K_0);
    mySM.initializeMeson(QCD::K_P);
}

double BR_Kp0nunu::computeThValue()
{
    double theta = asin(sqrt( (M_PI * SM.getAle() )/( sqrt(2) * SM.getGF() * 
                   SM.Mw_tree() * SM.Mw_tree()) ));
    
    return(SM.getOptionalParameter("IB_Kl") * (SM.getMesons(QCD::K_0).getLifetime() / HCUT / SM.getMesons(QCD::K_P).getLifetime() / HCUT)
           * 3. * SM.getAle() * SM.getAle() / (2.*M_PI*M_PI*pow(sin(theta),4.)) * SM.getOptionalParameter("Br_Kp_P0enu") *
           BRKp0nunu(NLO, NLO_QED).real());
}

gslpp::complex BR_Kp0nunu::BRKp0nunu(orders order, orders_qed order_qed)
{
    if (mySM.getFlavour().getHDS1().getCoeffDS1pnunu().getOrder() < order){
        std::stringstream out;
        out << order;
        throw std::runtime_error("BRKp0nunu::computeThValue(): requires cofficient of "
                                 "order" + out.str() + "not computed");
    }
    
    gslpp::vector<gslpp::complex> ** allcoeff = mySM.getFlavour().ComputeCoeffDS1pnunu();
    
    switch(order_qed) {
        case NLO_QED:
            return((*(allcoeff[LO]) + *(allcoeff[NLO]) + *(allcoeff[NLO_QED])) *
                   (*(allcoeff[LO]) + *(allcoeff[NLO]) + *(allcoeff[NLO_QED])));
        case LO_QED:
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
            out << order_qed;
            throw std::runtime_error("BRKp0nunu::BRKp0nunu(): order_qed " + out.str() + "not implemented");
    }
}


