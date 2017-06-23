/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "BR_Kppnunu.h"
#include "StandardModel.h"
#include "std_make_vector.h"
#include "HeffDS1.h"

BR_Kppnunu::BR_Kppnunu(StandardModel& SM_i)
: ThObservable(SM_i), mySM(SM_i), CKpnunu(SM)
{
    setParametersForObservable(make_vector<std::string>() << "Br_Kp_P0enu" << "DeltaP_cu" << "IB_Kp");
}

double BR_Kppnunu::computeThValue()
{
    double theta= asin(sqrt( (M_PI * SM.getAle() )/( sqrt(2) * SM.getGF() * 
                   SM.Mw_tree() * SM.Mw_tree()) ));
    
    return( SM.getOptionalParameter("IB_Kp") * 3.*SM.getAle()*SM.getAle()/(2.*M_PI*M_PI*pow(sin(theta),4.))
           * SM.getOptionalParameter("Br_Kp_P0enu") * BRKppnunu(NLO, NLO_QED).real());
}

gslpp::complex BR_Kppnunu::BRKppnunu(orders order, orders_qed order_qed)
{
    if (mySM.getFlavour().getHDS1().getCoeffDS1pnunu().getOrder() < order){
        std::stringstream out;
        out << order;
        throw std::runtime_error("BRKppnunu::computeThValue(): requires cofficient of "
                                 "order" + out.str() + "not computed");
    }
    
    gslpp::vector<gslpp::complex> ** allcoeff = mySM.getFlavour().ComputeCoeffDS1pnunu();
    
    switch(order_qed) {
        case NLO_QED:
            return((*(allcoeff[LO]) + *(allcoeff[NLO]) + *(allcoeff[NLO_QED])) *
                   (*(allcoeff[LO]) + *(allcoeff[NLO]) + *(allcoeff[NLO_QED])) +
                   CKpnunu.C_TOT(NNLO,NLO_QED)*CKpnunu.C_TOT(NNLO,NLO_QED));
        case LO_QED:
            switch(order) {
                case NLO:
                    return((*(allcoeff[LO]) + *(allcoeff[NLO])) *
                           (*(allcoeff[LO]) + *(allcoeff[NLO])) +
                           CKpnunu.C_TOT(NLO,LO_QED)*CKpnunu.C_TOT(NLO,LO_QED));
                case LO:
                    return((*(allcoeff[LO])) * (*(allcoeff[LO]) ) +
                           CKpnunu.C_TOT(LO,LO_QED)*CKpnunu.C_TOT(LO,LO_QED));
                default:
                    std::stringstream out;
                    out << order;
                    throw std::runtime_error("BRKppnunu::BRKppnunu(): order " + out.str() + "not implemented");
            }
        default:
            std::stringstream out;
            out << order_qed;
            throw std::runtime_error("BRKppnunu::BRKppnunu(): order_qed " + out.str() + "not implemented");
    }
}


