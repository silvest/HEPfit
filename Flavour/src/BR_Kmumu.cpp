/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "BR_Kmumu.h"
#include "StandardModel.h"
#include "std_make_vector.h"

BR_Kmumu::BR_Kmumu(StandardModel& SM_i)
: ThObservable(SM_i), mySM(SM_i), CPB(SM)
{
    setParametersForObservable(make_vector<std::string>() << "Br_Kp_munu" << "DeltaP_cu");
    mySM.initializeMeson(QCD::K_0);
    mySM.initializeMeson(QCD::K_P);
}

double BR_Kmumu::computeThValue()
{
    double theta= asin(sqrt( (M_PI * mySM.getAle() )/( sqrt(2) * mySM.getGF() * 
                   mySM.Mw_tree() * mySM.Mw_tree()) ));
    
    return((mySM.getMesons(QCD::K_0).getLifetime() / HCUT / mySM.getMesons(QCD::K_P).getLifetime() / HCUT)
           * mySM.getAle()*mySM.getAle()/(2.*M_PI*M_PI*pow(sin(theta),4.)) 
           * mySM.getOptionalParameter("Br_Kp_munu") * BRKmumu(NLO).real());
}

gslpp::complex BR_Kmumu::BRKmumu(orders order)
{
    if (mySM.getFlavour().getHDS1().getCoeffDS1mumu().getOrder() < order){
        std::stringstream out;
        out << order;
        throw std::runtime_error("BRKmumu::computeThValue(): requires cofficient of "
                                 "order" + out.str() + "not computed");
    }
    
    gslpp::vector<gslpp::complex> ** allcoeff = mySM.getFlavour().ComputeCoeffDS1mumu();
    
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