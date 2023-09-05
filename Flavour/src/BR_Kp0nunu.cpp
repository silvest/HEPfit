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

BR_Kp0nunu::BR_Kp0nunu(const StandardModel& SM_i)
: ThObservable(SM_i), mySM(SM_i), CKpnunu(SM_i)
{
    
    mySM.initializeMeson(QCD::K_0);
    mySM.initializeMeson(QCD::K_P);
    
    setParametersForObservable(make_vector<std::string>() << "Br_Kp_P0enu"  << "IB_Kl" << "tKl" << "EpsK" << "DeltaP_cu"); 
}

double BR_Kp0nunu::computeThValue()
{
    double K0 = k_zero();
    
    return( K0 * r_epsK() * BRKp0nunu(NLO, NLO_QED11) ); 
} 

double BR_Kp0nunu::BRKp0nunu(orders order, orders_qed order_qed)
{
    gslpp::vector<gslpp::complex> ** allcoeff = mySM.getFlavour().ComputeCoeffDS1pnunu();
    
    if (mySM.getFlavour().getHDS1().getCoeffDS1pnunu().getOrder() < order){
        std::stringstream out;
        out << order;
        throw std::runtime_error("BRKppnunu::computeThValue(): requires cofficient of "
                                 "order" + out.str() + "not computed");
    }
    
    if (mySM.getFlavour().getHDS1().getCoeffDS1pnunu().getOrder_qed() < order_qed){ // discrepancy between top and charm contribution
        std::stringstream out2;
        out2 << order_qed;
        throw std::runtime_error("BRKppnunu::computeThValue(): requires cofficient of "
                                 " qed order " + out2.str() + " not computed");
    }
    
    switch(order_qed) {
        case NLO_QED11:
            return( ((*(allcoeff[LO]) + *(allcoeff[NLO]) + *(allcoeff[NLO_QED11])) *
                    (*(allcoeff[LO]) + *(allcoeff[NLO]) + *(allcoeff[NLO_QED11]))   ).real()  );
        case LO_QED:
            switch(order) {
                case NLO:
                    return(((*(allcoeff[LO]) + *(allcoeff[NLO]) ) *
                            (*(allcoeff[LO]) + *(allcoeff[NLO]) )  ).real());
                case LO:
                    return( ((*(allcoeff[LO])) * (*(allcoeff[LO]))).real());
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


double BR_Kp0nunu::r_epsK()
{   
    
    double Delta_Pcu=mySM.getOptionalParameter("DeltaP_cu");
    //double ReXi_o_ImXi= - (1. + (CKpnunu.P_C(NNLO,NLO_QED21)+Delta_Pcu) / Xtop / mySM.getCKM().getA() / mySM.getCKM().getA() - mySM.getCKM().getRho() ) / mySM.getCKM().getEta() ;
    gslpp::complex Xi = (mySM.getCKM().computelamc()*pow(mySM.getCKM().getLambda(),4.)*(CKpnunu.P_C(NNLO,NLO_QED21)+Delta_Pcu)+mySM.getCKM().computelamt()*mySM.getMatching().Xt(NLO,NLO_QED11));
    double epsk_value = mySM.getOptionalParameter("EpsK");
    
    return ( 1. + sqrt(2.) * (Xi.real()/Xi.imag()) * epsk_value) ;
}

double BR_Kp0nunu::k_zero()
{   
    double sw4 = mySM.sW2_MSbar_Approx() * mySM.sW2_MSbar_Approx();
    double lambda8 = pow(mySM.getCKM().getLambda(),8.);
    double IB_KL=mySM.getOptionalParameter("IB_Kl"); //arxiv: 0603079v3
    double BR_Kp_P0nunu=mySM.getOptionalParameter("Br_Kp_P0enu");//BR_Kp_P0nunu= 0.0507 ; // mySM.getOptionalParameter("Br_Kp_P0enu")
    double KL_lifetime = mySM.getOptionalParameter("tKl") ; //mySM.getMesons(QCD::K_0).getLifetime()
    double KP_lifetime =mySM.getMesons(QCD::K_P).getLifetime() ;//KP_lifetime = 1.238e-8 ; //mySM.getMesons(QCD::K_P).getLifetime()
    
    return (IB_KL * KL_lifetime / KP_lifetime * //rk+
           3. * mySM.Ale(mySM.getMz(),FULLNLO) * mySM.Ale(mySM.getMz(),FULLNLO) / (2. * M_PI * M_PI * sw4) * BR_Kp_P0nunu * lambda8 );
}

