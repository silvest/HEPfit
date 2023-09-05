/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *ee doc/COPYING.
 */

#include "BR_Kppnunu.h"
#include "StandardModel.h"
#include "std_make_vector.h"
#include "StandardModelMatching.h"
#include "HeffDS1.h"

BR_Kppnunu::BR_Kppnunu(const StandardModel& SM_i)
: ThObservable(SM_i), mySM(SM_i), CKpnunu(SM_i)
{   
    
    mySM.initializeMeson(QCD::K_P);
    
    setParametersForObservable(make_vector<std::string>() <<  "DeltaP_cu"  << "IB_Kp" << "fKplus" << "tKl" << "PhSp_KP" << "Delta_EM"); 
   
}

double BR_Kppnunu::computeThValue()
{
    double Delta_EM = mySM.getOptionalParameter("Delta_EM"); // see arxiv: 2105.02868v1 and 0705.2025v2

    return(  k_plus() * (1. + Delta_EM) * BRKppnunu(NLO, NLO_QED21) ); 
}

double BR_Kppnunu::BRKppnunu(orders order, orders_qed order_qed)
{   
    if (mySM.getFlavour().getHDS1().getCoeffDS1pnunu().getOrder() < order){
        std::stringstream out;
        out << order;
        throw std::runtime_error("BRKppnunu::computeThValue(): requires cofficient of "
                                 "order" + out.str() + "not computed");
    }
    
    if (mySM.getFlavour().getHDS1().getCoeffDS1pnunu().getOrder_qed() + 1 < order_qed){ // discrepancy between top and charm contribution
        std::stringstream out2;
        out2 << order_qed;
        throw std::runtime_error("BRKppnunu::computeThValue(): requires cofficient of "
                                 " qed order " + out2.str() + " not computed");
    }
    gslpp::vector<gslpp::complex> ** allcoeff = mySM.getFlavour().ComputeCoeffDS1pnunu();
    //Charm_Kpnunu CKpnunu = Charm_Kpnunu(SM);

    
    switch(order_qed) {
        case NLO_QED21:
           return ( ((*(allcoeff[LO]) + *(allcoeff[NLO]) + *(allcoeff[NLO_QED11])) *
                     (*(allcoeff[LO]) + *(allcoeff[NLO]) + *(allcoeff[NLO_QED11]))).real() + 
                     CKpnunu.C_TOT(NNLO,NLO_QED21)*CKpnunu.C_TOT(NNLO,NLO_QED21) );
            break;
        case NLO_QED11:
            return( ((*(allcoeff[LO]) + *(allcoeff[NLO]) + *(allcoeff[NLO_QED11])) *
                    (*(allcoeff[LO]) + *(allcoeff[NLO]) + *(allcoeff[NLO_QED11]))).real() +
                    CKpnunu.C_TOT(NNLO,NLO_QED11)*CKpnunu.C_TOT(NNLO,NLO_QED11)  );
            break;
        case LO_QED:
            switch(order) {
                case NLO:
                    return( ((*(allcoeff[LO]) + *(allcoeff[NLO])) *
                            (*(allcoeff[LO]) + *(allcoeff[NLO]))).real() +
                            CKpnunu.C_TOT(NLO,LO_QED)*CKpnunu.C_TOT(NLO,LO_QED));
                    break;
                case LO:
                    return(  ((*(allcoeff[LO])) * (*(allcoeff[LO]))).real() +
                             CKpnunu.C_TOT(LO,LO_QED)*CKpnunu.C_TOT(LO,LO_QED) );
                    break;
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



double BR_Kppnunu::k_plus(){
    

    double lambda8 = pow(mySM.getCKM().getV_us().abs(),8.); 
    double sw4 = mySM.sW2_MSbar_Approx() * SM.sW2_MSbar_Approx();
    
    double prefactor = mySM.getGF() * mySM.getGF() * pow(mySM.getMesons(QCD::K_P).getMass(),5.) * 
                       mySM.Ale(mySM.getMz(),FULLNLO) * mySM.Ale(mySM.getMz(),FULLNLO) / 256. / pow(M_PI,5.) / sw4;
    
    return (prefactor * lambda8 * mySM.getMesons(QCD::K_P).getLifetime() / HCUT  * 
            pow(mySM.getOptionalParameter("IB_Kp") * mySM.getOptionalParameter("fKplus") * mySM.getCKM().getLambda() , 2. ) * mySM.getOptionalParameter("PhSp_KP")) ;
    
    
}

