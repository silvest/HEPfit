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
: ThObservable(SM_i), mySM(SM_i)
{
    
    mySM.initializeMeson(QCD::K_0);
    
    setParametersForObservable(make_vector<std::string>() << "PhSp_KL"  << "IB_Kl" << "IB_Kp" << "IB_K0p"  << "EpsK" << "DeltaP_cu" << "Vus_fpK0Pip"); 
}

double BR_Kp0nunu::computeThValue()
{    
    return( k_zero()  * BRKp0nunu(NLO, NLO_QED11) ); 
} 

double BR_Kp0nunu::BRKp0nunu(orders order, orders_qed order_qed)
{
    gslpp::vector<gslpp::complex>** c0_lamt_Xt = mySM.getFlavour().ComputeCoeffDS1pnunu(); 
    gslpp::vector<gslpp::complex> ** c0_lamc_Xc = mySM.getFlavour().ComputeCoeffDS1pnunuC();
    gslpp::complex eps(mySM.getOptionalParameter("EpsK")/sqrt(2.),mySM.getOptionalParameter("EpsK")/sqrt(2.));
    
     if (mySM.getFlavour().getHDS1().getCoeffDS1pnunu().getOrder() < order){
        std::stringstream out;
        out << order;
        throw std::runtime_error("BRKppnunu::computeThValue(): requires cofficient of "
                                 "order" + out.str() + "not computed");
    }
    
    if (mySM.getFlavour().getHDS1().getCoeffDS1pnunu().getOrder_qed()  < order_qed){ 
        std::stringstream out2;
        out2 << order_qed;
        throw std::runtime_error("BRKppnunu::computeThValue(): requires cofficient of "
                                 " qed order " + out2.str() + " not computed");
    }
    gslpp::complex c0_Xi;
    double br = 0.;
    switch(order_qed) {
        case NLO_QED11:
            c0_Xi =(  (*(c0_lamt_Xt[LO]))(0) + (*(c0_lamt_Xt[NLO]))(0) + (*(c0_lamt_Xt[NLO_QED11]))(0) +
                                     (*(c0_lamc_Xc[LO]))(0) + (*(c0_lamc_Xc[NLO]))(0) + (*(c0_lamc_Xc[NNLO]))(0) + (*(c0_lamc_Xc[LO_QED]))(0) + (*(c0_lamc_Xc[NLO_QED11]))(0) + LongDistance() );
            br = ( 0.5*(c0_Xi - c0_Xi.conjugate()*(1-eps)/(1+eps)).abs2() / (1+((1-eps)/(1+eps)).abs2()) );
            c0_Xi =(  (*(c0_lamt_Xt[LO]))(1) + (*(c0_lamt_Xt[NLO]))(1) + (*(c0_lamt_Xt[NLO_QED11]))(1) );
            br += ( 0.5*(c0_Xi - c0_Xi.conjugate()*(1-eps)/(1+eps)).abs2() / (1+((1-eps)/(1+eps)).abs2()) );
            return br;
        case LO_QED:{
            switch(order) {
                case NNLO:{
                    c0_Xi = (*(c0_lamt_Xt[LO]))(0) + (*(c0_lamt_Xt[NLO]))(0)  +
                                           (*(c0_lamc_Xc[LO]))(0) + (*(c0_lamc_Xc[NLO]))(0) + (*(c0_lamc_Xc[NNLO]))(0) + LongDistance()  ;
                    br = ( 0.5*(c0_Xi - c0_Xi.conjugate()*(1-eps)/(1+eps)).abs2() / (1+((1-eps)/(1+eps)).abs2()) );
                    c0_Xi = (*(c0_lamt_Xt[LO]))(1) + (*(c0_lamt_Xt[NLO]))(1)  ;
                    br += ( 0.5*(c0_Xi - c0_Xi.conjugate()*(1-eps)/(1+eps)).abs2() / (1+((1-eps)/(1+eps)).abs2()) );
                return br;
                }
                case NLO:{
                    c0_Xi = (*(c0_lamt_Xt[LO]))(0) + (*(c0_lamt_Xt[NLO]))(0)  +
                                           (*(c0_lamc_Xc[LO]))(0) + (*(c0_lamc_Xc[NLO]))(0) + LongDistance()  ;
                    br = ( 0.5*(c0_Xi - c0_Xi.conjugate()*(1-eps)/(1+eps)).abs2() / (1+((1-eps)/(1+eps)).abs2()) );
                    c0_Xi = (*(c0_lamt_Xt[LO]))(1) + (*(c0_lamt_Xt[NLO]))(1)  ;
                    br += ( 0.5*(c0_Xi - c0_Xi.conjugate()*(1-eps)/(1+eps)).abs2() / (1+((1-eps)/(1+eps)).abs2()) );
            return br;
                }
                case LO:{
                    c0_Xi = (*(c0_lamt_Xt[LO]))(0) +
                                           (*(c0_lamc_Xc[LO]))(0) + LongDistance()  ;
                    br = ( 0.5*(c0_Xi - c0_Xi.conjugate()*(1-eps)/(1+eps)).abs2() / (1+((1-eps)/(1+eps)).abs2()) );
                    c0_Xi = (*(c0_lamt_Xt[LO]))(1) ;
                    br += ( 0.5*(c0_Xi - c0_Xi.conjugate()*(1-eps)/(1+eps)).abs2() / (1+((1-eps)/(1+eps)).abs2()) );
            return br;
                }
                default:
                    std::stringstream out;
                    out << order;
                    throw std::runtime_error("BRKp0nunu::BRKp0nunu(): order " + out.str() + "not implemented");
            }
        }
        default:
            std::stringstream out;
            out << order_qed;
            throw std::runtime_error("BRKp0nunu::BRKp0nunu(): order_qed " + out.str() + "not implemented");
    }   
}


double BR_Kp0nunu::k_zero()
{   
    return (mySM.getOptionalParameter("IB_Kp") * mySM.getOptionalParameter("IB_Kp") * mySM.getOptionalParameter("IB_Kl") * mySM.getOptionalParameter("IB_Kl") / mySM.getOptionalParameter("IB_K0p") / mySM.getOptionalParameter("IB_K0p") *
            mySM.getOptionalParameter("Vus_fpK0Pip") * mySM.getOptionalParameter("Vus_fpK0Pip") * mySM.getOptionalParameter("PhSp_KL") *
            pow(mySM.getMesons(QCD::K_0).getMass(),5.) * mySM.getMesons(QCD::K_0).getLifetime() / HCUT / 512 / M_PI / M_PI / M_PI / mySM.getCKM().getLambda() / mySM.getCKM().getLambda());
    
}

gslpp::complex BR_Kp0nunu::LongDistance() // c0 * lambda_c * delta_Pc * lambda^4
{   
    return 4. * mySM.getGF() / sqrt(2.) * mySM.alphaMz() / 2. / M_PI / mySM.sW2_ND() * mySM.getCKM().computelamc() * mySM.getOptionalParameter("DeltaP_cu") * pow(mySM.getCKM().getLambda(),4.)  ;
}

