/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "NPBR_kpnunu.h"
#include "StandardModel.h"
#include "std_make_vector.h"
#include "StandardModelMatching.h"
#include "HeffDS1.h"

NPBR_kp0nunu::NPBR_kp0nunu(const StandardModel& SM_i)
: ThObservable(SM_i), mySM(SM_i), myBR_Kp0nunu(SM_i)

{
    
    setParametersForObservable(make_vector<std::string>()  << "PhSp_KL"  << "IB_Kl" << "IB_Kp" << "IB_K0p"  << "EpsK" << "DeltaP_cu" << "Vus_fpK0Pip" << "R_K" << "Theta_K"  );  
    
    mySM.initializeMeson(QCD::K_0);
    
}

double NPBR_kp0nunu::computeThValue()
{   
    return( myBR_Kp0nunu.k_zero() * BRKp0nunu_NP(NLO, NLO_QED11) ); 
}

double NPBR_kp0nunu::BRKp0nunu_NP(orders order, orders_qed order_qed)
{
    gslpp::vector<gslpp::complex>** c0_lamt_Xt = mySM.getFlavour().ComputeCoeffDS1pnunu(); 
    gslpp::vector<gslpp::complex> ** c0_lamc_Xc = mySM.getFlavour().ComputeCoeffDS1pnunuC();
    gslpp::complex eps(mySM.getOptionalParameter("EpsK")/sqrt(2.),mySM.getOptionalParameter("EpsK")/sqrt(2.));
    gslpp::complex NP_contribution(mySM.getOptionalParameter("R_K")*cos(mySM.getOptionalParameter("Theta_K")),-mySM.getOptionalParameter("R_K")*sin(mySM.getOptionalParameter("Theta_K")));
    
    
     if (mySM.getFlavour().getHDS1().getCoeffDS1pnunu().getOrder() < order){
        std::stringstream out;
        out << order;
        throw std::runtime_error("NPBR_kp0nunu::computeThValue(): requires cofficient of "
                                 "order" + out.str() + "not computed");
    }
    
    if (mySM.getFlavour().getHDS1().getCoeffDS1pnunu().getOrder_qed() < order_qed){ 
        std::stringstream out2;
        out2 << order_qed;
        throw std::runtime_error("NPBR_kp0nunu::computeThValue(): requires cofficient of "
                                 " qed order " + out2.str() + " not computed");
    }
    switch(order_qed) {
        case NLO_QED11:{
            gslpp::complex c0_Xi =( ( (*(c0_lamt_Xt[LO]))(0) + (*(c0_lamt_Xt[NLO]))(0) + (*(c0_lamt_Xt[NLO_QED11]))(0) ) * NP_contribution  +
                                     (*(c0_lamc_Xc[LO]))(0) + (*(c0_lamc_Xc[NLO]))(0) + (*(c0_lamc_Xc[NNLO]))(0) + (*(c0_lamc_Xc[LO_QED]))(0) + (*(c0_lamc_Xc[NLO_QED11]))(0) + myBR_Kp0nunu.LongDistance() );
            return ( 0.5*(c0_Xi - c0_Xi.conjugate()*(1-eps)/(1+eps)).abs2() / (1+((1-eps)/(1+eps)).abs2()) );
        }
        case LO_QED:{
            switch(order) {
                case NNLO:{
                    gslpp::complex c0_Xi =( ( (*(c0_lamt_Xt[LO]))(0) + (*(c0_lamt_Xt[NLO]))(0) ) * NP_contribution +
                                           (*(c0_lamc_Xc[LO]))(0) + (*(c0_lamc_Xc[NLO]))(0) + (*(c0_lamc_Xc[NNLO]))(0) + myBR_Kp0nunu.LongDistance() ) ;
                return ( 0.5*(c0_Xi - c0_Xi.conjugate()*(1-eps)/(1+eps)).abs2() / (1+((1-eps)/(1+eps)).abs2()) );
                }
                case NLO:{
                    gslpp::complex c0_Xi =( ( (*(c0_lamt_Xt[LO]))(0) + (*(c0_lamt_Xt[NLO]))(0) ) * NP_contribution +
                                           (*(c0_lamc_Xc[LO]))(0) + (*(c0_lamc_Xc[NLO]))(0) + myBR_Kp0nunu.LongDistance()  );
                return ( 0.5*(c0_Xi - c0_Xi.conjugate()*(1-eps)/(1+eps)).abs2() / (1+((1-eps)/(1+eps)).abs2()) );
                }
                case LO:{
                    gslpp::complex c0_Xi = ( (*(c0_lamt_Xt[LO]))(0) * NP_contribution +
                                           (*(c0_lamc_Xc[LO]))(0) + myBR_Kp0nunu.LongDistance() );
                return ( 0.5*(c0_Xi - c0_Xi.conjugate()*(1-eps)/(1+eps)).abs2() / (1+((1-eps)/(1+eps)).abs2()) );
                }
                default:
                    std::stringstream out;
                    out << order;
                    throw std::runtime_error("NPBR_kp0nunu::BRKp0nunu(): order " + out.str() + "not implemented");
            }
        }
        default:
            std::stringstream out;
            out << order_qed;
            throw std::runtime_error("NPBR_kp0nunu::BRKp0nunu(): order_qed " + out.str() + "not implemented");
    }   
}


// ############################################################################## //

NPBR_kppnunu::NPBR_kppnunu(const StandardModel& SM_i)
: ThObservable(SM_i), mySM(SM_i), myBR_Kppnunu(SM_i)
{
    setParametersForObservable(make_vector<std::string>() <<  "DeltaP_cu"  << "IB_Kp" << "Vus_fpK0Pip" << "PhSp_KP" << "Delta_EM"  << "R_K" << "Theta_K"); 
    
}

double NPBR_kppnunu::computeThValue()
{

    return(   myBR_Kppnunu.k_plus() * (1. + mySM.getOptionalParameter("Delta_EM")) * BRKppnunu_NP(NLO, NLO_QED11) ); 
}

double NPBR_kppnunu::BRKppnunu_NP(orders order, orders_qed order_qed)
{
    
     if (mySM.getFlavour().getHDS1().getCoeffDS1pnunu().getOrder() < order){
        std::stringstream out;
        out << order;
        throw std::runtime_error("NPBR_kppnunu::computeThValue(): requires cofficient of "
                                 "order" + out.str() + "not computed");
    }
    
    if (mySM.getFlavour().getHDS1().getCoeffDS1pnunu().getOrder_qed() < order_qed){ 
        std::stringstream out2;
        out2 << order_qed;
        throw std::runtime_error("NPBR_kppnunu::computeThValue(): requires cofficient of "
                                 " qed order " + out2.str() + " not computed");
    }
    
    gslpp::vector<gslpp::complex> ** c0_lamt_Xt = mySM.getFlavour().ComputeCoeffDS1pnunu(); 
    gslpp::vector<gslpp::complex> ** c0_lamc_Xc = mySM.getFlavour().ComputeCoeffDS1pnunuC();
    gslpp::complex NP_contribution(mySM.getOptionalParameter("R_K")*cos(mySM.getOptionalParameter("Theta_K")),-mySM.getOptionalParameter("R_K")*sin(mySM.getOptionalParameter("Theta_K")));
    
    switch(order_qed) {
        case NLO_QED11:
            //std::cout << "Pc: " << ( (*(c0_lamc_Xc[LO])).real() + (*(c0_lamc_Xc[NLO])).real() + (*(c0_lamc_Xc[NNLO])).real() + (*(c0_lamc_Xc[NLO_QED11])).real() + (*(c0_lamc_Xc[NLO_QED21])).real() ) / (4. * mySM.getGF() / sqrt(2.) * mySM.alphaMz() / 2. / M_PI / mySM.sW2_ND()) / pow(mySM.getCKM().getLambda(),4.) / mySM.getCKM().computelamc().real()  << std::endl;
            //std::cout << "deltaPC: " << myBR_Kppnunu.LongDistance().real() / (4. * mySM.getGF() / sqrt(2.) * mySM.alphaMz() / 2. / M_PI / mySM.sW2_ND()) / pow(mySM.getCKM().getLambda(),4.) / mySM.getCKM().computelamc().real() << std::endl;
            //std::cout << "Xt: " << ( ((*(c0_lamt_Xt[LO]))).imag() + ((*(c0_lamt_Xt[NLO]))).imag() + ((*(c0_lamt_Xt[NLO_QED11]))).imag() ) / (4. * mySM.getGF() / sqrt(2.) * mySM.alphaMz() / 2. / M_PI / mySM.sW2_ND())  / mySM.getCKM().computelamt().imag() << std::endl;
            return ( ((*(c0_lamt_Xt[LO]))(0) + (*(c0_lamt_Xt[NLO]))(0) + (*(c0_lamt_Xt[NLO_QED11]))(0))*NP_contribution + (*(c0_lamc_Xc[LO]))(0) + (*(c0_lamc_Xc[NLO]))(0) + (*(c0_lamc_Xc[NNLO]))(0) + (*(c0_lamc_Xc[LO_QED]))(0) + (*(c0_lamc_Xc[NLO_QED11]))(0) + myBR_Kppnunu.LongDistance()  ).abs2();
            break;
        case LO_QED:
            return( ((*(c0_lamt_Xt[LO]))(0) + (*(c0_lamt_Xt[NLO]))(0))*NP_contribution  + (*(c0_lamc_Xc[LO]))(0) + (*(c0_lamc_Xc[NLO]))(0) + (*(c0_lamc_Xc[NNLO]))(0) + (*(c0_lamc_Xc[LO_QED]))(0) + myBR_Kppnunu.LongDistance()  ).abs2()  ;
            break;
        case NO_QED:
            switch(order) {
                case NNLO:
                    ( ((*(c0_lamt_Xt[LO]))(0) + (*(c0_lamt_Xt[NLO]))(0))*NP_contribution  + (*(c0_lamc_Xc[LO]))(0) + (*(c0_lamc_Xc[NLO]))(0) + (*(c0_lamc_Xc[NNLO]))(0) + myBR_Kppnunu.LongDistance()  ).abs2()  ;
                    break;
                case NLO:
                    return ( ((*(c0_lamt_Xt[LO]))(0) + (*(c0_lamt_Xt[NLO]))(0))*NP_contribution  + (*(c0_lamc_Xc[LO]))(0) + (*(c0_lamc_Xc[NLO]))(0)  + (*(c0_lamc_Xc[LO_QED]))(0) + myBR_Kppnunu.LongDistance()  ).abs2()    ;
                    break;
                case LO:
                    return ( ((*(c0_lamt_Xt[LO]))(0))*NP_contribution + (*(c0_lamc_Xc[LO]))(0) + myBR_Kppnunu.LongDistance()  ).abs2()   ;
                    break;
                default:
                    std::stringstream out;
                    out << order;
                    throw std::runtime_error("NPBR_kppnunu::BRKppnunu_NP(): order " + out.str() + " not implemented");
            }
        default:
            std::stringstream out;
            out << order_qed;
            throw std::runtime_error("NPBR_kppnunu::BRKppnunu_NP(): order_qed " + out.str() + " not implemented");
    }    
}