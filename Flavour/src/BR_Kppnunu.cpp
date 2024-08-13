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
: ThObservable(SM_i), mySM(SM_i)
{   
    mySM.initializeMeson(QCD::K_P);
    
    setParametersForObservable(make_vector<std::string>() <<  "DeltaP_cu"  << "IB_Kp" << "Vus_fpK0Pip" << "PhSp_KP" << "Delta_EM"); 
   
}

double BR_Kppnunu::computeThValue()
{
    //std::cout <<" k+" << k_plus() * pow( 4. * mySM.getGF() / sqrt(2.) * mySM.alphaMz() / 2. / M_PI / mySM.sW2_ND(),2.) * pow(mySM.getCKM().getLambda() ,10. ) << std::endl;
    return(  k_plus() * (1. + mySM.getOptionalParameter("Delta_EM")) * BRKppnunu(NLO, NLO_QED11) ); 
}

double BR_Kppnunu::BRKppnunu(orders order, orders_qed order_qed)
{   
    if (mySM.getFlavour().getHDS1().getCoeffDS1pnunu().getOrder() < order){
        std::stringstream out;
        out << order;
        throw std::runtime_error("BRKppnunu::computeThValue(): requires cofficient of "
                                 "order" + out.str() + "not computed");
    }
    
    if (mySM.getFlavour().getHDS1().getCoeffDS1pnunu().getOrder_qed()< order_qed){ 
        std::stringstream out2;
        out2 << order_qed;
        throw std::runtime_error("BRKppnunu::computeThValue(): requires cofficient of "
                                 " qed order " + out2.str() + " not computed");
    }
    
    
    gslpp::vector<gslpp::complex> ** c0_lamt_Xt = mySM.getFlavour().ComputeCoeffDS1pnunu() ; 
    gslpp::vector<gslpp::complex> ** c0_lamc_Xc = mySM.getFlavour().ComputeCoeffDS1pnunuC(); 
    
    switch(order_qed) {
        case NLO_QED11:
            //std::cout << "Pc: " << ( (*(c0_lamc_Xc[LO]))(0) + (*(c0_lamc_Xc[NLO]))(0) + (*(c0_lamc_Xc[NNLO]))(0) + (*(c0_lamc_Xc[LO_QED]))(0) + (*(c0_lamc_Xc[NLO_QED11]))(0)).real() / (4. * mySM.getGF() / sqrt(2.) * mySM.alphaMz() / 2. / M_PI / mySM.sW2_ND()) / pow(mySM.getCKM().getLambda(),4.) / mySM.getCKM().computelamc().real()  << std::endl;
            //std::cout << "deltaPC: " << LongDistance().real() / (4. * mySM.getGF() / sqrt(2.) * mySM.alphaMz() / 2. / M_PI / mySM.sW2_ND()) / pow(mySM.getCKM().getLambda(),4.) / mySM.getCKM().computelamc().real() << std::endl;
            //std::cout << "Xt: " << ( ((*(c0_lamt_Xt[LO]))).imag() + ((*(c0_lamt_Xt[NLO]))).imag() + ((*(c0_lamt_Xt[NLO_QED11]))).imag() ) / (4. * mySM.getGF() / sqrt(2.) * mySM.alphaMz() / 2. / M_PI / mySM.sW2_ND())  / mySM.getCKM().computelamt().imag() << std::endl;
            return (( (*(c0_lamt_Xt[LO]))(0) + (*(c0_lamt_Xt[NLO]))(0) + (*(c0_lamt_Xt[NLO_QED11]))(0) + (*(c0_lamc_Xc[LO]))(0) + (*(c0_lamc_Xc[NLO]))(0) + (*(c0_lamc_Xc[NNLO]))(0) + (*(c0_lamc_Xc[LO_QED]))(0) + (*(c0_lamc_Xc[NLO_QED11]))(0) + LongDistance()  ).abs2()
            + ( (*(c0_lamt_Xt[LO]))(1) + (*(c0_lamt_Xt[NLO]))(1) + (*(c0_lamt_Xt[NLO_QED11]))(1)).abs2() );
            break;
        case LO_QED:
            return(( (*(c0_lamt_Xt[LO]))(0) + (*(c0_lamt_Xt[NLO]))(0)  + (*(c0_lamc_Xc[LO]))(0) + (*(c0_lamc_Xc[NLO]))(0) + (*(c0_lamc_Xc[NNLO]))(0) + (*(c0_lamc_Xc[LO_QED]))(0) + LongDistance()  ).abs2() 
            + ( (*(c0_lamt_Xt[LO]))(1) + (*(c0_lamt_Xt[NLO]))(1)  ).abs2())  ;
            break;
        case NO_QED:
            switch(order) {
                case NNLO:
                    return (( (*(c0_lamt_Xt[LO]))(0) + (*(c0_lamt_Xt[NLO]))(0)  + (*(c0_lamc_Xc[LO]))(0) + (*(c0_lamc_Xc[NLO]))(0) + (*(c0_lamc_Xc[NNLO]))(0) + LongDistance()  ).abs2() 
                    + ( (*(c0_lamt_Xt[LO]))(1) + (*(c0_lamt_Xt[NLO]))(1)  ).abs2());
                    break;
                case NLO:
                    return (( (*(c0_lamt_Xt[LO]))(0) + (*(c0_lamt_Xt[NLO]))(0)  + (*(c0_lamc_Xc[LO]))(0) + (*(c0_lamc_Xc[NLO]))(0)  + (*(c0_lamc_Xc[LO_QED]))(0) + LongDistance()  ).abs2()   
                    + ( (*(c0_lamt_Xt[LO]))(1) + (*(c0_lamt_Xt[NLO]))(1) ).abs2() );
                    break;
                case LO:
                    return (( (*(c0_lamt_Xt[LO]))(0) + (*(c0_lamc_Xc[LO]))(0) + LongDistance()  ).abs2()   
                    + ( (*(c0_lamt_Xt[LO]))(1) ).abs2() );
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
    
    return (mySM.getOptionalParameter("IB_Kp") * mySM.getOptionalParameter("IB_Kp") * mySM.getOptionalParameter("Vus_fpK0Pip") * mySM.getOptionalParameter("Vus_fpK0Pip") * mySM.getOptionalParameter("PhSp_KP") *
            pow(mySM.getMesons(QCD::K_P).getMass(),5.) * mySM.getMesons(QCD::K_P).getLifetime() / HCUT / 512 / M_PI / M_PI / M_PI / mySM.getCKM().getLambda() / mySM.getCKM().getLambda());
    
}

gslpp::complex BR_Kppnunu::LongDistance() // c0 * lambda_c * delta_Pc * lambda^4
{   
    return 4. * mySM.getGF() / sqrt(2.) * mySM.alphaMz() / 2. / M_PI / mySM.sW2_ND() * mySM.getCKM().computelamc() * mySM.getOptionalParameter("DeltaP_cu") * pow(mySM.getCKM().getLambda(),4.)  ;
}
