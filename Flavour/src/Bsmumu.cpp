/*
 * Copyright (C) 2012 SusyFit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Bsmumu.h"
#include "StandardModel.h"
#include "EvolBsmm.h"
#include "HeffDB1.h"

Bsmumu::Bsmumu(const StandardModel& SM_i, int obsFlag, QCD::lepton lep_i)
: ThObservable(SM_i),
  evolbsmm(*(new EvolBsmm(8, NDR, NNLO, NLO_QED22, SM)))
{
    lep = lep_i;
    if (obsFlag > 0 and obsFlag < 5) obs = obsFlag;
    else throw std::runtime_error("obsFlag in Bsmumu(myFlavour, obsFlag) called from ThFactory::ThFactory() can only be 1 (BR) or 2 (BRbar) or 3 (Amumu) or 4 (Smumu)");
    SM.initializeMeson(QCD::B_S);
};

double Bsmumu::computeThValue()
{
    computeObs(FULLNLO, FULLNLO_QED);
    double FBs = SM.getMesons(QCD::B_S).getDecayconst();
    
    double coupling = SM.getGF() * SM.getGF() * SM.Mw() * SM.Mw() /M_PI /M_PI ; 
 
    double PRF = pow(coupling, 2.) / M_PI /8. / SM.getMesons(QCD::B_S).computeWidth() * pow(FBs, 2.) * pow(mlep, 2.) * mBs * beta;
    double ys = SM.getMesons(QCD::B_S).getDgamma_gamma()/2.; // For now. To be explicitly calculated.
    timeInt = (1. + Amumu * ys) / (1. - ys * ys); // Note modification in form due to algorithm
     
    if (obs == 1) return( PRF * ampSq);
    if (obs == 2) return( PRF * ampSq * timeInt);
    if (obs == 3) return( Amumu );
    if (obs == 4) return( Smumu );
    
    throw std::runtime_error("Bsmumu::computeThValue(): Observable type not defined. Can be only any of (1,2,3,4)");
    return (EXIT_FAILURE);
}

void Bsmumu::computeObs(orders order, orders_qed order_qed)
{   
    double mu = SM.getMub();  
    
    mlep = SM.getLeptons(lep).getMass();
    mBs = SM.getMesons(QCD::B_S).getMass();
    mW = SM.Mw();
    mb = SM.getQuarks(QCD::BOTTOM).getMass();
    ms = SM.getQuarks(QCD::STRANGE).getMass();
    chiral = pow(mBs, 2.) / 2. / mlep * mb / (mb + ms);
    beta = sqrt(1. - pow(2. * mlep / mBs, 2.));
    computeAmpSq(order, order_qed, mu);
    Amumu = (absP * absP * cos(2. * argP - phiNP) -  absS * absS * cos(2. * argS - phiNP)) / (absP * absP + absS * absS);
    Smumu = (absP * absP * sin(2. * argP - phiNP) -  absS * absS * sin(2. * argS - phiNP)) / (absP * absP + absS * absS);
}

double Bsmumu::computeAmumu(orders order)
{
    computeObs(FULLNLO, FULLNLO_QED);
    return(Amumu);
}

double Bsmumu::computeSmumu(orders order)
{
    computeObs(FULLNLO, FULLNLO_QED);
    return(Smumu);
}

void Bsmumu::computeAmpSq(orders order, orders_qed order_qed, double mu)
{  
    if (SM.getFlavour().getHDB1().getCoeffsmumu().getOrder() < order % 3){
        std::stringstream out;
        out << order;
        throw std::runtime_error("Bsmumu::computeAmpSq(): required cofficient of "
                                 "order " + out.str() + " not computed");
    }
    gslpp::vector<gslpp::complex> ** allcoeff = SM.getFlavour().ComputeCoeffsmumu(mu, NDR);

    double alsmu = evolbsmm.alphatilde_s(mu);
    double alemu = evolbsmm.alphatilde_e(mu);
//    double alemu = SM.ale_OS(mu)/4./M_PI; // to be checked
    double sw = sqrt( (M_PI * SM.getAle() ) / ( sqrt(2.) * SM.getGF() * SM.Mw() * SM.Mw()) );
    
    bool WET_NP_btos = SM.getFlavour().getFlagWET_NP_btos();
    bool SMEFT_NP_btos = SM.getFlavour().getFlagSMEFT_NP_btos();
    
    if(WET_NP_btos){
        if(lep == StandardModel::ELECTRON){
            C_10_NP = SM.getOptionalParameter("C10_e");
            C_10p_NP = SM.getOptionalParameter("C10p_e");
            C_S_NP = SM.getOptionalParameter("CS_e");
            C_Sp_NP = SM.getOptionalParameter("CSp_e");
            C_P_NP = SM.getOptionalParameter("CP_e");
            C_Pp_NP = SM.getOptionalParameter("CPp_e");
        }
        else{
            C_10_NP = SM.getOptionalParameter("C10_mu");
            C_10p_NP = SM.getOptionalParameter("C10p_mu");
            C_S_NP = SM.getOptionalParameter("CS_mu");
            C_Sp_NP = SM.getOptionalParameter("CSp_mu");
            C_P_NP = SM.getOptionalParameter("CP_mu");
            C_Pp_NP = SM.getOptionalParameter("CPp_mu");     
        }
    }
    else if(SMEFT_NP_btos){
        // normalization to \Lambda = 1 TeV
        gslpp::complex SMEFT_factor = (M_PI/SM.getAle())*(SM.v()*1.e-3)*(SM.v()*1.e-3)/SM.computelamt_s();
        if(lep == StandardModel::ELECTRON){
            C_10_NP = SM.getOptionalParameter("CQe23_11");
            C_10_NP -= SM.getOptionalParameter("C1LQ11_23");
            C_10_NP -= SM.getOptionalParameter("C3LQ11_23");
            C_10_NP += SM.getOptionalParameter("C1HQ");
            C_10_NP += SM.getOptionalParameter("C3HQ");
            C_10_NP *= SMEFT_factor; 
            C_10p_NP = SM.getOptionalParameter("Ced11_23");
            C_10p_NP -= SM.getOptionalParameter("CLd11_23");
            C_10p_NP += SM.getOptionalParameter("CHd");
            C_10p_NP *= SMEFT_factor; 
            C_S_NP = SM.getOptionalParameter("CLedQ_11");
            C_S_NP *= SMEFT_factor; 
            C_Sp_NP = SM.getOptionalParameter("CpLedQ_11");
            C_Sp_NP *= SMEFT_factor;
            C_P_NP = -SM.getOptionalParameter("CLedQ_11");
            C_P_NP *= SMEFT_factor;
            C_Pp_NP = SM.getOptionalParameter("CpLedQ_11");
            C_Pp_NP *= SMEFT_factor;
        }
        else{
            C_10_NP = SM.getOptionalParameter("CQe23_22");
            C_10_NP -= SM.getOptionalParameter("C1LQ22_23");
            C_10_NP -= SM.getOptionalParameter("C3LQ22_23");
            C_10_NP += SM.getOptionalParameter("C1HQ");
            C_10_NP += SM.getOptionalParameter("C3HQ");
            C_10_NP *= SMEFT_factor; 
            C_10p_NP = SM.getOptionalParameter("Ced22_23");
            C_10p_NP -= SM.getOptionalParameter("CLd22_23");
            C_10p_NP += SM.getOptionalParameter("CHd");
            C_10p_NP *= SMEFT_factor; 
            C_S_NP = SM.getOptionalParameter("CLedQ_22");
            C_S_NP *= SMEFT_factor; 
            C_Sp_NP = SM.getOptionalParameter("CpLedQ_22");
            C_Sp_NP *= SMEFT_factor;
            C_P_NP = -SM.getOptionalParameter("CLedQ_22");
            C_P_NP *= SMEFT_factor;
            C_Pp_NP = SM.getOptionalParameter("CpLedQ_22");
            C_Pp_NP *= SMEFT_factor;        
        }
    }
    else{
        C_10_NP = 0.;
        C_10p_NP = 0.;
        C_S_NP = 0.;
        C_Sp_NP = 0.;
        C_P_NP = 0.;
        C_Pp_NP = 0.;
    }
   
    if((order == FULLNLO) && (order_qed == FULLNLO_QED)){
    
    switch(order_qed) {
        case FULLNLO_QED:
        {
            gslpp::complex C10_SM = (*(allcoeff[LO]))(7) /alemu  + (*(allcoeff[NLO]))(7) * alsmu/alemu 
                    + (*(allcoeff[NNLO]))(7) * alsmu * alsmu/alemu + (*(allcoeff[LO_QED ]))(7) /alsmu
                    + (*(allcoeff[NLO_QED11]))(7) + (*(allcoeff[NLO_QED02]))(7) * alemu /alsmu /alsmu 
                    + (*(allcoeff[NLO_QED21]))(7) * alsmu 
                    + (*(allcoeff[NLO_QED12]))(7) * alemu /alsmu+ (*(allcoeff[NLO_QED22]))(7) * alemu;
            
            gslpp::complex CC_P = C10_SM + SM.computelamt_s() * sw * sw * ( C_10_NP - C_10p_NP + mBs*mBs*mb / ( 2.*mlep*(mb+ms)*mW ) * (C_P_NP - C_Pp_NP) );
          
            absP = CC_P.abs(); //contains only SM contributions (P, P', S, S' not added)
            argP = CC_P.arg();
            
            gslpp::complex CC_S = SM.computelamt_s() * sw * sw * ( beta * mBs*mBs*mb / ( 2.*mlep*(mb+ms)*mW ) * (C_S_NP - C_Sp_NP) );

            absS = CC_S.abs();
            argS = CC_S.arg();
           
            phiNP = 0.;
            
            ampSq = absP * absP ;
                  
        }
        break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Bsmumu::computeAmpSq(): order " + out.str() + " not implemented");;
        }
    }
        
}