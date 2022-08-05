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
    gslpp::matrix<gslpp::complex> Vckm = SM.getVCKM();
    double sw = sqrt( (M_PI * SM.getAle() ) / ( sqrt(2.) * SM.getGF() * SM.Mw() * SM.Mw()) ) ;
    
    gslpp::complex C_10p = 0.;
    gslpp::complex C_S = 0.;
    gslpp::complex C_Sp = 0.;
    gslpp::complex C_P = 0.;
    gslpp::complex C_Pp = 0.;
   
    if((order == FULLNLO) && (order_qed == FULLNLO_QED)){
    
    switch(order_qed) {
        case FULLNLO_QED:
        {
            gslpp::complex C10_SM = (*(allcoeff[LO]))(7) /alemu  + (*(allcoeff[NLO]))(7) * alsmu/alemu 
                    + (*(allcoeff[NNLO]))(7) * alsmu * alsmu/alemu + (*(allcoeff[LO_QED ]))(7) /alsmu
                    + (*(allcoeff[NLO_QED11]))(7) + (*(allcoeff[NLO_QED02]))(7) * alemu /alsmu /alsmu 
                    + (*(allcoeff[NLO_QED21]))(7) * alsmu 
                    + (*(allcoeff[NLO_QED12]))(7) * alemu /alsmu+ (*(allcoeff[NLO_QED22]))(7) * alemu;
            
            gslpp::complex NPfactor = (Vckm(2,2).conjugate() * Vckm(2,1)) * sw * sw;

            gslpp::complex CC_P = C10_SM + NPfactor * ( /*C10_NP*/ - C_10p + mBs*mBs*mb / ( 2.*mlep*(mb+ms)*mW ) * (C_P - C_Pp) );

                
            absP = CC_P.abs(); //contains only SM contributions (P, P', S, S' not added)
            argP = CC_P.arg();
            
            gslpp::complex CC_S = NPfactor * ( beta * mBs*mBs*mb / ( 2.*mlep*(mb+ms)*mW ) * (C_S - C_Sp) );

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