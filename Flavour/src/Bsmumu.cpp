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
    /* Temporary usage of MVll class here below  */
    // gslpp::vector<gslpp::complex> ** allcoeffmumu = SM.getFlavour().ComputeCoeffsmumu(mu, NDR);
    gslpp::vector<gslpp::complex> ** allcoeff = SM.getFlavour().ComputeCoeffBMll(mu, lep); //check the mass scale, scheme fixed to NDR
    gslpp::vector<gslpp::complex> ** allcoeffprime = SM.getFlavour().ComputeCoeffprimeBMll(mu, lep);

    double alsmu = evolbsmm.alphatilde_s(mu);
    double alemu = evolbsmm.alphatilde_e(mu);
//    double alemu = SM.ale_OS(mu)/4./M_PI; // to be checked
    double sw = sqrt( (M_PI * SM.getAle() ) / ( sqrt(2.) * SM.getGF() * SM.Mw() * SM.Mw()) );


    C_10 = (*(allcoeff[LO]))(9) + (*(allcoeff[NLO]))(9);
    C_10p = (*(allcoeffprime[LO]))(9) + (*(allcoeffprime[NLO]))(9);
    C_S = (*(allcoeff[LO]))(10) + (*(allcoeff[NLO]))(10);
    C_Sp = (*(allcoeffprime[LO]))(10) + (*(allcoeffprime[NLO]))(10);
    C_P = (*(allcoeff[LO]))(11) + (*(allcoeff[NLO]))(11);
    C_Pp = (*(allcoeffprime[LO]))(11) + (*(allcoeffprime[NLO]))(11);
    
    if((order == FULLNLO) && (order_qed == FULLNLO_QED)){
    
    switch(order_qed) {
        case FULLNLO_QED:
        {
            /* Implementation to be corrected and updated with new EVO: At present better to use MVll!
            gslpp::complex C10_SM = (*(allcoeffmumu[LO]))(7) /alemu  + (*(allcoeffmumu[NLO]))(7) * alsmu/alemu 
                    + (*(allcoeffmumu[NNLO]))(7) * alsmu * alsmu/alemu + (*(allcoeffmumu[LO_QED ]))(7) /alsmu
                    + (*(allcoeffmumu[NLO_QED11]))(7) + (*(allcoeffmumu[NLO_QED02]))(7) * alemu /alsmu /alsmu 
                    + (*(allcoeffmumu[NLO_QED21]))(7) * alsmu 
                    + (*(allcoeffmumu[NLO_QED12]))(7) * alemu /alsmu+ (*(allcoeffmumu[NLO_QED22]))(7) * alemu; */
            /* Temporary usage of MVll result */
            //std::cout << " C10_SM " << C10_SM / sw / sw / SM.computelamt_s() << std::endl;
//            gslpp::complex C10_SM_plus_NP = SM.computelamt_s() * sw * sw * ((*(allcoeff[LO]))(9) + (*(allcoeff[NLO]))(9));
            //std::cout << " C10_SM_plus_NP " << C10_SM_plus_NP / sw / sw / SM.computelamt_s() << std::endl;
            gslpp::complex CC_P = SM.computelamt_s() * sw * sw * ( C_10 - C_10p + mBs*mBs / ( 2.*mlep*(mb+ms) ) * (C_P - C_Pp) );
          
            absP = CC_P.abs(); //contains only SM contributions (P, P', S, S' not added)
            argP = CC_P.arg();
            
            gslpp::complex CC_S = SM.computelamt_s() * sw * sw * ( beta * mBs*mBs / ( 2.*mlep*(mb+ms) ) * (C_S - C_Sp) );

            absS = CC_S.abs();
            argS = CC_S.arg();
           
            phiNP = 0.;
            
            ampSq = absP * absP + absS * absS;
                  
        }
        break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Bsmumu::computeAmpSq(): order " + out.str() + " not implemented");;
        }
    }
        
}
