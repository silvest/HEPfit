/*
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Bsmumu.h"

Bsmumu::Bsmumu(const StandardModel& SM_i, int obsFlag)
: ThObservable(SM_i),
  evolbsmm(8, NDR, NNLO, NLO_ewt4, SM)
{
    if (obsFlag > 0 and obsFlag < 5) obs = obsFlag;
    else throw std::runtime_error("obsFlag in Bsmumu(myFlavour, obsFlag) called from ThFactory::ThFactory() can only be 1 (BR) or 2 (BRbar) or 3 (Amumu) or 4 (Smumu)");
};

double Bsmumu::computeThValue()
{
    computeObs(FULLNLO, FULLNLO_ew);
    double FBs = SM.getMesons(QCD::B_S).getDecayconst();
    
    double coupling = SM.getGF() * SM.getGF() * SM.Mw() * SM.Mw() /M_PI /M_PI ; 
 
    double PRF = pow(coupling, 2.) / M_PI /8. / SM.getMesons(QCD::B_S).computeWidth() * pow(FBs, 2.) * pow(mmu, 2.) * mBs * beta;
    double ys = SM.getMesons(QCD::B_S).getDgamma_gamma()/2.; // For now. To be explicitly calculated.
    timeInt = (1. + Amumu * ys) / (1. - ys * ys); // Note modification in form due to algorithm
     
    if (obs == 1) return( PRF * ampSq);
    if (obs == 2) return( PRF * ampSq * timeInt);
    if (obs == 3) return( Amumu );
    if (obs == 4) return( Smumu );
    
    throw std::runtime_error("Bsmumu::computeThValue(): Observable type not defined. Can be only any of (1,2,3,4)");
    return (EXIT_FAILURE);
}

void Bsmumu::computeObs(orders order, orders_ew order_ew)
{   
    double mu = SM.getMub();  
    
    mmu = SM.getLeptons(StandardModel::MU).getMass();
    mBs = SM.getMesons(QCD::B_S).getMass();
    mb = SM.getQuarks(QCD::BOTTOM).getMass();
    ms = SM.getQuarks(QCD::STRANGE).getMass();
    chiral = pow(mBs, 2.) / 2. / mmu * mb / (mb + ms);
    beta = sqrt(1. - pow(2. * mmu / mBs, 2.));
    computeAmpSq(order, order_ew, mu);
    Amumu = (absP * absP * cos(2. * argP - phiNP) -  absS * absS * cos(2. * argS - phiNP)) / (absP * absP + absS * absS);
    Smumu = (absP * absP * sin(2. * argP - phiNP) -  absS * absS * sin(2. * argS - phiNP)) / (absP * absP + absS * absS);
}

double Bsmumu::computeAmumu(orders order)
{
    computeObs(FULLNLO, FULLNLO_ew);
    return(Amumu);
}

double Bsmumu::computeSmumu(orders order)
{
    computeObs(FULLNLO, FULLNLO_ew);
    return(Smumu);
}

void Bsmumu::computeAmpSq(orders order, orders_ew order_ew, double mu)
{  
    if (SM.getMyFlavour()->getHDB1().getCoeffsmumu().getOrder() < order % 3){
        std::stringstream out;
        out << order;
        throw std::runtime_error("Bsmumu::computeAmpSq(): required cofficient of "
                                 "order " + out.str() + " not computed");
    }
    gslpp::vector<gslpp::complex> ** allcoeff = SM.getMyFlavour()->ComputeCoeffsmumu(mu, NDR);

    double alsmu = evolbsmm.alphatilde_s(mu);
    double alemu = evolbsmm.alphatilde_e(mu);
   
    if((order == FULLNLO) && (order_ew == FULLNLO_ew)){
    
    switch(order_ew) {
        case FULLNLO_ew:
        {
            gslpp::complex CC = (*(allcoeff[LO]))(7) /alemu  + (*(allcoeff[NLO]))(7) * alsmu/alemu 
                    + (*(allcoeff[NNLO]))(7) * alsmu * alsmu/alemu + (*(allcoeff[LO_ew ]))(7) /alsmu
                    + (*(allcoeff[NLO_ew]))(7) + (*(allcoeff[NLO_ewt1]))(7) * alemu /alsmu /alsmu 
                    + (*(allcoeff[NLO_ewt2]))(7) * alsmu 
                    + (*(allcoeff[NLO_ewt3]))(7) * alemu /alsmu+ (*(allcoeff[NLO_ewt4]))(7) * alemu;
            absP = CC.abs(); //contains only SM contributions (P, P', S, S' not added)
            argP = CC.arg();
            
            absS = 0.;
            argS = 0.;
           
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