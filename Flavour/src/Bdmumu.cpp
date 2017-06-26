/*
 * Copyright (C) 2012 SusyFit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Bdmumu.h"
#include "StandardModel.h"
#include "EvolDF1.h"
#include "HeffDB1.h"

Bdmumu::Bdmumu(const StandardModel& SM_i, int obsFlag)
: ThObservable(SM_i),
  evolbdmm(*(new EvolDF1(8, "CPL", NDR, NNLO, SM_i)))
{  
    if (obsFlag > 0 and obsFlag < 5) obs = obsFlag;
    else throw std::runtime_error("obsFlag in Bdmumu(myFlavour, obsFlag) called from ThFactory::ThFactory() can only be 1 (BR) or 2 (BRbar) or 3 (Amumu) or 4 (Smumu)");
    SM.initializeMeson(QCD::B_D);
};

double Bdmumu::computeThValue()
{   
    computeObs(FULLNLO, FULLNLO_QED);
    double FBd = SM.getMesons(QCD::B_D).getDecayconst();
    
    double coupling = SM.getGF() * SM.getGF() * SM.Mw() * SM.Mw() /M_PI /M_PI ; 
    
    double PRF = pow(coupling, 2.) / M_PI /8. / SM.getMesons(QCD::B_D).computeWidth() * pow(FBd, 2.) * pow(mmu, 2.) * mBd * beta;
    yd = 0; // For now. To be explicitly calculated.
    timeInt = (1. + Amumu * yd) / (1. - yd * yd); // Note modification in form due to algorithm
    
    if (obs == 1) return( PRF * ampSq);
    if (obs == 2) return( PRF * ampSq * timeInt);
    if (obs == 3) return( Amumu );
    if (obs == 4) return( Smumu );

    throw std::runtime_error("Bdmumu::computeThValue(): Observable type not defined. Can be only any of (1,2,3,4)");
    return (EXIT_FAILURE);
}

void Bdmumu::computeObs(orders order, orders_qed order_qed)
{
    double mu = SM.getMub();  
        
    mmu = SM.getLeptons(StandardModel::MU).getMass();
    mBd = SM.getMesons(QCD::B_D).getMass();
    mb = SM.getQuarks(QCD::BOTTOM).getMass();
    md = SM.getQuarks(QCD::DOWN).getMass();
    chiral = pow(mBd, 2.) / 2. / mmu * mb / (mb + md);
    beta = sqrt(1. - pow(2. * mmu / mBd, 2.));
    computeAmpSq(order, order_qed, mu);
    Amumu = (absP * absP * cos(2. * argP - phiNP) -  absS * absS * cos(2. * argS - phiNP)) / (absP * absP + absS * absS);
    Smumu = (absP * absP * sin(2. * argP - phiNP) -  absS * absS * sin(2. * argS - phiNP)) / (absP * absP + absS * absS);
}

double Bdmumu::computeAmumu(orders order)
{
    computeObs(FULLNLO, FULLNLO_QED);
    return(Amumu);
}

double Bdmumu::computeSmumu(orders order)
{
    computeObs(FULLNLO, FULLNLO_QED);
    return(Smumu);
}

void Bdmumu::computeAmpSq(orders order, orders_qed order_qed, double mu)
{
    if (SM.getFlavour().getHDB1().getCoeffdmumu().getOrder() < order % 3){
        std::stringstream out;
        out << order;
        throw std::runtime_error("Bdmumu::computeAmpSq(): required cofficient of "
                                 "order " + out.str() + " not computed");
    }
    gslpp::vector<gslpp::complex> ** allcoeff = SM.getFlavour().ComputeCoeffdmumu(mu, NDR);
    
    double alsmu = evolbdmm.alphatilde_s(mu);
//    double alemu = evolbdmm.alphatilde_e(mu);
    double alemu = SM.ale_OS(mu)/4./M_PI; // to be checked
    
    if((order == FULLNLO) && (order_qed == FULLNLO_QED)){
    
    switch(order_qed) {
        case FULLNLO_QED:
        {
            gslpp::complex CC = (*(allcoeff[LO]))(7) /alemu  + (*(allcoeff[NLO]))(7) * alsmu/alemu 
                    + (*(allcoeff[NNLO]))(7) * alsmu * alsmu/alemu + (*(allcoeff[LO_QED ]))(7) /alsmu
                    + (*(allcoeff[NLO_QED11]))(7) + (*(allcoeff[NLO_QED02]))(7) * alemu /alsmu /alsmu 
                    + (*(allcoeff[NLO_QED21]))(7) * alsmu 
                    + (*(allcoeff[NLO_QED12]))(7) * alemu /alsmu+ (*(allcoeff[NLO_QED22]))(7) * alemu;
            absP = CC.abs();
            argP = CC.arg();
           
            phiNP = 0.;
            
            ampSq = absP * absP ;
                  
        }
        break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Bdmumu::computeAmpSq(): order " + out.str() + " not implemented");;
    }
    }
    
    
}