/*
 * Copyright (C) 2012 SusyFit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Bdmumu.h"
#include "StandardModel.h"
#include "EvolBsmm.h"
#include "HeffDB1.h"
#include "NPSMEFTd6GeneralMatching.h"

Bdmumu::Bdmumu(const StandardModel& SM_i, int obsFlag, QCD::lepton lep_i)
: ThObservable(SM_i),
  evolbdmm(new EvolBsmm(8, NDR, NNLO, NLO_QED22, SM))
{  
    lep = lep_i;
    if(lep == QCD::MU) leptonindex = 1;
    else if(lep == QCD::ELECTRON) leptonindex = 0;
    else if(lep == QCD::TAU) leptonindex = 2;
    else throw std::runtime_error("Bsmumu::Bsmumu(): lepton type not defined. Can be only MU or ELECTRON");
    if (obsFlag > 0 and obsFlag < 5) obs = obsFlag;
    else throw std::runtime_error("obsFlag in Bdmumu(myFlavour, obsFlag) called from ThFactory::ThFactory() can only be 1 (BR) or 2 (BRbar) or 3 (Amumu) or 4 (Smumu)");
    SM.initializeMeson(QCD::B_D);
};

double Bdmumu::computeThValue()
{   
    double FBd = SM.getMesons(QCD::B_D).getDecayconst();
    
    double Mw = SM.Mw();
    double GF = SM.getGF();
    coupling = GF * GF * Mw * Mw /M_PI /M_PI ; 

    computeObs(FULLNLO, FULLNLO_QED);
    

    double PRF = pow(coupling, 2.) / M_PI /8. / SM.getMesons(QCD::B_D).computeWidth() * pow(FBd, 2.) * pow(mlep, 2.) * mBd * beta;
    yd = 0.; // Use the experimental number here
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
        
    mlep = SM.getLeptons(lep).getMass();
    mBd = SM.getMesons(QCD::B_D).getMass();
    mW = SM.Mw();
    mb = SM.getQuarks(QCD::BOTTOM).getMass();
    md = SM.getQuarks(QCD::DOWN).getMass();
    chiral = pow(mBd, 2.) / 2. / mlep * mb / (mb + md);
    beta = sqrt(1. - pow(2. * mlep / mBd, 2.));
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
    allcoeff = SM.getFlavour().ComputeCoeffdmumu(mu, NDR);
    
    double alsmu = evolbdmm->alphatilde_s(mu);
    double alemu = evolbdmm->alphatilde_e(mu);
//    double alemu = SM.ale_OS(mu)/4./M_PI; // to be checked
    gslpp::matrix<gslpp::complex> Vckm = SM.getVCKM();
    
    if((order == FULLNLO) && (order_qed == FULLNLO_QED)){
    
    switch(order_qed) {
        case FULLNLO_QED:
        {
            C_10 = (*(allcoeff[LO]))(7) /alemu  + (*(allcoeff[NLO]))(7) * alsmu/alemu 
                    + (*(allcoeff[NNLO]))(7) * alsmu * alsmu/alemu + (*(allcoeff[LO_QED ]))(7) /alsmu
                    + (*(allcoeff[NLO_QED11]))(7) + (*(allcoeff[NLO_QED02]))(7) * alemu /alsmu /alsmu 
                    + (*(allcoeff[NLO_QED21]))(7) * alsmu 
                    + (*(allcoeff[NLO_QED12]))(7) * alemu /alsmu+ (*(allcoeff[NLO_QED22]))(7) * alemu;

            gslpp::complex NPfactor = (Vckm(2,2).conjugate() * Vckm(2,0)) * coupling;

            if(SM.getModelName().compare("NPSMEFTd6U2") == 0 || SM.getModelName().compare("NPSMEFTd6U3") == 0)
            {
                C_10 = C_10 + (dynamic_cast<const NPSMEFTd6GeneralMatching&>(SM.getMatching()).getCdeVLR(0,2,leptonindex,leptonindex) - 
                    dynamic_cast<const NPSMEFTd6GeneralMatching&>(SM.getMatching()).getCedVLL(leptonindex,leptonindex,0,2)) / NPfactor; 
            }



            absP = C_10.abs();
            argP = C_10.arg();
           
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
