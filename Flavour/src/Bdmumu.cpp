/*
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Bdmumu.h"

Bdmumu::Bdmumu(const StandardModel& SM_i, int obsFlag): ThObservable(SM_i), mySM(SM_i){
    if (obsFlag > 0 and obsFlag < 5) obs = obsFlag;
    else throw std::runtime_error("obsFlag in Bsmumu(myFlavour, obsFlag) called from ThFactory::ThFactory() can only be 1 (BR) or 2 (BRbar) or 3 (Amumu) or 4 (Smumu)");
};

double Bdmumu::computeThValue(){
    computeObs(FULLNLO);
    double FBd = mySM.getMesons(QCD::B_D).getDecayconst();
    double coupling = mySM.getGF() * mySM.alphaMz() / 4. / M_PI;
    double PRF = pow(coupling, 2.) / M_PI / mySM.getMesons(QCD::B_D).computeWidth() * pow(FBd, 2.) * pow(mmu, 2.) * mBd * beta;
    yd = 0; // For now. To be explicitly calculated.
    timeInt = (1. + Amumu * yd) / (1. - yd * yd); // Note modification in form due to algorithm
    /*    double theta = asin(sqrt( (M_PI * mySM.getAle() )/( sqrt(2) * mySM.getGF() *
     mySM.Mw_tree() * mySM.Mw_tree()) ));
     
     return( mySM.getMesons(QCD::B_S).getLifetime() / HCUT * mySM.getGF()*mySM.getGF()/M_PI
     * mySM.getAle()*mySM.getAle()/(16.*M_PI*M_PI*pow(sin(theta),4.))
     * mBs * mySM.getMesons(QCD::B_S).getDecayconst()
     * mmu * mmu * sqrt(1.-4.*mmu*mmu/mBs/mBs) * BRBsmumu(NLO).real());*/
    
    //std::cout << getAmumu(FULLNLO) << "  " << getSmumu(FULLNLO) << argP << std::endl;
    
    if (obs == 1) return( PRF * ampSq);
    if (obs == 2) return( PRF * ampSq * timeInt);
    if (obs == 3) return( Amumu );
    if (obs == 4) return( Smumu );

    throw std::runtime_error("Bdmumu::computeThValue(): Observable type not defined. Can be only any of (1,2,3,4)");
    return (EXIT_FAILURE);
}

void Bdmumu::computeObs(orders order){
    mmu = mySM.getLeptons(StandardModel::MU).getMass();
    mBd = mySM.getMesons(QCD::B_D).getMass();
    mb = mySM.getQuarks(QCD::BOTTOM).getMass();
    md = mySM.getQuarks(QCD::DOWN).getMass();
    chiral = pow(mBd, 2.) / 2. / mmu * mb / (mb + md);
    beta = sqrt(1. - pow(2. * mmu / mBd, 2.));
    computeAmpSq(order);
    Amumu = (absP * absP * cos(2. * argP - phiNP) -  absS * absS * cos(2. * argS - phiNP)) / (absP * absP + absS * absS);
    Smumu = (absP * absP * sin(2. * argP - phiNP) -  absS * absS * sin(2. * argS - phiNP)) / (absP * absP + absS * absS);
}

double Bdmumu::computeAmumu(orders order){
    computeObs(FULLNLO);
    return(Amumu);
}

double Bdmumu::computeSmumu(orders order){
    computeObs(FULLNLO);
    return(Smumu);
}

void Bdmumu::computeAmpSq(orders order){
    if (mySM.getMyFlavour()->getHDB1().getCoeffdmumu().getOrder() < order % 3){
        std::stringstream out;
        out << order;
        throw std::runtime_error("Bdmumu::computeAmpSq(): required cofficient of "
                                 "order " + out.str() + " not computed");
    }
    vector<complex> ** allcoeff = mySM.getMyFlavour()->ComputeCoeffdmumu();
    
    switch(order) {
        case FULLNLO:
        {
            complex PP = (*(allcoeff[LO]) + *(allcoeff[NLO]))(0) - (*(allcoeff[LO]) + *(allcoeff[NLO]))(1)
            + chiral * ((*(allcoeff[LO]) + *(allcoeff[NLO]))(2) - (*(allcoeff[LO]) + *(allcoeff[NLO]))(3));
            absP = PP.abs();
            argP = PP.arg();
            
            complex SS = beta * chiral * ((*(allcoeff[LO]) + *(allcoeff[NLO]))(4) - (*(allcoeff[LO]) + *(allcoeff[NLO]))(5));
            absS = SS.abs();
            argS = SS.arg();
            phiNP = 0.;
            
            ampSq = absP * absP + absS * absS;
        }
            break;
        case LO:
        {
            complex PP = (*(allcoeff[LO]))(0) - (*(allcoeff[LO]))(1)
            + chiral * ((*(allcoeff[LO]))(2) - (*(allcoeff[LO]))(3));
            absP = PP.abs();
            argP = PP.arg();
            
            complex SS = beta * chiral * ((*(allcoeff[LO]))(4) - (*(allcoeff[LO]))(5));
            absS = SS.abs();
            argS = SS.arg();
            phiNP = 0.;
            
            ampSq = absP * absP + absS * absS;
        }
            break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Bdmumu::computeAmpSq(): order " + out.str() + " not implemented");;
    }
}