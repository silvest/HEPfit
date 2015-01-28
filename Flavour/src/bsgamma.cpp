/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "bsgamma.h"

Bsgamma::Bsgamma(const StandardModel& SM_i, int obsFlag): ThObservable(SM_i), mySM(SM_i){
    if (obsFlag > 0 and obsFlag < 3) obs = obsFlag;
    else throw std::runtime_error("obsFlag in bsgamma can only be 1 (BR) or 2 (A_{CP})");
};

void Bsgamma::computeCoeff(orders order){
    allcoeff = mySM.getMyFlavour()->ComputeCoeffsgamma(mu_b);
    
    switch(order) {
        case FULLNNLO: //coeff = ;
        case FULLNLO: //coeff = ;
        case LO: //coeff = ;
        break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Bsmumu::computeAmpSq(): order " + out.str() + " not implemented");;
    }
}

void Bsgamma::computeGamma(orders order){
    GF = mySM.getGF();
    ale=mySM.getAle();
    Mb=mySM.getQuarks(QCD::BOTTOM).getMass();
    mu_b = mySM.getMub();
    width = mySM.getMesons(StandardModel::B_D).computeWidth();
    lambda_t=mySM.computelamt_s();
    
    computeCoeff(order);
    coeff=1.;
    
    Gamma = GF*GF*pow(Mb,5.)*ale*lambda_t.abs2()/(32. * pow(M_PI,4.)) * coeff;
    Gamma_conj = GF*GF*pow(Mb,5.)*ale*(lambda_t.conjugate()).abs2()/(32. * pow(M_PI,4.)) * coeff;
}

double Bsgamma::computeThValue(){
    computeGamma(FULLNNLO);
    
    if (obs == 1) return Gamma/width;
    if (obs == 2) return (Gamma - Gamma_conj) / (Gamma + Gamma_conj);
    
    throw std::runtime_error("Bsgamma::computeThValue(): Observable type not defined. Can be only any of (1,2)");
}