/* 
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "MVllObservables.h"
#include "Flavour.h"

/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/



P_1::P_1(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) : ThObservable(SM_i), mySM(SM_i) {
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double P_1::computeThValue() {
    

    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateSigma(4,q_min,q_max)/(2.* mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateSigma(3,q_min,q_max));
}


P_2::P_2(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) : ThObservable(SM_i), mySM(SM_i) {  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double P_2::computeThValue() {
    
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateSigma(7,q_min,q_max)/(8.*mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateSigma(3,q_min,q_max));
}


P_3::P_3(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) : ThObservable(SM_i), mySM(SM_i) {  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double P_3::computeThValue() {
    
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return -mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateSigma(11,q_min,q_max)/(4.*mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateSigma(3,q_min,q_max));
}


P_4Prime::P_4Prime(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) : ThObservable(SM_i), mySM(SM_i) {  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double P_4Prime::computeThValue() {
    
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateSigma(5,q_min,q_max)/sqrt(-mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateSigma(2,q_min,q_max)*mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateSigma(3,q_min,q_max));
   
}


P_5Prime::P_5Prime(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) : ThObservable(SM_i), mySM(SM_i) {  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double P_5Prime::computeThValue() {
    
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateSigma(6,q_min,q_max)/(2.*sqrt(-mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateSigma(2,q_min,q_max)*mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateSigma(3,q_min,q_max)));
}


P_6Prime::P_6Prime(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) : ThObservable(SM_i), mySM(SM_i) {  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double P_6Prime::computeThValue() {
    
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return -mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateSigma(9,q_min,q_max)/(2.*sqrt(-mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateSigma(2,q_min,q_max)*mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateSigma(3,q_min,q_max)));
 
}


GammaPrime::GammaPrime(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) : ThObservable(SM_i), mySM(SM_i) {  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double GammaPrime::computeGammaPrime(double qmin, double qmax){
    
    double q_min = qmin;
    double q_max = qmax;
    
    return ((3.*mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateSigma(0,q_min,q_max) - mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateSigma(2,q_min,q_max)) + 2.*(3.*mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateSigma(1,q_min,q_max) - mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateSigma(3,q_min,q_max)))/4.;
}

double GammaPrime::computeThValue(){

    double q_min = getBinMin();
    double q_max = getBinMax();

    return computeGammaPrime(q_min, q_max);
}


A_FB::A_FB(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) : GammaPrime(SM_i, meson_i, vector_i, lep_i), mySM(SM_i) {
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}


double A_FB::computeThValue() {

    double q_min = getBinMin();
    double q_max = getBinMax();

    return -3. * mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateSigma(7,q_min,q_max) / 4. / computeGammaPrime(q_min, q_max);
}


BR_BKstarll::BR_BKstarll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) : GammaPrime(SM_i, meson_i, vector_i, lep_i), mySM(SM_i) {  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double BR_BKstarll::computeThValue() {
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return computeGammaPrime(q_min, q_max)/mySM.getMyFlavour()->getMVll(meson, vectorM, lep).width_Bd;
}


ACP::ACP(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) : GammaPrime(SM_i, meson_i, vector_i, lep_i), mySM(SM_i) {
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double ACP::computeThValue() {
    
    double q_min = getBinMin();
    double q_max = getBinMax();
          
    return (3.*mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateDelta(0,q_min,q_max) - mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateDelta(2,q_min,q_max) + 2. * ( 3.*mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateDelta(1,q_min,q_max) - mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateDelta(3,q_min,q_max) ) )/(4.*computeGammaPrime(q_min, q_max));

}


P3CP::P3CP(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) : ThObservable(SM_i), mySM(SM_i) {  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double P3CP::computeThValue() {

    double q_min = getBinMin();
    double q_max = getBinMax();

    return - mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateDelta(11,q_min,q_max)/(4.*mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateSigma(3,q_min,q_max));

}


F_L::F_L(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) : ThObservable(SM_i), mySM(SM_i) {
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double F_L::computeThValue() {
    
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    double sigma0 = mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateSigma(0,q_min,q_max);
    double sigma2 = mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateSigma(2,q_min,q_max);
    
    return (3.*sigma0 - sigma2)/(((3.*sigma0 - sigma2) + 2.*(3.*mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateSigma(1,q_min,q_max) - mySM.getMyFlavour()->getMVll(meson, vectorM, lep).integrateSigma(3,q_min,q_max))));

}


M_1Prime::M_1Prime(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) : ThObservable(SM_i), mySM(SM_i) {  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double M_1Prime::computeThValue() {
    
    double q_min = getBinMin();
    
    return ( mySM.getMyFlavour()->getMVll(meson, vectorM, lep).H_V(1,q_min,0).abs2() + mySM.getMyFlavour()->getMVll(meson, vectorM, lep).H_V(2,q_min,0).abs2() - mySM.getMyFlavour()->getMVll(meson, vectorM, lep).H_A(1,q_min,0).abs2() - mySM.getMyFlavour()->getMVll(meson, vectorM, lep).H_A(2,q_min,0).abs2() )/
            ( 2.*( mySM.getMyFlavour()->getMVll(meson, vectorM, lep).H_V(1,q_min,0).abs2() + mySM.getMyFlavour()->getMVll(meson, vectorM, lep).H_V(2,q_min,0).abs2() + mySM.getMyFlavour()->getMVll(meson, vectorM, lep).H_A(1,q_min,0).abs2() + mySM.getMyFlavour()->getMVll(meson, vectorM, lep).H_A(2,q_min,0).abs2() ) );
}


M_2Prime::M_2Prime(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) : ThObservable(SM_i), mySM(SM_i) {  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double M_2Prime::computeThValue() {
    
    double q_min = getBinMin();
    
    return ( q_min/(2.*mySM.getMyFlavour()->getMVll(meson, vectorM, lep).Mlep*mySM.getMyFlavour()->getMVll(meson, vectorM, lep).Mlep)*( mySM.getMyFlavour()->getMVll(meson, vectorM, lep).H_P(q_min,0).abs2() + mySM.getMyFlavour()->getMVll(meson, vectorM, lep).beta(q_min)*mySM.getMyFlavour()->getMVll(meson, vectorM, lep).beta(q_min)*mySM.getMyFlavour()->getMVll(meson, vectorM, lep).H_S(q_min,0).abs2() ) + mySM.getMyFlavour()->getMVll(meson, vectorM, lep).H_V(0,q_min,0).abs2() - mySM.getMyFlavour()->getMVll(meson, vectorM, lep).H_A(0,q_min,0).abs2() )/
            ( mySM.getMyFlavour()->getMVll(meson, vectorM, lep).H_V(0,q_min,0).abs2() + mySM.getMyFlavour()->getMVll(meson, vectorM, lep).H_A(0,q_min,0).abs2() );  
}