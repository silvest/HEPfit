/* 
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Flavour.h"
#include "MVgamma.h"
#include <gslpp_complex.h>
#include <boost/bind.hpp>


MVgamma::MVgamma(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i) : ThObservable(SM_i), mySM(SM_i){
    meson = meson_i;
    vectorM = vector_i;
}


MVgamma::~MVgamma() {
/** Check to see if GSL pointers are released!!*/
}

void MVgamma::updateParameters(){
    GF = mySM.getGF();
    ale=mySM.getAle();
    MM=mySM.getMesons(meson).getMass();
    MV=mySM.getMesons(vectorM).getMass();
    Mb=mySM.getQuarks(QCD::BOTTOM).getMass();    // add the PS b mass
    Ms=mySM.getQuarks(QCD::STRANGE).getMass();
    MW=mySM.Mw();
    lambda_t=mySM.computelamt_s();
    mu_b = mySM.getMub();
    width = mySM.getMesons(meson).computeWidth();
    lambda = pow(MM,4.) + pow(MV,4.) -2.*MM*MM*MV*MV;
    
    switch(vectorM){
        case StandardModel::K_star :
            r_1T1=mySM.getr_1T1();
            r_2T1=mySM.getr_2T1();
            m_RT1=mySM.getm_RT1();
            m_fit2T1=mySM.getm_fit2T1();
            r_2T2=mySM.getr_2T2();
            m_fit2T2=mySM.getm_fit2T2();
            
            break;
        default:
            std::stringstream out;
            out << vectorM;
            throw std::runtime_error("MVgamma: vector " + out.str() + " not implemented");
    }
    
    
    h[0]=mySM.geth_p();    //h_plus
    h[1]=mySM.geth_m();  //h_minus
    
    allcoeff = mySM.getMyFlavour()->ComputeCoeffBKstarll(mu_b);   //check the mass scale, scheme fixed to NDR
    allcoeffprime = mySM.getMyFlavour()->ComputeCoeffprimeBKstarll(mu_b);   //check the mass scale, scheme fixed to NDR
    
    C_7 = (*(allcoeff[LO]))(6) + (*(allcoeff[NLO]))(6);
    C_7p = (*(allcoeffprime[LO]))(6) + (*(allcoeffprime[NLO]))(6);
    
}


/*******************************************************************************
 * Transverse Form Factors                                                     *
 * ****************************************************************************/
double MVgamma::T_1(){
    return mySM.getMyFlavour()->getMVll(meson, vectorM, StandardModel::MU)->LCSR_fit1(0., r_1T1, r_2T1, pow(m_RT1, 2.), m_fit2T1);
}



double MVgamma::T_2(){
    return mySM.getMyFlavour()->getMVll(meson, vectorM, StandardModel::MU)->LCSR_fit3(0., r_2T2, m_fit2T2);
}



double MVgamma::T_L(int i){
    switch (i){
        case 1:
            return (MM*MM - MV*MV) / ( 2.*MM*MM ) * T_2() - sqrt(lambda) / ( 2.*MM*MM ) * T_1();
        case 2:
            return (MM*MM - MV*MV) / ( 2.*MM*MM ) * T_2() + sqrt(lambda) / ( 2.*MM*MM ) * T_1();
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("T_L: index " + out.str() + " not implemented");
    }
}



double MVgamma::T_R(int i){
    return -T_L(3-i);
}



/*******************************************************************************
 * Helicity amplitudes                                                         *
 * ****************************************************************************/
double MVgamma::H_V_plus2() {
    return (C_7*T_L(1) + C_7p*T_R(1) - MM/(2*Mb)*16*M_PI*M_PI*h[0]).abs2();
}



double MVgamma::H_V_minus2() {
    return (C_7*T_L(2) + C_7p*T_R(2) - MM/(2*Mb)*16*M_PI*M_PI*h[1]).abs2();
}


/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/


BR_MVgamma::BR_MVgamma(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i) : mySM(SM_i), MVgamma(SM_i, meson_i, vector_i) {
    meson = meson_i;
    vectorM = vector_i;
}

double BR_MVgamma::computeThValue(){
    updateParameters();
    
    return ale*GF*GF*Mb*Mb*MM*sqrt(lambda)*lambda_t.abs2()/(32.*pow(M_PI,4.)*width)*(H_V_plus2() + H_V_minus2());
}