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
    ale = mySM.getAle();
    MM = mySM.getMesons(meson).getMass();
    MM2 = MM * MM;
    MV = mySM.getMesons(vectorM).getMass();
    Mb = mySM.getQuarks(QCD::BOTTOM).getMass();    // add the PS b mass
    Ms = mySM.getQuarks(QCD::STRANGE).getMass();
    MW = mySM.Mw();
    lambda_t = mySM.computelamt_s();
    mu_b = mySM.getMub();
    width = mySM.getMesons(meson).computeWidth();
    lambda = MM2 - pow(MV,2.);
    
    switch(vectorM){
        case StandardModel::K_star :
            r_1T1=mySM.getr_1T1();
            r_2T1=mySM.getr_2T1();
            m_RT1=mySM.getm_RT1();
            m_fit2T1=mySM.getm_fit2T1();
            
            break;
        case StandardModel::PHI :
            r_1T1=mySM.getr_1T1phi();
            r_2T1=mySM.getr_2T1phi();
            m_RT1=mySM.getm_RT1phi();
            m_fit2T1=mySM.getm_fit2T1phi();
            
            break;
        default:
            std::stringstream out;
            out << vectorM;
            throw std::runtime_error("MVgamma: vector " + out.str() + " not implemented");
    }
    
    
    h[0]=mySM.geth_p();    //h_plus
    h[1]=mySM.geth_m();    //h_minus
    
    allcoeff = mySM.getMyFlavour()->ComputeCoeffBKstarll(mu_b);   //check the mass scale, scheme fixed to NDR
    allcoeffprime = mySM.getMyFlavour()->ComputeCoeffprimeBKstarll(mu_b);   //check the mass scale, scheme fixed to NDR
    
    C_7 = (*(allcoeff[LO]))(6) + (*(allcoeff[NLO]))(6);
    C_7p = (*(allcoeffprime[LO]))(6) + (*(allcoeffprime[NLO]))(6);
    
}


/*******************************************************************************
 * Form Factor                                                     *
 * ****************************************************************************/
double MVgamma::T_1(){
    return mySM.getMyFlavour()->getMVll(meson, vectorM, StandardModel::MU)->LCSR_fit1(0., r_1T1, r_2T1, pow(m_RT1, 2.), m_fit2T1);
}



/*******************************************************************************
 * Helicity amplitudes                                                         *
 * ****************************************************************************/
complex MVgamma::H_V() {
    return lambda_t * ((C_7 - C_7p) *T_1() * lambda / MM2 - MM/(2*Mb)*16*M_PI*M_PI*(h[0] + h[1]));
}

complex MVgamma::H_V_bar() {
    return lambda_t.conjugate() * ((C_7 - C_7p) *T_1() * lambda / MM2 - MM/(2*Mb)*16*M_PI*M_PI*(h[0] + h[1]));
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
    
    return ale * pow(GF * Mb / (4 * M_PI * M_PI), 2.) * MM * lambda /(4. * width) * (H_V().abs2() + H_V_bar().abs2());
}

ACP_MVgamma::ACP_MVgamma(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i) : mySM(SM_i), MVgamma(SM_i, meson_i, vector_i) {
    meson = meson_i;
    vectorM = vector_i;
}

double ACP_MVgamma::computeThValue(){
    updateParameters();
    
    return (H_V().abs2() - H_V_bar().abs2()) / (H_V().abs2() + H_V_bar().abs2());
}