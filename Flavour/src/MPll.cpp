/* 
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "MPll.h"
#include <gslpp.h>
#include <gslpp_complex.h>
#include <gsl/gsl_math.h>
#include <boost/bind.hpp>





MPll::MPll(const StandardModel& SM_i) : ThObservable(SM_i), mySM(SM_i){
}


MPll::~MPll() {
}

void MPll::updateParameters(){
    GF = mySM.getGF();
    ale=mySM.getAle();
    Me=mySM.getLeptons(StandardModel::ELECTRON).getMass();
    Mmu=mySM.getLeptons(StandardModel::MU).getMass();
    MB=mySM.getMesons(QCD::B_D).getMass();
    MK=mySM.getMesons(QCD::K_0).getMass();
    Mb=mySM.getQuarks(QCD::BOTTOM).getMass();    // add the PS b mass
    Ms=mySM.getQuarks(QCD::STRANGE).getMass();
    MW=mySM.Mw();
    lambda_t=mySM.computelamt_s();
    mu_b = mySM.getMub();
    width_Bd = mySM.getMesons(QCD::B_D).computeWidth();
    
    r_1_fplus = mySM.getr_1_fplus();
    r_2_fplus = mySM.getr_2_fplus();
    m_fit2_fplus = mySM.getm_fit2_fplus();
    r_1_fT = mySM.getr_1_fT();
    r_2_fT = mySM.getr_2_fT();
    m_fit2_fT = mySM.getm_fit2_fT();
    r_2_f0 = mySM.getr_2_f0();
    m_fit2_f0 = mySM.getm_fit2_f0();
    
    h_0=mySM.getReh_0() + gslpp::complex::i()*mySM.getImh_0();
    
    b=1.;                           //please check
    
    allcoeff = mySM.getMyFlavour()->ComputeCoeffBKstarll(mu_b);   //check the mass scale, scheme fixed to NDR
    allcoeffprime = mySM.getMyFlavour()->ComputeCoeffprimeBKstarll(mu_b);   //check the mass scale, scheme fixed to NDR

    C_7 = (*(allcoeff[LO]) + *(allcoeff[NLO]))(6);
    C_9 = (*(allcoeff[LO]) + *(allcoeff[NLO]))(8);
    C_10 = (*(allcoeff[LO]) + *(allcoeff[NLO]))(9);
    C_S = (*(allcoeff[LO]) + *(allcoeff[NLO]))(10);
    C_P = (*(allcoeff[LO]) + *(allcoeff[NLO]))(11);
    
    C_7p = (*(allcoeffprime[LO]) + *(allcoeffprime[NLO]))(6);
    C_9p = (*(allcoeffprime[LO]) + *(allcoeffprime[NLO]))(8);
    C_10p = (*(allcoeffprime[LO]) + *(allcoeffprime[NLO]))(9);
    C_Sp = (*(allcoeffprime[LO]) + *(allcoeffprime[NLO]))(10);
    C_Pp = (*(allcoeffprime[LO]) + *(allcoeffprime[NLO]))(11);
    
}

/*******************************************************************************
 * Transverse Form Factors                                                     *
 * ****************************************************************************/
double MPll::LCSR_fit1(double q2, double r_1, double r_2, double m_fit2){
    return r_1/( 1. - q2/m_fit2 ) + r_2/pow( ( 1. - q2/m_fit2 ) , 2.) ;
}



double MPll::LCSR_fit2(double q2, double r_2, double m_fit2){
    return r_2/( 1. - q2/m_fit2 ) ; 
}



double MPll::f_plus(double q2){
    return LCSR_fit1(q2, r_1_fplus, r_2_fplus, m_fit2_fplus);
}



double MPll::f_T(double q2){
    return LCSR_fit1(q2, r_1_fT, r_2_fT, m_fit2_fT);
}



double MPll::f_0(double q2){
    return LCSR_fit2(q2, r_2_f0, m_fit2_f0);
}



gslpp::complex MPll::V_L(double q2){
    return gslpp::complex::i()*f_plus(q2);
}



gslpp::complex MPll::V_R(double q2){
    return -V_L(q2);
}



gslpp::complex MPll::T_L(double q2){
    return gslpp::complex::i()* (2. * MB)/(MB + MK) * f_T(q2);
}



gslpp::complex MPll::T_R(double q2){
    return -T_L(q2);
}



double MPll::S_L(double q2){
    return ( 1 + Ms/Mb )/( 1 - Ms/Mb ) * ( MB*MB - MK*MK )/( sqrt(lambda(q2)) ) * f_0(q2);
}



double MPll::S_R(double q2){
    return -S_L(q2);
}



/*******************************************************************************
 * Helicity amplitudes                                                         *
 * ****************************************************************************/
complex MPll::N(){
    return -(4.*GF*MB*ale*lambda_t)/(sqrt(2.)*4.*M_PI);
}



gslpp::complex MPll::H_V(double q2, int bar) {
    gslpp::complex n;
    switch(bar){
        case 0:
            n = N();
            break;
        case 1:
            n = N().conjugate();
            break;
        default:
            std::stringstream out;
            out << bar;
            throw std::runtime_error("H_V: index " + out.str() + " not allowed for an Angular Coefficient");
    }
                    
    return -gslpp::complex::i()*n*( C_9*V_L(q2) + C_9p*V_R(q2)
            + MB*MB/q2*( 2*Mb/MB*( C_7*T_L(q2) + C_7p*T_R(q2) ) - 16*M_PI*M_PI*h_0 ) );
}



gslpp::complex MPll::H_A(double q2, int bar) {
    gslpp::complex n;
    switch(bar){
        case 0:
            n = N();
            break;
        case 1:
            n = N().conjugate();
            break;
        default:
            std::stringstream out;
            out << bar;
            throw std::runtime_error("H_A: index " + out.str() + " not allowed for an Angular Coefficient");
    }
     
    return -gslpp::complex::i()*n*( C_10*V_L(q2) + C_10p*V_R(q2) );
}



gslpp::complex MPll::H_S(double q2, int bar) {
    gslpp::complex n;
    switch(bar){
        case 0:
            n = N();
            break;
        case 1:
            n = N().conjugate();
            break;
        default:
            std::stringstream out;
            out << bar;
            throw std::runtime_error("H_S: index " + out.str() + " not allowed for an Angular Coefficient");
    }
     
    return gslpp::complex::i()*n*Mb/MW*( C_S*S_L(q2) + C_Sp*S_R(q2) );
}



gslpp::complex MPll::H_P(double q2, int bar, double Mlep) {
    gslpp::complex n;
    switch(bar){
        case 0:
            n = N();
            break;
        case 1:
            n = N().conjugate();
            break;
        default:
            std::stringstream out;
            out << bar;
            throw std::runtime_error("H_S: index " + out.str() + " not allowed for an Angular Coefficient");
    }
     
    return gslpp::complex::i()*n*( Mb/MW*( C_P*S_L(q2) + C_Pp*S_R(q2) ) 
            + 2.*Mlep*Mb/q2*( C_10*( S_L(q2) - Ms/Mb*S_R(q2) ) + C_10p*( S_R(q2) - Ms/Mb*S_L(q2) ) ) );
}



/*******************************************************************************
 * Angular coefficients                                                         *
 * ****************************************************************************/
double MPll::k2(double q2) {
    return (pow(MB,4.) + q2*q2 + pow(MK,4.) -2.*MK*MK*q2 -2.*MB*MB*(q2 + MK*MK))/(4.*MB*MB);
}



double MPll::beta(double q2, double Mlep) {
    return sqrt(1-4.*Mlep*Mlep/q2);
}



double MPll::lambda(double q2) {
    return 4.*MB*MB*k2(q2);
}



double MPll::F(double q2, double b_i, double Mlep) {
    return sqrt(lambda(q2))*beta(q2,Mlep)*q2*b_i/(96.*M_PI*M_PI*M_PI*MB*MB*MB);
}



double MPll::I(int i, double q2, int bar, double Mlep) {

    double Mlep2 = Mlep*Mlep;
    double beta2 = beta(q2,Mlep)*beta(q2,Mlep);
    

    switch (i){
        case 0: // I1c
            return F(q2,b,Mlep)*( ( H_V(q2,bar).abs2() + H_A(q2,bar).abs2() )/2.  +  H_P(q2,bar,Mlep).abs2()  +  2.*Mlep2/q2*( H_V(q2,bar).abs2() 
                    - H_A(q2,bar).abs2() )  + beta2*H_S(q2,bar).abs2() );
        case 2: // I2c
            return -F(q2,b,Mlep)*beta2/2.*( H_V(q2,bar).abs2() + H_A(q2,bar).abs2() );
        case 8: // I6c
            return 2.*F(q2,b,Mlep)*beta(q2,Mlep)*Mlep/sqrt(q2)*( H_S(q2,bar).conjugate()*H_V(q2,bar) ).real();
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("I: index " + out.str() + " not implemented");
    }
}



double MPll::Sigma(int i, double q2, double Mlep) {
    return (I(i, q2, 0, Mlep) + I(i, q2, 1, Mlep))/2;
}



double MPll::Delta(int i, double q2, double Mlep) {
    return (I(i, q2, 0, Mlep) - I(i, q2, 1, Mlep))/2;
}

/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/


BR_MPll_e::BR_MPll_e(const StandardModel& SM_i) : MPll(SM_i) {  
}

double BR_MPll_e::computeBR_MPll_e(double qmin, double qmax) {
    
    double q_min = qmin;
    double q_max = qmax;
    
    F1 = convertToGslFunction( boost::bind( &MPll::getSigma2e, &(*this), _1 ) );
    gsl_integration_workspace * w_sigma2 = gsl_integration_workspace_alloc (50);
    gsl_integration_qags (&F1, q_min, q_max, 1.e-5, 1.e-3, 50, w_sigma2, &avaSigma2, &errSigma2);
    gsl_integration_workspace_free (w_sigma2);

    F2 = convertToGslFunction( boost::bind( &MPll::getSigma0e, &(*this), _1 ) );
    gsl_integration_workspace * w_sigma0 = gsl_integration_workspace_alloc (50);
    gsl_integration_qags (&F2, q_min, q_max, 1.e-5, 1.e-3, 50, w_sigma0, &avaSigma0, &errSigma0);
    gsl_integration_workspace_free (w_sigma0);
    
    return (3.*avaSigma0 - avaSigma2)/(4. * width_Bd);
}

double BR_MPll_e::computeThValue(){
    updateParameters();
    double q_min = getBinMin();
    double q_max = getBinMax();
    return computeBR_MPll_e(q_min, q_max);
}

BR_MPll_mu::BR_MPll_mu(const StandardModel& SM_i) : BR_MPll_e(SM_i) {  
}

double BR_MPll_mu::computeBR_MPll_mu(double qmin, double qmax) {
    
    double q_min = qmin;
    double q_max = qmax;
    
    F1 = convertToGslFunction( boost::bind( &MPll::getSigma2mu, &(*this), _1 ) );
    gsl_integration_workspace * w_sigma2 = gsl_integration_workspace_alloc (50);
    gsl_integration_qags (&F1, q_min, q_max, 1.e-5, 1.e-3, 50, w_sigma2, &avaSigma2, &errSigma2);
    gsl_integration_workspace_free (w_sigma2);

    F2 = convertToGslFunction( boost::bind( &MPll::getSigma0mu, &(*this), _1 ) );
    gsl_integration_workspace * w_sigma0 = gsl_integration_workspace_alloc (50);
    gsl_integration_qags (&F2, q_min, q_max, 1.e-5, 1.e-3, 50, w_sigma0, &avaSigma0, &errSigma0);
    gsl_integration_workspace_free (w_sigma0);
    
    return (3.*avaSigma0 - avaSigma2)/(4. * width_Bd);
}

double BR_MPll_mu::computeThValue(){
    updateParameters();
    double q_min = getBinMin();
    double q_max = getBinMax();
    return computeBR_MPll_mu(q_min, q_max);
}


R_MPll::R_MPll(const StandardModel& SM_i) : BR_MPll_mu(SM_i) {  
}

double R_MPll::computeThValue() {
    updateParameters();
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return computeBR_MPll_e(q_min, q_max)/computeBR_MPll_mu(q_min, q_max);
}

ACP_e::ACP_e(const StandardModel& SM_i) : BR_MPll_e(SM_i){
}

double ACP_e::computeThValue() {
    
    updateParameters();
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    F1 = convertToGslFunction( boost::bind( &BR_MPll_e::getDelta2e, &(*this), _1 ) );
    gsl_integration_workspace * w_delta2 = gsl_integration_workspace_alloc (50);
    gsl_integration_qags (&F1, q_min, q_max, 1.e-5, 1.e-3, 50, w_delta2, &avaDelta2, &errDelta2);
    gsl_integration_workspace_free (w_delta2);

    F2 = convertToGslFunction( boost::bind( &BR_MPll_e::getDelta0e, &(*this), _1 ) );
    gsl_integration_workspace * w_delta0 = gsl_integration_workspace_alloc (50);
    gsl_integration_qags (&F2, q_min, q_max, 1.e-5, 1.e-3, 50, w_delta0, &avaDelta0, &errDelta0);
    gsl_integration_workspace_free (w_delta0);

    return (3.*avaDelta0 - avaDelta2)/(4.*computeBR_MPll_e(q_min, q_max)* width_Bd);

}

ACP_mu::ACP_mu(const StandardModel& SM_i) : BR_MPll_mu(SM_i){
}

double ACP_mu::computeThValue() {
    
    updateParameters();
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    F1 = convertToGslFunction( boost::bind( &BR_MPll_mu::getDelta2e, &(*this), _1 ) );
    gsl_integration_workspace * w_delta2 = gsl_integration_workspace_alloc (50);
    gsl_integration_qags (&F1, q_min, q_max, 1.e-5, 1.e-3, 50, w_delta2, &avaDelta2, &errDelta2);
    gsl_integration_workspace_free (w_delta2);

    F2 = convertToGslFunction( boost::bind( &BR_MPll_mu::getDelta0e, &(*this), _1 ) );
    gsl_integration_workspace * w_delta0 = gsl_integration_workspace_alloc (50);
    gsl_integration_qags (&F2, q_min, q_max, 1.e-5, 1.e-3, 50, w_delta0, &avaDelta0, &errDelta0);
    gsl_integration_workspace_free (w_delta0);

    return (3.*avaDelta0 - avaDelta2)/(4.*computeBR_MPll_mu(q_min, q_max)* width_Bd);

}