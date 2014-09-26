/* 
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "BKstarll.h"
#include <gslpp.h>
#include <gslpp_complex.h>
#include <gsl/gsl_math.h>
#include <boost/bind.hpp>





BKstarll::BKstarll(const StandardModel& SM_i, StandardModel::lepton lep_i) : ThObservable(SM_i), mySM(SM_i) {
    lep = lep_i;
}


BKstarll::~BKstarll() {
}

void BKstarll::updateParameters(){
    GF = mySM.getGF();
    ale=mySM.getAle();
    Mm=mySM.getLeptons(lep).getMass();
    MB=mySM.getMesons(QCD::B_D).getMass();
    MKstar=mySM.getMesons(QCD::K_star).getMass();
    Mb=mySM.getQuarks(QCD::BOTTOM).getMass();    // add the PS b mass
    Ms=mySM.getQuarks(QCD::STRANGE).getMass();
    MW=mySM.Mw();
    lambda_t=mySM.computelamt_s();
    mu_b = mySM.getMub();
    width_Bd = mySM.getMesons(QCD::B_D).computeWidth();
    
    a_0V=mySM.geta_0V();
    a_1V=mySM.geta_1V();
    dmV=mySM.getdmV();
    a_0A0=mySM.geta_0A0();
    a_1A0=mySM.geta_1A0();
    dmA0=mySM.getdmA0();
    a_0A1=mySM.geta_0A1();
    a_1A1=mySM.geta_1A1();
    dmA1=mySM.getdmA1();
    a_0A12=mySM.geta_0A12();
    a_1A12=mySM.geta_1A12();
    dmA12=mySM.getdmA12();
    a_0T1=mySM.geta_0T1();
    a_1T1=mySM.geta_1T1();
    dmT1=mySM.getdmT1();
    a_0T2=mySM.geta_0T2();
    a_1T2=mySM.geta_1T2();
    dmT2=mySM.getdmT2();
    a_0T23=mySM.geta_0T23();
    a_1T23=mySM.geta_1T23();
    dmT23=mySM.getdmT23();
    
    r_1V=mySM.getr_1V();
    r_2V=mySM.getr_2V();
    m_RV=mySM.getm_RV();
    m_fit2V=mySM.getm_fit2V();
    r_1A0=mySM.getr_1A0();
    r_2A0=mySM.getr_2A0();
    m_RA0=mySM.getm_RA0();
    m_fit2A0=mySM.getm_fit2A0();
    r_2A1=mySM.getr_2A1();
    m_fit2A1=mySM.getm_fit2A1();
    r_1A2=mySM.getr_1A2();
    r_2A2=mySM.getr_2A2();
    m_fit2A2=mySM.getm_fit2A2();
    r_1T1=mySM.getr_1T1();
    r_2T1=mySM.getr_2T1();
    m_RT1=mySM.getm_RT1();
    m_fit2T1=mySM.getm_fit2T1();
    r_2T2=mySM.getr_2T2();
    m_fit2T2=mySM.getm_fit2T2();
    r_1T3t=mySM.getr_1T3t();
    r_2T3t=mySM.getr_2T3t();
    m_fit2T3t=mySM.getm_fit2T3t();
    
    h[0]=mySM.getabsh_0() * exp(gslpp::complex::i()*mySM.getargh_0());
    h[1]=mySM.getabsh_plus() * exp(gslpp::complex::i()*mySM.getargh_plus());
    h[2]=mySM.getabsh_minus() * exp(gslpp::complex::i()*mySM.getargh_minus());
    
    b=1.;                           //please check
    
    allcoeff = mySM.getMyFlavour()->ComputeCoeffBKstarll(mu_b);   //check the mass scale, scheme fixed to NDR
    allcoeffprime = mySM.getMyFlavour()->ComputeCoeffprimeBKstarll(mu_b);   //check the mass scale, scheme fixed to NDR

}


/*******************************************************************************
 * Transverse Form Factors                                                     *
 * ****************************************************************************/
double BKstarll::LCSR_fit1(double q2, double r_1, double r_2, double m_R2, double m_fit2){
    return r_1/( 1 - q2/m_R2 ) + r_2/( 1 - q2/m_fit2 ) ;
}



double BKstarll::LCSR_fit2(double q2, double r_1, double r_2, double m_fit2){
    return r_1/( 1 - q2/m_fit2 ) + r_2/pow( ( 1 - q2/m_fit2 ) ,2) ;

}



double BKstarll::LCSR_fit3(double q2, double r_2, double m_fit2){
    return r_2/( 1 - q2/m_fit2 ) ; 
}



double BKstarll::z(double q2){
    double t_0=12;
    double t_p=pow(MB + MKstar,2);
    return ( sqrt(t_p - q2) - sqrt(t_p - t_0) ) / ( sqrt(t_p - q2) + sqrt(t_p - t_0) );
}



double BKstarll::lat_fit(double q2, double a_0, double a_1, double dm){
    return 1 / (1 - q2/pow(MB + dm,2)) * ( a_0 + a_1*z(q2) );
}



double BKstarll::V(double q2){
    if (q2<CUTOFF)
        return LCSR_fit1(q2, r_1V, r_2V, pow(m_RV, 2), m_fit2V);
    else
        return lat_fit(q2, a_0V, a_1V, dmV);
}



double BKstarll::A_0(double q2){
    if (q2<CUTOFF)
        return LCSR_fit1(q2, r_1A0, r_2A0, pow(m_RA0, 2), m_fit2A0);
    else
        return lat_fit(q2, a_0A0, a_1A0, dmA0);
}



double BKstarll::A_1(double q2){
    if (q2<CUTOFF)
        return LCSR_fit3(q2, r_2A1, m_fit2A1);
    else
        return lat_fit(q2, a_0A1, a_1A1, dmA1);
}



double BKstarll::A_2(double q2){
    return LCSR_fit2(q2, r_1A2, r_2A2, m_fit2A2);
}



double BKstarll::T_1(double q2){
    if (q2<CUTOFF)
        return LCSR_fit1(q2, r_1T1, r_2T1, pow(m_RT1, 2), m_fit2T1);
    else
        return lat_fit(q2, a_0T1, a_1T1, dmT1);
}



double BKstarll::T_2(double q2){
    if (q2<CUTOFF)
        return LCSR_fit3(q2, r_2T2, m_fit2T2);
    else
        return lat_fit(q2, a_0T2, a_1T2, dmT2);
}



double BKstarll::T_3tilde(double q2){
    return LCSR_fit2(q2, r_1T3t, r_2T3t, m_fit2T3t);
}



double BKstarll::T_3(double q2){
    if (q2 < 2.) return (0.178168 - 0.202)/2. * q2 + 0.202;
    else return (MB*MB - MKstar*MKstar)/q2*(T_3tilde(q2) - T_2(q2));
}



double BKstarll::V_L(int i, double q2){
    switch (i){
        case 0:
            if (q2 < CUTOFF)
                return 1. / ( 4*MKstar*MB*(MB + MKstar)*sqrt(q2) ) * ( pow((MB + MKstar),2)*(MB*MB - q2 - MKstar*MKstar)*A_1(q2) - lambda(q2)*A_2(q2) );
            else
                return 4*MKstar/sqrt(q2)*lat_fit(q2, a_0A12, a_1A12, dmA12);
        case 1:
            return 1./2. * ( ( 1. + MKstar/MB)*A_1(q2) - sqrt(lambda(q2))/ ( MB* (MB + MKstar) ) * V(q2) );
        case 2:
            return 1./2. * ( ( 1. + MKstar/MB)*A_1(q2) + sqrt(lambda(q2))/ ( MB* (MB + MKstar) ) * V(q2) );
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("I: index " + out.str() + " not implemented");
    }
}



double BKstarll::V_R(int i, double q2){
    if (i != 0) i=3-i;
    return -V_L(i,q2);
}



double BKstarll::T_L(int i, double q2){
    switch (i){
        case 0:
            if (q2 < CUTOFF)
                return sqrt(q2)/(4.*MB*MB*MKstar) * ( ( MB*MB+ 3.*MKstar*MKstar - q2 ) * T_2(q2) - lambda(q2) / (MB*MB - MKstar*MKstar) * T_3(q2) );
            else
                return 2*sqrt(q2)*MKstar/(MB + MKstar)*lat_fit(q2, a_0T23, a_1T23, dmT23);
        case 1:
            return (MB*MB - MKstar*MKstar) / ( 2.*MB*MB ) * T_2(q2) - sqrt(lambda(q2)) / ( 2.*MB*MB ) * T_1(q2);
        case 2:
            return (MB*MB - MKstar*MKstar) / ( 2.*MB*MB ) * T_2(q2) + sqrt(lambda(q2)) / ( 2.*MB*MB ) * T_1(q2);
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("I: index " + out.str() + " not implemented");
    }
}



double BKstarll::T_R(int i, double q2){
    if (i != 0) i=3-i;
    return -T_L(i,q2);
}



double BKstarll::S_L(double q2){
    return -sqrt(lambda(q2))/ ( 2*MB*(Mb + Ms) ) *A_0(q2);
}



double BKstarll::S_R(double q2){
    return -S_L(q2);
}



/*******************************************************************************
 * Helicity amplitudes                                                         *
 * ****************************************************************************/
complex BKstarll::N(){
    return -(4*GF*MB*ale*lambda_t)/(sqrt(2)*4*M_PI);
}



gslpp::complex BKstarll::H_V(int i, double q2, int bar) {
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
                    
    return -gslpp::complex::i()*n*( (*(allcoeff[LO]) + *(allcoeff[NLO]))(8)*V_L(i,q2)
            + (*(allcoeffprime[LO]) + *(allcoeffprime[NLO]))(8)*V_R(i,q2)
            + MB*MB/q2*( 2*Mb/MB*( (*(allcoeff[LO]) + *(allcoeff[NLO]))(6)*T_L(i,q2)
            + (*(allcoeffprime[LO]) + *(allcoeffprime[NLO]))(6)*T_R(i,q2) ) - 16*M_PI*M_PI*h[i] ) );
}



gslpp::complex BKstarll::H_A(int i, double q2, int bar) {
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
     
    return -gslpp::complex::i()*n*( (*(allcoeff[LO]) + *(allcoeff[NLO]))(9)*V_L(i,q2) 
            + (*(allcoeffprime[LO]) + *(allcoeffprime[NLO]))(9)*V_R(i,q2) );
}



gslpp::complex BKstarll::H_S(double q2, int bar) {
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
     
    return gslpp::complex::i()*n*Mb/MW*( (*(allcoeff[LO]) + *(allcoeff[NLO]))(10)*S_L(q2) + 
            (*(allcoeffprime[LO]) + *(allcoeffprime[NLO]))(10)*S_R(q2) );
}



gslpp::complex BKstarll::H_P(double q2, int bar) {
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
     
    return gslpp::complex::i()*n*( Mb/MW*( (*(allcoeff[LO]) + *(allcoeff[NLO]))(11)*S_L(q2) 
            + (*(allcoeffprime[LO]) + *(allcoeffprime[NLO]))(11)*S_R(q2) ) 
            + 2*Mm*Mb/q2*( (*(allcoeff[LO]) + *(allcoeff[NLO]))(9)*( S_L(q2) - Ms/Mb*S_R(q2) ) 
            + (*(allcoeffprime[LO]) + *(allcoeffprime[NLO]))(9)*( S_R(q2) - Ms/Mb*S_L(q2) ) ) );
}



/*******************************************************************************
 * Angular coefficients                                                         *
 * ****************************************************************************/
double BKstarll::k2(double q2) {
    return (pow(MB,4) + q2*q2 + pow(MKstar,4) -2*MKstar*MKstar*q2 -2*MB*MB*(q2 + MKstar*MKstar))/(4*MB*MB);
}



double BKstarll::beta(double q2) {
    return sqrt(1-4*Mm*Mm/q2);
}



double BKstarll::lambda(double q2) {
    return 4*MB*MB*k2(q2);
}



double BKstarll::F(double q2, double b) {
    return sqrt(lambda(q2))*beta(q2)*q2*b/(96*M_PI*M_PI*M_PI*MB*MB*MB);   
}



double BKstarll::I(int i, double q2, int bar) {

    double Mm2 = Mm*Mm;
    double beta2 = beta(q2)*beta(q2);
    
    if (i < 10) {
        switch (i){
            case 0:
                return F(q2,b)*( ( H_V(0,q2,bar).abs2() + H_A(0,q2,bar).abs2() )/2  +  H_P(q2,bar).abs2()  +  2*Mm2/q2*( H_V(0,q2,bar).abs2() 
                        - H_A(0,q2,bar).abs2() )  + beta2*H_S(q2,bar).abs2() );
            case 1:
                return F(q2,b)*( (beta2 + 2.)/8.*( H_V(1,q2,bar).abs2() + H_V(2,q2,bar).abs2() + H_A(1,q2,bar).abs2() + H_A(2,q2,bar).abs2() )  +
                                Mm2/q2*( H_V(1,q2,bar).abs2() + H_V(2,q2,bar).abs2() - H_A(1,q2,bar).abs2() - H_A(2,q2,bar).abs2() ) );
            case 2:
                return -F(q2,b)*beta2/2*( H_V(0,q2,bar).abs2() + H_A(0,q2,bar).abs2() );
            case 3:
                return F(q2,b)*beta2/8*( H_V(1,q2,bar).abs2() + H_V(2,q2,bar).abs2()  +  H_A(1,q2,bar).abs2() + H_A(2,q2,bar).abs2() );
            case 4:
                return -F(q2,b)/2*( ( H_V(1,q2,bar)*H_V(2,q2,bar).conjugate() ).real()  +  ( H_A(1,q2,bar)*H_A(2,q2,bar).conjugate() ).real() );
            case 5:
                return F(q2,b)*beta2/4*( ( (H_V(2,q2,bar) + H_V(1,q2,bar))*H_V(0,q2,bar).conjugate() ).real()  +  ( (H_A(2,q2,bar) + H_A(1,q2,bar))*H_A(0,q2,bar).conjugate() ).real() );
            case 6:
                return F(q2,b)*( beta(q2)/2*( ( (H_V(2,q2,bar) - H_V(1,q2,bar))*H_A(0,q2,bar).conjugate() ).real()  +  ( (H_A(2,q2,bar) - H_A(1,q2,bar))*H_V(0,q2,bar).conjugate() ).real() )  -
                                beta(q2)*Mm/sqrt(q2)*( H_S(q2,bar).conjugate()*(H_V(1,q2,bar) + H_V(2,q2,bar)) ).real() );
            case 7:
                return F(q2,b)*beta(q2)*( H_V(2,q2,bar)*(H_A(2,q2,bar).conjugate()) - H_V(1,q2,bar)*(H_A(1,q2,bar).conjugate()) ).real();
            case 8:
                return 2*F(q2,b)*beta(q2)*Mm/sqrt(q2)*( H_S(q2,bar).conjugate()*H_V(0,q2,bar) ).real();
            case 9:
                return F(q2,b)*( beta(q2)/2*( ( (H_V(2,q2,bar) + H_V(1,q2,bar))*H_A(0,q2,bar).conjugate() ).imag()  +  ( (H_A(2,q2,bar) + H_A(1,q2,bar))*H_V(0,q2,bar).conjugate() ).imag() )  -
                                beta(q2)*Mm/sqrt(q2)*( H_S(q2,bar).conjugate()*(H_V(2,q2,bar) - H_V(1,q2,bar)) ).imag() );
            
            default:
                std::stringstream out;
                out << i;
                throw std::runtime_error("I: index " + out.str() + " not implemented");
        }
    } else {
        int j = i - 10;
        switch (j) {
            case 0:
                return F(q2,b)*beta2/4*( ( (H_V(2,q2,bar) - H_V(1,q2,bar))*H_V(0,q2,bar).conjugate() ).imag()  +  ( (H_A(2,q2,bar) - H_A(1,q2,bar))*H_A(0,q2,bar).conjugate() ).imag() );
            case 1:
                return F(q2,b)*beta2/2*( ( H_V(1,q2,bar)*H_V(2,q2,bar).conjugate() ).imag()  +  ( H_A(1,q2,bar)*H_A(2,q2,bar).conjugate() ).imag() );
            default:
                std::stringstream out;
                out << j + 10;
                throw std::runtime_error("I: index " + out.str() + " not implemented");
        }
    }
}



double BKstarll::Sigma(int i, double q2) {
    return (I(i, q2,0) + I(i, q2,1))/2;
}



double BKstarll::Delta(int i, double q2) {
    return (I(i, q2,0) - I(i, q2,1))/2;
}

/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/



P_1::P_1(const StandardModel& SM_i, StandardModel::lepton lep_i) : BKstarll(SM_i, lep_i) {
}

double P_1::computeThValue() {
    
    updateParameters();
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    
    gsl_function F1 = convertToGslFunction( boost::bind( &BKstarll::getSigma3, &(*this), _1 ) );
    gsl_function F2 = convertToGslFunction( boost::bind( &BKstarll::getSigma4, &(*this), _1 ) );
    
    double avaSigma3, errSigma3, avaSigma4, errSigma4;
    gsl_integration_workspace * w_sigma3 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F1, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma3, &avaSigma3, &errSigma3);
    gsl_integration_workspace_free (w_sigma3);

    
    gsl_integration_workspace * w_sigma4 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F2, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma4, &avaSigma4, &errSigma4);
    gsl_integration_workspace_free (w_sigma4);

    return avaSigma4/(2.*avaSigma3);
}


P_2::P_2(const StandardModel& SM_i, StandardModel::lepton lep_i) : BKstarll(SM_i, lep_i) {  
}

double P_2::computeThValue() {
    
    updateParameters();
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    
    gsl_function F1 = convertToGslFunction( boost::bind( &BKstarll::getSigma7, &(*this), _1 ) );
    gsl_function F2 = convertToGslFunction( boost::bind( &BKstarll::getSigma3, &(*this), _1 ) );
    
    double avaSigma7, errSigma7, avaSigma3, errSigma3;
    gsl_integration_workspace * w_sigma7 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F1, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma7, &avaSigma7, &errSigma7);
    gsl_integration_workspace_free (w_sigma7);

    gsl_integration_workspace * w_sigma3 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F2, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma3, &avaSigma3, &errSigma3);
    gsl_integration_workspace_free (w_sigma3);
    
    return avaSigma7/(8.*avaSigma3);
}


P_3::P_3(const StandardModel& SM_i, StandardModel::lepton lep_i) : BKstarll(SM_i, lep_i) {  
}

double P_3::computeThValue() {
    
    updateParameters();
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    
    gsl_function F1 = convertToGslFunction( boost::bind( &BKstarll::getSigma11, &(*this), _1 ) );
    gsl_function F2 = convertToGslFunction( boost::bind( &BKstarll::getSigma3, &(*this), _1 ) );
    
    double avaSigma11, errSigma11, avaSigma3, errSigma3;
    gsl_integration_workspace * w_sigma11 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F1, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma11, &avaSigma11, &errSigma11);
    gsl_integration_workspace_free (w_sigma11);
    
    gsl_integration_workspace * w_sigma3 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F2, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma3, &avaSigma3, &errSigma3);
    gsl_integration_workspace_free (w_sigma3);
    
    return -avaSigma11/(4.*avaSigma3);
}


P_4Prime::P_4Prime(const StandardModel& SM_i, StandardModel::lepton lep_i) : BKstarll(SM_i, lep_i) {  
}

double P_4Prime::computeThValue() {
    
    updateParameters();
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    
    gsl_function F1 = convertToGslFunction( boost::bind( &BKstarll::getSigma5, &(*this), _1 ) );
    gsl_function F2 = convertToGslFunction( boost::bind( &BKstarll::getSigma2, &(*this), _1 ) );
    gsl_function F3 = convertToGslFunction( boost::bind( &BKstarll::getSigma3, &(*this), _1 ) );
    
    double avaSigma5, errSigma5, avaSigma2, errSigma2, avaSigma3, errSigma3;
    gsl_integration_workspace * w_sigma5 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F1, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma5, &avaSigma5, &errSigma5);
    gsl_integration_workspace_free (w_sigma5);
    
    gsl_integration_workspace * w_sigma2 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F2, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma2, &avaSigma2, &errSigma2);
    gsl_integration_workspace_free (w_sigma2);
    
    gsl_integration_workspace * w_sigma3 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F3, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma3, &avaSigma3, &errSigma3);
    gsl_integration_workspace_free (w_sigma3);
    
    return avaSigma5/sqrt(-avaSigma2*avaSigma3);
   
}


P_5Prime::P_5Prime(const StandardModel& SM_i, StandardModel::lepton lep_i) : BKstarll(SM_i, lep_i) {  
}

double P_5Prime::computeThValue() {
    
    updateParameters();
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    
    gsl_function F1 = convertToGslFunction( boost::bind( &BKstarll::getSigma6, &(*this), _1 ) );
    gsl_function F2 = convertToGslFunction( boost::bind( &BKstarll::getSigma2, &(*this), _1 ) );
    gsl_function F3 = convertToGslFunction( boost::bind( &BKstarll::getSigma3, &(*this), _1 ) );
    
    double avaSigma6, errSigma6, avaSigma2, errSigma2, avaSigma3, errSigma3;
    gsl_integration_workspace * w_sigma6 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F1, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma6, &avaSigma6, &errSigma6);
    gsl_integration_workspace_free (w_sigma6);
    
    gsl_integration_workspace * w_sigma2 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F2, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma2, &avaSigma2, &errSigma2);
    gsl_integration_workspace_free (w_sigma2);
    
    gsl_integration_workspace * w_sigma3 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F3, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma3, &avaSigma3, &errSigma3);
    gsl_integration_workspace_free (w_sigma3);
    
    return avaSigma6/(2.*sqrt(-avaSigma2*avaSigma3));
}


P_6Prime::P_6Prime(const StandardModel& SM_i, StandardModel::lepton lep_i) : BKstarll(SM_i, lep_i) {  
}

double P_6Prime::computeThValue() {
    
    updateParameters();
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    
    gsl_function F1 = convertToGslFunction( boost::bind( &BKstarll::getSigma9, &(*this), _1 ) );
    gsl_function F2 = convertToGslFunction( boost::bind( &BKstarll::getSigma2, &(*this), _1 ) );
    gsl_function F3 = convertToGslFunction( boost::bind( &BKstarll::getSigma3, &(*this), _1 ) );
    
    double avaSigma9, errSigma9, avaSigma2, errSigma2, avaSigma3, errSigma3;
    gsl_integration_workspace * w_sigma9 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F1, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma9, &avaSigma9, &errSigma9);
    gsl_integration_workspace_free (w_sigma9);
    
    gsl_integration_workspace * w_sigma2 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F2, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma2, &avaSigma2, &errSigma2);
    gsl_integration_workspace_free (w_sigma2);
    
    gsl_integration_workspace * w_sigma3 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F3, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma3, &avaSigma3, &errSigma3);
    gsl_integration_workspace_free (w_sigma3);
    
    return -avaSigma9/(2.*sqrt(-avaSigma2*avaSigma3));
 
}


GammaPrime::GammaPrime(const StandardModel& SM_i, StandardModel::lepton lep_i) : BKstarll(SM_i, lep_i) {  
}

double GammaPrime::computeThValue() {
    
    updateParameters();
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    
    gsl_function F1 = convertToGslFunction( boost::bind( &BKstarll::getSigma0, &(*this), _1 ) );
    gsl_function F2 = convertToGslFunction( boost::bind( &BKstarll::getSigma1, &(*this), _1 ) );
    gsl_function F3 = convertToGslFunction( boost::bind( &BKstarll::getSigma2, &(*this), _1 ) );
    gsl_function F4 = convertToGslFunction( boost::bind( &BKstarll::getSigma3, &(*this), _1 ) );
    
    double avaSigma0, errSigma0, avaSigma1, errSigma1, avaSigma2, errSigma2, avaSigma3, errSigma3;
    gsl_integration_workspace * w_sigma0 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F1, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma0, &avaSigma0, &errSigma0);
    gsl_integration_workspace_free (w_sigma0);

    gsl_integration_workspace * w_sigma1 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F2, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma1, &avaSigma1, &errSigma1);
    gsl_integration_workspace_free (w_sigma1);
    
    gsl_integration_workspace * w_sigma2 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F3, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma2, &avaSigma2, &errSigma2);
    gsl_integration_workspace_free (w_sigma2);
    
    gsl_integration_workspace * w_sigma3 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F4, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma3, &avaSigma3, &errSigma3);
    gsl_integration_workspace_free (w_sigma3);
    
    return ((3.*avaSigma0 - avaSigma2) + 2.*(3.*avaSigma1 - avaSigma3))/4.;

}


BR_BKstarll::BR_BKstarll(const StandardModel& SM_i, StandardModel::lepton lep_i) : BKstarll(SM_i, lep_i) {  
}

double BR_BKstarll::computeThValue() {
    
    updateParameters();
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    
    gsl_function F1 = convertToGslFunction( boost::bind( &BKstarll::getSigma0, &(*this), _1 ) );
    gsl_function F2 = convertToGslFunction( boost::bind( &BKstarll::getSigma1, &(*this), _1 ) );
    gsl_function F3 = convertToGslFunction( boost::bind( &BKstarll::getSigma2, &(*this), _1 ) );
    gsl_function F4 = convertToGslFunction( boost::bind( &BKstarll::getSigma3, &(*this), _1 ) );
    
    double avaSigma0, errSigma0, avaSigma1, errSigma1, avaSigma2, errSigma2, avaSigma3, errSigma3;
    gsl_integration_workspace * w_sigma0 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F1, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma0, &avaSigma0, &errSigma0);
    gsl_integration_workspace_free (w_sigma0);

    gsl_integration_workspace * w_sigma1 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F2, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma1, &avaSigma1, &errSigma1);
    gsl_integration_workspace_free (w_sigma1);
    
    gsl_integration_workspace * w_sigma2 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F3, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma2, &avaSigma2, &errSigma2);
    gsl_integration_workspace_free (w_sigma2);
    
    gsl_integration_workspace * w_sigma3 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F4, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma3, &avaSigma3, &errSigma3);
    gsl_integration_workspace_free (w_sigma3);
    
    return ((3.*avaSigma0 - avaSigma2) + 2.*(3.*avaSigma1 - avaSigma3))/4./width_Bd;

}


ACP::ACP(const StandardModel& SM_i, StandardModel::lepton lep_i) : BKstarll(SM_i, lep_i), mySM(SM_i) {
    lep = lep_i;
}

double ACP::computeThValue() {
    
    updateParameters();
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    
    gsl_function F1 = convertToGslFunction( boost::bind( &BKstarll::getDelta0, &(*this), _1 ) );
    gsl_function F2 = convertToGslFunction( boost::bind( &BKstarll::getDelta1, &(*this), _1 ) );
    gsl_function F3 = convertToGslFunction( boost::bind( &BKstarll::getDelta2, &(*this), _1 ) );
    gsl_function F4 = convertToGslFunction( boost::bind( &BKstarll::getDelta3, &(*this), _1 ) );
    
    double avaDelta0, errDelta0, avaDelta1, errDelta1, avaDelta2, errDelta2, avaDelta3, errDelta3;
    gsl_integration_workspace * w_delta0 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F1, q_min, q_max, 1.e-5, 1.e-3, 1000, w_delta0, &avaDelta0, &errDelta0);
    gsl_integration_workspace_free (w_delta0);
    
    gsl_integration_workspace * w_delta1 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F2, q_min, q_max, 1.e-5, 1.e-3, 1000, w_delta1, &avaDelta1, &errDelta1);
    gsl_integration_workspace_free (w_delta1);
    
    gsl_integration_workspace * w_delta2 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F3, q_min, q_max, 1.e-5, 1.e-3, 1000, w_delta2, &avaDelta2, &errDelta2);
    gsl_integration_workspace_free (w_delta2);
    
    gsl_integration_workspace * w_delta3 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F4, q_min, q_max, 1.e-5, 1.e-3, 1000, w_delta3, &avaDelta3, &errDelta3);
    gsl_integration_workspace_free (w_delta3);
    
    GammaPrime myGammaPrime(mySM, lep);
    myGammaPrime.setBinMin(q_min);
    myGammaPrime.setBinMax(q_max);
            
    return (3.*avaDelta0 - avaDelta2 + 2. * ( 3*avaDelta1 - avaDelta3 ) )/(4.*myGammaPrime.computeThValue());

}


P3CP::P3CP(const StandardModel& SM_i, StandardModel::lepton lep_i) : BKstarll(SM_i, lep_i) {  
}

double P3CP::computeThValue() {
    
    updateParameters();
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    
    gsl_function F1 = convertToGslFunction( boost::bind( &BKstarll::getDelta11, &(*this), _1 ) );
    gsl_function F2 = convertToGslFunction( boost::bind( &BKstarll::getSigma3, &(*this), _1 ) );

    double avaDelta11, errDelta11, avaSigma3, errSigma3;
    gsl_integration_workspace * w_delta11 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F1, q_min, q_max, 1.e-5, 1.e-3, 1000, w_delta11, &avaDelta11, &errDelta11);   // integrational parameters tuned
    gsl_integration_workspace_free (w_delta11);

    gsl_integration_workspace * w_sigma3 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F2, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma3, &avaSigma3, &errSigma3);
    gsl_integration_workspace_free (w_sigma3);

    return - avaDelta11/(4.*avaSigma3);

}


F_L::F_L(const StandardModel& SM_i, StandardModel::lepton lep_i) : BKstarll(SM_i, lep_i) , mySM(SM_i) {
    lep = lep_i;
}

double F_L::computeThValue() {
    
    updateParameters();
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    
    gsl_function F1 = convertToGslFunction( boost::bind( &BKstarll::getSigma0, &(*this), _1 ) );
    gsl_function F2 = convertToGslFunction( boost::bind( &BKstarll::getSigma2, &(*this), _1 ) );
    
    double avaSigma0, errSigma0, avaSigma2, errSigma2;
    gsl_integration_workspace * w_sigma0 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F1, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma0, &avaSigma0, &errSigma0);
    gsl_integration_workspace_free (w_sigma0);
    
    gsl_integration_workspace * w_sigma2 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F2, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma2, &avaSigma2, &errSigma2);
    gsl_integration_workspace_free (w_sigma2);
    
    GammaPrime myGammaPrime(mySM, lep);
    myGammaPrime.setBinMin(q_min);
    myGammaPrime.setBinMax(q_max);
    
    return (3.*avaSigma0 - avaSigma2)/(4.*myGammaPrime.computeThValue());

}


M_1Prime::M_1Prime(const StandardModel& SM_i, StandardModel::lepton lep_i) : BKstarll(SM_i, lep_i) {  
}

double M_1Prime::computeThValue() {
    
    updateParameters();
    double q_min = getBinMin();
    /*double q_max = getBinMax();
    
    
    gsl_function F1 = convertToGslFunction( boost::bind( &BKstarll::getHV1_abs2, &(*this), _1 ) );
    gsl_function F2 = convertToGslFunction( boost::bind( &BKstarll::getHV2_abs2, &(*this), _1 ) );
    gsl_function F3 = convertToGslFunction( boost::bind( &BKstarll::getHA1_abs2, &(*this), _1 ) );
    gsl_function F4 = convertToGslFunction( boost::bind( &BKstarll::getHA2_abs2, &(*this), _1 ) );
    
    double avaHV1_abs2, errHV1_abs2, avaHV2_abs2, errHV2_abs2, avaHA1_abs2, errHA1_abs2, avaHA2_abs2, errHA2_abs2;
    gsl_integration_workspace * w_hv1_abs2 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F1, q_min, q_max, 0, 1e-7, 1000, w_hv1_abs2, &avaHV1_abs2, &errHV1_abs2);
    gsl_integration_workspace_free (w_hv1_abs2);
    
    gsl_integration_workspace * w_hv2_abs2 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F2, q_min, q_max, 0, 1e-7, 1000, w_hv1_abs2, &avaHV2_abs2, &errHV2_abs2);
    gsl_integration_workspace_free (w_hv2_abs2);
    
    gsl_integration_workspace * w_ha1_abs2 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F3, q_min, q_max, 0, 1e-7, 1000, w_ha1_abs2, &avaHA1_abs2, &errHA1_abs2);
    gsl_integration_workspace_free (w_ha1_abs2);
    
    gsl_integration_workspace * w_ha2_abs2 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F4, q_min, q_max, 0, 1e-7, 1000, w_ha2_abs2, &avaHA2_abs2, &errHA2_abs2);
    gsl_integration_workspace_free (w_ha2_abs2);*/
    
    return ( H_V(1,q_min,0).abs2() + H_V(2,q_min,0).abs2() - H_A(1,q_min,0).abs2() - H_A(2,q_min,0).abs2() )/
            ( 2*( H_V(1,q_min,0).abs2() + H_V(2,q_min,0).abs2() + H_A(1,q_min,0).abs2() + H_A(2,q_min,0).abs2() ) );
}


M_2Prime::M_2Prime(const StandardModel& SM_i, StandardModel::lepton lep_i) : BKstarll(SM_i, lep_i) {  
}

double M_2Prime::computeThValue() {
    
    updateParameters();
    double q_min = getBinMin();
    /*double q_max = getBinMax();
    
    
    gsl_function F1 = convertToGslFunction( boost::bind( &BKstarll::getHV0_abs2, &(*this), _1 ) );
    gsl_function F2 = convertToGslFunction( boost::bind( &BKstarll::getHA0_abs2, &(*this), _1 ) );
    gsl_function F3 = convertToGslFunction( boost::bind( &BKstarll::getHS_abs2, &(*this), _1 ) );
    gsl_function F4 = convertToGslFunction( boost::bind( &BKstarll::getHP_abs2, &(*this), _1 ) );
    gsl_function F5 = convertToGslFunction( boost::bind( &BKstarll::beta, &(*this), _1 ) );
    gsl_function F6 = convertToGslFunction( boost::bind( &BKstarll::getFactor, &(*this), _1 ) );


    
    double avaHV0_abs2, errHV0_abs2, avaHA0_abs2, errHA0_abs2, avaHS_abs2, errHS_abs2, avaHP_abs2, errHP_abs2;
    gsl_integration_workspace * w_hv0_abs2 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F1, q_min, q_max, 0, 1e-7, 1000, w_hv0_abs2, &avaHV0_abs2, &errHV0_abs2);
    gsl_integration_workspace_free (w_hv0_abs2);
    
    gsl_integration_workspace * w_ha0_abs2 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F2, q_min, q_max, 0, 1e-7, 1000, w_ha0_abs2, &avaHA0_abs2, &errHA0_abs2);
    gsl_integration_workspace_free (w_ha0_abs2);
    
    gsl_integration_workspace * w_hs_abs2 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F3, q_min, q_max, 0, 1e-7, 1000, w_hs_abs2, &avaHS_abs2, &errHS_abs2);
    gsl_integration_workspace_free (w_hs_abs2);
    
    gsl_integration_workspace * w_hp_abs2 = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F4, q_min, q_max, 0, 1e-7, 1000, w_hp_abs2, &avaHP_abs2, &errHP_abs2);
    gsl_integration_workspace_free (w_hp_abs2);
    
    double avabeta, errbeta, avafactor, errfactor;
    gsl_integration_workspace * w_beta = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F5, q_min, q_max, 0, 1e-7, 1000, w_beta, &avabeta, &errbeta);
    gsl_integration_workspace_free (w_beta);
    
    gsl_integration_workspace * w_factor = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F6, q_min, q_max, 0, 1e-7, 1000, w_factor, &avafactor, &errfactor);
    gsl_integration_workspace_free (w_factor);*/
    
    return ( q_min/(2*Mm*Mm)*( H_P(q_min,0).abs2() + beta(q_min)*beta(q_min)*H_S(q_min,0).abs2() ) + H_V(0,q_min,0).abs2() - H_A(0,q_min,0).abs2() )/
            ( H_V(0,q_min,0).abs2() + H_A(0,q_min,0).abs2() );  
}