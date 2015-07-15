/* 
 * Copyright (C) 2014 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Flavour.h"
#include "MPll.h"
#include <gslpp.h>
#include <gslpp_complex.h>
#include <gsl/gsl_math.h>
#include <boost/bind.hpp>





MPll::MPll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson pseudoscalar_i, StandardModel::lepton lep_i) 
:       mySM(SM_i),
        fplus_cache(2, 0.),
        fT_cache(2, 0.),
        k2_cache(2, 0.),
        SL_cache(2, 0.),
        N_cache(3, 0.),
        Ycache(3, 0.),
        H_V0cache(2, 0.),
        H_Scache(2, 0.),
        H_P_cache(4, 0.)
{
    lep = lep_i;
    meson = meson_i;
    pseudoscalar = pseudoscalar_i;
    I0_updated = 0;
    I2_updated = 0;
    I8_updated = 0;
    
    w_sigma0 = gsl_integration_workspace_alloc (50);
    w_sigma2 = gsl_integration_workspace_alloc (50);
    
    w_delta0 = gsl_integration_workspace_alloc (50);
    w_delta2 = gsl_integration_workspace_alloc (50);
}


MPll::~MPll() 
{}

void MPll::updateParameters()
{
    GF = mySM.getGF();
    ale=mySM.getAle();
    Mlep=mySM.getLeptons(lep).getMass();
    MM=mySM.getMesons(meson).getMass();
    MP=mySM.getMesons(pseudoscalar).getMass();
    Mb=mySM.getQuarks(QCD::BOTTOM).getMass();    // add the PS b mass
    Mc=mySM.getQuarks(QCD::CHARM).getMass();
    Ms=mySM.getQuarks(QCD::STRANGE).getMass();
    MW=mySM.Mw();
    lambda_t=mySM.computelamt_s();
    mu_b = mySM.getMub();
    width = mySM.getMesons(meson).computeWidth();
    
    switch(pseudoscalar){
        case StandardModel::K_P :
            r_1_fplus = mySM.getr_1_fplus();
            r_2_fplus = mySM.getr_2_fplus();
            m_fit2_fplus = mySM.getm_fit2_fplus();
            r_1_fT = mySM.getr_1_fT();
            r_2_fT = mySM.getr_2_fT();
            m_fit2_fT = mySM.getm_fit2_fT();
            r_2_f0 = mySM.getr_2_f0();
            m_fit2_f0 = mySM.getm_fit2_f0();
            break;
        default:
            std::stringstream out;
            out << pseudoscalar;
            throw std::runtime_error("MPll: pseudoscalar " + out.str() + " not implemented");
    }
    
    h_0 = mySM.geth_0_MP();
    h_0_1 = mySM.geth_0_1_MP();
    
    b = 1.;                           //please check
    
    allcoeff = mySM.getMyFlavour()->ComputeCoeffBMll(mu_b);   //check the mass scale, scheme fixed to NDR
    allcoeffprime = mySM.getMyFlavour()->ComputeCoeffprimeBMll(mu_b);   //check the mass scale, scheme fixed to NDR

    C_1 = (*(allcoeff[LO]))(0) + (*(allcoeff[NLO]))(0);
    C_2 =  ((*(allcoeff[LO]))(1) + (*(allcoeff[NLO]))(1));
    C_3 = ((*(allcoeff[LO]))(2) + (*(allcoeff[NLO]))(2));
    C_4 = (*(allcoeff[LO]))(3) + (*(allcoeff[NLO]))(3);
    C_5 =  ((*(allcoeff[LO]))(4) + (*(allcoeff[NLO]))(4));
    C_6 = ((*(allcoeff[LO]))(5) + (*(allcoeff[NLO]))(5));
    C_7 = (*(allcoeff[LO]))(6) + (*(allcoeff[NLO]))(6);
    C_9 = (*(allcoeff[LO]))(8) + (*(allcoeff[NLO]))(8);
    C_10 = (*(allcoeff[LO]))(9) + (*(allcoeff[NLO]))(9);
    C_S = (*(allcoeff[LO]))(10) + (*(allcoeff[NLO]))(10);
    C_P = (*(allcoeff[LO]))(11) + (*(allcoeff[NLO]))(11);
    
    C_7p = (*(allcoeffprime[LO]))(6) + (*(allcoeffprime[NLO]))(6);
    C_9p = (*(allcoeffprime[LO]))(8) + (*(allcoeffprime[NLO]))(8);
    C_10p = (*(allcoeffprime[LO]))(9) + (*(allcoeffprime[NLO]))(9);
    C_Sp = (*(allcoeffprime[LO]))(10) + (*(allcoeffprime[NLO]))(10);
    C_Pp = (*(allcoeffprime[LO]))(11) + (*(allcoeffprime[NLO]))(11);
    
    checkCache();
    
    std::map<std::pair<double, double>, unsigned int >::iterator it;
    
    if (I0_updated == 0) for (it = sigma0Cached.begin(); it != sigma0Cached.end(); ++it) it->second = 0;
    if (I2_updated == 0) for (it = sigma2Cached.begin(); it != sigma2Cached.end(); ++it) it->second = 0;
    
    if (I0_updated == 0) for (it = delta0Cached.begin(); it != delta0Cached.end(); ++it) it->second = 0;
    if (I2_updated == 0) for (it = delta2Cached.begin(); it != delta2Cached.end(); ++it) it->second = 0;
    
}

void MPll::checkCache()
{
    
    if(r_1_fplus == fplus_cache(0) && r_2_fplus == fplus_cache(1)) {
        fplus_updated = 1;
    } else {
        fplus_updated = 0;
        fplus_cache(0) = r_1_fplus;
        fplus_cache(1) = r_2_fplus;
    }
    
    if(r_1_fT == fT_cache(0) && r_2_fT == fT_cache(1)) {
        fT_updated = 1;
    } else {
        fT_updated = 0;
        fT_cache(0) = r_1_fT;
        fT_cache(1) = r_2_fT;
    }
    
    if(r_2_f0 == f0_cache) {
        f0_updated = 1;
    } else {
        f0_updated = 0;
        f0_cache = r_2_f0;
    }
    
    if (MM == k2_cache(0) && MP == k2_cache(1) ) {
        k2_updated = 1;
    } else {
        k2_updated = 0;
        k2_cache(0) = MM;
        k2_cache(1) = MP;
    }
    
    if (Mlep == beta_cache) {
        beta_updated = 1;
    } else {
        beta_updated = 0;
        beta_cache = Mlep;
    }
    
    if (MP == lambda_cache) {
        lambda_updated = k2_updated;
        F_updated = lambda_updated * beta_updated;
    } else {
        lambda_updated = 0;
        F_updated = 0;
        lambda_cache = MP;
    }
    
    VL_updated = k2_updated * fplus_updated;
    
    VR_updated = VL_updated;
    
    TL_updated = k2_updated * fT_updated;
    
    TR_updated = TL_updated;
    
    if (Mb == SL_cache(0) && Ms == SL_cache(1) ){
        SL_updated = k2_updated * f0_updated;
        SR_updated = SL_updated;
    } else {
        SL_updated = 0;
        SR_updated = SL_updated;
        SL_cache(0) = Mb;
        SL_cache(1) = Ms;
    }
    
    if (GF == N_cache(0) && ale == N_cache(1) && MM == N_cache(2) && lambda_t == Nc_cache) {
        N_updated = 1;
    } else {
        N_updated = 0;
        N_cache(0) = GF;
        N_cache(1) = ale;
        N_cache(2) = MM;
        Nc_cache = lambda_t;
    }
    
    if (C_1 == C_1_cache) {
        C_1_updated = 1;
    } else {
        C_1_updated = 0;
        C_1_cache = C_1;
    }
    
    if (C_2 == C_2_cache) {
        C_2_updated = 1;
    } else {
        C_2_updated = 0;
        C_2_cache = C_2;
    }
    
    if (C_3 == C_3_cache) {
        C_3_updated = 1;
    } else {
        C_3_updated = 0;
        C_3_cache = C_3;
    }
    
    if (C_4 == C_4_cache) {
        C_4_updated = 1;
    } else {
        C_4_updated = 0;
        C_4_cache = C_4;
    }
    
    if (C_5 == C_5_cache) {
        C_5_updated = 1;
    } else {
        C_5_updated = 0;
        C_5_cache = C_5;
    }
    
    if (C_6 == C_6_cache) {
        C_6_updated = 1;
    } else {
        C_6_updated = 0;
        C_6_cache = C_6;
    }
    
    if (C_7 == C_7_cache) {
        C_7_updated = 1;
    } else {
        C_7_updated = 0;
        C_7_cache = C_7;
    }
    
    if (C_9 == C_9_cache) {
        C_9_updated = 1;
    } else {
        C_9_updated = 0;
        C_9_cache = C_9;
    }
    
    if (C_10 == C_10_cache) {
        C_10_updated = 1;
    } else {
        C_10_updated = 0;
        C_10_cache = C_10;
    }
    
    if (C_S == C_S_cache) {
        C_S_updated = 1;
    } else {
        C_S_updated = 0;
        C_S_cache = C_S;
    }
    
    if (C_P == C_P_cache) {
        C_P_updated = 1;
    } else {
        C_P_updated = 0;
        C_P_cache = C_P;
    }
    
    if (C_7p == C_7p_cache) {
        C_7p_updated = 1;
    } else {
        C_7p_updated = 0;
        C_7p_cache = C_7p;
    }
    
    if (C_9p == C_9p_cache) {
        C_9p_updated = 1;
    } else {
        C_9p_updated = 0;
        C_9p_cache = C_9p;
    }
    
    if (C_10p == C_10p_cache) {
        C_10p_updated = 1;
    } else {
        C_10p_updated = 0;
        C_10p_cache = C_10p;
    }
    
    if (C_Sp == C_Sp_cache) {
        C_Sp_updated = 1;
    } else {
        C_Sp_updated = 0;
        C_Sp_cache = C_Sp;
    }
    
    if (C_Pp == C_Pp_cache) {
        C_Pp_updated = 1;
    } else {
        C_Pp_updated = 0;
        C_Pp_cache = C_Pp;
    }
    
    if (Mb == Ycache(0) && Mc == Ycache(1) ) {
        Yupdated = C_1_updated * C_2_updated * C_3_updated * C_4_updated * C_5_updated * C_6_updated;
    } else {
        Yupdated = 0;
        Ycache(0) = Mb;
        Ycache(1) = Mc;
    }
    
    if (MM == H_V0cache(0) && Mb == H_V0cache(1) && h_0 == H_V0Ccache[0] && h_0_1 == H_V0Ccache[1]) {
        H_V0updated = N_updated * C_9_updated * Yupdated * VL_updated * C_9p_updated * VR_updated * C_7_updated * TL_updated * C_7p_updated * TR_updated;
    } else {
        H_V0updated = 0;
        H_V0cache(0) = MM;
        H_V0cache(1) = Mb;
        H_V0Ccache[0] = h_0;
        H_V0Ccache[1] = h_0_1;
    }
    
    H_A0updated = N_updated * C_10_updated * VL_updated * C_10p_updated * VR_updated;
    
    if (Mb == H_Scache(0) && MW == H_Scache(1)) {
        H_Supdated = N_updated * C_S_updated * SL_updated * C_Sp_updated * SR_updated;
    } else {
        H_Supdated = 0;
        H_Scache(0) = Mb;
        H_Scache(1) = MW;
    }
    
    if (Mb == H_P_cache(0) && MW == H_P_cache(1) && Mlep == H_P_cache(2) && Ms == H_P_cache(3)) {
        H_P_updated = N_updated * C_P_updated * SL_updated * C_Pp_updated * SR_updated * C_10_updated * C_10p_updated;
    } else {
        H_P_updated = 0;
        H_P_cache(0) = Mb;
        H_P_cache(1) = MW;
        H_P_cache(2) = Mlep;
        H_P_cache(3) = Ms;   
    }
    
    I0_updated = F_updated * H_V0updated * H_A0updated * H_P_updated * beta_updated * H_Supdated;
    I2_updated = F_updated * beta_updated * H_V0updated * H_A0updated;
    I8_updated = F_updated * beta_updated * H_Supdated * H_V0updated;
    
}

/*******************************************************************************
 * Transverse Form Factors                                                     *
 * ****************************************************************************/
double MPll::LCSR_fit1(double q2, double r_1, double r_2, double m_fit2)
{
    return r_1/( 1. - q2/m_fit2 ) + r_2/pow( ( 1. - q2/m_fit2 ) , 2.) ;

}

double MPll::LCSR_fit2(double q2, double r_2, double m_fit2)
{
    return r_2/( 1. - q2/m_fit2 ) ; 
}

double MPll::f_plus(double q2)
{
    return LCSR_fit1(q2, r_1_fplus, r_2_fplus, m_fit2_fplus);
}

double MPll::f_T(double q2)
{
    return LCSR_fit1(q2, r_1_fT, r_2_fT, m_fit2_fT);
}

double MPll::f_0(double q2)
{
    return LCSR_fit2(q2, r_2_f0, m_fit2_f0);
}

gslpp::complex MPll::V_L(double q2)
{
    return gslpp::complex::i() * sqrt(lambda(q2)) / (2*MM*sqrt(q2)) * f_plus(q2);
}

gslpp::complex MPll::V_R(double q2)
{
    return -V_L(q2);
}

gslpp::complex MPll::T_L(double q2)
{
    return gslpp::complex::i()  * sqrt(lambda(q2)*q2) / (2.*MM*MM*(MM+MP)) * f_T(q2);
}

gslpp::complex MPll::T_R(double q2)
{
    return -T_L(q2);
}

double MPll::S_L(double q2)
{
    return -( MM*MM - MP*MP )/(2*MM*(Mb + Ms)) * ( 1 + Ms/Mb )/( 1 - Ms/Mb ) * f_0(q2);
}

double MPll::S_R(double q2)
{
    return -S_L(q2);
}



/*******************************************************************************
 * Helicity amplitudes                                                         *
 * ****************************************************************************/
gslpp::complex MPll::N()
{
    return -(4.*GF*MM*ale*lambda_t)/(sqrt(2.)*4.*M_PI);
}

gslpp::complex MPll::H(double q2, double m)
{
    double x = 4.*m*m/q2;
    gslpp::complex par;
    
    if (x>1.) par = sqrt(x - 1.) * atan( 1. / sqrt(x - 1.) );
    else par = sqrt(1. - x) * ( log( ( 1. + sqrt(1. - x) ) / sqrt(x) ) - gslpp::complex::i()*M_PI/2.);
    
    if (x == 0.) return 8. / 27. * (1. + 3. * gslpp::complex::i() * M_PI);
    else return - 4./9. * ( log( m*m / q2 ) - 2./3. - x ) - 4./9. * (-2. + x) * par;
}

gslpp::complex MPll::Y(double q2)
{
    return 4./3. * C_3 + 64./9. * C_5 + 64./27. * C_6 - 1./2. * H(q2,0.) * ( C_3 + 4./3.*C_4 + 16. * C_5 + 64./3.*C_6 )
            + H(q2, Mc) * ( 4./3.*C_1 + C_2 + 6.*C_3 + 60.*C_5 ) - 1./2. * H(q2, Mb) * ( 7.*C_3 + 4./3.*C_4 + 76.*C_5 + 64./3.*C_6 );
}

gslpp::complex MPll::H_V(double q2, int bar) 
{
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
                    
    return -gslpp::complex::i()*n*( ( C_9 + Y(q2) )*V_L(q2) 
            + C_9p*V_R(q2) + MM*MM/q2*( 2*Mb/MM*( C_7*T_L(q2) + C_7p*T_R(q2) ) - 16*M_PI*M_PI*(h_0 + h_0_1 * q2)) );
}

gslpp::complex MPll::H_A(double q2, int bar) 
{
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

gslpp::complex MPll::H_S(double q2, int bar) 
{
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

gslpp::complex MPll::H_P(double q2, int bar) 
{
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
double MPll::k2(double q2) 
{
    return (pow(MM,4.) + q2*q2 + pow(MP,4.) -2.*MP*MP*q2 -2.*MM*MM*(q2 + MP*MP))/(4.*MM*MM);
}

double MPll::beta(double q2) 
{
    return sqrt(1. - 4.*Mlep*Mlep/q2);
}

double MPll::lambda(double q2) 
{
    return 4.*MM*MM*k2(q2);
}

double MPll::F(double q2) 
{
    return sqrt(lambda(q2))*beta(q2)*q2/(96.*M_PI*M_PI*M_PI*MM*MM*MM);
}

double MPll::I(int i, double q2, int bar) 
{

    double Mlep2 = Mlep*Mlep;
    double beta2 = beta(q2)*beta(q2);
    

    switch (i){
        case 0: // I1c
            return F(q2)*( ( H_V(q2,bar).abs2() + H_A(q2,bar).abs2() )/2.  +  H_P(q2,bar).abs2()  +  2.*Mlep2/q2*( H_V(q2,bar).abs2() 
                    - H_A(q2,bar).abs2() )  + beta2*H_S(q2,bar).abs2() );
        case 2: // I2c
            return -F(q2)*beta2/2.*( H_V(q2,bar).abs2() + H_A(q2,bar).abs2() );
        case 8: // I6c
            return 2.*F(q2)*beta(q2)*Mlep/sqrt(q2)*( H_S(q2,bar).conjugate()*H_V(q2,bar) ).real();
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("I: index " + out.str() + " not implemented");
    }
}

double MPll::Sigma(int i, double q2) 
{
    return (I(i, q2, 0) + I(i, q2, 1))/2.;
}

double MPll::Delta(int i, double q2) 
{
    return (I(i, q2, 0) - I(i, q2, 1))/2.;
}

double MPll::integrateSigma(int i, double q_min, double q_max)
{
    
    if (mySM.getMyFlavour()->getUpdateFlag(meson, pseudoscalar, lep)){
        updateParameters();
        mySM.getMyFlavour()->setUpdateFlag(meson, pseudoscalar, lep, false);
    }
    
    std::pair<double, double > qbin = std::make_pair(q_min, q_max);
    
    switch(i){
        case 0:
            if (sigma0Cached[qbin] == 0) {
                FS0 = convertToGslFunction( boost::bind( &MPll::getSigma0, &(*this), _1 ) );
                gsl_integration_qags (&FS0, q_min, q_max, 1.e-5, 1.e-3, 50, w_sigma0, &avaSigma0, &errSigma0);
                cacheSigma0[qbin] = avaSigma0;
                sigma0Cached[qbin] = 1;
            }
            return cacheSigma0[qbin];
            break;
        case 2:
            if (sigma2Cached[qbin] == 0) {
                FS2 = convertToGslFunction( boost::bind( &MPll::getSigma2, &(*this), _1 ) );
                gsl_integration_qags (&FS2, q_min, q_max, 1.e-5, 1.e-3, 50, w_sigma2, &avaSigma2, &errSigma2);
                cacheSigma2[qbin] = avaSigma2;
                sigma2Cached[qbin] = 1;
            }
            return cacheSigma2[qbin];
            break;
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("MVll::integrateSigma_e: index " + out.str() + " not implemented");
    }
    
}

double MPll::integrateDelta(int i, double q_min, double q_max)
{
    
    if (mySM.getMyFlavour()->getUpdateFlag(meson, pseudoscalar, lep)){
        updateParameters();
        mySM.getMyFlavour()->setUpdateFlag(meson, pseudoscalar, lep, false);
    }
        
    std::pair<double, double > qbin = std::make_pair(q_min, q_max);
    switch(i){
        case 0:
            if (delta0Cached[qbin] == 0) {
                FD0 = convertToGslFunction( boost::bind( &MPll::getDelta0, &(*this), _1 ) );
                gsl_integration_qags (&FD0, q_min, q_max, 1.e-5, 1.e-3, 50, w_delta0, &avaDelta0, &errDelta0);
                cacheDelta0[qbin] = avaDelta0;
                delta0Cached[qbin] = 1;
            }
            return cacheDelta0[qbin];
            break;
        case 2:
            if (delta2Cached[qbin] == 0) {
                FD2 = convertToGslFunction( boost::bind( &MPll::getDelta2, &(*this), _1 ) );
                gsl_integration_qags (&FD2, q_min, q_max, 1.e-5, 1.e-3, 50, w_delta2, &avaDelta2, &errDelta2);
                cacheDelta2[qbin] = avaDelta2;
                delta2Cached[qbin] = 1;
            }
            return cacheDelta2[qbin];
            break;
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("integrateDelta_e: index " + out.str() + " not implemented"); 
    }
}