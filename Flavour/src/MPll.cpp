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





MPll::MPll(const StandardModel& SM_i) : ThObservable(SM_i), mySM(SM_i),
        N_cache(3, 0.),
        fplus_cache(2, 0.),
        fT_cache(2, 0.),
        k2_cache(2, 0.),
        beta_cache(2, 0.),
        SL_cache(2, 0.),
        H_V0cache(2, 0.),
        H_Scache(2, 0.),
        H_Pe_cache(4, 0.),
        H_Pmu_cache(4, 0.)
{
    I0e_updated = 0;
    I2e_updated = 0;
    I8e_updated = 0;
    I0mu_updated = 0;
    I2mu_updated = 0;
    I8mu_updated = 0;
    iter = 0;
    
    w_sigma0e = gsl_integration_workspace_alloc (50);
    w_sigma2e = gsl_integration_workspace_alloc (50);
    
    w_delta0e = gsl_integration_workspace_alloc (50);
    w_delta2e = gsl_integration_workspace_alloc (50);
    
    w_sigma0mu = gsl_integration_workspace_alloc (50);
    w_sigma2mu = gsl_integration_workspace_alloc (50);
    
    w_delta0mu = gsl_integration_workspace_alloc (50);
    w_delta2mu = gsl_integration_workspace_alloc (50);
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
    
    h_0 = mySM.geth_0();
    h_0_1 = mySM.geth_0_1();
    
    b = 1.;                           //please check
    
    allcoeff = mySM.getMyFlavour()->ComputeCoeffBKstarll(mu_b);   //check the mass scale, scheme fixed to NDR
    allcoeffprime = mySM.getMyFlavour()->ComputeCoeffprimeBKstarll(mu_b);   //check the mass scale, scheme fixed to NDR

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
    
    if (I0e_updated == 0) for (it = sigma0eCached.begin(); it != sigma0eCached.end(); ++it) it->second = 0;
    if (I2e_updated == 0) for (it = sigma2eCached.begin(); it != sigma2eCached.end(); ++it) it->second = 0;
    
    if (I0e_updated == 0) for (it = delta0eCached.begin(); it != delta0eCached.end(); ++it) it->second = 0;
    if (I2e_updated == 0) for (it = delta2eCached.begin(); it != delta2eCached.end(); ++it) it->second = 0;
    
    if (I0mu_updated == 0) for (it = sigma0muCached.begin(); it != sigma0muCached.end(); ++it) it->second = 0;
    if (I2mu_updated == 0) for (it = sigma2muCached.begin(); it != sigma2muCached.end(); ++it) it->second = 0;
    
    if (I0mu_updated == 0) for (it = delta0muCached.begin(); it != delta0muCached.end(); ++it) it->second = 0;
    if (I2mu_updated == 0) for (it = delta2muCached.begin(); it != delta2muCached.end(); ++it) it->second = 0;
    
}

void MPll::checkCache(){
    
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
    
    if (MB == k2_cache(0) && MK == k2_cache(1) ) {
        k2_updated = 1;
    } else {
        k2_updated = 0;
        k2_cache(0) = MB;
        k2_cache(1) = MK;
    }
    
    if (Mmu == beta_cache(0)) {
        beta_mu_updated = 1;
    } else {
        beta_mu_updated = 0;
        beta_cache(0) = Mmu;
    }
    
    if (Me == beta_cache(1)) {
        beta_e_updated = 1;
    } else {
        beta_e_updated = 0;
        beta_cache(1) = Me;
    }
    
    if (MK == lambda_cache) {
        lambda_updated = k2_updated;
        F_e_updated = lambda_updated * beta_e_updated;
        F_mu_updated = lambda_updated * beta_mu_updated;
    } else {
        lambda_updated = 0;
        F_e_updated = 0;
        F_mu_updated = 0;
        lambda_cache = MK;
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
    
    if (GF == N_cache(0) && ale == N_cache(1) && MB == N_cache(2) && lambda_t == Nc_cache) {
        N_updated = 1;
    } else {
        N_updated = 0;
        N_cache(0) = GF;
        N_cache(1) = ale;
        N_cache(2) = MB;
        Nc_cache = lambda_t;
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
    
    if (MB == H_V0cache(0) && Mb == H_V0cache(1) && h_0 == H_V0Ccache[0] && h_0_1 == H_V0Ccache[1]) {
        H_V0updated = N_updated * C_9_updated * VL_updated * C_9p_updated * VR_updated * C_7_updated * TL_updated * C_7p_updated * TR_updated;
    } else {
        H_V0updated = 0;
        H_V0cache(0) = MB;
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
    
    if (Mb == H_Pe_cache(0) && MW == H_Pe_cache(1) && Me == H_Pe_cache(2) && Ms == H_Pe_cache(3)) {
        H_Pe_updated = N_updated * C_P_updated * SL_updated * C_Pp_updated * SR_updated * C_10_updated * C_10p_updated;
    } else {
        H_Pe_updated = 0;
        H_Pe_cache(0) = Mb;
        H_Pe_cache(1) = MW;
        H_Pe_cache(2) = Me;
        H_Pe_cache(3) = Ms;   
    }
    
    if (Mb == H_Pmu_cache(0) && MW == H_Pmu_cache(1) && Mmu == H_Pmu_cache(2) && Ms == H_Pmu_cache(3)) {
        H_Pmu_updated = N_updated * C_P_updated * SL_updated * C_Pp_updated * SR_updated * C_10_updated * C_10p_updated;
    } else {
        H_Pmu_updated = 0;
        H_Pmu_cache(0) = Mb;
        H_Pmu_cache(1) = MW;
        H_Pmu_cache(2) = Mmu;
        H_Pmu_cache(3) = Ms;   
    }
    
    I0e_updated = F_e_updated * H_V0updated * H_A0updated * H_Pe_updated * beta_e_updated * H_Supdated;
    I2e_updated = F_e_updated * beta_e_updated * H_V0updated * H_A0updated;
    I8e_updated = F_e_updated * beta_e_updated * H_Supdated * H_V0updated;
    
    I0mu_updated = F_mu_updated * H_V0updated * H_A0updated * H_Pmu_updated * beta_mu_updated * H_Supdated;
    I2mu_updated = F_mu_updated * beta_mu_updated * H_V0updated * H_A0updated;
    I8mu_updated = F_mu_updated * beta_mu_updated * H_Supdated * H_V0updated;
    
    iter += 1 ;
    
}

/*******************************************************************************
 * Transverse Form Factors                                                     *
 * ****************************************************************************/
double MPll::f_plus(double q2){
    return mySM.getMyFlavour()->getMVll(QCD::B_D, QCD::K_star, StandardModel::MU)->LCSR_fit2(q2, r_1_fplus, r_2_fplus, m_fit2_fplus);
}

double MPll::f_T(double q2){
    return mySM.getMyFlavour()->getMVll(QCD::B_D, QCD::K_star, StandardModel::MU)->LCSR_fit2(q2, r_1_fT, r_2_fT, m_fit2_fT);
}

double MPll::f_0(double q2){
    return mySM.getMyFlavour()->getMVll(QCD::B_D, QCD::K_star, StandardModel::MU)->LCSR_fit3(q2, r_2_f0, m_fit2_f0);
}

gslpp::complex MPll::V_L(double q2){
    return gslpp::complex::i() * sqrt(lambda(q2)) / (2*MB*sqrt(q2)) * f_plus(q2);
}

gslpp::complex MPll::V_R(double q2){
    return -V_L(q2);
}

gslpp::complex MPll::T_L(double q2){
    return gslpp::complex::i()  * sqrt(lambda(q2)*q2) / (MB*MB*(MB+MK)) * f_T(q2);
}

gslpp::complex MPll::T_R(double q2){
    return -T_L(q2);
}

double MPll::S_L(double q2){
    return -( MB*MB - MK*MK )/(2*MB*(Mb + Ms)) * ( 1 + Ms/Mb )/( 1 - Ms/Mb ) * f_0(q2);
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
            + MB*MB/q2*( 2*Mb/MB*( C_7*T_L(q2) + C_7p*T_R(q2) ) - 16*M_PI*M_PI*(h_0 + h_0_1 * q2)) );
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
    return sqrt(1. - 4.*Mlep*Mlep/q2);
}

double MPll::lambda(double q2) {
    return 4.*MB*MB*k2(q2);
}

double MPll::F(double q2, double Mlep) {
    return sqrt(lambda(q2))*beta(q2,Mlep)*q2/(96.*M_PI*M_PI*M_PI*MB*MB*MB);
}

double MPll::I(int i, double q2, int bar, double Mlep) {

    double Mlep2 = Mlep*Mlep;
    double beta2 = beta(q2,Mlep)*beta(q2,Mlep);
    

    switch (i){
        case 0: // I1c
            return F(q2,Mlep)*( ( H_V(q2,bar).abs2() + H_A(q2,bar).abs2() )/2.  +  H_P(q2,bar,Mlep).abs2()  +  2.*Mlep2/q2*( H_V(q2,bar).abs2() 
                    - H_A(q2,bar).abs2() )  + beta2*H_S(q2,bar).abs2() );
        case 2: // I2c
            return -F(q2,Mlep)*beta2/2.*( H_V(q2,bar).abs2() + H_A(q2,bar).abs2() );
        case 8: // I6c
            return 2.*F(q2,Mlep)*beta(q2,Mlep)*Mlep/sqrt(q2)*( H_S(q2,bar).conjugate()*H_V(q2,bar) ).real();
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("I: index " + out.str() + " not implemented");
    }
}

double MPll::Sigma(int i, double q2, double Mlep) {
    return (I(i, q2, 0, Mlep) + I(i, q2, 1, Mlep))/2.;
}

double MPll::Delta(int i, double q2, double Mlep) {
    return (I(i, q2, 0, Mlep) - I(i, q2, 1, Mlep))/2.;
}

double MPll::integrateSigma_e(int i, double q_min, double q_max){
    
    if (mySM.getMyFlavour()->getUpdateFlag()){
        updateParameters();
        mySM.getMyFlavour()->setUpdateFlag(false);
    }
    
    std::pair<double, double > qbin = std::make_pair(q_min, q_max);
    
    switch(i){
        case 0:
            if (sigma0eCached[qbin] == 0) {
                FS0e = convertToGslFunction( boost::bind( &MPll::getSigma0e, &(*this), _1 ) );
                gsl_integration_qags (&FS0e, q_min, q_max, 1.e-5, 1.e-3, 50, w_sigma0e, &avaSigma0e, &errSigma0e);
                cacheSigma0e[qbin] = avaSigma0e;
                sigma0eCached[qbin] = 1;
            }
            return cacheSigma0e[qbin];
            break;
        case 2:
            if (sigma2eCached[qbin] == 0) {
                FS2e = convertToGslFunction( boost::bind( &MPll::getSigma2e, &(*this), _1 ) );
                gsl_integration_qags (&FS2e, q_min, q_max, 1.e-5, 1.e-3, 50, w_sigma2e, &avaSigma2e, &errSigma2e);
                cacheSigma2e[qbin] = avaSigma2e;
                sigma2eCached[qbin] = 1;
            }
            return cacheSigma2e[qbin];
            break;
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("MVll::integrateSigma_e: index " + out.str() + " not implemented");
    }
    
}

double MPll::integrateDelta_e(int i, double q_min, double q_max){
    
    if (mySM.getMyFlavour()->getUpdateFlag()){
        updateParameters();
        mySM.getMyFlavour()->setUpdateFlag(false);
    }
        
    std::pair<double, double > qbin = std::make_pair(q_min, q_max);
    switch(i){
        case 0:
            if (delta0eCached[qbin] == 0) {
                FD0e = convertToGslFunction( boost::bind( &MPll::getDelta0e, &(*this), _1 ) );
                gsl_integration_qags (&FD0e, q_min, q_max, 1.e-5, 1.e-3, 50, w_delta0e, &avaDelta0e, &errDelta0e);
                cacheDelta0e[qbin] = avaDelta0e;
                delta0eCached[qbin] = 1;
            }
            return cacheDelta0e[qbin];
            break;
        case 2:
            if (delta2eCached[qbin] == 0) {
                FD2e = convertToGslFunction( boost::bind( &MPll::getDelta2e, &(*this), _1 ) );
                gsl_integration_qags (&FD2e, q_min, q_max, 1.e-5, 1.e-3, 50, w_delta2e, &avaDelta2e, &errDelta2e);
                cacheDelta2e[qbin] = avaDelta2e;
                delta2eCached[qbin] = 1;
            }
            return cacheDelta2e[qbin];
            break;
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("integrateDelta_e: index " + out.str() + " not implemented"); 
    }
}

double MPll::integrateSigma_mu(int i, double q_min, double q_max){
    
    if (mySM.getMyFlavour()->getUpdateFlag()){
        updateParameters();
        mySM.getMyFlavour()->setUpdateFlag(false);
    }
    
    std::pair<double, double > qbin = std::make_pair(q_min, q_max);
    
    switch(i){
        case 0:
            if (sigma0muCached[qbin] == 0) {
                FS0mu = convertToGslFunction( boost::bind( &MPll::getSigma0mu, &(*this), _1 ) );
                gsl_integration_qags (&FS0mu, q_min, q_max, 1.e-5, 1.e-3, 50, w_sigma0mu, &avaSigma0mu, &errSigma0mu);
                cacheSigma0mu[qbin] = avaSigma0mu;
                sigma0muCached[qbin] = 1;
            }
            return cacheSigma0mu[qbin];
            break;
        case 2:
            if (sigma2muCached[qbin] == 0) {
                FS2mu = convertToGslFunction( boost::bind( &MPll::getSigma2mu, &(*this), _1 ) );
                gsl_integration_qags (&FS2mu, q_min, q_max, 1.e-5, 1.e-3, 50, w_sigma2mu, &avaSigma2mu, &errSigma2mu);
                cacheSigma2mu[qbin] = avaSigma2mu;
                sigma2muCached[qbin] = 1;
            }
            return cacheSigma2mu[qbin];
            break;
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("MVll::integrateSigma_mu: index " + out.str() + " not implemented");
    }
    
}

double MPll::integrateDelta_mu(int i, double q_min, double q_max){
    
    if (mySM.getMyFlavour()->getUpdateFlag()){
        updateParameters();
        mySM.getMyFlavour()->setUpdateFlag(false);
    }
        
    std::pair<double, double > qbin = std::make_pair(q_min, q_max);
    switch(i){
        case 0:
            if (delta0muCached[qbin] == 0) {
                FD0mu = convertToGslFunction( boost::bind( &MPll::getDelta0mu, &(*this), _1 ) );
                gsl_integration_qags (&FD0mu, q_min, q_max, 1.e-5, 1.e-3, 50, w_delta0mu, &avaDelta0mu, &errDelta0mu);
                cacheDelta0mu[qbin] = avaDelta0mu;
                delta0muCached[qbin] = 1;
            }
            return cacheDelta0mu[qbin];
            break;
        case 2:
            if (delta2muCached[qbin] == 0) {
                FD2mu = convertToGslFunction( boost::bind( &MPll::getDelta2mu, &(*this), _1 ) );
                gsl_integration_qags (&FD2mu, q_min, q_max, 1.e-5, 1.e-3, 50, w_delta2mu, &avaDelta2mu, &errDelta2mu);
                cacheDelta2mu[qbin] = avaDelta2mu;
                delta2muCached[qbin] = 1;
            }
            return cacheDelta2mu[qbin];
            break;
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("integrateDelta_mu: index " + out.str() + " not implemented"); 
    }
}
/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/


BR_MPll_e::BR_MPll_e(const StandardModel& SM_i) : MPll(SM_i) {  
}

double BR_MPll_e::computeBR_MPll_e(double qmin, double qmax) {
    
    double q_min = qmin;
    double q_max = qmax;
    
    return (3.*integrateSigma_e(0,q_min,q_max) - integrateSigma_e(2,q_min,q_max))/(4. * width_Bd);
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
    
    return (3.*integrateSigma_mu(0,q_min,q_max) - integrateSigma_mu(2,q_min,q_max))/(4. * width_Bd);
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
    
    return (3.*integrateDelta_e(0, q_min, q_max) - integrateDelta_e(2, q_min, q_max))/(4.*computeBR_MPll_e(q_min, q_max)* width_Bd);
}

ACP_mu::ACP_mu(const StandardModel& SM_i) : BR_MPll_mu(SM_i){
}

double ACP_mu::computeThValue() {
    
    updateParameters();
    
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return (3.*integrateDelta_mu(0, q_min, q_max) - integrateDelta_mu(2, q_min, q_max))/(4.*computeBR_MPll_e(q_min, q_max)* width_Bd);
}