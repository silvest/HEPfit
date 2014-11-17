/* 
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Flavour.h"
#include "MVll.h"
#include <gslpp_complex.h>
#include <boost/bind.hpp>


MVll::MVll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) : mySM(SM_i),
        H_V0cache(2, 0.),
        H_V1cache(2, 0.),
        H_V2cache(2, 0.),
        H_Scache(2, 0.),
        H_Pcache(4, 0.),
        k2_cache(2, 0.),
        N_cache(3, 0.),
        A0_cache(4, 0.),
        A1_cache(3, 0.),
        A2_cache(2, 0.),
        V_cache(4, 0.),
        T1_cache(4, 0.),
        T2_cache(3, 0.),
        T3t_cache(2, 0.),
        SL_cache(2, 0.),
        VL0_cache(3, 0.),
        TL0_cache(2, 0.)
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    I0_updated = 0;
    I1_updated = 0;
    I2_updated = 0;
    I3_updated = 0;
    I4_updated = 0;
    I5_updated = 0;
    I6_updated = 0;
    I7_updated = 0;
    I8_updated = 0;
    I9_updated = 0;
    I10_updated = 0;
    I11_updated = 0;
    it = 0;
    
    w_sigma0 = gsl_integration_workspace_alloc (50);
    w_sigma1 = gsl_integration_workspace_alloc (50);
    w_sigma2 = gsl_integration_workspace_alloc (50);
    w_sigma3 = gsl_integration_workspace_alloc (50);
    w_sigma4 = gsl_integration_workspace_alloc (50);
    w_sigma5 = gsl_integration_workspace_alloc (50);
    w_sigma6 = gsl_integration_workspace_alloc (50);
    w_sigma7 = gsl_integration_workspace_alloc (50);
    w_sigma9 = gsl_integration_workspace_alloc (50);
    w_sigma11 = gsl_integration_workspace_alloc (50);
    
    w_delta0 = gsl_integration_workspace_alloc (50);
    w_delta1 = gsl_integration_workspace_alloc (50);
    w_delta2 = gsl_integration_workspace_alloc (50);
    w_delta3 = gsl_integration_workspace_alloc (50);
    w_delta11 = gsl_integration_workspace_alloc (50);
}


MVll::~MVll() {
/** Check to see if GSL pointers are released!!*/
}

void MVll::updateParameters(){
    GF = mySM.getGF();
    ale=mySM.getAle();
    Mlep=mySM.getLeptons(lep).getMass();
    MB=mySM.getMesons(meson).getMass();
    MKstar=mySM.getMesons(vectorM).getMass();
    Mb=mySM.getQuarks(QCD::BOTTOM).getMass();    // add the PS b mass
    Ms=mySM.getQuarks(QCD::STRANGE).getMass();
    MW=mySM.Mw();
    lambda_t=mySM.computelamt_s();
    mu_b = mySM.getMub();
    width_Bd = mySM.getMesons(meson).computeWidth();
    
    switch(vectorM){
        case QCD::K_star : {
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
            
            b=1;
        }
        case QCD::PHI : {
            a_0V=mySM.geta_0Vphi();
            a_1V=mySM.geta_1Vphi();
            dmV=mySM.getdmVphi();
            a_0A0=mySM.geta_0A0phi();
            a_1A0=mySM.geta_1A0phi();
            dmA0=mySM.getdmA0phi();
            a_0A1=mySM.geta_0A1phi();
            a_1A1=mySM.geta_1A1phi();
            dmA1=mySM.getdmA1phi();
            a_0A12=mySM.geta_0A12phi();
            a_1A12=mySM.geta_1A12phi();
            dmA12=mySM.getdmA12phi();
            a_0T1=mySM.geta_0T1phi();
            a_1T1=mySM.geta_1T1phi();
            dmT1=mySM.getdmT1phi();
            a_0T2=mySM.geta_0T2phi();
            a_1T2=mySM.geta_1T2phi();
            dmT2=mySM.getdmT2phi();
            a_0T23=mySM.geta_0T23phi();
            a_1T23=mySM.geta_1T23phi();
            dmT23=mySM.getdmT23phi();

            r_1V=mySM.getr_1Vphi();
            r_2V=mySM.getr_2Vphi();
            m_RV=mySM.getm_RVphi();
            m_fit2V=mySM.getm_fit2Vphi();
            r_1A0=mySM.getr_1A0phi();
            r_2A0=mySM.getr_2A0phi();
            m_RA0=mySM.getm_RA0phi();
            m_fit2A0=mySM.getm_fit2A0phi();
            r_2A1=mySM.getr_2A1phi();
            m_fit2A1=mySM.getm_fit2A1phi();
            r_1A2=mySM.getr_1A2phi();
            r_2A2=mySM.getr_2A2phi();
            m_fit2A2=mySM.getm_fit2A2phi();
            r_1T1=mySM.getr_1T1phi();
            r_2T1=mySM.getr_2T1phi();
            m_RT1=mySM.getm_RT1phi();
            m_fit2T1=mySM.getm_fit2T1phi();
            r_2T2=mySM.getr_2T2phi();
            m_fit2T2=mySM.getm_fit2T2phi();
            r_1T3t=mySM.getr_1T3tphi();
            r_2T3t=mySM.getr_2T3tphi();
            m_fit2T3t=mySM.getm_fit2T3tphi();
            
            b= 0.489;
        }
        default:
            std::stringstream out;
            out << vectorM;
            throw std::runtime_error("MVll: vector " + out.str() + " not implemented");
    }
    
    
    h[0]=mySM.getReh_0() + gslpp::complex::i()*mySM.getImh_0();
    h[1]=mySM.getReh_plus() + gslpp::complex::i()*mySM.getImh_plus();
    h[2]=mySM.getReh_minus() + gslpp::complex::i()*mySM.getImh_minus();
    
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

void MVll::checkCache(double qmin, double qmax){
    
    if (MB == k2_cache(0) && MKstar == k2_cache(1) ) {
        k2_updated = 1;
        z_updated = 1;
    } else {
        k2_updated = 0;
        z_updated = 0;
        k2_cache(0) = MB;
        k2_cache(1) = MKstar;
    }
    
    if (Mlep == beta_cache) {
        beta_updated = 1;
    } else {
        beta_updated = 0;
        beta_cache = Mlep;
    }
    
    if (MB == lambda_cache) {
        lambda_updated = k2_updated;
        F_updated = lambda_updated * beta_updated;
    } else {
        lambda_updated = 0;
        F_updated = 0;
        lambda_cache = MB;
    }
    
    if (GF == N_cache(0) && ale == N_cache(1) && MB == N_cache(2) && lambda_t == Nc_cache ) {
        N_updated = 1;
    } else {
        N_updated = 0;
        N_cache(0) = GF;
        N_cache(1) = ale;
        N_cache(2) = MB;
        Nc_cache = lambda_t;
    }
    
    if (r_1A2 == A2_cache(0) && r_2A2 == A2_cache(1) ) {
        A2_updated = 1;
    } else {
        A2_updated = 0;
        A2_cache(0) = r_1A2;
        A2_cache(1) = r_2A2;
    }
    
    if (qmax < CUTOFF) {
        if (r_1V == V_cache(0) && r_2V == V_cache(1) ) {
            V_updated = 1;
        } else {
            V_updated = 0;
            V_cache(0) = r_1V;
            V_cache(1) = r_2V;
        }
        
        if (r_1A0 == A0_cache(0) && r_2A0 == A0_cache(1) ) {
            A0_updated = 1;
        } else {
            A0_updated = 0;
            A0_cache(0) = r_1A0;
            A0_cache(1) = r_2A0;
        }
        
        if ( r_2A1 == A1_cache(0) ) {
            A1_updated = 1;
        } else {
            A1_updated = 0;
            A1_cache(0) = r_2A1;
        }
        
        if (r_1T1 == T1_cache(0) && r_2T1 == T1_cache(1) ) {
            T1_updated = 1;
        } else {
            T1_updated = 0;
            T1_cache(0) = r_1T1;
            T1_cache(1) = r_2T1;
        }
        
        if ( r_2T2 == T2_cache(0) ) {
            T2_updated = 1;
        } else {
            T2_updated = 0;
            T2_cache(0) = r_2T2;
        }
    } else if (qmin >= CUTOFF){
        if (a_0V == V_cache(2) && a_1V == V_cache(3) ) {
            V_updated = z_updated;
        } else {
            V_updated = 0;
            V_cache(2) = a_0V;
            V_cache(3) = a_1V;
        }
        
        if (a_0A0 == A0_cache(2) && a_1A0 == A0_cache(3) ) {
            A0_updated = z_updated;
        } else {
            A0_updated = 0;
            A0_cache(2) = a_0A0;
            A0_cache(3) = a_1A0;
        }
        
        if (a_0A1 == A1_cache(1) && a_1A1 == A1_cache(2) ) {
            A1_updated = z_updated;
        } else {
            A1_updated = 0;
            A1_cache(1) = a_0A1;
            A1_cache(2) = a_1A1;
        }
        
        if (a_0T1 == T1_cache(2) && a_1T1 == T1_cache(3) ) {
            T1_updated = z_updated;
        } else {
            T1_updated = 0;
            T1_cache(2) = a_0T1;
            T1_cache(3) = a_1T1;
        }
        
        if (a_0T2 == T2_cache(1) && a_1T2 == T2_cache(2) ) {
            T2_updated = z_updated;
        } else {
            T2_updated = 0;
            T2_cache(1) = a_0T2;
            T2_cache(2) = a_1T2;
        }
    } else {
        V_updated = 0;
        A0_updated = 0;
        A1_updated = 0;
        T1_updated = 0;
        T2_updated = 0;
    }
    
    if (r_1T3t == T3t_cache(0) && r_2T3t == T3t_cache(1) ) {
        T3t_updated = 1;
    } else {
        T3t_updated = 0;
        T3t_cache(0) = r_1T3t;
        T3t_cache(1) = r_2T3t;
    }
    
    T3_updated = k2_updated * T3t_updated * T2_updated;
    
    VL1_updated = k2_updated * lambda_updated * A1_updated * V_updated;
    VL2_updated = VL1_updated;
    
    TL1_updated = k2_updated * lambda_updated * T1_updated * T2_updated;
    TL2_updated = TL1_updated;
    
    VR1_updated = VL2_updated;
    VR2_updated = VL1_updated;
    
    TR1_updated = TL2_updated;
    TR2_updated = TL1_updated;
    
    if (Mb == SL_cache(0) && Ms == SL_cache(1) ){
        SL_updated = lambda_updated * A0_updated;
        SR_updated = SL_updated;
    } else {
        SL_updated = 0;
        SR_updated = SL_updated;
        SL_cache(0) = Mb;
        SL_cache(1) = Ms;
    }
    
    if (qmax < CUTOFF) {
        VL0_updated = k2_updated * lambda_updated * A1_updated * A2_updated;
        VR0_updated = VL0_updated;
        TL0_updated = k2_updated * lambda_updated * T2_updated * T3_updated;
        TR0_updated = TL0_updated;
    } else if (qmin >= CUTOFF) {
        if (a_0A12 == VL0_cache(0) && a_1A12 == VL0_cache(1) && MKstar == VL0_cache(2) ){
            VL0_updated = z_updated;
            VR0_updated = VL0_updated;
        } else {
            VL0_updated = 0;
            VR0_updated = VL0_updated;
            VL0_cache(0) = a_0A12;
            VL0_cache(1) = a_1A12;
            VL0_cache(2) = MKstar;
        }
        if (a_0T23 == TL0_cache(0) && a_1T23 == TL0_cache(1) ){
            TL0_updated = k2_updated * z_updated;
            TR0_updated = TL0_updated;
        } else {
            TL0_updated = 0;
            TR0_updated = TL0_updated;
            TL0_cache(0) = a_0T23;
            TL0_cache(1) = a_1T23;
            VL0_cache(2) = MKstar;
        }
    } else {
        VL0_updated = 0;
        VR0_updated = 0;
        TL0_updated = 0;
        TR0_updated = 0;
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
    
    if (MB == H_V0cache(0) && Mb == H_V0cache(1) && h[0] == H_V0Ccache) {
        H_V0updated = N_updated * C_9_updated * VL0_updated * C_9p_updated * VR0_updated * C_7_updated * TL0_updated * C_7p_updated * TR0_updated;
    } else {
        H_V0updated = 0;
        H_V0cache(0) = MB;
        H_V0cache(1) = Mb;
        H_V0Ccache = h[0];
    }
    
    if (MB == H_V1cache(0) && Mb == H_V1cache(1) && h[1] == H_V1Ccache) {
        H_V1updated = N_updated * C_9_updated * VL1_updated * C_9p_updated * VR1_updated * C_7_updated * TL1_updated * C_7p_updated * TR1_updated;
    } else {
        H_V1updated = 0;
        H_V1cache(0) = MB;
        H_V1cache(1) = Mb;
        H_V1Ccache = h[1];
    }
    
    if (MB == H_V2cache(0) && Mb == H_V2cache(1) && h[2] == H_V2Ccache) {
        H_V2updated = N_updated * C_9_updated * VL2_updated * C_9p_updated * VR2_updated * C_7_updated * TL2_updated * C_7p_updated * TR2_updated;
    } else {
        H_V2updated = 0;
        H_V2cache(0) = MB;
        H_V2cache(1) = Mb;
        H_V2Ccache = h[2];
    }
    
    H_A0updated = N_updated * C_10_updated * VL0_updated * C_10p_updated * VR0_updated;
    H_A1updated = N_updated * C_10_updated * VL1_updated * C_10p_updated * VR1_updated;
    H_A2updated = N_updated * C_10_updated * VL2_updated * C_10p_updated * VR2_updated;
    
    if (Mb == H_Scache(0) && MW == H_Scache(1)) {
        H_Supdated = N_updated * C_S_updated * SL_updated * C_Sp_updated * SR_updated;
    } else {
        H_Supdated = 0;
        H_Scache(0) = Mb;
        H_Scache(1) = MW;
    }
    
    if (Mb == H_Pcache(0) && MW == H_Pcache(1) && Mlep == H_Pcache(2) && Ms == H_Pcache(3)) {
        H_Pupdated = N_updated * C_P_updated * SL_updated * C_Pp_updated * SR_updated * C_10_updated * C_10p_updated;
    } else {
        H_Pupdated = 0;
        H_Pcache(0) = Mb;
        H_Pcache(1) = MW;
        H_Pcache(2) = Mlep;
        H_Pcache(3) = Ms;
        
    }
    
    I0_updated = F_updated * H_V0updated * H_A0updated * H_Pupdated * beta_updated * H_Supdated;
    I1_updated = F_updated * beta_updated * H_V1updated * H_V2updated * H_A1updated * H_A2updated;
    I2_updated = F_updated * beta_updated * H_V0updated * H_A0updated;
    I3_updated = F_updated * H_V1updated * H_V2updated * H_A1updated * H_A2updated * beta_updated;
    I4_updated = F_updated * H_V1updated * H_V2updated * H_A1updated * H_A2updated;
    I5_updated = F_updated * H_V0updated * H_V1updated * H_V2updated * H_A0updated * H_A1updated * H_A2updated * beta_updated;
    I6_updated = F_updated * H_V1updated * H_V2updated * H_A0updated * H_A1updated * H_A2updated * H_V0updated * beta_updated * H_Supdated;
    I7_updated = I4_updated * beta_updated;
    I8_updated = F_updated * beta_updated * H_Supdated * H_V0updated;
    I9_updated = I6_updated;
    I10_updated = I5_updated;
    I11_updated = I7_updated;
    
    it += 1 ;
    
//    if (I0_updated == 1) std::cout << I0_updated << " I0 " << it << std::endl;
//    if (I1_updated == 1) std::cout << I1_updated << " I1 " << it << std::endl;
//    if (I2_updated == 1) std::cout << I2_updated << " I2 " << it << std::endl;
//    if (I3_updated == 0) std::cout << I3_updated << " I3 " << it << std::endl;
//    if (I4_updated == 1) std::cout << I4_updated << " I4 " << it << std::endl;
//    if (I5_updated == 1) std::cout << I5_updated << " I5 " << it << std::endl;
//    if (I6_updated == 1) std::cout << I6_updated << " I6 " << it << std::endl;
//    if (I7_updated == 1) std::cout << I7_updated << " I7 " << it << std::endl;
//    if (I8_updated == 1) std::cout << I8_updated << " I8 " << it << std::endl;
//    if (I9_updated == 1) std::cout << I9_updated << " I9 " << it << std::endl;
//    if (I10_updated == 1) std::cout << I10_updated << " I10 " << it << std::endl;
//    if (I11_updated == 1) std::cout << I11_updated << " I11 " << it << "\n" <<std::endl;
    
}

/*******************************************************************************
 * Transverse Form Factors                                                     *
 * ****************************************************************************/
double MVll::LCSR_fit1(double q2, double r_1, double r_2, double m_R2, double m_fit2){
    return r_1/( 1. - q2/m_R2 ) + r_2/( 1. - q2/m_fit2 ) ;
}



double MVll::LCSR_fit2(double q2, double r_1, double r_2, double m_fit2){
    return r_1/( 1. - q2/m_fit2 ) + r_2/pow( ( 1. - q2/m_fit2 ) , 2.) ;

}



double MVll::LCSR_fit3(double q2, double r_2, double m_fit2){
    return r_2/( 1. - q2/m_fit2 ) ; 
}



double MVll::z(double q2){
    double t_0 = 12.;
    double t_p=pow(MB + MKstar,2.);
    return ( sqrt(t_p - q2) - sqrt(t_p - t_0) ) / ( sqrt(t_p - q2) + sqrt(t_p - t_0) );
}



double MVll::lat_fit(double q2, double a_0, double a_1, double dm){
    return 1 / (1 - q2/pow(MB + dm,2.)) * ( a_0 + a_1*z(q2) );
}



double MVll::V(double q2){
    if (q2<CUTOFF)
        return LCSR_fit1(q2, r_1V, r_2V, pow(m_RV, 2.), m_fit2V);
    else
        return lat_fit(q2, a_0V, a_1V, dmV);
}



double MVll::A_0(double q2){
    if (q2<CUTOFF)
        return LCSR_fit1(q2, r_1A0, r_2A0, pow(m_RA0, 2.), m_fit2A0);
    else
        return lat_fit(q2, a_0A0, a_1A0, dmA0);
}



double MVll::A_1(double q2){
    if (q2<CUTOFF)
        return LCSR_fit3(q2, r_2A1, m_fit2A1);
    else
        return lat_fit(q2, a_0A1, a_1A1, dmA1);
}



double MVll::A_2(double q2){
    return LCSR_fit2(q2, r_1A2, r_2A2, m_fit2A2);
}



double MVll::T_1(double q2){
    if (q2<CUTOFF)
        return LCSR_fit1(q2, r_1T1, r_2T1, pow(m_RT1, 2.), m_fit2T1);
    else
        return lat_fit(q2, a_0T1, a_1T1, dmT1);
}



double MVll::T_2(double q2){
    if (q2<CUTOFF)
        return LCSR_fit3(q2, r_2T2, m_fit2T2);
    else
        return lat_fit(q2, a_0T2, a_1T2, dmT2);
}



double MVll::T_3tilde(double q2){
    return LCSR_fit2(q2, r_1T3t, r_2T3t, m_fit2T3t);
}



double MVll::T_3(double q2){
    if (q2 < 2.) {
        switch(vectorM){
            case QCD::K_star : return (0.178168 - 0.202)/2. * q2 + 0.202;
            case QCD::PHI : return (0.14667 - 0.175)/2. * q2 + 0.175;
            default:
            std::stringstream out;
            out << vectorM;
            throw std::runtime_error("MVll: vector " + out.str() + " not implemented");
                   
        }
    }
    else return (MB*MB - MKstar*MKstar)/q2*(T_3tilde(q2) - T_2(q2));
}



double MVll::V_L(int i, double q2){
    switch (i){
        case 0:
            if (q2 < CUTOFF)
                return 1. / ( 4.*MKstar*MB*(MB + MKstar)*sqrt(q2) ) * ( pow((MB + MKstar),2.)*(MB*MB - q2 - MKstar*MKstar)*A_1(q2) - lambda(q2)*A_2(q2) );
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



double MVll::V_R(int i, double q2){
    if (i != 0) i=3-i;
    return -V_L(i,q2);
}



double MVll::T_L(int i, double q2){
    switch (i){
        case 0:
            if (q2 < CUTOFF)
                return sqrt(q2)/(4.*MB*MB*MKstar) * ( ( MB*MB+ 3.*MKstar*MKstar - q2 ) * T_2(q2) - lambda(q2) / (MB*MB - MKstar*MKstar) * T_3(q2) );
            else
                return 2*sqrt(q2)*MKstar/MB/(MB + MKstar)*lat_fit(q2, a_0T23, a_1T23, dmT23);
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



double MVll::T_R(int i, double q2){
    if (i != 0) i=3-i;
    return -T_L(i,q2);
}



double MVll::S_L(double q2){
    return -sqrt(lambda(q2))/ ( 2*MB*(Mb + Ms) ) *A_0(q2);
}



double MVll::S_R(double q2){
    return -S_L(q2);
}



/*******************************************************************************
 * Helicity amplitudes                                                         *
 * ****************************************************************************/
complex MVll::N(){
    return -(4.*GF*MB*ale*lambda_t)/(sqrt(2.)*4.*M_PI);
}



gslpp::complex MVll::H_V(int i, double q2, int bar) {
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
                    
    return -gslpp::complex::i()*n*( C_9*V_L(i,q2)
            + C_9p*V_R(i,q2)
            + MB*MB/q2*( 2*Mb/MB*( C_7*T_L(i,q2)
            + C_7p*T_R(i,q2) ) - 16*M_PI*M_PI*h[i] ) );
}



gslpp::complex MVll::H_A(int i, double q2, int bar) {
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
     
    return -gslpp::complex::i()*n*( C_10*V_L(i,q2) 
            + C_10p*V_R(i,q2) );
}



gslpp::complex MVll::H_S(double q2, int bar) {
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
     
    return gslpp::complex::i()*n*Mb/MW*( C_S*S_L(q2) + 
            C_Sp*S_R(q2) );
}



gslpp::complex MVll::H_P(double q2, int bar) {
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
     
    return gslpp::complex::i()*n*( Mb/MW*( C_P*S_L(q2) 
            + C_Pp*S_R(q2) ) 
            + 2.*Mlep*Mb/q2*( C_10*( S_L(q2) - Ms/Mb*S_R(q2) ) 
            + C_10p*( S_R(q2) - Ms/Mb*S_L(q2) ) ) );
}



/*******************************************************************************
 * Angular coefficients                                                         *
 * ****************************************************************************/
double MVll::k2(double q2) {
    return (pow(MB,4.) + q2*q2 + pow(MKstar,4.) -2.*MKstar*MKstar*q2 -2.*MB*MB*(q2 + MKstar*MKstar))/(4.*MB*MB);
}



double MVll::beta(double q2) {
    return sqrt(1-4.*Mlep*Mlep/q2);
}



double MVll::lambda(double q2) {
    return 4.*MB*MB*k2(q2);
}



double MVll::F(double q2, double b_i) {
    return sqrt(lambda(q2))*beta(q2)*q2*b_i/(96.*M_PI*M_PI*M_PI*MB*MB*MB);
}



double MVll::I(int i, double q2, int bar) {

    double Mlep2 = Mlep*Mlep;
    double beta2 = beta(q2)*beta(q2);
    

    switch (i){
        case 0: // I1c
            return F(q2,b)*( ( H_V(0,q2,bar).abs2() + H_A(0,q2,bar).abs2() )/2.  +  H_P(q2,bar).abs2()  +  2.*Mlep2/q2*( H_V(0,q2,bar).abs2() 
                    - H_A(0,q2,bar).abs2() )  + beta2*H_S(q2,bar).abs2() );
        case 1: // I1s
            return F(q2,b)*( (beta2 + 2.)/8.*( H_V(1,q2,bar).abs2() + H_V(2,q2,bar).abs2() + H_A(1,q2,bar).abs2() + H_A(2,q2,bar).abs2() )  +
                            Mlep2/q2*( H_V(1,q2,bar).abs2() + H_V(2,q2,bar).abs2() - H_A(1,q2,bar).abs2() - H_A(2,q2,bar).abs2() ) );
        case 2: // I2c
            return -F(q2,b)*beta2/2.*( H_V(0,q2,bar).abs2() + H_A(0,q2,bar).abs2() );
        case 3: // I2s
            return F(q2,b)*beta2/8.*( H_V(1,q2,bar).abs2() + H_V(2,q2,bar).abs2()  +  H_A(1,q2,bar).abs2() + H_A(2,q2,bar).abs2() );
        case 4: // I3
            return -F(q2,b)/2.*( ( H_V(1,q2,bar)*H_V(2,q2,bar).conjugate() ).real()  +  ( H_A(1,q2,bar)*H_A(2,q2,bar).conjugate() ).real() );
        case 5: // I4
            return F(q2,b)*beta2/4.*( ( (H_V(2,q2,bar) + H_V(1,q2,bar))*H_V(0,q2,bar).conjugate() ).real()  +  ( (H_A(2,q2,bar) + H_A(1,q2,bar))*H_A(0,q2,bar).conjugate() ).real() );
        case 6: // I5
            return F(q2,b)*( beta(q2)/2.*( ( (H_V(2,q2,bar) - H_V(1,q2,bar))*H_A(0,q2,bar).conjugate() ).real()  +  ( (H_A(2,q2,bar) - H_A(1,q2,bar))*H_V(0,q2,bar).conjugate() ).real() )  -
                            beta(q2)*Mlep/sqrt(q2)*( H_S(q2,bar).conjugate()*(H_V(1,q2,bar) + H_V(2,q2,bar)) ).real() );
        case 7: // I6s
            return F(q2,b)*beta(q2)*( H_V(2,q2,bar)*(H_A(2,q2,bar).conjugate()) - H_V(1,q2,bar)*(H_A(1,q2,bar).conjugate()) ).real();
        case 8: // I6c
            return 2.*F(q2,b)*beta(q2)*Mlep/sqrt(q2)*( H_S(q2,bar).conjugate()*H_V(0,q2,bar) ).real();
        case 9: // I7
            return F(q2,b)*( beta(q2)/2.*( ( (H_V(2,q2,bar) + H_V(1,q2,bar))*H_A(0,q2,bar).conjugate() ).imag()  +  ( (H_A(2,q2,bar) + H_A(1,q2,bar))*H_V(0,q2,bar).conjugate() ).imag() )  -
                            beta(q2)*Mlep/sqrt(q2)*( H_S(q2,bar).conjugate()*(H_V(2,q2,bar) - H_V(1,q2,bar)) ).imag() );
        case 10: // I8
            return F(q2,b)*beta2/4.*( ( (H_V(2,q2,bar) - H_V(1,q2,bar))*H_V(0,q2,bar).conjugate() ).imag()  +  ( (H_A(2,q2,bar) - H_A(1,q2,bar))*H_A(0,q2,bar).conjugate() ).imag() );
        case 11: // I9
            return F(q2,b)*beta2/2.*( ( H_V(1,q2,bar)*H_V(2,q2,bar).conjugate() ).imag()  +  ( H_A(1,q2,bar)*H_A(2,q2,bar).conjugate() ).imag() );
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("I: index " + out.str() + " not implemented");
    }
}

double MVll::Sigma(int i, double q2) {
    return (I(i, q2,0) + I(i, q2,1))/2;
}

double MVll::Delta(int i, double q2) {
    return (I(i, q2,0) - I(i, q2,1))/2;
}

double MVll::integrateSigma(int i, double q_min, double q_max){

    checkCache(q_min, q_max);
    
    switch(i){
        case 0:
            if (I0_updated == 0) {
                FS0 = convertToGslFunction( boost::bind( &MVll::getSigma0, &(*this), _1 ) );
                gsl_integration_qags (&FS0, q_min, q_max, 1.e-5, 1.e-3, 50, w_sigma0, &avaSigma0, &errSigma0);
                cacheSigma0[std::make_pair(q_min,q_max)] = avaSigma0;
            }
            return cacheSigma0[std::make_pair(q_min,q_max)];
            break;
        case 1:
            if (I1_updated == 0) {
                FS1 = convertToGslFunction( boost::bind( &MVll::getSigma1, &(*this), _1 ) );
                gsl_integration_qags (&FS1, q_min, q_max, 1.e-5, 1.e-3, 50, w_sigma1, &avaSigma1, &errSigma1);
                cacheSigma1[std::make_pair(q_min,q_max)] = avaSigma1;
            }
            return cacheSigma1[std::make_pair(q_min,q_max)];
            break;
        case 2:
            if (I2_updated == 0) {
                FS2 = convertToGslFunction( boost::bind( &MVll::getSigma2, &(*this), _1 ) );
                gsl_integration_qags (&FS2, q_min, q_max, 1.e-5, 1.e-3, 50, w_sigma2, &avaSigma2, &errSigma2);
                cacheSigma2[std::make_pair(q_min,q_max)] = avaSigma2;
            }
            return cacheSigma2[std::make_pair(q_min,q_max)];
            break;
        case 3:
            if (I3_updated == 0) {
                FS3 = convertToGslFunction( boost::bind( &MVll::getSigma3, &(*this), _1 ) );
                gsl_integration_qags (&FS3, q_min, q_max, 1.e-5, 1.e-3, 50, w_sigma3, &avaSigma3, &errSigma3);
                cacheSigma3[std::make_pair(q_min,q_max)] = avaSigma3;
            } else {
            }
            return cacheSigma3[std::make_pair(q_min,q_max)];
            break;
        case 4:
            if (I4_updated == 0) {
                FS4 = convertToGslFunction( boost::bind( &MVll::getSigma4, &(*this), _1 ) );
                gsl_integration_qags (&FS4, q_min, q_max, 1.e-5, 1.e-3, 50, w_sigma4, &avaSigma4, &errSigma4);
                cacheSigma4[std::make_pair(q_min,q_max)] = avaSigma4;
            }
            return cacheSigma4[std::make_pair(q_min,q_max)];
            break;
        case 5:
            if (I5_updated == 0) {
                FS5 = convertToGslFunction( boost::bind( &MVll::getSigma5, &(*this), _1 ) );
                gsl_integration_qags (&FS5, q_min, q_max, 1.e-5, 1.e-3, 50, w_sigma5, &avaSigma5, &errSigma5);
                cacheSigma5[std::make_pair(q_min,q_max)] = avaSigma5;
            }
            return cacheSigma5[std::make_pair(q_min,q_max)];
            break;
        case 6:
            if (I6_updated == 0) {
                FS6 = convertToGslFunction( boost::bind( &MVll::getSigma6, &(*this), _1 ) );
                gsl_integration_qags (&FS6, q_min, q_max, 1.e-5, 1.e-3, 50, w_sigma6, &avaSigma6, &errSigma6);
                cacheSigma6[std::make_pair(q_min,q_max)] = avaSigma6;
            }
            return cacheSigma6[std::make_pair(q_min,q_max)];
            break;
        case 7:
            if (I7_updated == 0) {
                FS7 = convertToGslFunction( boost::bind( &MVll::getSigma7, &(*this), _1 ) );
                gsl_integration_qags (&FS7, q_min, q_max, 1.e-5, 1.e-3, 50, w_sigma7, &avaSigma7, &errSigma7);
                cacheSigma7[std::make_pair(q_min,q_max)] = avaSigma7;
            }
            return cacheSigma7[std::make_pair(q_min,q_max)];
            break;
        case 9:
            if (I9_updated == 0) {
                FS9 = convertToGslFunction( boost::bind( &MVll::getSigma9, &(*this), _1 ) );
                gsl_integration_qags (&FS9, q_min, q_max, 1.e-5, 1.e-3, 50, w_sigma9, &avaSigma9, &errSigma9);
                cacheSigma9[std::make_pair(q_min,q_max)] = avaSigma9;
            }
            return cacheSigma9[std::make_pair(q_min,q_max)];
            break;
        case 11:
            if (I11_updated == 0) {
                FS11 = convertToGslFunction( boost::bind( &MVll::getSigma11, &(*this), _1 ) );
                gsl_integration_qags (&FS11, q_min, q_max, 1.e-5, 1.e-3, 50, w_sigma11, &avaSigma11, &errSigma11);
                cacheSigma11[std::make_pair(q_min,q_max)] = avaSigma11;
            }
            return cacheSigma11[std::make_pair(q_min,q_max)];
            break;
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("MVll::integrateSigma: index " + out.str() + " not implemented");
    }
    
}

double MVll::integrateDelta(int i, double q_min, double q_max){
    
    checkCache(q_min, q_max);
    
    switch(i){
        case 0:
            if (I0_updated == 0) {
                FD0 = convertToGslFunction( boost::bind( &MVll::getDelta0, &(*this), _1 ) );
                gsl_integration_qags (&FD0, q_min, q_max, 1.e-5, 1.e-3, 50, w_delta0, &avaDelta0, &errDelta0);
                cacheDelta0[std::make_pair(q_min,q_max)] = avaDelta0;
            }
            return cacheDelta0[std::make_pair(q_min,q_max)];
            break;
        case 1:
            if (I1_updated == 0) {
                FD1 = convertToGslFunction( boost::bind( &MVll::getDelta1, &(*this), _1 ) );
                gsl_integration_qags (&FD1, q_min, q_max, 1.e-5, 1.e-3, 50, w_delta1, &avaDelta1, &errDelta1);
                cacheDelta1[std::make_pair(q_min,q_max)] = avaDelta1;
            }
            return cacheDelta1[std::make_pair(q_min,q_max)];
            break;
        case 2:
            if (I2_updated == 0) {
                FD2 = convertToGslFunction( boost::bind( &MVll::getDelta2, &(*this), _1 ) );
                gsl_integration_qags (&FD2, q_min, q_max, 1.e-5, 1.e-3, 50, w_delta2, &avaDelta2, &errDelta2);
                cacheDelta2[std::make_pair(q_min,q_max)] = avaDelta2;
            }
            return cacheDelta2[std::make_pair(q_min,q_max)];
            break;
        case 3:
            if (I3_updated == 0) {
                FD3 = convertToGslFunction( boost::bind( &MVll::getDelta3, &(*this), _1 ) );
                gsl_integration_qags (&FD3, q_min, q_max, 1.e-5, 1.e-3, 50, w_delta3, &avaDelta3, &errDelta3);
                cacheDelta3[std::make_pair(q_min,q_max)] = avaDelta3;
            }
            return cacheDelta3[std::make_pair(q_min,q_max)];
            break;
        case 11:
            if (I11_updated == 0) {
                FD11 = convertToGslFunction( boost::bind( &MVll::getDelta11, &(*this), _1 ) );
                gsl_integration_qags (&FD11, q_min, q_max, 1.e-5, 1.e-3, 50, w_delta11, &avaDelta11, &errDelta11);
                cacheDelta11[std::make_pair(q_min,q_max)] = avaDelta11;
            }
            return cacheDelta11[std::make_pair(q_min,q_max)];
            break;
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("integrateDelta: index " + out.str() + " not implemented"); 
    }
}