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





BKstarll::BKstarll(const StandardModel& SM_i, StandardModel::lepton lep_i) : ThObservable(SM_i), mySM(SM_i),
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

void BKstarll::checkCache(double qmin, double qmax){
    
    if (MB == k2_cache(0) && MKstar == k2_cache(1) ) {
        k2_updated = 1;
        z_updated = 1;
    } else {
        k2_updated = 0;
        z_updated = 0;
        k2_cache(0) = MB;
        k2_cache(1) = MKstar;
    }
    
    if (Mm == beta_cache) {
        beta_updated = 1;
    } else {
        beta_updated = 0;
        beta_cache = Mm;
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
    
    if ((*(allcoeff[LO]) + *(allcoeff[NLO]))(6) == C_7_cache) {
        C_7_updated = 1;
    } else {
        C_7_updated = 0;
        C_7_cache = (*(allcoeff[LO]) + *(allcoeff[NLO]))(6);
    }
    
    if ((*(allcoeff[LO]) + *(allcoeff[NLO]))(8) == C_9_cache) {
        C_9_updated = 1;
    } else {
        C_9_updated = 0;
        C_9_cache = (*(allcoeff[LO]) + *(allcoeff[NLO]))(8);
    }
    
    if ((*(allcoeff[LO]) + *(allcoeff[NLO]))(9) == C_10_cache) {
        C_10_updated = 1;
    } else {
        C_10_updated = 0;
        C_10_cache = (*(allcoeff[LO]) + *(allcoeff[NLO]))(9);
    }
    
    if ((*(allcoeff[LO]) + *(allcoeff[NLO]))(10) == C_S_cache) {
        C_S_updated = 1;
    } else {
        C_S_updated = 0;
        C_S_cache = (*(allcoeff[LO]) + *(allcoeff[NLO]))(10);
    }
    
    if ((*(allcoeff[LO]) + *(allcoeff[NLO]))(11) == C_P_cache) {
        C_P_updated = 1;
    } else {
        C_P_updated = 0;
        C_P_cache = (*(allcoeff[LO]) + *(allcoeff[NLO]))(11);
    }
    
    if ((*(allcoeffprime[LO]) + *(allcoeffprime[NLO]))(6) == C_7p_cache) {
        C_7p_updated = 1;
    } else {
        C_7p_updated = 0;
        C_7p_cache = (*(allcoeffprime[LO]) + *(allcoeffprime[NLO]))(6);
    }
    
    if ((*(allcoeffprime[LO]) + *(allcoeffprime[NLO]))(8) == C_9p_cache) {
        C_9p_updated = 1;
    } else {
        C_9p_updated = 0;
        C_9p_cache = (*(allcoeffprime[LO]) + *(allcoeffprime[NLO]))(8);
    }
    
    if ((*(allcoeffprime[LO]) + *(allcoeffprime[NLO]))(9) == C_10p_cache) {
        C_10p_updated = 1;
    } else {
        C_10p_updated = 0;
        C_10p_cache = (*(allcoeffprime[LO]) + *(allcoeffprime[NLO]))(9);
    }
    
    if ((*(allcoeffprime[LO]) + *(allcoeffprime[NLO]))(10) == C_Sp_cache) {
        C_Sp_updated = 1;
    } else {
        C_Sp_updated = 0;
        C_Sp_cache = (*(allcoeffprime[LO]) + *(allcoeffprime[NLO]))(10);
    }
    
    if ((*(allcoeffprime[LO]) + *(allcoeffprime[NLO]))(11) == C_Pp_cache) {
        C_Pp_updated = 1;
    } else {
        C_Pp_updated = 0;
        C_Pp_cache = (*(allcoeffprime[LO]) + *(allcoeffprime[NLO]))(11);
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
    
    if (Mb == H_Pcache(0) && MW == H_Pcache(1) && Mm == H_Pcache(2) && Ms == H_Pcache(3)) {
        H_Pupdated = N_updated * C_P_updated * SL_updated * C_Pp_updated * SR_updated * C_10_updated * C_10p_updated;
    } else {
        H_Pupdated = 0;
        H_Pcache(0) = Mb;
        H_Pcache(1) = MW;
        H_Pcache(2) = Mm;
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
//    if (I3_updated == 1) std::cout << I3_updated << " I3 " << it << std::endl;
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
    double t_0 = 12.;
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



double BKstarll::F(double q2, double b_i) {
    return sqrt(lambda(q2))*beta(q2)*q2*b_i/(96*M_PI*M_PI*M_PI*MB*MB*MB);
}



double BKstarll::I(int i, double q2, int bar) {

    double Mm2 = Mm*Mm;
    double beta2 = beta(q2)*beta(q2);
    

    switch (i){
        case 0: // I1c
            return F(q2,b)*( ( H_V(0,q2,bar).abs2() + H_A(0,q2,bar).abs2() )/2  +  H_P(q2,bar).abs2()  +  2*Mm2/q2*( H_V(0,q2,bar).abs2() 
                    - H_A(0,q2,bar).abs2() )  + beta2*H_S(q2,bar).abs2() );
        case 1: // I1s
            return F(q2,b)*( (beta2 + 2.)/8.*( H_V(1,q2,bar).abs2() + H_V(2,q2,bar).abs2() + H_A(1,q2,bar).abs2() + H_A(2,q2,bar).abs2() )  +
                            Mm2/q2*( H_V(1,q2,bar).abs2() + H_V(2,q2,bar).abs2() - H_A(1,q2,bar).abs2() - H_A(2,q2,bar).abs2() ) );
        case 2: // I2c
            return -F(q2,b)*beta2/2*( H_V(0,q2,bar).abs2() + H_A(0,q2,bar).abs2() );
        case 3: // I2s
            return F(q2,b)*beta2/8*( H_V(1,q2,bar).abs2() + H_V(2,q2,bar).abs2()  +  H_A(1,q2,bar).abs2() + H_A(2,q2,bar).abs2() );
        case 4: // I3
            return -F(q2,b)/2*( ( H_V(1,q2,bar)*H_V(2,q2,bar).conjugate() ).real()  +  ( H_A(1,q2,bar)*H_A(2,q2,bar).conjugate() ).real() );
        case 5: // I4
            return F(q2,b)*beta2/4*( ( (H_V(2,q2,bar) + H_V(1,q2,bar))*H_V(0,q2,bar).conjugate() ).real()  +  ( (H_A(2,q2,bar) + H_A(1,q2,bar))*H_A(0,q2,bar).conjugate() ).real() );
        case 6: // I5
            return F(q2,b)*( beta(q2)/2*( ( (H_V(2,q2,bar) - H_V(1,q2,bar))*H_A(0,q2,bar).conjugate() ).real()  +  ( (H_A(2,q2,bar) - H_A(1,q2,bar))*H_V(0,q2,bar).conjugate() ).real() )  -
                            beta(q2)*Mm/sqrt(q2)*( H_S(q2,bar).conjugate()*(H_V(1,q2,bar) + H_V(2,q2,bar)) ).real() );
        case 7: // I6s
            return F(q2,b)*beta(q2)*( H_V(2,q2,bar)*(H_A(2,q2,bar).conjugate()) - H_V(1,q2,bar)*(H_A(1,q2,bar).conjugate()) ).real();
        case 8: // I6c
            return  2*F(q2,b)*beta(q2)*Mm/sqrt(q2)*( H_S(q2,bar).conjugate()*H_V(0,q2,bar) ).real();
        case 9: // I7
            return F(q2,b)*( beta(q2)/2*( ( (H_V(2,q2,bar) + H_V(1,q2,bar))*H_A(0,q2,bar).conjugate() ).imag()  +  ( (H_A(2,q2,bar) + H_A(1,q2,bar))*H_V(0,q2,bar).conjugate() ).imag() )  -
                            beta(q2)*Mm/sqrt(q2)*( H_S(q2,bar).conjugate()*(H_V(2,q2,bar) - H_V(1,q2,bar)) ).imag() );
        case 10: // I8
            return F(q2,b)*beta2/4*( ( (H_V(2,q2,bar) - H_V(1,q2,bar))*H_V(0,q2,bar).conjugate() ).imag()  +  ( (H_A(2,q2,bar) - H_A(1,q2,bar))*H_A(0,q2,bar).conjugate() ).imag() );
        case 11: // I9
            return F(q2,b)*beta2/2*( ( H_V(1,q2,bar)*H_V(2,q2,bar).conjugate() ).imag()  +  ( H_A(1,q2,bar)*H_A(2,q2,bar).conjugate() ).imag() );
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("I: index " + out.str() + " not implemented");
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
    checkCache(q_min, q_max);
    
    if (I3_updated == 1) {
        avaSigma3 = sigma3_cache;
    } else {
        F1 = convertToGslFunction( boost::bind( &BKstarll::getSigma3, &(*this), _1 ) );
        gsl_integration_workspace * w_sigma3 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F1, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma3, &avaSigma3, &errSigma3);
        gsl_integration_workspace_free (w_sigma3);
        sigma3_cache = avaSigma3;
    }
    
    if (I4_updated == 1) {
        avaSigma4 = sigma4_cache;
    } else {
        F2 = convertToGslFunction( boost::bind( &BKstarll::getSigma4, &(*this), _1 ) );
        gsl_integration_workspace * w_sigma4 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F2, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma4, &avaSigma4, &errSigma4);
        gsl_integration_workspace_free (w_sigma4);
        sigma4_cache = avaSigma4;
    }
    
    return sigma4_cache/(2.* sigma3_cache);
}


P_2::P_2(const StandardModel& SM_i, StandardModel::lepton lep_i) : BKstarll(SM_i, lep_i) {  
}

double P_2::computeThValue() {
    
    updateParameters();
    double q_min = getBinMin();
    double q_max = getBinMax();
    checkCache(q_min, q_max);
    
    if (I3_updated == 1) {
        avaSigma3 = sigma3_cache;
    } else {
        F1 = convertToGslFunction( boost::bind( &BKstarll::getSigma3, &(*this), _1 ) );
        gsl_integration_workspace * w_sigma3 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F1, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma3, &avaSigma3, &errSigma3);
        gsl_integration_workspace_free (w_sigma3);
        sigma3_cache = avaSigma3;
    }
    
    if (I7_updated == 1) {
        avaSigma7 = sigma7_cache;
    } else {
        F2 = convertToGslFunction( boost::bind( &BKstarll::getSigma7, &(*this), _1 ) );
        gsl_integration_workspace * w_sigma7 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F2, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma7, &avaSigma7, &errSigma7);
        gsl_integration_workspace_free (w_sigma7);
        sigma7_cache = avaSigma7;
    }
    
    return sigma7_cache/(8.*sigma3_cache);
}


P_3::P_3(const StandardModel& SM_i, StandardModel::lepton lep_i) : BKstarll(SM_i, lep_i) {  
}

double P_3::computeThValue() {
    
    updateParameters();
    double q_min = getBinMin();
    double q_max = getBinMax();
    checkCache(q_min, q_max);
    
    if (I3_updated == 1) {
        avaSigma3 = sigma3_cache;
    } else {
        F1 = convertToGslFunction( boost::bind( &BKstarll::getSigma3, &(*this), _1 ) );
        gsl_integration_workspace * w_sigma3 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F1, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma3, &avaSigma3, &errSigma3);
        gsl_integration_workspace_free (w_sigma3);
        sigma3_cache = avaSigma3;
    }
    
    if (I11_updated == 1) {
        avaSigma11 = sigma11_cache;
    } else {
        F2 = convertToGslFunction( boost::bind( &BKstarll::getSigma11, &(*this), _1 ) );
        gsl_integration_workspace * w_sigma11 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F2, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma11, &avaSigma11, &errSigma11);
        gsl_integration_workspace_free (w_sigma11);
        sigma11_cache = avaSigma11;
    }
    
    return -sigma11_cache/(4.*sigma3_cache);
}


P_4Prime::P_4Prime(const StandardModel& SM_i, StandardModel::lepton lep_i) : BKstarll(SM_i, lep_i) {  
}

double P_4Prime::computeThValue() {
    
    updateParameters();
    double q_min = getBinMin();
    double q_max = getBinMax();
    checkCache(q_min, q_max);
    
    if (I2_updated == 1) {
        avaSigma2 = sigma2_cache;
    } else {
        F1 = convertToGslFunction( boost::bind( &BKstarll::getSigma2, &(*this), _1 ) );
        gsl_integration_workspace * w_sigma2 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F1, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma2, &avaSigma2, &errSigma2);
        gsl_integration_workspace_free (w_sigma2);
        sigma2_cache = avaSigma2;
    }
    if (I3_updated == 1) {
        avaSigma3 = sigma3_cache;
    } else {
        F2 = convertToGslFunction( boost::bind( &BKstarll::getSigma3, &(*this), _1 ) );
        gsl_integration_workspace * w_sigma3 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F2, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma3, &avaSigma3, &errSigma3);
        gsl_integration_workspace_free (w_sigma3);
        sigma3_cache = avaSigma3;
    }
    
    if (I5_updated == 1) {
        avaSigma5 = sigma5_cache;
    } else {
        gsl_function F3 = convertToGslFunction( boost::bind( &BKstarll::getSigma5, &(*this), _1 ) );
        gsl_integration_workspace * w_sigma5 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F3, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma5, &avaSigma5, &errSigma5);
        gsl_integration_workspace_free (w_sigma5);
        sigma5_cache = avaSigma5;
    }
    
    return sigma5_cache/sqrt(-sigma2_cache*sigma3_cache);
   
}


P_5Prime::P_5Prime(const StandardModel& SM_i, StandardModel::lepton lep_i) : BKstarll(SM_i, lep_i) {  
}

double P_5Prime::computeThValue() {
    
    updateParameters();
    double q_min = getBinMin();
    double q_max = getBinMax();
    checkCache(q_min, q_max);
    
    if (I2_updated == 1) {
        avaSigma2 = sigma2_cache;
    } else {
        F1 = convertToGslFunction( boost::bind( &BKstarll::getSigma2, &(*this), _1 ) );
        gsl_integration_workspace * w_sigma2 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F1, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma2, &avaSigma2, &errSigma2);
        gsl_integration_workspace_free (w_sigma2);
        sigma2_cache = avaSigma2;
    }
    if (I3_updated == 1) {
        avaSigma3 = sigma3_cache;
    } else {
        F2 = convertToGslFunction( boost::bind( &BKstarll::getSigma3, &(*this), _1 ) );
        gsl_integration_workspace * w_sigma3 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F2, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma3, &avaSigma3, &errSigma3);
        gsl_integration_workspace_free (w_sigma3);
        sigma3_cache = avaSigma3;
    }
    
    if (I6_updated == 1) {
        avaSigma6 = sigma6_cache;
    } else {
        F3 = convertToGslFunction( boost::bind( &BKstarll::getSigma6, &(*this), _1 ) );
        gsl_integration_workspace * w_sigma6 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F3, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma6, &avaSigma6, &errSigma6);
        gsl_integration_workspace_free (w_sigma6);
        sigma6_cache = avaSigma6;
    }
    
    return sigma6_cache/(2.*sqrt(-sigma2_cache*sigma3_cache));
}


P_6Prime::P_6Prime(const StandardModel& SM_i, StandardModel::lepton lep_i) : BKstarll(SM_i, lep_i) {  
}

double P_6Prime::computeThValue() {
    
    updateParameters();
    double q_min = getBinMin();
    double q_max = getBinMax();
    checkCache(q_min, q_max);
    
    if (I2_updated == 1) {
        avaSigma2 = sigma2_cache;
    } else {
        F1 = convertToGslFunction( boost::bind( &BKstarll::getSigma2, &(*this), _1 ) );
        gsl_integration_workspace * w_sigma2 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F1, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma2, &avaSigma2, &errSigma2);
        gsl_integration_workspace_free (w_sigma2);
        sigma2_cache = avaSigma2;
    }
    if (I3_updated == 1) {
        avaSigma3 = sigma3_cache;
    } else {
        F2 = convertToGslFunction( boost::bind( &BKstarll::getSigma3, &(*this), _1 ) );
        gsl_integration_workspace * w_sigma3 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F2, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma3, &avaSigma3, &errSigma3);
        gsl_integration_workspace_free (w_sigma3);
        sigma3_cache = avaSigma3;
    }
    
    if (I9_updated == 1) {
        avaSigma9 = sigma9_cache;
    } else {
        F3 = convertToGslFunction( boost::bind( &BKstarll::getSigma9, &(*this), _1 ) );
        gsl_integration_workspace * w_sigma9 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F3, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma9, &avaSigma9, &errSigma9);
        gsl_integration_workspace_free (w_sigma9);
        sigma9_cache = avaSigma9;
    }
    
    return -sigma9_cache/(2.*sqrt(-sigma2_cache*sigma3_cache));
 
}


GammaPrime::GammaPrime(const StandardModel& SM_i, StandardModel::lepton lep_i) : BKstarll(SM_i, lep_i) {  
}

double GammaPrime::computeGammaPrime(double qmin, double qmax) {
    
    double q_min = qmin;
    double q_max = qmax;
    
    if (I2_updated == 1) {
        avaSigma2 = sigma2_cache;
    } else {
        F1 = convertToGslFunction( boost::bind( &BKstarll::getSigma2, &(*this), _1 ) );
        gsl_integration_workspace * w_sigma2 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F1, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma2, &avaSigma2, &errSigma2);
        gsl_integration_workspace_free (w_sigma2);
        sigma2_cache = avaSigma2;
    }
    if (I3_updated == 1) {
        avaSigma3 = sigma3_cache;
    } else {
        F2 = convertToGslFunction( boost::bind( &BKstarll::getSigma3, &(*this), _1 ) );
        gsl_integration_workspace * w_sigma3 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F2, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma3, &avaSigma3, &errSigma3);
        gsl_integration_workspace_free (w_sigma3);
        sigma3_cache = avaSigma3;
    }
    
    if (I0_updated == 1) {
        avaSigma0 = sigma0_cache;
    } else {
        F3 = convertToGslFunction( boost::bind( &BKstarll::getSigma0, &(*this), _1 ) );
        gsl_integration_workspace * w_sigma0 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F3, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma0, &avaSigma0, &errSigma0);
        gsl_integration_workspace_free (w_sigma0);
        sigma0_cache = avaSigma0;
    }
    if (I1_updated == 1) {
        avaSigma1 = sigma1_cache;
    } else {
        F4 = convertToGslFunction( boost::bind( &BKstarll::getSigma1, &(*this), _1 ) );
        gsl_integration_workspace * w_sigma1 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F4, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma1, &avaSigma1, &errSigma1);
        gsl_integration_workspace_free (w_sigma1);
        sigma1_cache = avaSigma1;
    }
    return ((3.*sigma0_cache - sigma2_cache) + 2.*(3.*sigma1_cache - sigma3_cache))/4.;
}

double GammaPrime::computeThValue(){
    updateParameters();
    double q_min = getBinMin();
    double q_max = getBinMax();
    checkCache(q_min, q_max);
    return computeGammaPrime(q_min, q_max);
}


A_FB::A_FB(const StandardModel& SM_i, StandardModel::lepton lep_i) : GammaPrime(SM_i, lep_i){
}


double A_FB::computeThValue() {
    updateParameters();
    double q_min = getBinMin();
    double q_max = getBinMax();
    checkCache(q_min, q_max);
    
    if (I7_updated == 1) {
        avaSigma7 = sigma7_cache;
        
    } else {
        F5 = convertToGslFunction( boost::bind( &GammaPrime::getSigma7, &(*this), _1 ) );
        gsl_integration_workspace * w_sigma7 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F5, q_min, q_max, 1.e-20, 1.e-5, 1000, w_sigma7, &avaSigma7, &errSigma7);
        gsl_integration_workspace_free (w_sigma7);
        sigma7_cache = avaSigma7;
        
    }
    
    
    return -3. * sigma7_cache / 4. / computeGammaPrime(q_min, q_max);
}


BR_BKstarll::BR_BKstarll(const StandardModel& SM_i, StandardModel::lepton lep_i) : GammaPrime(SM_i, lep_i) {  
}

double BR_BKstarll::computeThValue() {
    updateParameters();
    double q_min = getBinMin();
    double q_max = getBinMax();
    checkCache(q_min, q_max);
    
    return computeGammaPrime(q_min, q_max)/width_Bd;
}


ACP::ACP(const StandardModel& SM_i, StandardModel::lepton lep_i) : GammaPrime(SM_i, lep_i){
}

double ACP::computeThValue() {
    
    updateParameters();
    double q_min = getBinMin();
    double q_max = getBinMax();
    checkCache(q_min, q_max);
    
    
    if (I2_updated == 1) {
        avaDelta2 = delta2_cache;
    } else {
        F1 = convertToGslFunction( boost::bind( &GammaPrime::getDelta2, &(*this), _1 ) );
        gsl_integration_workspace * w_delta2 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F1, q_min, q_max, 1.e-5, 1.e-3, 1000, w_delta2, &avaDelta2, &errDelta2);
        gsl_integration_workspace_free (w_delta2);
        delta2_cache = avaDelta2;
    }
    if (I3_updated == 1) {
        avaDelta3 = delta3_cache;
    } else {
        F2 = convertToGslFunction( boost::bind( &GammaPrime::getDelta3, &(*this), _1 ) );
        gsl_integration_workspace * w_delta3 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F2, q_min, q_max, 1.e-5, 1.e-3, 1000, w_delta3, &avaDelta3, &errDelta3);
        gsl_integration_workspace_free (w_delta3);
        delta3_cache = avaDelta3;
    }
    
    if (I0_updated == 1) {
        avaDelta0 = delta0_cache;
    } else {
        F3 = convertToGslFunction( boost::bind( &GammaPrime::getDelta0, &(*this), _1 ) );
        gsl_integration_workspace * w_delta0 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F3, q_min, q_max, 1.e-5, 1.e-3, 1000, w_delta0, &avaDelta0, &errDelta0);
        gsl_integration_workspace_free (w_delta0);
        delta0_cache = avaDelta0;
    }
    if (I1_updated == 1) {
        avaDelta1 = delta1_cache;
    } else {
        F4 = convertToGslFunction( boost::bind( &GammaPrime::getDelta1, &(*this), _1 ) );
        gsl_integration_workspace * w_delta1 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F4, q_min, q_max, 1.e-5, 1.e-3, 1000, w_delta1, &avaDelta1, &errDelta1);
        gsl_integration_workspace_free (w_delta1);
        delta1_cache = avaDelta1;
    }
    
            
    return (3.*delta0_cache - delta2_cache + 2. * ( 3*delta1_cache - delta3_cache ) )/(4.*computeGammaPrime(q_min, q_max));

}


P3CP::P3CP(const StandardModel& SM_i, StandardModel::lepton lep_i) : BKstarll(SM_i, lep_i) {  
}

double P3CP::computeThValue() {
    
    updateParameters();
    double q_min = getBinMin();
    double q_max = getBinMax();
    checkCache(q_min, q_max);
    
    if (I11_updated == 1) {
        avaDelta11 = delta11_cache;
    } else {
        F1 = convertToGslFunction( boost::bind( &BKstarll::getDelta11, &(*this), _1 ) );
        gsl_integration_workspace * w_delta11 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F1, q_min, q_max, 1.e-5, 1.e-3, 1000, w_delta11, &avaDelta11, &errDelta11);
        gsl_integration_workspace_free (w_delta11);
        delta11_cache = avaDelta11;
    }
    if (I3_updated == 1) {
        avaSigma3 = sigma3_cache;
    } else {
        F2 = convertToGslFunction( boost::bind( &BKstarll::getSigma3, &(*this), _1 ) );
        gsl_integration_workspace * w_sigma3 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F2, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma3, &avaSigma3, &errSigma3);
        gsl_integration_workspace_free (w_sigma3);
        sigma3_cache = avaSigma3;
    }

    return - delta11_cache/(4.*sigma3_cache);

}


F_L::F_L(const StandardModel& SM_i, StandardModel::lepton lep_i) : BKstarll(SM_i, lep_i) {
}

double F_L::computeThValue() {
    
    updateParameters();
    double q_min = getBinMin();
    double q_max = getBinMax();
    checkCache(q_min, q_max);
    
    if (I2_updated == 1) {
        avaSigma2 = sigma2_cache;
    } else {
        F1 = convertToGslFunction( boost::bind( &BKstarll::getSigma2, &(*this), _1 ) );
        gsl_integration_workspace * w_sigma2 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F1, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma2, &avaSigma2, &errSigma2);
        gsl_integration_workspace_free (w_sigma2);
        sigma2_cache = avaSigma2;
    }
    if (I3_updated == 1) {
        avaSigma3 = sigma3_cache;
    } else {
        F2 = convertToGslFunction( boost::bind( &BKstarll::getSigma3, &(*this), _1 ) );
        gsl_integration_workspace * w_sigma3 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F2, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma3, &avaSigma3, &errSigma3);
        gsl_integration_workspace_free (w_sigma3);
        sigma3_cache = avaSigma3;
    }
    
    if (I0_updated == 1) {
        avaSigma0 = sigma0_cache;
    } else {
        F3 = convertToGslFunction( boost::bind( &BKstarll::getSigma0, &(*this), _1 ) );
        gsl_integration_workspace * w_sigma0 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F3, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma0, &avaSigma0, &errSigma0);
        gsl_integration_workspace_free (w_sigma0);
        sigma0_cache = avaSigma0;
    }
    if (I1_updated == 1) {
        avaSigma1 = sigma1_cache;
    } else {
        F4 = convertToGslFunction( boost::bind( &BKstarll::getSigma1, &(*this), _1 ) );
        gsl_integration_workspace * w_sigma1 = gsl_integration_workspace_alloc (1000);
        gsl_integration_qags (&F4, q_min, q_max, 1.e-5, 1.e-3, 1000, w_sigma1, &avaSigma1, &errSigma1);
        gsl_integration_workspace_free (w_sigma1);
        sigma1_cache = avaSigma1;
    }
    
    return (3.*sigma0_cache - sigma2_cache)/(((3.*sigma0_cache - sigma2_cache) + 2.*(3.*sigma1_cache - sigma3_cache)));

}


M_1Prime::M_1Prime(const StandardModel& SM_i, StandardModel::lepton lep_i) : BKstarll(SM_i, lep_i) {  
}

double M_1Prime::computeThValue() {
    
    updateParameters();
    double q_min = getBinMin();
    
    return ( H_V(1,q_min,0).abs2() + H_V(2,q_min,0).abs2() - H_A(1,q_min,0).abs2() - H_A(2,q_min,0).abs2() )/
            ( 2*( H_V(1,q_min,0).abs2() + H_V(2,q_min,0).abs2() + H_A(1,q_min,0).abs2() + H_A(2,q_min,0).abs2() ) );
}


M_2Prime::M_2Prime(const StandardModel& SM_i, StandardModel::lepton lep_i) : BKstarll(SM_i, lep_i) {  
}

double M_2Prime::computeThValue() {
    
    updateParameters();
    double q_min = getBinMin();
    
    return ( q_min/(2*Mm*Mm)*( H_P(q_min,0).abs2() + beta(q_min)*beta(q_min)*H_S(q_min,0).abs2() ) + H_V(0,q_min,0).abs2() - H_A(0,q_min,0).abs2() )/
            ( H_V(0,q_min,0).abs2() + H_A(0,q_min,0).abs2() );  
}