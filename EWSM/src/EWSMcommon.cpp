/* 
 * File:   EWSMcommon.cpp
 * Author: mishima
 */

#include "EWSMcommon.h"
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <gsl/gsl_sf.h>


const int EWSMcommon::CacheSize;


EWSMcommon::EWSMcommon(const StandardModel& SM_i) : SM(SM_i) {
}


//////////////////////////////////////////////////////////////////////// 

void EWSMcommon::SetConstants() {
    Mz = SM.getMz();
    mh = SM.getMHl();    
    Mt = SM.getQuarks(SM.TOP).getMass();
    for (int i=0; i<6; i++) {
        ml[i] = SM.getLeptons((StandardModel::lepton) i).getMass();
        mq[i] = SM.getQuarks((StandardModel::quark) i).getMass();
    }
    
    /* X_t with G_F */
    Xt_GF = SM.getGF()*Mt*Mt/8.0/sqrt(2.0)/M_PI/M_PI;    
    
    /* TEST (should be modified later!!) */
    AlsMt = 0.1074432788;
    
    /* zeta functions */
    zeta2 = gsl_sf_zeta_int(2);
    zeta3 = gsl_sf_zeta_int(3);
    zeta4 = gsl_sf_zeta_int(4);
    zeta5 = gsl_sf_zeta_int(5);

    /* Constants for three-loop contribution */
    double Cl2_Pi_3 = Clausen.Cl2(M_PI/3.0);
    S2 = 4.0/9.0/sqrt(3.0)*Cl2_Pi_3;
    D3 = 6.0*zeta3 - 15.0/4.0*zeta4 - 6.0*Cl2_Pi_3*Cl2_Pi_3;
    //double Li4_1_2 = ??;
    //B4 = 16.0*Li4_1_2 - 4.0*zeta2*log2*log2 + 2.0/3.0*pow(log2,4.0) - 13.0/2.0*zeta4;
    B4 = - 1.76280008707377;
    
    /* Logarithms */
    log2 = log(2.0);
    logMZtoME = log( Mz/SM.getLeptons(SM.ELECTRON).getMass() );
    logMZtoMMU = log( Mz/SM.getLeptons(SM.MU).getMass() );    
    logMZtoMTAU = log( Mz/SM.getLeptons(SM.TAU).getMass() );    
    logMZtoMTOP = log( Mz/Mt );
    logMTOPtoMH = log( Mt/mh );
    
    /* Logarithms, Clausen functions, etc for two-loop QCD corrections */
    double r_QCD2 = Mz*Mz/4.0/Mt/Mt;
    Phi_QCD2 = asin(sqrt(r_QCD2));
    gamma_QCD2 = log(2.0*sqrt(r_QCD2));
    h_QCD2 = log(2.0*sqrt(1.0-r_QCD2));
    gsl_complex OneMinusE2Iphi = gsl_complex_rect(1.0-cos(2.0*Phi_QCD2), 
                                                  -sin(2.0*Phi_QCD2));
    gsl_complex OneMinusE4Iphi = gsl_complex_rect(1.0-cos(4.0*Phi_QCD2), 
                                                  -sin(4.0*Phi_QCD2));
    logV1primeAndA1prime = GSL_REAL(gsl_complex_log(OneMinusE2Iphi))
                           - 2.0*GSL_REAL(gsl_complex_log(OneMinusE4Iphi));
    double Phi= asin(Mz/2.0/Mt);            
    Cl3_2Phi = Clausen.Cl3(2.0*Phi);
    Cl3_4Phi = Clausen.Cl3(4.0*Phi);    
    Cl2_2Phi = Clausen.Cl2(2.0*Phi); 
    Cl2_4Phi = Clausen.Cl2(4.0*Phi);     
    
    /* One-loop functions */
    A0_Mz_Mz = PV.A0(Mz, Mz);
    A0_Mz_mh = PV.A0(Mz, SM.getMHl());    
    B0_Mz_Mz2_mh_Mz = PV.B0(Mz, Mz*Mz, mh, Mz);        
    B0p_Mz_Mz2_mh_Mz = PV.B0p(Mz, Mz*Mz, mh, Mz);    
    for (int i=0; i<6; i++) {
        B0_Mz_Mz2_ml_ml[i] = PV.B0(Mz, Mz*Mz, ml[i], ml[i]);
        B0_Mz_Mz2_mq_mq[i] = PV.B0(Mz, Mz*Mz, mq[i], mq[i]);
        Bf_Mz_Mz2_ml_ml[i] = PV.Bf(Mz, Mz*Mz, ml[i], ml[i]);
        Bf_Mz_Mz2_mq_mq[i] = PV.Bf(Mz, Mz*Mz, mq[i], mq[i]);
        if (ml[i]!=0.0) Bf_Mz_0_ml_ml[i] = PV.Bf(Mz, 0.0, ml[i], ml[i]);
        if (mq[i]!=0.0) Bf_Mz_0_mq_mq[i] = PV.Bf(Mz, 0.0, mq[i], mq[i]);
        B0p_Mz_Mz2_ml_ml[i] = PV.B0p(Mz, Mz*Mz, ml[i], ml[i]);
        B0p_Mz_Mz2_mq_mq[i] = PV.B0p(Mz, Mz*Mz, mq[i], mq[i]);
        Bfp_Mz_Mz2_ml_ml[i] = PV.Bfp(Mz, Mz*Mz, ml[i], ml[i]);
        Bfp_Mz_Mz2_mq_mq[i] = PV.Bfp(Mz, Mz*Mz, mq[i], mq[i]);
    }
    for (int gen=0; gen<3; gen++) {
        B1_Mz_0_ml_mlprime[gen] = PV.B1(Mz, 0.0, ml[2*gen], ml[2*gen+1]);
        B1_Mz_0_mq_mqprime[gen] = PV.B1(Mz, 0.0, mq[2*gen], mq[2*gen+1]);
        Bf_Mz_0_mlprime_ml[gen] = PV.Bf(Mz, 0.0, ml[2*gen+1], ml[2*gen]);
        Bf_Mz_0_mqprime_mq[gen] = PV.Bf(Mz, 0.0, mq[2*gen+1], mq[2*gen]);            
        B1_Mz_0_mlprime_ml[gen] = PV.B1(Mz, 0.0, ml[2*gen+1], ml[2*gen]);
        B1_Mz_0_mqprime_mq[gen] = PV.B1(Mz, 0.0, mq[2*gen+1], mq[2*gen]);
    } 
    
    /* One-loop function in vertex corrections */
    C0_Mz2_0_Mz_0 = PV.C0(Mz*Mz, 0.0, Mz, 0.0);     
}

void EWSMcommon::ComputeForCC(const double Mw_i) {
    SetConstants();
    
    Mw = Mw_i;
    cW2 = Mw*Mw/Mz/Mz;
    sW2 = 1.0 - cW2;    
    
    f_AlphaToGF = sqrt(2.0)*SM.getGF()*pow(Mz,2.0)*sW2*cW2/M_PI/SM.getAle();
    Xt_alpha = Xt_GF/f_AlphaToGF;

    /* Logarithms */
    log_cW2 = log(cW2);

    /* Dilogarithm and Trilogarithm for two-loop QCD corrections */
    Li2_MW2toMTOP2 = gsl_sf_dilog(Mw*Mw/Mt/Mt);
    Li3_MW2toMTOP2 = PolyLog.Li3(Mw*Mw/Mt/Mt);
    Li3_for_F1 = PolyLog.Li3(-Mw*Mw/Mt/Mt/(1.0 - Mw*Mw/Mt/Mt)); 
    
    /* One-loop functions in Delta r */
    A0_Mz_Mw = PV.A0(Mz, Mw);
    B0_Mz_Mw2_Mz_Mw = PV.B0(Mz, Mw*Mw, Mz, Mw);
    B0_Mz_Mw2_0_Mw = PV.B0(Mz, Mw*Mw, 0.0, Mw);
    B0_Mz_Mw2_mh_Mw = PV.B0(Mz, Mw*Mw, mh, Mw);
    B0_Mz_0_Mz_Mw = PV.B0(Mz, 0.0, Mz, Mw);
    B0_Mz_0_0_Mw = PV.B0(Mz, 0.0, 0.0, Mw);
    B0_Mz_0_mh_Mw = PV.B0(Mz, 0.0, mh, Mw);   
    B0_Mz_Mz2_Mw_Mw = PV.B0(Mz, Mz*Mz, Mw, Mw);
    B0p_Mz_0_Mz_Mw = PV.B0p(Mz, 0.0, Mz, Mw);
    B0p_Mz_0_mh_Mw = PV.B0p(Mz, 0.0, mh, Mw);
    for (int gen=0; gen<3; gen++) {
        B1_Mz_Mw2_ml_mlprime[gen] = PV.B1(Mz, Mw*Mw, ml[2*gen], ml[2*gen+1]);
        B1_Mz_Mw2_mq_mqprime[gen] = PV.B1(Mz, Mw*Mw, mq[2*gen], mq[2*gen+1]);
        Bf_Mz_Mw2_mlprime_ml[gen] = PV.Bf(Mz, Mw*Mw, ml[2*gen+1], ml[2*gen]);
        Bf_Mz_Mw2_mqprime_mq[gen] = PV.Bf(Mz, Mw*Mw, mq[2*gen+1], mq[2*gen]);            
        B1_Mz_Mw2_mlprime_ml[gen] = PV.B1(Mz, Mw*Mw, ml[2*gen+1], ml[2*gen]);
        B1_Mz_Mw2_mqprime_mq[gen] = PV.B1(Mz, Mw*Mw, mq[2*gen+1], mq[2*gen]);
    }     
}

void EWSMcommon::ComputeForNC(const double Mw_i) {    
    ComputeForCC(Mw_i);

    /* One-loop function in OneLoopEW::SigmaPrime_ZZ_bos_Mz2() */
    B0p_Mz_Mz2_Mw_Mw = PV.B0p(Mz, Mz*Mz, Mw, Mw);
    
    /* One-loop functions in vertex corrections */    
    B0_Mw_Mz2_Mw_Mw = PV.B0(Mw, Mz*Mz, Mw, Mw);
    B0_Mw_Mz2_Mt_Mt = PV.B0(Mw, Mz*Mz, Mt, Mt);
    C0_Mz2_0_Mw_0 = PV.C0(Mz*Mz, 0.0, Mw, 0.0);
    C0_Mz2_Mt_Mw_Mt = PV.C0(Mz*Mz, Mt, Mw, Mt);
    C0_Mz2_Mw_0_Mw = PV.C0(Mz*Mz, Mw, 0.0, Mw);   
    C0_Mz2_Mw_Mt_Mw = PV.C0(Mz*Mz, Mw, Mt, Mw);   
}

void EWSMcommon::ComputeForRhoWij(const double Mw_i) {
    ComputeForCC(Mw_i);    
    
    /* One-loop functions in vertex corrections to Gamma_W */
    A0_Mw_Mw = PV.A0(Mw, Mw);
    A0_Mw_Mz = PV.A0(Mw, Mz);
    A0_Mw_mh = PV.A0(Mw, mh);
    B0_Mw_Mw2_Mz_Mw = PV.B0(Mw, Mw*Mw, Mz, Mw);
    B0_Mw_Mw2_0_Mw = PV.B0(Mw, Mw*Mw, 0.0, Mw);
    B0_Mw_Mw2_mh_Mw = PV.B0(Mw, Mw*Mw, mh, Mw);
    B0p_Mw_Mw2_Mz_Mw = PV.B0p(Mw, Mw*Mw, Mz, Mw);
    B0p_Mw_Mw2_0_Mw = PV.B0p(Mw, Mw*Mw, 0.0, Mw);
    B0p_Mw_Mw2_mh_Mw = PV.B0p(Mw, Mw*Mw, mh, Mw);
    for (int gen=0; gen<3; gen++) {  
        Bf_Mw_Mw2_mlprime_ml[gen] = PV.Bf(Mw, Mw*Mw, ml[2*gen+1], ml[2*gen]);
        Bf_Mw_Mw2_mqprime_mq[gen] = PV.Bf(Mw, Mw*Mw, mq[2*gen+1], mq[2*gen]);               
        Bfp_Mw_Mw2_mlprime_ml[gen] = PV.Bfp(Mw, Mw*Mw, ml[2*gen+1], ml[2*gen]);
        Bfp_Mw_Mw2_mqprime_mq[gen] = PV.Bfp(Mw, Mw*Mw, mq[2*gen+1], mq[2*gen]);          
        B1p_Mw_Mw2_ml_mlprime[gen] = PV.B1p(Mw, Mw*Mw, ml[2*gen], ml[2*gen+1]);
        B1p_Mw_Mw2_mq_mqprime[gen] = PV.B1p(Mw, Mw*Mw, mq[2*gen], mq[2*gen+1]);
        B1p_Mw_Mw2_mlprime_ml[gen] = PV.B1p(Mw, Mw*Mw, ml[2*gen+1], ml[2*gen]);
        B1p_Mw_Mw2_mqprime_mq[gen] = PV.B1p(Mw, Mw*Mw, mq[2*gen+1], mq[2*gen]);
    }
    C0_Mw2_0_Mz_0 = PV.C0(Mw*Mw, 0.0, Mz, 0.0); 
    C0_Mw2_Mw_0_Mz = PV.C0(Mw*Mw, Mw, 0.0, Mz); 
}


//////////////////////////////////////////////////////////////////////// 

double EWSMcommon::Qf(const StandardModel::lepton l) const {
    switch(l) {
        case StandardModel::NEUTRINO_1:
        case StandardModel::NEUTRINO_2:
        case StandardModel::NEUTRINO_3:
            return (0.0);
        case StandardModel::ELECTRON:
        case StandardModel::MU:
        case StandardModel::TAU:
            return (-1.0);
        default:
            throw "Error in EWSMcommon::Qf()";  
    }
}

double EWSMcommon::Qf(const StandardModel::quark q) const {
    switch(q) {
        case StandardModel::UP:
        case StandardModel::CHARM:
        case StandardModel::TOP:
            return (2.0/3.0);
        case StandardModel::DOWN:
        case StandardModel::STRANGE:
        case StandardModel::BOTTOM:
            return (-1.0/3.0);
        default:
            throw "Error in EWSMcommon::Qf()";  
    }
}

double EWSMcommon::vf(const StandardModel::lepton l) const {
    return ( af(l) - 2.0*Qf(l)*sW2 );
}

double EWSMcommon::vf(const StandardModel::quark q) const {
    return ( af(q) - 2.0*Qf(q)*sW2 );    
}

double EWSMcommon::af(const StandardModel::lepton l) const {
    switch(l) {
        case StandardModel::NEUTRINO_1:
        case StandardModel::NEUTRINO_2:
        case StandardModel::NEUTRINO_3:
            return ( 1.0/2.0 );
        case StandardModel::ELECTRON:
        case StandardModel::MU:
        case StandardModel::TAU:
            return ( -1.0/2.0 );
        default:
            throw "Error in EWSMcommon::af()";  
    }        
}

double EWSMcommon::af(const StandardModel::quark q) const {
    switch(q) {
        case StandardModel::UP:
        case StandardModel::CHARM:
        case StandardModel::TOP:
            return ( 1.0/2.0 );
        case StandardModel::DOWN:
        case StandardModel::STRANGE:
        case StandardModel::BOTTOM:
            return ( -1.0/2.0 );
        default:
            throw "Error in EWSMcommon::af()";  
    }    
}

double EWSMcommon::sigmaf(const StandardModel::lepton l) const {
    return ( 1.0 - 2.0*fabs(Qf(l))*sW2 );
}

double EWSMcommon::sigmaf(const StandardModel::quark q) const {
    return ( 1.0 - 2.0*fabs(Qf(q))*sW2 );    
}

double EWSMcommon::deltaf(const StandardModel::lepton l) const {
    return ( - 2.0*Qf(l)*sW2 );   
}

double EWSMcommon::deltaf(const StandardModel::quark q) const {
    return ( - 2.0*Qf(q)*sW2 );       
}








