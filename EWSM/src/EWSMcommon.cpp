/* 
 * File:   EWSMcommon.cpp
 * Author: mishima
 * 
 * Created on August 31, 2011, 2:30 AM
 */

#include "EWSMcommon.h"
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <gsl/gsl_sf.h>


EWSMcommon::EWSMcommon(const StandardModel& SM_i) : SM(SM_i) {
}

//EWSMcommon::EWSMcommon(const EWSMcommon& orig) {
//}

EWSMcommon::~EWSMcommon() {
}


//////////////////////////////////////////////////////////////////////// 

void EWSMcommon::SetConstants() {
 
    /* zeta functions */
    zeta2 = gsl_sf_zeta_int(2);
    zeta3 = gsl_sf_zeta_int(3);
    zeta4 = gsl_sf_zeta_int(4);
    zeta5 = gsl_sf_zeta_int(5);
        
    /* logarithmic functions */
    log2 = log(2.0);
    logMZtoME = log( SM.getMz()/SM.getLeptons(SM.ELECTRON).getMass() );
    logMZtoMMU = log( SM.getMz()/SM.getLeptons(SM.MU).getMass() );    
    logMZtoMTAU = log( SM.getMz()/SM.getLeptons(SM.TAU).getMass() );    
    logMZtoMTOP = log( SM.getMz()/SM.getQuarks(SM.TOP).getMass() );
    logMTOPtoMH = log( SM.getQuarks(SM.TOP).getMass()/SM.getMHl() );
    
    /* X_t with G_F */
    Xt_GF = SM.getGF()*pow(SM.getQuarks(SM.TOP).getMass(), 2.0)
            /8.0/sqrt(2.0)/M_PI/M_PI;    
    
    /* TEST */
    AlsMt = 0.1074432788;

    double Cl2_Pi_3 = gsl_sf_clausen(M_PI/3.0); /* Clausen function */
    S2 = 4.0/9.0/sqrt(3)*Cl2_Pi_3;
    D3 = 6.0*zeta3 - 15.0/4.0*zeta4 - 6.0*Cl2_Pi_3*Cl2_Pi_3;
    //double Li4_1_2 = ;
    //B4 = 16.0*Li4_1_2 - 4.0*zeta2*log2*log2 + 2.0/3.0*pow(log2,4.0) - 13.0/2.0*zeta4;
    B4 = - 1.76280008707377;
    
}

void EWSMcommon::Compute(const double Mw_i) {
    
    Mw = Mw_i;
    double Mz = SM.getMz();
    cW2 = Mw*Mw/Mz/Mz;
    sW2 = 1.0 - cW2;    
    
    f_AlphaToGF = sqrt(2.0)*SM.getGF()*pow(Mz,2.0)*sW2*cW2/M_PI/SM.getAle();

    /* X_t with alpha(0) */
    Xt_alpha = Xt_GF/f_AlphaToGF;
    
    A0_Mz_Mw = PV.A0(Mz, Mw);
    A0_Mz_Mz = PV.A0(Mz, Mz);
    A0_Mz_mh = PV.A0(Mz, SM.getMHl());
    B0_Mz_Mw2_Mz_Mw = PV.B0(Mz, Mw*Mw, Mz, Mw);
    B0_Mz_Mw2_0_Mw = PV.B0(Mz, Mw*Mw, 0.0, Mw);
    B0_Mz_Mw2_mh_Mw = PV.B0(Mz, Mw*Mw, SM.getMHl(), Mw);
    B0_Mz_0_Mz_Mw = PV.B0(Mz, 0.0, Mz, Mw);
    B0_Mz_0_0_Mw = PV.B0(Mz, 0.0, 0.0, Mw);
    B0_Mz_0_mh_Mw = PV.B0(Mz, 0.0, SM.getMHl(), Mw);   
    B0_Mz_Mz2_Mw_Mw = PV.B0(Mz, Mz*Mz, Mw, Mw);
    B0_Mz_Mz2_mh_Mw = PV.B0(Mz, Mz*Mz, SM.getMHl(), Mw);
    B0_Mz_Mz2_mh_Mz = PV.B0(Mz, Mz*Mz, SM.getMHl(), Mz);    
    
    double ml[6], mq[6];
    for (int i=0; i<6; i++) {
        ml[i] = SM.getLeptons((StandardModel::lepton) i).getMass();
        mq[i] = SM.getQuarks((StandardModel::quark) i).getMass();
        B0_Mz_Mz2_ml_ml[i] = PV.B0(Mz, Mz*Mz, ml[i], ml[i]);
        B0_Mz_Mz2_mq_mq[i] = PV.B0(Mz, Mz*Mz, mq[i], mq[i]);
        Bf_Mz_Mz2_ml_ml[i] = PV.Bf(Mz, Mz*Mz, ml[i], ml[i]);
        Bf_Mz_Mz2_mq_mq[i] = PV.Bf(Mz, Mz*Mz, mq[i], mq[i]);
    }
    for (int gen=0; gen<3; gen++) {
        Bf_Mz_Mw2_ml_mlprime[gen] = PV.Bf(Mz, Mw*Mw, ml[2*gen], ml[2*gen+1]);
        Bf_Mz_Mw2_mq_mqprime[gen] = PV.Bf(Mz, Mw*Mw, mq[2*gen], mq[2*gen+1]);            
        B1_Mz_Mw2_ml_mlprime[gen] = PV.B1(Mz, Mw*Mw, ml[2*gen], ml[2*gen+1]);
        B1_Mz_Mw2_mq_mqprime[gen] = PV.B1(Mz, Mw*Mw, mq[2*gen], mq[2*gen+1]);
        Bf_Mz_Mw2_mlprime_ml[gen] = PV.Bf(Mz, Mw*Mw, ml[2*gen+1], ml[2*gen]);
        Bf_Mz_Mw2_mqprime_mq[gen] = PV.Bf(Mz, Mw*Mw, mq[2*gen+1], mq[2*gen]);            
        B1_Mz_Mw2_mlprime_ml[gen] = PV.B1(Mz, Mw*Mw, ml[2*gen+1], ml[2*gen]);
        B1_Mz_Mw2_mqprime_mq[gen] = PV.B1(Mz, Mw*Mw, mq[2*gen+1], mq[2*gen]);
    }    
    
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








