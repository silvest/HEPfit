/* 
 * File:   EWSMcache.cpp
 * Author: mishima
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <gsl/gsl_sf.h>
#include "EWSMcache.h"


EWSMcache::EWSMcache(const StandardModel& SM_i, bool bDebug_i) : SM(SM_i) {
    bDebug = bDebug_i;
    
    log2 = log(2.0);
    
    /* zeta functions */
    zeta2 = gsl_sf_zeta_int(2);
    zeta3 = gsl_sf_zeta_int(3);
    zeta4 = gsl_sf_zeta_int(4);
    zeta5 = gsl_sf_zeta_int(5);

    /* Constants for three-loop contribution */
    double Cl2_Pi_3 = Clausen.Cl2(M_PI/3.0);
    S2 = 4.0/9.0/sqrt(3.0)*Cl2_Pi_3;
    D3 = 6.0*zeta3 - 15.0/4.0*zeta4 - 6.0*Cl2_Pi_3*Cl2_Pi_3;
    B4 = - 1.76280008707377;
    //double Li4_1_2 = ??;
    //B4 = 16.0*Li4_1_2 - 4.0*zeta2*log2*log2 + 2.0/3.0*pow(log2,4.0) - 13.0/2.0*zeta4;

}


//////////////////////////////////////////////////////////////////////// 

int EWSMcache::CacheCheck(const double cache[][CacheSize], 
                          const int NumPar, const double params[]) const {
    bool bCache;
    for(int i=0; i<CacheSize; i++) {
        bCache = true;
        for(int j=0; j<NumPar; j++)
            bCache &= (params[j] == cache[j][i]);
        if (bCache) return i;
    }
    return -1;
}


int EWSMcache::CacheCheck(const complex cache[][CacheSize], 
                          const int NumPar, const double params[]) const {
    bool bCache;
    for(int i=0; i<CacheSize; i++) {
        bCache = true;
        for(int j=0; j<NumPar; j++)
            bCache &= (params[j] == cache[j][i].real());
        if (bCache) return i;
    }
    return -1;
}


void EWSMcache::CacheShift(double cache[][CacheSize], const int NumPar, 
                           const double params[], const double newResult) const {
    // shift old parameters and result
    for(int i=CacheSize-1; i>0; i--)
        for(int j=0; j<NumPar+1; j++)
            cache[j][i] = cache[j][i-1];

    // store new parameters and result
    for(int j=0; j<NumPar; j++) {
        cache[j][0] = params[j];
        cache[NumPar][0] = newResult;
    }
}


void EWSMcache::CacheShift(complex cache[][CacheSize], const int NumPar, 
                           const double params[], const complex newResult) const {
    // shift old parameters and result
    for(int i=CacheSize-1; i>0; i--)
        for(int j=0; j<NumPar+1; j++)
            cache[j][i] = cache[j][i-1];

    // store new parameters and result
    for(int j=0; j<NumPar; j++) {
        cache[j][0] = complex(params[j], 0.0, false);
        cache[NumPar][0] = newResult;
    }
}


//////////////////////////////////////////////////////////////////////// 

double EWSMcache::ml(const StandardModel::lepton l) const {
    return ( SM.getLeptons(l).getMass() );
}


double EWSMcache::mq(const StandardModel::quark q) const {
    if (q==StandardModel::TOP) 
        return Mt();
    else {
        if (!bDebug) {
            return ( SM.getQuarks(q).getMass() );            
        } else {
            double mq_fixed;
            switch(q) {
                case StandardModel::UP:
                    mq_fixed = 0.062;
                    break;
                case StandardModel::DOWN:
                    mq_fixed = 0.083;
                    break;
                case StandardModel::CHARM:
                    mq_fixed = 1.50;
                    break;
                case StandardModel::STRANGE:
                    mq_fixed = 0.215;
                    break;
                case StandardModel::BOTTOM:
                    mq_fixed = 4.70;
                    break;
                default:
                    throw "Error in EWSMcache::mq()";  
            }
            return mq_fixed; // for debug
        }
    }
}


double EWSMcache::Mt() const {
    return ( SM.getMtpole() );
}


double EWSMcache::alsMz() const {
    return ( SM.getAlsMz() );
}


double EWSMcache::GF() const {
    return ( SM.getGF() );
}


double EWSMcache::ale() const {
    return ( SM.getAle() );
}


double EWSMcache::dAle5Mz() const {
    return ( SM.getDAle5Mz() );
}


double EWSMcache::Mz() const {
    return ( SM.getMz() );
}


double EWSMcache::mh() const {
    return ( SM.getMHl() );
}


double EWSMcache::Mw(const double Mw_i) const {
    return ( Mw_i );
}


double EWSMcache::cW2(const double Mw_i) const {
    return ( Mw(Mw_i)*Mw(Mw_i)/Mz()/Mz() );
}


double EWSMcache::sW2(const double Mw_i) const {
    return ( 1.0 - cW2(Mw_i) );
}


double EWSMcache::Ql(const StandardModel::lepton l) const {
    return ( SM.getLeptons(l).getCharge() );
}        
 

double EWSMcache::Qq(const StandardModel::quark q) const {
    return ( SM.getQuarks(q).getCharge() );
}


double EWSMcache::vl(const StandardModel::lepton l, const double Mw_i) const {
    return ( al(l) - 2.0*Ql(l)*sW2(Mw_i) );
}


double EWSMcache::vq(const StandardModel::quark q, const double Mw_i) const {
    return ( aq(q) - 2.0*Qq(q)*sW2(Mw_i) );
}


double EWSMcache::al(const StandardModel::lepton l) const {
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
            throw "Error in EWSMcache::al()";  
    }    
}


double EWSMcache::aq(const StandardModel::quark q) const {
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
            throw "Error in EWSMcache::aq()";  
    }    
}


double EWSMcache::sigmal(const StandardModel::lepton l, const double Mw_i) const {
    return ( 1.0 - 2.0*fabs(Ql(l))*sW2(Mw_i) );
}


double EWSMcache::sigmaq(const StandardModel::quark q, const double Mw_i) const {
    return ( 1.0 - 2.0*fabs(Qq(q))*sW2(Mw_i) );
}


double EWSMcache::deltal(const StandardModel::lepton l, const double Mw_i) const {
    return ( - 2.0*Ql(l)*sW2(Mw_i) );   
}


double EWSMcache::deltaq(const StandardModel::quark q, const double Mw_i) const {
    return ( - 2.0*Qq(q)*sW2(Mw_i) );   
}


double EWSMcache::f_AlphaToGF(const double Mw_i) const {
    return ( sqrt(2.0)*GF()*pow(Mz(),2.0)*sW2(Mw_i)*cW2(Mw_i)/M_PI/ale() );
}


double EWSMcache::Xt_GF() const {
    return ( GF()*Mt()*Mt()/8.0/sqrt(2.0)/M_PI/M_PI );
}
 

double EWSMcache::Xt_alpha(const double Mw_i) const {
    return ( Xt_GF()/f_AlphaToGF(Mw_i) );
}


double EWSMcache::alsMt() const {
    if (!bDebug) {
        
        /* Write codes */

        
        
        
        return ( 0.1074432788 );
    } else 
        return ( 0.1074432788 );// for debug
}


double EWSMcache::Phi_QCD2() const {
    double r_QCD2 = Mz()*Mz()/4.0/Mt()/Mt();
    return ( asin(sqrt(r_QCD2)) );
}
    
double EWSMcache::gamma_QCD2() const {
    double r_QCD2 = Mz()*Mz()/4.0/Mt()/Mt();
    return ( log(2.0*sqrt(r_QCD2)) );
}


double EWSMcache::h_QCD2() const {
    double r_QCD2 = Mz()*Mz()/4.0/Mt()/Mt();
    return ( log(2.0*sqrt(1.0-r_QCD2)) );
}
    
    
double EWSMcache::logV1primeAndA1prime() const {
    gsl_complex OneMinusE2Iphi = gsl_complex_rect(1.0-cos(2.0*Phi_QCD2()), 
                                                  -sin(2.0*Phi_QCD2()));
    gsl_complex OneMinusE4Iphi = gsl_complex_rect(1.0-cos(4.0*Phi_QCD2()), 
                                                  -sin(4.0*Phi_QCD2()));
    return ( GSL_REAL(gsl_complex_log(OneMinusE2Iphi))
             - 2.0*GSL_REAL(gsl_complex_log(OneMinusE4Iphi)) );
}
    

double EWSMcache::Cl3_2Phi() const {
    double Phi= asin(Mz()/2.0/Mt());
    return ( Clausen.Cl3(2.0*Phi) );
}


double EWSMcache::Cl3_4Phi() const {
    double Phi= asin(Mz()/2.0/Mt());
    return ( Clausen.Cl3(4.0*Phi) );
}


double EWSMcache::Cl2_2Phi() const {
    double Phi= asin(Mz()/2.0/Mt());
    return ( Clausen.Cl2(2.0*Phi) );
}


double EWSMcache::Cl2_4Phi() const {
    double Phi= asin(Mz()/2.0/Mt());
    return ( Clausen.Cl2(4.0*Phi) );
}


//////////////////////////////////////////////////////////////////////// 

double EWSMcache::logMZtoME() const {
    int NumPar = 2;
    double params[] = {Mz(), ml(SM.ELECTRON)};

    int i = CacheCheck(logMZtoME_cache, NumPar, params);
    if (i>=0) {
        return ( logMZtoME_cache[NumPar][i] );
    } else {
        double newResult = log( Mz()/ml(SM.ELECTRON) );
        CacheShift(logMZtoME_cache, NumPar, params, newResult);
        return newResult;
    }
}


double EWSMcache::logMZtoMMU() const {
    int NumPar = 2;
    double params[] = {Mz(), ml(SM.MU)};

    int i = CacheCheck(logMZtoMMU_cache, NumPar, params);
    if (i>=0) {
        return ( logMZtoMMU_cache[NumPar][i] );
    } else {
        double newResult = log( Mz()/ml(SM.MU) );
        CacheShift(logMZtoMMU_cache, NumPar, params, newResult);
        return newResult;
    }
}
    
    
double EWSMcache::logMZtoMTAU() const {
    int NumPar = 2;
    double params[] = {Mz(), ml(SM.TAU)};

    int i = CacheCheck(logMZtoMTAU_cache, NumPar, params);
    if (i>=0) {
        return ( logMZtoMTAU_cache[NumPar][i] );
    } else {
        double newResult = log( Mz()/ml(SM.TAU) );
        CacheShift(logMZtoMTAU_cache, NumPar, params, newResult);
        return newResult;
    }
}
    

double EWSMcache::logMZtoMTOP() const {
    int NumPar = 2;
    double params[] = {Mz(), Mt()};

    int i = CacheCheck(logMZtoMTOP_cache, NumPar, params);
    if (i>=0) {
        return ( logMZtoMTOP_cache[NumPar][i] );
    } else {
        double newResult = log( Mz()/Mt() );
        CacheShift(logMZtoMTOP_cache, NumPar, params, newResult);
        return newResult;
    }
}
    
    
double EWSMcache::logMTOPtoMH() const {
    int NumPar = 2;
    double params[] = {Mt(), mh()};

    int i = CacheCheck(logMTOPtoMH_cache, NumPar, params);
    if (i>=0) {
        return ( logMTOPtoMH_cache[NumPar][i] );
    } else {
        double newResult = log( Mt()/mh() );
        CacheShift(logMTOPtoMH_cache, NumPar, params, newResult);
        return newResult;
    }
}

    
double EWSMcache::log_cW2(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    int i = CacheCheck(log_cW2_cache, NumPar, params);
    if (i>=0) {
        return ( log_cW2_cache[NumPar][i] );
    } else {
        double newResult = log(cW2(Mw_i));
        CacheShift(log_cW2_cache, NumPar, params, newResult);
        return newResult;
    }
}    


double EWSMcache::Li2_MW2toMTOP2(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mw(Mw_i), Mt()};

    int i = CacheCheck(Li2_MW2toMTOP2_cache, NumPar, params);
    if (i>=0) {
        return ( Li2_MW2toMTOP2_cache[NumPar][i] );
    } else {
        double newResult = gsl_sf_dilog(Mw(Mw_i)*Mw(Mw_i)/Mt()/Mt());
        CacheShift(Li2_MW2toMTOP2_cache, NumPar, params, newResult);
        return newResult;
    }
}


double EWSMcache::Li3_MW2toMTOP2(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mw(Mw_i), Mt()};

    int i = CacheCheck(Li3_MW2toMTOP2_cache, NumPar, params);
    if (i>=0) {
        return ( Li3_MW2toMTOP2_cache[NumPar][i] );
    } else {
        double newResult = PolyLog.Li3(Mw(Mw_i)*Mw(Mw_i)/Mt()/Mt());
        CacheShift(Li3_MW2toMTOP2_cache, NumPar, params, newResult);
        return newResult;
    }
}


double EWSMcache::Li3_for_F1(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mw(Mw_i), Mt()};

    int i = CacheCheck(Li3_for_F1_cache, NumPar, params);
    if (i>=0) {
        return ( Li3_for_F1_cache[NumPar][i] );
    } else {
        double tmp = Mw(Mw_i)*Mw(Mw_i)/Mt()/Mt();
        double newResult = PolyLog.Li3(-tmp/(1.0 - tmp));
        CacheShift(Li3_for_F1_cache, NumPar, params, newResult);
        return newResult;
    }
}


double EWSMcache::A0_Mz_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    int i = CacheCheck(A0_Mz_Mw_cache, NumPar, params);
    if (i>=0) {
        return ( A0_Mz_Mw_cache[NumPar][i] );
    } else {
        double newResult = PV.A0(Mz(), Mw(Mw_i));
        CacheShift(A0_Mz_Mw_cache, NumPar, params, newResult);
        return newResult;
    }    
}


double EWSMcache::A0_Mz_mh() const {
    int NumPar = 2;
    double params[] = {Mz(), mh()};

    int i = CacheCheck(A0_Mz_mh_cache, NumPar, params);
    if (i>=0) {
        return ( A0_Mz_mh_cache[NumPar][i] );
    } else {
        double newResult = PV.A0(Mz(), mh());
        CacheShift(A0_Mz_mh_cache, NumPar, params, newResult);
        return newResult;
    }      
}


double EWSMcache::A0_Mw_Mz(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mw(Mw_i), Mz()};

    int i = CacheCheck(A0_Mw_Mz_cache, NumPar, params);
    if (i>=0) {
        return ( A0_Mw_Mz_cache[NumPar][i] );
    } else {
        double newResult = PV.A0(Mw(Mw_i), Mz());
        CacheShift(A0_Mw_Mz_cache, NumPar, params, newResult);
        return newResult;
    }         
}


double EWSMcache::A0_Mw_mh(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mw(Mw_i), mh()};

    int i = CacheCheck(A0_Mw_mh_cache, NumPar, params);
    if (i>=0) {
        return ( A0_Mw_mh_cache[NumPar][i] );
    } else {
        double newResult = PV.A0(Mw(Mw_i), mh());
        CacheShift(A0_Mw_mh_cache, NumPar, params, newResult);
        return newResult;
    }     
}


double EWSMcache::A0_Mz_Mz() const {
    int NumPar = 1;
    double params[] = {Mz()};

    int i = CacheCheck(A0_Mz_Mz_cache, NumPar, params);
    if (i>=0) {
        return ( A0_Mz_Mz_cache[NumPar][i] );
    } else {
        double newResult = PV.A0(Mz(), Mz());
        CacheShift(A0_Mz_Mz_cache, NumPar, params, newResult);
        return newResult;
    }     
}


double EWSMcache::A0_Mw_Mw(const double Mw_i) const {
    int NumPar = 1;
    double params[] = {Mw(Mw_i)};

    int i = CacheCheck(A0_Mw_Mw_cache, NumPar, params);
    if (i>=0) {
        return ( A0_Mw_Mw_cache[NumPar][i] );
    } else {
        double newResult = PV.A0(Mw(Mw_i), Mw(Mw_i));
        CacheShift(A0_Mw_Mw_cache, NumPar, params, newResult);
        return newResult;
    }       
}


complex EWSMcache::B0_Mz_Mw2_mh_Mw(const double Mw_i) const {
    int NumPar = 3;
    double params[] = {Mz(), Mw(Mw_i), mh()};

    int i = CacheCheck(B0_Mz_Mw2_mh_Mw_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mz_Mw2_mh_Mw_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0(Mz(), Mw(Mw_i)*Mw(Mw_i), mh(), Mw(Mw_i));
        CacheShift(B0_Mz_Mw2_mh_Mw_cache, NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::B0_Mz_0_mh_Mw(const double Mw_i) const {
    int NumPar = 3;
    double params[] = {Mz(), mh(), Mw(Mw_i)};

    int i = CacheCheck(B0_Mz_0_mh_Mw_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mz_0_mh_Mw_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0(Mz(), 0.0, mh(), Mw(Mw_i));
        CacheShift(B0_Mz_0_mh_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mw_Mz2_Mt_Mt(const double Mw_i) const {
    int NumPar = 3;
    double params[] = {Mw(Mw_i), Mz(), Mt()};

    int i = CacheCheck(B0_Mw_Mz2_Mt_Mt_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mw_Mz2_Mt_Mt_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0(Mw(Mw_i), Mz()*Mz(), Mt(), Mt());
        CacheShift(B0_Mw_Mz2_Mt_Mt_cache, NumPar, params, newResult);
        return newResult;
    }    
}


complex EWSMcache::B0_Mz_Mz2_Mw_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    int i = CacheCheck(B0_Mz_Mz2_Mw_Mw_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mz_Mz2_Mw_Mw_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0(Mz(), Mz()*Mz(), Mw(Mw_i), Mw(Mw_i));
        CacheShift(B0_Mz_Mz2_Mw_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mz_Mz2_mh_Mz() const {
    int NumPar = 2;
    double params[] = {Mz(), mh()};

    int i = CacheCheck(B0_Mz_Mz2_mh_Mz_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mz_Mz2_mh_Mz_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0(Mz(), Mz()*Mz(), mh(), Mz());
        CacheShift(B0_Mz_Mz2_mh_Mz_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mz_Mw2_Mz_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    int i = CacheCheck(B0_Mz_Mw2_Mz_Mw_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mz_Mw2_Mz_Mw_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0(Mz(), Mw(Mw_i)*Mw(Mw_i), Mz(), Mw(Mw_i));
        CacheShift(B0_Mz_Mw2_Mz_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mz_Mw2_0_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    int i = CacheCheck(B0_Mz_Mw2_0_Mw_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mz_Mw2_0_Mw_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0(Mz(), Mw(Mw_i)*Mw(Mw_i), 0.0, Mw(Mw_i));
        CacheShift(B0_Mz_Mw2_0_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mz_0_Mz_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    int i = CacheCheck(B0_Mz_0_Mz_Mw_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mz_0_Mz_Mw_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0(Mz(), 0.0, Mz(), Mw(Mw_i));
        CacheShift(B0_Mz_0_Mz_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mz_0_0_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    int i = CacheCheck(B0_Mz_0_0_Mw_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mz_0_0_Mw_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0(Mz(), 0.0, 0.0, Mw(Mw_i));
        CacheShift(B0_Mz_0_0_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mw_Mz2_Mw_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mw(Mw_i), Mz()};

    int i = CacheCheck(B0_Mw_Mz2_Mw_Mw_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mw_Mz2_Mw_Mw_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0(Mw(Mw_i), Mz()*Mz(), Mw(Mw_i), Mw(Mw_i));
        CacheShift(B0_Mw_Mz2_Mw_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mw_Mw2_Mz_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mw(Mw_i), Mz()};

    int i = CacheCheck(B0_Mw_Mw2_Mz_Mw_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mw_Mw2_Mz_Mw_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), Mz(), Mw(Mw_i));
        CacheShift(B0_Mw_Mw2_Mz_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mw_Mw2_mh_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mw(Mw_i), mh()};

    int i = CacheCheck(B0_Mw_Mw2_mh_Mw_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mw_Mw2_mh_Mw_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), mh(), Mw(Mw_i));
        CacheShift(B0_Mw_Mw2_mh_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mw_Mw2_0_Mw(const double Mw_i) const {
    int NumPar = 1;
    double params[] = {Mw(Mw_i)};
    
    int i = CacheCheck(B0_Mw_Mw2_0_Mw_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mw_Mw2_0_Mw_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), 0.0, Mw(Mw_i));
        CacheShift(B0_Mw_Mw2_0_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mz_Mz2_ml_ml(const StandardModel::lepton l) const {
    int NumPar = 2;
    double params[] = {Mz(), ml(l)};

    int i = CacheCheck(B0_Mz_Mz2_ml_ml_cache[l], NumPar, params);
    if (i>=0) {
        return ( B0_Mz_Mz2_ml_ml_cache[l][NumPar][i] );
    } else {
        complex newResult = PV.B0(Mz(), Mz()*Mz(), ml(l), ml(l));
        CacheShift(B0_Mz_Mz2_ml_ml_cache[l], NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mz_Mz2_mq_mq(const StandardModel::quark q) const {
    int NumPar = 2;
    double params[] = {Mz(), mq(q)};

    int i = CacheCheck(B0_Mz_Mz2_mq_mq_cache[q], NumPar, params);
    if (i>=0) {
        return ( B0_Mz_Mz2_mq_mq_cache[q][NumPar][i] );
    } else {
        complex newResult = PV.B0(Mz(), Mz()*Mz(), mq(q), mq(q));
        CacheShift(B0_Mz_Mz2_mq_mq_cache[q], NumPar, params, newResult);
        return newResult;
    }    
}


complex EWSMcache::B0p_Mz_0_mh_Mw(const double Mw_i) const {
    int NumPar = 3;
    double params[] = {Mz(), mh(), Mw(Mw_i)};

    int i = CacheCheck(B0p_Mz_0_mh_Mw_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_Mz_0_mh_Mw_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0p(Mz(), 0.0, mh(), Mw(Mw_i));
        CacheShift(B0p_Mz_0_mh_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0p_Mz_Mz2_mh_Mz() const {
    int NumPar = 2;
    double params[] = {Mz(), mh()};

    int i = CacheCheck(B0p_Mz_Mz2_mh_Mz_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_Mz_Mz2_mh_Mz_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0p(Mz(), Mz()*Mz(), mh(), Mz());
        CacheShift(B0p_Mz_Mz2_mh_Mz_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0p_Mz_0_Mz_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    int i = CacheCheck(B0p_Mz_0_Mz_Mw_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_Mz_0_Mz_Mw_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0p(Mz(), 0.0, Mz(), Mw(Mw_i));
        CacheShift(B0p_Mz_0_Mz_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0p_Mz_Mz2_Mw_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    int i = CacheCheck(B0p_Mz_Mz2_Mw_Mw_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_Mz_Mz2_Mw_Mw_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0p(Mz(), Mz()*Mz(), Mw(Mw_i), Mw(Mw_i));
        CacheShift(B0p_Mz_Mz2_Mw_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0p_Mw_Mw2_Mz_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mw(Mw_i), Mz()};

    int i = CacheCheck(B0p_Mw_Mw2_Mz_Mw_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_Mw_Mw2_Mz_Mw_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0p(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), Mz(), Mw(Mw_i));
        CacheShift(B0p_Mw_Mw2_Mz_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0p_Mw_Mw2_mh_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mw(Mw_i), mh()};

    int i = CacheCheck(B0p_Mw_Mw2_mh_Mw_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_Mw_Mw2_mh_Mw_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0p(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), mh(), Mw(Mw_i));
        CacheShift(B0p_Mw_Mw2_mh_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0p_Mw_Mw2_0_Mw(const double Mw_i) const {
    int NumPar = 1;
    double params[] = {Mw(Mw_i)};

    int i = CacheCheck(B0p_Mw_Mw2_0_Mw_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_Mw_Mw2_0_Mw_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0p(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), 0.0, Mw(Mw_i));
        CacheShift(B0p_Mw_Mw2_0_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0p_Mz_Mz2_ml_ml(const StandardModel::lepton l) const {
    int NumPar = 2;
    double params[] = {Mz(), ml(l)};

    int i = CacheCheck(B0p_Mz_Mz2_ml_ml_cache[l], NumPar, params);
    if (i>=0) {
        return ( B0p_Mz_Mz2_ml_ml_cache[l][NumPar][i] );
    } else {
        complex newResult = PV.B0p(Mz(), Mz()*Mz(), ml(l), ml(l));
        CacheShift(B0p_Mz_Mz2_ml_ml_cache[l], NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0p_Mz_Mz2_mq_mq(const StandardModel::quark q) const {
    int NumPar = 2;
    double params[] = {Mz(), mq(q)};

    int i = CacheCheck(B0p_Mz_Mz2_mq_mq_cache[q], NumPar, params);
    if (i>=0) {
        return ( B0p_Mz_Mz2_mq_mq_cache[q][NumPar][i] );
    } else {
        complex newResult = PV.B0p(Mz(), Mz()*Mz(), mq(q), mq(q));
        CacheShift(B0p_Mz_Mz2_mq_mq_cache[q], NumPar, params, newResult);
        return newResult;
    }    
}

        
complex EWSMcache::B1_Mz_0_ml_mlprime(const int gen) const {
    int NumPar = 3;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mz(), mf, mfprime};

    int i = CacheCheck(B1_Mz_0_ml_mlprime_cache[gen], NumPar, params);
    if (i>=0) {
        return ( B1_Mz_0_ml_mlprime_cache[gen][NumPar][i] );
    } else {
        complex newResult = PV.B1(Mz(), 0.0, mf, mfprime);
        CacheShift(B1_Mz_0_ml_mlprime_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}   


complex EWSMcache::B1_Mz_0_mq_mqprime(const int gen) const {
    int NumPar = 3;
    double mf = mq((StandardModel::quark)(2*gen));    
    double mfprime = mq((StandardModel::quark)(2*gen+1));
    double params[] = {Mz(), mf, mfprime};

    int i = CacheCheck(B1_Mz_0_mq_mqprime_cache[gen], NumPar, params);
    if (i>=0) {
        return ( B1_Mz_0_mq_mqprime_cache[gen][NumPar][i] );
    } else {
        complex newResult = PV.B1(Mz(), 0.0, mf, mfprime);
        CacheShift(B1_Mz_0_mq_mqprime_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
} 


complex EWSMcache::B1_Mz_0_mlprime_ml(const int gen) const {
    int NumPar = 3;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mz(), mf, mfprime};

    int i = CacheCheck(B1_Mz_0_mlprime_ml_cache[gen], NumPar, params);
    if (i>=0) {
        return ( B1_Mz_0_mlprime_ml_cache[gen][NumPar][i] );
    } else {
        complex newResult = PV.B1(Mz(), 0.0, mfprime, mf);
        CacheShift(B1_Mz_0_mlprime_ml_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
} 


complex EWSMcache::B1_Mz_0_mqprime_mq(const int gen) const {
    int NumPar = 3;
    double mf = mq((StandardModel::quark)(2*gen));    
    double mfprime = mq((StandardModel::quark)(2*gen+1));
    double params[] = {Mz(), mf, mfprime};

    int i = CacheCheck(B1_Mz_0_mqprime_mq_cache[gen], NumPar, params);
    if (i>=0) {
        return ( B1_Mz_0_mqprime_mq_cache[gen][NumPar][i] );
    } else {
        complex newResult = PV.B1(Mz(), 0.0, mfprime, mf);
        CacheShift(B1_Mz_0_mqprime_mq_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::B1_Mz_Mw2_ml_mlprime(const int gen, const double Mw_i) const {
    int NumPar = 4;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mz(), Mw(Mw_i), mf, mfprime};

    int i = CacheCheck(B1_Mz_Mw2_ml_mlprime_cache[gen], NumPar, params);
    if (i>=0) {
        return ( B1_Mz_Mw2_ml_mlprime_cache[gen][NumPar][i] );
    } else {
        complex newResult = PV.B1(Mz(), Mw(Mw_i)*Mw(Mw_i), mf, mfprime);
        CacheShift(B1_Mz_Mw2_ml_mlprime_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}  


complex EWSMcache::B1_Mz_Mw2_mq_mqprime(const int gen, const double Mw_i) const {
    int NumPar = 4;
    double mf = mq((StandardModel::quark)(2*gen));    
    double mfprime = mq((StandardModel::quark)(2*gen+1));
    double params[] = {Mz(), Mw(Mw_i), mf, mfprime};

    int i = CacheCheck(B1_Mz_Mw2_mq_mqprime_cache[gen], NumPar, params);
    if (i>=0) {
        return ( B1_Mz_Mw2_mq_mqprime_cache[gen][NumPar][i] );
    } else {
        complex newResult = PV.B1(Mz(), Mw(Mw_i)*Mw(Mw_i), mf, mfprime);
        CacheShift(B1_Mz_Mw2_mq_mqprime_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
} 


complex EWSMcache::B1_Mz_Mw2_mlprime_ml(const int gen, const double Mw_i) const {
    int NumPar = 4;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mz(), Mw(Mw_i), mf, mfprime};

    int i = CacheCheck(B1_Mz_Mw2_mlprime_ml_cache[gen], NumPar, params);
    if (i>=0) {
        return ( B1_Mz_Mw2_mlprime_ml_cache[gen][NumPar][i] );
    } else {
        complex newResult = PV.B1(Mz(), Mw(Mw_i)*Mw(Mw_i), mfprime, mf);
        CacheShift(B1_Mz_Mw2_mlprime_ml_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}  


complex EWSMcache::B1_Mz_Mw2_mqprime_mq(const int gen, const double Mw_i) const {
    int NumPar = 4;
    double mf = mq((StandardModel::quark)(2*gen));    
    double mfprime = mq((StandardModel::quark)(2*gen+1));
    double params[] = {Mz(), Mw(Mw_i), mf, mfprime};

    int i = CacheCheck(B1_Mz_Mw2_mqprime_mq_cache[gen], NumPar, params);
    if (i>=0) {
        return ( B1_Mz_Mw2_mqprime_mq_cache[gen][NumPar][i] );
    } else {
        complex newResult = PV.B1(Mz(), Mw(Mw_i)*Mw(Mw_i), mfprime, mf);
        CacheShift(B1_Mz_Mw2_mqprime_mq_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
} 


complex EWSMcache::B1p_Mw_Mw2_ml_mlprime(const int gen, const double Mw_i) const {
    int NumPar = 3;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mw(Mw_i), mf, mfprime};

    int i = CacheCheck(B1p_Mw_Mw2_ml_mlprime_cache[gen], NumPar, params);
    if (i>=0) {
        return ( B1p_Mw_Mw2_ml_mlprime_cache[gen][NumPar][i] );
    } else {
        complex newResult = PV.B1p(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), mf, mfprime);
        CacheShift(B1p_Mw_Mw2_ml_mlprime_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::B1p_Mw_Mw2_mq_mqprime(const int gen, const double Mw_i) const {
    int NumPar = 3;
    double mf = mq((StandardModel::quark)(2*gen));    
    double mfprime = mq((StandardModel::quark)(2*gen+1));
    double params[] = {Mw(Mw_i), mf, mfprime};

    int i = CacheCheck(B1p_Mw_Mw2_mq_mqprime_cache[gen], NumPar, params);
    if (i>=0) {
        return ( B1p_Mw_Mw2_mq_mqprime_cache[gen][NumPar][i] );
    } else {
        complex newResult = PV.B1p(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), mf, mfprime);
        CacheShift(B1p_Mw_Mw2_mq_mqprime_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::B1p_Mw_Mw2_mlprime_ml(const int gen, const double Mw_i) const {
    int NumPar = 3;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mw(Mw_i), mf, mfprime};

    int i = CacheCheck(B1p_Mw_Mw2_mlprime_ml_cache[gen], NumPar, params);
    if (i>=0) {
        return ( B1p_Mw_Mw2_mlprime_ml_cache[gen][NumPar][i] );
    } else {
        complex newResult = PV.B1p(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), mfprime, mf);
        CacheShift(B1p_Mw_Mw2_mlprime_ml_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::B1p_Mw_Mw2_mqprime_mq(const int gen, const double Mw_i) const {
    int NumPar = 3;
    double mf = mq((StandardModel::quark)(2*gen));    
    double mfprime = mq((StandardModel::quark)(2*gen+1));
    double params[] = {Mw(Mw_i), mf, mfprime};

    int i = CacheCheck(B1p_Mw_Mw2_mqprime_mq_cache[gen], NumPar, params);
    if (i>=0) {
        return ( B1p_Mw_Mw2_mqprime_mq_cache[gen][NumPar][i] );
    } else {
        complex newResult = PV.B1p(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), mfprime, mf);
        CacheShift(B1p_Mw_Mw2_mqprime_mq_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::Bf_Mz_Mz2_ml_ml(const StandardModel::lepton l) const {
    int NumPar = 2;
    double params[] = {Mz(), ml(l)};

    int i = CacheCheck(Bf_Mz_Mz2_ml_ml_cache[l], NumPar, params);
    if (i>=0) {
        return ( Bf_Mz_Mz2_ml_ml_cache[l][NumPar][i] );
    } else {
        complex newResult = PV.Bf(Mz(), Mz()*Mz(), ml(l), ml(l));
        CacheShift(Bf_Mz_Mz2_ml_ml_cache[l], NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::Bf_Mz_Mz2_mq_mq(const StandardModel::quark q) const {
    int NumPar = 2;
    double params[] = {Mz(), mq(q)};

    int i = CacheCheck(Bf_Mz_Mz2_mq_mq_cache[q], NumPar, params);
    if (i>=0) {
        return ( Bf_Mz_Mz2_mq_mq_cache[q][NumPar][i] );
    } else {
        complex newResult = PV.Bf(Mz(), Mz()*Mz(), mq(q), mq(q));
        CacheShift(Bf_Mz_Mz2_mq_mq_cache[q], NumPar, params, newResult);
        return newResult;
    }    
}


complex EWSMcache::Bf_Mz_0_ml_ml(const StandardModel::lepton l) const {
    int NumPar = 2;
    double params[] = {Mz(), ml(l)};
    if (ml(l)==0.0)
        throw "Error in EWSMcache::Bf_Mz_0_ml_ml()";
    
    int i = CacheCheck(Bf_Mz_0_ml_ml_cache[l], NumPar, params);
    if (i>=0) {
        return ( Bf_Mz_0_ml_ml_cache[l][NumPar][i] );
    } else {
        complex newResult = PV.Bf(Mz(), 0.0, ml(l), ml(l));
            CacheShift(Bf_Mz_0_ml_ml_cache[l], NumPar, params, newResult);
            return newResult;
    }
}


complex EWSMcache::Bf_Mz_0_mq_mq(const StandardModel::quark q) const {
    int NumPar = 2;
    double params[] = {Mz(), mq(q)};
    if (mq(q)==0.0)
        throw "Error in EWSMcache::Bf_Mz_0_mq_mq()";
    
    int i = CacheCheck(Bf_Mz_0_mq_mq_cache[q], NumPar, params);
    if (i>=0) {
        return ( Bf_Mz_0_mq_mq_cache[q][NumPar][i] );
    } else {
        complex newResult = PV.Bf(Mz(), 0.0, mq(q), mq(q));
        CacheShift(Bf_Mz_0_mq_mq_cache[q], NumPar, params, newResult);
        return newResult;
    }    
}


complex EWSMcache::Bf_Mz_Mw2_mlprime_ml(const int gen, const double Mw_i) const {
    int NumPar = 4;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mz(), Mw(Mw_i), mf, mfprime};

    int i = CacheCheck(Bf_Mz_Mw2_mlprime_ml_cache[gen], NumPar, params);
    if (i>=0) {
        return ( Bf_Mz_Mw2_mlprime_ml_cache[gen][NumPar][i] );
    } else {
        complex newResult = PV.Bf(Mz(), Mw(Mw_i)*Mw(Mw_i), mfprime, mf);
        CacheShift(Bf_Mz_Mw2_mlprime_ml_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}        
        
        
complex EWSMcache::Bf_Mz_Mw2_mqprime_mq(const int gen, const double Mw_i) const {
    int NumPar = 4;
    double mf = mq((StandardModel::quark)(2*gen));    
    double mfprime = mq((StandardModel::quark)(2*gen+1));
    double params[] = {Mz(), Mw(Mw_i), mf, mfprime};

    int i = CacheCheck(Bf_Mz_Mw2_mqprime_mq_cache[gen], NumPar, params);
    if (i>=0) {
        return ( Bf_Mz_Mw2_mqprime_mq_cache[gen][NumPar][i] );
    } else {
        complex newResult = PV.Bf(Mz(), Mw(Mw_i)*Mw(Mw_i), mfprime, mf);
        CacheShift(Bf_Mz_Mw2_mqprime_mq_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::Bf_Mz_0_mlprime_ml(const int gen) const {
    int NumPar = 3;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mz(), mf, mfprime};

    int i = CacheCheck(Bf_Mz_0_mlprime_ml_cache[gen], NumPar, params);
    if (i>=0) {
        return ( Bf_Mz_0_mlprime_ml_cache[gen][NumPar][i] );
    } else {
        complex newResult = PV.Bf(Mz(), 0.0, mfprime, mf);
        CacheShift(Bf_Mz_0_mlprime_ml_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::Bf_Mz_0_mqprime_mq(const int gen) const {
    int NumPar = 3;
    double mf = mq((StandardModel::quark)(2*gen));    
    double mfprime = mq((StandardModel::quark)(2*gen+1));
    double params[] = {Mz(), mf, mfprime};

    int i = CacheCheck(Bf_Mz_0_mqprime_mq_cache[gen], NumPar, params);
    if (i>=0) {
        return ( Bf_Mz_0_mqprime_mq_cache[gen][NumPar][i] );
    } else {
        complex newResult = PV.Bf(Mz(), 0.0, mfprime, mf);
        CacheShift(Bf_Mz_0_mqprime_mq_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::Bf_Mw_Mw2_mlprime_ml(const int gen, const double Mw_i) const {
    int NumPar = 3;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mw(Mw_i), mf, mfprime};

    int i = CacheCheck(Bf_Mw_Mw2_mlprime_ml_cache[gen], NumPar, params);
    if (i>=0) {
        return ( Bf_Mw_Mw2_mlprime_ml_cache[gen][NumPar][i] );
    } else {
        complex newResult = PV.Bf(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), mfprime, mf);
        CacheShift(Bf_Mw_Mw2_mlprime_ml_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::Bf_Mw_Mw2_mqprime_mq(const int gen, const double Mw_i) const {
    int NumPar = 3;
    double mf = mq((StandardModel::quark)(2*gen));    
    double mfprime = mq((StandardModel::quark)(2*gen+1));
    double params[] = {Mw(Mw_i), mf, mfprime};

    int i = CacheCheck(Bf_Mw_Mw2_mqprime_mq_cache[gen], NumPar, params);
    if (i>=0) {
        return ( Bf_Mw_Mw2_mqprime_mq_cache[gen][NumPar][i] );
    } else {
        complex newResult = PV.Bf(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), mfprime, mf);
        CacheShift(Bf_Mw_Mw2_mqprime_mq_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::Bfp_Mz_Mz2_ml_ml(const StandardModel::lepton l) const {
    int NumPar = 2;
    double params[] = {Mz(), ml(l)};

    int i = CacheCheck(Bfp_Mz_Mz2_ml_ml_cache[l], NumPar, params);
    if (i>=0) {
        return ( Bfp_Mz_Mz2_ml_ml_cache[l][NumPar][i] );
    } else {
        complex newResult = PV.Bfp(Mz(), Mz()*Mz(), ml(l), ml(l));
        CacheShift(Bfp_Mz_Mz2_ml_ml_cache[l], NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::Bfp_Mz_Mz2_mq_mq(const StandardModel::quark q) const {
    int NumPar = 2;
    double params[] = {Mz(), mq(q)};

    int i = CacheCheck(Bfp_Mz_Mz2_mq_mq_cache[q], NumPar, params);
    if (i>=0) {
        return ( Bfp_Mz_Mz2_mq_mq_cache[q][NumPar][i] );
    } else {
        complex newResult = PV.Bfp(Mz(), Mz()*Mz(), mq(q), mq(q));
        CacheShift(Bfp_Mz_Mz2_mq_mq_cache[q], NumPar, params, newResult);
        return newResult;
    }    
}


complex EWSMcache::Bfp_Mw_Mw2_mlprime_ml(const int gen, const double Mw_i) const {
    int NumPar = 3;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mw(Mw_i), mf, mfprime};

    int i = CacheCheck(Bfp_Mw_Mw2_mlprime_ml_cache[gen], NumPar, params);
    if (i>=0) {
        return ( Bfp_Mw_Mw2_mlprime_ml_cache[gen][NumPar][i] );
    } else {
        complex newResult = PV.Bfp(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), mfprime, mf);
        CacheShift(Bfp_Mw_Mw2_mlprime_ml_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::Bfp_Mw_Mw2_mqprime_mq(const int gen, const double Mw_i) const {
    int NumPar = 3;
    double mf = mq((StandardModel::quark)(2*gen));    
    double mfprime = mq((StandardModel::quark)(2*gen+1));
    double params[] = {Mw(Mw_i), mf, mfprime};

    int i = CacheCheck(Bfp_Mw_Mw2_mqprime_mq_cache[gen], NumPar, params);
    if (i>=0) {
        return ( Bfp_Mw_Mw2_mqprime_mq_cache[gen][NumPar][i] );
    } else {
        complex newResult = PV.Bfp(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), mfprime, mf);
        CacheShift(Bfp_Mw_Mw2_mqprime_mq_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::C0_Mz2_Mw_Mt_Mw(const double Mw_i) const {
    int NumPar = 3;
    double params[] = {Mz(), Mw(Mw_i), Mt()};

    int i = CacheCheck(C0_Mz2_Mw_Mt_Mw_cache, NumPar, params);
    if (i>=0) {
        return ( C0_Mz2_Mw_Mt_Mw_cache[NumPar][i] );
    } else {
        complex newResult = PV.C0(Mz()*Mz(), Mw(Mw_i), Mt(), Mw(Mw_i));
        CacheShift(C0_Mz2_Mw_Mt_Mw_cache, NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::C0_Mz2_Mt_Mw_Mt(const double Mw_i) const {
    int NumPar = 3;
    double params[] = {Mz(), Mt(), Mw(Mw_i)};

    int i = CacheCheck(C0_Mz2_Mt_Mw_Mt_cache, NumPar, params);
    if (i>=0) {
        return ( C0_Mz2_Mt_Mw_Mt_cache[NumPar][i] );
    } else {
        complex newResult = PV.C0(Mz()*Mz(), Mt(), Mw(Mw_i), Mt());
        CacheShift(C0_Mz2_Mt_Mw_Mt_cache, NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::C0_Mz2_0_Mw_0(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    int i = CacheCheck(C0_Mz2_0_Mw_0_cache, NumPar, params);
    if (i>=0) {
        return ( C0_Mz2_0_Mw_0_cache[NumPar][i] );
    } else {
        complex newResult = PV.C0(Mz()*Mz(), 0.0, Mw(Mw_i), 0.0);
        CacheShift(C0_Mz2_0_Mw_0_cache, NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::C0_Mz2_Mw_0_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    int i = CacheCheck(C0_Mz2_Mw_0_Mw_cache, NumPar, params);
    if (i>=0) {
        return ( C0_Mz2_Mw_0_Mw_cache[NumPar][i] );
    } else {
        complex newResult = PV.C0(Mz()*Mz(), Mw(Mw_i), 0.0, Mw(Mw_i));
        CacheShift(C0_Mz2_Mw_0_Mw_cache, NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::C0_Mw2_Mw_0_Mz(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mw(Mw_i), Mz()};

    int i = CacheCheck(C0_Mw2_Mw_0_Mz_cache, NumPar, params);
    if (i>=0) {
        return ( C0_Mw2_Mw_0_Mz_cache[NumPar][i] );
    } else {
        complex newResult = PV.C0(Mw(Mw_i)*Mw(Mw_i), Mw(Mw_i), 0.0, Mz());
        CacheShift(C0_Mw2_Mw_0_Mz_cache, NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::C0_Mw2_0_Mz_0(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mw(Mw_i), Mz()};

    int i = CacheCheck(C0_Mw2_0_Mz_0_cache, NumPar, params);
    if (i>=0) {
        return ( C0_Mw2_0_Mz_0_cache[NumPar][i] );
    } else {
        complex newResult = PV.C0(Mw(Mw_i)*Mw(Mw_i), 0.0, Mz(), 0.0); 
        CacheShift(C0_Mw2_0_Mz_0_cache, NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::C0_Mz2_0_Mz_0() const {
    int NumPar = 1;
    double params[] = {Mz()};

    int i = CacheCheck(C0_Mz2_0_Mz_0_cache, NumPar, params);
    if (i>=0) {
        return ( C0_Mz2_0_Mz_0_cache[NumPar][i] );
    } else {
        complex newResult = PV.C0(Mz()*Mz(), 0.0, Mz(), 0.0);
        CacheShift(C0_Mz2_0_Mz_0_cache, NumPar, params, newResult);
        return newResult;
    } 
}




