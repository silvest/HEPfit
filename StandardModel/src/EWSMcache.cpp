/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <gsl/gsl_sf.h>
#include "EWSMcache.h"


EWSMcache::EWSMcache(const StandardModel& SM_i) : SM(SM_i) {
    bDebug = SM_i.isBDebug();
    
    bUseCacheEWSMcache = true;// use caches in the current class
    //bUseCacheEWSMcache = false;// do not use caches in the current class (for test)
    
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

double EWSMcache::ml(const StandardModel::lepton l) const {
    return SM.getLeptons(l).getMass();
}


double EWSMcache::mq(const StandardModel::quark q, const double mu, 
                     const orders order) const {
    switch(q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::STRANGE:
            if (bDebug)
                return SM.getQuarks(q).getMass(); // for debug
            else
                return SM.Mrun(mu, SM.getQuarks(q).getMass_scale(), 
                        SM.getQuarks(q).getMass(), order);                    
        case StandardModel::CHARM:
        case StandardModel::BOTTOM:
            if (bDebug)
                return SM.getQuarks(q).getMass(); // for debug
            else
                return SM.Mrun(mu, SM.getQuarks(q).getMass(), order);
        case StandardModel::TOP:
            return SM.getMtpole(); // the pole mass
        default:
            throw std::runtime_error("Error in EWSMcache::mq()"); 
    }
}


//////////////////////////////////////////////////////////////////////// 

double EWSMcache::logMZtoME() const {
    int NumPar = 2;
    double params[] = {Mz(), ml(SM.ELECTRON)};

    if ( CacheCheck(logMZtoME_cache, NumPar, params) )
        return logMZtoME_cache[NumPar];
    else {
        double newResult = log( Mz()/ml(SM.ELECTRON) );
        newCacheForDouble(logMZtoME_cache, NumPar, params, newResult);
        return newResult;
    }
}


double EWSMcache::logMZtoMMU() const {
    int NumPar = 2;
    double params[] = {Mz(), ml(SM.MU)};

    if ( CacheCheck(logMZtoMMU_cache, NumPar, params) )
        return logMZtoMMU_cache[NumPar];
    else {
        double newResult = log( Mz()/ml(SM.MU) );
        newCacheForDouble(logMZtoMMU_cache, NumPar, params, newResult);
        return newResult;
    }
}
    
    
double EWSMcache::logMZtoMTAU() const {
    int NumPar = 2;
    double params[] = {Mz(), ml(SM.TAU)};

    if ( CacheCheck(logMZtoMTAU_cache, NumPar, params) )
        return logMZtoMTAU_cache[NumPar];
    else {
        double newResult = log( Mz()/ml(SM.TAU) );
        newCacheForDouble(logMZtoMTAU_cache, NumPar, params, newResult);
        return newResult;
    }
}
    

double EWSMcache::logMZtoMTOP() const {
    int NumPar = 2;
    double params[] = {Mz(), Mt()};

    if ( CacheCheck(logMZtoMTOP_cache, NumPar, params) )
        return logMZtoMTOP_cache[NumPar];
    else {
        double newResult = log( Mz()/Mt() );
        newCacheForDouble(logMZtoMTOP_cache, NumPar, params, newResult);
        return newResult;
    }
}
    
    
double EWSMcache::logMTOPtoMH() const {
    int NumPar = 2;
    double params[] = {Mt(), mh()};

    if ( CacheCheck(logMTOPtoMH_cache, NumPar, params) )
        return logMTOPtoMH_cache[NumPar];
    else {
        double newResult = log( Mt()/mh() );
        newCacheForDouble(logMTOPtoMH_cache, NumPar, params, newResult);
        return newResult;
    }
}

    
double EWSMcache::log_cW2(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    if ( CacheCheck(log_cW2_cache, NumPar, params) )
        return log_cW2_cache[NumPar];
    else {
        double newResult = log(cW2(Mw_i));
        newCacheForDouble(log_cW2_cache, NumPar, params, newResult);
        return newResult;
    }
}    


double EWSMcache::Li2_MW2toMTOP2(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mw(Mw_i), Mt()};

    if ( CacheCheck(Li2_MW2toMTOP2_cache, NumPar, params) )
        return Li2_MW2toMTOP2_cache[NumPar];
    else {
        double newResult = PolyLog.Li2( Mw(Mw_i)*Mw(Mw_i)/Mt()/Mt() ).real();
        newCacheForDouble(Li2_MW2toMTOP2_cache, NumPar, params, newResult);
        return newResult;
    }
}


double EWSMcache::Li3_MW2toMTOP2(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mw(Mw_i), Mt()};

    if ( CacheCheck(Li3_MW2toMTOP2_cache, NumPar, params) )
        return Li3_MW2toMTOP2_cache[NumPar];
    else {
        double newResult = PolyLog.Li3(Mw(Mw_i)*Mw(Mw_i)/Mt()/Mt());
        newCacheForDouble(Li3_MW2toMTOP2_cache, NumPar, params, newResult);
        return newResult;
    }
}


double EWSMcache::Li3_for_F1(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mw(Mw_i), Mt()};

    if ( CacheCheck(Li3_for_F1_cache, NumPar, params) )
        return Li3_for_F1_cache[NumPar];
    else {
        double tmp = Mw(Mw_i)*Mw(Mw_i)/Mt()/Mt();
        double newResult = PolyLog.Li3(-tmp/(1.0 - tmp));
        newCacheForDouble(Li3_for_F1_cache, NumPar, params, newResult);
        return newResult;
    }
}


double EWSMcache::A0_Mz_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    if ( CacheCheck(A0_Mz_Mw_cache, NumPar, params) )
        return A0_Mz_Mw_cache[NumPar];
    else {
        double newResult = PV.A0(Mz(), Mw(Mw_i));
        newCacheForDouble(A0_Mz_Mw_cache, NumPar, params, newResult);
        return newResult;
    }    
}


double EWSMcache::A0_Mz_mh() const {
    int NumPar = 2;
    double params[] = {Mz(), mh()};

    if ( CacheCheck(A0_Mz_mh_cache, NumPar, params) )
        return A0_Mz_mh_cache[NumPar];
    else {
        double newResult = PV.A0(Mz(), mh());
        newCacheForDouble(A0_Mz_mh_cache, NumPar, params, newResult);
        return newResult;
    }      
}


double EWSMcache::A0_Mw_Mz(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mw(Mw_i), Mz()};

    if ( CacheCheck(A0_Mw_Mz_cache, NumPar, params) )
        return A0_Mw_Mz_cache[NumPar];
    else {
        double newResult = PV.A0(Mw(Mw_i), Mz());
        newCacheForDouble(A0_Mw_Mz_cache, NumPar, params, newResult);
        return newResult;
    }         
}


double EWSMcache::A0_Mw_mh(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mw(Mw_i), mh()};

    if ( CacheCheck(A0_Mw_mh_cache, NumPar, params) )
        return A0_Mw_mh_cache[NumPar];
    else {
        double newResult = PV.A0(Mw(Mw_i), mh());
        newCacheForDouble(A0_Mw_mh_cache, NumPar, params, newResult);
        return newResult;
    }     
}


double EWSMcache::A0_Mz_Mz() const {
    int NumPar = 1;
    double params[] = {Mz()};

    if ( CacheCheck(A0_Mz_Mz_cache, NumPar, params) )
        return A0_Mz_Mz_cache[NumPar];
    else {
        double newResult = PV.A0(Mz(), Mz());
        newCacheForDouble(A0_Mz_Mz_cache, NumPar, params, newResult);
        return newResult;
    }     
}


double EWSMcache::A0_Mw_Mw(const double Mw_i) const {
    int NumPar = 1;
    double params[] = {Mw(Mw_i)};

    if ( CacheCheck(A0_Mw_Mw_cache, NumPar, params) )
        return A0_Mw_Mw_cache[NumPar];
    else {
        double newResult = PV.A0(Mw(Mw_i), Mw(Mw_i));
        newCacheForDouble(A0_Mw_Mw_cache, NumPar, params, newResult);
        return newResult;
    }       
}


complex EWSMcache::B0_Mz_Mw2_mh_Mw(const double Mw_i) const {
    int NumPar = 3;
    double params[] = {Mz(), Mw(Mw_i), mh()};

    if ( CacheCheck(B0_Mz_Mw2_mh_Mw_cache, NumPar, params) )
        return complex(B0_Mz_Mw2_mh_Mw_cache[NumPar], 
                       B0_Mz_Mw2_mh_Mw_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0(Mz(), Mw(Mw_i)*Mw(Mw_i), mh(), Mw(Mw_i));
        newCacheForComplex(B0_Mz_Mw2_mh_Mw_cache, NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::B0_Mz_0_mh_Mw(const double Mw_i) const {
    int NumPar = 3;
    double params[] = {Mz(), mh(), Mw(Mw_i)};

    if ( CacheCheck(B0_Mz_0_mh_Mw_cache, NumPar, params) )
        return complex(B0_Mz_0_mh_Mw_cache[NumPar],
                       B0_Mz_0_mh_Mw_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0(Mz(), 0.0, mh(), Mw(Mw_i));
        newCacheForComplex(B0_Mz_0_mh_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mw_Mz2_Mt_Mt(const double Mw_i) const {
    int NumPar = 3;
    double params[] = {Mw(Mw_i), Mz(), Mt()};

    if ( CacheCheck(B0_Mw_Mz2_Mt_Mt_cache, NumPar, params) )
        return complex(B0_Mw_Mz2_Mt_Mt_cache[NumPar],
                       B0_Mw_Mz2_Mt_Mt_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0(Mw(Mw_i), Mz()*Mz(), Mt(), Mt());
        newCacheForComplex(B0_Mw_Mz2_Mt_Mt_cache, NumPar, params, newResult);
        return newResult;
    }    
}


complex EWSMcache::B0_Mz_Mz2_Mw_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    if ( CacheCheck(B0_Mz_Mz2_Mw_Mw_cache, NumPar, params) )
        return complex(B0_Mz_Mz2_Mw_Mw_cache[NumPar],
                       B0_Mz_Mz2_Mw_Mw_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0(Mz(), Mz()*Mz(), Mw(Mw_i), Mw(Mw_i));
        newCacheForComplex(B0_Mz_Mz2_Mw_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mz_Mz2_mh_Mz() const {
    int NumPar = 2;
    double params[] = {Mz(), mh()};

    if ( CacheCheck(B0_Mz_Mz2_mh_Mz_cache, NumPar, params) )
        return complex(B0_Mz_Mz2_mh_Mz_cache[NumPar],
                       B0_Mz_Mz2_mh_Mz_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0(Mz(), Mz()*Mz(), mh(), Mz());
        newCacheForComplex(B0_Mz_Mz2_mh_Mz_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mz_Mw2_Mz_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    if ( CacheCheck(B0_Mz_Mw2_Mz_Mw_cache, NumPar, params) )
        return complex(B0_Mz_Mw2_Mz_Mw_cache[NumPar],
                       B0_Mz_Mw2_Mz_Mw_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0(Mz(), Mw(Mw_i)*Mw(Mw_i), Mz(), Mw(Mw_i));
        newCacheForComplex(B0_Mz_Mw2_Mz_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mz_Mw2_0_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    if ( CacheCheck(B0_Mz_Mw2_0_Mw_cache, NumPar, params) )
        return complex(B0_Mz_Mw2_0_Mw_cache[NumPar], 
                       B0_Mz_Mw2_0_Mw_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0(Mz(), Mw(Mw_i)*Mw(Mw_i), 0.0, Mw(Mw_i));
        newCacheForComplex(B0_Mz_Mw2_0_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mz_0_Mz_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    if ( CacheCheck(B0_Mz_0_Mz_Mw_cache, NumPar, params) )
        return complex(B0_Mz_0_Mz_Mw_cache[NumPar],
                       B0_Mz_0_Mz_Mw_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0(Mz(), 0.0, Mz(), Mw(Mw_i));
        newCacheForComplex(B0_Mz_0_Mz_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mz_0_0_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    if ( CacheCheck(B0_Mz_0_0_Mw_cache, NumPar, params) )
        return complex(B0_Mz_0_0_Mw_cache[NumPar], 
                       B0_Mz_0_0_Mw_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0(Mz(), 0.0, 0.0, Mw(Mw_i));
        newCacheForComplex(B0_Mz_0_0_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mw_Mz2_Mw_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mw(Mw_i), Mz()};

    if ( CacheCheck(B0_Mw_Mz2_Mw_Mw_cache, NumPar, params) )
        return complex(B0_Mw_Mz2_Mw_Mw_cache[NumPar], 
                       B0_Mw_Mz2_Mw_Mw_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0(Mw(Mw_i), Mz()*Mz(), Mw(Mw_i), Mw(Mw_i));
        newCacheForComplex(B0_Mw_Mz2_Mw_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mw_Mw2_Mz_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mw(Mw_i), Mz()};

    if ( CacheCheck(B0_Mw_Mw2_Mz_Mw_cache, NumPar, params) )
        return complex(B0_Mw_Mw2_Mz_Mw_cache[NumPar], 
                       B0_Mw_Mw2_Mz_Mw_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), Mz(), Mw(Mw_i));
        newCacheForComplex(B0_Mw_Mw2_Mz_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mw_Mw2_mh_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mw(Mw_i), mh()};

    if ( CacheCheck(B0_Mw_Mw2_mh_Mw_cache, NumPar, params) )
        return complex(B0_Mw_Mw2_mh_Mw_cache[NumPar], 
                       B0_Mw_Mw2_mh_Mw_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), mh(), Mw(Mw_i));
        newCacheForComplex(B0_Mw_Mw2_mh_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mw_Mw2_0_Mw(const double Mw_i) const {
    int NumPar = 1;
    double params[] = {Mw(Mw_i)};
    
    if ( CacheCheck(B0_Mw_Mw2_0_Mw_cache, NumPar, params) )
        return complex(B0_Mw_Mw2_0_Mw_cache[NumPar], 
                       B0_Mw_Mw2_0_Mw_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), 0.0, Mw(Mw_i));
        newCacheForComplex(B0_Mw_Mw2_0_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mz_Mz2_ml_ml(const StandardModel::lepton l) const {
    int NumPar = 2;
    double params[] = {Mz(), ml(l)};

    if ( CacheCheck(B0_Mz_Mz2_ml_ml_cache[l], NumPar, params) )
        return complex(B0_Mz_Mz2_ml_ml_cache[l][NumPar], 
                       B0_Mz_Mz2_ml_ml_cache[l][NumPar+1], false);
    else {
        complex newResult = PV.B0(Mz(), Mz()*Mz(), ml(l), ml(l));
        newCacheForComplex(B0_Mz_Mz2_ml_ml_cache[l], NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mz_Mz2_mq_mq(const StandardModel::quark q) const {
    int NumPar = 2;
    double params[] = {Mz(), mq(q, Mz())};

    if ( CacheCheck(B0_Mz_Mz2_mq_mq_cache[q], NumPar, params) )
        return complex(B0_Mz_Mz2_mq_mq_cache[q][NumPar], 
                       B0_Mz_Mz2_mq_mq_cache[q][NumPar+1], false);
    else {
        complex newResult = PV.B0(Mz(), Mz()*Mz(), mq(q, Mz()), mq(q, Mz()));
        newCacheForComplex(B0_Mz_Mz2_mq_mq_cache[q], NumPar, params, newResult);
        return newResult;
    }    
}


complex EWSMcache::B0p_Mz_0_mh_Mw(const double Mw_i) const {
    int NumPar = 3;
    double params[] = {Mz(), mh(), Mw(Mw_i)};

    if ( CacheCheck(B0p_Mz_0_mh_Mw_cache, NumPar, params) )
        return complex(B0p_Mz_0_mh_Mw_cache[NumPar], 
                       B0p_Mz_0_mh_Mw_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0p(Mz(), 0.0, mh(), Mw(Mw_i));
        newCacheForComplex(B0p_Mz_0_mh_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0p_Mz_Mz2_mh_Mz() const {
    int NumPar = 2;
    double params[] = {Mz(), mh()};

    if ( CacheCheck(B0p_Mz_Mz2_mh_Mz_cache, NumPar, params) )
        return complex(B0p_Mz_Mz2_mh_Mz_cache[NumPar], 
                       B0p_Mz_Mz2_mh_Mz_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0p(Mz(), Mz()*Mz(), mh(), Mz());
        newCacheForComplex(B0p_Mz_Mz2_mh_Mz_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0p_Mz_0_Mz_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    if ( CacheCheck(B0p_Mz_0_Mz_Mw_cache, NumPar, params) )
        return complex(B0p_Mz_0_Mz_Mw_cache[NumPar], 
                       B0p_Mz_0_Mz_Mw_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0p(Mz(), 0.0, Mz(), Mw(Mw_i));
        newCacheForComplex(B0p_Mz_0_Mz_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0p_Mz_Mz2_Mw_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    if ( CacheCheck(B0p_Mz_Mz2_Mw_Mw_cache, NumPar, params) )
        return complex(B0p_Mz_Mz2_Mw_Mw_cache[NumPar], 
                       B0p_Mz_Mz2_Mw_Mw_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0p(Mz(), Mz()*Mz(), Mw(Mw_i), Mw(Mw_i));
        newCacheForComplex(B0p_Mz_Mz2_Mw_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0p_Mw_Mw2_Mz_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mw(Mw_i), Mz()};

    if ( CacheCheck(B0p_Mw_Mw2_Mz_Mw_cache, NumPar, params) )
        return complex(B0p_Mw_Mw2_Mz_Mw_cache[NumPar], 
                       B0p_Mw_Mw2_Mz_Mw_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0p(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), Mz(), Mw(Mw_i));
        newCacheForComplex(B0p_Mw_Mw2_Mz_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0p_Mw_Mw2_mh_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mw(Mw_i), mh()};

    if ( CacheCheck(B0p_Mw_Mw2_mh_Mw_cache, NumPar, params) )
        return complex(B0p_Mw_Mw2_mh_Mw_cache[NumPar], 
                       B0p_Mw_Mw2_mh_Mw_cache[NumPar+1], false); 
    else {
        complex newResult = PV.B0p(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), mh(), Mw(Mw_i));
        newCacheForComplex(B0p_Mw_Mw2_mh_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0p_Mw_Mw2_0_Mw(const double Mw_i) const {
    int NumPar = 1;
    double params[] = {Mw(Mw_i)};

    if ( CacheCheck(B0p_Mw_Mw2_0_Mw_cache, NumPar, params) )
        return complex(B0p_Mw_Mw2_0_Mw_cache[NumPar], 
                       B0p_Mw_Mw2_0_Mw_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0p(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), 0.0, Mw(Mw_i));
        newCacheForComplex(B0p_Mw_Mw2_0_Mw_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0p_Mz_Mz2_ml_ml(const StandardModel::lepton l) const {
    int NumPar = 2;
    double params[] = {Mz(), ml(l)};

    if ( CacheCheck(B0p_Mz_Mz2_ml_ml_cache[l], NumPar, params) )
        return complex(B0p_Mz_Mz2_ml_ml_cache[l][NumPar], 
                       B0p_Mz_Mz2_ml_ml_cache[l][NumPar+1], false);
    else {
        complex newResult = PV.B0p(Mz(), Mz()*Mz(), ml(l), ml(l));
        newCacheForComplex(B0p_Mz_Mz2_ml_ml_cache[l], NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0p_Mz_Mz2_mq_mq(const StandardModel::quark q) const {
    int NumPar = 2;
    double params[] = {Mz(), mq(q, Mz())};

    if ( CacheCheck(B0p_Mz_Mz2_mq_mq_cache[q], NumPar, params) )
        return complex(B0p_Mz_Mz2_mq_mq_cache[q][NumPar], 
                       B0p_Mz_Mz2_mq_mq_cache[q][NumPar+1], false);
    else {
        complex newResult = PV.B0p(Mz(), Mz()*Mz(), mq(q, Mz()), mq(q, Mz()));
        newCacheForComplex(B0p_Mz_Mz2_mq_mq_cache[q], NumPar, params, newResult);
        return newResult;
    }    
}

        
complex EWSMcache::B1_Mz_0_ml_mlprime(const int gen) const {
    int NumPar = 3;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mz(), mf, mfprime};

    if ( CacheCheck(B1_Mz_0_ml_mlprime_cache[gen], NumPar, params) )
        return complex(B1_Mz_0_ml_mlprime_cache[gen][NumPar], 
                       B1_Mz_0_ml_mlprime_cache[gen][NumPar+1], false);
    else {
        complex newResult = PV.B1(Mz(), 0.0, mf, mfprime);
        newCacheForComplex(B1_Mz_0_ml_mlprime_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}   


complex EWSMcache::B1_Mz_0_mq_mqprime(const int gen) const {
    int NumPar = 3;
    double mf = mq((StandardModel::quark)(2*gen), Mz());    
    double mfprime = mq((StandardModel::quark)(2*gen+1), Mz());
    double params[] = {Mz(), mf, mfprime};

    if ( CacheCheck(B1_Mz_0_mq_mqprime_cache[gen], NumPar, params) )
        return complex(B1_Mz_0_mq_mqprime_cache[gen][NumPar], 
                       B1_Mz_0_mq_mqprime_cache[gen][NumPar+1], false);
    else {
        complex newResult = PV.B1(Mz(), 0.0, mf, mfprime);
        newCacheForComplex(B1_Mz_0_mq_mqprime_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
} 


complex EWSMcache::B1_Mz_0_mlprime_ml(const int gen) const {
    int NumPar = 3;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mz(), mf, mfprime};

    if ( CacheCheck(B1_Mz_0_mlprime_ml_cache[gen], NumPar, params) )
        return complex(B1_Mz_0_mlprime_ml_cache[gen][NumPar], 
                       B1_Mz_0_mlprime_ml_cache[gen][NumPar+1], false);
    else {
        complex newResult = PV.B1(Mz(), 0.0, mfprime, mf);
        newCacheForComplex(B1_Mz_0_mlprime_ml_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
} 


complex EWSMcache::B1_Mz_0_mqprime_mq(const int gen) const {
    int NumPar = 3;
    double mf = mq((StandardModel::quark)(2*gen), Mz());    
    double mfprime = mq((StandardModel::quark)(2*gen+1), Mz());
    double params[] = {Mz(), mf, mfprime};

    if ( CacheCheck(B1_Mz_0_mqprime_mq_cache[gen], NumPar, params) )
        return complex(B1_Mz_0_mqprime_mq_cache[gen][NumPar], 
                       B1_Mz_0_mqprime_mq_cache[gen][NumPar+1], false);
    else {
        complex newResult = PV.B1(Mz(), 0.0, mfprime, mf);
        newCacheForComplex(B1_Mz_0_mqprime_mq_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::B1_Mz_Mw2_ml_mlprime(const int gen, const double Mw_i) const {
    int NumPar = 4;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mz(), Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(B1_Mz_Mw2_ml_mlprime_cache[gen], NumPar, params) )
        return complex(B1_Mz_Mw2_ml_mlprime_cache[gen][NumPar], 
                       B1_Mz_Mw2_ml_mlprime_cache[gen][NumPar+1], false);
    else {
        complex newResult = PV.B1(Mz(), Mw(Mw_i)*Mw(Mw_i), mf, mfprime);
        newCacheForComplex(B1_Mz_Mw2_ml_mlprime_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}  


complex EWSMcache::B1_Mz_Mw2_mq_mqprime(const int gen, const double Mw_i) const {
    int NumPar = 4;
    double mf = mq((StandardModel::quark)(2*gen), Mz());    
    double mfprime = mq((StandardModel::quark)(2*gen+1), Mz());
    double params[] = {Mz(), Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(B1_Mz_Mw2_mq_mqprime_cache[gen], NumPar, params) )
        return complex(B1_Mz_Mw2_mq_mqprime_cache[gen][NumPar], 
                       B1_Mz_Mw2_mq_mqprime_cache[gen][NumPar+1], false);
    else {
        complex newResult = PV.B1(Mz(), Mw(Mw_i)*Mw(Mw_i), mf, mfprime);
        newCacheForComplex(B1_Mz_Mw2_mq_mqprime_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
} 


complex EWSMcache::B1_Mz_Mw2_mlprime_ml(const int gen, const double Mw_i) const {
    int NumPar = 4;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mz(), Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(B1_Mz_Mw2_mlprime_ml_cache[gen], NumPar, params) )
        return complex(B1_Mz_Mw2_mlprime_ml_cache[gen][NumPar], 
                       B1_Mz_Mw2_mlprime_ml_cache[gen][NumPar+1], false);
    else {
        complex newResult = PV.B1(Mz(), Mw(Mw_i)*Mw(Mw_i), mfprime, mf);
        newCacheForComplex(B1_Mz_Mw2_mlprime_ml_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}  


complex EWSMcache::B1_Mz_Mw2_mqprime_mq(const int gen, const double Mw_i) const {
    int NumPar = 4;
    double mf = mq((StandardModel::quark)(2*gen), Mz());    
    double mfprime = mq((StandardModel::quark)(2*gen+1), Mz());
    double params[] = {Mz(), Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(B1_Mz_Mw2_mqprime_mq_cache[gen], NumPar, params) )
        return complex(B1_Mz_Mw2_mqprime_mq_cache[gen][NumPar], 
                       B1_Mz_Mw2_mqprime_mq_cache[gen][NumPar+1], false);
    else {
        complex newResult = PV.B1(Mz(), Mw(Mw_i)*Mw(Mw_i), mfprime, mf);
        newCacheForComplex(B1_Mz_Mw2_mqprime_mq_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
} 


complex EWSMcache::B1p_Mw_Mw2_ml_mlprime(const int gen, const double Mw_i) const {
    int NumPar = 3;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(B1p_Mw_Mw2_ml_mlprime_cache[gen], NumPar, params) )
        return complex(B1p_Mw_Mw2_ml_mlprime_cache[gen][NumPar], 
                       B1p_Mw_Mw2_ml_mlprime_cache[gen][NumPar+1], false);
    else {
        complex newResult = PV.B1p(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), mf, mfprime);
        newCacheForComplex(B1p_Mw_Mw2_ml_mlprime_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::B1p_Mw_Mw2_mq_mqprime(const int gen, const double Mw_i) const {
    int NumPar = 3;
    double mf = mq((StandardModel::quark)(2*gen), Mz());    
    double mfprime = mq((StandardModel::quark)(2*gen+1), Mz());
    double params[] = {Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(B1p_Mw_Mw2_mq_mqprime_cache[gen], NumPar, params) )
        return complex(B1p_Mw_Mw2_mq_mqprime_cache[gen][NumPar], 
                       B1p_Mw_Mw2_mq_mqprime_cache[gen][NumPar+1], false);
    else {
        complex newResult = PV.B1p(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), mf, mfprime);
        newCacheForComplex(B1p_Mw_Mw2_mq_mqprime_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::B1p_Mw_Mw2_mlprime_ml(const int gen, const double Mw_i) const {
    int NumPar = 3;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(B1p_Mw_Mw2_mlprime_ml_cache[gen], NumPar, params) )
        return complex(B1p_Mw_Mw2_mlprime_ml_cache[gen][NumPar], 
                       B1p_Mw_Mw2_mlprime_ml_cache[gen][NumPar+1], false);
    else {
        complex newResult = PV.B1p(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), mfprime, mf);
        newCacheForComplex(B1p_Mw_Mw2_mlprime_ml_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::B1p_Mw_Mw2_mqprime_mq(const int gen, const double Mw_i) const {
    int NumPar = 3;
    double mf = mq((StandardModel::quark)(2*gen), Mz());    
    double mfprime = mq((StandardModel::quark)(2*gen+1), Mz());
    double params[] = {Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(B1p_Mw_Mw2_mqprime_mq_cache[gen], NumPar, params) )
        return complex(B1p_Mw_Mw2_mqprime_mq_cache[gen][NumPar], 
                       B1p_Mw_Mw2_mqprime_mq_cache[gen][NumPar+1], false);
    else {
        complex newResult = PV.B1p(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), mfprime, mf);
        newCacheForComplex(B1p_Mw_Mw2_mqprime_mq_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::Bf_Mz_Mz2_ml_ml(const StandardModel::lepton l) const {
    int NumPar = 2;
    double params[] = {Mz(), ml(l)};

    if ( CacheCheck(Bf_Mz_Mz2_ml_ml_cache[l], NumPar, params) )
        return complex(Bf_Mz_Mz2_ml_ml_cache[l][NumPar], 
                       Bf_Mz_Mz2_ml_ml_cache[l][NumPar+1], false);
    else {
        complex newResult = PV.Bf(Mz(), Mz()*Mz(), ml(l), ml(l));
        newCacheForComplex(Bf_Mz_Mz2_ml_ml_cache[l], NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::Bf_Mz_Mz2_mq_mq(const StandardModel::quark q) const {
    int NumPar = 2;
    double params[] = {Mz(), mq(q, Mz())};

    if ( CacheCheck(Bf_Mz_Mz2_mq_mq_cache[q], NumPar, params) )
        return complex(Bf_Mz_Mz2_mq_mq_cache[q][NumPar], 
                       Bf_Mz_Mz2_mq_mq_cache[q][NumPar+1], false);
    else {
        complex newResult = PV.Bf(Mz(), Mz()*Mz(), mq(q, Mz()), mq(q, Mz()));
        newCacheForComplex(Bf_Mz_Mz2_mq_mq_cache[q], NumPar, params, newResult);
        return newResult;
    }    
}


complex EWSMcache::Bf_Mz_0_ml_ml(const StandardModel::lepton l) const {
    int NumPar = 2;
    double params[] = {Mz(), ml(l)};
    if (ml(l)==0.0)
        throw std::runtime_error("Error in EWSMcache::Bf_Mz_0_ml_ml()"); 
    
    if ( CacheCheck(Bf_Mz_0_ml_ml_cache[l], NumPar, params) )
        return complex(Bf_Mz_0_ml_ml_cache[l][NumPar], 
                       Bf_Mz_0_ml_ml_cache[l][NumPar+1], false);
    else {
        complex newResult = PV.Bf(Mz(), 0.0, ml(l), ml(l));
            newCacheForComplex(Bf_Mz_0_ml_ml_cache[l], NumPar, params, newResult);
            return newResult;
    }
}


complex EWSMcache::Bf_Mz_0_mq_mq(const StandardModel::quark q) const {
    int NumPar = 2;
    double params[] = {Mz(), mq(q, Mz())};
    if (mq(q, Mz())==0.0)
        throw std::runtime_error("Error in EWSMcache::Bf_Mz_0_mq_mq()"); 
    
    if ( CacheCheck(Bf_Mz_0_mq_mq_cache[q], NumPar, params) )
        return complex(Bf_Mz_0_mq_mq_cache[q][NumPar], 
                       Bf_Mz_0_mq_mq_cache[q][NumPar+1], false);
    else {
        complex newResult = PV.Bf(Mz(), 0.0, mq(q, Mz()), mq(q, Mz()));
        newCacheForComplex(Bf_Mz_0_mq_mq_cache[q], NumPar, params, newResult);
        return newResult;
    }    
}


complex EWSMcache::Bf_Mz_Mw2_mlprime_ml(const int gen, const double Mw_i) const {
    int NumPar = 4;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mz(), Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(Bf_Mz_Mw2_mlprime_ml_cache[gen], NumPar, params) )
        return complex(Bf_Mz_Mw2_mlprime_ml_cache[gen][NumPar], 
                       Bf_Mz_Mw2_mlprime_ml_cache[gen][NumPar+1], false);
    else {
        complex newResult = PV.Bf(Mz(), Mw(Mw_i)*Mw(Mw_i), mfprime, mf);
        newCacheForComplex(Bf_Mz_Mw2_mlprime_ml_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}        
        
        
complex EWSMcache::Bf_Mz_Mw2_mqprime_mq(const int gen, const double Mw_i) const {
    int NumPar = 4;
    double mf = mq((StandardModel::quark)(2*gen), Mz());    
    double mfprime = mq((StandardModel::quark)(2*gen+1), Mz());
    double params[] = {Mz(), Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(Bf_Mz_Mw2_mqprime_mq_cache[gen], NumPar, params) )
        return complex(Bf_Mz_Mw2_mqprime_mq_cache[gen][NumPar], 
                       Bf_Mz_Mw2_mqprime_mq_cache[gen][NumPar+1], false);
    else {
        complex newResult = PV.Bf(Mz(), Mw(Mw_i)*Mw(Mw_i), mfprime, mf);
        newCacheForComplex(Bf_Mz_Mw2_mqprime_mq_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::Bf_Mz_0_mlprime_ml(const int gen) const {
    int NumPar = 3;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mz(), mf, mfprime};

    if ( CacheCheck(Bf_Mz_0_mlprime_ml_cache[gen], NumPar, params) )
        return complex(Bf_Mz_0_mlprime_ml_cache[gen][NumPar], 
                       Bf_Mz_0_mlprime_ml_cache[gen][NumPar+1], false);
    else {
        complex newResult = PV.Bf(Mz(), 0.0, mfprime, mf);
        newCacheForComplex(Bf_Mz_0_mlprime_ml_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::Bf_Mz_0_mqprime_mq(const int gen) const {
    int NumPar = 3;
    double mf = mq((StandardModel::quark)(2*gen), Mz());    
    double mfprime = mq((StandardModel::quark)(2*gen+1), Mz());
    double params[] = {Mz(), mf, mfprime};

    if ( CacheCheck(Bf_Mz_0_mqprime_mq_cache[gen], NumPar, params) )
        return complex(Bf_Mz_0_mqprime_mq_cache[gen][NumPar], 
                       Bf_Mz_0_mqprime_mq_cache[gen][NumPar+1], false);
    else {
        complex newResult = PV.Bf(Mz(), 0.0, mfprime, mf);
        newCacheForComplex(Bf_Mz_0_mqprime_mq_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::Bf_Mw_Mw2_mlprime_ml(const int gen, const double Mw_i) const {
    int NumPar = 3;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(Bf_Mw_Mw2_mlprime_ml_cache[gen], NumPar, params) )
        return complex(Bf_Mw_Mw2_mlprime_ml_cache[gen][NumPar], 
                       Bf_Mw_Mw2_mlprime_ml_cache[gen][NumPar+1], false);
    else {
        complex newResult = PV.Bf(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), mfprime, mf);
        newCacheForComplex(Bf_Mw_Mw2_mlprime_ml_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::Bf_Mw_Mw2_mqprime_mq(const int gen, const double Mw_i) const {
    int NumPar = 3;
    double mf = mq((StandardModel::quark)(2*gen), Mz());    
    double mfprime = mq((StandardModel::quark)(2*gen+1), Mz());
    double params[] = {Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(Bf_Mw_Mw2_mqprime_mq_cache[gen], NumPar, params) )
        return complex(Bf_Mw_Mw2_mqprime_mq_cache[gen][NumPar], 
                       Bf_Mw_Mw2_mqprime_mq_cache[gen][NumPar+1], false);
    else {
        complex newResult = PV.Bf(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), mfprime, mf);
        newCacheForComplex(Bf_Mw_Mw2_mqprime_mq_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::Bfp_Mz_Mz2_ml_ml(const StandardModel::lepton l) const {
    int NumPar = 2;
    double params[] = {Mz(), ml(l)};

    if ( CacheCheck(Bfp_Mz_Mz2_ml_ml_cache[l], NumPar, params) )
        return complex(Bfp_Mz_Mz2_ml_ml_cache[l][NumPar], 
                       Bfp_Mz_Mz2_ml_ml_cache[l][NumPar+1], false);
    else {
        complex newResult = PV.Bfp(Mz(), Mz()*Mz(), ml(l), ml(l));
        newCacheForComplex(Bfp_Mz_Mz2_ml_ml_cache[l], NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::Bfp_Mz_Mz2_mq_mq(const StandardModel::quark q) const {
    int NumPar = 2;
    double params[] = {Mz(), mq(q, Mz())};

    if ( CacheCheck(Bfp_Mz_Mz2_mq_mq_cache[q], NumPar, params) )
        return complex(Bfp_Mz_Mz2_mq_mq_cache[q][NumPar], 
                       Bfp_Mz_Mz2_mq_mq_cache[q][NumPar+1], false);
    else {
        complex newResult = PV.Bfp(Mz(), Mz()*Mz(), mq(q, Mz()), mq(q, Mz()));
        newCacheForComplex(Bfp_Mz_Mz2_mq_mq_cache[q], NumPar, params, newResult);
        return newResult;
    }    
}


complex EWSMcache::Bfp_Mw_Mw2_mlprime_ml(const int gen, const double Mw_i) const {
    int NumPar = 3;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(Bfp_Mw_Mw2_mlprime_ml_cache[gen], NumPar, params) )
        return complex(Bfp_Mw_Mw2_mlprime_ml_cache[gen][NumPar], 
                       Bfp_Mw_Mw2_mlprime_ml_cache[gen][NumPar+1], false);
    else {
        complex newResult = PV.Bfp(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), mfprime, mf);
        newCacheForComplex(Bfp_Mw_Mw2_mlprime_ml_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::Bfp_Mw_Mw2_mqprime_mq(const int gen, const double Mw_i) const {
    int NumPar = 3;
    double mf = mq((StandardModel::quark)(2*gen), Mz());    
    double mfprime = mq((StandardModel::quark)(2*gen+1), Mz());
    double params[] = {Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(Bfp_Mw_Mw2_mqprime_mq_cache[gen], NumPar, params) )
        return complex(Bfp_Mw_Mw2_mqprime_mq_cache[gen][NumPar], 
                       Bfp_Mw_Mw2_mqprime_mq_cache[gen][NumPar+1], false);
    else {
        complex newResult = PV.Bfp(Mw(Mw_i), Mw(Mw_i)*Mw(Mw_i), mfprime, mf);
        newCacheForComplex(Bfp_Mw_Mw2_mqprime_mq_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::C0_Mz2_Mw_Mt_Mw(const double Mw_i) const {
    int NumPar = 3;
    double params[] = {Mz(), Mw(Mw_i), Mt()};

    if ( CacheCheck(C0_Mz2_Mw_Mt_Mw_cache, NumPar, params) )
        return complex(C0_Mz2_Mw_Mt_Mw_cache[NumPar], 
                       C0_Mz2_Mw_Mt_Mw_cache[NumPar+1], false);
    else {
        complex newResult = PV.C0(Mz()*Mz(), Mw(Mw_i), Mt(), Mw(Mw_i));
        newCacheForComplex(C0_Mz2_Mw_Mt_Mw_cache, NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::C0_Mz2_Mt_Mw_Mt(const double Mw_i) const {
    int NumPar = 3;
    double params[] = {Mz(), Mt(), Mw(Mw_i)};

    if ( CacheCheck(C0_Mz2_Mt_Mw_Mt_cache, NumPar, params) )
        return complex(C0_Mz2_Mt_Mw_Mt_cache[NumPar], 
                       C0_Mz2_Mt_Mw_Mt_cache[NumPar+1], false);
    else {
        complex newResult = PV.C0(Mz()*Mz(), Mt(), Mw(Mw_i), Mt());
        newCacheForComplex(C0_Mz2_Mt_Mw_Mt_cache, NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::C0_Mz2_0_Mw_0(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    if ( CacheCheck(C0_Mz2_0_Mw_0_cache, NumPar, params) )
        return complex(C0_Mz2_0_Mw_0_cache[NumPar], 
                       C0_Mz2_0_Mw_0_cache[NumPar+1], false);
    else {
        complex newResult = PV.C0(Mz()*Mz(), 0.0, Mw(Mw_i), 0.0);
        newCacheForComplex(C0_Mz2_0_Mw_0_cache, NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::C0_Mz2_Mw_0_Mw(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    if ( CacheCheck(C0_Mz2_Mw_0_Mw_cache, NumPar, params) )
        return complex(C0_Mz2_Mw_0_Mw_cache[NumPar], 
                       C0_Mz2_Mw_0_Mw_cache[NumPar+1], false);
    else {
        complex newResult = PV.C0(Mz()*Mz(), Mw(Mw_i), 0.0, Mw(Mw_i));
        newCacheForComplex(C0_Mz2_Mw_0_Mw_cache, NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::C0_Mw2_Mw_0_Mz(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mw(Mw_i), Mz()};

    if ( CacheCheck(C0_Mw2_Mw_0_Mz_cache, NumPar, params) )
        return complex(C0_Mw2_Mw_0_Mz_cache[NumPar], 
                       C0_Mw2_Mw_0_Mz_cache[NumPar+1], false);
    else {
        complex newResult = PV.C0(Mw(Mw_i)*Mw(Mw_i), Mw(Mw_i), 0.0, Mz());
        newCacheForComplex(C0_Mw2_Mw_0_Mz_cache, NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::C0_Mw2_0_Mz_0(const double Mw_i) const {
    int NumPar = 2;
    double params[] = {Mw(Mw_i), Mz()};

    if ( CacheCheck(C0_Mw2_0_Mz_0_cache, NumPar, params) )
        return complex(C0_Mw2_0_Mz_0_cache[NumPar], 
                       C0_Mw2_0_Mz_0_cache[NumPar+1], false);
    else {
        complex newResult = PV.C0(Mw(Mw_i)*Mw(Mw_i), 0.0, Mz(), 0.0); 
        newCacheForComplex(C0_Mw2_0_Mz_0_cache, NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::C0_Mz2_0_Mz_0() const {
    int NumPar = 1;
    double params[] = {Mz()};

    if ( CacheCheck(C0_Mz2_0_Mz_0_cache, NumPar, params) )
        return complex(C0_Mz2_0_Mz_0_cache[NumPar],
                       C0_Mz2_0_Mz_0_cache[NumPar+1], false);
    else {
        complex newResult = PV.C0(Mz()*Mz(), 0.0, Mz(), 0.0);
        newCacheForComplex(C0_Mz2_0_Mz_0_cache, NumPar, params, newResult);
        return newResult;
    } 
}




