/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
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


EWSMcache::EWSMcache(const StandardModel& SM_i) 
: SM(SM_i), PV(true)
{
    FlagDebug = false;
    FlagCacheInEWSMcache = true;// use caches in the current class
    //FlagCacheInEWSMcache = false;// do not use caches in the current class (for test)
    
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

double EWSMcache::ml(const StandardModel::lepton l) const 
{
    return SM.getLeptons(l).getMass();
}


double EWSMcache::mq(const QCD::quark q, const double mu, 
                     const orders order) const 
{
    switch(q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::STRANGE:
            if (FlagDebug)
                return SM.getQuarks(q).getMass();// for debug
            else
                return SM.Mrun(mu, SM.getQuarks(q).getMass_scale(),
                               SM.getQuarks(q).getMass(), order);

        case StandardModel::CHARM:
        case StandardModel::BOTTOM:
            if (FlagDebug)
                return SM.getQuarks(q).getMass();// for debug
            else
                return SM.Mrun(mu, SM.getQuarks(q).getMass(), order);
        case StandardModel::TOP:
            return SM.getMtpole(); // the pole mass
        default:
            throw std::runtime_error("Error in EWSMcache::mq()"); 
    }
}


//////////////////////////////////////////////////////////////////////// 

double EWSMcache::logMZtoME() const 
{
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


double EWSMcache::logMZtoMMU() const 
{
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
    
    
double EWSMcache::logMZtoMTAU() const 
{
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
    

double EWSMcache::logMZtoMTOP() const 
{
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
    
    
double EWSMcache::logMTOPtoMH() const
{
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

    
double EWSMcache::log_cW2(const double Mw_i) const 
{
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


double EWSMcache::Li2_MW2toMTOP2(const double Mw_i) const 
{
    int NumPar = 2;
    double params[] = {Mw(Mw_i), Mt()};

    if ( CacheCheck(Li2_MW2toMTOP2_cache, NumPar, params) )
        return Li2_MW2toMTOP2_cache[NumPar];
    else {
        double newResult = PolyLog.Li2( Mw2(Mw_i)/Mt()/Mt() ).real();
        newCacheForDouble(Li2_MW2toMTOP2_cache, NumPar, params, newResult);
        return newResult;
    }
}


double EWSMcache::Li3_MW2toMTOP2(const double Mw_i) const 
{
    int NumPar = 2;
    double params[] = {Mw(Mw_i), Mt()};

    if ( CacheCheck(Li3_MW2toMTOP2_cache, NumPar, params) )
        return Li3_MW2toMTOP2_cache[NumPar];
    else {
        double newResult = PolyLog.Li3(Mw2(Mw_i)/Mt()/Mt());
        newCacheForDouble(Li3_MW2toMTOP2_cache, NumPar, params, newResult);
        return newResult;
    }
}


double EWSMcache::Li3_for_F1(const double Mw_i) const 
{
    int NumPar = 2;
    double params[] = {Mw(Mw_i), Mt()};

    if ( CacheCheck(Li3_for_F1_cache, NumPar, params) )
        return Li3_for_F1_cache[NumPar];
    else {
        double tmp = Mw2(Mw_i)/Mt()/Mt();
        double newResult = PolyLog.Li3(-tmp/(1.0 - tmp));
        newCacheForDouble(Li3_for_F1_cache, NumPar, params, newResult);
        return newResult;
    }
}


double EWSMcache::A0_Mz2_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    if ( CacheCheck(A0_Mz2_Mw2_cache, NumPar, params) )
        return A0_Mz2_Mw2_cache[NumPar];
    else {
        double newResult = PV.A0(Mz2(), Mw2(Mw_i));
        newCacheForDouble(A0_Mz2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }    
}


double EWSMcache::A0_Mz2_mh2() const
{
    int NumPar = 2;
    double params[] = {Mz(), mh()};

    if ( CacheCheck(A0_Mz2_mh2_cache, NumPar, params) )
        return A0_Mz2_mh2_cache[NumPar];
    else {
        double newResult = PV.A0(Mz2(), mh2());
        newCacheForDouble(A0_Mz2_mh2_cache, NumPar, params, newResult);
        return newResult;
    }      
}


double EWSMcache::A0_Mw2_Mz2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mw(Mw_i), Mz()};

    if ( CacheCheck(A0_Mw2_Mz2_cache, NumPar, params) )
        return A0_Mw2_Mz2_cache[NumPar];
    else {
        double newResult = PV.A0(Mw2(Mw_i), Mz2());
        newCacheForDouble(A0_Mw2_Mz2_cache, NumPar, params, newResult);
        return newResult;
    }         
}


double EWSMcache::A0_Mw2_mh2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mw(Mw_i), mh()};

    if ( CacheCheck(A0_Mw2_mh2_cache, NumPar, params) )
        return A0_Mw2_mh2_cache[NumPar];
    else {
        double newResult = PV.A0(Mw2(Mw_i), mh2());
        newCacheForDouble(A0_Mw2_mh2_cache, NumPar, params, newResult);
        return newResult;
    }     
}


double EWSMcache::A0_Mz2_Mz2() const
{
    int NumPar = 1;
    double params[] = {Mz()};

    if ( CacheCheck(A0_Mz2_Mz2_cache, NumPar, params) )
        return A0_Mz2_Mz2_cache[NumPar];
    else {
        double newResult = PV.A0(Mz2(), Mz2());
        newCacheForDouble(A0_Mz2_Mz2_cache, NumPar, params, newResult);
        return newResult;
    }     
}


double EWSMcache::A0_Mw2_Mw2(const double Mw_i) const
{
    int NumPar = 1;
    double params[] = {Mw(Mw_i)};

    if ( CacheCheck(A0_Mw2_Mw2_cache, NumPar, params) )
        return A0_Mw2_Mw2_cache[NumPar];
    else {
        double newResult = PV.A0(Mw2(Mw_i), Mw2(Mw_i));
        newCacheForDouble(A0_Mw2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }       
}


complex EWSMcache::B0_Mz2_Mw2_mh2_Mw2(const double Mw_i) const
{
    int NumPar = 3;
    double params[] = {Mz(), Mw(Mw_i), mh()};

    if ( CacheCheck(B0_Mz2_Mw2_mh2_Mw2_cache, NumPar, params) )
        return complex(B0_Mz2_Mw2_mh2_Mw2_cache[NumPar],
                       B0_Mz2_Mw2_mh2_Mw2_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0(Mz2(), Mw2(Mw_i), mh2(), Mw2(Mw_i));
        newCacheForComplex(B0_Mz2_Mw2_mh2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::B0_Mz2_0_mh2_Mw2(const double Mw_i) const
{
    int NumPar = 3;
    double params[] = {Mz(), mh(), Mw(Mw_i)};

    if ( CacheCheck(B0_Mz2_0_mh2_Mw2_cache, NumPar, params) )
        return complex(B0_Mz2_0_mh2_Mw2_cache[NumPar],
                       B0_Mz2_0_mh2_Mw2_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0(Mz2(), 0.0, mh2(), Mw2(Mw_i));
        newCacheForComplex(B0_Mz2_0_mh2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mw2_Mz2_Mt2_Mt2(const double Mw_i) const
{
    int NumPar = 3;
    double params[] = {Mw(Mw_i), Mz(), Mt()};

    if ( CacheCheck(B0_Mw2_Mz2_Mt2_Mt2_cache, NumPar, params) )
        return complex(B0_Mw2_Mz2_Mt2_Mt2_cache[NumPar],
                       B0_Mw2_Mz2_Mt2_Mt2_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0(Mw2(Mw_i), Mz2(), Mt2(), Mt2());
        newCacheForComplex(B0_Mw2_Mz2_Mt2_Mt2_cache, NumPar, params, newResult);
        return newResult;
    }    
}


complex EWSMcache::B0_Mz2_Mz2_Mw2_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    if ( CacheCheck(B0_Mz2_Mz2_Mw2_Mw2_cache, NumPar, params) )
        return complex(B0_Mz2_Mz2_Mw2_Mw2_cache[NumPar],
                       B0_Mz2_Mz2_Mw2_Mw2_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0(Mz2(), Mz2(), Mw2(Mw_i), Mw2(Mw_i));
        newCacheForComplex(B0_Mz2_Mz2_Mw2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mz2_Mz2_mh2_Mz2() const
{
    int NumPar = 2;
    double params[] = {Mz(), mh()};

    if ( CacheCheck(B0_Mz2_Mz2_mh2_Mz2_cache, NumPar, params) )
        return complex(B0_Mz2_Mz2_mh2_Mz2_cache[NumPar],
                       B0_Mz2_Mz2_mh2_Mz2_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0(Mz2(), Mz2(), mh2(), Mz2());
        newCacheForComplex(B0_Mz2_Mz2_mh2_Mz2_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mz2_Mw2_Mz2_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    if ( CacheCheck(B0_Mz2_Mw2_Mz2_Mw2_cache, NumPar, params) )
        return complex(B0_Mz2_Mw2_Mz2_Mw2_cache[NumPar],
                       B0_Mz2_Mw2_Mz2_Mw2_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0(Mz2(), Mw2(Mw_i), Mz2(), Mw2(Mw_i));
        newCacheForComplex(B0_Mz2_Mw2_Mz2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mz2_Mw2_0_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    if ( CacheCheck(B0_Mz2_Mw2_0_Mw2_cache, NumPar, params) )
        return complex(B0_Mz2_Mw2_0_Mw2_cache[NumPar],
                       B0_Mz2_Mw2_0_Mw2_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0(Mz2(), Mw2(Mw_i), 0.0, Mw2(Mw_i));
        newCacheForComplex(B0_Mz2_Mw2_0_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mz2_0_Mz2_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    if ( CacheCheck(B0_Mz2_0_Mz2_Mw2_cache, NumPar, params) )
        return complex(B0_Mz2_0_Mz2_Mw2_cache[NumPar],
                       B0_Mz2_0_Mz2_Mw2_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0(Mz2(), 0.0, Mz2(), Mw2(Mw_i));
        newCacheForComplex(B0_Mz2_0_Mz2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mz2_0_0_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    if ( CacheCheck(B0_Mz2_0_0_Mw2_cache, NumPar, params) )
        return complex(B0_Mz2_0_0_Mw2_cache[NumPar],
                       B0_Mz2_0_0_Mw2_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0(Mz2(), 0.0, 0.0, Mw2(Mw_i));
        newCacheForComplex(B0_Mz2_0_0_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mw2_Mz2_Mw2_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mw(Mw_i), Mz()};

    if ( CacheCheck(B0_Mw2_Mz2_Mw2_Mw2_cache, NumPar, params) )
        return complex(B0_Mw2_Mz2_Mw2_Mw2_cache[NumPar],
                       B0_Mw2_Mz2_Mw2_Mw2_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0(Mw2(Mw_i), Mz2(), Mw2(Mw_i), Mw2(Mw_i));
        newCacheForComplex(B0_Mw2_Mz2_Mw2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mw2_Mw2_Mz2_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mw(Mw_i), Mz()};

    if ( CacheCheck(B0_Mw2_Mw2_Mz2_Mw2_cache, NumPar, params) )
        return complex(B0_Mw2_Mw2_Mz2_Mw2_cache[NumPar],
                       B0_Mw2_Mw2_Mz2_Mw2_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0(Mw2(Mw_i), Mw2(Mw_i), Mz2(), Mw2(Mw_i));
        newCacheForComplex(B0_Mw2_Mw2_Mz2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mw2_Mw2_mh2_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mw(Mw_i), mh()};

    if ( CacheCheck(B0_Mw2_Mw2_mh2_Mw2_cache, NumPar, params) )
        return complex(B0_Mw2_Mw2_mh2_Mw2_cache[NumPar],
                       B0_Mw2_Mw2_mh2_Mw2_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0(Mw2(Mw_i), Mw2(Mw_i), mh2(), Mw2(Mw_i));
        newCacheForComplex(B0_Mw2_Mw2_mh2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mw2_Mw2_0_Mw2(const double Mw_i) const
{
    int NumPar = 1;
    double params[] = {Mw(Mw_i)};
    
    if ( CacheCheck(B0_Mw2_Mw2_0_Mw2_cache, NumPar, params) )
        return complex(B0_Mw2_Mw2_0_Mw2_cache[NumPar],
                       B0_Mw2_Mw2_0_Mw2_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0(Mw2(Mw_i), Mw2(Mw_i), 0.0, Mw2(Mw_i));
        newCacheForComplex(B0_Mw2_Mw2_0_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mz2_Mz2_ml2_ml2(const StandardModel::lepton l) const
{
    int NumPar = 2;
    double params[] = {Mz(), ml(l)};

    if ( CacheCheck(B0_Mz2_Mz2_ml2_ml2_cache[l], NumPar, params) )
        return complex(B0_Mz2_Mz2_ml2_ml2_cache[l][NumPar],
                       B0_Mz2_Mz2_ml2_ml2_cache[l][NumPar+1], false);
    else {
        complex newResult = PV.B0(Mz2(), Mz2(), ml2(l), ml2(l));
        newCacheForComplex(B0_Mz2_Mz2_ml2_ml2_cache[l], NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0_Mz2_Mz2_mq2_mq2(const QCD::quark q) const
{
    int NumPar = 2;
    double params[] = {Mz(), mq(q, Mz())};

    if ( CacheCheck(B0_Mz2_Mz2_mq2_mq2_cache[q], NumPar, params) )
        return complex(B0_Mz2_Mz2_mq2_mq2_cache[q][NumPar],
                       B0_Mz2_Mz2_mq2_mq2_cache[q][NumPar+1], false);
    else {
        complex newResult = PV.B0(Mz2(), Mz2(), mq2(q, Mz()), mq2(q, Mz()));
        newCacheForComplex(B0_Mz2_Mz2_mq2_mq2_cache[q], NumPar, params, newResult);
        return newResult;
    }    
}


complex EWSMcache::B0p_Mz2_0_mh2_Mw2(const double Mw_i) const
{
    int NumPar = 3;
    double params[] = {Mz(), mh(), Mw(Mw_i)};

    if ( CacheCheck(B0p_Mz2_0_mh2_Mw2_cache, NumPar, params) )
        return complex(B0p_Mz2_0_mh2_Mw2_cache[NumPar],
                       B0p_Mz2_0_mh2_Mw2_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0p(Mz2(), 0.0, mh2(), Mw2(Mw_i));
        newCacheForComplex(B0p_Mz2_0_mh2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0p_Mz2_Mz2_mh2_Mz2() const
{
    int NumPar = 2;
    double params[] = {Mz(), mh()};

    if ( CacheCheck(B0p_Mz2_Mz2_mh2_Mz2_cache, NumPar, params) )
        return complex(B0p_Mz2_Mz2_mh2_Mz2_cache[NumPar],
                       B0p_Mz2_Mz2_mh2_Mz2_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0p(Mz2(), Mz2(), mh2(), Mz2());
        newCacheForComplex(B0p_Mz2_Mz2_mh2_Mz2_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0p_Mz2_0_Mz2_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    if ( CacheCheck(B0p_Mz2_0_Mz2_Mw2_cache, NumPar, params) )
        return complex(B0p_Mz2_0_Mz2_Mw2_cache[NumPar],
                       B0p_Mz2_0_Mz2_Mw2_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0p(Mz2(), 0.0, Mz2(), Mw2(Mw_i));
        newCacheForComplex(B0p_Mz2_0_Mz2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0p_Mz2_Mz2_Mw2_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    if ( CacheCheck(B0p_Mz2_Mz2_Mw2_Mw2_cache, NumPar, params) )
        return complex(B0p_Mz2_Mz2_Mw2_Mw2_cache[NumPar],
                       B0p_Mz2_Mz2_Mw2_Mw2_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0p(Mz2(), Mz2(), Mw2(Mw_i), Mw2(Mw_i));
        newCacheForComplex(B0p_Mz2_Mz2_Mw2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0p_Mw2_Mw2_Mz2_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mw(Mw_i), Mz()};

    if ( CacheCheck(B0p_Mw2_Mw2_Mz2_Mw2_cache, NumPar, params) )
        return complex(B0p_Mw2_Mw2_Mz2_Mw2_cache[NumPar],
                       B0p_Mw2_Mw2_Mz2_Mw2_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0p(Mw2(Mw_i), Mw2(Mw_i), Mz2(), Mw2(Mw_i));
        newCacheForComplex(B0p_Mw2_Mw2_Mz2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0p_Mw2_Mw2_mh2_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mw(Mw_i), mh()};

    if ( CacheCheck(B0p_Mw2_Mw2_mh2_Mw2_cache, NumPar, params) )
        return complex(B0p_Mw2_Mw2_mh2_Mw2_cache[NumPar],
                       B0p_Mw2_Mw2_mh2_Mw2_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0p(Mw2(Mw_i), Mw2(Mw_i), mh2(), Mw2(Mw_i));
        newCacheForComplex(B0p_Mw2_Mw2_mh2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0p_Mw2_Mw2_0_Mw2(const double Mw_i) const
{
    int NumPar = 1;
    double params[] = {Mw(Mw_i)};

    if ( CacheCheck(B0p_Mw2_Mw2_0_Mw2_cache, NumPar, params) )
        return complex(B0p_Mw2_Mw2_0_Mw2_cache[NumPar],
                       B0p_Mw2_Mw2_0_Mw2_cache[NumPar+1], false);
    else {
        complex newResult = PV.B0p(Mw2(Mw_i), Mw2(Mw_i), 0.0, Mw2(Mw_i));
        newCacheForComplex(B0p_Mw2_Mw2_0_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0p_Mz2_Mz2_ml2_ml2(const StandardModel::lepton l) const
{
    int NumPar = 2;
    double params[] = {Mz(), ml(l)};

    if ( CacheCheck(B0p_Mz2_Mz2_ml2_ml2_cache[l], NumPar, params) )
        return complex(B0p_Mz2_Mz2_ml2_ml2_cache[l][NumPar],
                       B0p_Mz2_Mz2_ml2_ml2_cache[l][NumPar+1], false);
    else {
        complex newResult = PV.B0p(Mz2(), Mz2(), ml2(l), ml2(l));
        newCacheForComplex(B0p_Mz2_Mz2_ml2_ml2_cache[l], NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::B0p_Mz2_Mz2_mq2_mq2(const QCD::quark q) const
{
    int NumPar = 2;
    double params[] = {Mz(), mq(q, Mz())};

    if ( CacheCheck(B0p_Mz2_Mz2_mq2_mq2_cache[q], NumPar, params) )
        return complex(B0p_Mz2_Mz2_mq2_mq2_cache[q][NumPar],
                       B0p_Mz2_Mz2_mq2_mq2_cache[q][NumPar+1], false);
    else {
        complex newResult = PV.B0p(Mz2(), Mz2(), mq2(q, Mz()), mq2(q, Mz()));
        newCacheForComplex(B0p_Mz2_Mz2_mq2_mq2_cache[q], NumPar, params, newResult);
        return newResult;
    }    
}

        
complex EWSMcache::B1_Mz2_0_ml2_mlprime2(const int gen) const
{
    int NumPar = 3;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mz(), mf, mfprime};

    if ( CacheCheck(B1_Mz2_0_ml2_mlprime2_cache[gen], NumPar, params) )
        return complex(B1_Mz2_0_ml2_mlprime2_cache[gen][NumPar],
                       B1_Mz2_0_ml2_mlprime2_cache[gen][NumPar+1], false);
    else {
        double mf2 = mf*mf;
        double mfprime2 = mfprime*mfprime;
        complex newResult = PV.B1(Mz2(), 0.0, mf2, mfprime2);
        newCacheForComplex(B1_Mz2_0_ml2_mlprime2_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}   


complex EWSMcache::B1_Mz2_0_mq2_mqprime2(const int gen) const
{
    int NumPar = 3;
    double mf = mq((QCD::quark)(2*gen), Mz());    
    double mfprime = mq((QCD::quark)(2*gen+1), Mz());
    double params[] = {Mz(), mf, mfprime};
    
    if ( CacheCheck(B1_Mz2_0_mq2_mqprime2_cache[gen], NumPar, params) )
        return complex(B1_Mz2_0_mq2_mqprime2_cache[gen][NumPar],
                       B1_Mz2_0_mq2_mqprime2_cache[gen][NumPar+1], false);
    else {
        double mf2 = mf*mf;
        double mfprime2 = mfprime*mfprime;
        complex newResult = PV.B1(Mz2(), 0.0, mf2, mfprime2);
        newCacheForComplex(B1_Mz2_0_mq2_mqprime2_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
} 


complex EWSMcache::B1_Mz2_0_mlprime2_ml2(const int gen) const
{
    int NumPar = 3;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mz(), mf, mfprime};

    if ( CacheCheck(B1_Mz2_0_mlprime2_ml2_cache[gen], NumPar, params) )
        return complex(B1_Mz2_0_mlprime2_ml2_cache[gen][NumPar],
                       B1_Mz2_0_mlprime2_ml2_cache[gen][NumPar+1], false);
    else {
        double mf2 = mf*mf;
        double mfprime2 = mfprime*mfprime;
        complex newResult = PV.B1(Mz2(), 0.0, mfprime2, mf2);
        newCacheForComplex(B1_Mz2_0_mlprime2_ml2_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
} 


complex EWSMcache::B1_Mz2_0_mqprime2_mq2(const int gen) const
{
    int NumPar = 3;
    double mf = mq((QCD::quark)(2*gen), Mz());    
    double mfprime = mq((QCD::quark)(2*gen+1), Mz());
    double params[] = {Mz(), mf, mfprime};

    if ( CacheCheck(B1_Mz2_0_mqprime2_mq2_cache[gen], NumPar, params) )
        return complex(B1_Mz2_0_mqprime2_mq2_cache[gen][NumPar],
                       B1_Mz2_0_mqprime2_mq2_cache[gen][NumPar+1], false);
    else {
        double mf2 = mf*mf;
        double mfprime2 = mfprime*mfprime;
        complex newResult = PV.B1(Mz2(), 0.0, mfprime2, mf2);
        newCacheForComplex(B1_Mz2_0_mqprime2_mq2_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::B1_Mz2_Mw2_ml2_mlprime2(const int gen, const double Mw_i) const
{
    int NumPar = 4;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mz(), Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(B1_Mz2_Mw2_ml2_mlprime2_cache[gen], NumPar, params) )
        return complex(B1_Mz2_Mw2_ml2_mlprime2_cache[gen][NumPar],
                       B1_Mz2_Mw2_ml2_mlprime2_cache[gen][NumPar+1], false);
    else {
        double mf2 = mf*mf;
        double mfprime2 = mfprime*mfprime;
        complex newResult = PV.B1(Mz2(), Mw2(Mw_i), mf2, mfprime2);
        newCacheForComplex(B1_Mz2_Mw2_ml2_mlprime2_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}  


complex EWSMcache::B1_Mz2_Mw2_mq2_mqprime2(const int gen, const double Mw_i) const
{
    int NumPar = 4;
    double mf = mq((QCD::quark)(2*gen), Mz());    
    double mfprime = mq((QCD::quark)(2*gen+1), Mz());
    double params[] = {Mz(), Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(B1_Mz2_Mw2_mq2_mqprime2_cache[gen], NumPar, params) )
        return complex(B1_Mz2_Mw2_mq2_mqprime2_cache[gen][NumPar],
                       B1_Mz2_Mw2_mq2_mqprime2_cache[gen][NumPar+1], false);
    else {
        double mf2 = mf*mf;
        double mfprime2 = mfprime*mfprime;
        complex newResult = PV.B1(Mz2(), Mw2(Mw_i), mf2, mfprime2);
        newCacheForComplex(B1_Mz2_Mw2_mq2_mqprime2_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
} 


complex EWSMcache::B1_Mz2_Mw2_mlprime2_ml2(const int gen, const double Mw_i) const
{
    int NumPar = 4;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mz(), Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(B1_Mz2_Mw2_mlprime2_ml2_cache[gen], NumPar, params) )
        return complex(B1_Mz2_Mw2_mlprime2_ml2_cache[gen][NumPar],
                       B1_Mz2_Mw2_mlprime2_ml2_cache[gen][NumPar+1], false);
    else {
        double mf2 = mf*mf;
        double mfprime2 = mfprime*mfprime;
        complex newResult = PV.B1(Mz2(), Mw2(Mw_i), mfprime2, mf2);
        newCacheForComplex(B1_Mz2_Mw2_mlprime2_ml2_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}  


complex EWSMcache::B1_Mz2_Mw2_mqprime2_mq2(const int gen, const double Mw_i) const
{
    int NumPar = 4;
    double mf = mq((QCD::quark)(2*gen), Mz());    
    double mfprime = mq((QCD::quark)(2*gen+1), Mz());
    double params[] = {Mz(), Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(B1_Mz2_Mw2_mqprime2_mq2_cache[gen], NumPar, params) )
        return complex(B1_Mz2_Mw2_mqprime2_mq2_cache[gen][NumPar],
                       B1_Mz2_Mw2_mqprime2_mq2_cache[gen][NumPar+1], false);
    else {
        double mf2 = mf*mf;
        double mfprime2 = mfprime*mfprime;
        complex newResult = PV.B1(Mz2(), Mw2(Mw_i), mfprime2, mf2);
        newCacheForComplex(B1_Mz2_Mw2_mqprime2_mq2_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
} 


complex EWSMcache::B1p_Mw2_Mw2_ml2_mlprime2(const int gen, const double Mw_i) const
{
    int NumPar = 3;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(B1p_Mw2_Mw2_ml2_mlprime2_cache[gen], NumPar, params) )
        return complex(B1p_Mw2_Mw2_ml2_mlprime2_cache[gen][NumPar],
                       B1p_Mw2_Mw2_ml2_mlprime2_cache[gen][NumPar+1], false);
    else {
        double mf2 = mf*mf;
        double mfprime2 = mfprime*mfprime;
        complex newResult = PV.B1p(Mw2(Mw_i), Mw2(Mw_i), mf2, mfprime2);
        newCacheForComplex(B1p_Mw2_Mw2_ml2_mlprime2_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::B1p_Mw2_Mw2_mq2_mqprime2(const int gen, const double Mw_i) const
{
    int NumPar = 3;
    double mf = mq((QCD::quark)(2*gen), Mz());    
    double mfprime = mq((QCD::quark)(2*gen+1), Mz());
    double params[] = {Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(B1p_Mw2_Mw2_mq2_mqprime2_cache[gen], NumPar, params) )
        return complex(B1p_Mw2_Mw2_mq2_mqprime2_cache[gen][NumPar],
                       B1p_Mw2_Mw2_mq2_mqprime2_cache[gen][NumPar+1], false);
    else {
        double mf2 = mf*mf;
        double mfprime2 = mfprime*mfprime;
        complex newResult = PV.B1p(Mw2(Mw_i), Mw2(Mw_i), mf2, mfprime2);
        newCacheForComplex(B1p_Mw2_Mw2_mq2_mqprime2_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::B1p_Mw2_Mw2_mlprime2_ml2(const int gen, const double Mw_i) const
{
    int NumPar = 3;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(B1p_Mw2_Mw2_mlprime2_ml2_cache[gen], NumPar, params) )
        return complex(B1p_Mw2_Mw2_mlprime2_ml2_cache[gen][NumPar],
                       B1p_Mw2_Mw2_mlprime2_ml2_cache[gen][NumPar+1], false);
    else {
        double mf2 = mf*mf;
        double mfprime2 = mfprime*mfprime;
        complex newResult = PV.B1p(Mw2(Mw_i), Mw2(Mw_i), mfprime2, mf2);
        newCacheForComplex(B1p_Mw2_Mw2_mlprime2_ml2_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::B1p_Mw2_Mw2_mqprime2_mq2(const int gen, const double Mw_i) const
{
    int NumPar = 3;
    double mf = mq((QCD::quark)(2*gen), Mz());    
    double mfprime = mq((QCD::quark)(2*gen+1), Mz());
    double params[] = {Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(B1p_Mw2_Mw2_mqprime2_mq2_cache[gen], NumPar, params) )
        return complex(B1p_Mw2_Mw2_mqprime2_mq2_cache[gen][NumPar],
                       B1p_Mw2_Mw2_mqprime2_mq2_cache[gen][NumPar+1], false);
    else {
        double mf2 = mf*mf;
        double mfprime2 = mfprime*mfprime;
        complex newResult = PV.B1p(Mw2(Mw_i), Mw2(Mw_i), mfprime2, mf2);
        newCacheForComplex(B1p_Mw2_Mw2_mqprime2_mq2_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::Bf_Mz2_Mz2_ml2_ml2(const StandardModel::lepton l) const
{
    int NumPar = 2;
    double params[] = {Mz(), ml(l)};

    if ( CacheCheck(Bf_Mz2_Mz2_ml2_ml2_cache[l], NumPar, params) )
        return complex(Bf_Mz2_Mz2_ml2_ml2_cache[l][NumPar],
                       Bf_Mz2_Mz2_ml2_ml2_cache[l][NumPar+1], false);
    else {
        complex newResult = PV.Bf(Mz2(), Mz2(), ml2(l), ml2(l));
        newCacheForComplex(Bf_Mz2_Mz2_ml2_ml2_cache[l], NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::Bf_Mz2_Mz2_mq2_mq2(const QCD::quark q) const
{
    int NumPar = 2;
    double params[] = {Mz(), mq(q, Mz())};

    if ( CacheCheck(Bf_Mz2_Mz2_mq2_mq2_cache[q], NumPar, params) )
        return complex(Bf_Mz2_Mz2_mq2_mq2_cache[q][NumPar],
                       Bf_Mz2_Mz2_mq2_mq2_cache[q][NumPar+1], false);
    else {
        complex newResult = PV.Bf(Mz2(), Mz2(), mq2(q, Mz()), mq2(q, Mz()));
        newCacheForComplex(Bf_Mz2_Mz2_mq2_mq2_cache[q], NumPar, params, newResult);
        return newResult;
    }    
}


complex EWSMcache::Bf_Mz2_0_ml2_ml2(const StandardModel::lepton l) const
{
    int NumPar = 2;
    double params[] = {Mz(), ml(l)};
    if (ml(l)==0.0)
        throw std::runtime_error("Error in EWSMcache::Bf_Mz_0_ml_ml()"); 
    
    if ( CacheCheck(Bf_Mz2_0_ml2_ml2_cache[l], NumPar, params) )
        return complex(Bf_Mz2_0_ml2_ml2_cache[l][NumPar],
                       Bf_Mz2_0_ml2_ml2_cache[l][NumPar+1], false);
    else {
        complex newResult = PV.Bf(Mz2(), 0.0, ml2(l), ml2(l));
            newCacheForComplex(Bf_Mz2_0_ml2_ml2_cache[l], NumPar, params, newResult);
            return newResult;
    }
}


complex EWSMcache::Bf_Mz2_0_mq2_mq2(const QCD::quark q) const
{
    int NumPar = 2;
    double params[] = {Mz(), mq(q, Mz())};
    if (mq(q, Mz())==0.0)
        throw std::runtime_error("Error in EWSMcache::Bf_Mz_0_mq_mq()"); 
    
    if ( CacheCheck(Bf_Mz2_0_mq2_mq2_cache[q], NumPar, params) )
        return complex(Bf_Mz2_0_mq2_mq2_cache[q][NumPar],
                       Bf_Mz2_0_mq2_mq2_cache[q][NumPar+1], false);
    else {
        complex newResult = PV.Bf(Mz2(), 0.0, mq2(q, Mz()), mq2(q, Mz()));
        newCacheForComplex(Bf_Mz2_0_mq2_mq2_cache[q], NumPar, params, newResult);
        return newResult;
    }    
}


complex EWSMcache::Bf_Mz2_Mw2_mlprime2_ml2(const int gen, const double Mw_i) const
{
    int NumPar = 4;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mz(), Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(Bf_Mz2_Mw2_mlprime2_ml2_cache[gen], NumPar, params) )
        return complex(Bf_Mz2_Mw2_mlprime2_ml2_cache[gen][NumPar],
                       Bf_Mz2_Mw2_mlprime2_ml2_cache[gen][NumPar+1], false);
    else {
        double mf2 = mf*mf;
        double mfprime2 = mfprime*mfprime;
        complex newResult = PV.Bf(Mz2(), Mw2(Mw_i), mfprime2, mf2);
        newCacheForComplex(Bf_Mz2_Mw2_mlprime2_ml2_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}        
        
        
complex EWSMcache::Bf_Mz2_Mw2_mqprime2_mq2(const int gen, const double Mw_i) const
{
    int NumPar = 4;
    double mf = mq((QCD::quark)(2*gen), Mz());    
    double mfprime = mq((QCD::quark)(2*gen+1), Mz());
    double params[] = {Mz(), Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(Bf_Mz2_Mw2_mqprime2_mq2_cache[gen], NumPar, params) )
        return complex(Bf_Mz2_Mw2_mqprime2_mq2_cache[gen][NumPar],
                       Bf_Mz2_Mw2_mqprime2_mq2_cache[gen][NumPar+1], false);
    else {
        double mf2 = mf*mf;
        double mfprime2 = mfprime*mfprime;
        complex newResult = PV.Bf(Mz2(), Mw2(Mw_i), mfprime2, mf2);
        newCacheForComplex(Bf_Mz2_Mw2_mqprime2_mq2_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::Bf_Mz2_0_mlprime2_ml2(const int gen) const
{
    int NumPar = 3;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mz(), mf, mfprime};

    if ( CacheCheck(Bf_Mz2_0_mlprime2_ml2_cache[gen], NumPar, params) )
        return complex(Bf_Mz2_0_mlprime2_ml2_cache[gen][NumPar],
                       Bf_Mz2_0_mlprime2_ml2_cache[gen][NumPar+1], false);
    else {
        double mf2 = mf*mf;
        double mfprime2 = mfprime*mfprime;
        complex newResult = PV.Bf(Mz2(), 0.0, mfprime2, mf2);
        newCacheForComplex(Bf_Mz2_0_mlprime2_ml2_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::Bf_Mz2_0_mqprime2_mq2(const int gen) const
{
    int NumPar = 3;
    double mf = mq((QCD::quark)(2*gen), Mz());    
    double mfprime = mq((QCD::quark)(2*gen+1), Mz());
    double params[] = {Mz(), mf, mfprime};

    if ( CacheCheck(Bf_Mz2_0_mqprime2_mq2_cache[gen], NumPar, params) )
        return complex(Bf_Mz2_0_mqprime2_mq2_cache[gen][NumPar],
                       Bf_Mz2_0_mqprime2_mq2_cache[gen][NumPar+1], false);
    else {
        double mf2 = mf*mf;
        double mfprime2 = mfprime*mfprime;
        complex newResult = PV.Bf(Mz2(), 0.0, mfprime2, mf2);
        newCacheForComplex(Bf_Mz2_0_mqprime2_mq2_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::Bf_Mw2_Mw2_mlprime2_ml2(const int gen, const double Mw_i) const
{
    int NumPar = 3;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(Bf_Mw2_Mw2_mlprime2_ml2_cache[gen], NumPar, params) )
        return complex(Bf_Mw2_Mw2_mlprime2_ml2_cache[gen][NumPar],
                       Bf_Mw2_Mw2_mlprime2_ml2_cache[gen][NumPar+1], false);
    else {
        double mf2 = mf*mf;
        double mfprime2 = mfprime*mfprime;
        complex newResult = PV.Bf(Mw2(Mw_i), Mw2(Mw_i), mfprime2, mf2);
        newCacheForComplex(Bf_Mw2_Mw2_mlprime2_ml2_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::Bf_Mw2_Mw2_mqprime2_mq2(const int gen, const double Mw_i) const
{
    int NumPar = 3;
    double mf = mq((QCD::quark)(2*gen), Mz());    
    double mfprime = mq((QCD::quark)(2*gen+1), Mz());
    double params[] = {Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(Bf_Mw2_Mw2_mqprime2_mq2_cache[gen], NumPar, params) )
        return complex(Bf_Mw2_Mw2_mqprime2_mq2_cache[gen][NumPar],
                       Bf_Mw2_Mw2_mqprime2_mq2_cache[gen][NumPar+1], false);
    else {
        double mf2 = mf*mf;
        double mfprime2 = mfprime*mfprime;
        complex newResult = PV.Bf(Mw2(Mw_i), Mw2(Mw_i), mfprime2, mf2);
        newCacheForComplex(Bf_Mw2_Mw2_mqprime2_mq2_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::Bfp_Mz2_Mz2_ml2_ml2(const StandardModel::lepton l) const
{
    int NumPar = 2;
    double params[] = {Mz(), ml(l)};

    if ( CacheCheck(Bfp_Mz2_Mz2_ml2_ml2_cache[l], NumPar, params) )
        return complex(Bfp_Mz2_Mz2_ml2_ml2_cache[l][NumPar],
                       Bfp_Mz2_Mz2_ml2_ml2_cache[l][NumPar+1], false);
    else {
        complex newResult = PV.Bfp(Mz2(), Mz2(), ml2(l), ml2(l));
        newCacheForComplex(Bfp_Mz2_Mz2_ml2_ml2_cache[l], NumPar, params, newResult);
        return newResult;
    }
}


complex EWSMcache::Bfp_Mz2_Mz2_mq2_mq2(const QCD::quark q) const
{
    int NumPar = 2;
    double params[] = {Mz(), mq(q, Mz())};

    if ( CacheCheck(Bfp_Mz2_Mz2_mq2_mq2_cache[q], NumPar, params) )
        return complex(Bfp_Mz2_Mz2_mq2_mq2_cache[q][NumPar],
                       Bfp_Mz2_Mz2_mq2_mq2_cache[q][NumPar+1], false);
    else {
        complex newResult = PV.Bfp(Mz2(), Mz2(), mq2(q, Mz()), mq2(q, Mz()));
        newCacheForComplex(Bfp_Mz2_Mz2_mq2_mq2_cache[q], NumPar, params, newResult);
        return newResult;
    }    
}


complex EWSMcache::Bfp_Mw2_Mw2_mlprime2_ml2(const int gen, const double Mw_i) const
{
    int NumPar = 3;
    double mf = ml((StandardModel::lepton)(2*gen));    
    double mfprime = ml((StandardModel::lepton)(2*gen+1));
    double params[] = {Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(Bfp_Mw2_Mw2_mlprime2_ml2_cache[gen], NumPar, params) )
        return complex(Bfp_Mw2_Mw2_mlprime2_ml2_cache[gen][NumPar],
                       Bfp_Mw2_Mw2_mlprime2_ml2_cache[gen][NumPar+1], false);
    else {
        double mf2 = mf*mf;
        double mfprime2 = mfprime*mfprime;
        complex newResult = PV.Bfp(Mw2(Mw_i), Mw2(Mw_i), mfprime2, mf2);
        newCacheForComplex(Bfp_Mw2_Mw2_mlprime2_ml2_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::Bfp_Mw2_Mw2_mqprime2_mq2(const int gen, const double Mw_i) const
{
    int NumPar = 3;
    double mf = mq((QCD::quark)(2*gen), Mz());    
    double mfprime = mq((QCD::quark)(2*gen+1), Mz());
    double params[] = {Mw(Mw_i), mf, mfprime};

    if ( CacheCheck(Bfp_Mw2_Mw2_mqprime2_mq2_cache[gen], NumPar, params) )
        return complex(Bfp_Mw2_Mw2_mqprime2_mq2_cache[gen][NumPar],
                       Bfp_Mw2_Mw2_mqprime2_mq2_cache[gen][NumPar+1], false);
    else {
        double mf2 = mf*mf;
        double mfprime2 = mfprime*mfprime;
        complex newResult = PV.Bfp(Mw2(Mw_i), Mw2(Mw_i), mfprime2, mf2);
        newCacheForComplex(Bfp_Mw2_Mw2_mqprime2_mq2_cache[gen], NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::C0_Mz2_Mw2_Mt2_Mw2(const double Mw_i) const
{
    int NumPar = 3;
    double params[] = {Mz(), Mw(Mw_i), Mt()};

    if ( CacheCheck(C0_Mz2_Mw2_Mt2_Mw2_cache, NumPar, params) )
        return complex(C0_Mz2_Mw2_Mt2_Mw2_cache[NumPar],
                       C0_Mz2_Mw2_Mt2_Mw2_cache[NumPar+1], false);
    else {
        complex newResult = PV.C0(Mz2(), Mw2(Mw_i), Mt2(), Mw2(Mw_i));
        newCacheForComplex(C0_Mz2_Mw2_Mt2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::C0_Mz2_Mt2_Mw2_Mt2(const double Mw_i) const
{
    int NumPar = 3;
    double params[] = {Mz(), Mt(), Mw(Mw_i)};

    if ( CacheCheck(C0_Mz2_Mt2_Mw2_Mt2_cache, NumPar, params) )
        return complex(C0_Mz2_Mt2_Mw2_Mt2_cache[NumPar],
                       C0_Mz2_Mt2_Mw2_Mt2_cache[NumPar+1], false);
    else {
        complex newResult = PV.C0(Mz2(), Mt2(), Mw2(Mw_i), Mt2());
        newCacheForComplex(C0_Mz2_Mt2_Mw2_Mt2_cache, NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::C0_Mz2_0_Mw2_0(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    if ( CacheCheck(C0_Mz2_0_Mw2_0_cache, NumPar, params) )
        return complex(C0_Mz2_0_Mw2_0_cache[NumPar],
                       C0_Mz2_0_Mw2_0_cache[NumPar+1], false);
    else {
        complex newResult = PV.C0(Mz2(), 0.0, Mw2(Mw_i), 0.0);
        newCacheForComplex(C0_Mz2_0_Mw2_0_cache, NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::C0_Mz2_Mw2_0_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mz(), Mw(Mw_i)};

    if ( CacheCheck(C0_Mz2_Mw2_0_Mw2_cache, NumPar, params) )
        return complex(C0_Mz2_Mw2_0_Mw2_cache[NumPar],
                       C0_Mz2_Mw2_0_Mw2_cache[NumPar+1], false);
    else {
        complex newResult = PV.C0(Mz2(), Mw2(Mw_i), 0.0, Mw2(Mw_i));
        newCacheForComplex(C0_Mz2_Mw2_0_Mw2_cache, NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::C0_Mw2_Mw2_0_Mz2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mw(Mw_i), Mz()};

    if ( CacheCheck(C0_Mw2_Mw2_0_Mz2_cache, NumPar, params) )
        return complex(C0_Mw2_Mw2_0_Mz2_cache[NumPar],
                       C0_Mw2_Mw2_0_Mz2_cache[NumPar+1], false);
    else {
        complex newResult = PV.C0(Mw2(Mw_i), Mw2(Mw_i), 0.0, Mz2());
        newCacheForComplex(C0_Mw2_Mw2_0_Mz2_cache, NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::C0_Mw2_0_Mz2_0(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mw(Mw_i), Mz()};

    if ( CacheCheck(C0_Mw2_0_Mz2_0_cache, NumPar, params) )
        return complex(C0_Mw2_0_Mz2_0_cache[NumPar],
                       C0_Mw2_0_Mz2_0_cache[NumPar+1], false);
    else {
        complex newResult = PV.C0(Mw2(Mw_i), 0.0, Mz2(), 0.0);
        newCacheForComplex(C0_Mw2_0_Mz2_0_cache, NumPar, params, newResult);
        return newResult;
    } 
}


complex EWSMcache::C0_Mz2_0_Mz2_0() const
{
    int NumPar = 1;
    double params[] = {Mz()};

    if ( CacheCheck(C0_Mz2_0_Mz2_0_cache, NumPar, params) )
        return complex(C0_Mz2_0_Mz2_0_cache[NumPar],
                       C0_Mz2_0_Mz2_0_cache[NumPar+1], false);
    else {
        complex newResult = PV.C0(Mz2(), 0.0, Mz2(), 0.0);
        newCacheForComplex(C0_Mz2_0_Mz2_0_cache, NumPar, params, newResult);
        return newResult;
    } 
}




