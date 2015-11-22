/* 
 * Copyright (C) 2012 HEPfit Collaboration
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
    FlagCacheInEWSMcache = true; // use caches in the current class
    //FlagCacheInEWSMcache = false;// do not use caches in the current class (for test)

    log2 = log(2.0);

    /* zeta functions */
    zeta2 = gsl_sf_zeta_int(2);
    zeta3 = gsl_sf_zeta_int(3);
    zeta4 = gsl_sf_zeta_int(4);
    zeta5 = gsl_sf_zeta_int(5);

    /* Constants for three-loop contribution */
    double Cl2_Pi_3 = Clausen.Cl2(M_PI / 3.0);
    S2 = 4.0 / 9.0 / sqrt(3.0) * Cl2_Pi_3;
    D3 = 6.0 * zeta3 - 15.0 / 4.0 * zeta4 - 6.0 * Cl2_Pi_3*Cl2_Pi_3;
    B4 = -1.76280008707377;
    //double Li4_1_2 = ??;
    //B4 = 16.0*Li4_1_2 - 4.0*zeta2*log2*log2 + 2.0/3.0*pow(log2,4.0) - 13.0/2.0*zeta4;

    // Initializations of the cache
    for (int i = 0; i < 12; ++i) {
        mf_atMz_cache[i] = 0.0;
        for (int j = 0; j < StandardModel::NumSMParamsForEWPO; ++j)
            mf_atMz_params_cache[i][j] = 0.0;
    }
}


////////////////////////////////////////////////////////////////////////

double EWSMcache::mf(const Particle f, const double mu, const orders order) const
{
    if (f.is("TOP"))
        return SM.getMtpole(); // the pole mass
    else if (f.is("QUARK") && !FlagDebug) {
        /* These codes are slow and not effective. */
        //if (mu == SM.getMz()) {
        //    if (FlagCacheInEWSMcache && order == FULLNNLO)
        //        if (SM.checkSMparams(mf_atMz_params_cache[f.getIndex()]))
        //            return mf_atMz_cache[f.getIndex()];
        //    mf_atMz_cache[f.getIndex()] = SM.Mrun(mu, f.getMass_scale(), f.getMass(), order);
        //    return mf_atMz_cache[f.getIndex()];
        //}
        return SM.Mrun(mu, f.getMass_scale(), f.getMass(), order);
    } else
        return f.getMass();
}


//////////////////////////////////////////////////////////////////////// 

double EWSMcache::logMZtoME() const
{
    int NumPar = 2;
    double params[] = {SM.getMz(), mf(SM.getLeptons(SM.ELECTRON))};

    if (CacheCheck(logMZtoME_cache, NumPar, params))
        return logMZtoME_cache[NumPar];
    else {
        double newResult = log(SM.getMz() / mf(SM.getLeptons(SM.ELECTRON)));
        newCacheForDouble(logMZtoME_cache, NumPar, params, newResult);
        return newResult;
    }
}

double EWSMcache::logMZtoMMU() const
{
    int NumPar = 2;
    double params[] = {SM.getMz(), mf(SM.getLeptons(SM.MU))};

    if (CacheCheck(logMZtoMMU_cache, NumPar, params))
        return logMZtoMMU_cache[NumPar];
    else {
        double newResult = log(SM.getMz() / mf(SM.getLeptons(SM.MU)));
        newCacheForDouble(logMZtoMMU_cache, NumPar, params, newResult);
        return newResult;
    }
}

double EWSMcache::logMZtoMTAU() const
{
    int NumPar = 2;
    double params[] = {SM.getMz(), mf(SM.getLeptons(SM.TAU))};

    if (CacheCheck(logMZtoMTAU_cache, NumPar, params))
        return logMZtoMTAU_cache[NumPar];
    else {
        double newResult = log(SM.getMz() / mf(SM.getLeptons(SM.TAU)));
        newCacheForDouble(logMZtoMTAU_cache, NumPar, params, newResult);
        return newResult;
    }
}

double EWSMcache::logMZtoMTOP() const
{
    int NumPar = 2;
    double params[] = {SM.getMz(), SM.getMtpole()};

    if (CacheCheck(logMZtoMTOP_cache, NumPar, params))
        return logMZtoMTOP_cache[NumPar];
    else {
        double newResult = log(SM.getMz() / SM.getMtpole());
        newCacheForDouble(logMZtoMTOP_cache, NumPar, params, newResult);
        return newResult;
    }
}

double EWSMcache::logMTOPtoMH() const
{
    int NumPar = 2;
    double params[] = {SM.getMtpole(), SM.getMHl()};

    if (CacheCheck(logMTOPtoMH_cache, NumPar, params))
        return logMTOPtoMH_cache[NumPar];
    else {
        double newResult = log(SM.getMtpole() / SM.getMHl());
        newCacheForDouble(logMTOPtoMH_cache, NumPar, params, newResult);
        return newResult;
    }
}

double EWSMcache::log_cW2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {SM.getMz(), Mw_i};

    if (CacheCheck(log_cW2_cache, NumPar, params))
        return log_cW2_cache[NumPar];
    else {
        double newResult = log(SM.cW2(Mw_i));
        newCacheForDouble(log_cW2_cache, NumPar, params, newResult);
        return newResult;
    }
}

double EWSMcache::Li2_MW2toMTOP2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mw_i, SM.getMtpole()};

    if (CacheCheck(Li2_MW2toMTOP2_cache, NumPar, params))
        return Li2_MW2toMTOP2_cache[NumPar];
    else {
        double newResult = PolyLog.Li2(Mw_i * Mw_i / SM.getMtpole() / SM.getMtpole()).real();
        newCacheForDouble(Li2_MW2toMTOP2_cache, NumPar, params, newResult);
        return newResult;
    }
}

double EWSMcache::Li3_MW2toMTOP2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mw_i, SM.getMtpole()};

    if (CacheCheck(Li3_MW2toMTOP2_cache, NumPar, params))
        return Li3_MW2toMTOP2_cache[NumPar];
    else {
        double newResult = PolyLog.Li3(Mw_i * Mw_i / SM.getMtpole() / SM.getMtpole());
        newCacheForDouble(Li3_MW2toMTOP2_cache, NumPar, params, newResult);
        return newResult;
    }
}

double EWSMcache::Li3_for_F1(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mw_i, SM.getMtpole()};

    if (CacheCheck(Li3_for_F1_cache, NumPar, params))
        return Li3_for_F1_cache[NumPar];
    else {
        double tmp = Mw_i * Mw_i / SM.getMtpole() / SM.getMtpole();
        double newResult = PolyLog.Li3(-tmp / (1.0 - tmp));
        newCacheForDouble(Li3_for_F1_cache, NumPar, params, newResult);
        return newResult;
    }
}

double EWSMcache::A0_Mz2_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {SM.getMz(), Mw_i};

    if (CacheCheck(A0_Mz2_Mw2_cache, NumPar, params))
        return A0_Mz2_Mw2_cache[NumPar];
    else {
        double newResult = PV.A0(SM.getMz() * SM.getMz(), Mw_i * Mw_i);
        newCacheForDouble(A0_Mz2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}

double EWSMcache::A0_Mz2_mh2() const
{
    int NumPar = 2;
    double params[] = {SM.getMz(), SM.getMHl()};

    if (CacheCheck(A0_Mz2_mh2_cache, NumPar, params))
        return A0_Mz2_mh2_cache[NumPar];
    else {
        double newResult = PV.A0(SM.getMz() * SM.getMz(), SM.getMHl() * SM.getMHl());
        newCacheForDouble(A0_Mz2_mh2_cache, NumPar, params, newResult);
        return newResult;
    }
}

double EWSMcache::A0_Mw2_Mz2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mw_i, SM.getMz()};

    if (CacheCheck(A0_Mw2_Mz2_cache, NumPar, params))
        return A0_Mw2_Mz2_cache[NumPar];
    else {
        double newResult = PV.A0(Mw_i*Mw_i, SM.getMz() * SM.getMz());
        newCacheForDouble(A0_Mw2_Mz2_cache, NumPar, params, newResult);
        return newResult;
    }
}

double EWSMcache::A0_Mw2_mh2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mw_i, SM.getMHl()};

    if (CacheCheck(A0_Mw2_mh2_cache, NumPar, params))
        return A0_Mw2_mh2_cache[NumPar];
    else {
        double newResult = PV.A0(Mw_i*Mw_i, SM.getMHl() * SM.getMHl());
        newCacheForDouble(A0_Mw2_mh2_cache, NumPar, params, newResult);
        return newResult;
    }
}

double EWSMcache::A0_Mz2_Mz2() const
{
    int NumPar = 1;
    double params[] = {SM.getMz()};

    if (CacheCheck(A0_Mz2_Mz2_cache, NumPar, params))
        return A0_Mz2_Mz2_cache[NumPar];
    else {
        double newResult = PV.A0(SM.getMz() * SM.getMz(), SM.getMz() * SM.getMz());
        newCacheForDouble(A0_Mz2_Mz2_cache, NumPar, params, newResult);
        return newResult;
    }
}

double EWSMcache::A0_Mw2_Mw2(const double Mw_i) const
{
    int NumPar = 1;
    double params[] = {Mw_i};

    if (CacheCheck(A0_Mw2_Mw2_cache, NumPar, params))
        return A0_Mw2_Mw2_cache[NumPar];
    else {
        double newResult = PV.A0(Mw_i*Mw_i, Mw_i * Mw_i);
        newCacheForDouble(A0_Mw2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B0_Mz2_Mw2_mh2_Mw2(const double Mw_i) const
{
    int NumPar = 3;
    double params[] = {SM.getMz(), Mw_i, SM.getMHl()};

    if (CacheCheck(B0_Mz2_Mw2_mh2_Mw2_cache, NumPar, params))
        return gslpp::complex(B0_Mz2_Mw2_mh2_Mw2_cache[NumPar],
            B0_Mz2_Mw2_mh2_Mw2_cache[NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.B0(SM.getMz() * SM.getMz(), Mw_i*Mw_i, SM.getMHl() * SM.getMHl(), Mw_i * Mw_i);
        newCacheForComplex(B0_Mz2_Mw2_mh2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B0_Mz2_0_mh2_Mw2(const double Mw_i) const
{
    int NumPar = 3;
    double params[] = {SM.getMz(), SM.getMHl(), Mw_i};

    if (CacheCheck(B0_Mz2_0_mh2_Mw2_cache, NumPar, params))
        return gslpp::complex(B0_Mz2_0_mh2_Mw2_cache[NumPar],
            B0_Mz2_0_mh2_Mw2_cache[NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.B0(SM.getMz() * SM.getMz(), 0.0, SM.getMHl() * SM.getMHl(), Mw_i * Mw_i);
        newCacheForComplex(B0_Mz2_0_mh2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B0_Mw2_Mz2_Mt2_Mt2(const double Mw_i) const
{
    int NumPar = 3;
    double params[] = {Mw_i, SM.getMz(), SM.getMtpole()};

    if (CacheCheck(B0_Mw2_Mz2_Mt2_Mt2_cache, NumPar, params))
        return gslpp::complex(B0_Mw2_Mz2_Mt2_Mt2_cache[NumPar],
            B0_Mw2_Mz2_Mt2_Mt2_cache[NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.B0(Mw_i*Mw_i, SM.getMz() * SM.getMz(), SM.getMtpole() * SM.getMtpole(), SM.getMtpole() * SM.getMtpole());
        newCacheForComplex(B0_Mw2_Mz2_Mt2_Mt2_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B0_Mz2_Mz2_Mw2_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {SM.getMz(), Mw_i};

    if (CacheCheck(B0_Mz2_Mz2_Mw2_Mw2_cache, NumPar, params))
        return gslpp::complex(B0_Mz2_Mz2_Mw2_Mw2_cache[NumPar],
            B0_Mz2_Mz2_Mw2_Mw2_cache[NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.B0(SM.getMz() * SM.getMz(), SM.getMz() * SM.getMz(), Mw_i*Mw_i, Mw_i * Mw_i);
        newCacheForComplex(B0_Mz2_Mz2_Mw2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B0_Mz2_Mz2_mh2_Mz2() const
{
    int NumPar = 2;
    double params[] = {SM.getMz(), SM.getMHl()};

    if (CacheCheck(B0_Mz2_Mz2_mh2_Mz2_cache, NumPar, params))
        return gslpp::complex(B0_Mz2_Mz2_mh2_Mz2_cache[NumPar],
            B0_Mz2_Mz2_mh2_Mz2_cache[NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.B0(SM.getMz() * SM.getMz(), SM.getMz() * SM.getMz(), SM.getMHl() * SM.getMHl(), SM.getMz() * SM.getMz());
        newCacheForComplex(B0_Mz2_Mz2_mh2_Mz2_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B0_Mz2_Mw2_Mz2_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {SM.getMz(), Mw_i};

    if (CacheCheck(B0_Mz2_Mw2_Mz2_Mw2_cache, NumPar, params))
        return gslpp::complex(B0_Mz2_Mw2_Mz2_Mw2_cache[NumPar],
            B0_Mz2_Mw2_Mz2_Mw2_cache[NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.B0(SM.getMz() * SM.getMz(), Mw_i*Mw_i, SM.getMz() * SM.getMz(), Mw_i * Mw_i);
        newCacheForComplex(B0_Mz2_Mw2_Mz2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B0_Mz2_Mw2_0_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {SM.getMz(), Mw_i};

    if (CacheCheck(B0_Mz2_Mw2_0_Mw2_cache, NumPar, params))
        return gslpp::complex(B0_Mz2_Mw2_0_Mw2_cache[NumPar],
            B0_Mz2_Mw2_0_Mw2_cache[NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.B0(SM.getMz() * SM.getMz(), Mw_i*Mw_i, 0.0, Mw_i * Mw_i);
        newCacheForComplex(B0_Mz2_Mw2_0_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B0_Mz2_0_Mz2_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {SM.getMz(), Mw_i};

    if (CacheCheck(B0_Mz2_0_Mz2_Mw2_cache, NumPar, params))
        return gslpp::complex(B0_Mz2_0_Mz2_Mw2_cache[NumPar],
            B0_Mz2_0_Mz2_Mw2_cache[NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.B0(SM.getMz() * SM.getMz(), 0.0, SM.getMz() * SM.getMz(), Mw_i * Mw_i);
        newCacheForComplex(B0_Mz2_0_Mz2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B0_Mz2_0_0_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {SM.getMz(), Mw_i};

    if (CacheCheck(B0_Mz2_0_0_Mw2_cache, NumPar, params))
        return gslpp::complex(B0_Mz2_0_0_Mw2_cache[NumPar],
            B0_Mz2_0_0_Mw2_cache[NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.B0(SM.getMz() * SM.getMz(), 0.0, 0.0, Mw_i * Mw_i);
        newCacheForComplex(B0_Mz2_0_0_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B0_Mw2_Mz2_Mw2_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mw_i, SM.getMz()};

    if (CacheCheck(B0_Mw2_Mz2_Mw2_Mw2_cache, NumPar, params))
        return gslpp::complex(B0_Mw2_Mz2_Mw2_Mw2_cache[NumPar],
            B0_Mw2_Mz2_Mw2_Mw2_cache[NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.B0(Mw_i*Mw_i, SM.getMz() * SM.getMz(), Mw_i*Mw_i, Mw_i * Mw_i);
        newCacheForComplex(B0_Mw2_Mz2_Mw2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B0_Mw2_Mw2_Mz2_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mw_i, SM.getMz()};

    if (CacheCheck(B0_Mw2_Mw2_Mz2_Mw2_cache, NumPar, params))
        return gslpp::complex(B0_Mw2_Mw2_Mz2_Mw2_cache[NumPar],
            B0_Mw2_Mw2_Mz2_Mw2_cache[NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.B0(Mw_i*Mw_i, Mw_i*Mw_i, SM.getMz() * SM.getMz(), Mw_i * Mw_i);
        newCacheForComplex(B0_Mw2_Mw2_Mz2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B0_Mw2_Mw2_mh2_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mw_i, SM.getMHl()};

    if (CacheCheck(B0_Mw2_Mw2_mh2_Mw2_cache, NumPar, params))
        return gslpp::complex(B0_Mw2_Mw2_mh2_Mw2_cache[NumPar],
            B0_Mw2_Mw2_mh2_Mw2_cache[NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.B0(Mw_i*Mw_i, Mw_i*Mw_i, SM.getMHl() * SM.getMHl(), Mw_i * Mw_i);
        newCacheForComplex(B0_Mw2_Mw2_mh2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B0_Mw2_Mw2_0_Mw2(const double Mw_i) const
{
    int NumPar = 1;
    double params[] = {Mw_i};

    if (CacheCheck(B0_Mw2_Mw2_0_Mw2_cache, NumPar, params))
        return gslpp::complex(B0_Mw2_Mw2_0_Mw2_cache[NumPar],
            B0_Mw2_Mw2_0_Mw2_cache[NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.B0(Mw_i*Mw_i, Mw_i*Mw_i, 0.0, Mw_i * Mw_i);
        newCacheForComplex(B0_Mw2_Mw2_0_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B0_Mz2_Mz2_mf2_mf2(const Particle f) const
{
    int NumPar = 2;
    double params[] = {SM.getMz(), mf(f, SM.getMz())};
    int ind = f.getIndex();
    if (CacheCheck(B0_Mz2_Mz2_mf2_mf2_cache[ind], NumPar, params))
        return gslpp::complex(B0_Mz2_Mz2_mf2_mf2_cache[ind][NumPar],
            B0_Mz2_Mz2_mf2_mf2_cache[ind][NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.B0(SM.getMz() * SM.getMz(), SM.getMz() * SM.getMz(), mf2(f, SM.getMz()), mf2(f, SM.getMz()));
        newCacheForComplex(B0_Mz2_Mz2_mf2_mf2_cache[ind], NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B0p_Mz2_0_mh2_Mw2(const double Mw_i) const
{
    int NumPar = 3;
    double params[] = {SM.getMz(), SM.getMHl(), Mw_i};

    if (CacheCheck(B0p_Mz2_0_mh2_Mw2_cache, NumPar, params))
        return gslpp::complex(B0p_Mz2_0_mh2_Mw2_cache[NumPar],
            B0p_Mz2_0_mh2_Mw2_cache[NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.B0p(SM.getMz() * SM.getMz(), 0.0, SM.getMHl() * SM.getMHl(), Mw_i * Mw_i);
        newCacheForComplex(B0p_Mz2_0_mh2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B0p_Mz2_Mz2_mh2_Mz2() const
{
    int NumPar = 2;
    double params[] = {SM.getMz(), SM.getMHl()};

    if (CacheCheck(B0p_Mz2_Mz2_mh2_Mz2_cache, NumPar, params))
        return gslpp::complex(B0p_Mz2_Mz2_mh2_Mz2_cache[NumPar],
            B0p_Mz2_Mz2_mh2_Mz2_cache[NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.B0p(SM.getMz() * SM.getMz(), SM.getMz() * SM.getMz(), SM.getMHl() * SM.getMHl(), SM.getMz() * SM.getMz());
        newCacheForComplex(B0p_Mz2_Mz2_mh2_Mz2_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B0p_Mz2_0_Mz2_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {SM.getMz(), Mw_i};

    if (CacheCheck(B0p_Mz2_0_Mz2_Mw2_cache, NumPar, params))
        return gslpp::complex(B0p_Mz2_0_Mz2_Mw2_cache[NumPar],
            B0p_Mz2_0_Mz2_Mw2_cache[NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.B0p(SM.getMz() * SM.getMz(), 0.0, SM.getMz() * SM.getMz(), Mw_i * Mw_i);
        newCacheForComplex(B0p_Mz2_0_Mz2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B0p_Mz2_Mz2_Mw2_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {SM.getMz(), Mw_i};

    if (CacheCheck(B0p_Mz2_Mz2_Mw2_Mw2_cache, NumPar, params))
        return gslpp::complex(B0p_Mz2_Mz2_Mw2_Mw2_cache[NumPar],
            B0p_Mz2_Mz2_Mw2_Mw2_cache[NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.B0p(SM.getMz() * SM.getMz(), SM.getMz() * SM.getMz(), Mw_i*Mw_i, Mw_i * Mw_i);
        newCacheForComplex(B0p_Mz2_Mz2_Mw2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B0p_Mw2_Mw2_Mz2_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mw_i, SM.getMz()};

    if (CacheCheck(B0p_Mw2_Mw2_Mz2_Mw2_cache, NumPar, params))
        return gslpp::complex(B0p_Mw2_Mw2_Mz2_Mw2_cache[NumPar],
            B0p_Mw2_Mw2_Mz2_Mw2_cache[NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.B0p(Mw_i*Mw_i, Mw_i*Mw_i, SM.getMz() * SM.getMz(), Mw_i * Mw_i);
        newCacheForComplex(B0p_Mw2_Mw2_Mz2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B0p_Mw2_Mw2_mh2_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mw_i, SM.getMHl()};

    if (CacheCheck(B0p_Mw2_Mw2_mh2_Mw2_cache, NumPar, params))
        return gslpp::complex(B0p_Mw2_Mw2_mh2_Mw2_cache[NumPar],
            B0p_Mw2_Mw2_mh2_Mw2_cache[NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.B0p(Mw_i*Mw_i, Mw_i*Mw_i, SM.getMHl() * SM.getMHl(), Mw_i * Mw_i);
        newCacheForComplex(B0p_Mw2_Mw2_mh2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B0p_Mw2_Mw2_0_Mw2(const double Mw_i) const
{
    int NumPar = 1;
    double params[] = {Mw_i};

    if (CacheCheck(B0p_Mw2_Mw2_0_Mw2_cache, NumPar, params))
        return gslpp::complex(B0p_Mw2_Mw2_0_Mw2_cache[NumPar],
            B0p_Mw2_Mw2_0_Mw2_cache[NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.B0p(Mw_i*Mw_i, Mw_i*Mw_i, 0.0, Mw_i * Mw_i);
        newCacheForComplex(B0p_Mw2_Mw2_0_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B0p_Mz2_Mz2_mf2_mf2(const Particle f) const
{
    int NumPar = 2;
    double params[] = {SM.getMz(), mf(f, SM.getMz())};
    int ind = f.getIndex();
    if (CacheCheck(B0p_Mz2_Mz2_mf2_mf2_cache[ind], NumPar, params))
        return gslpp::complex(B0p_Mz2_Mz2_mf2_mf2_cache[ind][NumPar],
            B0p_Mz2_Mz2_mf2_mf2_cache[ind][NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.B0p(SM.getMz() * SM.getMz(), SM.getMz() * SM.getMz(), mf2(f, SM.getMz()), mf2(f, SM.getMz()));
        newCacheForComplex(B0p_Mz2_Mz2_mf2_mf2_cache[ind], NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B1_Mz2_0_mf2_mfprime2(const int gen) const
{
    int NumPar = 3;
    double mymf, mymfprime;
    double Mz = SM.getMz();
    if (gen < 3) {
        mymf = mf(SM.getLeptons((StandardModel::lepton) (2 * gen)), Mz);
        mymfprime = mf(SM.getLeptons((StandardModel::lepton) (2 * gen + 1)), Mz);
    } else {
        int genq = gen - 3;
        mymf = mf(SM.getQuarks((QCD::quark) (2 * genq)), Mz);
        mymfprime = mf(SM.getQuarks((QCD::quark) (2 * genq + 1)), Mz);
    }
    double params[] = {Mz, mymf, mymfprime};

    if (CacheCheck(B1_Mz2_0_mf2_mfprime2_cache[gen], NumPar, params))
        return gslpp::complex(B1_Mz2_0_mf2_mfprime2_cache[gen][NumPar],
            B1_Mz2_0_mf2_mfprime2_cache[gen][NumPar + 1], false);
    else {
        double mf2 = mymf*mymf;
        double mfprime2 = mymfprime*mymfprime;
        gslpp::complex newResult = PV.B1(Mz * Mz, 0.0, mf2, mfprime2);
        newCacheForComplex(B1_Mz2_0_mf2_mfprime2_cache[gen], NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B1_Mz2_0_mfprime2_mf2(const int gen) const
{
    int NumPar = 3;
    double mymf, mymfprime;
    double Mz = SM.getMz();
    if (gen < 3) {
        mymf = mf(SM.getLeptons((StandardModel::lepton) (2 * gen)), Mz);
        mymfprime = mf(SM.getLeptons((StandardModel::lepton) (2 * gen + 1)), Mz);
    } else {
        int genq = gen - 3;
        mymf = mf(SM.getQuarks((QCD::quark) (2 * genq)), Mz);
        mymfprime = mf(SM.getQuarks((QCD::quark) (2 * genq + 1)), Mz);
    }
    double params[] = {Mz, mymf, mymfprime};

    if (CacheCheck(B1_Mz2_0_mfprime2_mf2_cache[gen], NumPar, params))
        return gslpp::complex(B1_Mz2_0_mfprime2_mf2_cache[gen][NumPar],
            B1_Mz2_0_mfprime2_mf2_cache[gen][NumPar + 1], false);
    else {
        double mf2 = mymf*mymf;
        double mfprime2 = mymfprime*mymfprime;
        gslpp::complex newResult = PV.B1(Mz * Mz, 0.0, mfprime2, mf2);
        newCacheForComplex(B1_Mz2_0_mfprime2_mf2_cache[gen], NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B1_Mz2_Mw2_mf2_mfprime2(const int gen, const double Mw_i) const
{
    int NumPar = 4;
    double mymf, mymfprime;
    double Mz = SM.getMz();
    if (gen < 3) {
        mymf = mf(SM.getLeptons((StandardModel::lepton) (2 * gen)), Mz);
        mymfprime = mf(SM.getLeptons((StandardModel::lepton) (2 * gen + 1)), Mz);
    } else {
        int genq = gen - 3;
        mymf = mf(SM.getQuarks((QCD::quark) (2 * genq)), Mz);
        mymfprime = mf(SM.getQuarks((QCD::quark) (2 * genq + 1)), Mz);
    }
    double params[] = {Mz, Mw_i, mymf, mymfprime};

    if (CacheCheck(B1_Mz2_Mw2_mf2_mfprime2_cache[gen], NumPar, params))
        return gslpp::complex(B1_Mz2_Mw2_mf2_mfprime2_cache[gen][NumPar],
            B1_Mz2_Mw2_mf2_mfprime2_cache[gen][NumPar + 1], false);
    else {
        double mf2 = mymf*mymf;
        double mfprime2 = mymfprime*mymfprime;
        gslpp::complex newResult = PV.B1(Mz * Mz, Mw_i*Mw_i, mf2, mfprime2);
        newCacheForComplex(B1_Mz2_Mw2_mf2_mfprime2_cache[gen], NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B1_Mz2_Mw2_mfprime2_mf2(const int gen, const double Mw_i) const
{
    int NumPar = 4;
    double mymf, mymfprime;
    double Mz = SM.getMz();
    if (gen < 3) {
        mymf = mf(SM.getLeptons((StandardModel::lepton) (2 * gen)), Mz);
        mymfprime = mf(SM.getLeptons((StandardModel::lepton) (2 * gen + 1)), Mz);
    } else {
        int genq = gen - 3;
        mymf = mf(SM.getQuarks((QCD::quark) (2 * genq)), Mz);
        mymfprime = mf(SM.getQuarks((QCD::quark) (2 * genq + 1)), Mz);
    }
    double params[] = {Mz, Mw_i, mymf, mymfprime};

    if (CacheCheck(B1_Mz2_Mw2_mfprime2_mf2_cache[gen], NumPar, params))
        return gslpp::complex(B1_Mz2_Mw2_mfprime2_mf2_cache[gen][NumPar],
            B1_Mz2_Mw2_mfprime2_mf2_cache[gen][NumPar + 1], false);
    else {
        double mf2 = mymf*mymf;
        double mfprime2 = mymfprime*mymfprime;
        gslpp::complex newResult = PV.B1(Mz * Mz, Mw_i*Mw_i, mfprime2, mf2);
        newCacheForComplex(B1_Mz2_Mw2_mfprime2_mf2_cache[gen], NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B1p_Mw2_Mw2_mf2_mfprime2(const int gen, const double Mw_i) const
{
    int NumPar = 3;
    double mymf, mymfprime;
    double Mz = SM.getMz();
    if (gen < 3) {
        mymf = mf(SM.getLeptons((StandardModel::lepton) (2 * gen)), Mz);
        mymfprime = mf(SM.getLeptons((StandardModel::lepton) (2 * gen + 1)), Mz);
    } else {
        int genq = gen - 3;
        mymf = mf(SM.getQuarks((QCD::quark) (2 * genq)), Mz);
        mymfprime = mf(SM.getQuarks((QCD::quark) (2 * genq + 1)), Mz);
    }
    double params[] = {Mw_i, mymf, mymfprime};

    if (CacheCheck(B1p_Mw2_Mw2_mf2_mfprime2_cache[gen], NumPar, params))
        return gslpp::complex(B1p_Mw2_Mw2_mf2_mfprime2_cache[gen][NumPar],
            B1p_Mw2_Mw2_mf2_mfprime2_cache[gen][NumPar + 1], false);
    else {
        double mf2 = mymf*mymf;
        double mfprime2 = mymfprime*mymfprime;
        gslpp::complex newResult = PV.B1p(Mw_i*Mw_i, Mw_i*Mw_i, mf2, mfprime2);
        newCacheForComplex(B1p_Mw2_Mw2_mf2_mfprime2_cache[gen], NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::B1p_Mw2_Mw2_mfprime2_mf2(const int gen, const double Mw_i) const
{
    int NumPar = 3;
    double mymf, mymfprime;
    double Mz = SM.getMz();
    if (gen < 3) {
        mymf = mf(SM.getLeptons((StandardModel::lepton) (2 * gen)), Mz);
        mymfprime = mf(SM.getLeptons((StandardModel::lepton) (2 * gen + 1)), Mz);
    } else {
        int genq = gen - 3;
        mymf = mf(SM.getQuarks((QCD::quark) (2 * genq)), Mz);
        mymfprime = mf(SM.getQuarks((QCD::quark) (2 * genq + 1)), Mz);
    }
    double params[] = {Mw_i, mymf, mymfprime};

    if (CacheCheck(B1p_Mw2_Mw2_mfprime2_mf2_cache[gen], NumPar, params))
        return gslpp::complex(B1p_Mw2_Mw2_mfprime2_mf2_cache[gen][NumPar],
            B1p_Mw2_Mw2_mfprime2_mf2_cache[gen][NumPar + 1], false);
    else {
        double mf2 = mymf*mymf;
        double mfprime2 = mymfprime*mymfprime;
        gslpp::complex newResult = PV.B1p(Mw_i*Mw_i, Mw_i*Mw_i, mfprime2, mf2);
        newCacheForComplex(B1p_Mw2_Mw2_mfprime2_mf2_cache[gen], NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::Bf_Mz2_Mz2_mf2_mf2(const Particle f) const
{
    int NumPar = 2;
    double params[] = {SM.getMz(), mf(f, SM.getMz())};
    int ind = f.getIndex();
    if (CacheCheck(Bf_Mz2_Mz2_mf2_mf2_cache[ind], NumPar, params))
        return gslpp::complex(Bf_Mz2_Mz2_mf2_mf2_cache[ind][NumPar],
            Bf_Mz2_Mz2_mf2_mf2_cache[ind][NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.Bf(SM.getMz() * SM.getMz(), SM.getMz() * SM.getMz(), mf2(f, SM.getMz()), mf2(f, SM.getMz()));
        newCacheForComplex(Bf_Mz2_Mz2_mf2_mf2_cache[ind], NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::Bf_Mz2_0_mf2_mf2(const Particle f) const
{
    int NumPar = 2;
    double params[] = {SM.getMz(), mf(f, SM.getMz())};
    if (params[1] == 0.0)
        throw std::runtime_error("Error in EWSMcache::Bf_Mz_0_mf_mf()");
    int ind = f.getIndex();
    if (CacheCheck(Bf_Mz2_0_mf2_mf2_cache[ind], NumPar, params))
        return gslpp::complex(Bf_Mz2_0_mf2_mf2_cache[ind][NumPar],
            Bf_Mz2_0_mf2_mf2_cache[ind][NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.Bf(SM.getMz() * SM.getMz(), 0.0, mf2(f, SM.getMz()), mf2(f, SM.getMz()));
        newCacheForComplex(Bf_Mz2_0_mf2_mf2_cache[ind], NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::Bf_Mz2_Mw2_mfprime2_mf2(const int gen, const double Mw_i) const
{
    int NumPar = 4;
    double mymf, mymfprime;
    double Mz = SM.getMz();
    if (gen < 3) {
        mymf = mf(SM.getLeptons((StandardModel::lepton) (2 * gen)), Mz);
        mymfprime = mf(SM.getLeptons((StandardModel::lepton) (2 * gen + 1)), Mz);
    } else {
        int genq = gen - 3;
        mymf = mf(SM.getQuarks((QCD::quark) (2 * genq)), Mz);
        mymfprime = mf(SM.getQuarks((QCD::quark) (2 * genq + 1)), Mz);
    }
    double params[] = {Mz, Mw_i, mymf, mymfprime};
    if (CacheCheck(Bf_Mz2_Mw2_mfprime2_mf2_cache[gen], NumPar, params))
        return gslpp::complex(Bf_Mz2_Mw2_mfprime2_mf2_cache[gen][NumPar],
            Bf_Mz2_Mw2_mfprime2_mf2_cache[gen][NumPar + 1], false);
    else {
        double mf2 = mymf*mymf;
        double mfprime2 = mymfprime*mymfprime;
        gslpp::complex newResult = PV.Bf(Mz * Mz, Mw_i*Mw_i, mfprime2, mf2);
        newCacheForComplex(Bf_Mz2_Mw2_mfprime2_mf2_cache[gen], NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::Bf_Mz2_0_mfprime2_mf2(const int gen) const
{
    int NumPar = 3;
    double mymf, mymfprime;
    double Mz = SM.getMz();
    if (gen < 3) {
        mymf = mf(SM.getLeptons((StandardModel::lepton) (2 * gen)), Mz);
        mymfprime = mf(SM.getLeptons((StandardModel::lepton) (2 * gen + 1)), Mz);
    } else {
        int genq = gen - 3;
        mymf = mf(SM.getQuarks((QCD::quark) (2 * genq)), Mz);
        mymfprime = mf(SM.getQuarks((QCD::quark) (2 * genq + 1)), Mz);
    }
    double params[] = {Mz, mymf, mymfprime};

    if (CacheCheck(Bf_Mz2_0_mfprime2_mf2_cache[gen], NumPar, params))
        return gslpp::complex(Bf_Mz2_0_mfprime2_mf2_cache[gen][NumPar],
            Bf_Mz2_0_mfprime2_mf2_cache[gen][NumPar + 1], false);
    else {
        double mf2 = mymf*mymf;
        double mfprime2 = mymfprime*mymfprime;
        gslpp::complex newResult = PV.Bf(Mz * Mz, 0.0, mfprime2, mf2);
        newCacheForComplex(Bf_Mz2_0_mfprime2_mf2_cache[gen], NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::Bf_Mw2_Mw2_mfprime2_mf2(const int gen, const double Mw_i) const
{
    int NumPar = 3;
    double mymf, mymfprime;
    double Mz = SM.getMz();
    if (gen < 3) {
        mymf = mf(SM.getLeptons((StandardModel::lepton) (2 * gen)), Mz);
        mymfprime = mf(SM.getLeptons((StandardModel::lepton) (2 * gen + 1)), Mz);
    } else {
        int genq = gen - 3;
        mymf = mf(SM.getQuarks((QCD::quark) (2 * genq)), Mz);
        mymfprime = mf(SM.getQuarks((QCD::quark) (2 * genq + 1)), Mz);
    }
    double params[] = {Mw_i, mymf, mymfprime};

    if (CacheCheck(Bf_Mw2_Mw2_mfprime2_mf2_cache[gen], NumPar, params))
        return gslpp::complex(Bf_Mw2_Mw2_mfprime2_mf2_cache[gen][NumPar],
            Bf_Mw2_Mw2_mfprime2_mf2_cache[gen][NumPar + 1], false);
    else {
        double mf2 = mymf*mymf;
        double mfprime2 = mymfprime*mymfprime;
        gslpp::complex newResult = PV.Bf(Mw_i*Mw_i, Mw_i*Mw_i, mfprime2, mf2);
        newCacheForComplex(Bf_Mw2_Mw2_mfprime2_mf2_cache[gen], NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::Bfp_Mz2_Mz2_mf2_mf2(const Particle f) const
{
    int NumPar = 2;
    double params[] = {SM.getMz(), mf(f, SM.getMz())};
    int ind = f.getIndex();
    if (CacheCheck(Bfp_Mz2_Mz2_mf2_mf2_cache[ind], NumPar, params))
        return gslpp::complex(Bfp_Mz2_Mz2_mf2_mf2_cache[ind][NumPar],
            Bfp_Mz2_Mz2_mf2_mf2_cache[ind][NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.Bfp(SM.getMz() * SM.getMz(), SM.getMz() * SM.getMz(), mf2(f, SM.getMz()), mf2(f, SM.getMz()));
        newCacheForComplex(Bfp_Mz2_Mz2_mf2_mf2_cache[ind], NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::Bfp_Mw2_Mw2_mfprime2_mf2(const int gen, const double Mw_i) const
{
    int NumPar = 3;
    double mymf, mymfprime;
    double Mz = SM.getMz();
    if (gen < 3) {
        mymf = mf(SM.getLeptons((StandardModel::lepton) (2 * gen)), Mz);
        mymfprime = mf(SM.getLeptons((StandardModel::lepton) (2 * gen + 1)), Mz);
    } else {
        int genq = gen - 3;
        mymf = mf(SM.getQuarks((QCD::quark) (2 * genq)), Mz);
        mymfprime = mf(SM.getQuarks((QCD::quark) (2 * genq + 1)), Mz);
    }
    double params[] = {Mw_i, mymf, mymfprime};

    if (CacheCheck(Bfp_Mw2_Mw2_mfprime2_mf2_cache[gen], NumPar, params))
        return gslpp::complex(Bfp_Mw2_Mw2_mfprime2_mf2_cache[gen][NumPar],
            Bfp_Mw2_Mw2_mfprime2_mf2_cache[gen][NumPar + 1], false);
    else {
        double mf2 = mymf*mymf;
        double mfprime2 = mymfprime*mymfprime;
        gslpp::complex newResult = PV.Bfp(Mw_i*Mw_i, Mw_i*Mw_i, mfprime2, mf2);
        newCacheForComplex(Bfp_Mw2_Mw2_mfprime2_mf2_cache[gen], NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::C0_Mz2_Mw2_Mt2_Mw2(const double Mw_i) const
{
    int NumPar = 3;
    double params[] = {SM.getMz(), Mw_i, SM.getMtpole()};

    if (CacheCheck(C0_Mz2_Mw2_Mt2_Mw2_cache, NumPar, params))
        return gslpp::complex(C0_Mz2_Mw2_Mt2_Mw2_cache[NumPar],
            C0_Mz2_Mw2_Mt2_Mw2_cache[NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.C0(SM.getMz() * SM.getMz(), Mw_i*Mw_i, SM.getMtpole() * SM.getMtpole(), Mw_i * Mw_i);
        newCacheForComplex(C0_Mz2_Mw2_Mt2_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::C0_Mz2_Mt2_Mw2_Mt2(const double Mw_i) const
{
    int NumPar = 3;
    double params[] = {SM.getMz(), SM.getMtpole(), Mw_i};

    if (CacheCheck(C0_Mz2_Mt2_Mw2_Mt2_cache, NumPar, params))
        return gslpp::complex(C0_Mz2_Mt2_Mw2_Mt2_cache[NumPar],
            C0_Mz2_Mt2_Mw2_Mt2_cache[NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.C0(SM.getMz() * SM.getMz(), SM.getMtpole() * SM.getMtpole(), Mw_i*Mw_i, SM.getMtpole() * SM.getMtpole());
        newCacheForComplex(C0_Mz2_Mt2_Mw2_Mt2_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::C0_Mz2_0_Mw2_0(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {SM.getMz(), Mw_i};

    if (CacheCheck(C0_Mz2_0_Mw2_0_cache, NumPar, params))
        return gslpp::complex(C0_Mz2_0_Mw2_0_cache[NumPar],
            C0_Mz2_0_Mw2_0_cache[NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.C0(SM.getMz() * SM.getMz(), 0.0, Mw_i*Mw_i, 0.0);
        newCacheForComplex(C0_Mz2_0_Mw2_0_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::C0_Mz2_Mw2_0_Mw2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {SM.getMz(), Mw_i};

    if (CacheCheck(C0_Mz2_Mw2_0_Mw2_cache, NumPar, params))
        return gslpp::complex(C0_Mz2_Mw2_0_Mw2_cache[NumPar],
            C0_Mz2_Mw2_0_Mw2_cache[NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.C0(SM.getMz() * SM.getMz(), Mw_i*Mw_i, 0.0, Mw_i * Mw_i);
        newCacheForComplex(C0_Mz2_Mw2_0_Mw2_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::C0_Mw2_Mw2_0_Mz2(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mw_i, SM.getMz()};

    if (CacheCheck(C0_Mw2_Mw2_0_Mz2_cache, NumPar, params))
        return gslpp::complex(C0_Mw2_Mw2_0_Mz2_cache[NumPar],
            C0_Mw2_Mw2_0_Mz2_cache[NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.C0(Mw_i*Mw_i, Mw_i*Mw_i, 0.0, SM.getMz() * SM.getMz());
        newCacheForComplex(C0_Mw2_Mw2_0_Mz2_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::C0_Mw2_0_Mz2_0(const double Mw_i) const
{
    int NumPar = 2;
    double params[] = {Mw_i, SM.getMz()};

    if (CacheCheck(C0_Mw2_0_Mz2_0_cache, NumPar, params))
        return gslpp::complex(C0_Mw2_0_Mz2_0_cache[NumPar],
            C0_Mw2_0_Mz2_0_cache[NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.C0(Mw_i*Mw_i, 0.0, SM.getMz() * SM.getMz(), 0.0);
        newCacheForComplex(C0_Mw2_0_Mz2_0_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex EWSMcache::C0_Mz2_0_Mz2_0() const
{
    int NumPar = 1;
    double params[] = {SM.getMz()};

    if (CacheCheck(C0_Mz2_0_Mz2_0_cache, NumPar, params))
        return gslpp::complex(C0_Mz2_0_Mz2_0_cache[NumPar],
            C0_Mz2_0_Mz2_0_cache[NumPar + 1], false);
    else {
        gslpp::complex newResult = PV.C0(SM.getMz() * SM.getMz(), 0.0, SM.getMz() * SM.getMz(), 0.0);
        newCacheForComplex(C0_Mz2_0_Mz2_0_cache, NumPar, params, newResult);
        return newResult;
    }
}




