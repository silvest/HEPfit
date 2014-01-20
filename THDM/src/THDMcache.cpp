/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "THDMcache.h"

THDMcache::THDMcache()
: PV(true)
{
}


/////////////////////////////////////////////////////////////////////////////////////////////////

int THDMcache::CacheCheck(const complex cache[][CacheSize], 
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

void THDMcache::CacheShift(complex cache[][CacheSize], const int NumPar, 
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

////////////////////////////////////////////////////////////////////////////////
/*One-loop functions*/


complex THDMcache::B0_Mz_Mz2_Mz_mH(const double Mz, const double mH) const {
    int NumPar = 2;
    double params[] = {Mz, mH};

    int i = CacheCheck(B0_Mz_Mz2_Mz_mH_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mz_Mz2_Mz_mH_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0(Mz*Mz, Mz*Mz, Mz*Mz, mH*mH);
        CacheShift(B0_Mz_Mz2_Mz_mH_cache, NumPar, params, newResult);
        return newResult;
    } 
}

complex THDMcache::B0_Mz_0_Mz_mH(const double Mz, const double mH) const {
    int NumPar = 2;
    double params[] = {Mz, mH};

    int i = CacheCheck(B0_Mz_0_Mz_mH_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mz_0_Mz_mH_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0(Mz*Mz, 0., Mz*Mz, mH*mH);
        CacheShift(B0_Mz_0_Mz_mH_cache, NumPar, params, newResult);
        return newResult;
    } 
}

complex THDMcache::B0_Mz_Mz2_Mz_mh(const double Mz, const double mh) const {
    int NumPar = 2;
    double params[] = {Mz, mh};

    int i = CacheCheck(B0_Mz_Mz2_Mz_mh_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mz_Mz2_Mz_mh_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0(Mz*Mz, Mz*Mz, Mz*Mz, mh*mh);
        CacheShift(B0_Mz_Mz2_Mz_mh_cache, NumPar, params, newResult);
        return newResult;
    } 
}

complex THDMcache::B0_Mz_0_Mz_mh(const double Mz, const double mh) const {
    int NumPar = 2;
    double params[] = {Mz, mh};

    int i = CacheCheck(B0_Mz_0_Mz_mh_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mz_0_Mz_mh_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0(Mz*Mz, 0., Mz*Mz, mh*mh);
        CacheShift(B0_Mz_0_Mz_mh_cache, NumPar, params, newResult);
        return newResult;
    } 
}

complex THDMcache::B0_Mz_Mw2_Mw_mH(const double Mz, const double Mw, const double mH) const {
    int NumPar = 3;
    double params[] = {Mz, Mw, mH};

    int i = CacheCheck(B0_Mz_Mw2_Mw_mH_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mz_Mw2_Mw_mH_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0(Mz*Mz, Mw*Mw, Mw*Mz, mH*mH);
        CacheShift(B0_Mz_Mw2_Mw_mH_cache, NumPar, params, newResult);
        return newResult;
    } 
}

complex THDMcache::B0_Mz_0_Mw_mH(const double Mz, const double Mw, const double mH) const {
    int NumPar = 3;
    double params[] = {Mz, Mw, mH};

    int i = CacheCheck(B0_Mz_0_Mw_mH_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mz_0_Mw_mH_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0(Mz*Mz, 0., Mw*Mw, mH*mH);
        CacheShift(B0_Mz_0_Mw_mH_cache, NumPar, params, newResult);
        return newResult;
    } 
}

complex THDMcache::B0_Mz_Mw2_Mw_mh(const double Mz, const double Mw, const double mh) const {
    int NumPar = 3;
    double params[] = {Mz, Mw, mh};

    int i = CacheCheck(B0_Mz_Mw2_Mw_mh_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mz_Mw2_Mw_mh_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0(Mz*Mz, Mw*Mw, Mw*Mw, mh*mh);
        CacheShift(B0_Mz_Mw2_Mw_mh_cache, NumPar, params, newResult);
        return newResult;
    } 
}

complex THDMcache::B0_Mz_0_Mw_mh(const double Mz, const double Mw, const double mh) const {
    int NumPar = 3;
    double params[] = {Mz, Mw, mh};

    int i = CacheCheck(B0_Mz_0_Mw_mh_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mz_0_Mw_mh_cache[NumPar][i] );
    } else {
        complex newResult = PV.B0(Mz*Mz, 0., Mw*Mw, mh*mh);
        CacheShift(B0_Mz_0_Mw_mh_cache, NumPar, params, newResult);
        return newResult;
    } 
}

complex THDMcache::B00_Mz_Mz2_mH_mA(const double Mz, const double mH, const double mA) const {
    int NumPar = 3;
    double params[] = {Mz, mH, mA};

    int i = CacheCheck(B00_Mz_Mz2_mH_mA_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_Mz2_mH_mA_cache[NumPar][i] );
    } else {
        complex newResult = PV.B00(Mz*Mz, Mz*Mz, mH*mH, mA*mA);
        CacheShift(B00_Mz_Mz2_mH_mA_cache, NumPar, params, newResult);
        return newResult;
    } 
}

complex THDMcache::B00_Mz_0_mH_mA(const double Mz, const double mH, const double mA) const {
    int NumPar = 3;
    double params[] = {Mz, mH, mA};

    int i = CacheCheck(B00_Mz_0_mH_mA_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_0_mH_mA_cache[NumPar][i] );
    } else {
        complex newResult = PV.B00(Mz*Mz, 0., mH*mH, mA*mA);
        CacheShift(B00_Mz_0_mH_mA_cache, NumPar, params, newResult);
        return newResult;
    } 
}

complex THDMcache::B00_Mz_Mz2_mHp_mHp(const double Mz, const double mHp) const {
    int NumPar = 2;
    double params[] = {Mz, mHp};

    int i = CacheCheck(B00_Mz_Mz2_mHp_mHp_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_Mz2_mHp_mHp_cache[NumPar][i] );
    } else {
        complex newResult = PV.B00(Mz*Mz, Mz*Mz, mHp*mHp, mHp*mHp);
        CacheShift(B00_Mz_Mz2_mHp_mHp_cache, NumPar, params, newResult);
        return newResult;
    } 
}

complex THDMcache::B00_Mz_0_mHp_mHp(const double Mz, const double mHp) const {
    int NumPar = 2;
    double params[] = {Mz, mHp};

    int i = CacheCheck(B00_Mz_0_mHp_mHp_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_0_mHp_mHp_cache[NumPar][i] );
    } else {
        complex newResult = PV.B00(Mz*Mz, 0., mHp*mHp, mHp*mHp);
        CacheShift(B00_Mz_0_mHp_mHp_cache, NumPar, params, newResult);
        return newResult;
    } 
}

complex THDMcache::B00_Mz_Mz2_mh_mA(const double Mz, const double mh, const double mA) const {
    int NumPar = 3;
    double params[] = {Mz, mh, mA};

    int i = CacheCheck(B00_Mz_Mz2_mh_mA_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_Mz2_mh_mA_cache[NumPar][i] );
    } else {
        complex newResult = PV.B00(Mz*Mz, Mz*Mz, mh*mh, mA*mA);
        CacheShift(B00_Mz_Mz2_mh_mA_cache, NumPar, params, newResult);
        return newResult;
    } 
}

complex THDMcache::B00_Mz_0_mh_mA(const double Mz, const double mh, const double mA) const {
    int NumPar = 3;
    double params[] = {Mz, mh, mA};

    int i = CacheCheck(B00_Mz_0_mh_mA_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_0_mh_mA_cache[NumPar][i] );
    } else {
        complex newResult = PV.B00(Mz*Mz, 0., mh*mh, mA*mA);
        CacheShift(B00_Mz_0_mh_mA_cache, NumPar, params, newResult);
        return newResult;
    } 
}

complex THDMcache::B00_Mz_Mz2_Mz_mH(const double Mz, const double mH) const {
    int NumPar = 2;
    double params[] = {Mz, mH};

    int i = CacheCheck(B00_Mz_Mz2_Mz_mH_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_Mz2_Mz_mH_cache[NumPar][i] );
    } else {
        complex newResult = PV.B00(Mz*Mz, Mz*Mz, Mz*Mz, mH*mH);
        CacheShift(B00_Mz_Mz2_Mz_mH_cache, NumPar, params, newResult);
        return newResult;
    } 
}

complex THDMcache::B00_Mz_0_Mz_mH(const double Mz, const double mH) const {
    int NumPar = 2;
    double params[] = {Mz, mH};

    int i = CacheCheck(B00_Mz_0_Mz_mH_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_0_Mz_mH_cache[NumPar][i] );
    } else {
        complex newResult = PV.B00(Mz*Mz, 0., Mz*Mz, mH*mH);
        CacheShift(B00_Mz_0_Mz_mH_cache, NumPar, params, newResult);
        return newResult;
    } 
}

complex THDMcache::B00_Mz_Mz2_Mz_mh(const double Mz, const double mh) const {
    int NumPar = 2;
    double params[] = {Mz, mh};

    int i = CacheCheck(B00_Mz_Mz2_Mz_mh_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_Mz2_Mz_mh_cache[NumPar][i] );
    } else {
        complex newResult = PV.B00(Mz*Mz, Mz*Mz, Mz*Mz, mh*mh);
        CacheShift(B00_Mz_Mz2_Mz_mh_cache, NumPar, params, newResult);
        return newResult;
    } 
}

complex THDMcache::B00_Mz_0_Mz_mh(const double Mz, const double mh) const {
    int NumPar = 2;
    double params[] = {Mz, mh};

    int i = CacheCheck(B00_Mz_0_Mz_mh_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_0_Mz_mh_cache[NumPar][i] );
    } else {
        complex newResult = PV.B00(Mz*Mz, 0., Mz*Mz, mh*mh);
        CacheShift(B00_Mz_0_Mz_mh_cache, NumPar, params, newResult);
        return newResult;
    } 
}

complex THDMcache::B00_Mz_Mw2_mA_mHp(const double Mz, const double Mw, const double mA, const double mHp) const {
    int NumPar = 4;
    double params[] = {Mz,Mw, mA, mHp};

    int i = CacheCheck(B00_Mz_Mw2_mA_mHp_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_Mw2_mA_mHp_cache[NumPar][i] );
    } else {
        complex newResult = PV.B00(Mz*Mz, Mw*Mw, mA*mA, mHp*mHp);
        CacheShift(B00_Mz_Mw2_mA_mHp_cache, NumPar, params, newResult);
        return newResult;
    } 
}

complex THDMcache::B00_Mz_0_mA_mHp(const double Mz, const double mA, const double mHp) const {
    int NumPar = 3;
    double params[] = {Mz, mA, mHp};

    int i = CacheCheck(B00_Mz_0_mA_mHp_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_0_mA_mHp_cache[NumPar][i] );
    } else {
        complex newResult = PV.B00(Mz*Mz, 0., mA*mA, mHp*mHp);
        CacheShift(B00_Mz_0_mA_mHp_cache, NumPar, params, newResult);
        return newResult;
    } 
}

complex THDMcache::B00_Mz_Mw2_mHp_mHp(const double Mz, const double Mw, const double mHp) const {
    int NumPar = 3;
    double params[] = {Mz,Mw, mHp};

    int i = CacheCheck(B00_Mz_Mw2_mHp_mHp_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_Mw2_mHp_mHp_cache[NumPar][i] );
    } else {
        complex newResult = PV.B00(Mz*Mz, Mw*Mw, mHp*mHp, mHp*mHp);
        CacheShift(B00_Mz_Mw2_mHp_mHp_cache, NumPar, params, newResult);
        return newResult;
    } 
}

complex THDMcache::B00_Mz_Mw2_Mw_mH(const double Mz, const double Mw, const double mH) const {
    int NumPar = 3;
    double params[] = {Mz,Mw, mH};

    int i = CacheCheck(B00_Mz_Mw2_Mw_mH_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_Mw2_Mw_mH_cache[NumPar][i] );
    } else {
        complex newResult = PV.B00(Mz*Mz, Mw*Mw, Mw*Mw, mH*mH);
        CacheShift(B00_Mz_Mw2_Mw_mH_cache, NumPar, params, newResult);
        return newResult;
    } 
}

complex THDMcache::B00_Mz_0_Mw_mH(const double Mz, const double Mw, const double mH) const {
    int NumPar = 3;
    double params[] = {Mz,Mw, mH};

    int i = CacheCheck(B00_Mz_0_Mw_mH_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_0_Mw_mH_cache[NumPar][i] );
    } else {
        complex newResult = PV.B00(Mz*Mz, Mw*Mw, Mw*Mw, mH*mH);
        CacheShift(B00_Mz_0_Mw_mH_cache, NumPar, params, newResult);
        return newResult;
    } 
}


complex THDMcache::B00_Mz_Mw2_Mw_mh(const double Mz, const double Mw, const double mh) const {
    int NumPar = 3;
    double params[] = {Mz,Mw, mh};

    int i = CacheCheck(B00_Mz_Mw2_Mw_mh_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_Mw2_Mw_mh_cache[NumPar][i] );
    } else {
        complex newResult = PV.B00(Mz*Mz, Mw*Mw, Mw*Mw, mh*mh);
        CacheShift(B00_Mz_Mw2_Mw_mh_cache, NumPar, params, newResult);
        return newResult;
    } 
}

complex THDMcache::B00_Mz_0_Mw_mh(const double Mz, const double Mw, const double mh) const {
    int NumPar = 3;
    double params[] = {Mz,Mw, mh};

    int i = CacheCheck(B00_Mz_0_Mw_mh_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_0_Mw_mh_cache[NumPar][i] );
    } else {
        complex newResult = PV.B00(Mz*Mz, 0., Mw*Mw, mh*mh);
        CacheShift(B00_Mz_0_Mw_mh_cache, NumPar, params, newResult);
        return newResult;
    } 
}

complex THDMcache::B00_Mz_Mw2_mH_mHp(const double Mz, const double Mw, const double mH, const double mHp) const {
    int NumPar = 4;
    double params[] = {Mz,Mw, mH, mHp};

    int i = CacheCheck(B00_Mz_Mw2_mH_mHp_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_Mw2_mH_mHp_cache[NumPar][i] );
    } else {
        complex newResult = PV.B00(Mz*Mz, Mw*Mw, mH*mH, mHp*mHp);
        CacheShift(B00_Mz_Mw2_mH_mHp_cache, NumPar, params, newResult);
        return newResult;
    } 
}

complex THDMcache::B00_Mz_0_mH_mHp(const double Mz, const double mH, const double mHp) const {
    int NumPar = 3;
    double params[] = {Mz, mH, mHp};

    int i = CacheCheck(B00_Mz_0_mH_mHp_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_0_mH_mHp_cache[NumPar][i] );
    } else {
        complex newResult = PV.B00(Mz*Mz, 0., mH*mH, mHp*mHp);
        CacheShift(B00_Mz_0_mH_mHp_cache, NumPar, params, newResult);
        return newResult;
    } 
}

complex THDMcache::B00_Mz_Mw2_mh_mHp(const double Mz, const double Mw, const double mh, const double mHp) const {
    int NumPar = 4;
    double params[] = {Mz,Mw, mh, mHp};

    int i = CacheCheck(B00_Mz_Mw2_mh_mHp_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_Mw2_mh_mHp_cache[NumPar][i] );
    } else {
        complex newResult = PV.B00(Mz*Mz, Mw*Mw, mh*mh, mHp*mHp);
        CacheShift(B00_Mz_Mw2_mh_mHp_cache, NumPar, params, newResult);
        return newResult;
    } 
}

complex THDMcache::B00_Mz_0_mh_mHp(const double Mz, const double mh, const double mHp) const {
    int NumPar = 3;
    double params[] = {Mz, mh, mHp};

    int i = CacheCheck(B00_Mz_0_mh_mHp_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_0_mh_mHp_cache[NumPar][i] );
    } else {
        complex newResult = PV.B00(Mz*Mz, 0., mh*mh, mHp*mHp);
        CacheShift(B00_Mz_0_mh_mHp_cache, NumPar, params, newResult);
        return newResult;
    } 
}
































