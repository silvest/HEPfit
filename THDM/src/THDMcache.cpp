/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "THDMcache.h"
#include <fstream>
#include "gslpp.h"
#include <sstream>
#include <string>

THDMcache::THDMcache()
        
        : PV(true),
        
        br_tt(19861, 2, 0.),
        br_bb(19861, 2, 0.),
        br_tautau(19861, 2, 0.),
        br_cc(19861, 2, 0.),
        br_mumu(19861, 2, 0.),
        br_ZZ(19861, 2, 0.),
        br_WW(19861, 2, 0.),
        pc_ggF(1971, 2, 0.),
        pc_VBF(1971, 2, 0.),
        pc_WH(1971, 2, 0.),
        pc_ZH(1971, 2, 0.),
        pc_ttH(1971, 2, 0.),
        GammaHtotSM(19861, 2, 0.),
        cs_ggH(186, 2, 0.),
        cs_ggH_tt(186, 2, 0.),
        cs_ggH_bb(186, 2, 0.),
        cs_ggA(186, 2, 0.),
        cs_ggA_tt(186, 2, 0.),
        cs_ggA_bb(186, 2, 0.),
        cs_bbFtoHP(185, 2, 0.),
        ATLAS_ggF_phi_gaga(991, 2, 0.),
        ATLAS_ggF_phi_tautau(496, 2, 0.),
        ATLAS_bbF_phi_tautau(496, 2, 0.),
        ATLAS_ggF_A_hZ_tautauZ(986, 2, 0.),
        ATLAS_ggF_A_hZ_bbZ(986, 2, 0.),
        ATLAS_pp_phi_tt(200, 2, 0.),
        ATLAS_ggF_H_WW(100,2,0.),
        ATLAS_VBF_H_WW(100,2,0.),
        ATLAS_ggF_H_hh(1000,2,0.),
        CMS_pp_H_ZZ(9851, 2, 0.),
        CMS_ggF_A_hZ_bbll(986, 2, 0.),
        CMS_pp_H_hh_gagabb(496, 2, 0.),
        CMS_pp_H_hh_bbbb(496, 2, 0.),
        CMS_bbF_phi_bb(1000, 2, 0.),
        CMS_ggF_phi_tautau(1000,2,0.),
        CMS_bbF_phi_tautau(1000,2,0.),
        CMS_ggF_phi_gaga(2000,2,0.),
        CMS_ggF_H_hh_bbtautau(1000,2,0.),
        CMS_ggF_A_hZ_tautaull(1000,2,0.),
        arraybsgamma(1111, 3, 0.)
{
  read();
}


/////////////////////////////////////////////////////////////////////////////////////////////////

int THDMcache::CacheCheck(const gslpp::complex cache[][CacheSize], 
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

void THDMcache::CacheShift(gslpp::complex cache[][CacheSize], const int NumPar, 
                           const double params[], const gslpp::complex newResult) const {
    // shift old parameters and result
    for(int i=CacheSize-1; i>0; i--)
        for(int j=0; j<NumPar+1; j++)
            cache[j][i] = cache[j][i-1];

    // store new parameters and result
    for(int j=0; j<NumPar; j++) {
        cache[j][0] = gslpp::complex(params[j], 0.0, false);
        cache[NumPar][0] = newResult;
    }
}

////////////////////////////////////////////////////////////////////////////////
/*One-loop functions*/
////////////////////////////////////////////////////////////////////////////////

gslpp::complex THDMcache::B0_MZ2_0_MW2_mHh2(const double MZ2, const double MW2, const double mHh2) const {
    int NumPar = 3;
    double params[] = {MZ2, MW2, mHh2};

    int i = CacheCheck(B0_MZ2_0_MW2_mHh2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_0_MW2_mHh2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, 0., MW2, mHh2);
        CacheShift(B0_MZ2_0_MW2_mHh2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B0_MZ2_0_MW2_mHl2(const double MZ2, const double MW2, const double mHl2) const {
    int NumPar = 3;
    double params[] = {MZ2, MW2, mHl2};

    int i = CacheCheck(B0_MZ2_0_MW2_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_0_MW2_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, 0., MW2, mHl2);
        CacheShift(B0_MZ2_0_MW2_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B0_MZ2_0_MZ2_mHh2(const double MZ2, const double mHh2) const {
    int NumPar = 2;
    double params[] = {MZ2, mHh2};

    int i = CacheCheck(B0_MZ2_0_MZ2_mHh2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_0_MZ2_mHh2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, 0., MZ2, mHh2);
        CacheShift(B0_MZ2_0_MZ2_mHh2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B0_MZ2_0_MZ2_mHl2(const double MZ2, const double mHl2) const {
    int NumPar = 2;
    double params[] = {MZ2, mHl2};

    int i = CacheCheck(B0_MZ2_0_MZ2_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_0_MZ2_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, 0., MZ2, mHl2);
        CacheShift(B0_MZ2_0_MZ2_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B0_MZ2_MW2_MW2_mHh2(const double MZ2, const double MW2, const double mHh2) const {
    int NumPar = 3;
    double params[] = {MZ2, MW2, mHh2};

    int i = CacheCheck(B0_MZ2_MW2_MW2_mHh2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_MW2_MW2_mHh2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, MW2, MW2, mHh2);
        CacheShift(B0_MZ2_MW2_MW2_mHh2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B0_MZ2_MW2_MW2_mHl2(const double MZ2, const double MW2, const double mHl2) const {
    int NumPar = 3;
    double params[] = {MZ2, MW2, mHl2};

    int i = CacheCheck(B0_MZ2_MW2_MW2_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_MW2_MW2_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, MW2, MW2, mHl2);
        CacheShift(B0_MZ2_MW2_MW2_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B0_MZ2_MZ2_MZ2_mHh2(const double MZ2, const double mHh2) const {
    int NumPar = 2;
    double params[] = {MZ2, mHh2};

    int i = CacheCheck(B0_MZ2_MZ2_MZ2_mHh2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_MZ2_MZ2_mHh2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, MZ2, MZ2, mHh2);
        CacheShift(B0_MZ2_MZ2_MZ2_mHh2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B0_MZ2_MZ2_MZ2_mHl2(const double MZ2, const double mHl2) const {
    int NumPar = 2;
    double params[] = {MZ2, mHl2};

    int i = CacheCheck(B0_MZ2_MZ2_MZ2_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_MZ2_MZ2_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, MZ2, MZ2, mHl2);
        CacheShift(B0_MZ2_MZ2_MZ2_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

///////////////////////////////////////////////////////////////////////////////////////////

gslpp::complex THDMcache::B00_MZ2_0_mA2_mHp2(const double MZ2, const double mA2, const double mHp2) const {
    int NumPar = 3;
    double params[] = {MZ2, mA2, mHp2};

    int i = CacheCheck(B00_MZ2_0_mA2_mHp2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_0_mA2_mHp2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, 0., mA2, mHp2);
        CacheShift(B00_MZ2_0_mA2_mHp2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_MZ2_0_mHh2_mA2(const double MZ2, const double mHh2, const double mA2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHh2, mA2};

    int i = CacheCheck(B00_MZ2_0_mHh2_mA2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_0_mHh2_mA2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, 0., mHh2, mA2);
        CacheShift(B00_MZ2_0_mHh2_mA2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_MZ2_0_mHh2_mHp2(const double MZ2, const double mHh2, const double mHp2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHh2, mHp2};

    int i = CacheCheck(B00_MZ2_0_mHh2_mHp2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_0_mHh2_mHp2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, 0., mHh2, mHp2);
        CacheShift(B00_MZ2_0_mHh2_mHp2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_MZ2_0_mHl2_mA2(const double MZ2, const double mHl2, const double mA2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHl2, mA2};

    int i = CacheCheck(B00_MZ2_0_mHl2_mA2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_0_mHl2_mA2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, 0., mHl2, mA2);
        CacheShift(B00_MZ2_0_mHl2_mA2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_MZ2_0_mHl2_mHp2(const double MZ2, const double mHl2, const double mHp2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHl2, mHp2};

    int i = CacheCheck(B00_MZ2_0_mHl2_mHp2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_0_mHl2_mHp2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, 0., mHl2, mHp2);
        CacheShift(B00_MZ2_0_mHl2_mHp2_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex THDMcache::B00_MZ2_0_mHp2_mHp2(const double MZ2, const double mHp2) const {
    int NumPar = 2;
    double params[] = {MZ2, mHp2};

    int i = CacheCheck(B00_MZ2_0_mHp2_mHp2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_0_mHp2_mHp2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, 0., mHp2, mHp2);
        CacheShift(B00_MZ2_0_mHp2_mHp2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_MZ2_0_MW2_mHh2(const double MZ2, const double MW2, const double mHh2) const {
    int NumPar = 3;
    double params[] = {MZ2, MW2, mHh2};

    int i = CacheCheck(B00_MZ2_0_MW2_mHh2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_0_MW2_mHh2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, MW2, MW2, mHh2);
        CacheShift(B00_MZ2_0_MW2_mHh2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_MZ2_0_MW2_mHl2(const double MZ2, const double MW2, const double mHl2) const {
    int NumPar = 3;
    double params[] = {MZ2, MW2, mHl2};

    int i = CacheCheck(B00_MZ2_0_MW2_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_0_MW2_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, 0., MW2, mHl2);
        CacheShift(B00_MZ2_0_MW2_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_MZ2_0_MZ2_mHh2(const double MZ2, const double mHh2) const {
    int NumPar = 2;
    double params[] = {MZ2, mHh2};

    int i = CacheCheck(B00_MZ2_0_MZ2_mHh2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_0_MZ2_mHh2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, 0., MZ2, mHh2);
        CacheShift(B00_MZ2_0_MZ2_mHh2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_MZ2_0_MZ2_mHl2(const double MZ2, const double mHl2) const {
    int NumPar = 2;
    double params[] = {MZ2, mHl2};

    int i = CacheCheck(B00_MZ2_0_MZ2_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_0_MZ2_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, 0., MZ2, mHl2);
        CacheShift(B00_MZ2_0_MZ2_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_MZ2_MW2_mA2_mHp2(const double MZ2, const double MW2, const double mA2, const double mHp2) const {
    int NumPar = 4;
    double params[] = {MZ2, MW2, mA2, mHp2};

    int i = CacheCheck(B00_MZ2_MW2_mA2_mHp2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_MW2_mA2_mHp2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, MW2, mA2, mHp2);
        CacheShift(B00_MZ2_MW2_mA2_mHp2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_MZ2_MW2_mHh2_mHp2(const double MZ2, const double MW2, const double mHh2, const double mHp2) const {
    int NumPar = 4;
    double params[] = {MZ2, MW2, mHh2, mHp2};

    int i = CacheCheck(B00_MZ2_MW2_mHh2_mHp2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_MW2_mHh2_mHp2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, MW2, mHh2, mHp2);
        CacheShift(B00_MZ2_MW2_mHh2_mHp2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_MZ2_MW2_mHl2_mHp2(const double MZ2, const double MW2, const double mHl2, const double mHp2) const {
    int NumPar = 4;
    double params[] = {MZ2, MW2, mHl2, mHp2};

    int i = CacheCheck(B00_MZ2_MW2_mHl2_mHp2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_MW2_mHl2_mHp2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, MW2, mHl2, mHp2);
        CacheShift(B00_MZ2_MW2_mHl2_mHp2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_MZ2_MW2_mHp2_mHp2(const double MZ2, const double MW2, const double mHp2) const {
    int NumPar = 3;
    double params[] = {MZ2, MW2, mHp2};

    int i = CacheCheck(B00_MZ2_MW2_mHp2_mHp2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_MW2_mHp2_mHp2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, MW2, mHp2, mHp2);
        CacheShift(B00_MZ2_MW2_mHp2_mHp2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_MZ2_MW2_MW2_mHh2(const double MZ2, const double MW2, const double mHh2) const {
    int NumPar = 3;
    double params[] = {MZ2, MW2, mHh2};

    int i = CacheCheck(B00_MZ2_MW2_MW2_mHh2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_MW2_MW2_mHh2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, MW2, MW2, mHh2);
        CacheShift(B00_MZ2_MW2_MW2_mHh2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_MZ2_MW2_MW2_mHl2(const double MZ2, const double MW2, const double mHl2) const {
    int NumPar = 3;
    double params[] = {MZ2, MW2, mHl2};

    int i = CacheCheck(B00_MZ2_MW2_MW2_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_MW2_MW2_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, MW2, MW2, mHl2);
        CacheShift(B00_MZ2_MW2_MW2_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_MZ2_MZ2_mHh2_mA2(const double MZ2, const double mHh2, const double mA2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHh2, mA2};

    int i = CacheCheck(B00_MZ2_MZ2_mHh2_mA2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_MZ2_mHh2_mA2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, MZ2, mHh2, mA2);
        CacheShift(B00_MZ2_MZ2_mHh2_mA2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_MZ2_MZ2_mHl2_mA2(const double MZ2, const double mHl2, const double mA2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHl2, mA2};

    int i = CacheCheck(B00_MZ2_MZ2_mHl2_mA2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_MZ2_mHl2_mA2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, MZ2, mHl2, mA2);
        CacheShift(B00_MZ2_MZ2_mHl2_mA2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_MZ2_MZ2_mHp2_mHp2(const double MZ2, const double mHp2) const {
    int NumPar = 2;
    double params[] = {MZ2, mHp2};

    int i = CacheCheck(B00_MZ2_MZ2_mHp2_mHp2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_MZ2_mHp2_mHp2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, MZ2, mHp2, mHp2);
        CacheShift(B00_MZ2_MZ2_mHp2_mHp2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_MZ2_MZ2_MZ2_mHh2(const double MZ2, const double mHh2) const {
    int NumPar = 2;
    double params[] = {MZ2, mHh2};

    int i = CacheCheck(B00_MZ2_MZ2_MZ2_mHh2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_MZ2_MZ2_mHh2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, MZ2, MZ2, mHh2);
        CacheShift(B00_MZ2_MZ2_MZ2_mHh2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_MZ2_MZ2_MZ2_mHl2(const double MZ2, const double mHl2) const {
    int NumPar = 2;
    double params[] = {MZ2, mHl2};

    int i = CacheCheck(B00_MZ2_MZ2_MZ2_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_MZ2_MZ2_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, MZ2, MZ2, mHl2);
        CacheShift(B00_MZ2_MZ2_MZ2_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

void THDMcache::read(){

    std::string tablepath="/Users/Roma1/Desktop/HEPfit/THDM/tabs/";
    std::stringstream br1,br2,br3,br4,br5,br6,br7;
    std::stringstream pc1,pc2,pc3,pc4,pc5;
    std::stringstream dw1;
    std::stringstream cs1,cs2,cs3,cs4,cs5,cs6,cs7;
    std::stringstream ex1,ex2,ex3,ex4,ex5,ex6,ex7,ex8,ex9,ex10,ex11,ex12,ex13,ex14,ex15,ex16,ex17,ex18,ex19;
    std::stringstream bsg1;

    std::cout<<"reading tables"<<std::endl;

    br1 << tablepath << "br1.dat";
    br_tt = readTable(br1.str(),19861,2);
    br2 << tablepath << "br2.dat";
    br_bb = readTable(br2.str(),19861,2);
    br3 << tablepath << "br3.dat";
    br_tautau = readTable(br3.str(),19861,2); 
    br4 << tablepath << "br4.dat";
    br_cc = readTable(br4.str(),19861,2);
    br5 << tablepath << "br5.dat";
    br_mumu = readTable(br5.str(),19861,2);
    br6 << tablepath << "br6.dat";
    br_ZZ = readTable(br6.str(),19861,2);
    br7 << tablepath << "br7.dat";
    br_WW = readTable(br7.str(),19861,2);
    pc1 << tablepath << "pc1.dat";
    pc_ggF = readTable(pc1.str(),1971,2);
    pc2 << tablepath << "pc2.dat";
    pc_VBF = readTable(pc2.str(),1971,2);
    pc3 << tablepath << "pc3.dat";
    pc_WH = readTable(pc3.str(),1971,2);
    pc4 << tablepath << "pc4.dat";
    pc_ZH = readTable(pc4.str(),1971,2);
    pc5 << tablepath << "pc5.dat";
    pc_ttH = readTable(pc5.str(),1971,2);
    dw1 << tablepath << "dw1.dat";
    GammaHtotSM = readTable(dw1.str(),19861,2);
    cs1 << tablepath << "cs1.dat";
    cs_ggH = readTable(cs1.str(),186,2);
    cs2 << tablepath << "cs2.dat";
    cs_ggH_tt = readTable(cs2.str(),186,2);
    cs3 << tablepath << "cs3.dat";
    cs_ggH_bb = readTable(cs3.str(),186,2);
    cs4 << tablepath << "cs4.dat";
    cs_ggA = readTable(cs4.str(),186,2);
    cs5 << tablepath << "cs5.dat";
    cs_ggA_tt = readTable(cs5.str(),186,2);
    cs6 << tablepath << "cs6.dat";
    cs_ggA_bb = readTable(cs6.str(),186,2);
    cs7 << tablepath << "cs7.dat";
    cs_bbFtoHP = readTable(cs7.str(),185,2);
    ex5 << tablepath << "ATLAS_ggF_phi_gaga.dat";
    ATLAS_ggF_phi_gaga = readTable(ex5.str(),991,2);
    ex6 << tablepath << "ATLAS_ggF_phi_tautau.dat";
    ATLAS_ggF_phi_tautau = readTable(ex6.str(),496,2);
    ex7 << tablepath << "ATLAS_bbF_phi_tautau.dat";
    ATLAS_bbF_phi_tautau = readTable(ex7.str(),496,2);
    ex8 << tablepath << "ATLAS_ggF_A_hZ_tautauZ.dat";
    ATLAS_ggF_A_hZ_tautauZ = readTable(ex8.str(),986,2);
    ex9 << tablepath << "ATLAS_ggF_A_hZ_bbZ.dat";
    ATLAS_ggF_A_hZ_bbZ = readTable(ex9.str(),986,2);
    ex11 << tablepath << "ATLAS_pp_phi_tt.dat";
    ATLAS_pp_phi_tt = readTable(ex11.str(),200,2);
    ex15 << tablepath << "ATLAS_ggF_H_WW.dat";
    ATLAS_ggF_H_WW = readTable(ex15.str(),100,2);
    ex16 << tablepath << "ATLAS_VBF_H_WW.dat";
    ATLAS_VBF_H_WW = readTable(ex16.str(),100,2);
    ex17 << tablepath << "ATLAS_ggF_H_hh.dat";
    ATLAS_ggF_H_hh = readTable(ex17.str(),1000,2);
    ex1 << tablepath << "CMS_pp_H_ZZ.dat";
    CMS_pp_H_ZZ = readTable(ex1.str(),9851,2);
    ex2 << tablepath << "CMS_ggF_A_hZ_bbll.dat";
    CMS_ggF_A_hZ_bbll = readTable(ex2.str(),986,2);
    ex3 << tablepath << "CMS_pp_H_hh_gagabb.dat";
    CMS_pp_H_hh_gagabb = readTable(ex3.str(),496,2);
    ex4 << tablepath << "CMS_pp_H_hh_bbbb.dat";
    CMS_pp_H_hh_bbbb = readTable(ex4.str(),496,2);
    ex10 << tablepath << "CMS_bbF_phi_bb.dat";
    CMS_bbF_phi_bb = readTable(ex10.str(),1000,2);
    ex12 << tablepath << "CMS_ggF_phi_tautau.dat";
    CMS_ggF_phi_tautau = readTable(ex12.str(),1000,2);
    ex13 << tablepath << "CMS_bbF_phi_tautau.dat";
    CMS_bbF_phi_tautau = readTable(ex13.str(),1000,2);
    ex14 << tablepath << "CMS_ggF_phi_gaga.dat";
    CMS_ggF_phi_gaga = readTable(ex14.str(),2000,2);
    ex18 << tablepath << "CMS_ggF_H_hh_bbtautau.dat";
    CMS_ggF_H_hh_bbtautau = readTable(ex18.str(),1000,2);
    ex19 << tablepath << "CMS_ggF_A_hZ_tautaull.dat";
    CMS_ggF_A_hZ_tautaull = readTable(ex19.str(),1000,2);
    bsg1 << tablepath << "bsgammatable.dat";
    arraybsgamma = readTable(bsg1.str(),1111,3);
}    
    


double THDMcache::ip_Br_HPtott(double mass){
    return pow(10.0,interpolate(br_tt,mass));
}



double THDMcache::ip_Br_HPtobb(double mass){
    return pow(10.0,interpolate(br_bb,mass));
}



double THDMcache::ip_Br_HPtotautau(double mass){
    return pow(10.0,interpolate(br_tautau,mass));
}



double THDMcache::ip_Br_HPtocc(double mass){
    return pow(10.0,interpolate(br_cc,mass));
}



double THDMcache::ip_Br_HPtomumu(double mass){
    return pow(10.0,interpolate(br_mumu,mass));
}



double THDMcache::ip_Br_HPtoZZ(double mass){
    return pow(10.0,interpolate(br_ZZ,mass));
}



double THDMcache::ip_Br_HPtoWW(double mass){
    return pow(10.0,interpolate(br_WW,mass));
}



double THDMcache::ip_pc_ggFtoHP(double mass){
    return interpolate(pc_ggF,mass);
}



double THDMcache::ip_pc_VBFtoHP(double mass){
    return interpolate(pc_VBF,mass);
}



double THDMcache::ip_pc_WHP_HP(double mass){
    return interpolate(pc_WH,mass);
}


double THDMcache::ip_pc_ZHP_HP(double mass){
    return interpolate(pc_ZH,mass);
}



double THDMcache::ip_pc_ttFtoHP(double mass){
    return interpolate(pc_ttH,mass);
}



double THDMcache::ip_GammaHPtotSM(double mass){
    return pow(10.0,interpolate(GammaHtotSM,mass));
}



double THDMcache::ip_cs_ggFtoHP(double mass){
    return pow(10.0,interpolate (cs_ggH,log10(mass)));
}


double THDMcache::ip_cs_ggHP_tt(double mass){
    return pow(10.0,interpolate (cs_ggH_tt,log10(mass)));
}



double THDMcache::ip_cs_ggHP_bb(double mass){
    return pow(10.0,interpolate (cs_ggH_bb,log10(mass)));
}



double THDMcache::ip_cs_ggA(double mass){
    return pow(10.0,interpolate (cs_ggA,log10(mass)));
}



double THDMcache::ip_cs_ggA_tt(double mass){
    return pow(10.0,interpolate (cs_ggA_tt,log10(mass)));
}



double THDMcache::ip_cs_ggA_bb(double mass){
    return pow(10.0,interpolate (cs_ggA_bb,log10(mass)));
}



double THDMcache::ip_cs_bbFtoHP(double mass){
    return pow(10.0,interpolate (cs_bbFtoHP,log10(mass)));
}



double THDMcache::ip_ex_ggF_phi_gaga_ATLAS(double mass){
    return interpolate(ATLAS_ggF_phi_gaga,mass);
}



double THDMcache::ip_ex_ggF_phi_tautau_ATLAS(double mass){
    return interpolate(ATLAS_ggF_phi_tautau,mass);
}



double THDMcache::ip_ex_bbF_phi_tautau_ATLAS(double mass){
    return interpolate(ATLAS_bbF_phi_tautau,mass);
}



double THDMcache::ip_ex_ggF_A_hZ_tautauZ_ATLAS(double mass){
    return interpolate(ATLAS_ggF_A_hZ_tautauZ,mass);
}



double THDMcache::ip_ex_ggF_A_hZ_bbZ_ATLAS(double mass){
    return interpolate(ATLAS_ggF_A_hZ_bbZ,mass);
}



double THDMcache::ip_ex_pp_phi_tt_ATLAS(double mass){
    return interpolate (ATLAS_pp_phi_tt,mass);
}



double THDMcache::ip_ex_ggF_H_WW_ATLAS(double mass){
    return interpolate (ATLAS_ggF_H_WW,mass);
}



double THDMcache::ip_ex_VBF_H_WW_ATLAS(double mass){
    return interpolate (ATLAS_VBF_H_WW,mass);
}



double THDMcache::ip_ex_ggF_H_hh_ATLAS(double mass){
    return interpolate (ATLAS_ggF_H_hh,mass);
}



double THDMcache::ip_ex_pp_H_ZZ_CMS(double mass){
    return interpolate(CMS_pp_H_ZZ,mass);
}



double THDMcache::ip_ex_ggF_A_hZ_bbll_CMS(double mass){
    return interpolate(CMS_ggF_A_hZ_bbll,mass);
}



double THDMcache::ip_ex_pp_phi_hh_gagabb_CMS(double mass){
    return interpolate(CMS_pp_H_hh_gagabb,mass);
}



double THDMcache::ip_ex_pp_phi_hh_bbbb_CMS(double mass){
    return interpolate(CMS_pp_H_hh_bbbb,mass);
}



double THDMcache::ip_ex_bbF_phi_bb_CMS(double mass){
    return interpolate (CMS_bbF_phi_bb,mass);
}



double THDMcache::ip_ex_ggF_phi_tautau_CMS(double mass){
    return interpolate (CMS_ggF_phi_tautau,mass);
}



double THDMcache::ip_ex_bbF_phi_tautau_CMS(double mass){
    return interpolate (CMS_bbF_phi_tautau,mass);
}



double THDMcache::ip_ex_ggF_phi_gaga_CMS(double mass){
    return interpolate (CMS_ggF_phi_gaga,mass);
}



double THDMcache::ip_ex_ggF_H_hh_bbtautau_CMS(double mass){
    return interpolate (CMS_ggF_H_hh_bbtautau,mass);
}



double THDMcache::ip_ex_ggF_A_hZ_tautaull_CMS(double mass){
    return interpolate (CMS_ggF_A_hZ_tautaull,mass);
}



double THDMcache::ip_ex_bsgamma(double logtb, double logmHp){
    return interpolate2D(arraybsgamma, logtb, logmHp);
}



gslpp::matrix<double> THDMcache::readTable(std::string filename, int rowN, int colN){

    std::ifstream INfile;
    std::string lineTab;
    INfile.open( filename.c_str() );
    if(INfile.fail()){
        std::cout<<"error: in HEAVY"<<" doesn't exist!"<<std::endl;
    }

    gslpp::matrix<double> arrayTab(rowN, colN, 0.);
    int a =0;
    int b=0;
    double v;

    while(INfile.good()){
        while(getline(INfile, lineTab)){
            if( lineTab[0]=='#' )continue;
            else{
            std::istringstream streamTab(lineTab);
            b=0;
            while(streamTab >>v){
                arrayTab.assign(a,b,v);
                b++;
            }
            a++;
            }
        }
    }
    
    INfile.close();
    
    return arrayTab;
}

//1D interpolation

double THDMcache::interpolate(gslpp::matrix<double> arrayTab, double x){

    int rowN=arrayTab.size_i();
    
    double xmin = arrayTab(0,0);
    double xmax = arrayTab(rowN-1,0);
    double interval = arrayTab(1,0)-arrayTab(0,0);
    int Nintervals = (x-xmin)/interval;
    double y = 0.0;
       
    if(x<xmin){
        std::cout<<"error: your table parameter value is smaller than the minimum allowed value"<<std::endl;
        return 0.;
    }
    else if(x>xmax){
        std::cout<<"error: your table parameter value is greater than the maximum allowed value"<<std::endl;
        return 0.;
    }
    else{
        y =(arrayTab(Nintervals+1,1)-arrayTab(Nintervals,1))/(arrayTab(Nintervals+1,0)
                   -arrayTab(Nintervals,0))*(x-arrayTab(Nintervals,0))+arrayTab(Nintervals,1);
        return y;
    }
}

//2D interpolation

double THDMcache::interpolate2D(gslpp::matrix<double> arrayTab, double x, double y){

    int rowN=arrayTab.size_i();
    
    double xmin = arrayTab(0,0);
    double xmax = arrayTab(rowN-1,0);
    double ymin = arrayTab(0,1);
    double ymax = arrayTab(rowN-1,1);
    double intervalx = arrayTab(1,0)-arrayTab(0,0);
    int i=1;
    do i++;
    while(arrayTab(i,1)-arrayTab(i-1,1)==0&&i<10000);
    double intervaly = arrayTab(i,1)-arrayTab(i-1,1);
    int Nintervalsx = (x-xmin)/intervalx;
    int Nintervalsy = (y-ymin)/intervaly;
    if(x<xmin||x>xmax||y<ymin||y>ymax){
        std::cout<<"error: the parameter point lies outside the table"<<std::endl;
        return 0.;
    }
    else{
    double x1=arrayTab(i*Nintervalsy+Nintervalsx,0);
    double x2=arrayTab(i*Nintervalsy+Nintervalsx+1,0);
    double y1=arrayTab(i*Nintervalsy+Nintervalsx,1);
    double y2=arrayTab(i*(Nintervalsy+1)+Nintervalsx,1);
    return (arrayTab(i*Nintervalsy+Nintervalsx,2) * (x2-x) * (y2-y)
            +arrayTab(i*Nintervalsy+Nintervalsx+1,2) * (x-x1) * (y2-y)
            +arrayTab(i*(Nintervalsy+1)+Nintervalsx,2) * (x2-x) * (y-y1)
            +arrayTab(i*(Nintervalsy+1)+Nintervalsx+1,2) * (x-x1) * (y-y1))
           /((x2-x1)*(y2-y1));
    }
}
