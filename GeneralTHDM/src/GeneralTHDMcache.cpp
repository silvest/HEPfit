/* 
 * Copyright (C) 2016 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMcache.h"
#include <fstream>
#include "gslpp.h"
#include <sstream>
#include <string>

//#include "log_cs_ggH_13.h"



GeneralTHDMcache::GeneralTHDMcache(const StandardModel& SM_i)
:       br_tt(19961, 2, 0.),
        br_bb(19961, 2, 0.),
        br_tautau(19961, 2, 0.),
        br_cc(19961, 2, 0.),
        br_mumu(19961, 2, 0.),
        br_ZZ(19961, 2, 0.),
        br_WW(19961, 2, 0.),
        GammaHtot_SM(19961, 2, 0.),
        log_cs_ggH_8(199, 2, 0.),
        log_cs_VBF_8(199, 2, 0.),
        log_cs_WH_8(199, 2, 0.),
        log_cs_ZH_8(199, 2, 0.),
        log_cs_ggH_13(199, 2, 0.),
        log_cs_VBF_13(199, 2, 0.),
        log_cs_WH_13(199, 2, 0.),
        log_cs_ZH_13(199, 2, 0.),
        log_cs_ttH_8(199, 2, 0.),
        log_cs_ttH_13(199, 2, 0.),
        log_cs_bbH_8(199, 2, 0.),
        log_cs_bbH_13(199, 2, 0.),
        log_cs_ggA_8(199, 2, 0.),
        log_cs_ttA_8(199, 2, 0.),
        log_cs_bbA_8(199, 2, 0.),
        log_cs_ggA_13(199, 2, 0.),
        log_cs_ttA_13(199, 2, 0.),
        log_cs_bbA_13(199, 2, 0.),
        log_cs_ggHp_8(744, 3, 0.),
        log_cs_ggHp_13(1104, 3, 0.),
        csrH_top_8(199, 2, 0.),
        csrH_bottom_8(199, 2, 0.),
        csrA_top_8(199, 2, 0.),
        csrA_bottom_8(199, 2, 0.),
        csrH_top_13(199, 2, 0.),
        csrH_bottom_13(199, 2, 0.),
        csrA_top_13(199, 2, 0.),
        csrA_bottom_13(199, 2, 0.),
        //
        ATLAS8_gg_phi_tautau(92, 2, 0.),
        ATLAS8_bb_phi_tautau(92, 2, 0.),
        ATLAS8_gg_phi_gaga(108, 2, 0.),
        ATLAS8_pp_phi_Zga_llga(141, 2, 0.),
        ATLAS8_gg_phi_ZZ(173, 2, 0.),
        ATLAS8_VV_phi_ZZ(173, 2, 0.),
        ATLAS8_gg_phi_WW(13, 2, 0.),
        ATLAS8_VV_phi_WW(13, 2, 0.),
        ATLAS8_gg_phi_phi1phi1(75, 2, 0.),
        ATLAS8_gg_phi_phi1Z_bbZ(79, 2, 0.),
        ATLAS8_gg_phi_phi1Z_tautauZ(79, 2, 0.),
        CMS8_mu_pp_phi_VV(172, 2, 0.),
        CMS8_bb_phi_bb(81, 2, 0.),
        CMS8_gg_phi_bb(88, 2, 0.),
        CMS8_gg_phi_tautau(92, 2, 0.),
        CMS8_bb_phi_tautau(92, 2, 0.),
        CMS8_pp_phi_Zga_llga(101, 2, 0.),
        CMS8_pp_phi_phi1phi1_bbbb(167, 2, 0.),
        CMS8_pp_phi_phi1phi1_bbgaga(85, 2, 0.),
        CMS8_gg_phi_phi1phi1_bbtautau(10, 2, 0.),
        CMS8_pp_phi_phi1phi1_bbtautau(71, 2, 0.),
        CMS8_gg_phi_phi1Z_bbll(16, 2, 0.),
        CMS8_gg_phi_phi1Z_tautaull(14, 2, 0.),
        CMS8_pp_phii_phijZ_bbll_1(28718, 3, 0.),
        CMS8_pp_phii_phijZ_bbll_2(29050, 3, 0.),
        CMS8_pp_phii_phijZ_tautaull_1(400, 3, 0.),
        CMS8_pp_phii_phijZ_tautaull_2(400, 3, 0.),
        ATLAS13_tt_phi_tt(61, 2, 0.),
        ATLAS13_bb_phi_tt(61, 2, 0.),
        ATLAS13_gg_phi_tautau(206, 2, 0.),
        ATLAS13_bb_phi_tautau(206, 2, 0.),
        ATLAS13_pp_phi_gaga(251, 2, 0.),
        ATLAS13_gg_phi_Zga_llga(216, 2, 0.),
        ATLAS13_gg_phi_Zga_qqga(581, 2, 0.),
        ATLAS13_gg_phi_ZZ_llllnunu(101, 2, 0.),
        ATLAS13_VV_phi_ZZ_llllnunu(101, 2, 0.),
        ATLAS13_gg_phi_ZZ_qqllnunu(271, 2, 0.),
        ATLAS13_VV_phi_ZZ_qqllnunu(271, 2, 0.),
        ATLAS13_gg_phi_WW_enumunu(381, 2, 0.),
        ATLAS13_VV_phi_WW_enumunu(281, 2, 0.),
        ATLAS13_gg_phi_WW_lnuqq(271, 2, 0.),
        ATLAS13_VV_phi_WW_lnuqq(271, 2, 0.),
        ATLAS13_pp_phi_VV_qqqq(181, 2, 0.),
        ATLAS13_pp_phi_phi1phi1_bbbb(275, 2, 0.),
        ATLAS13_pp_phi_phi1phi1_bbgaga(75, 2, 0.),
        ATLAS13_pp_phi_phi1phi1_bbtautau(75, 2, 0.),
        ATLAS13_pp_phi_phi1phi1_bbWW(51, 2, 0.),
        ATLAS13_gg_phi_phi1phi1_gagaWW(25, 2, 0.),
        ATLAS13_gg_phi_phi1Z_bbZ(181, 2, 0.),
        ATLAS13_bb_phi_phi1Z_bbZ(181, 2, 0.),
        ATLAS13_gg_phii_phijZ_bbZ(3364, 3, 0.),
        ATLAS13_bb_phii_phijZ_bbZ(3364, 3, 0.),
        CMS13_pp_phi_bb(66, 2, 0.),
        CMS13_bb_phi_bb(101, 2, 0.),
        CMS13_gg_phi_tautau(312, 2, 0.),
        CMS13_bb_phi_tautau(312, 2, 0.),
        CMS13_gg_phi_gaga(351, 2, 0.),
        CMS13_gg_phi_Zga(366, 2, 0.),
        CMS13_pp_phi_ZZ_llqqnunull(288, 2, 0.),
        CMS13_pp_phi_ZZ_qqnunu(301, 2, 0.),
        CMS13_ggVV_phi_WW_lnulnu(81, 2, 0.),
        CMS13_pp_phi_WW_lnuqq(341, 2, 0.),
        CMS13_pp_phi_phi1phi1_bbbb_1(95, 2, 0.),
        CMS13_pp_phi_phi1phi1_bbbb_2(181, 2, 0.),
        CMS13_pp_phi_phi1phi1_bbgaga(66, 2, 0.),
        CMS13_pp_phi_phi1phi1_bbtautau_1(66, 2, 0.),
        CMS13_pp_phi_phi1phi1_bbtautau_2(311, 2, 0.),
        CMS13_pp_phi_phi1phi1_bbVV(65, 2, 0.),
        CMS13_gg_phi_phi1Z_bbZ_1(79, 2, 0.),
        CMS13_gg_phi_phi1Z_bbZ_2(121, 2, 0.),
        CMS13_bb_phi_phi1Z_bbZ_1(79, 2, 0.),
        CMS13_bb_phi_phi1Z_bbZ_2(121, 2, 0.),
        ATLAS8_pp_Hpm_taunu(83, 2, 0.),
        ATLAS8_pp_Hpm_tb(41, 2, 0.),
        CMS8_pp_Hp_taunu(43, 2, 0.),
        CMS8_pp_Hp_tb(43, 2, 0.),
        ATLAS13_pp_Hpm_taunu(192, 2, 0.),
        ATLAS13_pp_Hpm_tb(181, 2, 0.),
        CMS13_pp_Hpm_taunu(283, 2, 0.),
        arraybsgamma(1111, 3, 0.),
        Mu_GTHDM(3,3,0.), Md_GTHDM(3,3,0.), Ml_GTHDM(3,3,0.),
        Nu_GTHDM(3,3,0.), Nd_GTHDM(3,3,0.), Nl_GTHDM(3,3,0.),
        Yu1_GTHDM(3,3,0.), Yu2_GTHDM(3,3,0.), Yd1_GTHDM(3,3,0.), Yd2_GTHDM(3,3,0.),
        Yl1_GTHDM(3,3,0.), Yl2_GTHDM(3,3,0.),
        //
        myGTHDM(static_cast<const GeneralTHDM*> (&SM_i)), 
        PV(true)
{
    read();
}
GeneralTHDMcache::~GeneralTHDMcache()
{}

/////////////////////////////////////////////////////////////////////////////////////////////////

int GeneralTHDMcache::CacheCheck(const gslpp::complex cache[][CacheSize], 
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

int GeneralTHDMcache::CacheCheckReal(const double cache[][CacheSize], 
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


void GeneralTHDMcache::CacheShift(gslpp::complex cache[][CacheSize], const int NumPar, 
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

void GeneralTHDMcache::CacheShiftReal(double cache[][CacheSize], const int NumPar,
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

////////////////////////////////////////////////////////////////////////////////
/*One-loop functions*/
////////////////////////////////////////////////////////////////////////////////

gslpp::complex GeneralTHDMcache::B0_MZ2_0_MW2_mHh2(const double MZ2, const double MW2, const double mHh2) const {
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

gslpp::complex GeneralTHDMcache::B0_MZ2_0_MW2_mHl2(const double MZ2, const double MW2, const double mHl2) const {
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

gslpp::complex GeneralTHDMcache::B0_MZ2_0_MZ2_mHh2(const double MZ2, const double mHh2) const {
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

gslpp::complex GeneralTHDMcache::B0_MZ2_0_MZ2_mHl2(const double MZ2, const double mHl2) const {
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

gslpp::complex GeneralTHDMcache::B0_MZ2_MW2_MW2_mHh2(const double MZ2, const double MW2, const double mHh2) const {
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

gslpp::complex GeneralTHDMcache::B0_MZ2_MW2_MW2_mHl2(const double MZ2, const double MW2, const double mHl2) const {
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

gslpp::complex GeneralTHDMcache::B0_MZ2_MZ2_MZ2_mHh2(const double MZ2, const double mHh2) const {
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

gslpp::complex GeneralTHDMcache::B0_MZ2_MZ2_MZ2_mHl2(const double MZ2, const double mHl2) const {
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

gslpp::complex GeneralTHDMcache::B0_MZ2_0_0_mHl2(const double MZ2, const double mHl2) const {
    int NumPar = 2;
    double params[] = {MZ2, mHl2};

    int i = CacheCheck(B0_MZ2_0_0_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_0_0_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, 0.0, 0.0, mHl2);
        CacheShift(B0_MZ2_0_0_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_0_0_mHh2(const double MZ2, const double mHh2) const {
    int NumPar = 2;
    double params[] = {MZ2, mHh2};

    int i = CacheCheck(B0_MZ2_0_0_mHh2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_0_0_mHh2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, 0.0, 0.0, mHh2);
        CacheShift(B0_MZ2_0_0_mHh2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_0_mHp2_mHl2(const double MZ2, const double mHp2, const double mHl2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHp2, mHl2};

    int i = CacheCheck(B0_MZ2_0_mHp2_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_0_mHp2_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, 0.0, mHp2, mHl2);
        CacheShift(B0_MZ2_0_mHp2_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_0_mHp2_mHh2(const double MZ2, const double mHp2, const double mHh2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHp2, mHh2};

    int i = CacheCheck(B0_MZ2_0_mHp2_mHh2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_0_mHp2_mHh2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, 0.0, mHp2, mHh2);
        CacheShift(B0_MZ2_0_mHp2_mHh2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_0_mA2_mHl2(const double MZ2, const double mA2, const double mHl2) const {
    int NumPar = 3;
    double params[] = {MZ2, mA2, mHl2};

    int i = CacheCheck(B0_MZ2_0_mA2_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_0_mA2_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, 0.0, mA2, mHl2);
        CacheShift(B0_MZ2_0_mA2_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_0_mA2_mHh2(const double MZ2, const double mA2, const double mHh2) const {
    int NumPar = 3;
    double params[] = {MZ2, mA2, mHh2};

    int i = CacheCheck(B0_MZ2_0_mA2_mHh2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_0_mA2_mHh2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, 0.0, mA2, mHh2);
        CacheShift(B0_MZ2_0_mA2_mHh2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_mHl2_0_0(const double MZ2, const double mHl2) const {
    int NumPar = 2;
    double params[] = {MZ2, mHl2};

    int i = CacheCheck(B0_MZ2_mHl2_0_0_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_mHl2_0_0_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, mHl2, 0.0, 0.0);
        CacheShift(B0_MZ2_mHl2_0_0_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_mHl2_0_mHp2(const double MZ2, const double mHl2, const double mHp2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHl2, mHp2};

    int i = CacheCheck(B0_MZ2_mHl2_0_mHp2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_mHl2_0_mHp2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, mHl2, 0.0, mHp2);
        CacheShift(B0_MZ2_mHl2_0_mHp2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_mHl2_0_mA2(const double MZ2, const double mHl2, const double mA2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHl2, mA2};

    int i = CacheCheck(B0_MZ2_mHl2_0_mA2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_mHl2_0_mA2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, mHl2, 0.0, mA2);
        CacheShift(B0_MZ2_mHl2_0_mA2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_mHl2_mHl2_mHl2(const double MZ2, const double mHl2) const {
    int NumPar = 2;
    double params[] = {MZ2, mHl2};

    int i = CacheCheck(B0_MZ2_mHl2_mHl2_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_mHl2_mHl2_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, mHl2, mHl2, mHl2);
        CacheShift(B0_MZ2_mHl2_mHl2_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_mHl2_mHh2_mHl2(const double MZ2, const double mHl2, const double mHh2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHl2, mHh2};

    int i = CacheCheck(B0_MZ2_mHl2_mHh2_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_mHl2_mHh2_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, mHl2, mHh2, mHl2);
        CacheShift(B0_MZ2_mHl2_mHh2_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_mHl2_mHh2_mHh2(const double MZ2, const double mHl2, const double mHh2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHl2, mHh2};

    int i = CacheCheck(B0_MZ2_mHl2_mHh2_mHh2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_mHl2_mHh2_mHh2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, mHl2, mHh2, mHh2);
        CacheShift(B0_MZ2_mHl2_mHh2_mHh2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_mHl2_mHp2_mHp2(const double MZ2, const double mHl2, const double mHp2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHl2, mHp2};

    int i = CacheCheck(B0_MZ2_mHl2_mHp2_mHp2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_mHl2_mHp2_mHp2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, mHl2, mHp2, mHp2);
        CacheShift(B0_MZ2_mHl2_mHp2_mHp2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_mHl2_mA2_mA2(const double MZ2, const double mHl2, const double mA2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHl2, mA2};

    int i = CacheCheck(B0_MZ2_mHl2_mA2_mA2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_mHl2_mA2_mA2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, mHl2, mA2, mA2);
        CacheShift(B0_MZ2_mHl2_mA2_mA2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_mHh2_0_0(const double MZ2, const double mHh2) const {
    int NumPar = 2;
    double params[] = {MZ2, mHh2};

    int i = CacheCheck(B0_MZ2_mHh2_0_0_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_mHh2_0_0_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, mHh2, 0.0, 0.0);
        CacheShift(B0_MZ2_mHh2_0_0_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_mHh2_0_mHp2(const double MZ2, const double mHh2, const double mHp2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHh2, mHp2};

    int i = CacheCheck(B0_MZ2_mHh2_0_mHp2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_mHh2_0_mHp2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, mHh2, 0.0, mHp2);
        CacheShift(B0_MZ2_mHh2_0_mHp2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_mHh2_0_mA2(const double MZ2, const double mHh2, const double mA2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHh2, mA2};

    int i = CacheCheck(B0_MZ2_mHh2_0_mA2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_mHh2_0_mA2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, mHh2, 0.0, mA2);
        CacheShift(B0_MZ2_mHh2_0_mA2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_mHh2_mHl2_mHl2(const double MZ2, const double mHh2, const double mHl2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHh2, mHl2};

    int i = CacheCheck(B0_MZ2_mHh2_mHl2_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_mHh2_mHl2_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, mHh2, mHl2, mHl2);
        CacheShift(B0_MZ2_mHh2_mHl2_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_mHh2_mHh2_mHl2(const double MZ2, const double mHh2, const double mHl2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHh2, mHl2};

    int i = CacheCheck(B0_MZ2_mHh2_mHh2_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_mHh2_mHh2_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, mHh2, mHh2, mHl2);
        CacheShift(B0_MZ2_mHh2_mHh2_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_mHh2_mHh2_mHh2(const double MZ2, const double mHh2) const {
    int NumPar = 2;
    double params[] = {MZ2, mHh2};

    int i = CacheCheck(B0_MZ2_mHh2_mHh2_mHh2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_mHh2_mHh2_mHh2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, mHh2, mHh2, mHh2);
        CacheShift(B0_MZ2_mHh2_mHh2_mHh2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_mHh2_mHp2_mHp2(const double MZ2, const double mHh2, const double mHp2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHh2, mHp2};

    int i = CacheCheck(B0_MZ2_mHh2_mHp2_mHp2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_mHh2_mHp2_mHp2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, mHh2, mHp2, mHp2);
        CacheShift(B0_MZ2_mHh2_mHp2_mHp2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_mHh2_mA2_mA2(const double MZ2, const double mHh2, const double mA2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHh2, mA2};

    int i = CacheCheck(B0_MZ2_mHh2_mA2_mA2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_mHh2_mA2_mA2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, mHh2, mA2, mA2);
        CacheShift(B0_MZ2_mHh2_mA2_mA2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_mHp2_0_mHl2(const double MZ2, const double mHp2, const double mHl2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHp2, mHl2};

    int i = CacheCheck(B0_MZ2_mHp2_0_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_mHp2_0_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, mHp2, 0.0, mHl2);
        CacheShift(B0_MZ2_mHp2_0_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_mHp2_0_mHh2(const double MZ2, const double mHp2, const double mHh2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHp2, mHh2};

    int i = CacheCheck(B0_MZ2_mHp2_0_mHh2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_mHp2_0_mHh2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, mHp2, 0.0, mHh2);
        CacheShift(B0_MZ2_mHp2_0_mHh2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_mHp2_mHp2_mHl2(const double MZ2, const double mHp2, const double mHl2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHp2, mHl2};

    int i = CacheCheck(B0_MZ2_mHp2_mHp2_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_mHp2_mHp2_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, mHp2, mHp2, mHl2);
        CacheShift(B0_MZ2_mHp2_mHp2_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_mHp2_mHp2_mHh2(const double MZ2, const double mHp2, const double mHh2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHp2, mHh2};

    int i = CacheCheck(B0_MZ2_mHp2_mHp2_mHh2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_mHp2_mHp2_mHh2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, mHp2, mHp2, mHh2);
        CacheShift(B0_MZ2_mHp2_mHp2_mHh2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_mA2_0_mHl2(const double MZ2, const double mA2, const double mHl2) const {
    int NumPar = 3;
    double params[] = {MZ2, mA2, mHl2};

    int i = CacheCheck(B0_MZ2_mA2_0_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_mA2_0_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, mA2, 0.0, mHl2);
        CacheShift(B0_MZ2_mA2_0_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_mA2_0_mHh2(const double MZ2, const double mA2, const double mHh2) const {
    int NumPar = 3;
    double params[] = {MZ2, mA2, mHh2};

    int i = CacheCheck(B0_MZ2_mA2_0_mHh2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_mA2_0_mHh2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, mA2, 0.0, mHh2);
        CacheShift(B0_MZ2_mA2_0_mHh2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_mA2_mA2_mHl2(const double MZ2, const double mA2, const double mHl2) const {
    int NumPar = 3;
    double params[] = {MZ2, mA2, mHl2};

    int i = CacheCheck(B0_MZ2_mA2_mA2_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_mA2_mA2_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, mA2, mA2, mHl2);
        CacheShift(B0_MZ2_mA2_mA2_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0_MZ2_mA2_mA2_mHh2(const double MZ2, const double mA2, const double mHh2) const {
    int NumPar = 3;
    double params[] = {MZ2, mA2, mHh2};

    int i = CacheCheck(B0_MZ2_mA2_mA2_mHh2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_mA2_mA2_mHh2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2, mA2, mA2, mHh2);
        CacheShift(B0_MZ2_mA2_mA2_mHh2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

///////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////

gslpp::complex GeneralTHDMcache::B0p_MZ2_0_0_mHl2(const double MZ2, const double mHl2) const {
    int NumPar = 2;
    double params[] = {MZ2, mHl2};

    int i = CacheCheck(B0p_MZ2_0_0_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_0_0_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, 0.0, 0.0, mHl2);
        CacheShift(B0p_MZ2_0_0_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_0_0_mHh2(const double MZ2, const double mHh2) const {
    int NumPar = 2;
    double params[] = {MZ2, mHh2};

    int i = CacheCheck(B0p_MZ2_0_0_mHh2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_0_0_mHh2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, 0.0, 0.0, mHh2);
        CacheShift(B0p_MZ2_0_0_mHh2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_0_mHp2_mHl2(const double MZ2, const double mHp2, const double mHl2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHp2, mHl2};

    int i = CacheCheck(B0p_MZ2_0_mHp2_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_0_mHp2_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, 0.0, mHp2, mHl2);
        CacheShift(B0p_MZ2_0_mHp2_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_0_mHp2_mHh2(const double MZ2, const double mHp2, const double mHh2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHp2, mHh2};

    int i = CacheCheck(B0p_MZ2_0_mHp2_mHh2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_0_mHp2_mHh2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, 0.0, mHp2, mHh2);
        CacheShift(B0p_MZ2_0_mHp2_mHh2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_0_mHp2_mA2(const double MZ2, const double mHp2, const double mA2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHp2, mA2};

    int i = CacheCheck(B0p_MZ2_0_mHp2_mA2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_0_mHp2_mA2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, 0.0, mHp2, mA2);
        CacheShift(B0p_MZ2_0_mHp2_mA2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_0_mA2_mHl2(const double MZ2, const double mA2, const double mHl2) const {
    int NumPar = 3;
    double params[] = {MZ2, mA2, mHl2};

    int i = CacheCheck(B0p_MZ2_0_mA2_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_0_mA2_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, 0.0, mA2, mHl2);
        CacheShift(B0p_MZ2_0_mA2_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_0_mA2_mHh2(const double MZ2, const double mA2, const double mHh2) const {
    int NumPar = 3;
    double params[] = {MZ2, mA2, mHh2};

    int i = CacheCheck(B0p_MZ2_0_mA2_mHh2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_0_mA2_mHh2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, 0.0, mA2, mHh2);
        CacheShift(B0p_MZ2_0_mA2_mHh2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_mHl2_0_0(const double MZ2, const double mHl2) const {
    int NumPar = 2;
    double params[] = {MZ2, mHl2};

    int i = CacheCheck(B0p_MZ2_0_0_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_0_0_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, mHl2, 0.0, 0.0);
        CacheShift(B0p_MZ2_0_0_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_mHl2_0_mHp2(const double MZ2, const double mHl2, const double mHp2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHl2, mHp2};

    int i = CacheCheck(B0p_MZ2_mHl2_0_mHp2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_mHl2_0_mHp2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, mHl2, 0.0, mHp2);
        CacheShift(B0p_MZ2_mHl2_0_mHp2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_mHl2_0_mA2(const double MZ2, const double mHl2, const double mA2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHl2, mA2};

    int i = CacheCheck(B0p_MZ2_mHl2_0_mA2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_mHl2_0_mA2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, mHl2, 0.0, mA2);
        CacheShift(B0p_MZ2_mHl2_0_mA2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_mHl2_mHl2_mHl2(const double MZ2, const double mHl2) const {
    int NumPar = 2;
    double params[] = {MZ2, mHl2};

    int i = CacheCheck(B0p_MZ2_mHl2_mHl2_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_mHl2_mHl2_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, mHl2, mHl2, mHl2);
        CacheShift(B0p_MZ2_mHl2_mHl2_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_mHl2_mHh2_mHl2(const double MZ2, const double mHl2, const double mHh2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHl2, mHh2};

    int i = CacheCheck(B0p_MZ2_mHl2_mHh2_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_mHl2_mHh2_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, mHl2, mHh2, mHl2);
        CacheShift(B0p_MZ2_mHl2_mHh2_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_mHl2_mHh2_mHh2(const double MZ2, const double mHl2, const double mHh2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHl2, mHh2};

    int i = CacheCheck(B0p_MZ2_mHl2_mHh2_mHh2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_mHl2_mHh2_mHh2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, mHl2, mHh2, mHh2);
        CacheShift(B0p_MZ2_mHl2_mHh2_mHh2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_mHl2_mHp2_mHp2(const double MZ2, const double mHl2, const double mHp2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHl2, mHp2};

    int i = CacheCheck(B0p_MZ2_mHl2_mHp2_mHp2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_mHl2_mHp2_mHp2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, mHl2, mHp2, mHp2);
        CacheShift(B0p_MZ2_mHl2_mHp2_mHp2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_mHl2_mA2_mA2(const double MZ2, const double mHl2, const double mA2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHl2, mA2};

    int i = CacheCheck(B0p_MZ2_mHl2_mA2_mA2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_mHl2_mA2_mA2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, mHl2, mA2, mA2);
        CacheShift(B0p_MZ2_mHl2_mA2_mA2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_mHh2_0_0(const double MZ2, const double mHh2) const {
    int NumPar = 2;
    double params[] = {MZ2, mHh2};

    int i = CacheCheck(B0p_MZ2_0_0_mHh2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_0_0_mHh2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, mHh2, 0.0, 0.0);
        CacheShift(B0p_MZ2_0_0_mHh2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_mHh2_0_mHp2(const double MZ2, const double mHh2, const double mHp2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHh2, mHp2};

    int i = CacheCheck(B0p_MZ2_mHh2_0_mHp2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_mHh2_0_mHp2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, mHh2, 0.0, mHp2);
        CacheShift(B0p_MZ2_mHh2_0_mHp2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_mHh2_0_mA2(const double MZ2, const double mHh2, const double mA2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHh2, mA2};

    int i = CacheCheck(B0p_MZ2_mHh2_0_mA2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_mHh2_0_mA2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, mHh2, 0.0, mA2);
        CacheShift(B0p_MZ2_mHh2_0_mA2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_mHh2_mHl2_mHl2(const double MZ2, const double mHh2, const double mHl2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHh2, mHl2};

    int i = CacheCheck(B0p_MZ2_mHh2_mHl2_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_mHh2_mHl2_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, mHh2, mHl2, mHl2);
        CacheShift(B0p_MZ2_mHh2_mHl2_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_mHh2_mHh2_mHl2(const double MZ2, const double mHh2, const double mHl2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHh2, mHl2};

    int i = CacheCheck(B0p_MZ2_mHh2_mHh2_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_mHh2_mHh2_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, mHh2, mHh2, mHl2);
        CacheShift(B0p_MZ2_mHh2_mHh2_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_mHh2_mHh2_mHh2(const double MZ2, const double mHh2) const {
    int NumPar = 2;
    double params[] = {MZ2, mHh2};

    int i = CacheCheck(B0p_MZ2_mHh2_mHh2_mHh2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_mHh2_mHh2_mHh2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, mHh2, mHh2, mHh2);
        CacheShift(B0p_MZ2_mHh2_mHh2_mHh2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_mHh2_mHp2_mHp2(const double MZ2, const double mHh2, const double mHp2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHh2, mHp2};

    int i = CacheCheck(B0p_MZ2_mHh2_mHp2_mHp2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_mHh2_mHp2_mHp2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, mHh2, mHp2, mHp2);
        CacheShift(B0p_MZ2_mHh2_mHp2_mHp2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_mHh2_mA2_mA2(const double MZ2, const double mHh2, const double mA2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHh2, mA2};

    int i = CacheCheck(B0p_MZ2_mHh2_mA2_mA2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_mHh2_mA2_mA2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, mHh2, mA2, mA2);
        CacheShift(B0p_MZ2_mHh2_mA2_mA2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_mHp2_0_mHl2(const double MZ2, const double mHp2, const double mHl2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHp2, mHl2};

    int i = CacheCheck(B0p_MZ2_mHp2_0_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_mHp2_0_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, mHp2, 0.0, mHl2);
        CacheShift(B0p_MZ2_mHp2_0_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_mHp2_0_mHh2(const double MZ2, const double mHp2, const double mHh2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHp2, mHh2};

    int i = CacheCheck(B0p_MZ2_mHp2_0_mHh2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_mHp2_0_mHh2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, mHp2, 0.0, mHh2);
        CacheShift(B0p_MZ2_mHp2_0_mHh2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_mHp2_0_mA2(const double MZ2, const double mHp2, const double mA2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHp2, mA2};

    int i = CacheCheck(B0p_MZ2_mHp2_0_mA2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_mHp2_0_mA2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, mHp2, 0.0, mA2);
        CacheShift(B0p_MZ2_mHp2_0_mA2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_mHp2_mHp2_mHl2(const double MZ2, const double mHp2, const double mHl2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHp2, mHl2};

    int i = CacheCheck(B0p_MZ2_mHp2_mHp2_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_mHp2_mHp2_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, mHp2, mHp2, mHl2);
        CacheShift(B0p_MZ2_mHp2_mHp2_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_mHp2_mHp2_mHh2(const double MZ2, const double mHp2, const double mHh2) const {
    int NumPar = 3;
    double params[] = {MZ2, mHp2, mHh2};

    int i = CacheCheck(B0p_MZ2_mHp2_mHp2_mHh2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_mHp2_mHp2_mHh2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, mHp2, mHp2, mHh2);
        CacheShift(B0p_MZ2_mHp2_mHp2_mHh2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_mA2_0_mHl2(const double MZ2, const double mA2, const double mHl2) const {
    int NumPar = 3;
    double params[] = {MZ2, mA2, mHl2};

    int i = CacheCheck(B0p_MZ2_mA2_0_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_mA2_0_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, mA2, 0.0, mHl2);
        CacheShift(B0p_MZ2_mA2_0_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_mA2_0_mHh2(const double MZ2, const double mA2, const double mHh2) const {
    int NumPar = 3;
    double params[] = {MZ2, mA2, mHh2};

    int i = CacheCheck(B0p_MZ2_mA2_0_mHh2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_mA2_0_mHh2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, mA2, 0.0, mHh2);
        CacheShift(B0p_MZ2_mA2_0_mHh2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_mA2_0_mHp2(const double MZ2, const double mA2, const double mHp2) const {
    int NumPar = 3;
    double params[] = {MZ2, mA2, mHp2};

    int i = CacheCheck(B0p_MZ2_mA2_0_mHp2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_mA2_0_mHp2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, mA2, 0.0, mHp2);
        CacheShift(B0p_MZ2_mA2_0_mHp2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_mA2_mA2_mHl2(const double MZ2, const double mA2, const double mHl2) const {
    int NumPar = 3;
    double params[] = {MZ2, mA2, mHl2};

    int i = CacheCheck(B0p_MZ2_mA2_mA2_mHl2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_mA2_mA2_mHl2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, mA2, mA2, mHl2);
        CacheShift(B0p_MZ2_mA2_mA2_mHl2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex GeneralTHDMcache::B0p_MZ2_mA2_mA2_mHh2(const double MZ2, const double mA2, const double mHh2) const {
    int NumPar = 3;
    double params[] = {MZ2, mA2, mHh2};

    int i = CacheCheck(B0p_MZ2_mA2_mA2_mHh2_cache, NumPar, params);
    if (i>=0) {
        return ( B0p_MZ2_mA2_mA2_mHh2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0p(MZ2, mA2, mA2, mHh2);
        CacheShift(B0p_MZ2_mA2_mA2_mHh2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

///////////////////////////////////////////////////////////////////////////////////////////

gslpp::complex GeneralTHDMcache::B00_MZ2_0_mA2_mHp2(const double MZ2, const double mA2, const double mHp2) const {
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

gslpp::complex GeneralTHDMcache::B00_MZ2_0_mHh2_mA2(const double MZ2, const double mHh2, const double mA2) const {
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

gslpp::complex GeneralTHDMcache::B00_MZ2_0_mHh2_mHp2(const double MZ2, const double mHh2, const double mHp2) const {
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

gslpp::complex GeneralTHDMcache::B00_MZ2_0_mHl2_mA2(const double MZ2, const double mHl2, const double mA2) const {
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

gslpp::complex GeneralTHDMcache::B00_MZ2_0_mHl2_mHp2(const double MZ2, const double mHl2, const double mHp2) const {
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

gslpp::complex GeneralTHDMcache::B00_MZ2_0_mHp2_mHp2(const double MZ2, const double mHp2) const {
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

gslpp::complex GeneralTHDMcache::B00_MZ2_0_MW2_mHh2(const double MZ2, const double MW2, const double mHh2) const {
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

gslpp::complex GeneralTHDMcache::B00_MZ2_0_MW2_mHl2(const double MZ2, const double MW2, const double mHl2) const {
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

gslpp::complex GeneralTHDMcache::B00_MZ2_0_MZ2_mHh2(const double MZ2, const double mHh2) const {
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

gslpp::complex GeneralTHDMcache::B00_MZ2_0_MZ2_mHl2(const double MZ2, const double mHl2) const {
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

gslpp::complex GeneralTHDMcache::B00_MZ2_MW2_mA2_mHp2(const double MZ2, const double MW2, const double mA2, const double mHp2) const {
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

gslpp::complex GeneralTHDMcache::B00_MZ2_MW2_mHh2_mHp2(const double MZ2, const double MW2, const double mHh2, const double mHp2) const {
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

gslpp::complex GeneralTHDMcache::B00_MZ2_MW2_mHl2_mHp2(const double MZ2, const double MW2, const double mHl2, const double mHp2) const {
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

gslpp::complex GeneralTHDMcache::B00_MZ2_MW2_mHp2_mHp2(const double MZ2, const double MW2, const double mHp2) const {
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

gslpp::complex GeneralTHDMcache::B00_MZ2_MW2_MW2_mHh2(const double MZ2, const double MW2, const double mHh2) const {
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

gslpp::complex GeneralTHDMcache::B00_MZ2_MW2_MW2_mHl2(const double MZ2, const double MW2, const double mHl2) const {
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

gslpp::complex GeneralTHDMcache::B00_MZ2_MZ2_mHh2_mA2(const double MZ2, const double mHh2, const double mA2) const {
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

gslpp::complex GeneralTHDMcache::B00_MZ2_MZ2_mHl2_mA2(const double MZ2, const double mHl2, const double mA2) const {
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

gslpp::complex GeneralTHDMcache::B00_MZ2_MZ2_mHp2_mHp2(const double MZ2, const double mHp2) const {
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

gslpp::complex GeneralTHDMcache::B00_MZ2_MZ2_MZ2_mHh2(const double MZ2, const double mHh2) const {
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

gslpp::complex GeneralTHDMcache::B00_MZ2_MZ2_MZ2_mHl2(const double MZ2, const double mHl2) const {
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


////////////////////////////////////////////////////////////////////////////////

gslpp::matrix<double> GeneralTHDMcache::readTable(std::string filename, int rowN, int colN){

    std::ifstream INfile;
    std::string lineTab;
    INfile.open( filename.c_str() );
    if(INfile.fail()){
        std::cout<<"error: in GeneralTHDMcache, table doesn't exist!"<<std::endl;
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

double GeneralTHDMcache::interpolate(gslpp::matrix<double> arrayTab, double x){

    int rowN=arrayTab.size_i();
    
    double xmin = arrayTab(0,0);
    double xmax = arrayTab(rowN-1,0);
    double interval = arrayTab(1,0)-arrayTab(0,0);
    int Nintervals = (x-xmin)/interval;
    double y = 0.0;
       
    if(x<xmin){
//        std::cout<<"warning: your table parameter value is smaller than the minimum allowed value"<<std::endl;
        return 0.;
    }
    else if(x>xmax){
//        std::cout<<"warning: your table parameter value is greater than the maximum allowed value"<<std::endl;
        return 0.;
    }
    else{
        y =(arrayTab(Nintervals+1,1)-arrayTab(Nintervals,1))/(arrayTab(Nintervals+1,0)
                   -arrayTab(Nintervals,0))*(x-arrayTab(Nintervals,0))+arrayTab(Nintervals,1);
        return y;
    }
}

//2D interpolation

double GeneralTHDMcache::interpolate2D(gslpp::matrix<double> arrayTab, double x, double y){

    int rowN=arrayTab.size_i();

    double xmin = arrayTab(0,0);
    double xmax = arrayTab(rowN-1,0);
    double ymin = arrayTab(0,1);
    double ymax = arrayTab(rowN-1,1);
    double intervalx = arrayTab(1,0)-arrayTab(0,0);
    int i=1;
    do i++;
    while(arrayTab(i,1)-arrayTab(i-1,1)==0&&i<30000);
    double intervaly = arrayTab(i,1)-arrayTab(i-1,1);
    int Nintervalsx = (x-xmin)/intervalx;
    int Nintervalsy = (y-ymin)/intervaly;
    if(x<xmin||x>xmax||y<ymin||y>ymax){
//        std::cout<<"warning: the parameter point lies outside the table"<<std::endl;
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

void GeneralTHDMcache::read(){
  std::stringstream br1,br2,br3,br4,br5,br6,br7;
    std::stringstream dw1;
    std::stringstream cs1,cs2,cs3,cs4,cs5,cs6,cs7,cs8,cs9;
    std::stringstream cs11,cs12,cs13,cs14,cs15,cs16,cs17,cs18,cs19;
    std::stringstream cs20,cs21;
    std::stringstream csr1,csr2,csr3,csr4;
    std::stringstream csr11,csr12,csr13,csr14;
    std::stringstream ex1,ex2,ex3,ex4,ex5,ex6,ex7,ex8,ex9,ex10,ex11,ex12,ex13,ex14,ex15,ex16,ex17,ex18,ex19,ex20,ex21,ex22,ex23;
    std::stringstream ex24,ex25,ex26,ex27,ex28,ex29,ex30,ex31,ex32,ex33,ex34,ex35,ex36,ex37,ex38,ex39,ex40,ex41,ex42,ex43,ex44,\
            ex45,ex46,ex46a,ex47,ex48,ex49,ex50,ex51,ex52,ex53,ex54,ex55,ex56;
    std::stringstream ex57,ex58,ex59,ex60,ex61,ex62,ex63,ex64,ex65,ex66,ex67,ex68,ex69,ex70,ex71,ex72,ex73,ex74,ex75,ex76,ex77,\
            ex78;//,ex79,ex80,ex81,ex82,ex83,ex84,ex85,ex86,ex87,ex88,ex89,ex90,ex91,ex92,ex93,ex94,ex95,ex96,ex97,ex98
    std::stringstream bsg1;

       std::cout<<"reading tables"<<std::endl;

//    std::cout << "HEPFITTABS = " << getenv("HEPFITPATH") << std::endl;
    std::stringstream path;
    path << getenv("HEPFITTABS") << "/THDM/tabs/";
    std::string tablepath=path.str();

    br1 << tablepath << "br1.dat";
    br_tt = readTable(br1.str(),19961,2);
    br2 << tablepath << "br2.dat";
    br_bb = readTable(br2.str(),19961,2);
    br3 << tablepath << "br3.dat";
    br_tautau = readTable(br3.str(),19961,2); 
    br4 << tablepath << "br4.dat";
    br_cc = readTable(br4.str(),19961,2);
    br5 << tablepath << "br5.dat";
    br_mumu = readTable(br5.str(),19961,2);
    br6 << tablepath << "br6.dat";
    br_ZZ = readTable(br6.str(),19961,2);
    br7 << tablepath << "br7.dat";
    br_WW = readTable(br7.str(),19961,2);
    dw1 << tablepath << "dw1.dat";
    GammaHtot_SM = readTable(dw1.str(),19961,2);
    cs1 << tablepath << "log_cs_ggH_8.dat";
    log_cs_ggH_8 = readTable(cs1.str(),199,2);
    cs11 << tablepath << "log_cs_ggH_13.dat";
    log_cs_ggH_13 = readTable(cs11.str(),199,2);
    cs2 << tablepath << "log_cs_VBF_8.dat";
    log_cs_VBF_8 = readTable(cs2.str(),199,2);
    cs12 << tablepath << "log_cs_VBF_13.dat";
    log_cs_VBF_13 = readTable(cs12.str(),199,2);
    cs3 << tablepath << "log_cs_WH_8.dat";
    log_cs_WH_8 = readTable(cs3.str(),199,2);
    cs13 << tablepath << "log_cs_WH_13.dat";
    log_cs_WH_13 = readTable(cs13.str(),199,2);
    cs4 << tablepath << "log_cs_ZH_8.dat";
    log_cs_ZH_8 = readTable(cs4.str(),199,2);
    cs14 << tablepath << "log_cs_ZH_13.dat";
    log_cs_ZH_13 = readTable(cs14.str(),199,2);
    cs5 << tablepath << "log_cs_ttH_8.dat";
    log_cs_ttH_8 = readTable(cs5.str(),199,2);
    cs15 << tablepath << "log_cs_ttH_13.dat";
    log_cs_ttH_13 = readTable(cs15.str(),199,2);
    cs6 << tablepath << "log_cs_bbH_8.dat";
    log_cs_bbH_8 = readTable(cs6.str(),199,2);
    cs16 << tablepath << "log_cs_bbH_13.dat";
    log_cs_bbH_13 = readTable(cs16.str(),199,2);
    cs7 << tablepath << "log_cs_ggA_8.dat";
    log_cs_ggA_8 = readTable(cs7.str(),199,2);
    cs17 << tablepath << "log_cs_ggA_13.dat";
    log_cs_ggA_13 = readTable(cs17.str(),199,2);
    cs8 << tablepath << "log_cs_ttA_8.dat";
    log_cs_ttA_8 = readTable(cs8.str(),199,2);
    cs18 << tablepath << "log_cs_ttA_13.dat";
    log_cs_ttA_13 = readTable(cs18.str(),199,2);
    cs9 << tablepath << "log_cs_bbA_8.dat";
    log_cs_bbA_8 = readTable(cs9.str(),199,2);
    cs19 << tablepath << "log_cs_bbA_13.dat";
    log_cs_bbA_13 = readTable(cs19.str(),199,2);
    cs20 << tablepath << "log_cs_ggHp_8.dat";
    log_cs_ggHp_8 = readTable(cs20.str(),744,3);
    cs21 << tablepath << "log_cs_ggHp_13.dat";
    log_cs_ggHp_13 = readTable(cs21.str(),1104,3);
    csr1 << tablepath << "csrH_top_8.dat";
    csrH_top_8 = readTable(csr1.str(),199,2);
    csr11 << tablepath << "csrH_top_13.dat";
    csrH_top_13 = readTable(csr11.str(),199,2);
    csr2 << tablepath << "csrH_bottom_8.dat";
    csrH_bottom_8 = readTable(csr2.str(),199,2);
    csr12 << tablepath << "csrH_bottom_13.dat";
    csrH_bottom_13 = readTable(csr12.str(),199,2);
    csr3 << tablepath << "csrA_top_8.dat";
    csrA_top_8 = readTable(csr3.str(),199,2);
    csr13 << tablepath << "csrA_top_13.dat";
    csrA_top_13 = readTable(csr13.str(),199,2);
    csr4 << tablepath << "csrA_bottom_8.dat";
    csrA_bottom_8 = readTable(csr4.str(),199,2);
    csr14 << tablepath << "csrA_bottom_13.dat";
    csrA_bottom_13 = readTable(csr14.str(),199,2);

    ex1 << tablepath << "14096064_a.dat";
    ATLAS8_gg_phi_tautau = readTable(ex1.str(),92,2);
    ex2 << tablepath << "14096064_b.dat";
    ATLAS8_bb_phi_tautau = readTable(ex2.str(),92,2);
    ex3 << tablepath << "14076583.dat";
    ATLAS8_gg_phi_gaga = readTable(ex3.str(),108,2);
    ex4 << tablepath << "14078150.dat";
    ATLAS8_pp_phi_Zga_llga = readTable(ex4.str(),141,2);
    ex5 << tablepath << "150705930_a.dat";
    ATLAS8_gg_phi_ZZ = readTable(ex5.str(),173,2);
    ex6 << tablepath << "150705930_b.dat";
    ATLAS8_VV_phi_ZZ = readTable(ex6.str(),173,2);
    ex7 << tablepath << "150900389_a.dat";
    ATLAS8_gg_phi_WW = readTable(ex7.str(),13,2);
    ex8 << tablepath << "150900389_b.dat";
    ATLAS8_VV_phi_WW = readTable(ex8.str(),13,2);
    ex9 << tablepath << "150904670.dat";
    ATLAS8_gg_phi_phi1phi1 = readTable(ex9.str(),75,2);
    ex10 << tablepath << "150204478_b.dat";
    ATLAS8_gg_phi_phi1Z_bbZ = readTable(ex10.str(),79,2);
    ex11 << tablepath << "150204478_a.dat";
    ATLAS8_gg_phi_phi1Z_tautauZ = readTable(ex11.str(),79,2);
    ex12 << tablepath << "150400936.dat";
    CMS8_mu_pp_phi_VV = readTable(ex12.str(),172,2);
    ex13 << tablepath << "150608329.dat";
    CMS8_bb_phi_bb = readTable(ex13.str(),81,2);
    ex14 << tablepath << "180206149.dat";
    CMS8_gg_phi_bb = readTable(ex14.str(),88,2);
    ex15 << tablepath << "CMS-PAS-HIG-14-029_a.dat";
    CMS8_gg_phi_tautau = readTable(ex15.str(),92,2);
    ex16 << tablepath << "CMS-PAS-HIG-14-029_b.dat";
    CMS8_bb_phi_tautau = readTable(ex16.str(),92,2);
    ex17 << tablepath << "CMS-PAS-HIG-16-014.dat";
    CMS8_pp_phi_Zga_llga = readTable(ex17.str(),101,2);
    ex18 << tablepath << "150304114.dat";
    CMS8_pp_phi_phi1phi1_bbbb = readTable(ex18.str(),167,2);
    ex19 << tablepath << "160306896.dat";
    CMS8_pp_phi_phi1phi1_bbgaga = readTable(ex19.str(),85,2);
    ex20 << tablepath << "151001181_a.dat";
    CMS8_gg_phi_phi1phi1_bbtautau = readTable(ex20.str(),10,2);
    ex21 << tablepath << "170700350.dat";
    CMS8_pp_phi_phi1phi1_bbtautau = readTable(ex21.str(),71,2);
    ex22 << tablepath << "150404710.dat";
    CMS8_gg_phi_phi1Z_bbll = readTable(ex22.str(),16,2);
    ex23 << tablepath << "151001181_b.dat";
    CMS8_gg_phi_phi1Z_tautaull = readTable(ex23.str(),14,2);

    ex24 << tablepath << "160302991_a.dat";
    CMS8_pp_phii_phijZ_bbll_1 = readTable(ex24.str(),28718,3);
    ex25 << tablepath << "160302991_b.dat";
    CMS8_pp_phii_phijZ_bbll_2 = readTable(ex25.str(),29050,3);
    ex26 << tablepath << "160302991_c.dat";
    CMS8_pp_phii_phijZ_tautaull_1 = readTable(ex26.str(),400,3);
    ex27 << tablepath << "160302991_d.dat";
    CMS8_pp_phii_phijZ_tautaull_2 = readTable(ex27.str(),400,3);

    ex28 << tablepath << "180711883.dat";
    ATLAS13_tt_phi_tt = readTable(ex28.str(),61,2);
    ex29 << tablepath << "ATLAS-CONF-2016-104_b.dat";
    ATLAS13_bb_phi_tt = readTable(ex29.str(),61,2);
    ex30 << tablepath << "170907242_a.dat";
    ATLAS13_gg_phi_tautau = readTable(ex30.str(),206,2);
    ex31 << tablepath << "170907242_b.dat";
    ATLAS13_bb_phi_tautau = readTable(ex31.str(),206,2);
    ex32 << tablepath << "170704147.dat";
    ATLAS13_pp_phi_gaga = readTable(ex32.str(),251,2);
    ex33 << tablepath << "170800212.dat";
    ATLAS13_gg_phi_Zga_llga = readTable(ex33.str(),216,2);
    ex34 << tablepath << "180501908.dat";
    ATLAS13_gg_phi_Zga_qqga = readTable(ex34.str(),581,2);
    ex35 << tablepath << "171206386_a.dat";
    ATLAS13_gg_phi_ZZ_llllnunu = readTable(ex35.str(),101,2);
    ex36 << tablepath << "171206386_b.dat";
    ATLAS13_VV_phi_ZZ_llllnunu = readTable(ex36.str(),101,2);
    ex37 << tablepath << "170809638_a.dat";
    ATLAS13_gg_phi_ZZ_qqllnunu = readTable(ex37.str(),271,2);
    ex38 << tablepath << "170809638_b.dat";
    ATLAS13_VV_phi_ZZ_qqllnunu = readTable(ex38.str(),271,2);
    ex39 << tablepath << "171001123_a.dat";
    ATLAS13_gg_phi_WW_enumunu = readTable(ex39.str(),381,2);
    ex40 << tablepath << "171001123_b.dat";
    ATLAS13_VV_phi_WW_enumunu = readTable(ex40.str(),281,2);
    ex41 << tablepath << "171007235_a.dat";
    ATLAS13_gg_phi_WW_lnuqq = readTable(ex41.str(),271,2);
    ex42 << tablepath << "171007235_b.dat";
    ATLAS13_VV_phi_WW_lnuqq = readTable(ex42.str(),271,2);
    ex43 << tablepath << "170804445.dat";
    ATLAS13_pp_phi_VV_qqqq = readTable(ex43.str(),181,2);
    ex44 << tablepath << "180406174.dat";
    ATLAS13_pp_phi_phi1phi1_bbbb = readTable(ex44.str(),275,2);
    ex45 << tablepath << "180704873.dat";
    ATLAS13_pp_phi_phi1phi1_bbgaga = readTable(ex45.str(),75,2);
    ex46 << tablepath << "180800336.dat";
    ATLAS13_pp_phi_phi1phi1_bbtautau = readTable(ex46.str(),75,2);
    ex46a << tablepath << "181104671.dat";
    ATLAS13_pp_phi_phi1phi1_bbWW = readTable(ex46a.str(),51,2);
    ex47 << tablepath << "180708567.dat";
    ATLAS13_gg_phi_phi1phi1_gagaWW = readTable(ex47.str(),25,2);
    ex48 << tablepath << "171206518_a.dat";
    ATLAS13_gg_phi_phi1Z_bbZ = readTable(ex48.str(),181,2);
    ex49 << tablepath << "171206518_b.dat";
    ATLAS13_bb_phi_phi1Z_bbZ = readTable(ex49.str(),181,2);

    ex50 << tablepath << "180401126_a.dat";
    ATLAS13_gg_phii_phijZ_bbZ = readTable(ex50.str(),3364,3);
    ex51 << tablepath << "180401126_b.dat";
    ATLAS13_bb_phii_phijZ_bbZ = readTable(ex51.str(),3364,3);

    ex52 << tablepath << "CMS-PAS-HIG-16-025.dat";
    CMS13_pp_phi_bb = readTable(ex52.str(),66,2);
    ex53 << tablepath << "180512191.dat";
    CMS13_bb_phi_bb = readTable(ex53.str(),101,2);
    ex54 << tablepath << "180306553_a.dat";
    CMS13_gg_phi_tautau = readTable(ex54.str(),312,2);
    ex55 << tablepath << "180306553_b.dat";
    CMS13_bb_phi_tautau = readTable(ex55.str(),312,2);
    ex56 << tablepath << "160902507.dat";
    CMS13_gg_phi_gaga = readTable(ex56.str(),351,2);
    ex57 << tablepath << "171203143.dat";
    CMS13_gg_phi_Zga = readTable(ex57.str(),366,2);
    ex58 << tablepath << "180401939_a.dat";
    CMS13_pp_phi_ZZ_llqqnunull = readTable(ex58.str(),288,2);
    ex59 << tablepath << "180303838.dat";
    CMS13_pp_phi_ZZ_qqnunu = readTable(ex59.str(),301,2);
    ex60 << tablepath << "CMS-PAS-HIG-16-023.dat";
    CMS13_ggVV_phi_WW_lnulnu = readTable(ex60.str(),81,2);
    ex61 << tablepath << "180209407.dat";
    CMS13_pp_phi_WW_lnuqq = readTable(ex61.str(),341,2);
    ex62 << tablepath << "180603548.dat";
    CMS13_pp_phi_phi1phi1_bbbb_1 = readTable(ex62.str(),95,2);
    ex63 << tablepath << "180801473.dat";
    CMS13_pp_phi_phi1phi1_bbbb_2 = readTable(ex63.str(),181,2);
    ex64 << tablepath << "180600408.dat";
    CMS13_pp_phi_phi1phi1_bbgaga = readTable(ex64.str(),66,2);
    ex65 << tablepath << "170702909.dat";
    CMS13_pp_phi_phi1phi1_bbtautau_1 = readTable(ex65.str(),66,2);
    ex66 << tablepath << "180801365.dat";
    CMS13_pp_phi_phi1phi1_bbtautau_2 = readTable(ex66.str(),311,2);
    ex67 << tablepath << "170804188.dat";
    CMS13_pp_phi_phi1phi1_bbVV = readTable(ex67.str(),65,2);
    ex68 << tablepath << "CMS-PAS-HIG-18-005_a.dat";
    CMS13_gg_phi_phi1Z_bbZ_1 = readTable(ex68.str(),79,2);
    ex69 << tablepath << "180702826_a.dat";
    CMS13_gg_phi_phi1Z_bbZ_2 = readTable(ex69.str(),121,2);
    ex70 << tablepath << "CMS-PAS-HIG-18-005_b.dat";
    CMS13_bb_phi_phi1Z_bbZ_1 = readTable(ex70.str(),79,2);
    ex71 << tablepath << "180702826_b.dat";
    CMS13_bb_phi_phi1Z_bbZ_2 = readTable(ex71.str(),121,2);

    ex72 << tablepath << "14126663.dat";
    ATLAS8_pp_Hpm_taunu = readTable(ex72.str(),83,2);
    ex73 << tablepath << "151203704.dat";
    ATLAS8_pp_Hpm_tb = readTable(ex73.str(),41,2);
    ex74 << tablepath << "150807774_a.dat";
    CMS8_pp_Hp_taunu = readTable(ex74.str(),43,2);
    ex75 << tablepath << "150807774_b.dat";
    CMS8_pp_Hp_tb = readTable(ex75.str(),43,2);
    ex76 << tablepath << "180707915.dat";
    ATLAS13_pp_Hpm_taunu = readTable(ex76.str(),192,2);
    ex77 << tablepath << "180803599.dat";
    ATLAS13_pp_Hpm_tb = readTable(ex77.str(),181,2);
    ex78 << tablepath << "CMS-PAS-HIG-16-031.dat";
    CMS13_pp_Hpm_taunu = readTable(ex78.str(),283,2);

    bsg1 << tablepath << "bsgammatable.dat";
    arraybsgamma = readTable(bsg1.str(),1111,3);
}



double GeneralTHDMcache::ip_Br_HPtott(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_Br_HPtott_cache, NumPar, params);
    if (i>=0) {
        return ( ip_Br_HPtott_cache[NumPar][i] );
    } else {
        double newResult = pow(10.0,interpolate(br_tt,mass));
        CacheShiftReal(ip_Br_HPtott_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_Br_HPtobb(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_Br_HPtobb_cache, NumPar, params);
    if (i>=0) {
        return ( ip_Br_HPtobb_cache[NumPar][i] );
    } else {
        double newResult = pow(10.0,interpolate(br_bb,mass));
        CacheShiftReal(ip_Br_HPtobb_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_Br_HPtotautau(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_Br_HPtotautau_cache, NumPar, params);
    if (i>=0) {
        return ( ip_Br_HPtotautau_cache[NumPar][i] );
    } else {
        double newResult = pow(10.0,interpolate(br_tautau,mass));
        CacheShiftReal(ip_Br_HPtotautau_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_Br_HPtocc(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_Br_HPtocc_cache, NumPar, params);
    if (i>=0) {
        return ( ip_Br_HPtocc_cache[NumPar][i] );
    } else {
        double newResult = pow(10.0,interpolate(br_cc,mass));
        CacheShiftReal(ip_Br_HPtocc_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_Br_HPtomumu(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_Br_HPtomumu_cache, NumPar, params);
    if (i>=0) {
        return ( ip_Br_HPtomumu_cache[NumPar][i] );
    } else {
        double newResult = pow(10.0,interpolate(br_mumu,mass));
        CacheShiftReal(ip_Br_HPtomumu_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_Br_HPtoZZ(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_Br_HPtoZZ_cache, NumPar, params);
    if (i>=0) {
        return ( ip_Br_HPtoZZ_cache[NumPar][i] );
    } else {
        double newResult = pow(10.0,interpolate(br_ZZ,mass));
        CacheShiftReal(ip_Br_HPtoZZ_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_Br_HPtoWW(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_Br_HPtoWW_cache, NumPar, params);
    if (i>=0) {
        return ( ip_Br_HPtoWW_cache[NumPar][i] );
    } else {
        double newResult = pow(10.0,interpolate(br_WW,mass));
        CacheShiftReal(ip_Br_HPtoWW_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_GammaHPtotSM(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_GammaHPtotSM_cache, NumPar, params);
    if (i>=0) {
        return ( ip_GammaHPtotSM_cache[NumPar][i] );
    } else {
        double newResult = pow(10.0,interpolate(GammaHtot_SM,mass));
        CacheShiftReal(ip_GammaHPtotSM_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_cs_ggtoH_8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_ggtoH_8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_ggtoH_8_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_ggH_8,mass));
        }
        CacheShiftReal(ip_cs_ggtoH_8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_cs_ggtoH_13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_ggtoH_13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_ggtoH_13_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_ggH_13,mass));
        }
        CacheShiftReal(ip_cs_ggtoH_13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_cs_VBFtoH_8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_VBFtoH_8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_VBFtoH_8_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_VBF_8,mass));
        }
        CacheShiftReal(ip_cs_VBFtoH_8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_cs_VBFtoH_13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_VBFtoH_13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_VBFtoH_13_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_VBF_13,mass));
        }
        CacheShiftReal(ip_cs_VBFtoH_13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_cs_WtoWH_8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_WtoWH_8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_WtoWH_8_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_WH_8,mass));
        }
        CacheShiftReal(ip_cs_WtoWH_8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_cs_WtoWH_13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_WtoWH_13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_WtoWH_13_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_WH_13,mass));
        }
        CacheShiftReal(ip_cs_WtoWH_13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_cs_ZtoZH_8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_ZtoZH_8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_ZtoZH_8_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_ZH_8,mass));
        }
        CacheShiftReal(ip_cs_ZtoZH_8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_cs_ZtoZH_13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_ZtoZH_13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_ZtoZH_13_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_ZH_13,mass));
        }
        CacheShiftReal(ip_cs_ZtoZH_13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_cs_pptottH_8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_pptottH_8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_pptottH_8_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_ttH_8,mass));
        }
        CacheShiftReal(ip_cs_pptottH_8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_cs_pptottH_13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_pptottH_13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_pptottH_13_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_ttH_13,mass));
        }
        CacheShiftReal(ip_cs_pptottH_13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_cs_pptobbH_8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_pptobbH_8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_pptobbH_8_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_bbH_8,mass));
        }
        CacheShiftReal(ip_cs_pptobbH_8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_cs_pptobbH_13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_pptobbH_13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_pptobbH_13_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_bbH_13,mass));
        }
        CacheShiftReal(ip_cs_pptobbH_13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_cs_ggtoA_8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_ggtoA_8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_ggtoA_8_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_ggA_8,mass));
        }
        CacheShiftReal(ip_cs_ggtoA_8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_cs_ggtoA_13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_ggtoA_13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_ggtoA_13_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_ggA_13,mass));
        }
        CacheShiftReal(ip_cs_ggtoA_13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_cs_pptottA_8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_pptottA_8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_pptottA_8_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_ttA_8,mass));
        }
        CacheShiftReal(ip_cs_pptottA_8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_cs_pptottA_13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_pptottA_13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_pptottA_13_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_ttA_13,mass));
        }
        CacheShiftReal(ip_cs_pptottA_13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_cs_pptobbA_8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_pptobbA_8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_pptobbA_8_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_bbA_8,mass));
        }
        CacheShiftReal(ip_cs_pptobbA_8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_cs_pptobbA_13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_pptobbA_13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_pptobbA_13_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_bbA_13,mass));
        }
        CacheShiftReal(ip_cs_pptobbA_13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_cs_ggtoHp_8(double mHp, double logtb){
    int NumPar = 2;
    double params[] = {mHp, logtb};

    int i = CacheCheckReal(ip_cs_ggtoHp_8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_ggtoHp_8_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mHp>=180. && mHp <=1400. && logtb>=-1. && logtb<=1.75) {
            newResult = pow(10.0,interpolate2D(log_cs_ggHp_8, logtb, mHp));
        }
        CacheShiftReal(ip_cs_ggtoHp_8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_cs_ggtoHp_13(double mHp, double logtb){
    int NumPar = 2;
    double params[] = {mHp, logtb};

    int i = CacheCheckReal(ip_cs_ggtoHp_13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_ggtoHp_13_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mHp>=180. && mHp <=2000. && logtb>=-1. && logtb<=1.75) {
            newResult = pow(10.0,interpolate2D(log_cs_ggHp_13, logtb, mHp));
        }
        CacheShiftReal(ip_cs_ggtoHp_13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_csr_ggH_t_8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_csr_ggH_t_8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_csr_ggH_t_8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (csrH_top_8,mass);
        CacheShiftReal(ip_csr_ggH_t_8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_csr_ggH_t_13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_csr_ggH_t_13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_csr_ggH_t_13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (csrH_top_13,mass);
        CacheShiftReal(ip_csr_ggH_t_13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_csr_ggH_b_8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_csr_ggH_b_8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_csr_ggH_b_8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (csrH_bottom_8,mass);
        CacheShiftReal(ip_csr_ggH_b_8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_csr_ggH_b_13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_csr_ggH_b_13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_csr_ggH_b_13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (csrH_bottom_13,mass);
        CacheShiftReal(ip_csr_ggH_b_13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_csr_ggA_t_8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_csr_ggA_t_8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_csr_ggA_t_8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (csrA_top_8,mass);
        CacheShiftReal(ip_csr_ggA_t_8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_csr_ggA_t_13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_csr_ggA_t_13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_csr_ggA_t_13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (csrA_top_13,mass);
        CacheShiftReal(ip_csr_ggA_t_13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_csr_ggA_b_8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_csr_ggA_b_8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_csr_ggA_b_8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (csrA_bottom_8,mass);
        CacheShiftReal(ip_csr_ggA_b_8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_csr_ggA_b_13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_csr_ggA_b_13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_csr_ggA_b_13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (csrA_bottom_13,mass);
        CacheShiftReal(ip_csr_ggA_b_13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_tt_phi_tt_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_tt_phi_tt_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_tt_phi_tt_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS13_tt_phi_tt,mass);
        CacheShiftReal(ip_ex_tt_phi_tt_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_bb_phi_tt_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_tt_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phi_tt_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS13_bb_phi_tt,mass);
        CacheShiftReal(ip_ex_bb_phi_tt_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_bb_phi_bb_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_bb_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phi_bb_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS8_bb_phi_bb,mass);
        CacheShiftReal(ip_ex_bb_phi_bb_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_gg_phi_bb_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_bb_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_bb_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS8_gg_phi_bb,mass);
        CacheShiftReal(ip_ex_gg_phi_bb_CMS8_cache, NumPar, params, newResult);
        return newResult;      
    }
}

double GeneralTHDMcache::ip_ex_pp_phi_bb_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_bb_CMS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_phi_bb_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS13_pp_phi_bb,mass);
        CacheShiftReal(ip_ex_pp_phi_bb_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_bb_phi_bb_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_bb_CMS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phi_bb_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS13_bb_phi_bb,mass);
        CacheShiftReal(ip_ex_bb_phi_bb_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_gg_phi_tautau_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_tautau_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_tautau_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS8_gg_phi_tautau,mass);
        CacheShiftReal(ip_ex_gg_phi_tautau_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_gg_phi_tautau_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_tautau_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_tautau_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS8_gg_phi_tautau,mass);
        CacheShiftReal(ip_ex_gg_phi_tautau_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_bb_phi_tautau_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_tautau_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phi_tautau_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS8_bb_phi_tautau,mass);
        CacheShiftReal(ip_ex_bb_phi_tautau_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_bb_phi_tautau_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_tautau_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phi_tautau_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS8_bb_phi_tautau,mass);
        CacheShiftReal(ip_ex_bb_phi_tautau_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_gg_phi_tautau_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_tautau_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_tautau_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS13_gg_phi_tautau,mass);
        CacheShiftReal(ip_ex_gg_phi_tautau_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_gg_phi_tautau_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_tautau_CMS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_tautau_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS13_gg_phi_tautau,mass);
        CacheShiftReal(ip_ex_gg_phi_tautau_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_bb_phi_tautau_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_tautau_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phi_tautau_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS13_bb_phi_tautau,mass);
        CacheShiftReal(ip_ex_bb_phi_tautau_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_bb_phi_tautau_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_tautau_CMS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phi_tautau_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS13_bb_phi_tautau,mass);
        CacheShiftReal(ip_ex_bb_phi_tautau_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_gg_phi_gaga_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_gaga_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_gaga_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS8_gg_phi_gaga,mass);
        CacheShiftReal(ip_ex_gg_phi_gaga_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_pp_phi_gaga_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_gaga_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_phi_gaga_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS13_pp_phi_gaga,mass);
        CacheShiftReal(ip_ex_pp_phi_gaga_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_gg_phi_gaga_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_gaga_CMS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_gaga_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS13_gg_phi_gaga,mass);
        CacheShiftReal(ip_ex_gg_phi_gaga_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_pp_phi_Zga_llga_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_Zga_llga_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_phi_Zga_llga_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS8_pp_phi_Zga_llga,mass);
        CacheShiftReal(ip_ex_pp_phi_Zga_llga_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_pp_phi_Zga_llga_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_Zga_llga_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_phi_Zga_llga_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS8_pp_phi_Zga_llga,mass);
        CacheShiftReal(ip_ex_pp_phi_Zga_llga_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_gg_phi_Zga_llga_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_Zga_llga_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_Zga_llga_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS13_gg_phi_Zga_llga,mass);
        CacheShiftReal(ip_ex_gg_phi_Zga_llga_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_gg_phi_Zga_qqga_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_Zga_qqga_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_Zga_qqga_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS13_gg_phi_Zga_qqga,mass);
        CacheShiftReal(ip_ex_gg_phi_Zga_qqga_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_gg_phi_Zga_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_Zga_CMS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_Zga_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS13_gg_phi_Zga,mass);
        CacheShiftReal(ip_ex_gg_phi_Zga_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_gg_phi_ZZ_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_ZZ_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_ZZ_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS8_gg_phi_ZZ,mass);
        CacheShiftReal(ip_ex_gg_phi_ZZ_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_VV_phi_ZZ_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VV_phi_ZZ_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_VV_phi_ZZ_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS8_VV_phi_ZZ,mass);
        CacheShiftReal(ip_ex_VV_phi_ZZ_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_gg_phi_ZZ_llllnunu_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_ZZ_llllnunu_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_ZZ_llllnunu_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS13_gg_phi_ZZ_llllnunu,mass);
        CacheShiftReal(ip_ex_gg_phi_ZZ_llllnunu_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_VV_phi_ZZ_llllnunu_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VV_phi_ZZ_llllnunu_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_VV_phi_ZZ_llllnunu_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS13_VV_phi_ZZ_llllnunu,mass);
        CacheShiftReal(ip_ex_VV_phi_ZZ_llllnunu_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_gg_phi_ZZ_qqllnunu_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_ZZ_qqllnunu_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_ZZ_qqllnunu_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS13_gg_phi_ZZ_qqllnunu,mass);
        CacheShiftReal(ip_ex_gg_phi_ZZ_qqllnunu_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_VV_phi_ZZ_qqllnunu_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VV_phi_ZZ_qqllnunu_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_VV_phi_ZZ_qqllnunu_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS13_VV_phi_ZZ_qqllnunu,mass);
        CacheShiftReal(ip_ex_VV_phi_ZZ_qqllnunu_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_pp_phi_ZZ_llqqnunull_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_ZZ_llqqnunull_CMS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_phi_ZZ_llqqnunull_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS13_pp_phi_ZZ_llqqnunull,mass);
        CacheShiftReal(ip_ex_pp_phi_ZZ_llqqnunull_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_pp_phi_ZZ_qqnunu_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_ZZ_qqnunu_CMS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_phi_ZZ_qqnunu_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS13_pp_phi_ZZ_qqnunu,mass);
        CacheShiftReal(ip_ex_pp_phi_ZZ_qqnunu_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_gg_phi_WW_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_WW_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_WW_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS8_gg_phi_WW,mass);
        CacheShiftReal(ip_ex_gg_phi_WW_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_VV_phi_WW_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VV_phi_WW_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_VV_phi_WW_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS8_VV_phi_WW,mass);
        CacheShiftReal(ip_ex_VV_phi_WW_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_gg_phi_WW_enumunu_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_WW_enumunu_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_WW_enumunu_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS13_gg_phi_WW_enumunu,mass);
        CacheShiftReal(ip_ex_gg_phi_WW_enumunu_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_VV_phi_WW_enumunu_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VV_phi_WW_enumunu_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_VV_phi_WW_enumunu_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS13_VV_phi_WW_enumunu,mass);
        CacheShiftReal(ip_ex_VV_phi_WW_enumunu_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_ggVV_phi_WW_lnulnu_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_ggVV_phi_WW_lnulnu_CMS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_ggVV_phi_WW_lnulnu_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS13_ggVV_phi_WW_lnulnu,mass);
        CacheShiftReal(ip_ex_ggVV_phi_WW_lnulnu_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_gg_phi_WW_lnuqq_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_WW_lnuqq_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_WW_lnuqq_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS13_gg_phi_WW_lnuqq,mass);
        CacheShiftReal(ip_ex_gg_phi_WW_lnuqq_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_VV_phi_WW_lnuqq_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VV_phi_WW_lnuqq_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_VV_phi_WW_lnuqq_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS13_VV_phi_WW_lnuqq,mass);
        CacheShiftReal(ip_ex_VV_phi_WW_lnuqq_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_pp_phi_WW_lnuqq_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_WW_lnuqq_CMS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_phi_WW_lnuqq_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS13_pp_phi_WW_lnuqq,mass);
        CacheShiftReal(ip_ex_pp_phi_WW_lnuqq_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_mu_pp_phi_VV_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_mu_pp_phi_VV_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_mu_pp_phi_VV_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS8_mu_pp_phi_VV,mass);
        CacheShiftReal(ip_ex_mu_pp_phi_VV_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_pp_phi_VV_qqqq_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_VV_qqqq_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_phi_VV_qqqq_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS13_pp_phi_VV_qqqq,mass);
        CacheShiftReal(ip_ex_pp_phi_VV_qqqq_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_gg_phi_phi1phi1_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_phi1phi1_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_phi_phi1phi1_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS8_gg_phi_phi1phi1,mass);
        CacheShiftReal(ip_ex_gg_phi_phi1phi1_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_pp_phi_phi1phi1_bbbb_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_phi1phi1_bbbb_CMS8_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_phi1phi1_bbbb_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_pp_phi_phi1phi1_bbbb,mass);
        CacheShiftReal(ip_ex_pp_phi_phi1phi1_bbbb_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_pp_phi_phi1phi1_bbgaga_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_phi1phi1_bbgaga_CMS8_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_phi1phi1_bbgaga_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_pp_phi_phi1phi1_bbgaga,mass);
        CacheShiftReal(ip_ex_pp_phi_phi1phi1_bbgaga_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_gg_phi_phi1phi1_bbtautau_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_phi1phi1_bbtautau_CMS8_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_phi_phi1phi1_bbtautau_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_gg_phi_phi1phi1_bbtautau,mass);
        CacheShiftReal(ip_ex_gg_phi_phi1phi1_bbtautau_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_pp_phi_phi1phi1_bbtautau_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_phi1phi1_bbtautau_CMS8_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_phi1phi1_bbtautau_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_pp_phi_phi1phi1_bbtautau,mass);
        CacheShiftReal(ip_ex_pp_phi_phi1phi1_bbtautau_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_pp_phi_phi1phi1_bbbb_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_phi1phi1_bbbb_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_phi1phi1_bbbb_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_phi_phi1phi1_bbbb,mass);
        CacheShiftReal(ip_ex_pp_phi_phi1phi1_bbbb_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_pp_phi_phi1phi1_bbbb_1_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_phi1phi1_bbbb_1_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_phi1phi1_bbbb_1_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_phi_phi1phi1_bbbb_1,mass);
        CacheShiftReal(ip_ex_pp_phi_phi1phi1_bbbb_1_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_pp_phi_phi1phi1_bbbb_2_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_phi1phi1_bbbb_2_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_phi1phi1_bbbb_2_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_phi_phi1phi1_bbbb_2,mass);
        CacheShiftReal(ip_ex_pp_phi_phi1phi1_bbbb_2_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_pp_phi_phi1phi1_bbgaga_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_phi1phi1_bbgaga_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_phi1phi1_bbgaga_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_phi_phi1phi1_bbgaga,mass);
        CacheShiftReal(ip_ex_pp_phi_phi1phi1_bbgaga_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_pp_phi_phi1phi1_bbgaga_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_phi1phi1_bbgaga_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_phi1phi1_bbgaga_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_phi_phi1phi1_bbgaga,mass);
        CacheShiftReal(ip_ex_pp_phi_phi1phi1_bbgaga_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_pp_phi_phi1phi1_bbtautau_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_phi1phi1_bbtautau_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_phi1phi1_bbtautau_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_phi_phi1phi1_bbtautau,mass);
        CacheShiftReal(ip_ex_pp_phi_phi1phi1_bbtautau_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_pp_phi_phi1phi1_bbtautau_1_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_phi1phi1_bbtautau_1_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_phi1phi1_bbtautau_1_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_phi_phi1phi1_bbtautau_1,mass);
        CacheShiftReal(ip_ex_pp_phi_phi1phi1_bbtautau_1_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_pp_phi_phi1phi1_bbtautau_2_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_phi1phi1_bbtautau_2_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_phi1phi1_bbtautau_2_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_phi_phi1phi1_bbtautau_2,mass);
        CacheShiftReal(ip_ex_pp_phi_phi1phi1_bbtautau_2_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_pp_phi_phi1phi1_bbVV_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_phi1phi1_bbVV_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_phi1phi1_bbVV_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_phi_phi1phi1_bbVV,mass);
        CacheShiftReal(ip_ex_pp_phi_phi1phi1_bbVV_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}


double GeneralTHDMcache::ip_ex_pp_phi_phi1phi1_bbWW_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_phi1phi1_bbWW_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_phi1phi1_bbWW_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_phi_phi1phi1_bbWW,mass);
        CacheShiftReal(ip_ex_pp_phi_phi1phi1_bbWW_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_gg_phi_phi1phi1_gagaWW_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_phi1phi1_gagaWW_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_phi_phi1phi1_gagaWW_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_gg_phi_phi1phi1_gagaWW,mass);
        CacheShiftReal(ip_ex_gg_phi_phi1phi1_gagaWW_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_gg_phi_phi1Z_bbZ_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_phi1Z_bbZ_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_phi1Z_bbZ_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS8_gg_phi_phi1Z_bbZ,mass);
        CacheShiftReal(ip_ex_gg_phi_phi1Z_bbZ_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_gg_phi_phi1Z_bbll_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_phi1Z_bbll_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_phi1Z_bbll_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS8_gg_phi_phi1Z_bbll,mass);
        CacheShiftReal(ip_ex_gg_phi_phi1Z_bbll_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_gg_phi_phi1Z_tautauZ_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_phi1Z_tautauZ_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_phi1Z_tautauZ_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS8_gg_phi_phi1Z_tautauZ,mass);
        CacheShiftReal(ip_ex_gg_phi_phi1Z_tautauZ_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_gg_phi_phi1Z_tautaull_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_phi1Z_tautaull_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_phi1Z_tautaull_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS8_gg_phi_phi1Z_tautaull,mass);
        CacheShiftReal(ip_ex_gg_phi_phi1Z_tautaull_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_gg_phi_phi1Z_bbZ_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_phi1Z_bbZ_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_phi1Z_bbZ_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS13_gg_phi_phi1Z_bbZ,mass);
        CacheShiftReal(ip_ex_gg_phi_phi1Z_bbZ_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_gg_phi_phi1Z_bbZ_1_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_phi1Z_bbZ_1_CMS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_phi1Z_bbZ_1_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS13_gg_phi_phi1Z_bbZ_1,mass);
        CacheShiftReal(ip_ex_gg_phi_phi1Z_bbZ_1_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_gg_phi_phi1Z_bbZ_2_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_phi1Z_bbZ_2_CMS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_phi1Z_bbZ_2_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS13_gg_phi_phi1Z_bbZ_2,mass);
        CacheShiftReal(ip_ex_gg_phi_phi1Z_bbZ_2_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_bb_phi_phi1Z_bbZ_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_phi1Z_bbZ_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phi_phi1Z_bbZ_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS13_bb_phi_phi1Z_bbZ,mass);
        CacheShiftReal(ip_ex_bb_phi_phi1Z_bbZ_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_bb_phi_phi1Z_bbZ_1_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_phi1Z_bbZ_1_CMS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phi_phi1Z_bbZ_1_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS13_bb_phi_phi1Z_bbZ_1,mass);
        CacheShiftReal(ip_ex_bb_phi_phi1Z_bbZ_1_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_bb_phi_phi1Z_bbZ_2_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_phi1Z_bbZ_2_CMS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phi_phi1Z_bbZ_2_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS13_bb_phi_phi1Z_bbZ_2,mass);
        CacheShiftReal(ip_ex_bb_phi_phi1Z_bbZ_2_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_pp_phii_phijZ_bbll_1_CMS8(double mi, double mj){
    int NumPar = 2;
    double params[] = {mi, mj};

    int i = CacheCheckReal(ip_ex_pp_phii_phijZ_bbll_1_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_phii_phijZ_bbll_1_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate2D(CMS8_pp_phii_phijZ_bbll_1, mi, mj);
        CacheShiftReal(ip_ex_pp_phii_phijZ_bbll_1_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_ex_pp_phii_phijZ_bbll_2_CMS8(double mi, double mj){
    int NumPar = 2;
    double params[] = {mi, mj};

    int i = CacheCheckReal(ip_ex_pp_phii_phijZ_bbll_2_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_phii_phijZ_bbll_2_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate2D(CMS8_pp_phii_phijZ_bbll_2,mi, mj);
        CacheShiftReal(ip_ex_pp_phii_phijZ_bbll_2_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}


double GeneralTHDMcache::ip_ex_pp_phii_phijZ_tautaull_1_CMS8(double mi, double mj){
    int NumPar = 2;
    double params[] = {mi, mj};

    int i = CacheCheckReal(ip_ex_pp_phii_phijZ_tautaull_1_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_phii_phijZ_tautaull_1_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate2D(CMS8_pp_phii_phijZ_tautaull_1,mi, mj);
        CacheShiftReal(ip_ex_pp_phii_phijZ_tautaull_1_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_ex_pp_phii_phijZ_tautaull_2_CMS8(double mi, double mj){
    int NumPar = 2;
    double params[] = {mi, mj};

    int i = CacheCheckReal(ip_ex_pp_phii_phijZ_tautaull_2_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_phii_phijZ_tautaull_2_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate2D(CMS8_pp_phii_phijZ_tautaull_2, mi, mj);
        CacheShiftReal(ip_ex_pp_phii_phijZ_tautaull_2_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_ex_gg_phii_phijZ_bbZ_ATLAS13(double mj, double mi){
    int NumPar = 2;
    double params[] = {mj, mi};

    int i = CacheCheckReal(ip_ex_gg_phii_phijZ_bbZ_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phii_phijZ_bbZ_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate2D(ATLAS13_gg_phii_phijZ_bbZ, mj, mi);
        CacheShiftReal(ip_ex_gg_phii_phijZ_bbZ_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}


double GeneralTHDMcache::ip_ex_bb_phii_phijZ_bbZ_ATLAS13(double mj, double mi){
    int NumPar = 2;
    double params[] = {mj, mi};

    int i = CacheCheckReal(ip_ex_bb_phii_phijZ_bbZ_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phii_phijZ_bbZ_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate2D(ATLAS13_bb_phii_phijZ_bbZ, mj, mi);
        CacheShiftReal(ip_ex_bb_phii_phijZ_bbZ_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}


//double GeneralTHDMcache::ip_ex_pp_phi_gaga_ATLAS8(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_pp_phi_gaga_ATLAS8_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_pp_phi_gaga_ATLAS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate(ATLAS8_pp_phi_gaga,mass);
//        CacheShiftReal(ip_ex_pp_phi_gaga_ATLAS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_pp_phi_gaga_ATLAS8_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_pp_phi_gaga_ATLAS8_cache_e, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_pp_phi_gaga_ATLAS8_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate(ATLAS8_pp_phi_gaga_e,mass);
//        CacheShiftReal(ip_ex_pp_phi_gaga_ATLAS8_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//double GeneralTHDMcache::ip_ex_pp_phi_Zga_llga_ATLAS8(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_pp_phi_Zga_llga_ATLAS8_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_pp_phi_Zga_llga_ATLAS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS8_pp_phi_Zga_llga,mass);
//        CacheShiftReal(ip_ex_pp_phi_Zga_llga_ATLAS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_pp_phi_Zga_llga_ATLAS8_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_pp_phi_Zga_llga_ATLAS8_cache_e, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_pp_phi_Zga_llga_ATLAS8_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS8_pp_phi_Zga_llga_e,mass);
//        CacheShiftReal(ip_ex_pp_phi_Zga_llga_ATLAS8_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_phi_tautau_ATLAS8(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_phi_tautau_ATLAS8_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_phi_tautau_ATLAS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate(ATLAS8_gg_phi_tautau,mass);
//        CacheShiftReal(ip_ex_gg_phi_tautau_ATLAS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_phi_tautau_ATLAS8_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_phi_tautau_ATLAS8_cache_e, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_phi_tautau_ATLAS8_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate(ATLAS8_gg_phi_tautau_e,mass);
//        CacheShiftReal(ip_ex_gg_phi_tautau_ATLAS8_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_bb_phi_tautau_ATLAS8(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_bb_phi_tautau_ATLAS8_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_bb_phi_tautau_ATLAS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate(ATLAS8_bb_phi_tautau,mass);
//        CacheShiftReal(ip_ex_bb_phi_tautau_ATLAS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_bb_phi_tautau_ATLAS8_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_bb_phi_tautau_ATLAS8_cache_e, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_bb_phi_tautau_ATLAS8_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate(ATLAS8_bb_phi_tautau_e,mass);
//        CacheShiftReal(ip_ex_bb_phi_tautau_ATLAS8_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_A_hZ_tautauZ_ATLAS8(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_A_hZ_tautauZ_ATLAS8_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_A_hZ_tautauZ_ATLAS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate(ATLAS8_gg_A_hZ_tautauZ,mass);
//        CacheShiftReal(ip_ex_gg_A_hZ_tautauZ_ATLAS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_A_hZ_tautauZ_ATLAS8_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_A_hZ_tautauZ_ATLAS8_cache_e, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_A_hZ_tautauZ_ATLAS8_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate(ATLAS8_gg_A_hZ_tautauZ_e,mass);
//        CacheShiftReal(ip_ex_gg_A_hZ_tautauZ_ATLAS8_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_A_hZ_bbZ_ATLAS8(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_A_hZ_bbZ_ATLAS8_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_A_hZ_bbZ_ATLAS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate(ATLAS8_gg_A_hZ_bbZ,mass);
//        CacheShiftReal(ip_ex_gg_A_hZ_bbZ_ATLAS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//
//
//double GeneralTHDMcache::ip_ex_gg_A_hZ_bbZ_ATLAS8_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_A_hZ_bbZ_ATLAS8_cache_e, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_A_hZ_bbZ_ATLAS8_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate(ATLAS8_gg_A_hZ_bbZ_e,mass);
//        CacheShiftReal(ip_ex_gg_A_hZ_bbZ_ATLAS8_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_phi_tt_ATLAS8(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_phi_tt_ATLAS8_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_phi_tt_ATLAS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS8_gg_phi_tt,mass);
//        CacheShiftReal(ip_ex_gg_phi_tt_ATLAS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_phi_tt_ATLAS8_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_phi_tt_ATLAS8_cache_e, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_phi_tt_ATLAS8_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS8_gg_phi_tt_e,mass);
//        CacheShiftReal(ip_ex_gg_phi_tt_ATLAS8_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_H_WW_ATLAS8(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_H_WW_ATLAS8_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_H_WW_ATLAS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS8_gg_H_WW,mass);
//        CacheShiftReal(ip_ex_gg_H_WW_ATLAS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_H_WW_ATLAS8_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_H_WW_ATLAS8_cache_e, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_H_WW_ATLAS8_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS8_gg_H_WW_e,mass);
//        CacheShiftReal(ip_ex_gg_H_WW_ATLAS8_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_VBF_H_WW_ATLAS8(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_VBF_H_WW_ATLAS8_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_VBF_H_WW_ATLAS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS8_VBF_H_WW,mass);
//        CacheShiftReal(ip_ex_VBF_H_WW_ATLAS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_VBF_H_WW_ATLAS8_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_VBF_H_WW_ATLAS8_cache_e, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_VBF_H_WW_ATLAS8_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS8_VBF_H_WW_e,mass);
//        CacheShiftReal(ip_ex_VBF_H_WW_ATLAS8_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_H_ZZ_ATLAS8(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_H_ZZ_ATLAS8_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_H_ZZ_ATLAS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS8_gg_H_ZZ,mass);
//        CacheShiftReal(ip_ex_gg_H_ZZ_ATLAS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_H_ZZ_ATLAS8_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_H_ZZ_ATLAS8_cache_e, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_H_ZZ_ATLAS8_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS8_gg_H_ZZ_e,mass);
//        CacheShiftReal(ip_ex_gg_H_ZZ_ATLAS8_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_VBF_H_ZZ_ATLAS8(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_VBF_H_ZZ_ATLAS8_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_VBF_H_ZZ_ATLAS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS8_VBF_H_ZZ,mass);
//        CacheShiftReal(ip_ex_VBF_H_ZZ_ATLAS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_VBF_H_ZZ_ATLAS8_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_VBF_H_ZZ_ATLAS8_cache_e, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_VBF_H_ZZ_ATLAS8_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS8_VBF_H_ZZ_e,mass);
//        CacheShiftReal(ip_ex_VBF_H_ZZ_ATLAS8_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_H_hh_ATLAS8(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_H_hh_ATLAS8_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_H_hh_ATLAS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS8_gg_H_hh,mass);
//        CacheShiftReal(ip_ex_gg_H_hh_ATLAS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_H_hh_ATLAS8_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_H_hh_ATLAS8_cache_e, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_H_hh_ATLAS8_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS8_gg_H_hh_e,mass);
//        CacheShiftReal(ip_ex_gg_H_hh_ATLAS8_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_mu_pp_H_VV_CMS8(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_mu_pp_H_VV_CMS8_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_mu_pp_H_VV_CMS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate(CMS8_mu_pp_H_VV,mass);
//        CacheShiftReal(ip_ex_mu_pp_H_VV_CMS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_mu_pp_H_VV_CMS8_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_mu_pp_H_VV_CMS8_cache_e, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_mu_pp_H_VV_CMS8_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate(CMS8_mu_pp_H_VV_e,mass);
//        CacheShiftReal(ip_ex_mu_pp_H_VV_CMS8_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_A_hZ_bbll_CMS8(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_A_hZ_bbll_CMS8_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_A_hZ_bbll_CMS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate(CMS8_gg_A_hZ_bbll,mass);
//        CacheShiftReal(ip_ex_gg_A_hZ_bbll_CMS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_A_hZ_bbll_CMS8_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_A_hZ_bbll_CMS8_cache_e, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_A_hZ_bbll_CMS8_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate(CMS8_gg_A_hZ_bbll_e,mass);
//        CacheShiftReal(ip_ex_gg_A_hZ_bbll_CMS8_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_pp_H_hh_CMS8(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_pp_H_hh_CMS8_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_pp_H_hh_CMS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate(CMS8_pp_H_hh,mass);
//        CacheShiftReal(ip_ex_pp_H_hh_CMS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_pp_H_hh_CMS8_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_pp_H_hh_CMS8_cache_e, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_pp_H_hh_CMS8_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate(CMS8_pp_H_hh_e,mass);
//        CacheShiftReal(ip_ex_pp_H_hh_CMS8_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_pp_phi_hh_gagabb_CMS8(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_pp_phi_hh_gagabb_CMS8_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_pp_phi_hh_gagabb_CMS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate(CMS8_pp_H_hh_gagabb,mass);
//        CacheShiftReal(ip_ex_pp_phi_hh_gagabb_CMS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_pp_phi_hh_gagabb_CMS8_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_pp_phi_hh_gagabb_CMS8_cache_e, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_pp_phi_hh_gagabb_CMS8_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate(CMS8_pp_H_hh_gagabb_e,mass);
//        CacheShiftReal(ip_ex_pp_phi_hh_gagabb_CMS8_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_pp_phi_hh_bbbb_CMS8(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_pp_phi_hh_bbbb_CMS8_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_pp_phi_hh_bbbb_CMS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate(CMS8_pp_H_hh_bbbb,mass);
//        CacheShiftReal(ip_ex_pp_phi_hh_bbbb_CMS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_pp_phi_hh_bbbb_CMS8_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_pp_phi_hh_bbbb_CMS8_cache_e, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_pp_phi_hh_bbbb_CMS8_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate(CMS8_pp_H_hh_bbbb_e,mass);
//        CacheShiftReal(ip_ex_pp_phi_hh_bbbb_CMS8_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_bb_phi_bb_CMS8(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_bb_phi_bb_CMS8_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_bb_phi_bb_CMS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (CMS8_bb_phi_bb,mass);
//        CacheShiftReal(ip_ex_bb_phi_bb_CMS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_bb_phi_bb_CMS8_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_bb_phi_bb_CMS8_cache_e, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_bb_phi_bb_CMS8_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (CMS8_bb_phi_bb_e,mass);
//        CacheShiftReal(ip_ex_bb_phi_bb_CMS8_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_phi_tautau_CMS8(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_phi_tautau_CMS8_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_phi_tautau_CMS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (CMS8_gg_phi_tautau,mass);
//        CacheShiftReal(ip_ex_gg_phi_tautau_CMS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_phi_tautau_CMS8_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_phi_tautau_CMS8_cache_e, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_phi_tautau_CMS8_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (CMS8_gg_phi_tautau_e,mass);
//        CacheShiftReal(ip_ex_gg_phi_tautau_CMS8_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_bb_phi_tautau_CMS8(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_bb_phi_tautau_CMS8_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_bb_phi_tautau_CMS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (CMS8_bb_phi_tautau,mass);
//        CacheShiftReal(ip_ex_bb_phi_tautau_CMS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_bb_phi_tautau_CMS8_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_bb_phi_tautau_CMS8_cache_e, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_bb_phi_tautau_CMS8_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (CMS8_bb_phi_tautau_e,mass);
//        CacheShiftReal(ip_ex_bb_phi_tautau_CMS8_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_phi_gaga_CMS8(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_phi_gaga_CMS8_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_phi_gaga_CMS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (CMS8_gg_phi_gaga,mass);
//        CacheShiftReal(ip_ex_gg_phi_gaga_CMS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_phi_gaga_CMS8_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_phi_gaga_CMS8_cache_e, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_phi_gaga_CMS8_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (CMS8_gg_phi_gaga_e,mass);
//        CacheShiftReal(ip_ex_gg_phi_gaga_CMS8_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
////
////double GeneralTHDMcache::ip_ex_gg_phi_gaga_CMS_ep1(double mass){
////    int NumPar = 1;
////    double params[] = {mass};
////
////    int i = CacheCheckReal(ip_ex_gg_phi_gaga_CMS_cache_ep1, NumPar, params);
////    if (i>=0) {
////        return ( ip_ex_gg_phi_gaga_CMS_cache_ep1[NumPar][i] );
////    } else {
////        double newResult = interpolate (CMS_ggF_phi_gaga_ep1,mass);
////        CacheShiftReal(ip_ex_gg_phi_gaga_CMS_cache_ep1, NumPar, params, newResult);
////        return newResult;
////    }
////}
//
////double GeneralTHDMcache::ip_ex_gg_phi_gaga_CMS_ep2(double mass){
////    int NumPar = 1;
////    double params[] = {mass};
////
////    int i = CacheCheckReal(ip_ex_gg_phi_gaga_CMS_cache_ep2, NumPar, params);
////    if (i>=0) {
////        return ( ip_ex_gg_phi_gaga_CMS_cache_ep2[NumPar][i] );
////    } else {
////        double newResult = interpolate (CMS_ggF_phi_gaga_ep2,mass);
////        CacheShiftReal(ip_ex_gg_phi_gaga_CMS_cache_ep2, NumPar, params, newResult);
////        return newResult;
////    }
////}
////
////double GeneralTHDMcache::ip_ex_gg_phi_gaga_CMS_em1(double mass){
////    int NumPar = 1;
////    double params[] = {mass};
////
////    int i = CacheCheckReal(ip_ex_gg_phi_gaga_CMS_cache_em1, NumPar, params);
////    if (i>=0) {
////        return ( ip_ex_gg_phi_gaga_CMS_cache_em1[NumPar][i] );
////    } else {
////        double newResult = interpolate (CMS_ggF_phi_gaga_em1,mass);
////        CacheShiftReal(ip_ex_gg_phi_gaga_CMS_cache_em1, NumPar, params, newResult);
////        return newResult;
////    }
////}
////
////double GeneralTHDMcache::ip_ex_gg_phi_gaga_CMS_em2(double mass){
////    int NumPar = 1;
////    double params[] = {mass};
////
////    int i = CacheCheckReal(ip_ex_gg_phi_gaga_CMS_cache_em2, NumPar, params);
////    if (i>=0) {
////        return ( ip_ex_gg_phi_gaga_CMS_cache_em2[NumPar][i] );
////    } else {
////        double newResult = interpolate (CMS_ggF_phi_gaga_em2,mass);
////        CacheShiftReal(ip_ex_gg_phi_gaga_CMS_cache_em2, NumPar, params, newResult);
////        return newResult;
////    }
////}
////
//
//
//double GeneralTHDMcache::ip_ex_pp_A_Zga_llga_CMS8(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_pp_A_Zga_llga_CMS8_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_pp_A_Zga_llga_CMS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (CMS8_pp_A_Zga_llga,mass);
//        CacheShiftReal(ip_ex_pp_A_Zga_llga_CMS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_pp_A_Zga_llga_CMS8_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_pp_A_Zga_llga_CMS8_cache_e, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_pp_A_Zga_llga_CMS8_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (CMS8_pp_A_Zga_llga_e,mass);
//        CacheShiftReal(ip_ex_pp_A_Zga_llga_CMS8_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_H_hh_bbtautau_CMS8(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_H_hh_bbtautau_CMS8_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_H_hh_bbtautau_CMS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (CMS8_gg_H_hh_bbtautau,mass);
//        CacheShiftReal(ip_ex_gg_H_hh_bbtautau_CMS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_H_hh_bbtautau_CMS8_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_H_hh_bbtautau_CMS8_cache_e, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_H_hh_bbtautau_CMS8_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (CMS8_gg_H_hh_bbtautau_e,mass);
//        CacheShiftReal(ip_ex_gg_H_hh_bbtautau_CMS8_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_A_hZ_tautaull_CMS8(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_A_hZ_tautaull_CMS8_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_A_hZ_tautaull_CMS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (CMS8_gg_A_hZ_tautaull,mass);
//        CacheShiftReal(ip_ex_gg_A_hZ_tautaull_CMS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_A_hZ_tautaull_CMS8_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_A_hZ_tautaull_CMS8_cache_e, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_A_hZ_tautaull_CMS8_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (CMS8_gg_A_hZ_tautaull_e,mass);
//        CacheShiftReal(ip_ex_gg_A_hZ_tautaull_CMS8_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_pp_A_HZ_bbll_CMS8(double m1, double m2){
//    int NumPar = 2;
//    double params[] = {m1, m2};
//    int i = CacheCheckReal(ip_ex_pp_A_HZ_bbll_CMS8_cache, NumPar, params);
//   
//    if (i>=0) {
//        
//        return ( ip_ex_pp_A_HZ_bbll_CMS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate2D(CMS8_pp_A_HZ_bbll, m1, m2);
//        
//        
//        CacheShiftReal(ip_ex_pp_A_HZ_bbll_CMS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_pp_H_AZ_bbll_CMS8(double mA, double mH){
//    int NumPar = 2;
//    double params[] = {mA, mH};
//
//    int i = CacheCheckReal(ip_ex_pp_H_AZ_bbll_CMS8_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_pp_H_AZ_bbll_CMS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate2D(CMS8_pp_H_AZ_bbll, mA, mH);
//        CacheShiftReal(ip_ex_pp_H_AZ_bbll_CMS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//double GeneralTHDMcache::ip_ex_pp_A_HZ_tautaull_CMS8(double mA, double mH){
//    int NumPar = 2;
//    double params[] = {mA, mH};
//
//    int i = CacheCheckReal(ip_ex_pp_A_HZ_tautaull_CMS8_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_pp_A_HZ_tautaull_CMS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate2D(CMS8_pp_A_HZ_tautaull, mA, mH);
//        CacheShiftReal(ip_ex_pp_A_HZ_tautaull_CMS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_pp_H_AZ_tautaull_CMS8(double mA, double mH){
//    int NumPar = 2;
//    double params[] = {mA, mH};
//
//    int i = CacheCheckReal(ip_ex_pp_H_AZ_tautaull_CMS8_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_pp_H_AZ_tautaull_CMS8_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate2D(CMS8_pp_H_AZ_tautaull, mA, mH);
//        CacheShiftReal(ip_ex_pp_H_AZ_tautaull_CMS8_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_bb_phi_tt_ATLAS13(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_bb_phi_tt_ATLAS13_cache, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_bb_phi_tt_ATLAS13_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_bb_phi_tt,mass);
//        CacheShiftReal(ip_ex_bb_phi_tt_ATLAS13_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_bb_phi_tt_ATLAS13_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_bb_phi_tt_ATLAS13_cache_e, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_bb_phi_tt_ATLAS13_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_bb_phi_tt_e,mass);
//        CacheShiftReal(ip_ex_bb_phi_tt_ATLAS13_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_tt_phi_tt_ATLAS13(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_tt_phi_tt_ATLAS13_cache, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_tt_phi_tt_ATLAS13_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_tt_phi_tt,mass);
//        CacheShiftReal(ip_ex_tt_phi_tt_ATLAS13_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_tt_phi_tt_ATLAS13_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_tt_phi_tt_ATLAS13_cache_e, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_tt_phi_tt_ATLAS13_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_tt_phi_tt_e,mass);
//        CacheShiftReal(ip_ex_tt_phi_tt_ATLAS13_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_phi_tautau_ATLAS13(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_phi_tautau_ATLAS13_cache, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_gg_phi_tautau_ATLAS13_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_gg_phi_tautau,mass);
//        CacheShiftReal(ip_ex_gg_phi_tautau_ATLAS13_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_phi_tautau_ATLAS13_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_phi_tautau_ATLAS13_cache_e, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_gg_phi_tautau_ATLAS13_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_gg_phi_tautau_e,mass);
//        CacheShiftReal(ip_ex_gg_phi_tautau_ATLAS13_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_bb_phi_tautau_ATLAS13(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_bb_phi_tautau_ATLAS13_cache, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_bb_phi_tautau_ATLAS13_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_bb_phi_tautau,mass);
//        CacheShiftReal(ip_ex_bb_phi_tautau_ATLAS13_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_bb_phi_tautau_ATLAS13_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_bb_phi_tautau_ATLAS13_cache_e, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_bb_phi_tautau_ATLAS13_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_bb_phi_tautau_e,mass);
//        CacheShiftReal(ip_ex_bb_phi_tautau_ATLAS13_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_pp_phi_gaga_ATLAS13(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_pp_phi_gaga_ATLAS13_cache, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_pp_phi_gaga_ATLAS13_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_pp_phi_gaga,mass);
//        CacheShiftReal(ip_ex_pp_phi_gaga_ATLAS13_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_pp_phi_gaga_ATLAS13_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_pp_phi_gaga_ATLAS13_cache_e, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_pp_phi_gaga_ATLAS13_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_pp_phi_gaga_e,mass);
//        CacheShiftReal(ip_ex_pp_phi_gaga_ATLAS13_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_pp_phi_Zga_llga_ATLAS13(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_pp_phi_Zga_llga_ATLAS13_cache, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_pp_phi_Zga_llga_ATLAS13_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_pp_phi_Zga,mass);
//        CacheShiftReal(ip_ex_pp_phi_Zga_llga_ATLAS13_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_pp_phi_Zga_llga_ATLAS13_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_pp_phi_Zga_llga_ATLAS13_cache_e, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_pp_phi_Zga_llga_ATLAS13_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_pp_phi_Zga_e,mass);
//        CacheShiftReal(ip_ex_pp_phi_Zga_llga_ATLAS13_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_phi_Zga_llga_ATLAS13(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_phi_Zga_llga_ATLAS13_cache, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_gg_phi_Zga_llga_ATLAS13_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_gg_phi_Zga_llga,mass);
//        CacheShiftReal(ip_ex_gg_phi_Zga_llga_ATLAS13_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_phi_Zga_llga_ATLAS13_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_phi_Zga_llga_ATLAS13_cache_e, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_gg_phi_Zga_llga_ATLAS13_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_gg_phi_Zga_llga_e,mass);
//        CacheShiftReal(ip_ex_gg_phi_Zga_llga_ATLAS13_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_H_ZZ_llllnunu_ATLAS13(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_H_ZZ_llllnunu_ATLAS13_cache, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_gg_H_ZZ_llllnunu_ATLAS13_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_gg_H_ZZ_llllnunu,mass);
//        CacheShiftReal(ip_ex_gg_H_ZZ_llllnunu_ATLAS13_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_H_ZZ_llllnunu_ATLAS13_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_H_ZZ_llllnunu_ATLAS13_cache_e, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_gg_H_ZZ_llllnunu_ATLAS13_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_gg_H_ZZ_llllnunu_e,mass);
//        CacheShiftReal(ip_ex_gg_H_ZZ_llllnunu_ATLAS13_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_cache, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_VBF_H_ZZ_llllnunu,mass);
//        CacheShiftReal(ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_cache_e, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_VBF_H_ZZ_llllnunu_e,mass);
//        CacheShiftReal(ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_H_ZZ_llnunu_ATLAS13(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_H_ZZ_llnunu_ATLAS13_cache, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_gg_H_ZZ_llnunu_ATLAS13_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_gg_H_ZZ_llnunu,mass);
//        CacheShiftReal(ip_ex_gg_H_ZZ_llnunu_ATLAS13_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_H_ZZ_llnunu_ATLAS13_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_H_ZZ_llnunu_ATLAS13_cache_e, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_gg_H_ZZ_llnunu_ATLAS13_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_gg_H_ZZ_llnunu_e,mass);
//        CacheShiftReal(ip_ex_gg_H_ZZ_llnunu_ATLAS13_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_H_ZZ_llll_ATLAS13(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_H_ZZ_llll_ATLAS13_cache, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_gg_H_ZZ_llll_ATLAS13_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_gg_H_ZZ_llll,mass);
//        CacheShiftReal(ip_ex_gg_H_ZZ_llll_ATLAS13_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_H_ZZ_llll_ATLAS13_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_H_ZZ_llll_ATLAS13_cache_e, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_gg_H_ZZ_llll_ATLAS13_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_gg_H_ZZ_llll_e,mass);
//        CacheShiftReal(ip_ex_gg_H_ZZ_llll_ATLAS13_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_VBF_H_ZZ_llll_ATLAS13(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_VBF_H_ZZ_llll_ATLAS13_cache, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_VBF_H_ZZ_llll_ATLAS13_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_VBF_H_ZZ_llll,mass);
//        CacheShiftReal(ip_ex_VBF_H_ZZ_llll_ATLAS13_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_VBF_H_ZZ_llll_ATLAS13_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_VBF_H_ZZ_llll_ATLAS13_cache_e, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_VBF_H_ZZ_llll_ATLAS13_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_VBF_H_ZZ_llll_e,mass);
//        CacheShiftReal(ip_ex_VBF_H_ZZ_llll_ATLAS13_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_cache, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_gg_H_ZZ_qqllnunu,mass);
//        CacheShiftReal(ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_cache_e, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_gg_H_ZZ_qqllnunu_e,mass);
//        CacheShiftReal(ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_cache, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_VBF_H_ZZ_qqllnunu,mass);
//        CacheShiftReal(ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_cache_e, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_VBF_H_ZZ_qqllnunu_e,mass);
//        CacheShiftReal(ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_H_ZZ_llqq_ATLAS13(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_H_ZZ_llqq_ATLAS13_cache, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_gg_H_ZZ_llqq_ATLAS13_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_gg_H_ZZ_llqq,mass);
//        CacheShiftReal(ip_ex_gg_H_ZZ_llqq_ATLAS13_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_H_ZZ_llqq_ATLAS13_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_H_ZZ_llqq_ATLAS13_cache_e, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_gg_H_ZZ_llqq_ATLAS13_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_gg_H_ZZ_llqq_e,mass);
//        CacheShiftReal(ip_ex_gg_H_ZZ_llqq_ATLAS13_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_VBF_H_ZZ_llqq_ATLAS13(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_VBF_H_ZZ_llqq_ATLAS13_cache, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_VBF_H_ZZ_llqq_ATLAS13_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_VBF_H_ZZ_llqq,mass);
//        CacheShiftReal(ip_ex_VBF_H_ZZ_llqq_ATLAS13_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_VBF_H_ZZ_llqq_ATLAS13_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_VBF_H_ZZ_llqq_ATLAS13_cache_e, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_VBF_H_ZZ_llqq_ATLAS13_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_VBF_H_ZZ_llqq_e,mass);
//        CacheShiftReal(ip_ex_VBF_H_ZZ_llqq_ATLAS13_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_H_ZZ_nunuqq_ATLAS13(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_H_ZZ_nunuqq_ATLAS13_cache, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_gg_H_ZZ_nunuqq_ATLAS13_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_gg_H_ZZ_nunuqq,mass);
//        CacheShiftReal(ip_ex_gg_H_ZZ_nunuqq_ATLAS13_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_H_ZZ_nunuqq_ATLAS13_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_H_ZZ_nunuqq_ATLAS13_cache_e, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_gg_H_ZZ_nunuqq_ATLAS13_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_gg_H_ZZ_nunuqq_e,mass);
//        CacheShiftReal(ip_ex_gg_H_ZZ_nunuqq_ATLAS13_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_H_WW_enumunu_ATLAS13(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_H_WW_enumunu_ATLAS13_cache, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_gg_H_WW_enumunu_ATLAS13_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_gg_H_WW_enumumu,mass);
//        CacheShiftReal(ip_ex_gg_H_WW_enumunu_ATLAS13_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_H_WW_enumunu_ATLAS13_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_H_WW_enumunu_ATLAS13_cache_e, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_gg_H_WW_enumunu_ATLAS13_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_gg_H_WW_enumumu_e,mass);
//        CacheShiftReal(ip_ex_gg_H_WW_enumunu_ATLAS13_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_VBF_H_WW_enumunu_ATLAS13(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_VBF_H_WW_enumunu_ATLAS13_cache, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_VBF_H_WW_enumunu_ATLAS13_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_VBF_H_WW_enumumu,mass);
//        CacheShiftReal(ip_ex_VBF_H_WW_enumunu_ATLAS13_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_VBF_H_WW_enumunu_ATLAS13_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_VBF_H_WW_enumunu_ATLAS13_cache_e, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_VBF_H_WW_enumunu_ATLAS13_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_VBF_H_WW_enumumu_e,mass);
//        CacheShiftReal(ip_ex_VBF_H_WW_enumunu_ATLAS13_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_H_WW_lnuqq_ATLAS13(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_H_WW_lnuqq_ATLAS13_cache, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_gg_H_WW_lnuqq_ATLAS13_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_gg_H_WW_lnuqq,mass);
//        CacheShiftReal(ip_ex_gg_H_WW_lnuqq_ATLAS13_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_gg_H_WW_lnuqq_ATLAS13_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_H_WW_lnuqq_ATLAS13_cache_e, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_gg_H_WW_lnuqq_ATLAS13_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_gg_H_WW_lnuqq_e,mass);
//        CacheShiftReal(ip_ex_gg_H_WW_lnuqq_ATLAS13_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_VBF_H_WW_lnuqq_ATLAS13(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_VBF_H_WW_lnuqq_ATLAS13_cache, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_VBF_H_WW_lnuqq_ATLAS13_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_VBF_H_WW_lnuqq,mass);
//        CacheShiftReal(ip_ex_VBF_H_WW_lnuqq_ATLAS13_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_VBF_H_WW_lnuqq_ATLAS13_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_VBF_H_WW_lnuqq_ATLAS13_cache_e, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_VBF_H_WW_lnuqq_ATLAS13_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_VBF_H_WW_lnuqq_e,mass);
//        CacheShiftReal(ip_ex_VBF_H_WW_lnuqq_ATLAS13_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_pp_H_VV_qqqq_ATLAS13(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_pp_H_VV_qqqq_ATLAS13_cache, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_pp_H_VV_qqqq_ATLAS13_cache[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_pp_H_VV_qqqq,mass);
//        CacheShiftReal(ip_ex_pp_H_VV_qqqq_ATLAS13_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//
//double GeneralTHDMcache::ip_ex_pp_H_VV_qqqq_ATLAS13_e(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_pp_H_VV_qqqq_ATLAS13_cache_e, NumPar, params);
//    if (i>=0) {
//        return(ip_ex_pp_H_VV_qqqq_ATLAS13_cache_e[NumPar][i] );
//    } else {
//        double newResult = interpolate (ATLAS13_pp_H_VV_qqqq_e,mass);
//        CacheShiftReal(ip_ex_pp_H_VV_qqqq_ATLAS13_cache_e, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//
//



double GeneralTHDMcache::ip_ex_pp_Hpm_taunu_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_Hpm_taunu_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_Hpm_taunu_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS8_pp_Hpm_taunu,mass);
        CacheShiftReal(ip_ex_pp_Hpm_taunu_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_ex_pp_Hp_taunu_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_Hp_taunu_CMS8_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_Hp_taunu_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_pp_Hp_taunu,mass);
        CacheShiftReal(ip_ex_pp_Hp_taunu_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_ex_pp_Hpm_taunu_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_Hpm_taunu_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_Hpm_taunu_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_Hpm_taunu,mass);
        CacheShiftReal(ip_ex_pp_Hpm_taunu_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_ex_pp_Hpm_taunu_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_Hpm_taunu_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_Hpm_taunu_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_Hpm_taunu,mass);
        CacheShiftReal(ip_ex_pp_Hpm_taunu_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_ex_pp_Hpm_tb_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_Hpm_tb_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_Hpm_tb_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS8_pp_Hpm_tb,mass);
        CacheShiftReal(ip_ex_pp_Hpm_tb_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_ex_pp_Hp_tb_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_Hp_tb_CMS8_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_Hp_tb_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_pp_Hp_tb,mass);
        CacheShiftReal(ip_ex_pp_Hp_tb_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_pp_Hpm_tb_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_Hpm_tb_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_Hpm_tb_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_Hpm_tb,mass);
        CacheShiftReal(ip_ex_pp_Hpm_tb_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_ex_bsgamma(double logtb, double logmHp){
    int NumPar = 2;
    double params[] = {logtb, logmHp};

    int i = CacheCheckReal(ip_ex_bsgamma_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bsgamma_cache[NumPar][i] );
    } else {
        double newResult = interpolate2D(arraybsgamma, logtb, logmHp);
        CacheShiftReal(ip_ex_bsgamma_cache, NumPar, params, newResult);
        return newResult;
    }
}

// Function needed to calculate some loop functions

gslpp::complex GeneralTHDMcache::f_func(const double x) const{
    if(x<1) {
    gslpp::complex z = -gslpp::complex::i()*M_PI;
    return -pow(log((1.0+sqrt(1.0-x))/(1.0-sqrt(1.0-x)))+z,2)/4.0;
    }
    else {
        return pow(asin(sqrt(1.0/x)),2);
    }
}

gslpp::complex GeneralTHDMcache::g_func(const double x) const{
    if(x<1) {
    gslpp::complex z = -gslpp::complex::i()*M_PI;
    gslpp::complex gs1 = sqrt(1.0-x)*(log((1.0+sqrt(1.0-x))/(1.0-sqrt(1.0-x)))+z)/2.0;
    return gs1;
    }
    else {
        gslpp::complex gg1 = sqrt(x-1.0)*asin(sqrt(1.0/x));
        return gg1;
    }
}

gslpp::complex GeneralTHDMcache::Int1(const double tau, const double lambda) const{
    return tau*lambda/(tau-lambda)/2.0+tau*tau*lambda*lambda/((tau-lambda)
           *(tau-lambda))/2.0*(f_func(tau)-f_func(lambda))+tau*tau*lambda/((tau-lambda)
           *(tau-lambda))*(g_func(tau)-g_func(lambda));
}

gslpp::complex GeneralTHDMcache::Int2(const double tau, const double lambda) const{
    return -tau*lambda/(tau-lambda)/2.0*(f_func(tau)-f_func(lambda));
}

 //Loop functions needed for decays and cross sections


gslpp::complex GeneralTHDMcache::I_h_U(const double mHl2, const double Mu, const double Mc, const double Mt) const {
    int NumPar = 4;
    double params[] = {mHl2, Mu, Mc, Mt};

    int i = CacheCheck(I_h_U_cache, NumPar, params);
    if (i>=0) {
        return ( I_h_U_cache[NumPar][i] );
    } else {
    	double TAUu=4.0*Mu*Mu/mHl2;
    	double TAUc=4.0*Mc*Mc/mHl2;
    	double TAUt=4.0*Mt*Mt/mHl2;
        gslpp::complex newResult = -(8./3.)*(TAUu*(1.0+(1.0-TAUu)*f_func(TAUu))
                         +TAUc*(1.0+(1.0-TAUc)*f_func(TAUc))+TAUt*(1.0+(1.0-TAUt)*f_func(TAUt)));
        CacheShift(I_h_U_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex GeneralTHDMcache::I_HH_U(const double mHh2, const double Mc, const double Mt) const {
    int NumPar = 3;
    double params[] = {mHh2, Mc, Mt};

    int i = CacheCheck(I_HH_U_cache, NumPar, params);
    if (i>=0) {
        return ( I_HH_U_cache[NumPar][i] );
    } else {
    	double TAUc=4.0*Mc*Mc/mHh2;
    	double TAUt=4.0*Mt*Mt/mHh2;
        gslpp::complex newResult = -(8./3.)*(TAUc*(1.0+(1.0-TAUc)*f_func(TAUc))
                      +TAUt*(1.0+(1.0-TAUt)*f_func(TAUt)));
        CacheShift(I_HH_U_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex GeneralTHDMcache::I_A_U(const double mA2, const double Mc, const double Mt) const {
    int NumPar = 3;
    double params[] = {mA2, Mc, Mt};

    int i = CacheCheck(I_A_U_cache, NumPar, params);
    if (i>=0) {
        return ( I_A_U_cache[NumPar][i] );
    } else {
    	double TAUc=4.0*Mc*Mc/mA2;
    	double TAUt=4.0*Mt*Mt/mA2;
        gslpp::complex newResult = -(8./3.)*(TAUc*f_func(TAUc)+TAUt*f_func(TAUt));
        CacheShift(I_A_U_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex GeneralTHDMcache::I_h_D(const double mHl2, const double Md, const double Ms, const double Mb) const {
    int NumPar = 4;
    double params[] = {mHl2, Md, Ms, Mb};

    int i = CacheCheck(I_h_D_cache, NumPar, params);
    if (i>=0) {
        return ( I_h_D_cache[NumPar][i] );
    } else {
    	double TAUd=4.0*Md*Md/mHl2;
    	double TAUs=4.0*Ms*Ms/mHl2;
    	double TAUb=4.0*Mb*Mb/mHl2;
        gslpp::complex newResult = -(2./3.)*(TAUd*(1.0+(1.0-TAUd)*f_func(TAUd))
                         +TAUs*(1.0+(1.0-TAUs)*f_func(TAUs))+TAUb*(1.0+(1.0-TAUb)*f_func(TAUb)));
        CacheShift(I_h_D_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex GeneralTHDMcache::I_HH_D(const double mHh2, const double Ms, const double Mb) const {
    int NumPar = 3;
    double params[] = {mHh2, Ms, Mb};

    int i = CacheCheck(I_HH_D_cache, NumPar, params);
    if (i>=0) {
        return ( I_HH_D_cache[NumPar][i] );
    } else {
    	double TAUs=4.0*Ms*Ms/mHh2;
    	double TAUb=4.0*Mb*Mb/mHh2;
        gslpp::complex newResult = -(2./3.)*(TAUs*(1.0+(1.0-TAUs)*f_func(TAUs))
                      +TAUb*(1.0+(1.0-TAUb)*f_func(TAUb)));
        CacheShift(I_HH_D_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex GeneralTHDMcache::I_A_D(const double mA2, const double Ms, const double Mb) const {
    int NumPar = 3;
    double params[] = {mA2, Ms, Mb};

    int i = CacheCheck(I_A_D_cache, NumPar, params);
    if (i>=0) {
        return ( I_A_D_cache[NumPar][i] );
    } else {
    	double TAUs=4.0*Ms*Ms/mA2;
    	double TAUb=4.0*Mb*Mb/mA2;
        gslpp::complex newResult = -(2./3.)*(TAUs*f_func(TAUs)+TAUb*f_func(TAUb));
        CacheShift(I_A_D_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex GeneralTHDMcache::I_h_L(const double mHl2, const double Me, const double Mmu, const double Mtau) const {
    int NumPar = 4;
    double params[] = {mHl2, Me, Mmu, Mtau};

    int i = CacheCheck(I_h_L_cache, NumPar, params);
    if (i>=0) {
        return ( I_h_L_cache[NumPar][i] );
    } else {
    	double TAUe=4.0*Me*Me/mHl2;
    	double TAUmu=4.0*Mmu*Mmu/mHl2;
    	double TAUtau=4.0*Mtau*Mtau/mHl2;
        gslpp::complex newResult = -2.0*(TAUe*(1.0+(1.0-TAUe)*f_func(TAUe))
                         +TAUmu*(1.0+(1.0-TAUmu)*f_func(TAUmu))
                         +TAUtau*(1.0+(1.0-TAUtau)*f_func(TAUtau)));
        CacheShift(I_h_L_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex GeneralTHDMcache::I_HH_L(const double mHh2, const double Mmu, const double Mtau) const {
    int NumPar = 3;
    double params[] = {mHh2, Mmu, Mtau};

    int i = CacheCheck(I_HH_L_cache, NumPar, params);
    if (i>=0) {
        return ( I_HH_L_cache[NumPar][i] );
    } else {
    	double TAUmu=4.0*Mmu*Mmu/mHh2;
    	double TAUtau=4.0*Mtau*Mtau/mHh2;
        gslpp::complex newResult = -2.0*(TAUmu*(1.0+(1.0-TAUmu)*f_func(TAUmu))+
                           TAUtau*(1.0+(1.0-TAUtau)*f_func(TAUtau)));
        CacheShift(I_HH_L_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex GeneralTHDMcache::I_A_L(const double mA2, const double Mmu, const double Mtau) const {
    int NumPar = 3;
    double params[] = {mA2, Mmu, Mtau};

    int i = CacheCheck(I_A_L_cache, NumPar, params);
    if (i>=0) {
        return ( I_A_L_cache[NumPar][i] );
    } else {
    	double TAUmu=4.0*Mmu*Mmu/mA2;
    	double TAUtau=4.0*Mtau*Mtau/mA2;
        gslpp::complex newResult = -2.0*(TAUmu*f_func(TAUmu)+TAUtau*f_func(TAUtau));
        CacheShift(I_A_L_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex GeneralTHDMcache::I_H_W(const double mH, const double MW) const {
    int NumPar = 2;
    double params[] = {mH, MW};

    int i = CacheCheck(I_H_W_cache, NumPar, params);
    if (i>=0) {
        return ( I_H_W_cache[NumPar][i] );
    } else {
        double TAUw=4.0*MW*MW/(mH*mH);
        gslpp::complex newResult = 2.0 + 3.0*TAUw + 3.0*TAUw*(2.0-TAUw)*f_func(TAUw);
        CacheShift(I_H_W_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex GeneralTHDMcache::I_H_Hp(const double mHp2, const double mH) const {
    int NumPar = 2;
    double params[] = {mHp2, mH};

    int i = CacheCheck(I_H_Hp_cache, NumPar, params);
    if (i>=0) {
        return ( I_H_Hp_cache[NumPar][i] );
    } else {
        double TAUhp=4.0*mHp2/(mH*mH);
        gslpp::complex newResult = -TAUhp*(1.0-TAUhp*f_func(TAUhp));
        CacheShift(I_H_Hp_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex GeneralTHDMcache::A_h_U(const double mHl2, const double cW2, const double Mu, const double Mc, const double Mt, const double MZ) const {
    int NumPar = 6;
    double params[] = {mHl2, cW2, Mu, Mc, Mt, MZ};

    int i = CacheCheck(A_h_U_cache, NumPar, params);
    if (i>=0) {
        return ( A_h_U_cache[NumPar][i] );
    } else {
    	double TAUu=4.0*Mu*Mu/mHl2;
    	double TAUc=4.0*Mc*Mc/mHl2;
    	double TAUt=4.0*Mt*Mt/mHl2;
    	double LAMu=4.0*Mu*Mu/(MZ*MZ);
    	double LAMc=4.0*Mc*Mc/(MZ*MZ);
    	double LAMt=4.0*Mt*Mt/(MZ*MZ);
	double sW2=1.0-cW2;
        gslpp::complex newResult = -4.0*(1.0/2.0-4.0/3.0*sW2)*(Int1(TAUu,LAMu)+Int1(TAUc,LAMc)
                           +Int1(TAUt,LAMt)-Int2(TAUu,LAMu)-Int2(TAUc,LAMc)-Int2(TAUt,LAMt));
        CacheShift(A_h_U_cache, NumPar, params, newResult);
        return newResult;
    }

}

gslpp::complex GeneralTHDMcache::A_HH_U(const double mHh2, const double cW2, const double Mc, const double Mt, const double MZ) const {
    int NumPar = 5;
    double params[] = {mHh2, cW2, Mc, Mt, MZ};

    int i = CacheCheck(A_HH_U_cache, NumPar, params);
    if (i>=0) {
        return ( A_HH_U_cache[NumPar][i] );
    } else {
    	double TAUc=4.0*Mc*Mc/mHh2;
    	double TAUt=4.0*Mt*Mt/mHh2;
    	double LAMc=4.0*Mc*Mc/(MZ*MZ);
    	double LAMt=4.0*Mt*Mt/(MZ*MZ);
	double sW2=1.0-cW2;
        gslpp::complex newResult = -4.0*(1.0/2.0-4.0/3.0*sW2)*(Int1(TAUc,LAMc)-Int2(TAUc,LAMc)
                                         +Int1(TAUt,LAMt)-Int2(TAUt,LAMt));
        CacheShift(A_HH_U_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex GeneralTHDMcache::A_A_U(const double mA2, const double cW2, const double Mc, const double Mt, const double MZ) const {
    int NumPar = 5;
    double params[] = {mA2, cW2, Mc, Mt, MZ};

    int i = CacheCheck(A_A_U_cache, NumPar, params);
    if (i>=0) {
        return ( A_A_U_cache[NumPar][i] );
    } else {
    	double TAUc=4.0*Mc*Mc/mA2;
    	double TAUt=4.0*Mt*Mt/mA2;
    	double LAMc=4.0*Mc*Mc/(MZ*MZ);
    	double LAMt=4.0*Mt*Mt/(MZ*MZ);
	double sW2=1.0-cW2;
        gslpp::complex newResult = -4.0*(1.0/2.0-4.0/3.0*sW2)*(-Int2(TAUc,LAMc)-Int2(TAUt,LAMt))/sqrt(sW2*cW2);
        CacheShift(A_A_U_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex GeneralTHDMcache::A_h_D(const double mHl2, const double cW2, const double Md, const double Ms, const double Mb, const double MZ) const {
    int NumPar = 6;
    double params[] = {mHl2, cW2, Md, Ms, Mb, MZ};

    int i = CacheCheck(A_h_D_cache, NumPar, params);
    if (i>=0) {
        return ( A_h_D_cache[NumPar][i] );
    } else {
    	double TAUd=4.0*Md*Md/mHl2;
    	double TAUs=4.0*Ms*Ms/mHl2;
    	double TAUb=4.0*Mb*Mb/mHl2;
    	double LAMd=4.0*Md*Md/(MZ*MZ);
    	double LAMs=4.0*Ms*Ms/(MZ*MZ);
	double LAMb=4.0*Mb*Mb/(MZ*MZ);
	double sW2=1.0-cW2;
        gslpp::complex newResult = 2.0*(-1.0/2.0+2.0/3.0*sW2)*(Int1(TAUd,LAMd)+Int1(TAUs,LAMs)
                           +Int1(TAUb,LAMb)-Int2(TAUd,LAMd)-Int2(TAUs,LAMs)-Int2(TAUb,LAMb));
        CacheShift(A_h_D_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex GeneralTHDMcache::A_HH_D(const double mHh2, const double cW2, const double Ms, const double Mb, const double MZ) const {
    int NumPar = 5;
    double params[] = {mHh2, cW2, Ms, Mb, MZ};

    int i = CacheCheck(A_HH_D_cache, NumPar, params);
    if (i>=0) {
        return ( A_HH_D_cache[NumPar][i] );
    } else {
    	double TAUs=4.0*Ms*Ms/mHh2;
    	double TAUb=4.0*Mb*Mb/mHh2;
    	double LAMs=4.0*Ms*Ms/(MZ*MZ);
	double LAMb=4.0*Mb*Mb/(MZ*MZ);
	double sW2=1.0-cW2;
        gslpp::complex newResult = 2.0*(-1.0/2.0+2.0/3.0*sW2)*(Int1(TAUs,LAMs)-Int2(TAUs,LAMs)
                                          +Int1(TAUb,LAMb)-Int2(TAUb,LAMb));
        CacheShift(A_HH_D_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex GeneralTHDMcache::A_A_D(const double mA2, const double cW2, const double Ms, const double Mb, const double MZ) const {
    int NumPar = 5;
    double params[] = {mA2, cW2, Ms, Mb, MZ};

    int i = CacheCheck(A_A_D_cache, NumPar, params);
    if (i>=0) {
        return ( A_A_D_cache[NumPar][i] );
    } else {
    	double TAUs=4.0*Ms*Ms/mA2;
    	double TAUb=4.0*Mb*Mb/mA2;
    	double LAMs=4.0*Ms*Ms/(MZ*MZ);
	double LAMb=4.0*Mb*Mb/(MZ*MZ);
	double sW2=1.0-cW2;
        gslpp::complex newResult = 2.0*(-1.0/2.0+2.0/3.0*sW2)*(-Int2(TAUs,LAMs)-Int2(TAUb,LAMb))/sqrt(sW2*cW2);
        CacheShift(A_A_D_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex GeneralTHDMcache::A_h_L(const double mHl2, const double cW2, const double Me, const double Mmu, const double Mtau, const double MZ) const {
    int NumPar = 6;
    double params[] = {mHl2, cW2, Me, Mmu, Mtau, MZ};

    int i = CacheCheck(A_h_L_cache, NumPar, params);
    if (i>=0) {
        return ( A_h_L_cache[NumPar][i] );
    } else {
    	double TAUe=4.0*Me*Me/mHl2;
    	double TAUmu=4.0*Mmu*Mmu/mHl2;
    	double TAUtau=4.0*Mtau*Mtau/mHl2;
    	double LAMe=4.0*Me*Me/(MZ*MZ);
    	double LAMmu=4.0*Mmu*Mmu/(MZ*MZ);
	double LAMtau=4.0*Mtau*Mtau/(MZ*MZ);
	double sW2=1.0-cW2;
        gslpp::complex newResult = 2.0*(-1.0/2.0+2.0*sW2)*(Int1(TAUe,LAMe)+Int1(TAUmu,LAMmu)
                            +Int1(TAUtau,LAMtau)-Int2(TAUe,LAMe)-Int2(TAUmu,LAMmu)
                            -Int2(TAUtau,LAMtau));
        CacheShift(A_h_L_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex GeneralTHDMcache::A_HH_L(const double mHh2, const double cW2, const double Mmu, const double Mtau, const double MZ) const {
    int NumPar = 5;
    double params[] = {mHh2, cW2, Mmu, Mtau, MZ};

    int i = CacheCheck(A_HH_L_cache, NumPar, params);
    if (i>=0) {
        return ( A_HH_L_cache[NumPar][i] );
    } else {
    	double TAUmu=4.0*Mmu*Mmu/mHh2;
    	double TAUtau=4.0*Mtau*Mtau/mHh2;
    	double LAMmu=4.0*Mmu*Mmu/(MZ*MZ);
	double LAMtau=4.0*Mtau*Mtau/(MZ*MZ);
	double sW2=1.0-cW2;
        gslpp::complex newResult = 2.0*(-1.0/2.0+2.0*sW2)*(Int1(TAUmu,LAMmu)-Int2(TAUmu,LAMmu)
                                      +Int1(TAUtau,LAMtau)-Int2(TAUtau,LAMtau));
        CacheShift(A_HH_L_cache, NumPar, params, newResult);
        return newResult;
    }
}   
    


gslpp::complex GeneralTHDMcache::A_A_L(const double mA2, const double cW2, const double Mmu, const double Mtau, const double MZ) const {
    int NumPar = 5;
    double params[] = {mA2, cW2, Mmu, Mtau, MZ};

    int i = CacheCheck(A_A_L_cache, NumPar, params);
    if (i>=0) {
        return ( A_A_L_cache[NumPar][i] );
    } else {
    	double TAUmu=4.0*Mmu*Mmu/mA2;
    	double TAUtau=4.0*Mtau*Mtau/mA2;
    	double LAMmu=4.0*Mmu*Mmu/(MZ*MZ);
	double LAMtau=4.0*Mtau*Mtau/(MZ*MZ);
	double sW2=1.0-cW2;
        gslpp::complex newResult = 2.0*(-1.0/2.0+2.0*sW2)*(-Int2(TAUmu,LAMmu)-Int2(TAUtau,LAMtau))/sqrt(sW2*cW2);
        CacheShift(A_A_L_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex GeneralTHDMcache::A_H_W(const double mH, const double cW2, const double MW, const double MZ) const {
    int NumPar = 4;
    double params[] = {mH, cW2, MW, MZ};

    int i = CacheCheck(A_H_W_cache, NumPar, params);
    if (i>=0) {
        return ( A_H_W_cache[NumPar][i] );
    } else {
        double TAUw=4.0*MW*MW/(mH*mH);
        double LAMw=4.0*MW*MW/(MZ*MZ);
	double sW2=1.0-cW2;
        gslpp::complex newResult = -sqrt(cW2/sW2)*(4.0*(3.0-sW2/cW2)*Int2(TAUw,LAMw)
                            +((1.0+2.0/TAUw)*sW2/cW2-(5.0+2.0/TAUw))*Int1(TAUw,LAMw));
        CacheShift(A_H_W_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex GeneralTHDMcache::A_H_Hp(const double mHp2, const double mH, const double cW2, const double MZ) const {
    int NumPar = 4;
    double params[] = {mHp2, mH, cW2, MZ};

    int i = CacheCheck(A_H_Hp_cache, NumPar, params);
    if (i>=0) {
        return ( A_H_Hp_cache[NumPar][i] );
    } else {
        double TAUhp=4.0*mHp2/(mH*mH);
        double LAMhp=4.0*mHp2/(MZ*MZ);
	double sW2=1.0-cW2;
        gslpp::complex newResult = (1.0-2.0*sW2)/sqrt(cW2*sW2)*Int1(TAUhp,LAMhp);
        CacheShift(A_H_Hp_cache, NumPar, params, newResult);
        return newResult;
    }
}

int GeneralTHDMcache::HSTheta (const double x) const{
    if(x<0) return 0.0;
    else return 1.0;
}


double GeneralTHDMcache::KaellenFunction(const double a2, const double b2, const double c2) const{
    int NumPar = 3;
    double params[] = {a2, b2, c2};

    int i = CacheCheckReal(KaellenFunction_cache, NumPar, params);
    if (i>=0) {
        return ( KaellenFunction_cache[NumPar][i] );
    }
    else {
        double newResult = 0.0;
        double x = (a2-b2-c2)*(a2-b2-c2)-4.0*b2*c2;
        if(x>0) newResult = sqrt(std::fabs(x/a2))/2.0;        
        CacheShiftReal(KaellenFunction_cache, NumPar, params, newResult);
        return newResult;
    }
}

    
    double GeneralTHDMcache::cW2GTHDM(const double c02) const{
    return c02;
    }



    double GeneralTHDMcache::MWGTHDM(const double MW) const{
    return MW;
    }



    double GeneralTHDMcache::beta(const double mf, const double m_2) const
    {
            return sqrt(1.0-4.0*mf*mf/m_2);
    }
    
    
    double GeneralTHDMcache::beta_mt_sq(const double mt,const double m_2) const
    {
        if (4.0*mt*mt/m_2 < 1 )
        {    
            return 1.0-4.0*mt*mt/m_2;
        }
        else
        {    
            return -(1.0-4.0*mt*mt/m_2);
        }
    }
  

     double GeneralTHDMcache::lambdaijk(const double Ri1,const double Ri2,const double Ri3,const double Rj1,const double Rj2,const double Rj3, const double Rk1,const double Rk2,const double Rk3, const double lambda1H, const double lambda3H, const double lambda4H, const double Relambda5H, const double Imlambda5H, const double Relambda6H, const double Imlambda6H, const double Relambda7H, const double Imlambda7H) const
    {
        return (1.0/2.0)*vev*(Imlambda7H*(-Ri3*Rj3*Rk3 - Ri2*Rj2*Rk3) - 3.0*Imlambda6H*Ri1*Rj1*Rk3 
                + lambda1H*Ri1*Rj1*Rk1 + Relambda7H*Ri2*Rj2*Rk2 + 3.0*Relambda6H*Ri1*Rj1*Rk2
                +(Relambda5H + lambda3H + lambda4H)*Ri1*Rj2*Rk2 - (2.0*Relambda5H - lambda3H - lambda4H)*Ri1*Rj3*Rk3
                + Relambda7H*Ri2*Rj3*Rk3 - Imlambda5H*Ri1*Rj2*Rk3);
    }
    
    double GeneralTHDMcache::lambdaipm(const double Ri1,const double Ri2, const double Ri3) const
    {
        return -vev*(lambda3H*Ri1 + Relambda7H*Ri2 - Imlambda7H*Ri3);
    }
       
void GeneralTHDMcache::computeSignalStrengths()
{
    
    m2_2 = mH2sq;
    m2 = sqrt(m2_2);
    m3_2 = mH3sq;
    m3 = sqrt(m3_2);
    
    
    double GF = 1/(sqrt(2.0)*vev*vev);
    double sW2=1.0-cW2;

      //FLAG to select only the model in which all the couplings are the same (by families)

    if (!myGTHDM->getATHDMflag())
    {
        throw std::runtime_error("Signal strengths are only available in the A2HDM.");
    }
  
        /*complex i */
    gslpp::complex i = gslpp::complex::i();
     
    Mt=myGTHDM->getQuarks(QCD::TOP).getMass();
    Mb=myGTHDM->getQuarks(QCD::BOTTOM).getMass(); 
    Mtau=myGTHDM->getLeptons(StandardModel::TAU).getMass();
    Mc=myGTHDM->getQuarks(QCD::CHARM).getMass();
    Ms=myGTHDM->getQuarks(QCD::STRANGE).getMass();
    Mmu=myGTHDM->getLeptons(StandardModel::MU).getMass();
    Mu=myGTHDM->getQuarks(QCD::UP).getMass();
    Md=myGTHDM->getQuarks(QCD::DOWN).getMass();
    Me=myGTHDM->getLeptons(StandardModel::ELECTRON).getMass();
    MW=MWGTHDM(myGTHDM->Mw_tree());
    cW2=cW2GTHDM(myGTHDM->c02());
//    MZ=myGTHDM->getMz();

   
    
    //The 125 GeV is always defined as the one of m_1, so we don't use the mass ordering. 
   // For the SM_Higgs flag it does not matter
  
    double m1_2 = mH1sq;
    double m1 = mHl;
  
    //fermionic couplings for phi1
    
    gslpp::complex yu1 = myGTHDM->getyu1();
    gslpp::complex yd1 = myGTHDM->getyd1();
    gslpp::complex yl1 = myGTHDM->getyl1();
    
    yu1R = myGTHDM->getyu1R();
    yd1R = myGTHDM->getyd1R();
    yl1R = myGTHDM->getyl1R();
   
        
     //The Standard Model h branching ratios

    BrSM_htobb = 5.77e-1;
    BrSM_htotautau = 6.32e-2;
    BrSM_htogaga = 2.28e-3;
    double BrSM_htoWW = 2.15e-1;
    double BrSM_htoZZ = 2.64e-2;
    double BrSM_htogg = 8.57e-2;
    double BrSM_htoZga = 1.54e-3;
    double BrSM_htocc = 2.91e-2;

    //The ggH cross section in the SM at 8 TeV
    double SigmaggF8 = myGTHDM->computeSigmaggH(8.0);
    //The ggH cross section in the SM at 13 TeV.
    double SigmaggF13 = myGTHDM->computeSigmaggH(13.0);
    //The square of the top-quark contribution to the ggH cross section in the SM at 8 TeV
    double Sigmaggh_tt8 = myGTHDM->computeSigmaggH_tt(8.0);
    //The square of the top-quark contribution to the ggH cross section in the SM at 13 TeV
//    double Sigmaggh_tt13 = myTHDM->computeSigmaggH_tt(13.0);
    //The square of the bottom-quark contribution to the ggH cross section in the SM at 8 TeV
    double Sigmaggh_bb8 = myGTHDM->computeSigmaggH_bb(8.0);
    //The square of the bottom-quark contribution to the ggH cross section in the SM at 13 TeV
//    double Sigmaggh_bb13 = myTHDM->computeSigmaggH_bb(13.0);
    //The ttH production cross section in the SM at 8 TeV
    double Sigmatth8 = myGTHDM->computeSigmattH(8.0);
    //The ttH production cross section in the SM at 13 TeV
    double Sigmatth13 = myGTHDM->computeSigmattH(13.0);
    //The bbH production cross section in the SM at 13 TeV
    double Sigmabbh13 = ip_cs_pptobbH_13(mHl);
    //The VBF plus Vh production cross section in the SM at 13 TeV
    double SigmaVBFVh13 = (myGTHDM->computeSigmaVBF(13.0)+myGTHDM->computeSigmaWH(13.0)+myGTHDM->computeSigmaZH(13.0));

    
 //////gg
    
    
    //gg -> A (phi odd) production cross section ratio at 8 TeV, top loop only over total
   // double rSigmagghO_t8 = ip_csr_ggA_t_8(m1);
    //gg -> A (phi odd) production cross section ratio at 8 TeV, bottom loop only over total
    double rSigmagghO_b8 = ip_csr_ggA_b_8(m1);
    
    //gg -> A (phiodd) production cross section at 8 TeV, total
    double SigmagghO_8 = ip_cs_ggtoA_8(m1);
    
     
    //gg -> A (phiodd) production cross section at 13 TeV, total
  //  double SigmagghO_13 = ip_cs_ggtoA_13(m1);

  //  beta_h_t = beta(Mt, m1_2);
    beta_h_t = sqrt(-1.0+4.0*Mt*Mt/m1_2);
    beta_h_b = beta(Mb, m1_2);
    beta_h_tau = beta(Mtau, m1_2);
    beta_h_c = beta(Mc, m1_2);
    beta_h_mu = beta(Mmu, m1_2);
     
    /* r_ii is the ratio of the squared 2HDM vertex coupling of h to
     * the particle i and the respective squared SM coupling. Where separated
       E means the even part (coming from the CP-even scalar) and O
       the odd part (coming from the CP-odd scalar) */
                 
    rh_QuQuE= yu1.real()*yu1.real(); 
    rh_QuQuO= yu1.imag()*yu1.imag(); 
    rh_QdQdE= yd1.real()*yd1.real(); 
    rh_QdQdO= yd1.imag()*yd1.imag(); 
    rh_QlQlE= yl1.real()*yl1.real(); 
    rh_QlQlO= yl1.imag()*yl1.imag(); 
    rh_ggE = yu1.real()*yd1.real() + (yu1.real()*yu1.real() - yu1.real()*yd1.real())*(Sigmaggh_tt8/SigmaggF8)  + (yd1.real()*yd1.real() - yu1.real()*yd1.real())*(Sigmaggh_bb8/SigmaggF8);
    rh_ggO = yu1.imag()*yu1.imag() + (yu1.imag()*yu1.imag() - yu1.imag()*yd1.imag())*rSigmagghO_b8  + (yd1.imag()*yd1.imag() - yu1.imag()*yd1.imag())*rSigmagghO_b8;
    rh_gg = rh_ggE+rh_ggO*(SigmagghO_8/SigmaggF8);
    rh_VV=0.0;
    if(myGTHDM->getSMHiggs()){
            rh_VV=R11_GTHDM*R11_GTHDM;
    }
    else{
            rh_VV=R21_GTHDM*R21_GTHDM;
    }        
            
            
    
    /*Loop functions needed to rh_gaga and rh_Zga ...*/
    
    
    gslpp::complex fermU=I_h_U(m1_2,Mu,Mc,Mt);
    gslpp::complex fermD=I_h_D(m1_2,Md,Ms,Mb);
    gslpp::complex fermL=I_h_L(m1_2,Me,Mmu,Mtau);
    gslpp::complex I_hSM_W=I_H_W(mHl,MW);
    gslpp::complex I_h_W=0.0;
      if(myGTHDM->getSMHiggs()){
            I_h_W=R11_GTHDM*I_hSM_W;
    }
    else{
            I_h_W=R21_GTHDM*I_hSM_W;
    }        
     
  
    gslpp::complex I_hSM_F= fermU+fermD+fermL;
    gslpp::complex I_hE_F= yu1.real()*fermU+ yd1.real()*fermD+yl1.real()*fermL;
    
                                                                               
    /*Coupling between h and two charged Higgs*/
    
       
    double lambdahHpHm = 0.0; 
    if(myGTHDM->getSMHiggs()){
         lambdahHpHm = lambdaipm(R11_GTHDM, R12_GTHDM, R13_GTHDM);
    }
    else{
         lambdahHpHm = lambdaipm(R21_GTHDM, R22_GTHDM, R23_GTHDM);
    }
       
    gslpp::complex I_h_Hp=(vev)/(2.0*mHp2)*I_H_Hp(mHp2,m1)*lambdahHpHm;
       
    /*CP ODD */

    gslpp::complex I_h_Ux=I_A_U(m1_2,Mc,Mt);
    gslpp::complex I_h_Dx=I_A_D(m1_2,Ms,Mb);
    gslpp::complex I_h_Lx=I_A_L(m1_2,Mmu,Mtau);

    gslpp::complex I_hO_F = yu1.imag()*I_h_Ux + yd1.imag()*I_h_Dx + yl1.imag()*I_h_Lx;

   // double Gamma_hgaga=(GF*Ale*Ale*m1*m1*m1/(sqrt(2.0)*128.0*M_PI*M_PI*M_PI))*((I_hE_F+I_h_W+I_h_Hp).abs2()+ (I_hO_F).abs2());
    rh_gaga = ((I_hE_F+I_h_W+I_h_Hp).abs2()+ (I_hO_F).abs2())/(I_hSM_F +I_hSM_W).abs2();    
    /*Decay to Z gamma
    CP-EVEN PART*/

    gslpp::complex A_hE_Ux = A_h_U(m1_2,cW2,Mu,Mc,Mt,MZ);
    gslpp::complex A_hE_Dx = A_h_D(m1_2,cW2,Md,Ms,Mb,MZ);
    gslpp::complex A_hE_Lx = A_h_L(m1_2,cW2,Me,Mmu,Mtau,MZ);
    gslpp::complex A_hSM_W = A_H_W(m1,cW2,MW,MZ);
    gslpp::complex A_h_W = 0.0;
    if(myGTHDM->getSMHiggs()){
         A_h_W = R11_GTHDM*A_hSM_W;
    }
    else{
         A_h_W = R21_GTHDM*A_hSM_W;
    }
    
    

    gslpp::complex A_hSM_F = (A_hE_Ux+ A_hE_Dx+ A_hE_Lx)/sqrt(sW2*cW2);
    gslpp::complex A_hE_F = (yu1.real()*A_hE_Ux+ yd1.real()*A_hE_Dx+ yl1.real()*A_hE_Lx)/sqrt(sW2*cW2);

    gslpp::complex A_h_Hp =(vev)/(2.0*mHp2)*A_H_Hp(mHp2,m1,cW2,MZ)*(lambdahHpHm);

    /*CP-ODD PART*/

    gslpp::complex A_hO_Ux = A_A_U(m1_2,cW2,Mc,Mt,MZ);
    gslpp::complex A_hO_Dx = A_A_D(m1_2,cW2,Ms,Mb,MZ);
    gslpp::complex A_hO_Lx = A_A_L(m1_2,cW2,Mmu,Mtau,MZ);

    gslpp::complex A_hO_F=yu1.imag()*A_hO_Ux + yd1.imag()*A_hO_Dx + yl1.imag()*A_hO_Lx;

   // double Gamma_hZga=HSTheta(m1-MZ)*GF*Ale*Ale*m1*m1*m1/(sqrt(2.0)*64.0*M_PI*M_PI*M_PI)*(1.0-MZ*MZ/(m1*m1))*(1.0-MZ*MZ/(m1*m1))*(1.0-MZ*MZ/(m1*m1))*((A_hE_F+A_h_W+A_h_Hp).abs2()+ A_hO_F.abs2());
    rh_Zga = ((A_hE_F+A_h_W+A_h_Hp).abs2()+ A_hO_F.abs2())/(A_hSM_F +A_hSM_W ).abs2();

    /*Decay to gluons*/

    double Gamma_hggSM=GF*Als*Als*m1*m1*m1/(sqrt(2.0)*16.0*M_PI*M_PI*M_PI)*(9.0/4.0)*(fermU/4.0+fermD).abs2();
    
    double Gamma_hgg=rh_gg*GF*Als*Als*m1*m1*m1/(sqrt(2.0)*16.0*M_PI*M_PI*M_PI)*(9.0/4.0)*(fermU/4.0+fermD).abs2();
    double lambda122 = 0;
   if(myGTHDM->getSMHiggs()){
         lambda122 = (2.0)*(lambdaijk(R11, R21, R31, R12, R22, R32, R12, R22, R32, lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H) +  lambdaijk(R12, R22, R32, R11, R21, R31, R12, R22, R32,  lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H) + lambdaijk(R12, R22, R32, R12, R22, R32, R11, R21, R31,  lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H) );
    }
    else{
         lambda122 = (2.0)*(lambdaijk(R12, R22, R32, R11, R21, R31, R11, R21, R31, lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H) +  lambdaijk(R11, R21, R31, R12, R22, R32, R11, R21, R31,  lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H) + lambdaijk(R11, R21, R31, R11, R21, R31, R12, R22, R32,  lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H) );
    }
    
    
    
    
    double Gamma_hHH = HSTheta(m1 - 2.0*m2)*(KaellenFunction(m1_2,m2_2,m2_2)*lambda122*lambda122)/(32.0*M_PI);
     
   double lambda133 = 0;
    if(myGTHDM->getSMHiggs()){
     lambda133 = (vev)*(Relambda7H*R21 - (2.0*Relambda5H - lambda3H - lambda4H)*R11);
    }
    else{
     lambda133 = (vev)*(Relambda7H*R22 - (2.0*Relambda5H - lambda3H - lambda4H)*R12);
    }
   
   
   
   double Gamma_hAA = HSTheta(m1 - 2.0*m3)*(KaellenFunction(m1_2,m3_2,m3_2)*lambda133*lambda133)/(32.0*M_PI);;

   
   //  /* ggF_tth8 is the ratio of the THDM and SM cross sections for ggF or tth production at 8 TeV*/
  //  ggF_tth8 = (SigmaggF8*rh_ggE + SigmagghO_8*rh_ggO + Sigmatth8*(rh_QuQuE + rh_QuQuO/(beta(Mc, m1_2)*beta(Mc, m1_2))))/(SigmaggF8 + Sigmatth8);
  //  /* ggF_tth13 is the ratio of the THDM and SM cross sections for ggF or tth production at 13 TeV */
  //  ggF_tth13 = (SigmaggF13*rh_ggE + SigmagghO_13*rh_ggO + Sigmatth8*(rh_QuQuE + rh_QuQuO/(beta(Mc, m1_2)*beta(Mc, m1_2))))/(SigmaggF13 + Sigmatth13);
  //  /* pph13 is the ratio of the THDM and SM cross sections for an h production at 13 TeV */
  //  pph13 = (SigmaggF13*rh_ggE + SigmagghO_13*rh_ggO+ SigmaVBFVh13*rh_VV + Sigmatth13*(rh_QuQuE + rh_QuQuO/(beta(Mc, m1_2)*beta(Mc, m1_2))) + Sigmabbh13*(rh_QdQdE + rh_QdQdE/(beta(Mb, m1_2)*beta(Mb, m1_2))))/(SigmaggF13 + SigmaVBFVh13 + Sigmatth13 + Sigmabbh13);
  //  /* VBF_Vh is the ratio of the THDM and SM cross sections for VBF or Vh production */
  //  VBF_Vh = rh_VV;

    /* ggF_tth8 is the ratio of the THDM and SM cross sections for ggF or tth production at 8 TeV*/
    ggF_tth8 = (SigmaggF8*rh_gg + Sigmatth8*(rh_QuQuE + rh_QuQuO/(beta(Mc, m1_2)*beta(Mc, m1_2))))/(SigmaggF8 + Sigmatth8);
    /* ggF_tth13 is the ratio of the THDM and SM cross sections for ggF or tth production at 13 TeV */
    ggF_tth13 = (SigmaggF13*rh_gg + Sigmatth8*(rh_QuQuE + rh_QuQuO/(beta(Mc, m1_2)*beta(Mc, m1_2))))/(SigmaggF13 + Sigmatth13);
    /* pph13 is the ratio of the THDM and SM cross sections for an h production at 13 TeV */
    pph13 = (SigmaggF13*rh_gg+ SigmaVBFVh13*rh_VV + Sigmatth13*(rh_QuQuE + rh_QuQuO/(beta(Mc, m1_2)*beta(Mc, m1_2))) + Sigmabbh13*(rh_QdQdE + rh_QdQdE/(beta(Mb, m1_2)*beta(Mb, m1_2))))/(SigmaggF13 + SigmaVBFVh13 + Sigmatth13 + Sigmabbh13);
    /* VBF_Vh is the ratio of the THDM and SM cross sections for VBF or Vh production */
    VBF_Vh = rh_VV;

    sumModBRs = BrSM_htobb*(rh_QdQdE + rh_QdQdO/(beta(Mb, m1_2)*beta(Mb, m1_2))) 
            + (BrSM_htoWW+BrSM_htoZZ)*rh_VV
            + BrSM_htotautau*(rh_QlQlE + rh_QlQlO/(beta(Mtau, m1_2)*beta(Mtau, m1_2))) 
            + BrSM_htogaga*rh_gaga
            + BrSM_htogg*rh_gg
            + BrSM_htoZga*rh_Zga
            + BrSM_htocc*(rh_QuQuE + rh_QuQuO/(beta(Mc, m1_2)*beta(Mc, m1_2)));
    Gamma_h = sumModBRs*myGTHDM->computeGammaHTotal() + Gamma_hHH + Gamma_hAA;
   
    GTHDM_BR_h_bb=(rh_QdQdE + rh_QdQdO/(beta(Mb, m1_2)*beta(Mb, m1_2)))*BrSM_htobb/sumModBRs;
    GTHDM_BR_h_WW = rh_VV*BrSM_htoWW/sumModBRs;
    GTHDM_BR_h_ZZ = rh_VV*BrSM_htoZZ/sumModBRs;
    GTHDM_BR_h_tautau = BrSM_htotautau*(rh_QlQlE + rh_QlQlO/(beta(Mtau, m1_2)*beta(Mtau, m1_2)))/sumModBRs;
    GTHDM_BR_h_cc =(rh_QuQuE + rh_QuQuO/(beta(Mc, m1_2)*beta(Mc, m1_2)))*BrSM_htocc/sumModBRs;
    GTHDM_BR_h_gaga = rh_gaga*BrSM_htogaga/sumModBRs;
    GTHDM_BR_h_gg = (Gamma_hgg/Gamma_hggSM)*BrSM_htogg/sumModBRs;
    
}

double GeneralTHDMcache::computephi2quantities()
{
    
    m2_2 = mH2sq;
    m2 = sqrt(m2_2);
    m3_2 = mH3sq;
    m3 = sqrt(m3_2);
      
    double GF=1/(sqrt(2.0)*vev*vev);
    double sW2=1.0-cW2;

    //FLAG to select only the model in which all the couplings are the same (by families)
    if (!myGTHDM->getATHDMflag())
    {
        throw std::runtime_error("Direct Searches are only aviable in the A2HDM.");
    }
  
    /*complex i */
     gslpp::complex i = gslpp::complex::i();
 

        //fermionic couplings for phi2
    
    gslpp::complex yu2 = R21 + (R22 - i*R23)*su.conjugate();
    gslpp::complex yd2 = R21 + (R22 + i*R23)*sd;
    gslpp::complex yl2 = R21 + (R22 + i*R23)*sl;
    
    
    yu2R = R21_GTHDM + (R22_GTHDM)*su.real();
    yd2R = R21_GTHDM + (R22_GTHDM)*sd.real();
    yl2R = R21_GTHDM + (R22_GTHDM)*sl.real();
   
     

    //These cross sections ratios are necessary for rphi2_gg
    //At 8 TeV
    
    //SM gg -> H (phi even) production cross section ratio at 8 TeV, top loop only over total
    double rSigmaggphi2E_t8 = ip_csr_ggH_t_8(m2);
    //SM gg -> H (phi even) production cross section ratio at 8 TeV, bottom loop only over total
    double rSigmaggphi2E_b8 = ip_csr_ggH_b_8(m2);
    //gg -> H (phieven) production cross section at 8 TeV, total
   //  double Sigmaggphi2E_8 = ip_cs_ggtoH_8(m2);

    
    //gg -> A (phi odd) production cross section ratio at 8 TeV, top loop only over total
    double rSigmaggphi2O_t8 = ip_csr_ggA_t_8(m2);
    //gg -> A (phi odd) production cross section ratio at 8 TeV, bottom loop only over total
    double rSigmaggphi2O_b8 = ip_csr_ggA_b_8(m2);
    
    //gg -> A (phiodd) production cross section at 8 TeV, total
  //  double Sigmaggphi2O_8 = ip_cs_ggtoA_8(m2);
    
    
    /* r_ii is the ratio of the squared 2HDM vertex coupling of phi2
     * to the particle phi2 and the respective squared SM coupling.
     * phi2 is fixed to be the non-SM and lightests (phi2), but can be translated*/

    double rphi2_QuQuE= yu2.real()*yu2.real(); 
    double rphi2_QuQuO= yu2.imag()*yu2.imag(); 
    double rphi2_QdQdE= yd2.real()*yd2.real(); 
    double rphi2_QdQdO= yd2.imag()*yd2.imag(); 
    double rphi2_QlQlE= yl2.real()*yl2.real(); 
    double rphi2_QlQlO= yl2.imag()*yl2.imag(); 
    rphi2_ggE = yu2.real()*yd2.real() + (yu2.real()*yu2.real() - yu2.real()*yd2.real())*rSigmaggphi2E_t8  + (yd2.real()*yd2.real() - yu2.real()*yd2.real())*rSigmaggphi2E_b8;
    rphi2_ggO = yu2.imag()*yu2.imag() + (yu2.imag()*yu2.imag() - yu2.imag()*yd2.imag())*rSigmaggphi2O_t8  + (yd2.imag()*yd2.imag() - yu2.imag()*yd2.imag())*rSigmaggphi2O_b8;
    rphi2_VV=R21*R21;

    

    /*Gamma_phi2gaga and Gamma_phi2Zga expressions ...*/
    
    /*Decay to photons. The fermionic contribution has a CP-even part (HH) and a CP-odd (A)*/
    /*CP EVEN*/
   
                                                                                  
    gslpp::complex I_HH2_Ux=I_HH_U(m2_2,Mc,Mt);
    gslpp::complex I_HH2_Dx=I_HH_D(m2_2,Ms,Mb);
    gslpp::complex I_HH2_Lx=I_HH_L(m2_2,Mmu,Mtau);
    gslpp::complex I_phi2E_F= yu2.real()*I_HH2_Ux+ yd2.real()*I_HH2_Dx+yl2.real()*I_HH2_Lx;
                                                                               
    gslpp::complex I_phi2_W=R21*I_H_W(m2,MW);

    
    double lambdaphi2HpHm = lambdaipm(R21, R22, R32);                                 
    gslpp::complex I_phi2_Hp=(vev*vev)/(2.0*mHp2)*I_H_Hp(mHp2,m2)*(lambdaphi2HpHm);
    
    
    /*CP ODD */
            
    gslpp::complex I_A2_Ux=I_A_U(m2_2,Mc,Mt);
    gslpp::complex I_A2_Dx=I_A_D(m2_2,Ms,Mb);
    gslpp::complex I_A2_Lx=I_A_L(m2_2,Mmu,Mtau);
                                                                               
    gslpp::complex I_phi2O_F = yu2.imag()*I_A2_Ux + yd2.imag()*I_A2_Dx + yl2.imag()*I_A2_Lx;
                                                                             
    double Gamma_phi2gaga=(GF*Ale*Ale*m2*m2*m2/(sqrt(2.0)*128.0*M_PI*M_PI*M_PI))*((I_phi2E_F+I_phi2_W+I_phi2_Hp).abs2()
                        + (I_phi2O_F).abs2());
            
    /*Decay to Z gamma
    CP-EVEN PART*/

    gslpp::complex A_HH2_Ux = A_HH_U(m2_2,cW2,Mc,Mt,MZ);
    gslpp::complex A_HH2_Dx = A_HH_D(m2_2,cW2,Ms,Mb,MZ);
    gslpp::complex A_HH2_Lx = A_HH_L(m2_2,cW2,Mmu,Mtau,MZ);
                                                                               
    gslpp::complex A_phi2E_F = (yu2.real()*A_HH2_Ux+ yd2.real()*A_HH2_Dx+ yl2.real()*A_HH2_Lx)/sqrt(sW2*cW2);
    gslpp::complex A_phi2_W = R21*A_H_W(m2,cW2,MW,MZ);
  
    
    gslpp::complex A_phi2_Hp = (vev*vev)/(2.0*mHp2)*A_H_Hp(mHp2,m2,cW2,MZ)*(lambdaphi2HpHm);

    /*CP-ODD PART*/
                                                                               
    gslpp::complex A_A2_Ux = A_A_U(m2_2,cW2,Mc,Mt,MZ);
    gslpp::complex A_A2_Dx = A_A_D(m2_2,cW2,Ms,Mb,MZ);
    gslpp::complex A_A2_Lx = A_A_L(m2_2,cW2,Mmu,Mtau,MZ);
                                                                               
    gslpp::complex A_phi2O_F=yu2.imag()*A_A2_Ux + yd2.imag()*A_A2_Dx + yl2.imag()*A_A2_Lx;
    
                                                                             
                                                                               
    double Gamma_phi2Zga=HSTheta(m2-MZ)*GF*Ale*Ale*m2*m2*m2/(sqrt(2.0)*64.0*M_PI*M_PI*M_PI)*(1.0-MZ*MZ/(m2*m2))*(1.0-MZ*MZ/(m2*m2))*(1.0-MZ*MZ/(m2*m2))*((A_phi2E_F+A_phi2_W+A_phi2_Hp).abs2()
                        + A_phi2O_F.abs2());
       
    /*Decay to gluons*/
          
    
   double Gamma_phi2gg=(rphi2_ggE)*GF*Als*Als*m2*m2*m2/(sqrt(2.0)*16.0*M_PI*M_PI*M_PI)*(9.0/4.0)*(I_HH2_Ux/4.0+I_HH2_Dx).abs2()
                        +rphi2_ggO*GF*Als*Als*m2*m2*m2/(sqrt(2.0)*16.0*M_PI*M_PI*M_PI)*(9.0/4.0)*(I_A2_Ux/4.0+I_A2_Dx).abs2(); 
       
 
    //Cross-sections of ggF, bbF and VBF at 8 TeV Sigmaxx_H8 = Sigmaxx_H8SM*rphi2_xx
    /*
     SigmaggF_phi2_8=ip_cs_ggtoH_8(m2)*rphi2_gg;
     SigmabbF_phi2_8=ip_cs_pptobbH_8(m2)*rphi2_QbQb;
     SigmaVBF_phi2_8=ip_cs_VBFtoH_8(m2)*rphi2_VV;
     SigmattF_phi2_8=ip_cs_pptottH_8(m2)*rphi2_QtQt;
     SigmaVH_phi2_8=(ip_cs_WtoWH_8(m2)+ip_cs_ZtoZH_8(m2))*rphi2_VV;*/


    SigmaggF_phi2_8=ip_cs_ggtoH_8(m2)* rphi2_ggE + ip_cs_ggtoA_8(m2)* rphi2_ggO;
    SigmabbF_phi2_8=ip_cs_pptobbH_8(m2)*rphi2_QdQdE + ip_cs_pptobbA_8(m2)*rphi2_QdQdO;
    SigmaVBF_phi2_8=ip_cs_VBFtoH_8(m2)*rphi2_VV;
    SigmattF_phi2_8=ip_cs_pptottH_8(m2)*rphi2_QuQuE + ip_cs_pptottA_8(m2)*rphi2_QuQuO;
    SigmaVH_phi2_8=(ip_cs_WtoWH_8(m2)+ip_cs_ZtoZH_8(m2))*rphi2_VV;
        

    //SM PREDICTIONS

    SigmaTotSM_phi2_8 = 1.0e-15;

    if (m2>=20. && m2 <=2000.) {
        SigmaTotSM_phi2_8=ip_cs_ggtoH_8(m2)+ip_cs_VBFtoH_8(m2)+ip_cs_WtoWH_8(m2)+ip_cs_ZtoZH_8(m2)+ip_cs_pptottH_8(m2)+ip_cs_pptobbH_8(m2);
    }
     SigmaSumphi2_8 = SigmaggF_phi2_8 + SigmaVBF_phi2_8 + SigmaVH_phi2_8 + SigmattF_phi2_8 + SigmabbF_phi2_8;

     /*   SigmaggF_phi2_13=ip_cs_ggtoH_13(m2)*rphi2_gg;
        SigmabbF_phi2_13=ip_cs_pptobbH_13(m2)*rphi2_QbQb;
        SigmaVBF_phi2_13=ip_cs_VBFtoH_13(m2)*rphi2_VV;
        SigmattF_phi2_13=ip_cs_pptottH_13(m2)*rphi2_QtQt;
        SigmaVH_phi2_13=(ip_cs_WtoWH_13(m2)+ip_cs_ZtoZH_13(m2))*rphi2_VV;*/
 
 
  
    SigmaggF_phi2_13=ip_cs_ggtoH_13(m2)*rphi2_ggE + ip_cs_ggtoA_13(m2)*rphi2_ggO;
    SigmabbF_phi2_13=ip_cs_pptobbH_13(m2)*rphi2_QdQdE + ip_cs_pptobbA_13(m2)*rphi2_QdQdO;
    SigmaVBF_phi2_13=ip_cs_VBFtoH_13(m2)*rphi2_VV;
    SigmattF_phi2_13=ip_cs_pptottH_13(m2)*rphi2_QuQuE + ip_cs_pptottA_13(m2)*rphi2_QuQuO;
    SigmaVH_phi2_13=(ip_cs_WtoWH_13(m2)+ip_cs_ZtoZH_13(m2))*rphi2_VV;


 
//    double SigmaTotSM_H13 = 1.0e-15;
//    if (mHh>=20. && mHh <=2000.) {
//            SigmaTotSM_H13=ip_cs_ggtoH_13(mHh)+ip_cs_VBFtoH_13(mHh)+ip_cs_WtoWH_13(mHh)+ip_cs_ZtoZH_13(mHh)+ip_cs_pptottH_13(mHh)+ip_cs_pptobbH_13(mHh);
//    }
    SigmaSumphi2_13 = SigmaggF_phi2_13 + SigmaVBF_phi2_13 + SigmaVH_phi2_13 + SigmattF_phi2_13 + SigmabbF_phi2_13;
     
    double BrSM_phi2tott=ip_Br_HPtott(m2);
    double BrSM_phi2tocc=ip_Br_HPtocc(m2);
    double BrSM_phi2tobb=ip_Br_HPtobb(m2);
    double BrSM_phi2totautau=ip_Br_HPtotautau(m2);
   // double BrSM_phi2tomumu=ip_Br_HPtomumu(m2);
    double BrSM_phi2toWW =ip_Br_HPtoWW(m2);
    double BrSM_phi2toZZ =ip_Br_HPtoZZ(m2);

Gammaphi2totSM=ip_GammaHPtotSM(m2);
    
 /*Decay of phi3 to the others scalars*/
double lambda132 = lambdaijk(R11, R21, R31, R13, R23, R33, R12, R22, R32, lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H)  + lambdaijk(R12, R22, R32, R11, R22, R32, R12, R22, R33, lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H)   + lambdaijk(R12, R22, R32, R13, R23, R33, R11, R21, R33, lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H)   + lambdaijk(R11, R21, R31, R13, R23, R33, R12, R22, R32, lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H)  + lambdaijk(R13, R23,R33, R11, R21, R31, R12, R22, R32,lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H)   + lambdaijk(R13, R23,R33, R12, R22, R32, R11, R21, R31,lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H)  ;
double lambda332 = lambdaijk(R13, R23, R33, R13, R23, R33, R12, R22, R32, lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H)  + lambdaijk(R12, R22, R32, R13, R23, R33, R12, R22, R32, lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H)  + lambdaijk(R13, R23, R33, R12, R22, R32, R12, R22, R32, lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H) ;
//phi2 -> phi1 phi3
double Gammaphi2_phi1phi3=HSTheta(m2 - (m1+m3))*KaellenFunction(m2_2,m1_2,m3_2)*lambda132*lambda132/(8.0*m3_2*M_PI);
double Gammaphi2_phi3phi3=HSTheta(m2 - 2.0*m3)*sqrt(std::fabs(1.0 - (4.0*m3_2)/m2_2))*lambda332*lambda332/(32.0*m2*M_PI);

    
    /*Decay of phi2 to the others scalars*/
   
double lambda112 = (2.0)*(lambdaijk(R11, R21, R31, R12, R22, R32, R11, R21, R31, lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H) + lambdaijk(R11, R21, R31, R11, R21, R31, R12, R22, R32, lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H) + lambdaijk(R12, R22, R32, R11, R21, R31, R11, R21, R31,  lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H) );

   
//phi2 -> phi1phi1
double Gammaphi2_phi1phi1=HSTheta(m2 - 2.0*m1)*sqrt(std::fabs(1.0 - (4.0*m1_2)/m2_2))*lambda112*lambda112/(32.0*m2*M_PI);

 //phi2 ->H+H-
double Gammaphi2_HpHm=HSTheta(m2 - 2.0*sqrt(mHp2))*sqrt(std::fabs(1.0 - (4.0*mHp2)/m2_2))*lambdaphi2HpHm*lambdaphi2HpHm/(8.0*m2*M_PI);
//phi2 -> phi1 Z
double Gammaphi2_phi1Z=HSTheta(m2-(m1+MZ))*pow(KaellenFunction(m2_2,MZ*MZ,m1_2),3)*(R23*R12 + R22*R13)*(R23*R12 + R22*R13)/(2.0*M_PI*vev*vev);

//phi2 -> phi3 Z
double Gammaphi2_phi3Z=HSTheta(m2-(m3+MZ))*pow(KaellenFunction(m2_2,MZ*MZ,m3_2),3)*(R33*R22 + R32*R23)*(R33*R22 + R32*R23)/(2.0*M_PI*vev*vev);


/* phi2 -> H+W- */
double Gammaphi2_HpW=HSTheta(m2-sqrt(mHp2)-MW)*pow(KaellenFunction(m2_2,MW*MW,mHp2),3)*(R23-i*R22).abs2()/(M_PI*vev*vev);

Gammaphi2tot = 1.e-10;

Gammaphi2tot= Gammaphi2tot + (BrSM_phi2tott*(rphi2_QuQuE + rphi2_QuQuO/(beta_mt_sq(Mt, m2_2)))
         +BrSM_phi2tocc*(rphi2_QuQuE + rphi2_QuQuO/(beta(Mc, m2_2)*beta(Mc, m2_2)))
         +BrSM_phi2tobb*(rphi2_QdQdE + rphi2_QdQdO/(beta(Mb, m2_2)*beta(Mb, m2_2)))
         +BrSM_phi2totautau*(rphi2_QlQlE + rphi2_QlQlO/(beta(Mtau, m2_2)*beta(Mtau, m2_2)))
         +(BrSM_phi2toWW+BrSM_phi2toZZ)*rphi2_VV)*Gammaphi2totSM+Gamma_phi2gaga+Gamma_phi2Zga+Gamma_phi2gg 
         +Gammaphi2_phi1phi3+ Gammaphi2_phi1phi1+Br_phi2tophi3phi3+Gammaphi2_HpHm+Gammaphi2_phi1Z+Gammaphi2_phi3Z
         +Gammaphi2_HpW;
  

    Br_phi2tott=BrSM_phi2tott*(rphi2_QuQuE + rphi2_QuQuO/(beta_mt_sq(Mt, m2_2)))*Gammaphi2totSM/Gammaphi2tot;
    Br_phi2tobb=BrSM_phi2tobb*(rphi2_QdQdE + rphi2_QdQdO/(beta(Mb, m2_2)*beta(Mb, m2_2)))*Gammaphi2totSM/Gammaphi2tot;
    Br_phi2totautau=BrSM_phi2totautau*(rphi2_QlQlE + rphi2_QlQlO/(beta(Mtau, m2_2)*beta(Mtau, m2_2)))*Gammaphi2totSM/Gammaphi2tot;

    Br_phi2toWW=BrSM_phi2toWW*rphi2_VV*Gammaphi2totSM/Gammaphi2tot;
    Br_phi2toZZ=BrSM_phi2toZZ*rphi2_VV*Gammaphi2totSM/Gammaphi2tot;
    Br_phi2togaga=Gamma_phi2gaga/Gammaphi2tot;
    Br_phi2toZga=Gamma_phi2Zga/Gammaphi2tot;
    Br_phi2tophi1phi1=Gammaphi2_phi1phi1/Gammaphi2tot;
    Br_phi2tophi3phi3=Gammaphi2_phi3phi3/Gammaphi2tot;
    Br_phi2tophi1phi3=Gammaphi2_phi1phi3/Gammaphi2tot;
    Br_phi2toHpHm=Gammaphi2_HpHm/Gammaphi2tot;
    Br_phi2tophi1Z=Gammaphi2_phi1Z/Gammaphi2tot;
    Br_phi2tophi3Z=Gammaphi2_phi3Z/Gammaphi2tot;
    Br_phi2toHpW=Gammaphi2_HpW/Gammaphi2tot;
    
 return 0.;  
}


    
double GeneralTHDMcache::computephi3quantities()
{
    
  
    m2_2 = mH2sq;
    m2 = sqrt(m2_2);
    m3_2 = mH3sq;
    m3 = sqrt(m3_2);
    
    
    double GF=1/(sqrt(2.0)*vev*vev);
    double sW2=1.0-cW2;

      //FLAG to select only the model in which all the couplings are the same (by families)
    
    if (!myGTHDM->getATHDMflag())
    {
        throw std::runtime_error("Direct Searches are only aviable in the A2HDM.");
    }
  
        /*complex i */
    
     gslpp::complex i = gslpp::complex::i();

     //fermionic couplings for phi3
    
    gslpp::complex yu3 = R31 + (R32 - i*R33)*su.conjugate();
    gslpp::complex yd3 = R31 + (R32 + i*R33)*sd;
    gslpp::complex yl3 = R31 + (R32 + i*R33)*sl;
    
    yu3R = R31_GTHDM + (R32_GTHDM)*su.real();
    yd3R = R31_GTHDM + (R32_GTHDM)*sd.real();
    yl3R = R31_GTHDM + (R32_GTHDM)*sl.real();
   
    
    
    //These cross sections ratios are necessary for rphi3_gg
    //At 8 TeV
    
    //SM gg -> H (phi even) production cross section ratio at 8 TeV, top loop only over total
    double rSigmaggphi3E_t8 = ip_csr_ggH_t_8(m3);
    //SM gg -> H (phi even) production cross section ratio at 8 TeV, bottom loop only over total
    double rSigmaggphi3E_b8 = ip_csr_ggH_b_8(m3);
    //gg -> H (phieven) production cross section at 8 TeV, total
   // double Sigmaggphi3E_8 = ip_cs_ggtoH_8(m3);
    
    
    //gg -> A (phi odd) production cross section ratio at 8 TeV, top loop only over total
    double rSigmaggphi3O_t8 = ip_csr_ggA_t_8(m3);
    //gg -> A (phi odd) production cross section ratio at 8 TeV, bottom loop only over total
    double rSigmaggphi3O_b8 = ip_csr_ggA_b_8(m3);
    
    //gg -> A (phiodd) production cross section at 8 TeV, total
   // double Sigmaggphi3O_8 = ip_cs_ggtoA_8(m3);
    
   

       /* r_ii is the ratio of the squared 2HDM vertex coupling of phi3
     * to the particle phi3 and the respective squared SM coupling.
     * phi is fixed to be the non-SM and heaviest (phi3), but can be translated*/

   
    double rphi3_QuQuE= yu3.real()*yu3.real();   
    double rphi3_QuQuO= yu3.imag()*yu3.imag(); 
    double rphi3_QdQdE= yd3.real()*yd3.real();  
    double rphi3_QdQdO= yd3.imag()*yd3.imag(); 
    double rphi3_QlQlE= yl3.real()*yl3.real(); 
    double rphi3_QlQlO= yl3.imag()*yl3.imag(); 
     rphi3_ggE = yu3.real()*yd3.real() + (yu3.real()*yu3.real() - yu3.real()*yd3.real())*rSigmaggphi3E_t8  + (yd3.real()*yd3.real() - yu3.real()*yd3.real())*rSigmaggphi3E_b8;
     rphi3_ggO = yu3.imag()*yu3.imag() + (yu3.imag()*yu3.imag() - yu3.imag()*yd3.imag())*rSigmaggphi3O_t8  + (yd3.imag()*yd3.imag() - yu3.imag()*yd3.imag())*rSigmaggphi3O_b8;
     rphi3_VV=R31*R31;
     

    /*Gamma_phi3gaga and Gamma_phi3Zga expressions ...*/
    
    /*Decay to photons. The fermionic contribution has a CP-even part (HH) and a CP-odd (A)*/
    /*CP EVEN*/
                                                                                   
    gslpp::complex I_HH_Ux=I_HH_U(m3_2,Mc,Mt);
    gslpp::complex I_HH_Dx=I_HH_D(m3_2,Ms,Mb);
    gslpp::complex I_HH_Lx=I_HH_L(m3_2,Mmu,Mtau);
    gslpp::complex I_phi3E_F= yu3.real()*I_HH_Ux+ yd3.real()*I_HH_Dx+yl3.real()*I_HH_Lx;
                                                                               
    gslpp::complex I_phi3_W=R31*I_H_W(m3,MW);


                          
    double lambdaphi3HpHm =  lambdaipm(R31, R32, R33);
    gslpp::complex I_phi3_Hp=(vev*vev)/(2.0*mHp2)*I_H_Hp(mHp2,m3)*(lambdaphi3HpHm);
    
    
    /*CP ODD */
            
    gslpp::complex I_A_Ux=I_A_U(m3_2,Mc,Mt);
    gslpp::complex I_A_Dx=I_A_D(m3_2,Ms,Mb);
    gslpp::complex I_A_Lx=I_A_L(m3_2,Mmu,Mtau);
                                                                               
    gslpp::complex I_phi3O_F = yu3.imag()*I_A_Ux + yd3.imag()*I_A_Dx + yl3.imag()*I_A_Lx;
                                                                             
    double Gamma_phi3gaga=(GF*Ale*Ale*m3*m3*m3/(sqrt(2.0)*128.0*M_PI*M_PI*M_PI))*((I_phi3E_F+I_phi3_W+I_phi3_Hp).abs2()
    + (I_phi3O_F).abs2());
    
    
                                                                              
    /*Decay to Z gamma
    CP-EVEN PART*/

    gslpp::complex A_HH_Ux = A_HH_U(m3_2,cW2,Mc,Mt,MZ);
    gslpp::complex A_HH_Dx = A_HH_D(m3_2,cW2,Ms,Mb,MZ);
    gslpp::complex A_HH_Lx = A_HH_L(m3_2,cW2,Mmu,Mtau,MZ);
                                                                               
    gslpp::complex A_phi3E_F = (yu3.real()*A_HH_Ux+ yd3.real()*A_HH_Dx+ yl3.real()*A_HH_Lx)/sqrt(sW2*cW2);
    gslpp::complex A_phi3_W = R31*A_H_W(m3,cW2,MW,MZ);
  
    gslpp::complex A_phi3_Hp = (vev*vev)/(2.0*mHp2)*A_H_Hp(mHp2,m3,cW2,MZ)*(lambdaphi3HpHm);

    /*CP-ODD PART*/
                                                                               
    gslpp::complex A_A_Ux = A_A_U(m3_2,cW2,Mc,Mt,MZ);
    gslpp::complex A_A_Dx = A_A_D(m3_2,cW2,Ms,Mb,MZ);
    gslpp::complex A_A_Lx = A_A_L(m3_2,cW2,Mmu,Mtau,MZ);
                                                                               
    gslpp::complex A_phi3O_F=yu3.imag()*A_A_Ux + yd3.imag()*A_A_Dx + yl3.imag()*A_A_Lx;
                                                                                 
    double Gamma_phi3Zga=HSTheta(m3-MZ)*GF*Ale*Ale*m3*m3*m3/(sqrt(2.0)*64.0*M_PI*M_PI*M_PI)*(1.0-MZ*MZ/(m3*m3))*(1.0-MZ*MZ/(m3*m3))*(1.0-MZ*MZ/(m3*m3))*((A_phi3E_F+A_phi3_W+A_phi3_Hp).abs2()
                        + A_phi3O_F.abs2());

    
    /*Decay to gluons */
                                                                               
   double Gamma_phi3gg=(rphi3_ggE)*GF*Als*Als*m3*m3*m3/(sqrt(2.0)*16.0*M_PI*M_PI*M_PI)*(9.0/4.0)*(I_HH_Ux/4.0+I_HH_Dx).abs2()
                        +rphi3_ggO*GF*Als*Als*m3*m3*m3/(sqrt(2.0)*16.0*M_PI*M_PI*M_PI)*(9.0/4.0)*(I_A_Ux/4.0+I_A_Dx).abs2();
    

    //Cross-sections of ggF, bbF and VBF at 8 TeV Sigmaxx_H8 = Sigmaxx_H8SM*rphi3_xx

     /*SigmaggF_phi3_8=ip_cs_ggtoH_8(m3)*rphi3_gg;
     SigmabbF_phi3_8=ip_cs_pptobbH_8(m3)*rphi3_QbQb;
     SigmattF_phi3_8=ip_cs_pptottH_8(m3)*rphi3_QtQt;
    */

    SigmaggF_phi3_8=ip_cs_ggtoH_8(m3)*rphi3_ggE + ip_cs_ggtoA_8(m3)*rphi3_ggO;
    SigmabbF_phi3_8=ip_cs_pptobbH_8(m3)*rphi3_QdQdE + ip_cs_pptobbA_8(m3)*rphi3_QdQdO;
    SigmaVBF_phi3_8=ip_cs_VBFtoH_8(m3)*rphi3_VV;
    SigmattF_phi3_8=ip_cs_pptottH_8(m3)*rphi3_QuQuE + ip_cs_pptottA_8(m3)*rphi3_QuQuO;
    SigmaVH_phi3_8=(ip_cs_WtoWH_8(m3)+ip_cs_ZtoZH_8(m3))*rphi3_VV;
      
    //SM PREDICTIONS
     SigmaTotSM_phi3_8 = 1.0e-15;

    if (m3>=20. && m3 <=2000.) {
        SigmaTotSM_phi3_8=ip_cs_ggtoH_8(m3)+ip_cs_VBFtoH_8(m3)+ip_cs_WtoWH_8(m3)+ip_cs_ZtoZH_8(m3)+ip_cs_pptottH_8(m3)+ip_cs_pptobbH_8(m3);
    }
     SigmaSumphi3_8 = SigmaggF_phi3_8 + SigmaVBF_phi3_8 + SigmaVH_phi3_8 + SigmattF_phi3_8 + SigmabbF_phi3_8;

  /*  SigmaggF_phi3_13=ip_cs_ggtoH_13(m3)*rphi3_gg;
    SigmabbF_phi3_13=ip_cs_pptobbH_13(m3)*rphi3_QbQb;
    SigmattF_phi3_13=ip_cs_pptottH_13(m3)*rphi3_QtQt;*/
    
    SigmaggF_phi3_13=ip_cs_ggtoH_13(m3)* rphi3_ggE + ip_cs_ggtoA_13(m3)* rphi3_ggO;
    SigmabbF_phi3_13=ip_cs_pptobbH_13(m3)*rphi3_QdQdE + ip_cs_pptobbA_13(m3)*rphi3_QdQdO;
    SigmaVBF_phi3_13=ip_cs_VBFtoH_13(m3)*rphi3_VV;
    SigmattF_phi3_13=ip_cs_pptottH_13(m3)*rphi3_QuQuE + ip_cs_pptottA_13(m3)*rphi3_QuQuO;
    SigmaVH_phi3_13=(ip_cs_WtoWH_13(m3)+ip_cs_ZtoZH_13(m3))*rphi3_VV;

  
//    double SigmaTotSM_H13 = 1.0e-15;
//    if (mHh>=20. && mHh <=2000.) {
//            SigmaTotSM_H13=ip_cs_ggtoH_13(mHh)+ip_cs_VBFtoH_13(mHh)+ip_cs_WtoWH_13(mHh)+ip_cs_ZtoZH_13(mHh)+ip_cs_pptottH_13(mHh)+ip_cs_pptobbH_13(mHh);
//    }
    SigmaSumphi3_13 = SigmaggF_phi3_13 + SigmaVBF_phi3_13 + SigmaVH_phi3_13 + SigmattF_phi3_13 + SigmabbF_phi3_13;

    
double BrSM_phi3tott=ip_Br_HPtott(m3);
double BrSM_phi3tocc=ip_Br_HPtocc(m3);
double BrSM_phi3tobb=ip_Br_HPtobb(m3);
double BrSM_phi3totautau=ip_Br_HPtotautau(m3);
//double BrSM_phi3tomumu=ip_Br_HPtomumu(m3);
double BrSM_phi3toWW =ip_Br_HPtoWW(m3);
double BrSM_phi3toZZ =ip_Br_HPtoZZ(m3);

Gammaphi3totSM=ip_GammaHPtotSM(m3);

 
 
/*Decay of phi3 to the others scalars*/
double lambda123 = lambdaijk(R11, R21, R13, R12, R22, R32, R13, R23, R33, lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H)  + lambdaijk(R12, R22, R32, R11, R22, R32, R12, R22, R33, lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H)   + lambdaijk(R12, R22, R32, R13, R23, R33, R11, R21, R33, lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H)   + lambdaijk(R11, R21, R31, R13, R23, R33, R12, R22, R32, lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H)  + lambdaijk(R13, R23,R33, R11, R21, R31, R12, R22, R32,lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H)   + lambdaijk(R13, R23,R33, R12, R22, R32, R11, R21, R31,lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H)  ;
double lambda223 = lambdaijk(R11, R22, R32, R12, R22, R32, R13, R23, R33, lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H)  + lambdaijk(R12, R22, R32, R13, R23, R33, R12, R22, R32, lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H)  + lambdaijk(R13, R23, R33, R12, R22, R32, R12, R22, R32, lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H) ;
double lambda113 = lambdaijk(R11, R21, R31, R11, R21, R31, R13, R23, R33, lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H)  + lambdaijk(R11, R21, R31, R13, R23, R33, R11, R21, R31,lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H)  + lambdaijk(R13, R23, R33, R11, R21, R31, R11, R21, R31,lambda1H, lambda3H, lambda4H, Relambda5H, Imlambda5H, Relambda6H, Imlambda6H, Relambda7H, Imlambda7H) ;
//phi3 -> phi1 phi2
double Gammaphi3_phi1phi2=HSTheta(m3 - (m1+m2))*KaellenFunction(m3_2,m1_2,m2_2)*lambda123*lambda123/(8.0*m3_2*M_PI);
// phi3 -> phi2 phi2
double Gammaphi3_phi2phi2=HSTheta(m3 - 2.0*m2)*sqrt(std::fabs(1.0 - (4.0*m2_2)/m3_2))*lambda223*lambda223/(32.0*m3*M_PI);
//phi3 -> phi1phi1
double Gammaphi3_phi1phi1=HSTheta(m3 - 2.0*m1)*sqrt(std::fabs(1.0 - (4.0*m1_2)/m3_2))*lambda113*lambda113/(32.0*m3*M_PI);
//phi3 ->H+H-
double Gammaphi3_HpHm=HSTheta(m3 - 2.0*sqrt(mHp2))*sqrt(std::fabs(1.0 - (4.0*mHp2)/m3_2))*lambdaphi3HpHm*lambdaphi3HpHm/(32.0*m3*M_PI);
//phi3 -> phi1 Z
double Gammaphi3_phi1Z=HSTheta(m3-(m1+MZ))*pow(KaellenFunction(m3_2,MZ*MZ,m1_2),3)*(R33*R12 + R32*R13)*(R33*R12 + R32*R13)/(2.0*M_PI*vev*vev);
//phi3 -> phi2 Z
double Gammaphi3_phi2Z=HSTheta(m3-(m2+MZ))*pow(KaellenFunction(m3_2,MZ*MZ,m2_2),3)*(R33*R22 + R32*R23)*(R33*R22 + R32*R23)/(2.0*M_PI*vev*vev);
/* phi3 -> H+W- */
double Gammaphi3_HpW=HSTheta(m3-sqrt(mHp2)-MW)*pow(KaellenFunction(m3_2,MW*MW,mHp2),3)*(R33-i*R32).abs2()/(M_PI*vev*vev);
                                                                              
Gammaphi3tot  = 1.e-10;

 Gammaphi3tot= Gammaphi3tot + (BrSM_phi3tott*(rphi3_QuQuE + rphi3_QuQuO/(beta_mt_sq(Mt, m3_2)))
         +BrSM_phi3tocc*(rphi3_QuQuE + rphi3_QuQuO/(beta(Mc, m3_2)*beta(Mc, m3_2)))
         +BrSM_phi3tobb*(rphi3_QdQdE + rphi3_QdQdO/(beta(Mb, m3_2)*beta(Mb, m3_2)))
         +BrSM_phi3totautau*(rphi3_QlQlE + rphi3_QlQlO/(beta(Mtau, m3_2)*beta(Mtau, m3_2)))
         +(BrSM_phi3toWW+BrSM_phi3toZZ)*rphi3_VV)*Gammaphi3totSM+Gamma_phi3gaga+Gamma_phi3Zga+Gamma_phi3gg 
         + Gammaphi3_phi1phi1+Gammaphi3_phi2phi2+Gammaphi3_phi1phi2+Gammaphi3_HpHm+Gammaphi3_phi1Z+Gammaphi3_phi2Z
         +Gammaphi3_HpW;

    Br_phi3tott=BrSM_phi3tott*(rphi3_QuQuE + rphi3_QuQuO/(beta_mt_sq(Mt, m3_2)))*Gammaphi3totSM/Gammaphi3tot;
    Br_phi3tobb=BrSM_phi3tobb*(rphi3_QdQdE + rphi3_QdQdO/(beta(Mb, m3_2)*beta(Mb, m3_2)))*Gammaphi3totSM/Gammaphi3tot;
    
   
    Br_phi3totautau=BrSM_phi3totautau*(rphi3_QlQlE + rphi3_QlQlO/(beta(Mtau, m3_2)*beta(Mtau, m3_2)))*Gammaphi3totSM/Gammaphi3tot;
    Br_phi3toWW=BrSM_phi3toWW*rphi3_VV*Gammaphi3totSM/Gammaphi3tot;
    Br_phi3toZZ=BrSM_phi3toZZ*rphi3_VV*Gammaphi3totSM/Gammaphi3tot;
    Br_phi3togaga=Gamma_phi3gaga/Gammaphi3tot;
    Br_phi3toZga=Gamma_phi3Zga/Gammaphi3tot;
    Br_phi3tophi1phi1=Gammaphi3_phi1phi1/Gammaphi3tot;
    Br_phi3tophi2phi2=Gammaphi3_phi2phi2/Gammaphi3tot;
    
    Br_phi3tophi1phi2=Gammaphi3_phi1phi2/Gammaphi3tot;
    Br_phi3toHpHm=Gammaphi3_HpHm/Gammaphi3tot;
    Br_phi3tophi1Z=Gammaphi3_phi1Z/Gammaphi3tot;
    Br_phi3tophi2Z=Gammaphi3_phi2Z/Gammaphi3tot;
    Br_phi3toHpW=Gammaphi3_HpW/Gammaphi3tot;
 
   
 return 0.;


}

    
double GeneralTHDMcache::computeHpquantities()
{
    
    m2_2 = mH2sq;
    m2 = sqrt(m2_2);
    m3_2 = mH3sq;
    m3 = sqrt(m3_2);
    
//    double GF=1/(sqrt(2.0)*vev*vev);
//    double sW2=1.0-cW2;
    double Mb2 = Mb*Mb;
    double Mt2 = Mt*Mt;
    double Mtau2=Mtau*Mtau;
    double MW2 = MW*MW;
    double Vtb=myGTHDM->getCKM().getV_tb().abs();
    
    m2 = sqrt(m2_2);
    m3 = sqrt(m3_2);

      //FLAG to select only the model in which all the couplings are the same (by families)
    
    if (!myGTHDM->getATHDMflag())
    {
        throw std::runtime_error("Direct Searches are only aviable in the A2HDM.");
    }
  
        /*complex i */
    
     gslpp::complex i = gslpp::complex::i();
     
    SigmaHp8=0.0;
    SigmaHpm13=0.0;
     
    double GammaHptaunu=HSTheta(mHp-Mtau)*(Mtau2*(mHp2-Mtau2)*(mHp2-Mtau2)*sl.abs2())/(8.0*mHp*mHp2*M_PI*vev*vev);
    
    
    
    double GammaHptb=HSTheta(mHp-Mt-Mb)*(Vtb*Vtb/(8.0*mHp*M_PI*vev*vev))*3.0*(-4.0*(su*sd).real()*Mb2*Mt2
                        -sd.abs2()*Mb2*(Mb2-mHp2+Mt2)-su.abs2()*Mt2*(Mb2-mHp2+Mt2))
                      *sqrt((Mb2*Mb2+(mHp2-Mt2)*(mHp2-Mt2)-2.0*Mb2*(mHp2+Mt2))/(mHp2*mHp2));
    
   
    //decay to SM Higgs
    double GammaHpHlW=KaellenFunction(1.0,mHl/mHp,MW/mHp)*KaellenFunction(1.0,mHp/MW,mHl/MW)*KaellenFunction(1.0,mHp/MW,mHl/MW)
                      *MW2*MW2/mHp*(R21-R31)*(R21-R31)/(2.0*M_PI*vev*vev);
    
    //decay to phi2
    double GammaHpphi2W=KaellenFunction(1.0,m2/mHp,MW/mHp)*KaellenFunction(1.0,mHp/MW,m2/MW)*KaellenFunction(1.0,mHp/MW,m2/MW)
                      *MW2*MW2/mHp*(R22-R32)*(R22-R32)/(2.0*M_PI*vev*vev);
    
    //decay to phi3
    double GammaHpphi3W=KaellenFunction(1.0,m3/mHp,MW/mHp)*KaellenFunction(1.0,mHp/MW,m3/MW)*KaellenFunction(1.0,mHp/MW,m3/MW)
                      *MW2*MW2/mHp*(R23-R33)*(R23-R33)/(2.0*M_PI*vev*vev);
    
    GammaHptot = 1.e-10;
    
    GammaHptot= GammaHptot + GammaHptaunu + GammaHptb + GammaHpHlW + GammaHpphi2W + GammaHpphi3W;


    Br_Hptotaunu=GammaHptaunu/GammaHptot;
    Br_Hptotb=GammaHptb/GammaHptot;

  
     return 0;
}

//Higgs direct searches

double GeneralTHDMcache::ComputeHeavyHiggs()
{
    
    m2_2 = mH2sq;
    m2 = sqrt(m2_2);
    m3_2 = mH3sq;
    m3 = sqrt(m3_2);
    
    
    //FLAG to select only the model in which all the couplings are the same (by families)
    if (!myGTHDM->getATHDMflag())
    {
        throw std::runtime_error("Direct Searches are only aviable in the A2HDM.");
    }
  
    /*complex i */
    gslpp::complex i = gslpp::complex::i();

    double Br_Ztoee=0.03363; //C. Patrignani et al.(Particle Data Group), Chin. Phys. C, 40, 100001 (2016)
    double Br_Ztomumu=0.03366; //C. Patrignani et al.(Particle Data Group), Chin. Phys. C, 40, 100001 (2016)
    double Br_Ztotautau=0.0337; //C. Patrignani et al.(Particle Data Group), Chin. Phys. C, 40, 100001 (2016)
    double Br_Ztoinv=0.2; //C. Patrignani et al.(Particle Data Group), Chin. Phys. C, 40, 100001 (2016)
    double Br_Wtoenu=0.1071; //C. Patrignani et al.(Particle Data Group), Chin. Phys. C, 40, 100001 (2016)
    double Br_Wtomunu=0.1063; //C. Patrignani et al.(Particle Data Group), Chin. Phys. C, 40, 100001 (2016)
    double Br_Wtotaunu=0.1138; //C. Patrignani et al.(Particle Data Group), Chin. Phys. C, 40, 100001 (2016)

    THoEX_bb_phi2_tt_ATLAS13=0.0;
    THoEX_bb_phi3_tt_ATLAS13=0.0;
    THoEX_bb_phi2_bb_CMS8=0.0;
    THoEX_bb_phi3_bb_CMS8=0.0;
    THoEX_gg_phi2_bb_CMS8=0.0;
    THoEX_gg_phi3_bb_CMS8=0.0;
    THoEX_pp_phi2_bb_CMS13=0.0;
    THoEX_pp_phi3_bb_CMS13=0.0;
    THoEX_bb_phi2_bb_CMS13=0.0;
    THoEX_bb_phi3_bb_CMS13=0.0;
    THoEX_gg_phi2_tautau_ATLAS8=0.0;
    THoEX_gg_phi3_tautau_ATLAS8=0.0;
    THoEX_gg_phi2_tautau_CMS8=0.0;
    THoEX_gg_phi3_tautau_CMS8=0.0;
    THoEX_bb_phi2_tautau_ATLAS8=0.0;
    THoEX_bb_phi3_tautau_ATLAS8=0.0;
    THoEX_bb_phi2_tautau_CMS8=0.0;
    THoEX_bb_phi3_tautau_CMS8=0.0;
    THoEX_gg_phi2_tautau_ATLAS13=0.0;
    THoEX_gg_phi3_tautau_ATLAS13=0.0;
    THoEX_gg_phi2_tautau_CMS13=0.0;
    THoEX_gg_phi3_tautau_CMS13=0.0;
    THoEX_bb_phi2_tautau_ATLAS13=0.0;
    THoEX_bb_phi3_tautau_ATLAS13=0.0;
    THoEX_bb_phi2_tautau_CMS13=0.0;
    THoEX_bb_phi3_tautau_CMS13=0.0;
    THoEX_gg_phi2_gaga_ATLAS8=0.0;
    THoEX_gg_phi3_gaga_ATLAS8=0.0;
    THoEX_pp_phi2_gaga_ATLAS13=0.0;
    THoEX_pp_phi3_gaga_ATLAS13=0.0;
    THoEX_gg_phi2_gaga_CMS13=0.0;
    THoEX_gg_phi3_gaga_CMS13=0.0;
    THoEX_pp_phi2_Zga_llga_ATLAS8=0.0;
    THoEX_pp_phi3_Zga_llga_ATLAS8=0.0;
    THoEX_pp_phi2_Zga_llga_CMS8=0.0;
    THoEX_pp_phi3_Zga_llga_CMS8=0.0;
    THoEX_gg_phi2_Zga_llga_ATLAS13=0.0;
    THoEX_gg_phi3_Zga_llga_ATLAS13=0.0;
    THoEX_gg_phi2_Zga_qqga_ATLAS13=0.0;
    THoEX_gg_phi3_Zga_qqga_ATLAS13=0.0;
    THoEX_gg_phi2_Zga_CMS13=0.0;
    THoEX_gg_phi3_Zga_CMS13=0.0;
    THoEX_gg_phi2_ZZ_ATLAS8=0.0;
    THoEX_gg_phi3_ZZ_ATLAS8=0.0;
    THoEX_VV_phi2_ZZ_ATLAS8=0.0;
    THoEX_VV_phi3_ZZ_ATLAS8=0.0;
    THoEX_gg_phi2_ZZ_llllnunu_ATLAS13=0.0;
    THoEX_gg_phi3_ZZ_llllnunu_ATLAS13=0.0;
    THoEX_VV_phi2_ZZ_llllnunu_ATLAS13=0.0;
    THoEX_VV_phi3_ZZ_llllnunu_ATLAS13=0.0;
    THoEX_gg_phi2_ZZ_qqllnunu_ATLAS13=0.0;
    THoEX_gg_phi3_ZZ_qqllnunu_ATLAS13=0.0;
    THoEX_VV_phi2_ZZ_qqllnunu_ATLAS13=0.0;
    THoEX_VV_phi3_ZZ_qqllnunu_ATLAS13=0.0;
    THoEX_pp_phi2_ZZ_llqqnunull_CMS13=0.0;
    THoEX_pp_phi3_ZZ_llqqnunull_CMS13=0.0;
    THoEX_pp_phi2_ZZ_qqnunu_CMS13=0.0;
    THoEX_pp_phi3_ZZ_qqnunu_CMS13=0.0;
    THoEX_gg_phi2_WW_ATLAS8=0.0;
    THoEX_gg_phi3_WW_ATLAS8=0.0;
    THoEX_VV_phi2_WW_ATLAS8=0.0;
    THoEX_VV_phi3_WW_ATLAS8=0.0;
    THoEX_gg_phi2_WW_enumunu_ATLAS13=0.0;
    THoEX_gg_phi3_WW_enumunu_ATLAS13=0.0;
    THoEX_VV_phi2_WW_enumunu_ATLAS13=0.0;
    THoEX_VV_phi3_WW_enumunu_ATLAS13=0.0;
    THoEX_ggVV_phi2_WW_lnulnu_CMS13=0.0;
    THoEX_ggVV_phi3_WW_lnulnu_CMS13=0.0;
    THoEX_gg_phi2_WW_lnuqq_ATLAS13=0.0;
    THoEX_gg_phi3_WW_lnuqq_ATLAS13=0.0;
    THoEX_VV_phi2_WW_lnuqq_ATLAS13=0.0;
    THoEX_VV_phi3_WW_lnuqq_ATLAS13=0.0;
    THoEX_pp_phi2_WW_lnuqq_CMS13=0.0;
    THoEX_pp_phi3_WW_lnuqq_CMS13=0.0;
    THoEX_mu_pp_phi2_VV_CMS8=0.0;
    THoEX_mu_pp_phi3_VV_CMS8=0.0;
    THoEX_pp_phi2_VV_qqqq_ATLAS13=0.0;
    THoEX_pp_phi3_VV_qqqq_ATLAS13=0.0;
    THoEX_gg_phi2_phi1phi1_ATLAS8=0.0;
    THoEX_gg_phi3_phi1phi1_ATLAS8=0.0;
    THoEX_pp_phi2_phi1phi1_bbbb_CMS8=0.0;
    THoEX_pp_phi3_phi1phi1_bbbb_CMS8=0.0;
    THoEX_pp_phi2_phi1phi1_bbgaga_CMS8=0.0;
    THoEX_pp_phi3_phi1phi1_bbgaga_CMS8=0.0;
    THoEX_gg_phi2_phi1phi1_bbtautau_CMS8=0.0;
    THoEX_gg_phi3_phi1phi1_bbtautau_CMS8=0.0;
    THoEX_pp_phi2_phi1phi1_bbtautau_CMS8=0.0;
    THoEX_pp_phi3_phi1phi1_bbtautau_CMS8=0.0;
    THoEX_pp_phi2_phi1phi1_bbbb_ATLAS13=0.0;
    THoEX_pp_phi3_phi1phi1_bbbb_ATLAS13=0.0;
    THoEX_pp_phi2_phi1phi1_bbbb_1_CMS13=0.0;
    THoEX_pp_phi3_phi1phi1_bbbb_1_CMS13=0.0;
    THoEX_pp_phi2_phi1phi1_bbbb_2_CMS13=0.0;
    THoEX_pp_phi3_phi1phi1_bbbb_2_CMS13=0.0;
    THoEX_pp_phi2_phi1phi1_bbgaga_ATLAS13=0.0;
    THoEX_pp_phi3_phi1phi1_bbgaga_ATLAS13=0.0;
    THoEX_pp_phi2_phi1phi1_bbgaga_CMS13=0.0;
    THoEX_pp_phi3_phi1phi1_bbgaga_CMS13=0.0;
    THoEX_pp_phi2_phi1phi1_bbtautau_ATLAS13=0.0;
    THoEX_pp_phi3_phi1phi1_bbtautau_ATLAS13=0.0;
    THoEX_pp_phi2_phi1phi1_bbtautau_1_CMS13=0.0;
    THoEX_pp_phi3_phi1phi1_bbtautau_1_CMS13=0.0;
    THoEX_pp_phi2_phi1phi1_bbtautau_2_CMS13=0.0;
    THoEX_pp_phi3_phi1phi1_bbtautau_2_CMS13=0.0;
    THoEX_pp_phi2_phi1phi1_bbVV_CMS13=0.0;
    THoEX_pp_phi3_phi1phi1_bbVV_CMS13=0.0;
    THoEX_pp_phi2_phi1phi1_bbWW_ATLAS13=0.0;
    THoEX_pp_phi3_phi1phi1_bbWW_ATLAS13=0.0;       
    THoEX_gg_phi2_phi1phi1_gagaWW_ATLAS13=0.0;
    THoEX_gg_phi3_phi1phi1_gagaWW_ATLAS13=0.0;
    THoEX_gg_phi2_phi1Z_bbZ_ATLAS8=0.0;
    THoEX_gg_phi3_phi1Z_bbZ_ATLAS8=0.0;
    THoEX_gg_phi2_phi1Z_bbll_CMS8=0.0;
    THoEX_gg_phi3_phi1Z_bbll_CMS8=0.0;
    THoEX_gg_phi2_phi1Z_tautauZ_ATLAS8=0.0;
    THoEX_gg_phi3_phi1Z_tautauZ_ATLAS8=0.0;
    THoEX_gg_phi2_phi1Z_tautaull_CMS8=0.0;
    THoEX_gg_phi3_phi1Z_tautaull_CMS8=0.0;
    THoEX_gg_phi2_phi1Z_bbZ_ATLAS13=0.0;
    THoEX_gg_phi3_phi1Z_bbZ_ATLAS13=0.0;
    THoEX_gg_phi2_phi1Z_bbZ_1_CMS13=0.0;
    THoEX_gg_phi3_phi1Z_bbZ_1_CMS13=0.0;
    THoEX_gg_phi2_phi1Z_bbZ_2_CMS13=0.0;
    THoEX_gg_phi3_phi1Z_bbZ_2_CMS13=0.0;
    THoEX_bb_phi2_phi1Z_bbZ_ATLAS13=0.0;
    THoEX_bb_phi3_phi1Z_bbZ_ATLAS13=0.0;
    THoEX_bb_phi2_phi1Z_bbZ_1_CMS13=0.0;
    THoEX_bb_phi3_phi1Z_bbZ_1_CMS13=0.0;
    THoEX_bb_phi2_phi1Z_bbZ_2_CMS13=0.0;
    THoEX_bb_phi3_phi1Z_bbZ_2_CMS13=0.0;

    THoEX_pp_phi3_phi2Z_bbll_1_CMS8=0.0;
    THoEX_pp_phi2_phi3Z_bbll_1_CMS8=0.0;   
    THoEX_pp_phi3_phi2Z_bbll_2_CMS8=0.0;
    THoEX_pp_phi2_phi3Z_bbll_2_CMS8=0.0;    
    THoEX_pp_phi3_phi2Z_tautaull_1_CMS8=0.0;
    THoEX_pp_phi2_phi3Z_tautaull_1_CMS8=0.0;    
    THoEX_pp_phi3_phi2Z_tautaull_2_CMS8=0.0;
    THoEX_pp_phi2_phi3Z_tautaull_2_CMS8=0.0;    
    THoEX_gg_phi3_phi2Z_bbZ_ATLAS13=0.0;
    THoEX_gg_phi2_phi3Z_bbZ_ATLAS13=0.0;  
    THoEX_bb_phi3_phi2Z_bbZ_ATLAS13=0.0;
    THoEX_bb_phi2_phi3Z_bbZ_ATLAS13=0.0;
    
    THoEX_pp_Hpm_taunu_ATLAS8=0.0;
    THoEX_pp_Hp_taunu_CMS8=0.0;
    THoEX_pp_Hpm_taunu_ATLAS13=0.0;
    THoEX_pp_Hpm_taunu_CMS13=0.0;
    THoEX_pp_Hpm_tb_ATLAS8=0.0;
    THoEX_pp_Hp_tb_CMS8=0.0;
    THoEX_pp_Hpm_tb_ATLAS13=0.0;

   //Theoretical expressions for the Heavy Higgs cross sections times branching ratios

    tt_phi2_tt_TH13=SigmattF_phi2_13*Br_phi2tott;
    tt_phi3_tt_TH13=SigmattF_phi3_13*Br_phi3tott;
    bb_phi2_tt_TH13=SigmabbF_phi2_13*Br_phi2tott;
    bb_phi3_tt_TH13=SigmabbF_phi3_13*Br_phi3tott;

    bb_phi2_bb_TH8=SigmabbF_phi2_8*Br_phi2tobb;
    bb_phi3_bb_TH8=SigmabbF_phi3_8*Br_phi3tobb;
    gg_phi2_bb_TH8=SigmaggF_phi2_8*Br_phi2tobb;
    gg_phi3_bb_TH8=SigmaggF_phi3_8*Br_phi3tobb;
    pp_phi2_bb_TH13=SigmaSumphi2_13*Br_phi2tobb;
    pp_phi3_bb_TH13=SigmaSumphi3_13*Br_phi3tobb;
    bb_phi2_bb_TH13=SigmabbF_phi2_13*Br_phi2tobb;
    bb_phi3_bb_TH13=SigmabbF_phi3_13*Br_phi3tobb;

    gg_phi2_tautau_TH8=SigmaggF_phi2_8*Br_phi2totautau;
    gg_phi3_tautau_TH8=SigmaggF_phi3_8*Br_phi3totautau;
    bb_phi2_tautau_TH8=SigmabbF_phi2_8*Br_phi2totautau;
    bb_phi3_tautau_TH8=SigmabbF_phi3_8*Br_phi3totautau;
    gg_phi2_tautau_TH13=SigmaggF_phi2_13*Br_phi2totautau;
    gg_phi3_tautau_TH13=SigmaggF_phi3_13*Br_phi3totautau;
    bb_phi2_tautau_TH13=SigmabbF_phi2_13*Br_phi2totautau;
    bb_phi3_tautau_TH13=SigmabbF_phi3_13*Br_phi3totautau;

    gg_phi2_gaga_TH8=SigmaggF_phi2_8*Br_phi2togaga;
    gg_phi3_gaga_TH8=SigmaggF_phi3_8*Br_phi3togaga;
    pp_phi2_gaga_TH13=SigmaSumphi2_13*Br_phi2togaga;
    pp_phi3_gaga_TH13=SigmaSumphi3_13*Br_phi3togaga;
    gg_phi2_gaga_TH13=SigmaggF_phi2_13*Br_phi2togaga;
    gg_phi3_gaga_TH13=SigmaggF_phi3_13*Br_phi3togaga;

    pp_phi2_Zga_llga_TH8=SigmaSumphi2_8*Br_phi2toZga*(Br_Ztoee+Br_Ztomumu);
    pp_phi3_Zga_llga_TH8=SigmaSumphi3_8*Br_phi3toZga*(Br_Ztoee+Br_Ztomumu);
    gg_phi2_Zga_TH13=SigmaggF_phi2_13*Br_phi2toZga;
    gg_phi3_Zga_TH13=SigmaggF_phi3_13*Br_phi3toZga;

    gg_phi2_ZZ_TH8=SigmaggF_phi2_8*Br_phi2toZZ;
    gg_phi3_ZZ_TH8=SigmaggF_phi3_8*Br_phi3toZZ;
    VV_phi2_ZZ_TH8=SigmaVBF_phi2_8*Br_phi2toZZ;
    VV_phi3_ZZ_TH8=SigmaVBF_phi3_8*Br_phi3toZZ;
    gg_phi2_ZZ_TH13=SigmaggF_phi2_13*Br_phi2toZZ;
    gg_phi3_ZZ_TH13=SigmaggF_phi3_13*Br_phi3toZZ;
    VV_phi2_ZZ_TH13=SigmaVBF_phi2_13*Br_phi2toZZ;
    VV_phi3_ZZ_TH13=SigmaVBF_phi3_13*Br_phi3toZZ;
    pp_phi2_ZZ_TH13=SigmaSumphi2_13*Br_phi2toZZ;
    pp_phi3_ZZ_TH13=SigmaSumphi3_13*Br_phi3toZZ;

    gg_phi2_WW_TH8=SigmaggF_phi2_8*Br_phi2toWW;
    gg_phi3_WW_TH8=SigmaggF_phi3_8*Br_phi3toWW;   
    VV_phi2_WW_TH8=SigmaVBF_phi2_8*Br_phi2toWW;
    VV_phi3_WW_TH8=SigmaVBF_phi3_8*Br_phi3toWW;
    gg_phi2_WW_TH13=SigmaggF_phi2_13*Br_phi2toWW;
    gg_phi3_WW_TH13=SigmaggF_phi3_13*Br_phi3toWW;
    VV_phi2_WW_TH13=SigmaVBF_phi2_13*Br_phi2toWW;
    VV_phi3_WW_TH13=SigmaVBF_phi3_13*Br_phi3toWW;
    ggVV_phi2_WW_lnulnu_TH13=(SigmaggF_phi2_13+SigmaVBF_phi2_13)*Br_phi2toWW*(Br_Wtoenu+Br_Wtomunu)*(Br_Wtoenu+Br_Wtomunu);
    ggVV_phi3_WW_lnulnu_TH13=(SigmaggF_phi3_13+SigmaVBF_phi3_13)*Br_phi3toWW*(Br_Wtoenu+Br_Wtomunu)*(Br_Wtoenu+Br_Wtomunu);
    pp_phi2_WW_TH13=SigmaSumphi2_13*Br_phi2toWW;
    pp_phi3_WW_TH13=SigmaSumphi3_13*Br_phi3toWW;

    mu_pp_phi2_VV_TH8=SigmaSumphi2_8/SigmaTotSM_phi2_8*rphi2_VV*Gammaphi2totSM/Gammaphi2tot;
    mu_pp_phi3_VV_TH8=SigmaSumphi3_8/SigmaTotSM_phi3_8*rphi3_VV*Gammaphi3totSM/Gammaphi3tot;
    pp_phi2_VV_TH13=SigmaSumphi2_13*(Br_phi2toZZ+Br_phi2toWW);
    pp_phi3_VV_TH13=SigmaSumphi3_13*(Br_phi3toZZ+Br_phi3toWW);

    gg_phi2_phi1phi1_TH8=SigmaggF_phi2_8*Br_phi2tophi1phi1; 
    gg_phi3_phi1phi1_TH8=SigmaggF_phi3_8*Br_phi3tophi1phi1;
    pp_phi2_phi1phi1_bbbb_TH8=SigmaSumphi2_8*Br_phi2tophi1phi1*GTHDM_BR_h_bb*GTHDM_BR_h_bb;
    pp_phi3_phi1phi1_bbbb_TH8=SigmaSumphi3_8*Br_phi3tophi1phi1*GTHDM_BR_h_bb*GTHDM_BR_h_bb;
    pp_phi2_phi1phi1_bbgaga_TH8=SigmaSumphi2_8*Br_phi2tophi1phi1*GTHDM_BR_h_bb*GTHDM_BR_h_gaga*2.0;
    pp_phi3_phi1phi1_bbgaga_TH8=SigmaSumphi3_8*Br_phi3tophi1phi1*GTHDM_BR_h_bb*GTHDM_BR_h_gaga*2.0;
    gg_phi2_phi1phi1_bbtautau_TH8=SigmaggF_phi2_8*Br_phi2tophi1phi1*GTHDM_BR_h_bb*GTHDM_BR_h_tautau*2.0;
    gg_phi3_phi1phi1_bbtautau_TH8=SigmaggF_phi3_8*Br_phi3tophi1phi1*GTHDM_BR_h_bb*GTHDM_BR_h_tautau*2.0;
    pp_phi2_phi1phi1_TH8=SigmaSumphi2_8*Br_phi2tophi1phi1;
    pp_phi3_phi1phi1_TH8=SigmaSumphi3_8*Br_phi3tophi1phi1;
    pp_phi2_phi1phi1_bbbb_TH13=SigmaSumphi2_13*Br_phi2tophi1phi1*GTHDM_BR_h_bb*GTHDM_BR_h_bb;
    pp_phi3_phi1phi1_bbbb_TH13=SigmaSumphi3_13*Br_phi3tophi1phi1*GTHDM_BR_h_bb*GTHDM_BR_h_bb;
    pp_phi2_phi1phi1_TH13=SigmaSumphi2_13*Br_phi2tophi1phi1;
    pp_phi3_phi1phi1_TH13=SigmaSumphi3_13*Br_phi3tophi1phi1;
    pp_phi2_phi1phi1_bbgaga_TH13=SigmaSumphi2_13*Br_phi2tophi1phi1*GTHDM_BR_h_bb*GTHDM_BR_h_gaga*2.0;
    pp_phi3_phi1phi1_bbgaga_TH13=SigmaSumphi3_13*Br_phi3tophi1phi1*GTHDM_BR_h_bb*GTHDM_BR_h_gaga*2.0;
    pp_phi2_phi1phi1_bbtautau_TH13=SigmaSumphi2_13*Br_phi2tophi1phi1*GTHDM_BR_h_bb*GTHDM_BR_h_tautau*2.0;
    pp_phi3_phi1phi1_bbtautau_TH13=SigmaSumphi3_13*Br_phi3tophi1phi1*GTHDM_BR_h_bb*GTHDM_BR_h_tautau*2.0;
    pp_phi2_phi1phi1_bbVV_TH13=SigmaSumphi2_13*Br_phi2tophi1phi1*2.0*GTHDM_BR_h_bb*(GTHDM_BR_h_WW*pow(Br_Wtoenu+Br_Wtomunu+Br_Wtotaunu*0.352,2)
                                                            +GTHDM_BR_h_ZZ*2.0*Br_Ztoinv*(Br_Ztoee+Br_Ztomumu+Br_Ztotautau*0.124));
    pp_phi3_phi1phi1_bbVV_TH13=SigmaSumphi3_13*Br_phi3tophi1phi1*2.0*GTHDM_BR_h_bb*(GTHDM_BR_h_WW*pow(Br_Wtoenu+Br_Wtomunu+Br_Wtotaunu*0.352,2)
                                                            +GTHDM_BR_h_ZZ*2.0*Br_Ztoinv*(Br_Ztoee+Br_Ztomumu+Br_Ztotautau*0.124));  
    pp_phi2_phi1phi1_bbWW_TH13=SigmaSumphi2_13*Br_phi2tophi1phi1*2.0*5.77e-1*2.15e-1; //SM Br of hh assumed
    pp_phi3_phi1phi1_bbWW_TH13=SigmaSumphi3_13*Br_phi3tophi1phi1*2.0*5.77e-1*2.15e-1; //SM Br of hh assumed
    gg_phi2_phi1phi1_gagaWW_TH13=SigmaggF_phi2_13*Br_phi2tophi1phi1*GTHDM_BR_h_gaga*GTHDM_BR_h_WW*2.0;
    gg_phi3_phi1phi1_gagaWW_TH13=SigmaggF_phi3_13*Br_phi3tophi1phi1*GTHDM_BR_h_gaga*GTHDM_BR_h_WW*2.0;

    gg_phi2_phi1Z_bbZ_TH8=SigmaggF_phi2_8*Br_phi2tophi1Z*GTHDM_BR_h_bb;
    gg_phi3_phi1Z_bbZ_TH8=SigmaggF_phi3_8*Br_phi3tophi1Z*GTHDM_BR_h_bb;
    gg_phi2_phi1Z_bbll_TH8=SigmaggF_phi2_8*Br_phi2tophi1Z*GTHDM_BR_h_bb*(Br_Ztoee+Br_Ztomumu);
    gg_phi3_phi1Z_bbll_TH8=SigmaggF_phi3_8*Br_phi3tophi1Z*GTHDM_BR_h_bb*(Br_Ztoee+Br_Ztomumu);
    gg_phi2_phi1Z_tautauZ_TH8=SigmaggF_phi2_8*Br_phi2tophi1Z*GTHDM_BR_h_tautau;
    gg_phi3_phi1Z_tautauZ_TH8=SigmaggF_phi3_8*Br_phi3tophi1Z*GTHDM_BR_h_tautau;
    gg_phi2_phi1Z_tautaull_TH8=SigmaggF_phi2_8*Br_phi2tophi1Z*GTHDM_BR_h_tautau*(Br_Ztoee+Br_Ztomumu);
    gg_phi3_phi1Z_tautaull_TH8=SigmaggF_phi3_8*Br_phi3tophi1Z*GTHDM_BR_h_tautau*(Br_Ztoee+Br_Ztomumu);
    gg_phi2_phi1Z_bbZ_TH13=SigmaggF_phi2_13*Br_phi2tophi1Z*GTHDM_BR_h_bb;
    gg_phi3_phi1Z_bbZ_TH13=SigmaggF_phi3_13*Br_phi3tophi1Z*GTHDM_BR_h_bb;
    bb_phi2_phi1Z_bbZ_TH13=SigmabbF_phi2_13*Br_phi2tophi1Z*GTHDM_BR_h_bb;
    bb_phi3_phi1Z_bbZ_TH13=SigmabbF_phi3_13*Br_phi3tophi1Z*GTHDM_BR_h_bb;
  
    pp_phi3_phi2Z_bbll_TH8=SigmaSumphi3_8*Br_phi3tophi2Z*Br_phi2tobb*(Br_Ztoee+Br_Ztomumu);
    pp_phi2_phi3Z_bbll_TH8=SigmaSumphi2_8*Br_phi2tophi3Z*Br_phi3tobb*(Br_Ztoee+Br_Ztomumu);     
    pp_phi3_phi2Z_tautaull_TH8=SigmaSumphi3_8*Br_phi3tophi2Z*Br_phi2totautau*(Br_Ztoee+Br_Ztomumu);
    pp_phi2_phi3Z_tautaull_TH8=SigmaSumphi2_8*Br_phi2tophi3Z*Br_phi3totautau*(Br_Ztoee+Br_Ztomumu);
    gg_phi3_phi2Z_bbZ_TH13=SigmaggF_phi3_13*Br_phi3tophi2Z*Br_phi2tobb;
    gg_phi2_phi3Z_bbZ_TH13=SigmaggF_phi2_13*Br_phi2tophi3Z*Br_phi3tobb;
    bb_phi3_phi2Z_bbZ_TH13=SigmabbF_phi3_13*Br_phi3tophi2Z*Br_phi2tobb;
    bb_phi2_phi3Z_bbZ_TH13=SigmabbF_phi2_13*Br_phi2tophi3Z*Br_phi3tobb;
    
    pp_Hpm_taunu_TH8=2.0*SigmaHp8*Br_Hptotaunu;
    pp_Hp_taunu_TH8=SigmaHp8*Br_Hptotaunu;
    pp_Hpm_taunu_TH13=2.0*SigmaHpm13*Br_Hptotaunu;
    pp_Hpm_tb_TH8=2.0*SigmaHp8*Br_Hptotb;
    pp_Hp_tb_TH8=SigmaHp8*Br_Hptotb;
    pp_Hpm_tb_TH13=2.0*SigmaHpm13*Br_Hptotb;

    //95% to 1 sigma conversion factor, roughly sqrt(3.84)
//    double nftos=1.95996398454;

    if(m2>= 400.0 && m2<1000.0) THoEX_tt_phi2_tt_ATLAS13=tt_phi2_tt_TH13/ip_ex_tt_phi_tt_ATLAS13(m2);
    if(m3>= 400.0 && m3<1000.0) THoEX_tt_phi3_tt_ATLAS13=tt_phi3_tt_TH13/ip_ex_tt_phi_tt_ATLAS13(m3);
    if(m2>= 400.0 && m2<1000.0) THoEX_bb_phi2_tt_ATLAS13=bb_phi2_tt_TH13/ip_ex_bb_phi_tt_ATLAS13(m2);
    if(m3>= 400.0 && m3<1000.0) THoEX_bb_phi3_tt_ATLAS13=bb_phi3_tt_TH13/ip_ex_bb_phi_tt_ATLAS13(m3);
    if(m2>= 100.0 && m2< 900.0) THoEX_bb_phi2_bb_CMS8=bb_phi2_bb_TH8/ip_ex_bb_phi_bb_CMS8(m2);
    if(m3>= 100.0 && m3< 900.0) THoEX_bb_phi3_bb_CMS8=bb_phi3_bb_TH8/ip_ex_bb_phi_bb_CMS8(m3);
    if(m2>= 330.0 && m2<1200.0) THoEX_gg_phi2_bb_CMS8=gg_phi2_bb_TH8/ip_ex_gg_phi_bb_CMS8(m2);
    if(m3>= 330.0 && m3<1200.0) THoEX_gg_phi3_bb_CMS8=gg_phi3_bb_TH8/ip_ex_gg_phi_bb_CMS8(m3);   
    if(m2>= 550.0 && m2<1200.0) THoEX_pp_phi2_bb_CMS13=pp_phi2_bb_TH13/ip_ex_pp_phi_bb_CMS13(m2);
    if(m3>= 550.0 && m3<1200.0) THoEX_pp_phi3_bb_CMS13=pp_phi3_bb_TH13/ip_ex_pp_phi_bb_CMS13(m3);
    if(m2>= 300.0 && m2<1300.0) THoEX_bb_phi2_bb_CMS13=bb_phi2_bb_TH13/ip_ex_bb_phi_bb_CMS13(m2);
    if(m3>= 300.0 && m3<1300.0) THoEX_bb_phi3_bb_CMS13=bb_phi3_bb_TH13/ip_ex_bb_phi_bb_CMS13(m3);
    if(m2>=  90.0 && m2<1000.0) THoEX_gg_phi2_tautau_ATLAS8=gg_phi2_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(m2);
    if(m3>=  90.0 && m3<1000.0) THoEX_gg_phi3_tautau_ATLAS8=gg_phi3_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(m3);
    if(m2>=  90.0 && m2<1000.0) THoEX_gg_phi2_tautau_CMS8=gg_phi2_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(m2);
    if(m3>=  90.0 && m3<1000.0) THoEX_gg_phi3_tautau_CMS8=gg_phi3_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(m3);
    if(m2>=  90.0 && m2<1000.0) THoEX_bb_phi2_tautau_ATLAS8=bb_phi2_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(m2);
    if(m3>=  90.0 && m3<1000.0) THoEX_bb_phi3_tautau_ATLAS8=bb_phi3_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(m3);
    if(m2>=  90.0 && m2<1000.0) THoEX_bb_phi2_tautau_CMS8=bb_phi2_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(m2);
    if(m3>=  90.0 && m3<1000.0) THoEX_bb_phi3_tautau_CMS8=bb_phi3_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(m3);
    if(m2>= 200.0 && m2<2250.0) THoEX_gg_phi2_tautau_ATLAS13=gg_phi2_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(m2);
    if(m3>= 200.0 && m3<2250.0) THoEX_gg_phi3_tautau_ATLAS13=gg_phi3_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(m3);
    if(m2>=  90.0 && m2<3200.0) THoEX_gg_phi2_tautau_CMS13=gg_phi2_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(m2);
    if(m3>=  90.0 && m3<3200.0) THoEX_gg_phi3_tautau_CMS13=gg_phi3_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(m3);
    if(m2>= 200.0 && m2<2250.0) THoEX_bb_phi2_tautau_ATLAS13=bb_phi2_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(m2);
    if(m3>= 200.0 && m3<2250.0) THoEX_bb_phi3_tautau_ATLAS13=bb_phi3_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(m3);
    if(m2>=  90.0 && m2<3200.0) THoEX_bb_phi2_tautau_CMS13=bb_phi2_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(m2);
    if(m3>=  90.0 && m3<3200.0) THoEX_bb_phi3_tautau_CMS13=bb_phi3_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(m3);
    if(m2>=  65.0 && m2< 600.0) THoEX_gg_phi2_gaga_ATLAS8=gg_phi2_gaga_TH8/ip_ex_gg_phi_gaga_ATLAS8(m2);
    if(m3>=  65.0 && m3< 600.0) THoEX_gg_phi3_gaga_ATLAS8=gg_phi3_gaga_TH8/ip_ex_gg_phi_gaga_ATLAS8(m3);
    if(m2>= 200.0 && m2<2700.0) THoEX_pp_phi2_gaga_ATLAS13=pp_phi2_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(m2);
    if(m3>= 200.0 && m3<2700.0) THoEX_pp_phi3_gaga_ATLAS13=pp_phi3_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(m3);   
    if(m2>= 500.0 && m2<4000.0) THoEX_gg_phi2_gaga_CMS13=gg_phi2_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(m2);    
    if(m3>= 500.0 && m3<4000.0) THoEX_gg_phi3_gaga_CMS13=gg_phi3_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(m3);
    if(m2>= 200.0 && m2<1600.0) THoEX_pp_phi2_Zga_llga_ATLAS8=pp_phi2_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(m2);
    if(m3>= 200.0 && m3<1600.0) THoEX_pp_phi3_Zga_llga_ATLAS8=pp_phi3_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(m3);
    if(m2>= 200.0 && m2<1200.0) THoEX_pp_phi2_Zga_llga_CMS8=pp_phi2_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_CMS8(m2);
    if(m3>= 200.0 && m3<1200.0) THoEX_pp_phi3_Zga_llga_CMS8=pp_phi3_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_CMS8(m3);
    if(m2>= 250.0 && m2<2400.0) THoEX_gg_phi2_Zga_llga_ATLAS13=gg_phi2_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(m2);
    if(m3>= 250.0 && m3<2400.0) THoEX_gg_phi3_Zga_llga_ATLAS13=gg_phi3_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(m3);
    if(m2>=1000.0 && m2<6800.0) THoEX_gg_phi2_Zga_qqga_ATLAS13=gg_phi2_Zga_TH13/ip_ex_gg_phi_Zga_qqga_ATLAS13(m2);
    if(m3>=1000.0 && m3<6800.0) THoEX_gg_phi3_Zga_qqga_ATLAS13=gg_phi3_Zga_TH13/ip_ex_gg_phi_Zga_qqga_ATLAS13(m3);
    if(m2>= 350.0 && m2<4000.0) THoEX_gg_phi2_Zga_CMS13=gg_phi2_Zga_TH13/ip_ex_gg_phi_Zga_CMS13(m2);
    if(m3>= 350.0 && m3<4000.0) THoEX_gg_phi3_Zga_CMS13=gg_phi3_Zga_TH13/ip_ex_gg_phi_Zga_CMS13(m3);
    if(m2>= 140.0 && m2<1000.0) THoEX_gg_phi2_ZZ_ATLAS8=gg_phi2_ZZ_TH8/ip_ex_gg_phi_ZZ_ATLAS8(m2);
    if(m3>= 140.0 && m3<1000.0) THoEX_gg_phi3_ZZ_ATLAS8=gg_phi3_ZZ_TH8/ip_ex_gg_phi_ZZ_ATLAS8(m3);
    if(m2>= 140.0 && m2<1000.0) THoEX_VV_phi2_ZZ_ATLAS8=VV_phi2_ZZ_TH8/ip_ex_VV_phi_ZZ_ATLAS8(m2);
    if(m3>= 140.0 && m3<1000.0) THoEX_VV_phi3_ZZ_ATLAS8=VV_phi3_ZZ_TH8/ip_ex_VV_phi_ZZ_ATLAS8(m3);
    if(m2>= 200.0 && m2<1200.0) THoEX_gg_phi2_ZZ_llllnunu_ATLAS13=gg_phi2_ZZ_TH13/ip_ex_gg_phi_ZZ_llllnunu_ATLAS13(m2);
    if(m3>= 200.0 && m3<1200.0) THoEX_gg_phi3_ZZ_llllnunu_ATLAS13=gg_phi3_ZZ_TH13/ip_ex_gg_phi_ZZ_llllnunu_ATLAS13(m3);
    if(m2>= 200.0 && m2<1200.0) THoEX_VV_phi2_ZZ_llllnunu_ATLAS13=VV_phi2_ZZ_TH13/ip_ex_VV_phi_ZZ_llllnunu_ATLAS13(m2);
    if(m3>= 200.0 && m3<1200.0) THoEX_VV_phi3_ZZ_llllnunu_ATLAS13=VV_phi3_ZZ_TH13/ip_ex_VV_phi_ZZ_llllnunu_ATLAS13(m3);
    if(m2>= 300.0 && m2<3000.0) THoEX_gg_phi2_ZZ_qqllnunu_ATLAS13=gg_phi2_ZZ_TH13/ip_ex_gg_phi_ZZ_qqllnunu_ATLAS13(m2);
    if(m3>= 300.0 && m3<3000.0) THoEX_gg_phi3_ZZ_qqllnunu_ATLAS13=gg_phi3_ZZ_TH13/ip_ex_gg_phi_ZZ_qqllnunu_ATLAS13(m3);
    if(m2>= 300.0 && m2<3000.0) THoEX_VV_phi2_ZZ_qqllnunu_ATLAS13=VV_phi2_ZZ_TH13/ip_ex_VV_phi_ZZ_qqllnunu_ATLAS13(m2);
    if(m3>= 300.0 && m3<3000.0) THoEX_VV_phi3_ZZ_qqllnunu_ATLAS13=VV_phi3_ZZ_TH13/ip_ex_VV_phi_ZZ_qqllnunu_ATLAS13(m3);
    if(m2>= 130.0 && m2<3000.0) THoEX_pp_phi2_ZZ_llqqnunull_CMS13=pp_phi2_ZZ_TH13/ip_ex_pp_phi_ZZ_llqqnunull_CMS13(m2);
    if(m3>= 130.0 && m3<3000.0) THoEX_pp_phi3_ZZ_llqqnunull_CMS13=pp_phi3_ZZ_TH13/ip_ex_pp_phi_ZZ_llqqnunull_CMS13(m3);
    if(m2>=1000.0 && m2<4000.0) THoEX_pp_phi2_ZZ_qqnunu_CMS13=pp_phi2_ZZ_TH13/ip_ex_pp_phi_ZZ_qqnunu_CMS13(m2);
    if(m3>=1000.0 && m3<4000.0) THoEX_pp_phi3_ZZ_qqnunu_CMS13=pp_phi3_ZZ_TH13/ip_ex_pp_phi_ZZ_qqnunu_CMS13(m3);
    if(m2>= 300.0 && m2<1500.0) THoEX_gg_phi2_WW_ATLAS8=gg_phi2_WW_TH8/ip_ex_gg_phi_WW_ATLAS8(m2);
    if(m3>= 300.0 && m3<1500.0) THoEX_gg_phi3_WW_ATLAS8=gg_phi3_WW_TH8/ip_ex_gg_phi_WW_ATLAS8(m3);  
    if(m2>= 300.0 && m2<1500.0) THoEX_VV_phi2_WW_ATLAS8=VV_phi2_WW_TH8/ip_ex_VV_phi_WW_ATLAS8(m2);
    if(m3>= 300.0 && m3<1500.0) THoEX_VV_phi3_WW_ATLAS8=VV_phi3_WW_TH8/ip_ex_VV_phi_WW_ATLAS8(m3);
    if(m2>= 250.0 && m2<4000.0) THoEX_gg_phi2_WW_enumunu_ATLAS13=gg_phi2_WW_TH13/ip_ex_gg_phi_WW_enumunu_ATLAS13(m2);
    if(m3>= 250.0 && m3<4000.0) THoEX_gg_phi3_WW_enumunu_ATLAS13=gg_phi3_WW_TH13/ip_ex_gg_phi_WW_enumunu_ATLAS13(m3);
    if(m2>= 250.0 && m2<3000.0) THoEX_VV_phi2_WW_enumunu_ATLAS13=VV_phi2_WW_TH13/ip_ex_VV_phi_WW_enumunu_ATLAS13(m2);
    if(m3>= 250.0 && m3<3000.0) THoEX_VV_phi3_WW_enumunu_ATLAS13=VV_phi3_WW_TH13/ip_ex_VV_phi_WW_enumunu_ATLAS13(m3);
    if(m2>= 200.0 && m2<1000.0) THoEX_ggVV_phi2_WW_lnulnu_CMS13=ggVV_phi2_WW_lnulnu_TH13/ip_ex_ggVV_phi_WW_lnulnu_CMS13(m2);
    if(m3>= 200.0 && m3<1000.0) THoEX_ggVV_phi3_WW_lnulnu_CMS13=ggVV_phi3_WW_lnulnu_TH13/ip_ex_ggVV_phi_WW_lnulnu_CMS13(m3);
    if(m2>= 300.0 && m2<3000.0) THoEX_gg_phi2_WW_lnuqq_ATLAS13=gg_phi2_WW_TH13/ip_ex_gg_phi_WW_lnuqq_ATLAS13(m2);
    if(m3>= 300.0 && m3<3000.0) THoEX_gg_phi3_WW_lnuqq_ATLAS13=gg_phi3_WW_TH13/ip_ex_gg_phi_WW_lnuqq_ATLAS13(m3);
    if(m2>= 300.0 && m2<3000.0) THoEX_VV_phi2_WW_lnuqq_ATLAS13=VV_phi2_WW_TH13/ip_ex_VV_phi_WW_lnuqq_ATLAS13(m2);
    if(m3>= 300.0 && m3<3000.0) THoEX_VV_phi3_WW_lnuqq_ATLAS13=VV_phi3_WW_TH13/ip_ex_VV_phi_WW_lnuqq_ATLAS13(m3);
    if(m2>=1000.0 && m2<4400.0) THoEX_pp_phi2_WW_lnuqq_CMS13=pp_phi2_WW_TH13/ip_ex_pp_phi_WW_lnuqq_CMS13(m2);
    if(m3>=1000.0 && m3<4400.0) THoEX_pp_phi3_WW_lnuqq_CMS13=pp_phi3_WW_TH13/ip_ex_pp_phi_WW_lnuqq_CMS13(m3);
    if(m2>= 145.0 && m2<1000.0) THoEX_mu_pp_phi2_VV_CMS8=mu_pp_phi2_VV_TH8/ip_ex_mu_pp_phi_VV_CMS8(m2);
    if(m3>= 145.0 && m3<1000.0) THoEX_mu_pp_phi3_VV_CMS8=mu_pp_phi3_VV_TH8/ip_ex_mu_pp_phi_VV_CMS8(m3);
    if(m2>=1200.0 && m2<3000.0) THoEX_pp_phi2_VV_qqqq_ATLAS13=pp_phi2_VV_TH13/ip_ex_pp_phi_VV_qqqq_ATLAS13(m2);
    if(m3>=1200.0 && m3<3000.0) THoEX_pp_phi3_VV_qqqq_ATLAS13=pp_phi3_VV_TH13/ip_ex_pp_phi_VV_qqqq_ATLAS13(m3);
    if(m2>= 260.0 && m2<1000.0) THoEX_gg_phi2_phi1phi1_ATLAS8=gg_phi2_phi1phi1_TH8/ip_ex_gg_phi_phi1phi1_ATLAS8(m2);   
    if(m3>= 260.0 && m3<1000.0) THoEX_gg_phi3_phi1phi1_ATLAS8=gg_phi3_phi1phi1_TH8/ip_ex_gg_phi_phi1phi1_ATLAS8(m3);
    if(m2>= 270.0 && m2<1100.0) THoEX_pp_phi2_phi1phi1_bbbb_CMS8=pp_phi2_phi1phi1_bbbb_TH8/ip_ex_pp_phi_phi1phi1_bbbb_CMS8(m2);
    if(m3>= 270.0 && m3<1100.0) THoEX_pp_phi3_phi1phi1_bbbb_CMS8=pp_phi3_phi1phi1_bbbb_TH8/ip_ex_pp_phi_phi1phi1_bbbb_CMS8(m3);
    if(m2>= 260.0 && m2<1100.0) THoEX_pp_phi2_phi1phi1_bbgaga_CMS8=pp_phi2_phi1phi1_bbgaga_TH8/ip_ex_pp_phi_phi1phi1_bbgaga_CMS8(m2);
    if(m3>= 260.0 && m3<1100.0) THoEX_pp_phi3_phi1phi1_bbgaga_CMS8=pp_phi3_phi1phi1_bbgaga_TH8/ip_ex_pp_phi_phi1phi1_bbgaga_CMS8(m3);
    if(m2>= 260.0 && m2< 350.0) THoEX_gg_phi2_phi1phi1_bbtautau_CMS8=gg_phi2_phi1phi1_bbtautau_TH8/ip_ex_gg_phi_phi1phi1_bbtautau_CMS8(m2);
    if(m3>= 260.0 && m3< 350.0) THoEX_gg_phi3_phi1phi1_bbtautau_CMS8=gg_phi3_phi1phi1_bbtautau_TH8/ip_ex_gg_phi_phi1phi1_bbtautau_CMS8(m3);
    if(m2>= 350.0 && m2<1000.0) THoEX_pp_phi2_phi1phi1_bbtautau_CMS8=pp_phi2_phi1phi1_TH8/ip_ex_pp_phi_phi1phi1_bbtautau_CMS8(m2);
    if(m3>= 350.0 && m3<1000.0) THoEX_pp_phi3_phi1phi1_bbtautau_CMS8=pp_phi3_phi1phi1_TH8/ip_ex_pp_phi_phi1phi1_bbtautau_CMS8(m3);
    if(m2>= 260.0 && m2<3000.0) THoEX_pp_phi2_phi1phi1_bbbb_ATLAS13=pp_phi2_phi1phi1_bbbb_TH13/ip_ex_pp_phi_phi1phi1_bbbb_ATLAS13(m2);
    if(m3>= 260.0 && m3<3000.0) THoEX_pp_phi3_phi1phi1_bbbb_ATLAS13=pp_phi3_phi1phi1_bbbb_TH13/ip_ex_pp_phi_phi1phi1_bbbb_ATLAS13(m3);
    if(m2>= 260.0 && m2<1200.0) THoEX_pp_phi2_phi1phi1_bbbb_1_CMS13=pp_phi2_phi1phi1_bbbb_TH13/ip_ex_pp_phi_phi1phi1_bbbb_1_CMS13(m2);
    if(m3>= 260.0 && m3<1200.0) THoEX_pp_phi3_phi1phi1_bbbb_1_CMS13=pp_phi3_phi1phi1_bbbb_TH13/ip_ex_pp_phi_phi1phi1_bbbb_1_CMS13(m3);
    if(m2>=1200.0 && m2<3000.0) THoEX_pp_phi2_phi1phi1_bbbb_2_CMS13=pp_phi2_phi1phi1_bbbb_TH13/ip_ex_pp_phi_phi1phi1_bbbb_2_CMS13(m2);
    if(m3>=1200.0 && m3<3000.0) THoEX_pp_phi3_phi1phi1_bbbb_2_CMS13=pp_phi3_phi1phi1_bbbb_TH13/ip_ex_pp_phi_phi1phi1_bbbb_2_CMS13(m3);
    if(m2>= 260.0 && m2<1000.0) THoEX_pp_phi2_phi1phi1_bbgaga_ATLAS13=pp_phi2_phi1phi1_TH13/ip_ex_pp_phi_phi1phi1_bbgaga_ATLAS13(m2);
    if(m3>= 260.0 && m3<1000.0) THoEX_pp_phi3_phi1phi1_bbgaga_ATLAS13=pp_phi3_phi1phi1_TH13/ip_ex_pp_phi_phi1phi1_bbgaga_ATLAS13(m3);
    if(m2>= 250.0 && m2< 900.0) THoEX_pp_phi2_phi1phi1_bbgaga_CMS13=pp_phi2_phi1phi1_bbgaga_TH13/ip_ex_pp_phi_phi1phi1_bbgaga_CMS13(m2);
    if(m3>= 250.0 && m3< 900.0) THoEX_pp_phi3_phi1phi1_bbgaga_CMS13=pp_phi3_phi1phi1_bbgaga_TH13/ip_ex_pp_phi_phi1phi1_bbgaga_CMS13(m3);
    if(m2>= 260.0 && m2<1000.0) THoEX_pp_phi2_phi1phi1_bbtautau_ATLAS13=pp_phi2_phi1phi1_bbtautau_TH13/ip_ex_pp_phi_phi1phi1_bbtautau_ATLAS13(m2);
    if(m3>= 260.0 && m3<1000.0) THoEX_pp_phi3_phi1phi1_bbtautau_ATLAS13=pp_phi3_phi1phi1_bbtautau_TH13/ip_ex_pp_phi_phi1phi1_bbtautau_ATLAS13(m3);
    if(m2>= 250.0 && m2< 900.0) THoEX_pp_phi2_phi1phi1_bbtautau_1_CMS13=pp_phi2_phi1phi1_bbtautau_TH13/ip_ex_pp_phi_phi1phi1_bbtautau_1_CMS13(m2);
    if(m3>= 250.0 && m3< 900.0) THoEX_pp_phi3_phi1phi1_bbtautau_1_CMS13=pp_phi3_phi1phi1_bbtautau_TH13/ip_ex_pp_phi_phi1phi1_bbtautau_1_CMS13(m3);
    if(m2>= 900.0 && m2<4000.0) THoEX_pp_phi2_phi1phi1_bbtautau_2_CMS13=pp_phi2_phi1phi1_TH13/ip_ex_pp_phi_phi1phi1_bbtautau_2_CMS13(m2);
    if(m3>= 900.0 && m3<4000.0) THoEX_pp_phi3_phi1phi1_bbtautau_2_CMS13=pp_phi3_phi1phi1_TH13/ip_ex_pp_phi_phi1phi1_bbtautau_2_CMS13(m3);
    if(m2>= 260.0 && m2< 900.0) THoEX_pp_phi2_phi1phi1_bbVV_CMS13=pp_phi2_phi1phi1_bbVV_TH13/ip_ex_pp_phi_phi1phi1_bbVV_CMS13(m2);
    if(m3>= 260.0 && m3< 900.0) THoEX_pp_phi3_phi1phi1_bbVV_CMS13=pp_phi3_phi1phi1_bbVV_TH13/ip_ex_pp_phi_phi1phi1_bbVV_CMS13(m3); 
    if(m2>= 500.0 && m2< 3000.0) THoEX_pp_phi2_phi1phi1_bbWW_ATLAS13=pp_phi2_phi1phi1_bbWW_TH13/ip_ex_pp_phi_phi1phi1_bbWW_ATLAS13(m2);
    if(m3>= 500.0 && m3< 3000.0) THoEX_pp_phi3_phi1phi1_bbWW_ATLAS13=pp_phi3_phi1phi1_bbWW_TH13/ip_ex_pp_phi_phi1phi1_bbWW_ATLAS13(m3);  
    if(m2>= 260.0 && m2< 500.0) THoEX_gg_phi2_phi1phi1_gagaWW_ATLAS13=gg_phi2_phi1phi1_gagaWW_TH13/ip_ex_gg_phi_phi1phi1_gagaWW_ATLAS13(m2);
    if(m3>= 260.0 && m3< 500.0) THoEX_gg_phi3_phi1phi1_gagaWW_ATLAS13=gg_phi3_phi1phi1_gagaWW_TH13/ip_ex_gg_phi_phi1phi1_gagaWW_ATLAS13(m3);
    
    
    if(m2>= 220.0 && m2<1000.0) THoEX_gg_phi2_phi1Z_bbZ_ATLAS8=gg_phi2_phi1Z_bbZ_TH8/ip_ex_gg_phi_phi1Z_bbZ_ATLAS8(m2);
    if(m3>= 220.0 && m3<1000.0) THoEX_gg_phi3_phi1Z_bbZ_ATLAS8=gg_phi3_phi1Z_bbZ_TH8/ip_ex_gg_phi_phi1Z_bbZ_ATLAS8(m3);
    if(m2>= 225.0 && m2< 600.0) THoEX_gg_phi2_phi1Z_bbll_CMS8=gg_phi2_phi1Z_bbll_TH8/ip_ex_gg_phi_phi1Z_bbll_CMS8(m2);
    if(m3>= 225.0 && m3< 600.0) THoEX_gg_phi3_phi1Z_bbll_CMS8=gg_phi3_phi1Z_bbll_TH8/ip_ex_gg_phi_phi1Z_bbll_CMS8(m3);
    if(m2>= 220.0 && m2<1000.0) THoEX_gg_phi2_phi1Z_tautauZ_ATLAS8=gg_phi2_phi1Z_tautauZ_TH8/ip_ex_gg_phi_phi1Z_tautauZ_ATLAS8(m2);
    if(m3>= 220.0 && m3<1000.0) THoEX_gg_phi3_phi1Z_tautauZ_ATLAS8=gg_phi3_phi1Z_tautauZ_TH8/ip_ex_gg_phi_phi1Z_tautauZ_ATLAS8(m3);
    if(m2>= 220.0 && m2< 350.0) THoEX_gg_phi2_phi1Z_tautaull_CMS8=gg_phi2_phi1Z_tautaull_TH8/ip_ex_gg_phi_phi1Z_tautaull_CMS8(m2);
    if(m3>= 220.0 && m3< 350.0) THoEX_gg_phi3_phi1Z_tautaull_CMS8=gg_phi3_phi1Z_tautaull_TH8/ip_ex_gg_phi_phi1Z_tautaull_CMS8(m3);
    if(m2>= 200.0 && m2<2000.0) THoEX_gg_phi2_phi1Z_bbZ_ATLAS13=gg_phi2_phi1Z_bbZ_TH13/ip_ex_gg_phi_phi1Z_bbZ_ATLAS13(m2);
    if(m3>= 200.0 && m3<2000.0) THoEX_gg_phi3_phi1Z_bbZ_ATLAS13=gg_phi3_phi1Z_bbZ_TH13/ip_ex_gg_phi_phi1Z_bbZ_ATLAS13(m3);
    if(m2>= 220.0 && m2< 800.0) THoEX_gg_phi2_phi1Z_bbZ_1_CMS13=gg_phi2_phi1Z_bbZ_TH13/ip_ex_gg_phi_phi1Z_bbZ_1_CMS13(m2);
    if(m3>= 220.0 && m3< 800.0) THoEX_gg_phi3_phi1Z_bbZ_1_CMS13=gg_phi3_phi1Z_bbZ_TH13/ip_ex_gg_phi_phi1Z_bbZ_1_CMS13(m3);
    if(m2>= 800.0 && m2<2000.0) THoEX_gg_phi2_phi1Z_bbZ_2_CMS13=gg_phi2_phi1Z_bbZ_TH13/ip_ex_gg_phi_phi1Z_bbZ_2_CMS13(m2);
    if(m3>= 800.0 && m3<2000.0) THoEX_gg_phi3_phi1Z_bbZ_2_CMS13=gg_phi3_phi1Z_bbZ_TH13/ip_ex_gg_phi_phi1Z_bbZ_2_CMS13(m3);
    if(m2>= 200.0 && m2<2000.0) THoEX_bb_phi2_phi1Z_bbZ_ATLAS13=bb_phi2_phi1Z_bbZ_TH13/ip_ex_bb_phi_phi1Z_bbZ_ATLAS13(m2);
    if(m3>= 200.0 && m3<2000.0) THoEX_bb_phi3_phi1Z_bbZ_ATLAS13=bb_phi3_phi1Z_bbZ_TH13/ip_ex_bb_phi_phi1Z_bbZ_ATLAS13(m3);
    if(m2>= 220.0 && m2< 800.0) THoEX_bb_phi2_phi1Z_bbZ_1_CMS13=bb_phi2_phi1Z_bbZ_TH13/ip_ex_bb_phi_phi1Z_bbZ_1_CMS13(m2);
    if(m3>= 220.0 && m3< 800.0) THoEX_bb_phi3_phi1Z_bbZ_1_CMS13=bb_phi3_phi1Z_bbZ_TH13/ip_ex_bb_phi_phi1Z_bbZ_1_CMS13(m3);
    if(m2>= 800.0 && m2<2000.0) THoEX_bb_phi2_phi1Z_bbZ_2_CMS13=bb_phi2_phi1Z_bbZ_TH13/ip_ex_bb_phi_phi1Z_bbZ_2_CMS13(m2);
    if(m3>= 800.0 && m3<2000.0) THoEX_bb_phi3_phi1Z_bbZ_2_CMS13=bb_phi3_phi1Z_bbZ_TH13/ip_ex_bb_phi_phi1Z_bbZ_2_CMS13(m3);
    if(m3>= 175.0 && m3<1000.0 && m2 >=50.0 && m2 <910.0) THoEX_pp_phi3_phi2Z_bbll_1_CMS8=pp_phi3_phi2Z_bbll_TH8/ip_ex_pp_phii_phijZ_bbll_1_CMS8(m3,m2); //mA=m3, mH=m2
    if(m2>= 175.0 && m2<1000.0 && m3 >=50.0 && m3 <910.0) THoEX_pp_phi2_phi3Z_bbll_1_CMS8=pp_phi2_phi3Z_bbll_TH8/ip_ex_pp_phii_phijZ_bbll_1_CMS8(m2,m3); //mA=m3, mH=m2
    if(m3>=  50.0 && m3<1000.0 && m2 >=50.0 && m2 <1000.0) THoEX_pp_phi3_phi2Z_tautaull_1_CMS8=pp_phi3_phi2Z_tautaull_TH8/ip_ex_pp_phii_phijZ_tautaull_1_CMS8(m3,m2); //mA=m3, mH=m2
    if(m2>=  50.0 && m2<1000.0 && m3 >=50.0 && m3 <1000.0) THoEX_pp_phi2_phi3Z_tautaull_1_CMS8=pp_phi2_phi3Z_tautaull_TH8/ip_ex_pp_phii_phijZ_tautaull_1_CMS8(m2,m3); //mA=m3, mH=m2
    if(m3>=  50.0 && m3<1000.0 && m2 >=50.0 && m2 <1000.0) THoEX_pp_phi3_phi2Z_tautaull_2_CMS8=pp_phi3_phi2Z_tautaull_TH8/ip_ex_pp_phii_phijZ_tautaull_2_CMS8(m2,m3); //mA=m2, mH=m3
    if(m2>=  50.0 && m2<1000.0 && m3 >=50.0 && m3 <1000.0) THoEX_pp_phi2_phi3Z_tautaull_2_CMS8=pp_phi2_phi3Z_tautaull_TH8/ip_ex_pp_phii_phijZ_tautaull_2_CMS8(m3,m2); //mA=m2, mH=m3
    if(m3 >= 230.0 && m3 <800.0 && m2>=130.0 && m2<700.0) THoEX_gg_phi3_phi2Z_bbZ_ATLAS13=gg_phi3_phi2Z_bbZ_TH13/ip_ex_gg_phii_phijZ_bbZ_ATLAS13(m3,m2);
    if(m2 >= 230.0 && m2 <800.0 && m3>=130.0 && m3<700.0) THoEX_gg_phi2_phi3Z_bbZ_ATLAS13=gg_phi2_phi3Z_bbZ_TH13/ip_ex_gg_phii_phijZ_bbZ_ATLAS13(m2,m3);
    if(m3 >= 230.0 && m3 <800.0 && m2>=130.0 && m2<700.0) THoEX_bb_phi3_phi2Z_bbZ_ATLAS13=bb_phi3_phi2Z_bbZ_TH13/ip_ex_bb_phii_phijZ_bbZ_ATLAS13(m3,m2);
    if(m2 >= 230.0 && m2 <800.0 && m3>=130.0 && m3<700.0) THoEX_bb_phi2_phi3Z_bbZ_ATLAS13=bb_phi2_phi3Z_bbZ_TH13/ip_ex_bb_phii_phijZ_bbZ_ATLAS13(m2,m3);
      
    if(mHp>= 180.0 && mHp<1000.0) THoEX_pp_Hpm_taunu_ATLAS8=pp_Hpm_taunu_TH8/ip_ex_pp_Hpm_taunu_ATLAS8(mHp);
    if(mHp>= 180.0 && mHp< 600.0) THoEX_pp_Hp_taunu_CMS8=pp_Hp_taunu_TH8/ip_ex_pp_Hp_taunu_CMS8(mHp);
    if(mHp>= 150.0 && mHp<2000.0) THoEX_pp_Hpm_taunu_ATLAS13=pp_Hpm_taunu_TH13/ip_ex_pp_Hpm_taunu_ATLAS13(mHp);
    if(mHp>= 180.0 && mHp<3000.0) THoEX_pp_Hpm_taunu_CMS13=pp_Hpm_taunu_TH13/ip_ex_pp_Hpm_taunu_CMS13(mHp);
    if(mHp>= 200.0 && mHp< 600.0) THoEX_pp_Hpm_tb_ATLAS8=pp_Hpm_tb_TH8/ip_ex_pp_Hpm_tb_ATLAS8(mHp);
    if(mHp>= 180.0 && mHp< 600.0) THoEX_pp_Hp_tb_CMS8=pp_Hp_tb_TH8/ip_ex_pp_Hp_tb_CMS8(mHp);
    if(mHp>= 200.0 && mHp<2000.0) THoEX_pp_Hpm_tb_ATLAS13=pp_Hpm_tb_TH13/ip_ex_pp_Hpm_tb_ATLAS13(mHp);

    return 0.;
}

void GeneralTHDMcache::runGeneralTHDMparameters()
{

    std::string RGEorder=myGTHDM->getRGEorderflag();
    //flag will be used to transport information about model and RGEorder to the Runner:
    //flag=0 for LO (1 for approxNLO and 2 for NLO - not implemented yet)
    int flag;
    if( RGEorder == "LO" ) flag=0;
//    else if( RGEorder == "approxNLO" ) flag=1;
//    else if( RGEorder == "NLO" ) flag=2;
    else {
        throw std::runtime_error("RGEorder can be only \"LO\" at the moment");//any of \"LO\", \"approxNLO\" or \"NLO\"
    }

    if(myGTHDM->getATHDMflag() && myGTHDM->getCPconservationflag())
    {
        double g1_at_MZ=sqrt(4.0*M_PI*Ale/cW2);
        double g2_at_MZ=sqrt(4.0*M_PI*Ale/(1-cW2));
        double g3_at_MZ=sqrt(4.0*M_PI*Als);
        double v1_at_MZ=0.;
        double v2_at_MZ=0.;
        double etaU1_at_MZ=0.;
        double etaU2_at_MZ=0.;
        double etaD1_at_MZ=0.;
        double etaD2_at_MZ=0.;
        double etaL1_at_MZ=0.;
        double etaL2_at_MZ=0.;
        double m11sq_at_MZ=0.;
        double m12sq_at_MZ=0.;
        double m22sq_at_MZ=0.;
        double lambda1_at_MZ=lambda1;
        double lambda2_at_MZ=lambda2;
        double lambda3_at_MZ=lambda3;
        double lambda4_at_MZ=lambda4;
        double lambda5_at_MZ=Relambda5;
        double lambda6_at_MZ=Relambda6;
        double lambda7_at_MZ=Relambda7;
        double RpepsGTHDM=myGTHDM->getRpepsGTHDM();
        double NLOuniscale=myGTHDM->getNLOuniscaleGTHDM();

        if(fabs(Q_GTHDM-log10(MZ))<0.005)   //at MZ scale
        {
            Q_cutoff=log10(MZ);

            g1_at_Q = g1_at_MZ;
            g2_at_Q = g2_at_MZ;
            g3_at_Q = g3_at_MZ;
            v1_at_Q = v1_at_MZ;
            v2_at_Q = v2_at_MZ;
            etaU1_at_Q = etaU1_at_MZ;
            etaU2_at_Q = etaU2_at_MZ;
            etaD1_at_Q = etaD1_at_MZ;
            etaD2_at_Q = etaD2_at_MZ;
            etaL1_at_Q = etaL1_at_MZ;
            etaL2_at_Q = etaL2_at_MZ;
            m11sq_at_Q = m11sq_at_MZ;
            m12sq_at_Q = m12sq_at_MZ;
            m22sq_at_Q = m22sq_at_MZ;
            lambda1_at_Q = lambda1_at_MZ;
            lambda2_at_Q = lambda2_at_MZ;
            lambda3_at_Q = lambda3_at_MZ;
            lambda4_at_Q = lambda4_at_MZ;
            Relambda5_at_Q = lambda5_at_MZ;
            Relambda6_at_Q = lambda6_at_MZ;
            Relambda7_at_Q = lambda7_at_MZ;
        }
        else   //at some other scale
        {
            double InitVals[21];
            InitVals[0]=g1_at_MZ;
            InitVals[1]=g2_at_MZ;
            InitVals[2]=g3_at_MZ;
            InitVals[3]=v1_at_MZ;
            InitVals[4]=v2_at_MZ;
            InitVals[5]=etaU1_at_MZ;
            InitVals[6]=etaU2_at_MZ;
            InitVals[7]=etaD1_at_MZ;
            InitVals[8]=etaD2_at_MZ;
            InitVals[9]=etaL1_at_MZ;
            InitVals[10]=etaL2_at_MZ;
            InitVals[11]=m11sq_at_MZ;
            InitVals[12]=m12sq_at_MZ;
            InitVals[13]=m22sq_at_MZ;
            InitVals[14]=lambda1_at_MZ;
            InitVals[15]=lambda2_at_MZ;
            InitVals[16]=lambda3_at_MZ;
            InitVals[17]=lambda4_at_MZ;
            InitVals[18]=lambda5_at_MZ;
            InitVals[19]=lambda6_at_MZ;
            InitVals[20]=lambda7_at_MZ;

            Q_cutoff=myRunnerGTHDM->RGEGeneralTHDMRunner(InitVals, 21, log10(MZ), Q_GTHDM, flag, RpepsGTHDM, NLOuniscale);  //Running up to Q_cutoff<=Q_GTHDM

            g1_at_Q = InitVals[0];
            g2_at_Q = InitVals[1];
            g3_at_Q = InitVals[2];
            v1_at_Q = InitVals[3];
            v2_at_Q = InitVals[4];
            etaU1_at_Q = InitVals[5];
            etaU2_at_Q = InitVals[6];
            etaD1_at_Q = InitVals[7];
            etaD2_at_Q = InitVals[8];
            etaL1_at_Q = InitVals[9];
            etaL2_at_Q = InitVals[10];
            m11sq_at_Q = InitVals[11];
            m12sq_at_Q = InitVals[12];
            m22sq_at_Q = InitVals[13];
            lambda1_at_Q = InitVals[14];
            lambda2_at_Q = InitVals[15];
            lambda3_at_Q = InitVals[16];
            lambda4_at_Q = InitVals[17];
            Relambda5_at_Q = InitVals[18];
            Relambda6_at_Q = InitVals[19];
            Relambda7_at_Q = InitVals[20];
        }
    }//End ATHDM case
    else {
        throw std::runtime_error("RGE are only defined for the CP conserving ATHDM at the moment");
    }
}

double GeneralTHDMcache::updateCache()
{      
    mHl=myGTHDM->getMHl();
    mH1sq=mHl*mHl;
    mH2sq=myGTHDM->getmH2sq();
    mH3sq=myGTHDM->getmH3sq();
    mHp2=myGTHDM->getmHp2();
    vev=myGTHDM->v();
    tanb=myGTHDM->gettanb();
    cosb=myGTHDM->getcosb();
    sinb=myGTHDM->getsinb();
    cosa1=myGTHDM->getcosalpha1();
    sina1=myGTHDM->getsinalpha1();
    cosa2=myGTHDM->getcosalpha2();
    sina2=myGTHDM->getsinalpha2();
    cosa3=myGTHDM->getcosalpha3();
    sina3=myGTHDM->getsinalpha3();
    Relambda5=myGTHDM->getRelambda5();
    Imlambda5=myGTHDM->getImlambda5();
    Relambda6=myGTHDM->getRelambda6();
    Relambda7=myGTHDM->getRelambda7();
    
    R11_GTHDM = cosa1*cosa2;
    R12_GTHDM = sina1*cosa2;
    R13_GTHDM = -sina2;
    R21_GTHDM = cosa1*sina2*sina3 - sina1*cosa3;
    R22_GTHDM = sina1*sina2*sina3 + cosa1*cosa3;
    R23_GTHDM = cosa2*sina3;
    R31_GTHDM = cosa2*sina2*cosa3 + sina1*sina3;
    R32_GTHDM = sina1*sina2*cosa3 - cosa1*sina3;
    R33_GTHDM = cosa2*cosa3;
    
    R11 = cosa1*cosa2;
    R12 = sina1*cosa2;
    R13 = -sina2;
    R21 = cosa1*sina2*sina3 - sina1*cosa3;
    R22 = sina1*sina2*sina3 + cosa1*cosa3;
    R23 = cosa2*sina3;
    R31 = cosa2*sina2*cosa3 + sina1*sina3;
    R32 = sina1*sina2*cosa3 - cosa1*sina3;
    R33 = cosa2*cosa3;

    /*The Mij_2 are defined such that Msqdiag = -2*RT*M_2*R with the rotation Matrix R
     * and Msqdiag containing the physical mass squares on the diagonal. */

    M11_2 = -0.5*(mH1sq*cosa1*cosa1*cosa2*cosa2
                  +mH2sq*sina1*sina1*cosa2*cosa2 + mH3sq*sina2*sina2);
    M12_2 = 0.5*cosa2*((mH1sq-mH2sq)*cosa1*cosa3*sina1
                       +(-mH3sq+mH1sq*cosa1*cosa1+mH2sq*sina1*sina1)*sina2*sina3);    
    M13_2 = 0.5*cosa2*(cosa3*(-mH3sq+mH1sq*cosa1*cosa1+mH2sq*sina1*sina1)*sina2
                       +(mH2sq-mH1sq)*cosa1*sina1*sina3);
    M22_2 = -0.5*(mH3sq*cosa2*cosa2*sina3*sina3
                  +mH1sq*(cosa3*sina1+cosa1*sina2*sina3)*(cosa3*sina1+cosa1*sina2*sina3)
                  +mH2sq*(cosa1*cosa3-sina1*sina2*sina3)*(cosa1*cosa3-sina1*sina2*sina3));
    M23_2 = 0.5*((mH2sq-mH1sq)*cosa1*(cosa3*cosa3-sina3*sina3)*sina1*sina2
                 +cosa1*cosa1*cosa3*(mH2sq-mH1sq*sina2*sina2)*sina3
                 -cosa3*sina3*(mH3sq*cosa2*cosa2+sina1*sina1*(-mH1sq+mH2sq*sina2*sina2)));
    M33_2 = -0.5*(mH3sq*cosa2*cosa2*cosa3*cosa3
                  +mH2sq*(cosa3*sina1*sina2+cosa1*sina3)*(cosa3*sina1*sina2+cosa1*sina3)
                  +mH1sq*(cosa1*cosa3*sina2-sina1*sina3)*(cosa1*cosa3*sina2-sina1*sina3));

    //Remaining general potential parameters
    m11sq     = M11_2 - M33_2 - M12_2*tanb + 0.5*Relambda5*vev*vev
                + (M33_2-0.5*Relambda5*vev*vev)*(cosb*cosb-sinb*sinb)
                + 0.5*vev*vev*((Relambda6-Relambda7)*sinb*cosb+Relambda7*tanb);

    m22sq     = M11_2 - M33_2 + M12_2/tanb + 0.5*Relambda5*vev*vev
                - (M33_2-0.5*Relambda5*vev*vev)*(cosb*cosb-sinb*sinb)
                + 0.25*vev*vev*(Relambda6+Relambda7+(Relambda6-Relambda7)*(cosb*cosb-sinb*sinb))/tanb;

    Rem12sq   = 0.25*vev*vev*(Relambda6+Relambda7+(Relambda6-Relambda7)*(cosb*cosb-sinb*sinb))
                - (2.0*M33_2-Relambda5*vev*vev)*sinb*cosb;

    Imm12sq   = M13_2;

    lambda1   = (-2.0*(M11_2-M22_2+M33_2) + Relambda5*vev*vev
                 - (2.0*M22_2-2.0*M33_2+Relambda5*vev*vev)/(cosb*cosb)
                 + (4.0*M12_2-2.0*Relambda6*vev*vev)*tanb)/(vev*vev);

    lambda2   = (-2.0*(M11_2-M22_2+M33_2) + Relambda5*vev*vev
                 - (2.0*M22_2-2.0*M33_2+Relambda5*vev*vev)/(sinb*sinb)
                 - (4.0*M12_2+2.0*Relambda7*vev*vev)/tanb)/(vev*vev);

    lambda3   = -(2.0*(M11_2-M22_2-M33_2-mHp2) + Relambda5*vev*vev
                  + (2.0*M12_2+Relambda6*vev*vev)/tanb
                  - (2.0*M12_2-Relambda7*vev*vev)*tanb)/(vev*vev);

    lambda4   = Relambda5 - (2.0*mHp2+4.0*M33_2)/(vev*vev);

    Imlambda6 = (2.0*M13_2-(2.0*M23_2+0.5*Imlambda5*vev*vev)*tanb)/(vev*vev);

    Imlambda7 = 2.0*M13_2/(vev*vev) + (-0.5*Imlambda5+(2.0*M23_2)/(vev*vev))/tanb;

    //Higgs potential parameters

    m11sqH     = M11_2;
    m22sqH     = M11_2-2.0*M33_2+Relambda5*vev*vev
                 +(M12_2+0.5*(Relambda6*vev*vev))/tanb
                 -(M12_2-0.5*(Relambda7*vev*vev))*tanb;
    Rem12sqH   = -M12_2;
    Imm12sqH   = M13_2;
    lambda1H   = -2.0*M11_2/(vev*vev);
    lambda2H   = -2.0*((2.0*M12_2+Relambda7*vev*vev)/tanb
                       +M11_2-4.0*M22_2+4.0*M33_2-2.0*Relambda5*vev*vev
                       +(M22_2-M33_2+0.5*Relambda5*vev*vev)/(sinb*sinb*cosb*cosb)
                       -2.0*M12_2*tanb+Relambda6*vev*vev*tanb)/(vev*vev);
    lambda3H   = -((2.0*(M11_2-2.0*M33_2-mHp2+Relambda5*vev*vev)
                    +(2.0*M12_2+Relambda6*vev*vev)/tanb
                    -(2.0*M12_2-Relambda7*vev*vev)*tanb)/(vev*vev));
    lambda4H   = -2.0*(M22_2+M33_2+mHp2)/(vev*vev);
    Relambda5H = -2.0*(M22_2-M33_2)/(vev*vev);
    Imlambda5H = 4.0*M23_2/(vev*vev);
    Relambda6H = -2.0*M12_2/(vev*vev);
    Imlambda6H = 2.0*M13_2/(vev*vev);
    Relambda7H = (-2.0*M12_2+(Relambda6-Relambda7)*vev*vev
                  +(2.0*M22_2-2.0*M33_2+Relambda5*vev*vev)*(tanb-1.0/tanb))/(vev*vev);
    Imlambda7H = 2.0*(M13_2-M23_2*(tanb-1.0/tanb))/(vev*vev)-0.5*Imlambda5/(sinb*cosb);


    M2 = Rem12sq/(sinb*cosb);

   
//    R11_GTHDM = cosalpha1*cosalpha2;
//    R12_GTHDM = sinalpha1*cosalpha2;
//    R13_GTHDM = -sinalpha2;
//    R21_GTHDM = cosalpha1*sinalpha2*sinalpha3 - sinalpha1*cosalpha3;
//    R22_GTHDM = sinalpha1*sinalpha2*sinalpha3 + cosalpha1*cosalpha3;
//    R23_GTHDM = cosalpha2*sinalpha3;
//    R31_GTHDM = cosalpha1*sinalpha2*cosalpha3 + sinalpha1*sinalpha3;
//    R32_GTHDM = sinalpha1*sinalpha2*cosalpha3 - cosalpha1*sinalpha3;
//    R33_GTHDM = cosalpha2*cosalpha3;

//    M13_2 = -vev*vev*(sinb*cosb*Imlambda5 + cosb*cosb*Imlambda6 + sinb*sinb*Imlambda7);
//    M23_2 = -vev*vev*((cosb*cosb - sinb*sinb)*Imlambda5 + 2.*sinb*cosb*(Imlambda7 - Imlambda6))/2.;

//        std::cout<<"mH1sq before ordering = "<<mH1sq<<std::endl;
//        std::cout<<"mH2sq before ordering = "<<mH2sq<<std::endl;
//        std::cout<<"mH3sq before ordering = "<<mH3sq<<std::endl;

    if(mH1sq<mH3sq && mH3sq<mH2sq)
    {
        //1<3<2 swap 2 and 3
        mHlight_2  = mH1sq;
        mHmedium_2 = mH3sq;
        mHheavy_2  = mH2sq;
    }
    else if(mH3sq<mH2sq && mH2sq<mH1sq)
    {
        //3<2<1 swap 1 and 3
        mHlight_2  = mH3sq;
        mHmedium_2 = mH2sq;
        mHheavy_2  = mH1sq;
    }
    else if(mH2sq<mH1sq && mH1sq<mH3sq)
    {
        //2<1<3 swap 1 and 2
        mHlight_2  = mH2sq;
        mHmedium_2 = mH1sq;
        mHheavy_2  = mH3sq;
    }
    else if(mH2sq<mH3sq && mH3sq<mH1sq)
    {
        //2<3<1: 3->2, 1->3, 2->1
        mHlight_2  = mH2sq;
        mHmedium_2 = mH3sq;
        mHheavy_2  = mH1sq;
    }
    else if(mH3sq<mH1sq && mH1sq<mH2sq)
    {
        //3<1<2 3->1, 1->2, 2->3
        mHlight_2  = mH3sq;
        mHmedium_2 = mH1sq;
        mHheavy_2  = mH2sq;
    }
    
    else if(mH1sq<mH2sq && mH2sq<mH3sq)
    {
        //1<2<3 ok
        mHlight_2  = mH1sq;
        mHmedium_2 = mH2sq;
        mHheavy_2  = mH3sq;
    }

//    M11_2 = (mH1sq*cosalpha1*cosalpha1*cosalpha2*cosalpha2 + mH2sq*sinalpha1*sinalpha1*cosalpha2*cosalpha2 + mH3sq*sinalpha2*sinalpha2);
//
//    M12_2 = (mH1sq*cosalpha1*cosalpha2*(cosalpha1*sinalpha2*sinalpha3 - cosalpha3*sinalpha1)
//             + mH2sq*cosalpha2*sinalpha1*(cosalpha1*cosalpha3 + sinalpha1*sinalpha2*sinalpha3)
//             - mH3sq*cosalpha2*sinalpha2*sinalpha3);
//
////    M13_2 = ?
//
//    M22_2 = (mH1sq*(cosalpha1*sinalpha2*sinalpha3 - cosalpha3*sinalpha1)*(cosalpha1*sinalpha2*sinalpha3 - cosalpha3*sinalpha1)
//             + mH2sq*(cosalpha1*cosalpha3 + sinalpha1*sinalpha2*sinalpha3)*(cosalpha1*cosalpha3 + sinalpha1*sinalpha2*sinalpha3)
//             + mH3sq*cosalpha2*cosalpha2*sinalpha3*sinalpha3);
//
////    M23_2 = ?
//
//    M33_2 = (mH1sq*(cosalpha1*cosalpha3*sinalpha2 + sinalpha1*sinalpha3)*(cosalpha1*cosalpha3*sinalpha2 + sinalpha1*sinalpha3)
//             + mH2sq*(cosalpha3*sinalpha1*sinalpha2 - cosalpha1*sinalpha3)*(cosalpha3*sinalpha1*sinalpha2 - cosalpha1*sinalpha3)
//             + mH3sq*cosalpha2*cosalpha2*cosalpha3*cosalpha3);
//
//    m11_2_GTHDM = M2_GTHDM*(1. - cosb*cosb + 3.*sinb*sinb)/4. + (M12_2*tanb - M11_2)/2.;
//    m22_2_GTHDM = M2_GTHDM*(1. + 3.*cosb*cosb - sinb*sinb)/4. - (M12_2/tanb + M11_2)/2.;
//    Imm12_2_GTHDM = 0.5*(cosb*sinb*Imlambda5 + cosb*cosb*Imlambda6 + sinb*sinb*Imlambda7)*vev*vev;
//    lambda1_GTHDM = (M11_2 + tanb*tanb*(M22_2-M2_GTHDM) - 2.0*tanb*M12_2)/(vev*vev) + tanb*(tanb*tanb*Relambda7 - 3.0*Relambda6)/2.0;
//    lambda2_GTHDM = (M11_2 + (M22_2-M2_GTHDM)/(tanb*tanb) + 2.0*M12_2/tanb)/(vev*vev) + (0.5*Relambda6/(tanb*tanb) - 1.5*Relambda7)/tanb;
//    lambda3_GTHDM = (M11_2 - M22_2 - M2_GTHDM + (1.0/tanb - tanb)*M12_2 + 2.0*mHp2)/(vev*vev) - (Relambda6/tanb + tanb*Relambda7)/2.0;
//    lambda4_GTHDM = (M2_GTHDM + M33_2 - 2.0*mHp2)/(vev*vev) - 0.5*(Relambda6/tanb + tanb*Relambda7);
//    Relambda5_GTHDM = (M2_GTHDM - M33_2)/(vev*vev) - 0.5*(Relambda6/tanb + tanb*Relambda7);
//    
    

    Mu_GTHDM.assign(0,0, myGTHDM->getQuarks(QCD::UP).getMass());
    Mu_GTHDM.assign(1,1, myGTHDM->getQuarks(QCD::CHARM).getMass());
    Mu_GTHDM.assign(2,2, myGTHDM->getQuarks(QCD::TOP).getMass());
    
    Md_GTHDM.assign(0,0, myGTHDM->getQuarks(QCD::DOWN).getMass());
    Md_GTHDM.assign(1,1, myGTHDM->getQuarks(QCD::STRANGE).getMass());
    Md_GTHDM.assign(2,2, myGTHDM->getQuarks(QCD::BOTTOM).getMass());
    
    Ml_GTHDM.assign(0,0, myGTHDM->getLeptons(StandardModel::ELECTRON).getMass());
    Ml_GTHDM.assign(1,1, myGTHDM->getLeptons(StandardModel::MU).getMass());
    Ml_GTHDM.assign(2,2, myGTHDM->getLeptons(StandardModel::TAU).getMass());
    
    if(myGTHDM->getATHDMflag() == true)
    {
        sigmau_ATHDM = myGTHDM->getNu_11()/myGTHDM->getQuarks(StandardModel::TOP).getMass();
        sigmad_ATHDM = myGTHDM->getNd_11()/myGTHDM->getQuarks(StandardModel::DOWN).getMass();
        sigmal_ATHDM = myGTHDM->getNl_11()/myGTHDM->getLeptons(StandardModel::TAU).getMass();
        
        Nu_GTHDM.assign(0,0, sigmau_ATHDM*Mu_GTHDM(0,0));
        Nu_GTHDM.assign(1,1, sigmau_ATHDM*Mu_GTHDM(1,1));
        Nu_GTHDM.assign(2,2, sigmau_ATHDM*Mu_GTHDM(2,2));
        
        Nd_GTHDM.assign(0,0, sigmad_ATHDM*Mu_GTHDM(0,0));
        Nd_GTHDM.assign(1,1, sigmad_ATHDM*Md_GTHDM(1,1));
        Nd_GTHDM.assign(2,2, sigmad_ATHDM*Mu_GTHDM(2,2));
        
        Nl_GTHDM.assign(0,0, sigmal_ATHDM*Ml_GTHDM(0,0));
        Nl_GTHDM.assign(1,1, sigmal_ATHDM*Ml_GTHDM(1,1));
        Nl_GTHDM.assign(2,2, sigmal_ATHDM*Ml_GTHDM(2,2));
    }
    else
    {
        Nu_GTHDM.assign(0,0, myGTHDM->getNu_11());
        Nu_GTHDM.assign(0,1, myGTHDM->getNu_12());
        Nu_GTHDM.assign(0,2, myGTHDM->getNu_13());
        Nu_GTHDM.assign(1,0, myGTHDM->getNu_21());
        Nu_GTHDM.assign(1,1, myGTHDM->getNu_22());
        Nu_GTHDM.assign(1,2, myGTHDM->getNu_23());
        Nu_GTHDM.assign(2,0, myGTHDM->getNu_31());
        Nu_GTHDM.assign(2,1, myGTHDM->getNu_32());
        Nu_GTHDM.assign(2,2, myGTHDM->getNu_33());
    
        Nd_GTHDM.assign(0,0, myGTHDM->getNd_11());
        Nd_GTHDM.assign(0,1, myGTHDM->getNd_12());
        Nd_GTHDM.assign(0,2, myGTHDM->getNd_13());
        Nd_GTHDM.assign(1,0, myGTHDM->getNd_21());
        Nd_GTHDM.assign(1,1, myGTHDM->getNd_22());
        Nd_GTHDM.assign(1,2, myGTHDM->getNd_23());
        Nd_GTHDM.assign(2,0, myGTHDM->getNd_31());
        Nd_GTHDM.assign(2,1, myGTHDM->getNd_32());
        Nd_GTHDM.assign(2,2, myGTHDM->getNd_33());
    
        Nl_GTHDM.assign(0,0, myGTHDM->getNl_11());
        Nl_GTHDM.assign(0,1, myGTHDM->getNl_12());
        Nl_GTHDM.assign(0,2, myGTHDM->getNl_13());
        Nl_GTHDM.assign(1,0, myGTHDM->getNl_21());
        Nl_GTHDM.assign(1,1, myGTHDM->getNl_22());
        Nl_GTHDM.assign(1,2, myGTHDM->getNl_23());
        Nl_GTHDM.assign(2,0, myGTHDM->getNl_31());
        Nl_GTHDM.assign(2,1, myGTHDM->getNl_32());
        Nl_GTHDM.assign(2,2, myGTHDM->getNl_33());
    }
       
    //Definition of Yukawa matrices
    Yu1_GTHDM = (cosb*Mu_GTHDM - sinb*Nu_GTHDM)*sqrt(2.)/vev;
    Yu2_GTHDM = (cosb*Nu_GTHDM + sinb*Mu_GTHDM)*sqrt(2.)/vev;
    Yd1_GTHDM = (cosb*Md_GTHDM - sinb*Nd_GTHDM)*sqrt(2.)/vev;
    Yd2_GTHDM = (cosb*Nd_GTHDM + sinb*Md_GTHDM)*sqrt(2.)/vev;
    Yl1_GTHDM = (cosb*Ml_GTHDM - sinb*Nl_GTHDM)*sqrt(2.)/vev;
    Yl2_GTHDM = (cosb*Nl_GTHDM + sinb*Ml_GTHDM)*sqrt(2.)/vev;
    
    /*up, down and leptonic couplings */
    su = myGTHDM->getNu_11();
    sd = myGTHDM->getNd_11();
    sl = myGTHDM->getNl_11();
    
  /*  std::cout << "su = " << su << std::endl;
   std::cout << "yu1R_GTHDM = " << myGTHDM->getyu1() << std::endl;
    std::cout << "cosa1 = " << cosa1 << std::endl;
    std::cout << "sina1 = " << sina1 << std::endl;*/

    
    
        
    //If to always set 1 as the SM Higgs and 3 as the heaviest
    /*This loop is different as the previous. Before it was ordered by mass. Now we 
    want 1 as the SM Higgs and 3 as the heaviest*/
   /* if(mH1sq == mHl*mHl)
    {
         m1_2 = mH1sq;
          
         R11 = R11_GTHDM;
         R12 = R12_GTHDM;
         R13 = R13_GTHDM;
         
        if(mH2sq<mH3sq)
        {   
            m2_2 = mH2sq;
            m3_2 = mH3sq;
            R21 = R21_GTHDM;
            R22 = R22_GTHDM;
            R23 = R23_GTHDM;
            R31 = R31_GTHDM;
            R32 = R32_GTHDM;
            R33 = R33_GTHDM;
        }
         else if(mH3sq<mH2sq)
        {   
            m2_2 = mH3sq;
            m3_2 = mH2sq;
            R21 = R31_GTHDM;
            R22 = R32_GTHDM;
            R23 = R33_GTHDM;
            R31 = R21_GTHDM;
            R32 = R22_GTHDM;
            R33 = R23_GTHDM;
        }
    }
    else if(mH2sq == mHl*mHl )
    {
         m1_2 = mH2sq;
         R11 = R21_GTHDM;
         R12 = R22_GTHDM;
         R13 = R23_GTHDM;
         
         if(mH1sq<mH3sq)
        {   
            m2_2 = mH1sq;
            m3_2 = mH3sq;
            R21 = R11_GTHDM;
            R22 = R12_GTHDM;
            R23 = R13_GTHDM;
            R31 = R31_GTHDM;
            R32 = R32_GTHDM;
            R33 = R33_GTHDM;

        }
         else if(mH3sq<mH1sq)
        {   
            m2_2 = mH3sq;
            m3_2 = mH1sq;
            R21 = R31_GTHDM;
            R22 = R32_GTHDM;
            R23 = R33_GTHDM;
            R31 = R11_GTHDM;
            R32 = R12_GTHDM;
            R33 = R13_GTHDM;
        }
  
    }
    else if(mH3sq == mHl*mHl )
    {
         m1_2 = mH3sq;

         R11 = R31_GTHDM;
         R12 = R32_GTHDM;
         R13 = R33_GTHDM;
         
          if(mH1sq<mH2sq)
        {   
            m2_2 = mH1sq;
            m3_2 = mH2sq;
            R21 = R11_GTHDM;
            R22 = R12_GTHDM;
            R23 = R13_GTHDM;
            R31 = R21_GTHDM;
            R32 = R22_GTHDM;
            R33 = R23_GTHDM;

        }
         else if(mH2sq<mH1sq)
        {   
            m2_2 = mH2sq;
            m3_2 = mH1sq;
            R21 = R21_GTHDM;
            R22 = R22_GTHDM;
            R23 = R23_GTHDM;
            R31 = R11_GTHDM;
            R32 = R12_GTHDM;
            R33 = R13_GTHDM;
        }
         
    }
    
      if (m1_2 < 0 || m2_2 < 0 || m3_2 < 0) 
                return std::numeric_limits<double>::quiet_NaN();
     
     */
    
    
    
    m1_2=mH1sq;
    
    
    R11 = R11_GTHDM;
    R12 = R12_GTHDM;
    R13 = R13_GTHDM;
    
   /* COMMENTED ON 08/01/2019- TEST
    *  if(mH2sq<=mH3sq)
    {   
        
            m2_2 = mH2sq;
            m3_2 = mH3sq;
            R21 = R21_GTHDM;
            R22 = R22_GTHDM;
            R23 = R23_GTHDM;
            R31 = R31_GTHDM;
            R32 = R32_GTHDM;
            R33 = R33_GTHDM;
        }
         else
        {   
            m2_2 = mH3sq;
            m3_2 = mH2sq;
            R21 = R31_GTHDM;
            R22 = R32_GTHDM;
            R23 = R33_GTHDM;THoEX_gg_phi3_bb_CMS8
            R31 = R21_GTHDM;
            R32 = R22_GTHDM;
            R33 = R23_GTHDM;
        }
    */
    if (m1_2 < 0 || m2_2 < 0 || m3_2 < 0) 
                return std::numeric_limits<double>::quiet_NaN();
    

        
    Q_GTHDM=myGTHDM->getQ_GTHDM();
    Ale=myGTHDM->getAle();
    Als=myGTHDM->getAlsMz();
//    Mt=myTHDM->getQuarks(QCD::TOP).getMass();
//    Mb=myTHDM->getQuarks(QCD::BOTTOM).getMass();   
//    Mtau=myTHDM->getLeptons(StandardModel::TAU).getMass();
//    Mc=myTHDM->getQuarks(QCD::CHARM).getMass();
//    Ms=myTHDM->getQuarks(QCD::STRANGE).getMass();
//    Mmu=myTHDM->getLeptons(StandardModel::MU).getMass();
//    Mu=myTHDM->getQuarks(QCD::UP).getMass();
//    Md=myTHDM->getQuarks(QCD::DOWN).getMass();
//    Me=myTHDM->getLeptons(StandardModel::ELECTRON).getMass();
    MZ=myGTHDM->getMz();
//    modelflag=myTHDM->getModelTypeflag();
//    WFRflag=myTHDM->getWFRflag();

    runGeneralTHDMparameters();
    computeSignalStrengths();
    computephi2quantities();
    computephi3quantities();
    ComputeHeavyHiggs();
    computeHpquantities();
 
    return mH1sq;
}
