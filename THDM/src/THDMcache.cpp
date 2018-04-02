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

#include "log_cs_ggH_13.h"

THDMcache::THDMcache(const StandardModel& SM_i)

        :br_tt(19961, 2, 0.),
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
        ATLAS8_pp_phi_gaga(108, 2, 0.),
        ATLAS8_pp_phi_Zga_llga(141, 2, 0.),
        ATLAS8_gg_phi_tautau(92, 2, 0.),
        ATLAS8_bb_phi_tautau(92, 2, 0.),
        ATLAS8_gg_A_hZ_tautauZ(79, 2, 0.),
        ATLAS8_gg_A_hZ_bbZ(79, 2, 0.),
        ATLAS8_gg_phi_tt(53, 2, 0.),
        ATLAS8_gg_H_WW(13,2,0.),
        ATLAS8_VBF_H_WW(13,2,0.),
        ATLAS8_gg_H_ZZ(173,2,0.),
        ATLAS8_VBF_H_ZZ(173,2,0.),
        ATLAS8_gg_H_hh(75,2,0.),
        ATLAS8_pp_phi_gaga_e(108, 2, 0.),
        ATLAS8_pp_phi_Zga_llga_e(141, 2, 0.),
        ATLAS8_gg_phi_tautau_e(92, 2, 0.),
        ATLAS8_bb_phi_tautau_e(92, 2, 0.),
        ATLAS8_gg_A_hZ_tautauZ_e(79, 2, 0.),
        ATLAS8_gg_A_hZ_bbZ_e(79, 2, 0.),
        ATLAS8_gg_phi_tt_e(53, 2, 0.),
        ATLAS8_gg_H_WW_e(13,2,0.),
        ATLAS8_VBF_H_WW_e(13,2,0.),
        ATLAS8_gg_H_ZZ_e(173,2,0.),
        ATLAS8_VBF_H_ZZ_e(173,2,0.),
        ATLAS8_gg_H_hh_e(75,2,0.),
        CMS8_mu_pp_H_VV(172, 2, 0.),
        CMS8_mu_pp_H_VV_e(172, 2, 0.),
        CMS8_gg_A_hZ_bbll(16, 2, 0.),
        CMS8_pp_H_hh(71, 2, 0.),
        CMS8_pp_H_hh_gagabb(85, 2, 0.),
        CMS8_pp_H_hh_bbbb(167, 2, 0.),
        CMS8_bb_phi_bb(81, 2, 0.),
        CMS8_gg_phi_tautau(92,2,0.),
        CMS8_bb_phi_tautau(92,2,0.),
        CMS8_gg_phi_gaga(141,2,0.),
        CMS8_pp_A_Zga_llga(101,2,0.),
        CMS8_pp_phi_Zga(101,2,0.),
        CMS8_gg_H_hh_bbtautau(10,2,0.),
        CMS8_gg_A_hZ_tautaull(14,2,0.),
        CMS8_pp_A_HZ_bbll(28718, 3, 0.),
        CMS8_pp_H_AZ_bbll(29050, 3, 0.),
        CMS8_pp_A_HZ_tautaull(400, 3, 0.),
        CMS8_pp_H_AZ_tautaull(400, 3, 0.),
        CMS8_gg_A_hZ_bbll_e(16, 2, 0.),
        CMS8_pp_H_hh_e(71, 2, 0.),
        CMS8_pp_H_hh_gagabb_e(85, 2, 0.),
        CMS8_pp_H_hh_bbbb_e(167, 2, 0.),
        CMS8_bb_phi_bb_e(81, 2, 0.),
        CMS8_gg_phi_tautau_e(92,2,0.),
        CMS8_bb_phi_tautau_e(92,2,0.),
        CMS8_gg_phi_gaga_e(141,2,0.),
        CMS8_pp_A_Zga_llga_e(101,2,0.),
        CMS8_gg_H_hh_bbtautau_e(10,2,0.),
        CMS8_gg_A_hZ_tautaull_e(14,2,0.),
//        CMS_ggF_phi_gaga_ep1(141,2,0.),
//        CMS_gg_phi_gaga_ep2(141,2,0.),
//        CMS_ggF_phi_gaga_em1(141,2,0.),
//        CMS_ggF_phi_gaga_em2(141,2,0.),
        ATLAS13_bb_phi_tt(61,2,0.),
        ATLAS13_tt_phi_tt(61,2,0.),
        ATLAS13_gg_phi_tautau(206,2,0.),
        ATLAS13_bb_phi_tautau(206,2,0.),
        ATLAS13_pp_phi_gaga(251,2,0.),
        ATLAS13_pp_phi_Zga(216,2,0.),
        ATLAS13_gg_phi_Zga_llga(216,2,0.),
        ATLAS13_gg_H_ZZ_llllnunu(101,2,0.),
        ATLAS13_VBF_H_ZZ_llllnunu(101,2,0.),
        ATLAS13_gg_H_ZZ_llnunu(71,2,0.),
        ATLAS13_gg_H_ZZ_llll(81,2,0.),
        ATLAS13_VBF_H_ZZ_llll(81,2,0.),
        ATLAS13_gg_H_ZZ_qqllnunu(271,2,0.),
        ATLAS13_VBF_H_ZZ_qqllnunu(271,2,0.),
        ATLAS13_gg_H_ZZ_llqq(271,2,0.),
        ATLAS13_VBF_H_ZZ_llqq(271,2,0.),
        ATLAS13_gg_H_ZZ_nunuqq(251,2,0.),
        ATLAS13_gg_H_WW_enumumu(381,2,0.),
        ATLAS13_VBF_H_WW_enumumu(281,2,0.),
        ATLAS13_gg_H_WW_lnuqq(271,2,0.),
        ATLAS13_VBF_H_WW_lnuqq(271,2,0.),
        ATLAS13_pp_H_VV_qqqq(181,2,0.),
        ATLAS13_pp_H_hh_bbbb(271,2,0.),
        ATLAS13_pp_H_hh_gagabb(26,2,0.),
        ATLAS13_pp_H_hh_gagaWW(25,2,0.),
        ATLAS13_gg_A_Zh_Zbb(181,2,0.),
        ATLAS13_bb_A_Zh_Zbb(181,2,0.),
        ATLAS13_bb_phi_tt_e(61,2,0.),
        ATLAS13_tt_phi_tt_e(61,2,0.),
        ATLAS13_gg_phi_tautau_e(206,2,0.),
        ATLAS13_bb_phi_tautau_e(206,2,0.),
        ATLAS13_pp_phi_gaga_e(251,2,0.),
        ATLAS13_pp_phi_Zga_e(216,2,0.),
        ATLAS13_gg_phi_Zga_llga_e(216,2,0.),
        ATLAS13_gg_H_ZZ_llllnunu_e(101,2,0.),
        ATLAS13_VBF_H_ZZ_llllnunu_e(101,2,0.),
        ATLAS13_gg_H_ZZ_llnunu_e(71,2,0.),
        ATLAS13_gg_H_ZZ_llll_e(81,2,0.),
        ATLAS13_VBF_H_ZZ_llll_e(81,2,0.),
        ATLAS13_gg_H_ZZ_qqllnunu_e(271,2,0.),
        ATLAS13_VBF_H_ZZ_qqllnunu_e(271,2,0.),
        ATLAS13_gg_H_ZZ_llqq_e(271,2,0.),
        ATLAS13_VBF_H_ZZ_llqq_e(271,2,0.),
        ATLAS13_gg_H_ZZ_nunuqq_e(251,2,0.),
        ATLAS13_gg_H_WW_enumumu_e(381,2,0.),
        ATLAS13_VBF_H_WW_enumumu_e(281,2,0.),
        ATLAS13_gg_H_WW_lnuqq_e(271,2,0.),
        ATLAS13_VBF_H_WW_lnuqq_e(271,2,0.),
        ATLAS13_pp_H_VV_qqqq_e(181,2,0.),
        ATLAS13_pp_H_hh_bbbb_e(271,2,0.),
        ATLAS13_pp_H_hh_gagabb_e(26,2,0.),
        ATLAS13_pp_H_hh_gagaWW_e(25,2,0.),
        ATLAS13_gg_A_Zh_Zbb_e(181,2,0.),
        ATLAS13_bb_A_Zh_Zbb_e(181,2,0.),
        CMS13_pp_phi_bb(66,2,0.),
        CMS13_gg_phi_tautau(312,2,0.),
        CMS13_bb_phi_tautau(312,2,0.),
        CMS13_gg_phi_gaga(351,2,0.),
        CMS13_pp_phi_Zga_llga(171,2,0.),
        CMS13_pp_phi_Zga_qqga(236,2,0.),
        CMS13_ggF_phi_Zga(366,2,0.),
        CMS13_pp_H_ZZ_llnunu(191,2,0.),
        CMS13_gg_H_ZZ_llnunu(131,2,0.),
        CMS13_VBF_H_ZZ_llnunu(131,2,0.),
        CMS13_pp_H_ZZ_llll(241,2,0.),
        CMS13_VBFVH_H_ZZ_llll(241,2,0.),
        CMS13_pp_H_ZZ_llqq(151,2,0.),
        CMS13_ggFVBF_H_WW_lnulnu(81,2,0.),
        CMS13_pp_H_hh_bbbb(95,2,0.),
        CMS13_ggF_H_hh_bbbb(226,2,0.),
        CMS13_pp_H_hh_gagabb(66,2,0.),
        CMS13_pp_H_hh_bbtautau(66,2,0.),
        CMS13_pp_H_hh_bbtautau1(66,2,0.),
        CMS13_pp_H_hh_bblnulnu(65,2,0.),
        CMS13_pp_H_hh_bbVV(65,2,0.),
        CMS13_pp_phi_bb_e(66,2,0.),
        CMS13_gg_phi_tautau_e(312,2,0.),
        CMS13_bb_phi_tautau_e(312,2,0.),
        CMS13_gg_phi_gaga_e(351,2,0.),
        CMS13_pp_phi_Zga_llga_e(171,2,0.),
        CMS13_pp_phi_Zga_qqga_e(236,2,0.),
        CMS13_ggF_phi_Zga_e(366,2,0.),
        CMS13_pp_H_ZZ_llnunu_e(191,2,0.),
        CMS13_gg_H_ZZ_llnunu_e(131,2,0.),
        CMS13_VBF_H_ZZ_llnunu_e(131,2,0.),
        CMS13_pp_H_ZZ_llll_e(241,2,0.),
        CMS13_VBFVH_H_ZZ_llll_e(241,2,0.),
        CMS13_pp_H_ZZ_llqq_e(151,2,0.),
        CMS13_ggFVBF_H_WW_lnulnu_e(81,2,0.),
        CMS13_pp_H_hh_bbbb_e(95,2,0.),
        CMS13_ggF_H_hh_bbbb_e(226,2,0.),
        CMS13_pp_H_hh_gagabb_e(66,2,0.),
        CMS13_pp_H_hh_bbtautau_e(66,2,0.),
        CMS13_pp_H_hh_bbtautau1_e(66,2,0.),
        CMS13_pp_H_hh_bblnulnu_e(65,2,0.),
        CMS13_pp_H_hh_bbVV_e(65,2,0.),
        temp1(1,2,0.), temp2(1,2,0.), temp3(1,2,0.), temp4(1,2,0.), temp5(1,2,0.), temp6(1,2,0.), temp7(1,2,0.), temp8(1,2,0.), temp9(1,2,0.), temp10(1,2,0.),
        temp11(1,2,0.), temp12(1,2,0.), temp13(1,2,0.), temp14(1,2,0.), temp15(1,2,0.), temp16(1,2,0.), temp17(1,2,0.), temp18(1,2,0.), temp19(1,2,0.), temp20(1,2,0.),
        temp21(1,2,0.), temp22(1,2,0.), temp23(1,2,0.), temp24(1,2,0.), temp25(1,2,0.), temp26(1,2,0.), temp27(1,2,0.), temp28(1,2,0.), temp29(1,2,0.), temp30(1,2,0.),
        temp31(1,2,0.), temp32(1,2,0.), temp33(1,2,0.), temp34(1,2,0.), temp35(1,2,0.), temp36(1,2,0.), temp37(1,2,0.), temp38(1,2,0.), temp39(1,2,0.), temp40(1,2,0.),
        temp1e(1,2,0.), temp2e(1,2,0.), temp3e(1,2,0.), temp4e(1,2,0.), temp5e(1,2,0.), temp6e(1,2,0.), temp7e(1,2,0.), temp8e(1,2,0.), temp9e(1,2,0.), temp10e(1,2,0.),
        temp11e(1,2,0.), temp12e(1,2,0.), temp13e(1,2,0.), temp14e(1,2,0.), temp15e(1,2,0.), temp16e(1,2,0.), temp17e(1,2,0.), temp18e(1,2,0.), temp19e(1,2,0.), temp20e(1,2,0.),
        temp21e(1,2,0.), temp22e(1,2,0.), temp23e(1,2,0.), temp24e(1,2,0.), temp25e(1,2,0.), temp26e(1,2,0.), temp27e(1,2,0.), temp28e(1,2,0.), temp29e(1,2,0.), temp30e(1,2,0.),
        temp31e(1,2,0.), temp32e(1,2,0.), temp33e(1,2,0.), temp34e(1,2,0.), temp35e(1,2,0.), temp36e(1,2,0.), temp37e(1,2,0.), temp38e(1,2,0.), temp39e(1,2,0.), temp40e(1,2,0.),
        ATLAS8_pp_Hpm_taunu(83,2,0.),
        ATLAS8_pp_Hpm_tb(41,2,0.),
        ATLAS8_pp_Hpm_taunu_e(83,2,0.),
        ATLAS8_pp_Hpm_tb_e(41,2,0.),
        CMS8_pp_Hp_taunu(43,2,0.),
        CMS8_pp_Hp_tb(43,2,0.),
        CMS8_pp_Hp_taunu_e(43,2,0.),
        CMS8_pp_Hp_tb_e(43,2,0.),
        ATLAS13_pp_Hpm_taunu(181,2,0.),
        ATLAS13_pp_Hp_tb1(71,2,0.),
        ATLAS13_pp_Hp_tb2(181,2,0.),
        ATLAS13_pp_Hpm_taunu_e(181,2,0.),
        ATLAS13_pp_Hp_tb1_e(71,2,0.),
        ATLAS13_pp_Hp_tb2_e(181,2,0.),
        CMS13_pp_Hpm_taunu(283,2,0.),
        CMS13_pp_Hpm_taunu_e(283,2,0.),
        arraybsgamma(1111, 3, 0.),
        myTHDM(static_cast<const THDM*> (&SM_i)),
        PV(true)
{
    mym11_2=new m11_2(SM_i);
    mym22_2=new m22_2(SM_i);
    mylambda1=new lambda1(SM_i);
    mylambda2=new lambda2(SM_i);
    mylambda3=new lambda3(SM_i);
    mylambda4=new lambda4(SM_i);
    mylambda5=new lambda5(SM_i);
    myRunner=new Runner(SM_i);
  read();
}

THDMcache::~THDMcache()
{
  delete mym11_2;
  delete mym22_2;
  delete mylambda1;
  delete mylambda2;
  delete mylambda3;
  delete mylambda4;
  delete mylambda5;
  delete myRunner;
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

int THDMcache::CacheCheckReal(const double cache[][CacheSize], 
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

void THDMcache::CacheShiftReal(double cache[][CacheSize], const int NumPar,
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

gslpp::complex THDMcache::B0_MZ2_0_0_mHl2(const double MZ2, const double mHl2) const {
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

gslpp::complex THDMcache::B0_MZ2_0_0_mHh2(const double MZ2, const double mHh2) const {
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

gslpp::complex THDMcache::B0_MZ2_0_mHp2_mHl2(const double MZ2, const double mHp2, const double mHl2) const {
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

gslpp::complex THDMcache::B0_MZ2_0_mHp2_mHh2(const double MZ2, const double mHp2, const double mHh2) const {
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

gslpp::complex THDMcache::B0_MZ2_0_mA2_mHl2(const double MZ2, const double mA2, const double mHl2) const {
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

gslpp::complex THDMcache::B0_MZ2_0_mA2_mHh2(const double MZ2, const double mA2, const double mHh2) const {
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

gslpp::complex THDMcache::B0_MZ2_mHl2_0_0(const double MZ2, const double mHl2) const {
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

gslpp::complex THDMcache::B0_MZ2_mHl2_0_mHp2(const double MZ2, const double mHl2, const double mHp2) const {
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

gslpp::complex THDMcache::B0_MZ2_mHl2_0_mA2(const double MZ2, const double mHl2, const double mA2) const {
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

gslpp::complex THDMcache::B0_MZ2_mHl2_mHl2_mHl2(const double MZ2, const double mHl2) const {
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

gslpp::complex THDMcache::B0_MZ2_mHl2_mHh2_mHl2(const double MZ2, const double mHl2, const double mHh2) const {
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

gslpp::complex THDMcache::B0_MZ2_mHl2_mHh2_mHh2(const double MZ2, const double mHl2, const double mHh2) const {
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

gslpp::complex THDMcache::B0_MZ2_mHl2_mHp2_mHp2(const double MZ2, const double mHl2, const double mHp2) const {
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

gslpp::complex THDMcache::B0_MZ2_mHl2_mA2_mA2(const double MZ2, const double mHl2, const double mA2) const {
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

gslpp::complex THDMcache::B0_MZ2_mHh2_0_0(const double MZ2, const double mHh2) const {
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

gslpp::complex THDMcache::B0_MZ2_mHh2_0_mHp2(const double MZ2, const double mHh2, const double mHp2) const {
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

gslpp::complex THDMcache::B0_MZ2_mHh2_0_mA2(const double MZ2, const double mHh2, const double mA2) const {
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

gslpp::complex THDMcache::B0_MZ2_mHh2_mHl2_mHl2(const double MZ2, const double mHh2, const double mHl2) const {
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

gslpp::complex THDMcache::B0_MZ2_mHh2_mHh2_mHl2(const double MZ2, const double mHh2, const double mHl2) const {
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

gslpp::complex THDMcache::B0_MZ2_mHh2_mHh2_mHh2(const double MZ2, const double mHh2) const {
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

gslpp::complex THDMcache::B0_MZ2_mHh2_mHp2_mHp2(const double MZ2, const double mHh2, const double mHp2) const {
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

gslpp::complex THDMcache::B0_MZ2_mHh2_mA2_mA2(const double MZ2, const double mHh2, const double mA2) const {
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

gslpp::complex THDMcache::B0_MZ2_mHp2_0_mHl2(const double MZ2, const double mHp2, const double mHl2) const {
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

gslpp::complex THDMcache::B0_MZ2_mHp2_0_mHh2(const double MZ2, const double mHp2, const double mHh2) const {
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

gslpp::complex THDMcache::B0_MZ2_mHp2_mHp2_mHl2(const double MZ2, const double mHp2, const double mHl2) const {
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

gslpp::complex THDMcache::B0_MZ2_mHp2_mHp2_mHh2(const double MZ2, const double mHp2, const double mHh2) const {
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

gslpp::complex THDMcache::B0_MZ2_mA2_0_mHl2(const double MZ2, const double mA2, const double mHl2) const {
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

gslpp::complex THDMcache::B0_MZ2_mA2_0_mHh2(const double MZ2, const double mA2, const double mHh2) const {
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

gslpp::complex THDMcache::B0_MZ2_mA2_mA2_mHl2(const double MZ2, const double mA2, const double mHl2) const {
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

gslpp::complex THDMcache::B0_MZ2_mA2_mA2_mHh2(const double MZ2, const double mA2, const double mHh2) const {
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

gslpp::complex THDMcache::B0p_MZ2_0_0_mHl2(const double MZ2, const double mHl2) const {
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

gslpp::complex THDMcache::B0p_MZ2_0_0_mHh2(const double MZ2, const double mHh2) const {
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

gslpp::complex THDMcache::B0p_MZ2_0_mHp2_mHl2(const double MZ2, const double mHp2, const double mHl2) const {
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

gslpp::complex THDMcache::B0p_MZ2_0_mHp2_mHh2(const double MZ2, const double mHp2, const double mHh2) const {
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

gslpp::complex THDMcache::B0p_MZ2_0_mHp2_mA2(const double MZ2, const double mHp2, const double mA2) const {
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

gslpp::complex THDMcache::B0p_MZ2_0_mA2_mHl2(const double MZ2, const double mA2, const double mHl2) const {
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

gslpp::complex THDMcache::B0p_MZ2_0_mA2_mHh2(const double MZ2, const double mA2, const double mHh2) const {
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

gslpp::complex THDMcache::B0p_MZ2_mHl2_0_0(const double MZ2, const double mHl2) const {
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

gslpp::complex THDMcache::B0p_MZ2_mHl2_0_mHp2(const double MZ2, const double mHl2, const double mHp2) const {
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

gslpp::complex THDMcache::B0p_MZ2_mHl2_0_mA2(const double MZ2, const double mHl2, const double mA2) const {
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

gslpp::complex THDMcache::B0p_MZ2_mHl2_mHl2_mHl2(const double MZ2, const double mHl2) const {
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

gslpp::complex THDMcache::B0p_MZ2_mHl2_mHh2_mHl2(const double MZ2, const double mHl2, const double mHh2) const {
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

gslpp::complex THDMcache::B0p_MZ2_mHl2_mHh2_mHh2(const double MZ2, const double mHl2, const double mHh2) const {
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

gslpp::complex THDMcache::B0p_MZ2_mHl2_mHp2_mHp2(const double MZ2, const double mHl2, const double mHp2) const {
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

gslpp::complex THDMcache::B0p_MZ2_mHl2_mA2_mA2(const double MZ2, const double mHl2, const double mA2) const {
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

gslpp::complex THDMcache::B0p_MZ2_mHh2_0_0(const double MZ2, const double mHh2) const {
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

gslpp::complex THDMcache::B0p_MZ2_mHh2_0_mHp2(const double MZ2, const double mHh2, const double mHp2) const {
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

gslpp::complex THDMcache::B0p_MZ2_mHh2_0_mA2(const double MZ2, const double mHh2, const double mA2) const {
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

gslpp::complex THDMcache::B0p_MZ2_mHh2_mHl2_mHl2(const double MZ2, const double mHh2, const double mHl2) const {
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

gslpp::complex THDMcache::B0p_MZ2_mHh2_mHh2_mHl2(const double MZ2, const double mHh2, const double mHl2) const {
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

gslpp::complex THDMcache::B0p_MZ2_mHh2_mHh2_mHh2(const double MZ2, const double mHh2) const {
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

gslpp::complex THDMcache::B0p_MZ2_mHh2_mHp2_mHp2(const double MZ2, const double mHh2, const double mHp2) const {
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

gslpp::complex THDMcache::B0p_MZ2_mHh2_mA2_mA2(const double MZ2, const double mHh2, const double mA2) const {
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

gslpp::complex THDMcache::B0p_MZ2_mHp2_0_mHl2(const double MZ2, const double mHp2, const double mHl2) const {
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

gslpp::complex THDMcache::B0p_MZ2_mHp2_0_mHh2(const double MZ2, const double mHp2, const double mHh2) const {
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

gslpp::complex THDMcache::B0p_MZ2_mHp2_0_mA2(const double MZ2, const double mHp2, const double mA2) const {
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

gslpp::complex THDMcache::B0p_MZ2_mHp2_mHp2_mHl2(const double MZ2, const double mHp2, const double mHl2) const {
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

gslpp::complex THDMcache::B0p_MZ2_mHp2_mHp2_mHh2(const double MZ2, const double mHp2, const double mHh2) const {
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

gslpp::complex THDMcache::B0p_MZ2_mA2_0_mHl2(const double MZ2, const double mA2, const double mHl2) const {
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

gslpp::complex THDMcache::B0p_MZ2_mA2_0_mHh2(const double MZ2, const double mA2, const double mHh2) const {
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

gslpp::complex THDMcache::B0p_MZ2_mA2_0_mHp2(const double MZ2, const double mA2, const double mHp2) const {
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

gslpp::complex THDMcache::B0p_MZ2_mA2_mA2_mHl2(const double MZ2, const double mA2, const double mHl2) const {
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

gslpp::complex THDMcache::B0p_MZ2_mA2_mA2_mHh2(const double MZ2, const double mA2, const double mHh2) const {
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

    std::stringstream br1,br2,br3,br4,br5,br6,br7;
    std::stringstream dw1;
    std::stringstream cs1,cs2,cs3,cs4,cs5,cs6,cs7,cs8,cs9;
    std::stringstream cs11,cs12,cs13,cs14,cs15,cs16,cs17,cs18,cs19;
    std::stringstream cs20,cs21;
    std::stringstream csr1,csr2,csr3,csr4;
    std::stringstream csr11,csr12,csr13,csr14;
    std::stringstream ex1,ex2,ex3,ex4,ex5,ex6,ex7,ex8,ex9,ex10,ex11,ex12,ex13,ex14,ex15,ex16,ex17,ex18,ex19,ex20,ex21,ex22,ex23;
    std::stringstream ex1e,ex2e,ex3e,ex4e,ex5e,ex6e,ex7e,ex8e,ex9e,ex10e,ex11e,ex12e,ex13e,ex14e,ex15e,ex16e,ex17e,ex18e,ex19e,ex20e,ex21e,ex22e,ex23e;
//    std::stringstream ex14ep2,ex14em2;
    std::stringstream ex24,ex25,ex26,ex27,ex28,ex29,ex30,ex31,ex32,ex33,ex34,ex35,ex36,ex37,ex38,ex39,ex40,ex41,ex42,ex43,ex44,\
            ex45,ex46,ex47,ex48,ex49,ex50,ex51,ex52,ex53,ex54,ex55,ex56;
    std::stringstream ex24e,ex25e,ex26e,ex27e,ex28e,ex29e,ex30e,ex31e,ex32e,ex33e,ex34e,ex35e,ex36e,ex37e,ex38e,ex39e,ex40e,ex41e,ex42e,ex43e,ex44e,\
            ex45e,ex46e,ex47e,ex48e,ex49e,ex50e,ex51e,ex52e,ex53e,ex54e,ex55e,ex56e;
    std::stringstream ex57,ex58,ex59,ex60,ex61,ex62,ex63,ex64,ex65,ex66,ex67,ex68,ex69,ex70,ex71,ex72,ex73,ex74,ex75,ex76,ex77,\
            ex78,ex79,ex80,ex81,ex82,ex83,ex84,ex85,ex86,ex87,ex88,ex89,ex90,ex91,ex92,ex93,ex94,ex95,ex96,ex97,ex98;
    std::stringstream ex57e,ex58e,ex59e,ex60e,ex61e,ex62e,ex63e,ex64e,ex65e,ex66e,ex67e,ex68e,ex69e,ex70e,ex71e,ex72e,ex73e,ex74e,ex75e,ex76e,ex77e,\
            ex78e,ex79e,ex80e,ex81e,ex82e,ex83e,ex84e,ex85e,ex86e,ex87e,ex88e,ex89e,ex90e,ex91e,ex92e,ex93e,ex94e,ex95e,ex96e,ex97e,ex98e;
    std::stringstream ex99,ex100,ex101,ex102,ex103,ex104,ex105,ex106,ex107,ex108,ex109,ex110,ex111,ex112,ex113,ex114,ex115,ex116,ex117,ex118,ex119,ex120;
    std::stringstream ex99e,ex100e,ex101e,ex102e,ex103e,ex104e,ex105e,ex106e,ex107e,ex108e,ex109e,ex110e,ex111e,ex112e,ex113e,ex114e,ex115e,ex116e,ex117e,ex118e,ex119e,ex120e;
    std::stringstream ex121,ex122,ex123,ex124;
    std::stringstream bsg1;

    std::cout<<"reading tables"<<std::endl;

//    std::cout << "HEPFITTABS = " << getenv("HEPFITPATH") << std::endl;
    std::stringstream path;
    path << getenv("HEPFITTABS") << "/THDM/tabs/";
    std::string tablepath=path.str();

    br1 << tablepath << "br1.dat";
    br_tt = readTable(br1.str(),19961,2);

//// THIS IS FOR THE FUTURE IMPLEMENTATION INTO HEADERS:
//    std::cout<<"br_tt="<<br_tt<<std::endl;
//    double brtt1[4][2];
//    brtt1[0][1]=1;
//        gslpp::matrix<double> brtt1(19861,2,0.);
//    std::stringstream br1x;
//    br1x << "log_cs_ggH_13.h";
//      //brtt1(2)=(3.,4.);
//      brtt1=readTable(br1x.str(),20,2);
//    std::cout<<"brtt1="<<bla1<<std::endl;


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
    ex1 << tablepath << "150400936.dat";
    CMS8_mu_pp_H_VV = readTable(ex1.str(),172,2);
    ex1e << tablepath << "150400936_e.dat";
    CMS8_mu_pp_H_VV_e = readTable(ex1e.str(),172,2);
    ex2 << tablepath << "150404710.dat";
    CMS8_gg_A_hZ_bbll = readTable(ex2.str(),16,2);
    ex2e << tablepath << "150404710_e.dat";
    CMS8_gg_A_hZ_bbll_e = readTable(ex2e.str(),16,2);
    ex3 << tablepath << "160306896.dat";
    CMS8_pp_H_hh_gagabb = readTable(ex3.str(),85,2);
    ex3e << tablepath << "160306896_e.dat";
    CMS8_pp_H_hh_gagabb_e = readTable(ex3e.str(),85,2);
    ex4 << tablepath << "150304114.dat";
    CMS8_pp_H_hh_bbbb = readTable(ex4.str(),167,2);
    ex4e << tablepath << "150304114_e.dat";
    CMS8_pp_H_hh_bbbb_e = readTable(ex4e.str(),167,2);
    ex5 << tablepath << "14076583.dat";
    ATLAS8_pp_phi_gaga = readTable(ex5.str(),108,2);
    ex5e << tablepath << "14076583_e.dat";
    ATLAS8_pp_phi_gaga_e = readTable(ex5e.str(),108,2);
    ex6 << tablepath << "14096064_a.dat";
    ATLAS8_gg_phi_tautau = readTable(ex6.str(),92,2);
    ex6e << tablepath << "14096064_a_e.dat";
    ATLAS8_gg_phi_tautau_e = readTable(ex6e.str(),92,2);
    ex7 << tablepath << "14096064_b.dat";
    ATLAS8_bb_phi_tautau = readTable(ex7.str(),92,2);
    ex7e << tablepath << "14096064_b_e.dat";
    ATLAS8_bb_phi_tautau_e = readTable(ex7e.str(),92,2);
    ex8 << tablepath << "150204478_a.dat";
    ATLAS8_gg_A_hZ_tautauZ = readTable(ex8.str(),79,2);
    ex8e << tablepath << "150204478_a_e.dat";
    ATLAS8_gg_A_hZ_tautauZ_e = readTable(ex8e.str(),79,2);
    ex9 << tablepath << "150204478_b.dat";
    ATLAS8_gg_A_hZ_bbZ = readTable(ex9.str(),79,2);
    ex9e << tablepath << "150204478_b_e.dat";
    ATLAS8_gg_A_hZ_bbZ_e = readTable(ex9e.str(),79,2);
    ex10 << tablepath << "150608329.dat";
    CMS8_bb_phi_bb = readTable(ex10.str(),81,2);
    ex10e << tablepath << "150608329_e.dat";
    CMS8_bb_phi_bb_e = readTable(ex10e.str(),81,2);
    ex11 << tablepath << "150507018.dat";
    ATLAS8_gg_phi_tt = readTable(ex11.str(),53,2);
    ex11e << tablepath << "150507018_e.dat";
    ATLAS8_gg_phi_tt_e = readTable(ex11e.str(),53,2);
    ex12 << tablepath << "CMS-PAS-HIG-14-029_a.dat";
    CMS8_gg_phi_tautau = readTable(ex12.str(),92,2);
    ex12e << tablepath << "CMS-PAS-HIG-14-029_a_e.dat";
    CMS8_gg_phi_tautau_e = readTable(ex12e.str(),92,2);
    ex13 << tablepath << "CMS-PAS-HIG-14-029_b.dat";
    CMS8_bb_phi_tautau = readTable(ex13.str(),92,2);
    ex13e << tablepath << "CMS-PAS-HIG-14-029_b_e.dat";
    CMS8_bb_phi_tautau_e = readTable(ex13e.str(),92,2);
    ex14 << tablepath << "150602301.dat";
    CMS8_gg_phi_gaga = readTable(ex14.str(),141,2);
    ex14e << tablepath << "150602301_e.dat";
    CMS8_gg_phi_gaga_e = readTable(ex14e.str(),141,2);

//    ex14ep1 << tablepath << "150602301_ep1.dat";
//    CMS_ggF_phi_gaga_ep1 = readTable(ex14ep1.str(),141,2);
    //CHANGE THIS DEFINITION!
//    ex14ep2 << tablepath << "150602301_e.dat";
//    CMS_ggF_phi_gaga_ep2 = readTable(ex14ep2.str(),141,2);
//    ex14em1 << tablepath << "150602301_em1.dat";
//    CMS_ggF_phi_gaga_em1 = readTable(ex14em1.str(),141,2);
    //CHANGE THIS DEFINITION!
//    ex14em2 << tablepath << "150602301_e.dat";
//    CMS_ggF_phi_gaga_em2 = readTable(ex14em2.str(),141,2);

    ex15 << tablepath << "150900389_a.dat";
    ATLAS8_gg_H_WW = readTable(ex15.str(),13,2);
    ex15e << tablepath << "150900389_a_e.dat";
    ATLAS8_gg_H_WW_e = readTable(ex15e.str(),13,2);
    ex16 << tablepath << "150900389_b.dat";
    ATLAS8_VBF_H_WW = readTable(ex16.str(),13,2);
    ex16e << tablepath << "150900389_b_e.dat";
    ATLAS8_VBF_H_WW_e = readTable(ex16e.str(),13,2);
    ex17 << tablepath << "150904670.dat";
    ATLAS8_gg_H_hh = readTable(ex17.str(),75,2);
    ex17e << tablepath << "150904670_e.dat";
    ATLAS8_gg_H_hh_e = readTable(ex17e.str(),75,2);
    ex18 << tablepath << "151001181_a.dat";
    CMS8_gg_H_hh_bbtautau = readTable(ex18.str(),10,2);
    ex18e << tablepath << "151001181_a_e.dat";
    CMS8_gg_H_hh_bbtautau_e = readTable(ex18e.str(),10,2);
    ex19 << tablepath << "151001181_b.dat";
    CMS8_gg_A_hZ_tautaull = readTable(ex19.str(),14,2);
    ex19e << tablepath << "151001181_b_e.dat";
    CMS8_gg_A_hZ_tautaull_e = readTable(ex19e.str(),14,2);
    ex20 << tablepath << "150705930_a.dat";
    ATLAS8_gg_H_ZZ = readTable(ex20.str(),173,2);
    ex20e << tablepath << "150705930_a_e.dat";
    ATLAS8_gg_H_ZZ_e = readTable(ex20e.str(),173,2);
    ex21 << tablepath << "150705930_b.dat";
    ATLAS8_VBF_H_ZZ = readTable(ex21.str(),173,2);
    ex21e << tablepath << "150705930_b_e.dat";
    ATLAS8_VBF_H_ZZ_e = readTable(ex21e.str(),173,2);
    ex22 << tablepath << "CMS-PAS-HIG-15-013.dat";
    CMS8_pp_H_hh = readTable(ex22.str(),71,2);
    ex22e << tablepath << "CMS-PAS-HIG-15-013_e.dat";
    CMS8_pp_H_hh_e = readTable(ex22e.str(),71,2);
    ex23 << tablepath << "CMS-PAS-HIG-16-014.dat";
    CMS8_pp_A_Zga_llga = readTable(ex23.str(),101,2);
    ex23e << tablepath << "CMS-PAS-HIG-16-014_e.dat";
    CMS8_pp_A_Zga_llga_e = readTable(ex23e.str(),101,2);
    
    ex24 << tablepath << "ATLAS-CONF-2016-104_b.dat";
    ATLAS13_bb_phi_tt = readTable(ex24.str(),61,2);
    ex24e << tablepath << "ATLAS-CONF-2016-104_b_e.dat";
    ATLAS13_bb_phi_tt_e = readTable(ex24e.str(),61,2);
    ex25 << tablepath << "ATLAS-CONF-2016-104_a.dat";
    ATLAS13_tt_phi_tt = readTable(ex25.str(),61,2);
    ex25e << tablepath << "ATLAS-CONF-2016-104_a_e.dat";
    ATLAS13_tt_phi_tt_e = readTable(ex25e.str(),61,2);
    ex26 << tablepath << "ATLAS-CONF-2017-050_a.dat";
    ATLAS13_gg_phi_tautau = readTable(ex26.str(),206,2);
    ex26e << tablepath << "ATLAS-CONF-2017-050_a_e.dat";
    ATLAS13_gg_phi_tautau_e = readTable(ex26e.str(),206,2);
    ex27 << tablepath << "ATLAS-CONF-2017-050_b.dat";
    ATLAS13_bb_phi_tautau = readTable(ex27.str(),206,2);
    ex27e << tablepath << "ATLAS-CONF-2017-050_b_e.dat";
    ATLAS13_bb_phi_tautau_e = readTable(ex27e.str(),206,2);
    ex28 << tablepath << "170704147.dat";
    ATLAS13_pp_phi_gaga = readTable(ex28.str(),251,2);
    ex28e << tablepath << "170704147_e.dat";
    ATLAS13_pp_phi_gaga_e = readTable(ex28e.str(),251,2);
    ex29 << tablepath << "ATLAS-CONF-2016-044.dat";
    ATLAS13_pp_phi_Zga = readTable(ex29.str(),216,2);
    ex29e << tablepath << "ATLAS-CONF-2016-044_e.dat";
    ATLAS13_pp_phi_Zga_e = readTable(ex29e.str(),216,2);
    ex30 << tablepath << "ATLAS-CONF-2016-056.dat";
    ATLAS13_gg_H_ZZ_llnunu = readTable(ex30.str(),71,2);
    ex30e << tablepath << "ATLAS-CONF-2016-056_e.dat";
    ATLAS13_gg_H_ZZ_llnunu_e = readTable(ex30e.str(),71,2);
    ex31 << tablepath << "ATLAS-CONF-2016-079_a.dat";
    ATLAS13_gg_H_ZZ_llll = readTable(ex31.str(),81,2);
    ex31e << tablepath << "ATLAS-CONF-2016-079_a_e.dat";
    ATLAS13_gg_H_ZZ_llll_e = readTable(ex31e.str(),81,2);
    ex32 << tablepath << "ATLAS-CONF-2016-079_b.dat";
    ATLAS13_VBF_H_ZZ_llll = readTable(ex32.str(),81,2);
    ex32e << tablepath << "ATLAS-CONF-2016-079_b_e.dat";
    ATLAS13_VBF_H_ZZ_llll_e = readTable(ex32e.str(),81,2);
    ex33 << tablepath << "ATLAS-CONF-2016-082_a.dat";
    ATLAS13_gg_H_ZZ_llqq = readTable(ex33.str(),271,2);
    ex33e << tablepath << "ATLAS-CONF-2016-082_a_e.dat";
    ATLAS13_gg_H_ZZ_llqq_e = readTable(ex33e.str(),271,2);
    ex34 << tablepath << "ATLAS-CONF-2016-082_b.dat";
    ATLAS13_VBF_H_ZZ_llqq = readTable(ex34.str(),271,2);
    ex34e << tablepath << "ATLAS-CONF-2016-082_b_e.dat";
    ATLAS13_VBF_H_ZZ_llqq_e = readTable(ex34e.str(),271,2);
    ex35 << tablepath << "ATLAS-CONF-2016-082_c.dat";
    ATLAS13_gg_H_ZZ_nunuqq = readTable(ex35.str(),251,2);
    ex35e << tablepath << "ATLAS-CONF-2016-082_c_e.dat";
    ATLAS13_gg_H_ZZ_nunuqq_e = readTable(ex35e.str(),251,2);
    ex36 << tablepath << "171001123_a.dat";
    ATLAS13_gg_H_WW_enumumu = readTable(ex36.str(),381,2);
    ex36e << tablepath << "171001123_a_e.dat";
    ATLAS13_gg_H_WW_enumumu_e = readTable(ex36e.str(),381,2);
    ex37 << tablepath << "171001123_b.dat";
    ATLAS13_VBF_H_WW_enumumu = readTable(ex37.str(),281,2);
    ex37e << tablepath << "171001123_b_e.dat";
    ATLAS13_VBF_H_WW_enumumu_e = readTable(ex37e.str(),281,2);
    ex38 << tablepath << "171007235_a.dat";
    ATLAS13_gg_H_WW_lnuqq = readTable(ex38.str(),271,2);
    ex38e << tablepath << "171007235_a_e.dat";
    ATLAS13_gg_H_WW_lnuqq_e = readTable(ex38e.str(),271,2);
    ex39 << tablepath << "ATLAS-CONF-2016-049.dat";
    ATLAS13_pp_H_hh_bbbb = readTable(ex39.str(),271,2);
    ex39e << tablepath << "ATLAS-CONF-2016-049_e.dat";
    ATLAS13_pp_H_hh_bbbb_e = readTable(ex39e.str(),271,2);
    ex40 << tablepath << "ATLAS-CONF-2016-004.dat";
    ATLAS13_pp_H_hh_gagabb = readTable(ex40.str(),26,2);
    ex40e << tablepath << "ATLAS-CONF-2016-004_e.dat";
    ATLAS13_pp_H_hh_gagabb_e = readTable(ex40e.str(),26,2);
    ex41 << tablepath << "ATLAS-CONF-2016-071.dat";
    ATLAS13_pp_H_hh_gagaWW = readTable(ex41.str(),25,2);
    ex41e << tablepath << "ATLAS-CONF-2016-071_e.dat";
    ATLAS13_pp_H_hh_gagaWW_e = readTable(ex41e.str(),25,2);
    ex42 << tablepath << "ATLAS-CONF-2017-055_a.dat";
    ATLAS13_gg_A_Zh_Zbb = readTable(ex42.str(),181,2);
    ex42e << tablepath << "ATLAS-CONF-2017-055_a_e.dat";
    ATLAS13_gg_A_Zh_Zbb_e = readTable(ex42e.str(),181,2);
    ex43 << tablepath << "ATLAS-CONF-2017-055_b.dat";
    ATLAS13_bb_A_Zh_Zbb = readTable(ex43.str(),181,2);
    ex43e << tablepath << "ATLAS-CONF-2017-055_b_e.dat";
    ATLAS13_bb_A_Zh_Zbb_e = readTable(ex43e.str(),181,2);
    ex44 << tablepath << "CMS-PAS-HIG-16-025.dat";
    CMS13_pp_phi_bb = readTable(ex44.str(),66,2);
    ex44e << tablepath << "CMS-PAS-HIG-16-025_e.dat";
    CMS13_pp_phi_bb_e = readTable(ex44e.str(),66,2);
    ex45 << tablepath << "CMS-PAS-HIG-16-037_a.dat";
    CMS13_gg_phi_tautau = readTable(ex45.str(),312,2);
    ex45e << tablepath << "CMS-PAS-HIG-16-037_a_e.dat";
    CMS13_gg_phi_tautau_e = readTable(ex45e.str(),312,2);
    ex46 << tablepath << "CMS-PAS-HIG-16-037_b.dat";
    CMS13_bb_phi_tautau = readTable(ex46.str(),312,2);
    ex46e << tablepath << "CMS-PAS-HIG-16-037_b_e.dat";
    CMS13_bb_phi_tautau_e = readTable(ex46e.str(),312,2);
    ex47 << tablepath << "CMS-PAS-EXO-16-027.dat";
    CMS13_gg_phi_gaga = readTable(ex47.str(),351,2);
    ex47e << tablepath << "CMS-PAS-EXO-16-027_e.dat";
    CMS13_gg_phi_gaga_e = readTable(ex47e.str(),351,2);
    ex48 << tablepath << "CMS-PAS-EXO-16-034.dat";
    CMS13_pp_phi_Zga_llga = readTable(ex48.str(),171,2);
    ex48e << tablepath << "CMS-PAS-EXO-16-034_e.dat";
    CMS13_pp_phi_Zga_llga_e = readTable(ex48e.str(),171,2);
    ex49 << tablepath << "CMS-PAS-EXO-16-035.dat";
    CMS13_pp_phi_Zga_qqga = readTable(ex49.str(),236,2);
    ex49e << tablepath << "CMS-PAS-EXO-16-035_e.dat";
    CMS13_pp_phi_Zga_qqga_e = readTable(ex49e.str(),236,2);
    ex50 << tablepath << "CMS-PAS-HIG-16-033_a.dat";
    CMS13_pp_H_ZZ_llll = readTable(ex50.str(),241,2);
    ex50e << tablepath << "CMS-PAS-HIG-16-033_a_e.dat";
    CMS13_pp_H_ZZ_llll_e = readTable(ex50e.str(),241,2);
    ex51 << tablepath << "CMS-PAS-HIG-16-033_b.dat";
    CMS13_VBFVH_H_ZZ_llll = readTable(ex51.str(),241,2);
    ex51e << tablepath << "CMS-PAS-HIG-16-033_b_e.dat";
    CMS13_VBFVH_H_ZZ_llll_e = readTable(ex51e.str(),241,2);
    ex52 << tablepath << "CMS-PAS-HIG-16-023.dat";
    CMS13_ggFVBF_H_WW_lnulnu = readTable(ex52.str(),81,2);
    ex52e << tablepath << "CMS-PAS-HIG-16-023_e.dat";
    CMS13_ggFVBF_H_WW_lnulnu_e = readTable(ex52e.str(),81,2);
    ex53 << tablepath << "CMS-PAS-HIG-17-009.dat";
    CMS13_pp_H_hh_bbbb = readTable(ex53.str(),95,2);
    ex53e << tablepath << "CMS-PAS-HIG-17-009_e.dat";
    CMS13_pp_H_hh_bbbb_e = readTable(ex53e.str(),95,2);
    ex54 << tablepath << "CMS-PAS-HIG-17-008.dat";
    CMS13_pp_H_hh_gagabb = readTable(ex54.str(),66,2);
    ex54e << tablepath << "CMS-PAS-HIG-17-008_e.dat";
    CMS13_pp_H_hh_gagabb_e = readTable(ex54e.str(),66,2);
    ex55 << tablepath << "CMS-PAS-HIG-16-029.dat";
    CMS13_pp_H_hh_bbtautau = readTable(ex55.str(),66,2);
    ex55e << tablepath << "CMS-PAS-HIG-16-029_e.dat";
    CMS13_pp_H_hh_bbtautau_e = readTable(ex55e.str(),66,2);
    ex56 << tablepath << "CMS-PAS-HIG-16-011.dat";
    CMS13_pp_H_hh_bblnulnu = readTable(ex56.str(),65,2);
    ex56e << tablepath << "CMS-PAS-HIG-16-011_e.dat";
    CMS13_pp_H_hh_bblnulnu_e = readTable(ex56e.str(),65,2);

    ex57 << tablepath << "t1.dat";
    temp1 = readTable(ex57.str(),1,2);
    ex57e << tablepath << "t1_e.dat";
    temp1 = readTable(ex57e.str(),1,2);
    ex58 << tablepath << "t2.dat";
    temp2 = readTable(ex58.str(),1,2);
    ex58e << tablepath << "t2_e.dat";
    temp2 = readTable(ex58e.str(),1,2);
    ex59 << tablepath << "t3.dat";
    temp3 = readTable(ex59.str(),1,2);
    ex59e << tablepath << "t3_e.dat";
    temp3 = readTable(ex59e.str(),1,2);
    ex60 << tablepath << "t4.dat";
    temp4 = readTable(ex60.str(),1,2);
    ex60e << tablepath << "t4_e.dat";
    temp4 = readTable(ex60e.str(),1,2);
    ex61 << tablepath << "t5.dat";
    temp5 = readTable(ex61.str(),1,2);
    ex61e << tablepath << "t5_e.dat";
    temp5 = readTable(ex61e.str(),1,2);
    ex62 << tablepath << "t6.dat";
    temp6 = readTable(ex62.str(),1,2);
    ex62e << tablepath << "t6_e.dat";
    temp6 = readTable(ex62e.str(),1,2);
    ex63 << tablepath << "t7.dat";
    temp7 = readTable(ex63.str(),1,2);
    ex63e << tablepath << "t7_e.dat";
    temp7 = readTable(ex63e.str(),1,2);
    ex64 << tablepath << "t8.dat";
    temp8 = readTable(ex64.str(),1,2);
    ex64e << tablepath << "t8_e.dat";
    temp8 = readTable(ex64e.str(),1,2);
    ex65 << tablepath << "t9.dat";
    temp9 = readTable(ex65.str(),1,2);
    ex65e << tablepath << "t9_e.dat";
    temp9 = readTable(ex65e.str(),1,2);
    ex66 << tablepath << "t10.dat";
    temp10 = readTable(ex66.str(),1,2);
    ex66e << tablepath << "t10_e.dat";
    temp10 = readTable(ex66e.str(),1,2);
    ex67 << tablepath << "t11.dat";
    temp11 = readTable(ex67.str(),1,2);
    ex67e << tablepath << "t11_e.dat";
    temp11 = readTable(ex67e.str(),1,2);
    ex68 << tablepath << "t12.dat";
    temp12 = readTable(ex68.str(),1,2);
    ex68e << tablepath << "t12_e.dat";
    temp12 = readTable(ex68e.str(),1,2);
    ex69 << tablepath << "t13.dat";
    temp13 = readTable(ex69.str(),1,2);
    ex69e << tablepath << "t13_e.dat";
    temp13 = readTable(ex69e.str(),1,2);
    ex70 << tablepath << "t14.dat";
    temp14 = readTable(ex70.str(),1,2);
    ex70e << tablepath << "t14_e.dat";
    temp14 = readTable(ex70e.str(),1,2);
    ex71 << tablepath << "t15.dat";
    temp15 = readTable(ex71.str(),1,2);
    ex71e << tablepath << "t15_e.dat";
    temp15 = readTable(ex71e.str(),1,2);
    ex72 << tablepath << "t16.dat";
    temp16 = readTable(ex72.str(),1,2);
    ex72e << tablepath << "t16_e.dat";
    temp16 = readTable(ex72e.str(),1,2);
    ex73 << tablepath << "t17.dat";
    temp17 = readTable(ex73.str(),1,2);
    ex73e << tablepath << "t17_e.dat";
    temp17 = readTable(ex73e.str(),1,2);
    ex74 << tablepath << "t18.dat";
    temp18 = readTable(ex74.str(),1,2);
    ex74e << tablepath << "t18_e.dat";
    temp18 = readTable(ex74e.str(),1,2);
    ex75 << tablepath << "t19.dat";
    temp19 = readTable(ex75.str(),1,2);
    ex75e << tablepath << "t19_e.dat";
    temp19 = readTable(ex75e.str(),1,2);
    ex76 << tablepath << "t20.dat";
    temp20 = readTable(ex76.str(),1,2);
    ex76e << tablepath << "t20_e.dat";
    temp20 = readTable(ex76e.str(),1,2);
    ex77 << tablepath << "t21.dat";
    temp21 = readTable(ex77.str(),1,2);
    ex77e << tablepath << "t21_e.dat";
    temp21 = readTable(ex77e.str(),1,2);
    ex78 << tablepath << "t22.dat";
    temp22 = readTable(ex78.str(),1,2);
    ex78e << tablepath << "t22_e.dat";
    temp22 = readTable(ex78e.str(),1,2);
    ex79 << tablepath << "t23.dat";
    temp23 = readTable(ex79.str(),1,2);
    ex79e << tablepath << "t23_e.dat";
    temp23 = readTable(ex79e.str(),1,2);
    ex80 << tablepath << "t24.dat";
    temp24 = readTable(ex80.str(),1,2);
    ex80e << tablepath << "t24_e.dat";
    temp24 = readTable(ex80e.str(),1,2);
    ex81 << tablepath << "t25.dat";
    temp25 = readTable(ex81.str(),1,2);
    ex81e << tablepath << "t25_e.dat";
    temp25 = readTable(ex81e.str(),1,2);
    ex82 << tablepath << "t26.dat";
    temp26 = readTable(ex82.str(),1,2);
    ex82e << tablepath << "t26_e.dat";
    temp26 = readTable(ex82e.str(),1,2);
    ex83 << tablepath << "t27.dat";
    temp27 = readTable(ex83.str(),1,2);
    ex83e << tablepath << "t27_e.dat";
    temp27 = readTable(ex83e.str(),1,2);
    ex84 << tablepath << "t28.dat";
    temp28 = readTable(ex84.str(),1,2);
    ex84e << tablepath << "t28_e.dat";
    temp28 = readTable(ex84e.str(),1,2);
    ex85 << tablepath << "t29.dat";
    temp29 = readTable(ex85.str(),1,2);
    ex85e << tablepath << "t29_e.dat";
    temp29 = readTable(ex85e.str(),1,2);
    ex86 << tablepath << "t30.dat";
    temp30 = readTable(ex86.str(),1,2);
    ex86e << tablepath << "t30_e.dat";
    temp30e = readTable(ex86e.str(),1,2);
    ex87 << tablepath << "t31.dat";
    temp31 = readTable(ex87.str(),1,2);
    ex87e << tablepath << "t31_e.dat";
    temp31e = readTable(ex87e.str(),1,2);
    ex88 << tablepath << "t32.dat";
    temp32 = readTable(ex88.str(),1,2);
    ex88e << tablepath << "t32_e.dat";
    temp32e = readTable(ex88e.str(),1,2);
    ex89 << tablepath << "t33.dat";
    temp33 = readTable(ex89.str(),1,2);
    ex89e << tablepath << "t33_e.dat";
    temp33e = readTable(ex89e.str(),1,2);
    ex90 << tablepath << "t34.dat";
    temp34 = readTable(ex90.str(),1,2);
    ex90e << tablepath << "t34_e.dat";
    temp34e = readTable(ex90e.str(),1,2);
    ex91 << tablepath << "t35.dat";
    temp35 = readTable(ex91.str(),1,2);
    ex91e << tablepath << "t35_e.dat";
    temp35e = readTable(ex91e.str(),1,2);
    ex92 << tablepath << "t36.dat";
    temp36 = readTable(ex92.str(),1,2);
    ex92e << tablepath << "t36_e.dat";
    temp36e = readTable(ex92e.str(),1,2);
    ex93 << tablepath << "t37.dat";
    temp37 = readTable(ex93.str(),1,2);
    ex93e << tablepath << "t37_e.dat";
    temp37e = readTable(ex93e.str(),1,2);
    ex94 << tablepath << "t38.dat";
    temp38 = readTable(ex94.str(),1,2);
    ex94e << tablepath << "t38_e.dat";
    temp38e = readTable(ex94e.str(),1,2);
    ex95 << tablepath << "t39.dat";
    temp39 = readTable(ex95.str(),1,2);
    ex95e << tablepath << "t39_e.dat";
    temp39e = readTable(ex95e.str(),1,2);
    ex96 << tablepath << "t40.dat";
    temp40 = readTable(ex96.str(),1,2);
    ex96e << tablepath << "t40_e.dat";
    temp40e = readTable(ex96e.str(),1,2);

    ex97 << tablepath << "CMS-PAS-HIG-16-034.dat";
    CMS13_pp_H_ZZ_llqq = readTable(ex97.str(),151,2);
    ex97e << tablepath << "CMS-PAS-HIG-16-034_e.dat";
    CMS13_pp_H_ZZ_llqq_e = readTable(ex97e.str(),151,2);
    ex98 << tablepath << "14078150.dat";
    ATLAS8_pp_phi_Zga_llga = readTable(ex98.str(),141,2);
    ex98e << tablepath << "14078150.dat";
    ATLAS8_pp_phi_Zga_llga_e = readTable(ex98e.str(),141,2);

    ex99 << tablepath << "14126663.dat";
    ATLAS8_pp_Hpm_taunu = readTable(ex99.str(),83,2);
    ex99e << tablepath << "14126663_e.dat";
    ATLAS8_pp_Hpm_taunu_e = readTable(ex99e.str(),83,2);
    ex100 << tablepath << "151203704.dat";
    ATLAS8_pp_Hpm_tb = readTable(ex100.str(),41,2);
    ex100e << tablepath << "151203704_e.dat";
    ATLAS8_pp_Hpm_tb_e = readTable(ex100e.str(),41,2);
    ex101 << tablepath << "150807774_a.dat";
    CMS8_pp_Hp_taunu = readTable(ex101.str(),43,2);
    ex101e << tablepath << "150807774_a_e.dat";
    CMS8_pp_Hp_taunu_e = readTable(ex101e.str(),43,2);
    ex102 << tablepath << "150807774_b.dat";
    CMS8_pp_Hp_tb = readTable(ex102.str(),43,2);
    ex102e << tablepath << "150807774_b_e.dat";
    CMS8_pp_Hp_tb_e = readTable(ex102e.str(),43,2);
    ex103 << tablepath << "ATLAS-CONF-2016-088.dat";
    ATLAS13_pp_Hpm_taunu = readTable(ex103.str(),181,2);
    ex103e << tablepath << "ATLAS-CONF-2016-088_e.dat";
    ATLAS13_pp_Hpm_taunu_e = readTable(ex103e.str(),181,2);
    ex104 << tablepath << "ATLAS-CONF-2016-089.dat";
    ATLAS13_pp_Hp_tb1 = readTable(ex104.str(),71,2);
    ex104e << tablepath << "ATLAS-CONF-2016-089_e.dat";
    ATLAS13_pp_Hp_tb1_e = readTable(ex104e.str(),71,2);
    ex105 << tablepath << "ATLAS-CONF-2016-104_c.dat";
    ATLAS13_pp_Hp_tb2 = readTable(ex105.str(),181,2);
    ex105e << tablepath << "ATLAS-CONF-2016-104_c_e.dat";
    ATLAS13_pp_Hp_tb2_e = readTable(ex105e.str(),181,2);
    ex106 << tablepath << "CMS-PAS-HIG-16-031.dat";
    CMS13_pp_Hpm_taunu = readTable(ex106.str(),283,2);
    ex106e << tablepath << "CMS-PAS-HIG-16-031_e.dat";
    CMS13_pp_Hpm_taunu_e = readTable(ex106e.str(),283,2);

    ex107 << tablepath << "CMS-PAS-HIG-17-002.dat";
    CMS13_pp_H_hh_bbtautau1 = readTable(ex107.str(),66,2);
    ex107e << tablepath << "CMS-PAS-HIG-17-002_e.dat";
    CMS13_pp_H_hh_bbtautau1_e = readTable(ex107e.str(),66,2);
    ex108 << tablepath << "170804188.dat";
    CMS13_pp_H_hh_bbVV = readTable(ex108.str(),65,2);
    ex108e << tablepath << "170804188_e.dat";
    CMS13_pp_H_hh_bbVV_e = readTable(ex108e.str(),65,2);
    ex109 << tablepath << "CMS-PAS-EXO-17-005.dat";
    CMS13_ggF_phi_Zga = readTable(ex109.str(),366,2);
    ex109e << tablepath << "CMS-PAS-EXO-17-005_e.dat";
    CMS13_ggF_phi_Zga_e = readTable(ex109e.str(),366,2);
    ex110 << tablepath << "171004960.dat";
    CMS13_ggF_H_hh_bbbb = readTable(ex110.str(),226,2);
    ex110e << tablepath << "171004960_e.dat";
    CMS13_ggF_H_hh_bbbb_e = readTable(ex110e.str(),226,2);

    ex111 << tablepath << "171007235_b.dat";
    ATLAS13_VBF_H_WW_lnuqq = readTable(ex111.str(),271,2);
    ex111e << tablepath << "171007235_b_e.dat";
    ATLAS13_VBF_H_WW_lnuqq_e = readTable(ex111e.str(),271,2);
    ex112 << tablepath << "170800212.dat";
    ATLAS13_gg_phi_Zga_llga = readTable(ex112.str(),216,2);
    ex112e << tablepath << "170800212_e.dat";
    ATLAS13_gg_phi_Zga_llga_e = readTable(ex112e.str(),216,2);
    ex113 << tablepath << "ATLAS-CONF-2017-058_a.dat";
    ATLAS13_gg_H_ZZ_llllnunu = readTable(ex113.str(),101,2);
    ex113e << tablepath << "ATLAS-CONF-2017-058_a_e.dat";
    ATLAS13_gg_H_ZZ_llllnunu_e = readTable(ex113e.str(),101,2);
    ex114 << tablepath << "ATLAS-CONF-2017-058_b.dat";
    ATLAS13_VBF_H_ZZ_llllnunu = readTable(ex114.str(),101,2);
    ex114e << tablepath << "ATLAS-CONF-2017-058_b_e.dat";
    ATLAS13_VBF_H_ZZ_llllnunu_e = readTable(ex114e.str(),101,2);
    ex115 << tablepath << "170809638_a.dat";
    ATLAS13_gg_H_ZZ_qqllnunu = readTable(ex115.str(),271,2);
    ex115e << tablepath << "170809638_a_e.dat";
    ATLAS13_gg_H_ZZ_qqllnunu_e = readTable(ex115e.str(),271,2);
    ex116 << tablepath << "170809638_b.dat";
    ATLAS13_VBF_H_ZZ_qqllnunu = readTable(ex116.str(),271,2);
    ex116e << tablepath << "170809638_b_e.dat";
    ATLAS13_VBF_H_ZZ_qqllnunu_e = readTable(ex116e.str(),271,2);
    ex117 << tablepath << "170804445.dat";
    ATLAS13_pp_H_VV_qqqq = readTable(ex117.str(),181,2);
    ex117e << tablepath << "170804445_e.dat";
    ATLAS13_pp_H_VV_qqqq_e = readTable(ex117e.str(),181,2);
    ex118 << tablepath << "CMS-PAS-B2G-16-023.dat";
    CMS13_pp_H_ZZ_llnunu = readTable(ex118.str(),191,2);
    ex118e << tablepath << "CMS-PAS-B2G-16-023_e.dat";
    CMS13_pp_H_ZZ_llnunu_e = readTable(ex118e.str(),191,2);
    ex119 << tablepath << "CMS-PAS-HIG-16-001_a.dat";
    CMS13_gg_H_ZZ_llnunu = readTable(ex119.str(),131,2);
    ex119e << tablepath << "CMS-PAS-HIG-16-001_a_e.dat";
    CMS13_gg_H_ZZ_llnunu_e = readTable(ex119e.str(),131,2);
    ex120 << tablepath << "CMS-PAS-HIG-16-001_b.dat";
    CMS13_VBF_H_ZZ_llnunu = readTable(ex120.str(),131,2);
    ex120e << tablepath << "CMS-PAS-HIG-16-001_b_e.dat";
    CMS13_VBF_H_ZZ_llnunu_e = readTable(ex120e.str(),131,2);

    ex121 << tablepath << "160302991_a.dat";
    CMS8_pp_A_HZ_bbll = readTable(ex121.str(),28718,3);
    ex122 << tablepath << "160302991_b.dat";
    CMS8_pp_H_AZ_bbll = readTable(ex122.str(),29050,3);
    ex123 << tablepath << "160302991_c.dat";
    CMS8_pp_A_HZ_tautaull = readTable(ex123.str(),400,3);
    ex124 << tablepath << "160302991_d.dat";
    CMS8_pp_H_AZ_tautaull = readTable(ex124.str(),400,3);

    bsg1 << tablepath << "bsgammatable.dat";
    arraybsgamma = readTable(bsg1.str(),1111,3);
}    
    


double THDMcache::ip_Br_HPtott(double mass){
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



double THDMcache::ip_Br_HPtobb(double mass){
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



double THDMcache::ip_Br_HPtotautau(double mass){
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



double THDMcache::ip_Br_HPtocc(double mass){
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



double THDMcache::ip_Br_HPtomumu(double mass){
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



double THDMcache::ip_Br_HPtoZZ(double mass){
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



double THDMcache::ip_Br_HPtoWW(double mass){
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



double THDMcache::ip_GammaHPtotSM(double mass){
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



double THDMcache::ip_cs_ggtoH_8(double mass){
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



double THDMcache::ip_cs_ggtoH_13(double mass){
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



double THDMcache::ip_cs_VBFtoH_8(double mass){
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



double THDMcache::ip_cs_VBFtoH_13(double mass){
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



double THDMcache::ip_cs_WtoWH_8(double mass){
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



double THDMcache::ip_cs_WtoWH_13(double mass){
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



double THDMcache::ip_cs_ZtoZH_8(double mass){
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



double THDMcache::ip_cs_ZtoZH_13(double mass){
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



double THDMcache::ip_cs_pptottH_8(double mass){
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



double THDMcache::ip_cs_pptottH_13(double mass){
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



double THDMcache::ip_cs_pptobbH_8(double mass){
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



double THDMcache::ip_cs_pptobbH_13(double mass){
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



double THDMcache::ip_cs_ggtoA_8(double mass){
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



double THDMcache::ip_cs_ggtoA_13(double mass){
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



double THDMcache::ip_cs_pptottA_8(double mass){
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



double THDMcache::ip_cs_pptottA_13(double mass){
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



double THDMcache::ip_cs_pptobbA_8(double mass){
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



double THDMcache::ip_cs_pptobbA_13(double mass){
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



double THDMcache::ip_cs_ggtoHp_8(double mHp, double logtb){
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



double THDMcache::ip_cs_ggtoHp_13(double mHp, double logtb){
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



double THDMcache::ip_csr_ggH_t_8(double mass){
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



double THDMcache::ip_csr_ggH_t_13(double mass){
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



double THDMcache::ip_csr_ggH_b_8(double mass){
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



double THDMcache::ip_csr_ggH_b_13(double mass){
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



double THDMcache::ip_csr_ggA_t_8(double mass){
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



double THDMcache::ip_csr_ggA_t_13(double mass){
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



double THDMcache::ip_csr_ggA_b_8(double mass){
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



double THDMcache::ip_csr_ggA_b_13(double mass){
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



double THDMcache::ip_ex_pp_phi_gaga_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_gaga_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_phi_gaga_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS8_pp_phi_gaga,mass);
        CacheShiftReal(ip_ex_pp_phi_gaga_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_phi_gaga_ATLAS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_gaga_ATLAS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_phi_gaga_ATLAS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS8_pp_phi_gaga_e,mass);
        CacheShiftReal(ip_ex_pp_phi_gaga_ATLAS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_phi_Zga_llga_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_Zga_llga_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_phi_Zga_llga_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS8_pp_phi_Zga_llga,mass);
        CacheShiftReal(ip_ex_pp_phi_Zga_llga_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_phi_Zga_llga_ATLAS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_Zga_llga_ATLAS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_phi_Zga_llga_ATLAS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS8_pp_phi_Zga_llga_e,mass);
        CacheShiftReal(ip_ex_pp_phi_Zga_llga_ATLAS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_phi_tautau_ATLAS8(double mass){
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



double THDMcache::ip_ex_gg_phi_tautau_ATLAS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_tautau_ATLAS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_tautau_ATLAS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS8_gg_phi_tautau_e,mass);
        CacheShiftReal(ip_ex_gg_phi_tautau_ATLAS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_bb_phi_tautau_ATLAS8(double mass){
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



double THDMcache::ip_ex_bb_phi_tautau_ATLAS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_tautau_ATLAS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phi_tautau_ATLAS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS8_bb_phi_tautau_e,mass);
        CacheShiftReal(ip_ex_bb_phi_tautau_ATLAS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_A_hZ_tautauZ_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_A_hZ_tautauZ_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_A_hZ_tautauZ_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS8_gg_A_hZ_tautauZ,mass);
        CacheShiftReal(ip_ex_gg_A_hZ_tautauZ_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_A_hZ_tautauZ_ATLAS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_A_hZ_tautauZ_ATLAS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_A_hZ_tautauZ_ATLAS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS8_gg_A_hZ_tautauZ_e,mass);
        CacheShiftReal(ip_ex_gg_A_hZ_tautauZ_ATLAS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_A_hZ_bbZ_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_A_hZ_bbZ_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_A_hZ_bbZ_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS8_gg_A_hZ_bbZ,mass);
        CacheShiftReal(ip_ex_gg_A_hZ_bbZ_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}





double THDMcache::ip_ex_gg_A_hZ_bbZ_ATLAS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_A_hZ_bbZ_ATLAS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_A_hZ_bbZ_ATLAS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS8_gg_A_hZ_bbZ_e,mass);
        CacheShiftReal(ip_ex_gg_A_hZ_bbZ_ATLAS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_phi_tt_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_tt_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_tt_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS8_gg_phi_tt,mass);
        CacheShiftReal(ip_ex_gg_phi_tt_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_phi_tt_ATLAS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_tt_ATLAS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_tt_ATLAS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS8_gg_phi_tt_e,mass);
        CacheShiftReal(ip_ex_gg_phi_tt_ATLAS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_H_WW_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_H_WW_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_H_WW_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS8_gg_H_WW,mass);
        CacheShiftReal(ip_ex_gg_H_WW_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_H_WW_ATLAS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_H_WW_ATLAS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_H_WW_ATLAS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS8_gg_H_WW_e,mass);
        CacheShiftReal(ip_ex_gg_H_WW_ATLAS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_VBF_H_WW_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VBF_H_WW_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_VBF_H_WW_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS8_VBF_H_WW,mass);
        CacheShiftReal(ip_ex_VBF_H_WW_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_VBF_H_WW_ATLAS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VBF_H_WW_ATLAS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_VBF_H_WW_ATLAS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS8_VBF_H_WW_e,mass);
        CacheShiftReal(ip_ex_VBF_H_WW_ATLAS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_H_ZZ_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_H_ZZ_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_H_ZZ_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS8_gg_H_ZZ,mass);
        CacheShiftReal(ip_ex_gg_H_ZZ_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_H_ZZ_ATLAS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_H_ZZ_ATLAS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_H_ZZ_ATLAS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS8_gg_H_ZZ_e,mass);
        CacheShiftReal(ip_ex_gg_H_ZZ_ATLAS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_VBF_H_ZZ_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VBF_H_ZZ_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_VBF_H_ZZ_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS8_VBF_H_ZZ,mass);
        CacheShiftReal(ip_ex_VBF_H_ZZ_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_VBF_H_ZZ_ATLAS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VBF_H_ZZ_ATLAS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_VBF_H_ZZ_ATLAS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS8_VBF_H_ZZ_e,mass);
        CacheShiftReal(ip_ex_VBF_H_ZZ_ATLAS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_H_hh_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_H_hh_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_H_hh_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS8_gg_H_hh,mass);
        CacheShiftReal(ip_ex_gg_H_hh_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_H_hh_ATLAS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_H_hh_ATLAS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_H_hh_ATLAS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS8_gg_H_hh_e,mass);
        CacheShiftReal(ip_ex_gg_H_hh_ATLAS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_mu_pp_H_VV_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_mu_pp_H_VV_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_mu_pp_H_VV_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS8_mu_pp_H_VV,mass);
        CacheShiftReal(ip_ex_mu_pp_H_VV_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_mu_pp_H_VV_CMS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_mu_pp_H_VV_CMS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_mu_pp_H_VV_CMS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate(CMS8_mu_pp_H_VV_e,mass);
        CacheShiftReal(ip_ex_mu_pp_H_VV_CMS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_A_hZ_bbll_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_A_hZ_bbll_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_A_hZ_bbll_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS8_gg_A_hZ_bbll,mass);
        CacheShiftReal(ip_ex_gg_A_hZ_bbll_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_A_hZ_bbll_CMS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_A_hZ_bbll_CMS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_A_hZ_bbll_CMS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate(CMS8_gg_A_hZ_bbll_e,mass);
        CacheShiftReal(ip_ex_gg_A_hZ_bbll_CMS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_hh_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_hh_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_H_hh_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS8_pp_H_hh,mass);
        CacheShiftReal(ip_ex_pp_H_hh_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_hh_CMS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_hh_CMS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_H_hh_CMS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate(CMS8_pp_H_hh_e,mass);
        CacheShiftReal(ip_ex_pp_H_hh_CMS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_phi_hh_gagabb_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_hh_gagabb_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_phi_hh_gagabb_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS8_pp_H_hh_gagabb,mass);
        CacheShiftReal(ip_ex_pp_phi_hh_gagabb_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_phi_hh_gagabb_CMS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_hh_gagabb_CMS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_phi_hh_gagabb_CMS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate(CMS8_pp_H_hh_gagabb_e,mass);
        CacheShiftReal(ip_ex_pp_phi_hh_gagabb_CMS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_phi_hh_bbbb_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_hh_bbbb_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_phi_hh_bbbb_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS8_pp_H_hh_bbbb,mass);
        CacheShiftReal(ip_ex_pp_phi_hh_bbbb_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_phi_hh_bbbb_CMS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_hh_bbbb_CMS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_phi_hh_bbbb_CMS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate(CMS8_pp_H_hh_bbbb_e,mass);
        CacheShiftReal(ip_ex_pp_phi_hh_bbbb_CMS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_bb_phi_bb_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_bb_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phi_bb_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_bb_phi_bb,mass);
        CacheShiftReal(ip_ex_bb_phi_bb_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_bb_phi_bb_CMS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_bb_CMS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phi_bb_CMS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_bb_phi_bb_e,mass);
        CacheShiftReal(ip_ex_bb_phi_bb_CMS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_phi_tautau_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_tautau_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_tautau_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_gg_phi_tautau,mass);
        CacheShiftReal(ip_ex_gg_phi_tautau_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_phi_tautau_CMS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_tautau_CMS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_tautau_CMS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_gg_phi_tautau_e,mass);
        CacheShiftReal(ip_ex_gg_phi_tautau_CMS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_bb_phi_tautau_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_tautau_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phi_tautau_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_bb_phi_tautau,mass);
        CacheShiftReal(ip_ex_bb_phi_tautau_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_bb_phi_tautau_CMS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_tautau_CMS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phi_tautau_CMS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_bb_phi_tautau_e,mass);
        CacheShiftReal(ip_ex_bb_phi_tautau_CMS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_phi_gaga_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_gaga_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_gaga_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_gg_phi_gaga,mass);
        CacheShiftReal(ip_ex_gg_phi_gaga_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_phi_gaga_CMS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_gaga_CMS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_gaga_CMS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_gg_phi_gaga_e,mass);
        CacheShiftReal(ip_ex_gg_phi_gaga_CMS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



//
//double THDMcache::ip_ex_gg_phi_gaga_CMS_ep1(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_phi_gaga_CMS_cache_ep1, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_phi_gaga_CMS_cache_ep1[NumPar][i] );
//    } else {
//        double newResult = interpolate (CMS_ggF_phi_gaga_ep1,mass);
//        CacheShiftReal(ip_ex_gg_phi_gaga_CMS_cache_ep1, NumPar, params, newResult);
//        return newResult;
//    }
//}

//double THDMcache::ip_ex_gg_phi_gaga_CMS_ep2(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_phi_gaga_CMS_cache_ep2, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_phi_gaga_CMS_cache_ep2[NumPar][i] );
//    } else {
//        double newResult = interpolate (CMS_ggF_phi_gaga_ep2,mass);
//        CacheShiftReal(ip_ex_gg_phi_gaga_CMS_cache_ep2, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//double THDMcache::ip_ex_gg_phi_gaga_CMS_em1(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_phi_gaga_CMS_cache_em1, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_phi_gaga_CMS_cache_em1[NumPar][i] );
//    } else {
//        double newResult = interpolate (CMS_ggF_phi_gaga_em1,mass);
//        CacheShiftReal(ip_ex_gg_phi_gaga_CMS_cache_em1, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//double THDMcache::ip_ex_gg_phi_gaga_CMS_em2(double mass){
//    int NumPar = 1;
//    double params[] = {mass};
//
//    int i = CacheCheckReal(ip_ex_gg_phi_gaga_CMS_cache_em2, NumPar, params);
//    if (i>=0) {
//        return ( ip_ex_gg_phi_gaga_CMS_cache_em2[NumPar][i] );
//    } else {
//        double newResult = interpolate (CMS_ggF_phi_gaga_em2,mass);
//        CacheShiftReal(ip_ex_gg_phi_gaga_CMS_cache_em2, NumPar, params, newResult);
//        return newResult;
//    }
//}
//


double THDMcache::ip_ex_pp_A_Zga_llga_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_A_Zga_llga_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_A_Zga_llga_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_pp_A_Zga_llga,mass);
        CacheShiftReal(ip_ex_pp_A_Zga_llga_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_A_Zga_llga_CMS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_A_Zga_llga_CMS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_A_Zga_llga_CMS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_pp_A_Zga_llga_e,mass);
        CacheShiftReal(ip_ex_pp_A_Zga_llga_CMS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_H_hh_bbtautau_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_H_hh_bbtautau_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_H_hh_bbtautau_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_gg_H_hh_bbtautau,mass);
        CacheShiftReal(ip_ex_gg_H_hh_bbtautau_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_H_hh_bbtautau_CMS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_H_hh_bbtautau_CMS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_H_hh_bbtautau_CMS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_gg_H_hh_bbtautau_e,mass);
        CacheShiftReal(ip_ex_gg_H_hh_bbtautau_CMS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_A_hZ_tautaull_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_A_hZ_tautaull_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_A_hZ_tautaull_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_gg_A_hZ_tautaull,mass);
        CacheShiftReal(ip_ex_gg_A_hZ_tautaull_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_A_hZ_tautaull_CMS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_A_hZ_tautaull_CMS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_A_hZ_tautaull_CMS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_gg_A_hZ_tautaull_e,mass);
        CacheShiftReal(ip_ex_gg_A_hZ_tautaull_CMS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_A_HZ_bbll_CMS8(double mA, double mH){
    int NumPar = 2;
    double params[] = {mA, mH};

    int i = CacheCheckReal(ip_ex_pp_A_HZ_bbll_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_A_HZ_bbll_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate2D(CMS8_pp_A_HZ_bbll, mA, mH);
        CacheShiftReal(ip_ex_pp_A_HZ_bbll_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_AZ_bbll_CMS8(double mA, double mH){
    int NumPar = 2;
    double params[] = {mA, mH};

    int i = CacheCheckReal(ip_ex_pp_H_AZ_bbll_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_H_AZ_bbll_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate2D(CMS8_pp_H_AZ_bbll, mA, mH);
        CacheShiftReal(ip_ex_pp_H_AZ_bbll_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_A_HZ_tautaull_CMS8(double mA, double mH){
    int NumPar = 2;
    double params[] = {mA, mH};

    int i = CacheCheckReal(ip_ex_pp_A_HZ_tautaull_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_A_HZ_tautaull_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate2D(CMS8_pp_A_HZ_tautaull, mA, mH);
        CacheShiftReal(ip_ex_pp_A_HZ_tautaull_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_AZ_tautaull_CMS8(double mA, double mH){
    int NumPar = 2;
    double params[] = {mA, mH};

    int i = CacheCheckReal(ip_ex_pp_H_AZ_tautaull_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_H_AZ_tautaull_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate2D(CMS8_pp_H_AZ_tautaull, mA, mH);
        CacheShiftReal(ip_ex_pp_H_AZ_tautaull_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_bb_phi_tt_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_tt_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phi_tt_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_bb_phi_tt,mass);
        CacheShiftReal(ip_ex_bb_phi_tt_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_bb_phi_tt_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_tt_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phi_tt_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_bb_phi_tt_e,mass);
        CacheShiftReal(ip_ex_bb_phi_tt_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_tt_phi_tt_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_tt_phi_tt_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_tt_phi_tt_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_tt_phi_tt,mass);
        CacheShiftReal(ip_ex_tt_phi_tt_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_tt_phi_tt_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_tt_phi_tt_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_tt_phi_tt_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_tt_phi_tt_e,mass);
        CacheShiftReal(ip_ex_tt_phi_tt_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_phi_tautau_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_tautau_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_phi_tautau_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_gg_phi_tautau,mass);
        CacheShiftReal(ip_ex_gg_phi_tautau_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_phi_tautau_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_tautau_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_phi_tautau_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_gg_phi_tautau_e,mass);
        CacheShiftReal(ip_ex_gg_phi_tautau_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_bb_phi_tautau_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_tautau_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_bb_phi_tautau_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_bb_phi_tautau,mass);
        CacheShiftReal(ip_ex_bb_phi_tautau_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_bb_phi_tautau_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_tautau_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_bb_phi_tautau_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_bb_phi_tautau_e,mass);
        CacheShiftReal(ip_ex_bb_phi_tautau_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_phi_gaga_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_gaga_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_gaga_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_phi_gaga,mass);
        CacheShiftReal(ip_ex_pp_phi_gaga_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_phi_gaga_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_gaga_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_gaga_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_phi_gaga_e,mass);
        CacheShiftReal(ip_ex_pp_phi_gaga_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_phi_Zga_llga_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_Zga_llga_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_Zga_llga_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_phi_Zga,mass);
        CacheShiftReal(ip_ex_pp_phi_Zga_llga_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_phi_Zga_llga_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_Zga_llga_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_Zga_llga_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_phi_Zga_e,mass);
        CacheShiftReal(ip_ex_pp_phi_Zga_llga_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_phi_Zga_llga_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_Zga_llga_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_phi_Zga_llga_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_gg_phi_Zga_llga,mass);
        CacheShiftReal(ip_ex_gg_phi_Zga_llga_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_phi_Zga_llga_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_Zga_llga_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_phi_Zga_llga_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_gg_phi_Zga_llga_e,mass);
        CacheShiftReal(ip_ex_gg_phi_Zga_llga_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_H_ZZ_llllnunu_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_H_ZZ_llllnunu_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_H_ZZ_llllnunu_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_gg_H_ZZ_llllnunu,mass);
        CacheShiftReal(ip_ex_gg_H_ZZ_llllnunu_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_H_ZZ_llllnunu_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_H_ZZ_llllnunu_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_H_ZZ_llllnunu_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_gg_H_ZZ_llllnunu_e,mass);
        CacheShiftReal(ip_ex_gg_H_ZZ_llllnunu_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_VBF_H_ZZ_llllnunu,mass);
        CacheShiftReal(ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_VBF_H_ZZ_llllnunu_e,mass);
        CacheShiftReal(ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_H_ZZ_llnunu_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_H_ZZ_llnunu_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_H_ZZ_llnunu_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_gg_H_ZZ_llnunu,mass);
        CacheShiftReal(ip_ex_gg_H_ZZ_llnunu_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_H_ZZ_llnunu_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_H_ZZ_llnunu_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_H_ZZ_llnunu_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_gg_H_ZZ_llnunu_e,mass);
        CacheShiftReal(ip_ex_gg_H_ZZ_llnunu_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_H_ZZ_llll_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_H_ZZ_llll_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_H_ZZ_llll_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_gg_H_ZZ_llll,mass);
        CacheShiftReal(ip_ex_gg_H_ZZ_llll_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_H_ZZ_llll_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_H_ZZ_llll_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_H_ZZ_llll_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_gg_H_ZZ_llll_e,mass);
        CacheShiftReal(ip_ex_gg_H_ZZ_llll_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_VBF_H_ZZ_llll_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VBF_H_ZZ_llll_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_VBF_H_ZZ_llll_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_VBF_H_ZZ_llll,mass);
        CacheShiftReal(ip_ex_VBF_H_ZZ_llll_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_VBF_H_ZZ_llll_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VBF_H_ZZ_llll_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_VBF_H_ZZ_llll_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_VBF_H_ZZ_llll_e,mass);
        CacheShiftReal(ip_ex_VBF_H_ZZ_llll_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_gg_H_ZZ_qqllnunu,mass);
        CacheShiftReal(ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_gg_H_ZZ_qqllnunu_e,mass);
        CacheShiftReal(ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_VBF_H_ZZ_qqllnunu,mass);
        CacheShiftReal(ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_VBF_H_ZZ_qqllnunu_e,mass);
        CacheShiftReal(ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_H_ZZ_llqq_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_H_ZZ_llqq_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_H_ZZ_llqq_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_gg_H_ZZ_llqq,mass);
        CacheShiftReal(ip_ex_gg_H_ZZ_llqq_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_H_ZZ_llqq_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_H_ZZ_llqq_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_H_ZZ_llqq_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_gg_H_ZZ_llqq_e,mass);
        CacheShiftReal(ip_ex_gg_H_ZZ_llqq_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_VBF_H_ZZ_llqq_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VBF_H_ZZ_llqq_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_VBF_H_ZZ_llqq_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_VBF_H_ZZ_llqq,mass);
        CacheShiftReal(ip_ex_VBF_H_ZZ_llqq_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_VBF_H_ZZ_llqq_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VBF_H_ZZ_llqq_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_VBF_H_ZZ_llqq_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_VBF_H_ZZ_llqq_e,mass);
        CacheShiftReal(ip_ex_VBF_H_ZZ_llqq_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_H_ZZ_nunuqq_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_H_ZZ_nunuqq_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_H_ZZ_nunuqq_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_gg_H_ZZ_nunuqq,mass);
        CacheShiftReal(ip_ex_gg_H_ZZ_nunuqq_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_H_ZZ_nunuqq_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_H_ZZ_nunuqq_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_H_ZZ_nunuqq_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_gg_H_ZZ_nunuqq_e,mass);
        CacheShiftReal(ip_ex_gg_H_ZZ_nunuqq_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_H_WW_enumunu_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_H_WW_enumunu_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_H_WW_enumunu_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_gg_H_WW_enumumu,mass);
        CacheShiftReal(ip_ex_gg_H_WW_enumunu_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_H_WW_enumunu_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_H_WW_enumunu_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_H_WW_enumunu_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_gg_H_WW_enumumu_e,mass);
        CacheShiftReal(ip_ex_gg_H_WW_enumunu_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_VBF_H_WW_enumunu_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VBF_H_WW_enumunu_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_VBF_H_WW_enumunu_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_VBF_H_WW_enumumu,mass);
        CacheShiftReal(ip_ex_VBF_H_WW_enumunu_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_VBF_H_WW_enumunu_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VBF_H_WW_enumunu_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_VBF_H_WW_enumunu_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_VBF_H_WW_enumumu_e,mass);
        CacheShiftReal(ip_ex_VBF_H_WW_enumunu_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_H_WW_lnuqq_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_H_WW_lnuqq_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_H_WW_lnuqq_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_gg_H_WW_lnuqq,mass);
        CacheShiftReal(ip_ex_gg_H_WW_lnuqq_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_H_WW_lnuqq_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_H_WW_lnuqq_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_H_WW_lnuqq_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_gg_H_WW_lnuqq_e,mass);
        CacheShiftReal(ip_ex_gg_H_WW_lnuqq_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_VBF_H_WW_lnuqq_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VBF_H_WW_lnuqq_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_VBF_H_WW_lnuqq_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_VBF_H_WW_lnuqq,mass);
        CacheShiftReal(ip_ex_VBF_H_WW_lnuqq_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_VBF_H_WW_lnuqq_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VBF_H_WW_lnuqq_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_VBF_H_WW_lnuqq_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_VBF_H_WW_lnuqq_e,mass);
        CacheShiftReal(ip_ex_VBF_H_WW_lnuqq_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_VV_qqqq_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_VV_qqqq_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_VV_qqqq_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_H_VV_qqqq,mass);
        CacheShiftReal(ip_ex_pp_H_VV_qqqq_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_VV_qqqq_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_VV_qqqq_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_VV_qqqq_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_H_VV_qqqq_e,mass);
        CacheShiftReal(ip_ex_pp_H_VV_qqqq_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_hh_bbbb_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_hh_bbbb_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_hh_bbbb_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_H_hh_bbbb,mass);
        CacheShiftReal(ip_ex_pp_H_hh_bbbb_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_hh_bbbb_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_hh_bbbb_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_hh_bbbb_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_H_hh_bbbb_e,mass);
        CacheShiftReal(ip_ex_pp_H_hh_bbbb_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_hh_gagabb_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_hh_gagabb_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_hh_gagabb_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_H_hh_gagabb,mass);
        CacheShiftReal(ip_ex_pp_H_hh_gagabb_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_hh_gagabb_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_hh_gagabb_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_hh_gagabb_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_H_hh_gagabb_e,mass);
        CacheShiftReal(ip_ex_pp_H_hh_gagabb_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_hh_gagaWW_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_hh_gagaWW_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_hh_gagaWW_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_H_hh_gagaWW,mass);
        CacheShiftReal(ip_ex_pp_H_hh_gagaWW_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_hh_gagaWW_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_hh_gagaWW_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_hh_gagaWW_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_H_hh_gagaWW_e,mass);
        CacheShiftReal(ip_ex_pp_H_hh_gagaWW_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_A_Zh_Zbb_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_A_Zh_Zbb_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_A_Zh_Zbb_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_gg_A_Zh_Zbb,mass);
        CacheShiftReal(ip_ex_gg_A_Zh_Zbb_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_A_Zh_Zbb_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_A_Zh_Zbb_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_A_Zh_Zbb_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_gg_A_Zh_Zbb_e,mass);
        CacheShiftReal(ip_ex_gg_A_Zh_Zbb_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_bb_A_Zh_Zbb_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_A_Zh_Zbb_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_bb_A_Zh_Zbb_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_bb_A_Zh_Zbb,mass);
        CacheShiftReal(ip_ex_bb_A_Zh_Zbb_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_bb_A_Zh_Zbb_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_A_Zh_Zbb_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_bb_A_Zh_Zbb_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_bb_A_Zh_Zbb_e,mass);
        CacheShiftReal(ip_ex_bb_A_Zh_Zbb_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_phi_bb_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_bb_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_bb_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_phi_bb,mass);
        CacheShiftReal(ip_ex_pp_phi_bb_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_phi_bb_CMS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_bb_CMS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_bb_CMS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_phi_bb_e,mass);
        CacheShiftReal(ip_ex_pp_phi_bb_CMS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_phi_tautau_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_tautau_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_phi_tautau_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_gg_phi_tautau,mass);
        CacheShiftReal(ip_ex_gg_phi_tautau_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_phi_tautau_CMS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_tautau_CMS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_phi_tautau_CMS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_gg_phi_tautau_e,mass);
        CacheShiftReal(ip_ex_gg_phi_tautau_CMS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_bb_phi_tautau_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_tautau_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_bb_phi_tautau_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_bb_phi_tautau,mass);
        CacheShiftReal(ip_ex_bb_phi_tautau_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_bb_phi_tautau_CMS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_tautau_CMS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_bb_phi_tautau_CMS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_bb_phi_tautau_e,mass);
        CacheShiftReal(ip_ex_bb_phi_tautau_CMS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_phi_gaga_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_gaga_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_phi_gaga_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_gg_phi_gaga,mass);
        CacheShiftReal(ip_ex_gg_phi_gaga_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_phi_gaga_CMS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_gaga_CMS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_phi_gaga_CMS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_gg_phi_gaga_e,mass);
        CacheShiftReal(ip_ex_gg_phi_gaga_CMS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_phi_Zga_llga_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_Zga_llga_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_Zga_llga_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_phi_Zga_llga,mass);
        CacheShiftReal(ip_ex_pp_phi_Zga_llga_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_phi_Zga_llga_CMS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_Zga_llga_CMS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_Zga_llga_CMS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_phi_Zga_llga_e,mass);
        CacheShiftReal(ip_ex_pp_phi_Zga_llga_CMS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_phi_Zga_qqga_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_Zga_qqga_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_Zga_qqga_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_phi_Zga_qqga,mass);
        CacheShiftReal(ip_ex_pp_phi_Zga_qqga_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_phi_Zga_qqga_CMS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_Zga_qqga_CMS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_Zga_qqga_CMS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_phi_Zga_qqga_e,mass);
        CacheShiftReal(ip_ex_pp_phi_Zga_qqga_CMS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_ggF_phi_Zga_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_ggF_phi_Zga_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_ggF_phi_Zga_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_ggF_phi_Zga,mass);
        CacheShiftReal(ip_ex_ggF_phi_Zga_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_ggF_phi_Zga_CMS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_ggF_phi_Zga_CMS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_ggF_phi_Zga_CMS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_ggF_phi_Zga_e,mass);
        CacheShiftReal(ip_ex_ggF_phi_Zga_CMS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_ZZ_llnunu_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_ZZ_llnunu_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_ZZ_llnunu_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_H_ZZ_llnunu,mass);
        CacheShiftReal(ip_ex_pp_H_ZZ_llnunu_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_ZZ_llnunu_CMS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_ZZ_llnunu_CMS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_ZZ_llnunu_CMS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_H_ZZ_llnunu_e,mass);
        CacheShiftReal(ip_ex_pp_H_ZZ_llnunu_CMS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_H_ZZ_llnunu_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_H_ZZ_llnunu_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_H_ZZ_llnunu_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_gg_H_ZZ_llnunu,mass);
        CacheShiftReal(ip_ex_gg_H_ZZ_llnunu_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_gg_H_ZZ_llnunu_CMS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_H_ZZ_llnunu_CMS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_gg_H_ZZ_llnunu_CMS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_gg_H_ZZ_llnunu_e,mass);
        CacheShiftReal(ip_ex_gg_H_ZZ_llnunu_CMS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_VBF_H_ZZ_llnunu_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VBF_H_ZZ_llnunu_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_VBF_H_ZZ_llnunu_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_VBF_H_ZZ_llnunu,mass);
        CacheShiftReal(ip_ex_VBF_H_ZZ_llnunu_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_VBF_H_ZZ_llnunu_CMS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VBF_H_ZZ_llnunu_CMS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_VBF_H_ZZ_llnunu_CMS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_VBF_H_ZZ_llnunu_e,mass);
        CacheShiftReal(ip_ex_VBF_H_ZZ_llnunu_CMS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_ZZ_llll_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_ZZ_llll_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_ZZ_llll_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_H_ZZ_llll,mass);
        CacheShiftReal(ip_ex_pp_H_ZZ_llll_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_ZZ_llll_CMS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_ZZ_llll_CMS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_ZZ_llll_CMS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_H_ZZ_llll_e,mass);
        CacheShiftReal(ip_ex_pp_H_ZZ_llll_CMS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_VBF_VH_H_ZZ_llll_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VBF_VH_H_ZZ_llll_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_VBF_VH_H_ZZ_llll_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_VBFVH_H_ZZ_llll,mass);
        CacheShiftReal(ip_ex_VBF_VH_H_ZZ_llll_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VBF_VH_H_ZZ_llll_CMS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_VBF_VH_H_ZZ_llll_CMS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_VBFVH_H_ZZ_llll_e,mass);
        CacheShiftReal(ip_ex_VBF_VH_H_ZZ_llll_CMS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_ZZ_llqq_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_ZZ_llqq_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_ZZ_llqq_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_H_ZZ_llqq,mass);
        CacheShiftReal(ip_ex_pp_H_ZZ_llqq_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_ZZ_llqq_CMS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_ZZ_llqq_CMS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_ZZ_llqq_CMS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_H_ZZ_llqq_e,mass);
        CacheShiftReal(ip_ex_pp_H_ZZ_llqq_CMS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_ggVV_H_WW_lnulnu_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_ggVV_H_WW_lnulnu_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_ggVV_H_WW_lnulnu_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_ggFVBF_H_WW_lnulnu,mass);
        CacheShiftReal(ip_ex_ggVV_H_WW_lnulnu_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_ggVV_H_WW_lnulnu_CMS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_ggVV_H_WW_lnulnu_CMS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_ggVV_H_WW_lnulnu_CMS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_ggFVBF_H_WW_lnulnu_e,mass);
        CacheShiftReal(ip_ex_ggVV_H_WW_lnulnu_CMS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_hh_bbbb_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_hh_bbbb_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_hh_bbbb_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_H_hh_bbbb,mass);
        CacheShiftReal(ip_ex_pp_H_hh_bbbb_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_hh_bbbb_CMS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_hh_bbbb_CMS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_hh_bbbb_CMS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_H_hh_bbbb_e,mass);
        CacheShiftReal(ip_ex_pp_H_hh_bbbb_CMS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_ggF_H_hh_bbbb_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_ggF_H_hh_bbbb_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_ggF_H_hh_bbbb_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_ggF_H_hh_bbbb,mass);
        CacheShiftReal(ip_ex_ggF_H_hh_bbbb_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_ggF_H_hh_bbbb_CMS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_ggF_H_hh_bbbb_CMS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_ggF_H_hh_bbbb_CMS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_ggF_H_hh_bbbb_e,mass);
        CacheShiftReal(ip_ex_ggF_H_hh_bbbb_CMS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_hh_gagabb_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_hh_gagabb_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_hh_gagabb_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_H_hh_gagabb,mass);
        CacheShiftReal(ip_ex_pp_H_hh_gagabb_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_hh_gagabb_CMS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_hh_gagabb_CMS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_hh_gagabb_CMS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_H_hh_gagabb_e,mass);
        CacheShiftReal(ip_ex_pp_H_hh_gagabb_CMS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_hh_bbtautau_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_hh_bbtautau_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_hh_bbtautau_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_H_hh_bbtautau,mass);
        CacheShiftReal(ip_ex_pp_H_hh_bbtautau_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_hh_bbtautau_CMS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_hh_bbtautau_CMS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_hh_bbtautau_CMS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_H_hh_bbtautau_e,mass);
        CacheShiftReal(ip_ex_pp_H_hh_bbtautau_CMS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_hh_bbtautau1_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_hh_bbtautau1_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_hh_bbtautau1_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_H_hh_bbtautau1,mass);
        CacheShiftReal(ip_ex_pp_H_hh_bbtautau1_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_hh_bbtautau1_CMS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_hh_bbtautau1_CMS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_hh_bbtautau1_CMS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_H_hh_bbtautau1_e,mass);
        CacheShiftReal(ip_ex_pp_H_hh_bbtautau1_CMS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_hh_bblnulnu_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_hh_bblnulnu_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_hh_bblnulnu_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_H_hh_bblnulnu,mass);
        CacheShiftReal(ip_ex_pp_H_hh_bblnulnu_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_hh_bblnulnu_CMS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_hh_bblnulnu_CMS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_hh_bblnulnu_CMS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_H_hh_bblnulnu_e,mass);
        CacheShiftReal(ip_ex_pp_H_hh_bblnulnu_CMS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_hh_bbVV_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_hh_bbVV_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_hh_bbVV_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_H_hh_bbVV,mass);
        CacheShiftReal(ip_ex_pp_H_hh_bbVV_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_H_hh_bbVV_CMS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_hh_bbVV_CMS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_hh_bbVV_CMS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_H_hh_bbVV_e,mass);
        CacheShiftReal(ip_ex_pp_H_hh_bbVV_CMS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie1(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie1c, NumPar, params);
    if (i>=0) {
        return(ie1c[NumPar][i] );
    } else {
        double newResult = interpolate (temp1,mass);
        CacheShiftReal(ie1c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie1e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie1ec, NumPar, params);
    if (i>=0) {
        return(ie1ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp1e,mass);
        CacheShiftReal(ie1ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie2(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie2c, NumPar, params);
    if (i>=0) {
        return(ie2c[NumPar][i] );
    } else {
        double newResult = interpolate (temp2,mass);
        CacheShiftReal(ie2c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie2e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie2ec, NumPar, params);
    if (i>=0) {
        return(ie2ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp2e,mass);
        CacheShiftReal(ie2ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie3(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie3c, NumPar, params);
    if (i>=0) {
        return(ie3c[NumPar][i] );
    } else {
        double newResult = interpolate (temp3,mass);
        CacheShiftReal(ie3c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie3e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie3ec, NumPar, params);
    if (i>=0) {
        return(ie3ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp3e,mass);
        CacheShiftReal(ie3ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie4(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie4c, NumPar, params);
    if (i>=0) {
        return(ie4c[NumPar][i] );
    } else {
        double newResult = interpolate (temp4,mass);
        CacheShiftReal(ie4c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie4e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie4ec, NumPar, params);
    if (i>=0) {
        return(ie4ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp4e,mass);
        CacheShiftReal(ie4ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie5(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie5c, NumPar, params);
    if (i>=0) {
        return(ie5c[NumPar][i] );
    } else {
        double newResult = interpolate (temp5,mass);
        CacheShiftReal(ie5c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie5e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie5ec, NumPar, params);
    if (i>=0) {
        return(ie5ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp5e,mass);
        CacheShiftReal(ie5ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie6(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie6c, NumPar, params);
    if (i>=0) {
        return(ie6c[NumPar][i] );
    } else {
        double newResult = interpolate (temp6,mass);
        CacheShiftReal(ie6c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie6e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie6ec, NumPar, params);
    if (i>=0) {
        return(ie6ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp6e,mass);
        CacheShiftReal(ie6ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie7(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie7c, NumPar, params);
    if (i>=0) {
        return(ie7c[NumPar][i] );
    } else {
        double newResult = interpolate (temp7,mass);
        CacheShiftReal(ie7c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie7e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie7ec, NumPar, params);
    if (i>=0) {
        return(ie7ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp7e,mass);
        CacheShiftReal(ie7ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie8c, NumPar, params);
    if (i>=0) {
        return(ie8c[NumPar][i] );
    } else {
        double newResult = interpolate (temp8,mass);
        CacheShiftReal(ie8c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie8e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie8ec, NumPar, params);
    if (i>=0) {
        return(ie8ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp8e,mass);
        CacheShiftReal(ie8ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie9(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie9c, NumPar, params);
    if (i>=0) {
        return(ie9c[NumPar][i] );
    } else {
        double newResult = interpolate (temp9,mass);
        CacheShiftReal(ie9c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie9e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie9ec, NumPar, params);
    if (i>=0) {
        return(ie9ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp9e,mass);
        CacheShiftReal(ie9ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie10(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie10c, NumPar, params);
    if (i>=0) {
        return(ie10c[NumPar][i] );
    } else {
        double newResult = interpolate (temp10,mass);
        CacheShiftReal(ie10c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie10e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie10ec, NumPar, params);
    if (i>=0) {
        return(ie10ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp10e,mass);
        CacheShiftReal(ie10ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie11(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie11c, NumPar, params);
    if (i>=0) {
        return(ie11c[NumPar][i] );
    } else {
        double newResult = interpolate (temp11,mass);
        CacheShiftReal(ie11c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie11e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie11ec, NumPar, params);
    if (i>=0) {
        return(ie11ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp11e,mass);
        CacheShiftReal(ie11ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie12(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie12c, NumPar, params);
    if (i>=0) {
        return(ie12c[NumPar][i] );
    } else {
        double newResult = interpolate (temp12,mass);
        CacheShiftReal(ie12c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie12e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie12ec, NumPar, params);
    if (i>=0) {
        return(ie12ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp12e,mass);
        CacheShiftReal(ie12ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie13c, NumPar, params);
    if (i>=0) {
        return(ie13c[NumPar][i] );
    } else {
        double newResult = interpolate (temp13,mass);
        CacheShiftReal(ie13c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie13e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie13ec, NumPar, params);
    if (i>=0) {
        return(ie13ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp13e,mass);
        CacheShiftReal(ie13ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie14(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie14c, NumPar, params);
    if (i>=0) {
        return(ie14c[NumPar][i] );
    } else {
        double newResult = interpolate (temp14,mass);
        CacheShiftReal(ie14c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie14e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie14ec, NumPar, params);
    if (i>=0) {
        return(ie14ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp14e,mass);
        CacheShiftReal(ie14ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie15(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie15c, NumPar, params);
    if (i>=0) {
        return(ie15c[NumPar][i] );
    } else {
        double newResult = interpolate (temp15,mass);
        CacheShiftReal(ie15c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie15e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie15ec, NumPar, params);
    if (i>=0) {
        return(ie15ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp15e,mass);
        CacheShiftReal(ie15ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie16(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie16c, NumPar, params);
    if (i>=0) {
        return(ie16c[NumPar][i] );
    } else {
        double newResult = interpolate (temp16,mass);
        CacheShiftReal(ie16c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie16e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie16ec, NumPar, params);
    if (i>=0) {
        return(ie16ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp16e,mass);
        CacheShiftReal(ie16ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie17(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie17c, NumPar, params);
    if (i>=0) {
        return(ie17c[NumPar][i] );
    } else {
        double newResult = interpolate (temp17,mass);
        CacheShiftReal(ie17c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie17e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie17ec, NumPar, params);
    if (i>=0) {
        return(ie17ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp17e,mass);
        CacheShiftReal(ie17ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie18(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie18c, NumPar, params);
    if (i>=0) {
        return(ie18c[NumPar][i] );
    } else {
        double newResult = interpolate (temp18,mass);
        CacheShiftReal(ie18c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie18e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie18ec, NumPar, params);
    if (i>=0) {
        return(ie18ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp18e,mass);
        CacheShiftReal(ie18ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie19(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie19c, NumPar, params);
    if (i>=0) {
        return(ie19c[NumPar][i] );
    } else {
        double newResult = interpolate (temp19,mass);
        CacheShiftReal(ie19c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie19e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie19ec, NumPar, params);
    if (i>=0) {
        return(ie19ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp19e,mass);
        CacheShiftReal(ie19ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie20(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie20c, NumPar, params);
    if (i>=0) {
        return(ie20c[NumPar][i] );
    } else {
        double newResult = interpolate (temp20,mass);
        CacheShiftReal(ie20c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie20e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie20ec, NumPar, params);
    if (i>=0) {
        return(ie20ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp20e,mass);
        CacheShiftReal(ie20ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie21(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie21c, NumPar, params);
    if (i>=0) {
        return(ie21c[NumPar][i] );
    } else {
        double newResult = interpolate (temp21,mass);
        CacheShiftReal(ie21c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie21e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie21ec, NumPar, params);
    if (i>=0) {
        return(ie21ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp21e,mass);
        CacheShiftReal(ie21ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie22(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie22c, NumPar, params);
    if (i>=0) {
        return(ie22c[NumPar][i] );
    } else {
        double newResult = interpolate (temp22,mass);
        CacheShiftReal(ie22c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie22e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie22ec, NumPar, params);
    if (i>=0) {
        return(ie22ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp22e,mass);
        CacheShiftReal(ie22ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie23(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie23c, NumPar, params);
    if (i>=0) {
        return(ie23c[NumPar][i] );
    } else {
        double newResult = interpolate (temp23,mass);
        CacheShiftReal(ie23c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie23e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie23ec, NumPar, params);
    if (i>=0) {
        return(ie23ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp23e,mass);
        CacheShiftReal(ie23ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie24(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie24c, NumPar, params);
    if (i>=0) {
        return(ie24c[NumPar][i] );
    } else {
        double newResult = interpolate (temp24,mass);
        CacheShiftReal(ie24c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie24e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie24ec, NumPar, params);
    if (i>=0) {
        return(ie24ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp24e,mass);
        CacheShiftReal(ie24ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie25(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie25c, NumPar, params);
    if (i>=0) {
        return(ie25c[NumPar][i] );
    } else {
        double newResult = interpolate (temp25,mass);
        CacheShiftReal(ie25c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie25e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie25ec, NumPar, params);
    if (i>=0) {
        return(ie25ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp25e,mass);
        CacheShiftReal(ie25ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie26(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie26c, NumPar, params);
    if (i>=0) {
        return(ie26c[NumPar][i] );
    } else {
        double newResult = interpolate (temp26,mass);
        CacheShiftReal(ie26c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie26e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie26ec, NumPar, params);
    if (i>=0) {
        return(ie26ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp26e,mass);
        CacheShiftReal(ie26ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie27(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie27c, NumPar, params);
    if (i>=0) {
        return(ie27c[NumPar][i] );
    } else {
        double newResult = interpolate (temp27,mass);
        CacheShiftReal(ie27c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie27e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie27ec, NumPar, params);
    if (i>=0) {
        return(ie27ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp27e,mass);
        CacheShiftReal(ie27ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie28(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie28c, NumPar, params);
    if (i>=0) {
        return(ie28c[NumPar][i] );
    } else {
        double newResult = interpolate (temp28,mass);
        CacheShiftReal(ie28c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie28e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie28ec, NumPar, params);
    if (i>=0) {
        return(ie28ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp28e,mass);
        CacheShiftReal(ie28ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie29(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie29c, NumPar, params);
    if (i>=0) {
        return(ie29c[NumPar][i] );
    } else {
        double newResult = interpolate (temp29,mass);
        CacheShiftReal(ie29c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie29e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie29ec, NumPar, params);
    if (i>=0) {
        return(ie29ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp29e,mass);
        CacheShiftReal(ie29ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie30(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie30c, NumPar, params);
    if (i>=0) {
        return(ie30c[NumPar][i] );
    } else {
        double newResult = interpolate (temp30,mass);
        CacheShiftReal(ie30c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie30e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie30ec, NumPar, params);
    if (i>=0) {
        return(ie30ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp30e,mass);
        CacheShiftReal(ie30ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie31(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie31c, NumPar, params);
    if (i>=0) {
        return(ie31c[NumPar][i] );
    } else {
        double newResult = interpolate (temp31,mass);
        CacheShiftReal(ie31c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie31e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie31ec, NumPar, params);
    if (i>=0) {
        return(ie31ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp31e,mass);
        CacheShiftReal(ie31ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie32(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie32c, NumPar, params);
    if (i>=0) {
        return(ie32c[NumPar][i] );
    } else {
        double newResult = interpolate (temp32,mass);
        CacheShiftReal(ie32c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie32e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie32ec, NumPar, params);
    if (i>=0) {
        return(ie32ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp32e,mass);
        CacheShiftReal(ie32ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie33(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie33c, NumPar, params);
    if (i>=0) {
        return(ie33c[NumPar][i] );
    } else {
        double newResult = interpolate (temp33,mass);
        CacheShiftReal(ie33c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie33e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie33ec, NumPar, params);
    if (i>=0) {
        return(ie33ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp33e,mass);
        CacheShiftReal(ie33ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie34(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie34c, NumPar, params);
    if (i>=0) {
        return(ie34c[NumPar][i] );
    } else {
        double newResult = interpolate (temp34,mass);
        CacheShiftReal(ie34c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie34e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie34ec, NumPar, params);
    if (i>=0) {
        return(ie34ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp34e,mass);
        CacheShiftReal(ie34ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie35(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie35c, NumPar, params);
    if (i>=0) {
        return(ie35c[NumPar][i] );
    } else {
        double newResult = interpolate (temp35,mass);
        CacheShiftReal(ie35c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie35e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie35ec, NumPar, params);
    if (i>=0) {
        return(ie35ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp35e,mass);
        CacheShiftReal(ie35ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie36(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie36c, NumPar, params);
    if (i>=0) {
        return(ie36c[NumPar][i] );
    } else {
        double newResult = interpolate (temp36,mass);
        CacheShiftReal(ie36c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie36e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie36ec, NumPar, params);
    if (i>=0) {
        return(ie36ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp36e,mass);
        CacheShiftReal(ie36ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie37(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie37c, NumPar, params);
    if (i>=0) {
        return(ie37c[NumPar][i] );
    } else {
        double newResult = interpolate (temp37,mass);
        CacheShiftReal(ie37c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie37e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie37ec, NumPar, params);
    if (i>=0) {
        return(ie37ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp37e,mass);
        CacheShiftReal(ie37ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie38(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie38c, NumPar, params);
    if (i>=0) {
        return(ie38c[NumPar][i] );
    } else {
        double newResult = interpolate (temp38,mass);
        CacheShiftReal(ie38c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie38e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie38ec, NumPar, params);
    if (i>=0) {
        return(ie38ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp38e,mass);
        CacheShiftReal(ie38ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie39(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie39c, NumPar, params);
    if (i>=0) {
        return(ie39c[NumPar][i] );
    } else {
        double newResult = interpolate (temp39,mass);
        CacheShiftReal(ie39c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie39e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie39ec, NumPar, params);
    if (i>=0) {
        return(ie39ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp39e,mass);
        CacheShiftReal(ie39ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie40(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie40c, NumPar, params);
    if (i>=0) {
        return(ie40c[NumPar][i] );
    } else {
        double newResult = interpolate (temp40,mass);
        CacheShiftReal(ie40c, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ie40e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ie40ec, NumPar, params);
    if (i>=0) {
        return(ie40ec[NumPar][i] );
    } else {
        double newResult = interpolate (temp40e,mass);
        CacheShiftReal(ie40ec, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_Hpm_taunu_ATLAS8(double mass){
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



double THDMcache::ip_ex_pp_Hpm_taunu_ATLAS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_Hpm_taunu_ATLAS8_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_Hpm_taunu_ATLAS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS8_pp_Hpm_taunu_e,mass);
        CacheShiftReal(ip_ex_pp_Hpm_taunu_ATLAS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_Hp_taunu_CMS8(double mass){
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



double THDMcache::ip_ex_pp_Hp_taunu_CMS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_Hp_taunu_CMS8_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_Hp_taunu_CMS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_pp_Hp_taunu_e,mass);
        CacheShiftReal(ip_ex_pp_Hp_taunu_CMS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_Hpm_tb_ATLAS8(double mass){
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



double THDMcache::ip_ex_pp_Hpm_tb_ATLAS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_Hpm_tb_ATLAS8_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_Hpm_tb_ATLAS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS8_pp_Hpm_tb_e,mass);
        CacheShiftReal(ip_ex_pp_Hpm_tb_ATLAS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_Hp_tb_CMS8(double mass){
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



double THDMcache::ip_ex_pp_Hp_tb_CMS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_Hp_tb_CMS8_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_Hp_tb_CMS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_pp_Hp_tb_e,mass);
        CacheShiftReal(ip_ex_pp_Hp_tb_CMS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_Hpm_taunu_ATLAS13(double mass){
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



double THDMcache::ip_ex_pp_Hpm_taunu_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_Hpm_taunu_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_Hpm_taunu_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_Hpm_taunu_e,mass);
        CacheShiftReal(ip_ex_pp_Hpm_taunu_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_Hpm_taunu_CMS13(double mass){
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



double THDMcache::ip_ex_pp_Hpm_taunu_CMS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_Hpm_taunu_CMS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_Hpm_taunu_CMS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_Hpm_taunu_e,mass);
        CacheShiftReal(ip_ex_pp_Hpm_taunu_CMS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_Hp_tb_ATLAS13_1(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_Hp_tb_ATLAS13_1_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_Hp_tb_ATLAS13_1_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_Hp_tb1,mass);
        CacheShiftReal(ip_ex_pp_Hp_tb_ATLAS13_1_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_Hp_tb_ATLAS13_1_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_Hp_tb_ATLAS13_1_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_Hp_tb_ATLAS13_1_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_Hp_tb1_e,mass);
        CacheShiftReal(ip_ex_pp_Hp_tb_ATLAS13_1_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_Hp_tb_ATLAS13_2(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_Hp_tb_ATLAS13_2_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_Hp_tb_ATLAS13_2_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_Hp_tb2,mass);
        CacheShiftReal(ip_ex_pp_Hp_tb_ATLAS13_2_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_pp_Hp_tb_ATLAS13_2_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_Hp_tb_ATLAS13_2_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_Hp_tb_ATLAS13_2_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_Hp_tb2_e,mass);
        CacheShiftReal(ip_ex_pp_Hp_tb_ATLAS13_2_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMcache::ip_ex_bsgamma(double logtb, double logmHp){
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



gslpp::matrix<double> THDMcache::readTable(std::string filename, int rowN, int colN){

    std::ifstream INfile;
    std::string lineTab;
    INfile.open( filename.c_str() );
    if(INfile.fail()){
        std::cout<<"error: in THDMcache, table doesn't exist!"<<std::endl;
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

double THDMcache::interpolate2D(gslpp::matrix<double> arrayTab, double x, double y){

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

double THDMcache::ghHpHm(const double mHp2, const double tanb, const double m12_2, const double bma, const double mHl2, const double vev) const {
    int NumPar = 6;
    double params[] = {mHp2, tanb, m12_2, bma, mHl2, vev};

    int i = CacheCheckReal(ghHpHm_cache, NumPar, params);
    if (i>=0) {
        return ( ghHpHm_cache[NumPar][i] );
    } else {
        double newResult = ((cos(bma)*mHl2*(tanb*tanb-1.0))/tanb 
                                    -(mHl2+2.0*mHp2)*sin(bma) 
                                    +(m12_2*(cos(bma)*(1.0-tanb*tanb)+2.0*sin(bma)*tanb)*(1.0+tanb*tanb))/(tanb*tanb))/vev;
        CacheShiftReal(ghHpHm_cache, NumPar, params, newResult);
        return newResult;
    }
}

double THDMcache::g_HH_HpHm(const double mHp2, const double mHh2, const double tanb, const double m12_2, const double bma, const double vev) const {
    int NumPar = 6;
    double params[] = {mHp2, mHh2, tanb, m12_2, bma, vev};

    int i = CacheCheckReal(g_HH_HpHm_cache, NumPar, params);
    if (i>=0) {
        return ( g_HH_HpHm_cache[NumPar][i] );
    } else {
        double newResult = (cos(bma)*(mHh2-2.0*mHp2)
                                    +((m12_2-mHh2*tanb+m12_2*tanb*tanb)
                                      *(2.0*cos(bma)*tanb+sin(bma)*(tanb*tanb-1.0)))/(tanb*tanb))/vev;
        CacheShiftReal(g_HH_HpHm_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex THDMcache::I_h_U(const double mHl2, const double Mu, const double Mc, const double Mt) const {
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

gslpp::complex THDMcache::I_HH_U(const double mHh2, const double Mc, const double Mt) const {
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

gslpp::complex THDMcache::I_A_U(const double mA2, const double Mc, const double Mt) const {
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

gslpp::complex THDMcache::I_h_D(const double mHl2, const double Md, const double Ms, const double Mb) const {
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

gslpp::complex THDMcache::I_HH_D(const double mHh2, const double Ms, const double Mb) const {
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

gslpp::complex THDMcache::I_A_D(const double mA2, const double Ms, const double Mb) const {
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

gslpp::complex THDMcache::I_h_L(const double mHl2, const double Me, const double Mmu, const double Mtau) const {
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

gslpp::complex THDMcache::I_HH_L(const double mHh2, const double Mmu, const double Mtau) const {
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

gslpp::complex THDMcache::I_A_L(const double mA2, const double Mmu, const double Mtau) const {
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

gslpp::complex THDMcache::I_H_W(const double mH, const double MW) const {
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

gslpp::complex THDMcache::I_H_Hp(const double mHp2, const double mH) const {
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

gslpp::complex THDMcache::A_h_U(const double mHl2, const double cW2, const double Mu, const double Mc, const double Mt, const double MZ) const {
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

gslpp::complex THDMcache::A_HH_U(const double mHh2, const double cW2, const double Mc, const double Mt, const double MZ) const {
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

gslpp::complex THDMcache::A_A_U(const double mA2, const double cW2, const double Mc, const double Mt, const double MZ) const {
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

gslpp::complex THDMcache::A_h_D(const double mHl2, const double cW2, const double Md, const double Ms, const double Mb, const double MZ) const {
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

gslpp::complex THDMcache::A_HH_D(const double mHh2, const double cW2, const double Ms, const double Mb, const double MZ) const {
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

gslpp::complex THDMcache::A_A_D(const double mA2, const double cW2, const double Ms, const double Mb, const double MZ) const {
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

gslpp::complex THDMcache::A_h_L(const double mHl2, const double cW2, const double Me, const double Mmu, const double Mtau, const double MZ) const {
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

gslpp::complex THDMcache::A_HH_L(const double mHh2, const double cW2, const double Mmu, const double Mtau, const double MZ) const {
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

gslpp::complex THDMcache::A_A_L(const double mA2, const double cW2, const double Mmu, const double Mtau, const double MZ) const {
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

gslpp::complex THDMcache::A_H_W(const double mH, const double cW2, const double MW, const double MZ) const {
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

gslpp::complex THDMcache::A_H_Hp(const double mHp2, const double mH, const double cW2, const double MZ) const {
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

gslpp::complex THDMcache::f_func(const double x) const{
    if(x<1) {
    gslpp::complex z = -gslpp::complex::i()*M_PI;
    return -pow(log((1.0+sqrt(1.0-x))/(1.0-sqrt(1.0-x)))+z,2)/4.0;
    }
    else {
        return pow(asin(sqrt(1.0/x)),2);
    }
}



gslpp::complex THDMcache::g_func(const double x) const{
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



gslpp::complex THDMcache::Int1(const double tau, const double lambda) const{
    return tau*lambda/(tau-lambda)/2.0+tau*tau*lambda*lambda/((tau-lambda)
           *(tau-lambda))/2.0*(f_func(tau)-f_func(lambda))+tau*tau*lambda/((tau-lambda)
           *(tau-lambda))*(g_func(tau)-g_func(lambda));
}



gslpp::complex THDMcache::Int2(const double tau, const double lambda) const{
    return -tau*lambda/(tau-lambda)/2.0*(f_func(tau)-f_func(lambda));
}


int THDMcache::HSTheta (const double x) const{
    if(x<0) return 0.0;
    else return 1.0;
}


double THDMcache::KaellenFunction(const double a2, const double b2, const double c2) const{
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


double THDMcache::cW2THDM(const double c02) const{
    return c02;
}



double THDMcache::MWTHDM(const double MW) const{
    return MW;
}



void THDMcache::computeSignalStrengthQuantities()
{
    double sW2=1.0-cW2;
    double sin_ba=sin(bma);
    double sinb=tanb/sqrt(1.0+tanb*tanb);
    double cosb=1.0/sqrt(1.0+tanb*tanb);
    double sina=sinb*cos(bma)-cosb*sin(bma);
    double cosa=cosb*cos(bma)+sinb*sin(bma);

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
    double SigmaggF8 = myTHDM->computeSigmaggH(8.0);
    //The ggH cross section in the SM at 13 TeV.
    double SigmaggF13 = myTHDM->computeSigmaggH(13.0);
    //The square of the top-quark contribution to the ggH cross section in the SM at 8 TeV
    double Sigmaggh_tt8 = myTHDM->computeSigmaggH_tt(8.0);
    //The square of the top-quark contribution to the ggH cross section in the SM at 13 TeV
//    double Sigmaggh_tt13 = myTHDM->computeSigmaggH_tt(13.0);
    //The square of the bottom-quark contribution to the ggH cross section in the SM at 8 TeV
    double Sigmaggh_bb8 = myTHDM->computeSigmaggH_bb(8.0);
    //The square of the bottom-quark contribution to the ggH cross section in the SM at 13 TeV
//    double Sigmaggh_bb13 = myTHDM->computeSigmaggH_bb(13.0);
    //The ttH production cross section in the SM at 8 TeV
    double Sigmatth8 = myTHDM->computeSigmattH(8.0);
    //The ttH production cross section in the SM at 13 TeV
    double Sigmatth13 = myTHDM->computeSigmattH(13.0);
    //The bbH production cross section in the SM at 13 TeV
    double Sigmabbh13 = ip_cs_pptobbH_13(mHl);
    //The VBF plus Vh production cross section in the SM at 13 TeV
    double SigmaVBFVh13 = (myTHDM->computeSigmaVBF(13.0)+myTHDM->computeSigmaWH(13.0)+myTHDM->computeSigmaZH(13.0));

    /* r_ii is the ratio of the squared 2HDM vertex coupling of h to
     * the particle i and the respective squared SM coupling.*/
    rh_QuQu=cosa*cosa/(sinb*sinb);
    rh_VV=sin_ba*sin_ba;
    rh_QdQd=0.0;//It depends on the modelType
    rh_ll=0.0;//It depends on the modelType
    rh_gg=0.0;//It depends on the modelType 

    //Calulation of rh_gg, rh_QdQd, rh_ll, rh_gaga, rh_Zga (depending on the model type): START

    //rh_gaga formula = abs(I_h_F+I_h_W+I_h_Hp)^2 / abs(I_hSM_F+I_hSM_W)^2

    gslpp::complex I_h_F=0.0;//It depends on the modelType
    gslpp::complex fermU=I_h_U(mHl*mHl,Mu,Mc,Mt);
    gslpp::complex fermD=I_h_D(mHl*mHl,Md,Ms,Mb);
    gslpp::complex fermL=I_h_L(mHl*mHl,Me,Mmu,Mtau);
    gslpp::complex I_hSM_W=I_H_W(mHl*mHl,MW);
    gslpp::complex I_h_W=sin_ba*I_hSM_W;
    gslpp::complex I_h_Hp=I_H_Hp(mHp2,mHl*mHl)*ghHpHm(mHp2,tanb,m12_2,bma,mHl*mHl,vev)*vev/(2.0*mHp2);

    double ABSgagaTHDM=0.0;
    double ABSgagaSM=0.0;

    //rh_Zga formula = abs(A_h_F+A_h_W+A_h_Hp)^2 / abs(A_hSM_F+A_hSM_W)^2

    gslpp::complex A_h_F = 0.0;//It depends on the modelType
    gslpp::complex A_h_Ux = A_h_U(mHl*mHl,cW2,Mu,Mc,Mt,MZ);
    gslpp::complex A_h_Dx = A_h_D(mHl*mHl,cW2,Md,Ms,Mb,MZ);
    gslpp::complex A_h_Lx  = A_h_L(mHl*mHl,cW2,Me,Mmu,Mtau,MZ);
    gslpp::complex A_hSM_W = A_H_W(mHl*mHl,cW2,MW,MZ);
    gslpp::complex A_h_W = sin_ba*A_hSM_W;
    gslpp::complex A_h_Hp = A_H_Hp(mHp2,mHl*mHl,cW2,MZ)*ghHpHm(mHp2,tanb,m12_2,bma,mHl*mHl,vev)*vev/(2.0*mHp2);

    double ABSZgaTHDM=0.0;
    double ABSZgaSM=0.0;

    if( modelflag == "type1" ) {
        rh_gg=cosa/sinb*cosa/sinb;
        rh_QdQd=cosa/sinb*cosa/sinb;
        rh_ll=cosa/sinb*cosa/sinb;
        I_h_F=cosa/sinb*(fermU+fermD+fermL);
        A_h_F = cosa/sinb*(A_h_Ux+A_h_Dx+A_h_Lx)/sqrt(sW2*cW2);
    }
    else if( modelflag == "type2" ) {
        rh_gg=-cosa/sinb*sina/cosb+(cosa/sinb+sina/cosb)
             *(Sigmaggh_tt8*cosa/sinb+Sigmaggh_bb8*sina/cosb)/SigmaggF8;
        rh_QdQd=sina/cosb*sina/cosb;
        rh_ll=sina/cosb*sina/cosb;
        I_h_F=cosa/sinb*fermU -sina/cosb*(fermD+fermL);
        A_h_F = (cosa/sinb*A_h_Ux-sina/cosb*(A_h_Dx+A_h_Lx))/sqrt(sW2*cW2);
    }
    else if( modelflag == "typeX" ) {
        rh_gg=cosa/sinb*cosa/sinb;
        rh_QdQd=cosa/sinb*cosa/sinb;
        rh_ll=sina/cosb*sina/cosb;
        I_h_F = cosa/sinb*(fermU+fermD) -sina/cosb*fermL;
        A_h_F = (cosa/sinb*(A_h_Ux+A_h_Dx)-sina/cosb*A_h_Lx)/sqrt(sW2*cW2);
    }
    else if( modelflag == "typeY" ) {
        rh_gg=-cosa/sinb*sina/cosb+(cosa/sinb+sina/cosb)
             *(Sigmaggh_tt8*cosa/sinb+Sigmaggh_bb8*sina/cosb)/SigmaggF8;
        rh_QdQd=sina/cosb*sina/cosb;
        rh_ll=cosa/sinb*cosa/sinb;
        I_h_F = cosa/sinb*(fermU+fermL) -sina/cosb*fermD;
        A_h_F = (cosa/sinb*(A_h_Ux+A_h_Lx)-sina/cosb*A_h_Dx)/sqrt(sW2*cW2);
    }
    else {
        throw std::runtime_error("modelflag can be only any of \"type1\", \"type2\", \"typeX\" or \"typeY\"");
    }

    ABSgagaTHDM=(I_h_F+I_h_W+I_h_Hp).abs2();
    ABSgagaSM=(fermU+fermL+fermD+I_hSM_W).abs2();
    rh_gaga=ABSgagaTHDM/ABSgagaSM;

    ABSZgaTHDM=(A_h_F+A_h_W+A_h_Hp).abs2();
    ABSZgaSM=((A_h_Ux+A_h_Lx+A_h_Dx)/sqrt(sW2*cW2)+A_hSM_W).abs2();
    rh_Zga=ABSZgaTHDM/ABSZgaSM;
    //Calulation of rh_gg, rh_QdQd, rh_ll, rh_gaga, rh_Zga (they depend on the model type): END

    /* ggF_tth8 is the ratio of the THDM and SM cross sections for ggF or tth production at 8 TeV*/
    ggF_tth8 = (SigmaggF8*rh_gg + Sigmatth8*rh_QuQu)/(SigmaggF8 + Sigmatth8);
    /* ggF_tth13 is the ratio of the THDM and SM cross sections for ggF or tth production at 13 TeV */
    ggF_tth13 = (SigmaggF13*rh_gg + Sigmatth13*rh_QuQu)/(SigmaggF13 + Sigmatth13);
    /* pph13 is the ratio of the THDM and SM cross sections for an h production at 13 TeV */
    pph13 = (SigmaggF13*rh_gg + SigmaVBFVh13*rh_VV + Sigmatth13*rh_QuQu + Sigmabbh13*rh_QdQd)/(SigmaggF13 + SigmaVBFVh13 + Sigmatth13 + Sigmabbh13);
    /* VBF_Vh is the ratio of the THDM and SM cross sections for VBF or Vh production */
    VBF_Vh = rh_VV;

    sumModBRs = rh_QdQd*BrSM_htobb + rh_VV*(BrSM_htoWW+BrSM_htoZZ) + rh_ll*BrSM_htotautau +
          rh_gaga*BrSM_htogaga + rh_gg*BrSM_htogg + rh_Zga*BrSM_htoZga + rh_QuQu*BrSM_htocc;

    Gamma_h = sumModBRs*myTHDM->computeGammaHTotal();
    
    THDM_BR_h_bb = rh_QdQd*BrSM_htobb/sumModBRs;
    THDM_BR_h_gaga = rh_gaga*BrSM_htogaga/sumModBRs;
    THDM_BR_h_tautau = rh_ll*BrSM_htotautau/sumModBRs;
    THDM_BR_h_WW = rh_VV*BrSM_htoWW/sumModBRs;
    THDM_BR_h_ZZ = rh_VV*BrSM_htoZZ/sumModBRs;
    THDM_BR_h_gg = rh_gg*BrSM_htogg/sumModBRs;
    THDM_BR_h_cc = rh_QuQu*BrSM_htocc/sumModBRs;
}

void THDMcache::computeHHquantities()
{
    double GF=1/(sqrt(2.0)*vev*vev);
    double sW2=1.0-cW2;
    double mHh=sqrt(mHh2);
    double sin_ba=sin(bma);
    double sinb=tanb/sqrt(1.0+tanb*tanb);
    double cosb=1.0/sqrt(1.0+tanb*tanb);
    double sina=sinb*cos(bma)-cosb*sin(bma);
    double cosa=cosb*cos(bma)+sinb*sin(bma);
    double cos_2b=cosb*cosb-sinb*sinb;
    double cos_ba=cos(bma);

    //These cross sections ratios are necessary for rHH_gg
    //SM gg -> H production cross section ratio at 8 TeV, top loop only over total
    double rSigmaggH_t8 = ip_csr_ggH_t_8(mHh);
    //SM gg -> H production cross section ratio at 8 TeV, bottom loop only over total
    double rSigmaggH_b8 = ip_csr_ggH_b_8(mHh);
    //SM gg -> H production cross section ratio at 13 TeV, top loop only over total
//    double rSigmaggH_t13 = ip_csr_ggH_t_13(mHh);
    //SM gg -> H production cross section ratio at 13 TeV, bottom loop only over total
//    double rSigmaggH_b13 = ip_csr_ggH_b_13(mHh);

    /* r_ii is the ratio between the squared 2HDM vertex coupling of the heavy Higgs to
     * the particle i and the corresponding coupling of the SM Higgs boson.*/
    double rHH_QuQu=sina/sinb*sina/sinb;
    double rHH_QdQd=0.0;//It depends on the modelType
    double rHH_ll=0.0;//It depends on the modelType
    rHH_gg=0.0;//It depends on the modelType
    rHH_VV=cos_ba*cos_ba;

    /*Calulation of rHH_QdQd, rHH_ll, rHH_gg, Gamma_Hgaga, Gamma_HZga, Gamma_Hgg
     * (they depend on the model type): START*/

    /*Gamma_Hgaga and Gamma_HZga expressions can be found in
     "The Higgs Hunter's Guide", Appendix C and in arXiv:0902.4665v3, Appendix A;
     *Gamma_Hgg expression can be found in arXiv:0902.4665v3, Appendix A*/

    /*I_HH_F, I_HH_W and I_HH_Hp are needed for Gamma_Hgaga;
     * their expressions can be found in "The Higgs Hunter's Guide", Appendix C, C.4*/
    gslpp::complex I_HH_F=0.0;//It depends on the modelType
    gslpp::complex I_HH_Ux=I_HH_U(mHh2,Mc,Mt);
    gslpp::complex I_HH_Dx=I_HH_D(mHh2,Ms,Mb);
    gslpp::complex I_HH_Lx=I_HH_L(mHh2,Mmu,Mtau);
    gslpp::complex I_HH_W=cos_ba*I_H_W(mHh,MW);
    /* g_HH_HpHm is the coupling of the heavy Higgs boson to Hp and Hm; its
     * expression can be found in arXiv:1403.1264v2, formula 5*/
    gslpp::complex I_HH_Hp=I_H_Hp(mHp2,mHh)*g_HH_HpHm(mHp2,mHh2,tanb,m12_2,bma,vev)*vev/(2.0*mHp2);

    /*A_HH_F, A_HH_W and A_HH_Hp are needed for Gamma_HZga;
     * their expressions can be found in "The Higgs Hunter's Guide", Appendix C, C.12*/
    gslpp::complex A_HH_F = 0.0;//It depends on the modelType
    gslpp::complex A_HH_Ux = A_HH_U(mHh2,cW2,Mc,Mt,MZ);
    gslpp::complex A_HH_Dx = A_HH_D(mHh2,cW2,Ms,Mb,MZ);
    gslpp::complex A_HH_Lx = A_HH_L(mHh2,cW2,Mmu,Mtau,MZ);
    /*A_HH_W expression can be found in "The Higgs Hunter's Guide", Appendix C, C.13*/
    gslpp::complex A_HH_W = cos_ba*A_H_W(mHh,cW2,MW,MZ);
    /*A_HH_Hp expression can be found in "The Higgs Hunter's Guide", Appendix C, C.14*/
    gslpp::complex A_HH_Hp= A_H_Hp(mHp2,mHh,cW2,MZ)*g_HH_HpHm(mHp2,mHh2,tanb,m12_2,bma,vev)*vev/(2.0*mHp2);

    if( modelflag == "type1" ) {
        rHH_gg=sina/sinb*sina/sinb;
        rHH_QdQd=sina/sinb*sina/sinb;
        rHH_ll=sina/sinb*sina/sinb;
        I_HH_F=sina/sinb*(I_HH_Ux+I_HH_Dx+I_HH_Lx);
        A_HH_F = sina/sinb*(A_HH_Ux+A_HH_Dx+A_HH_Lx)/sqrt(sW2*cW2);
    }
    else if( modelflag == "type2" ) {
        rHH_gg=sina/sinb*cosa/cosb+rSigmaggH_t8*sina/sinb*(sina/sinb-cosa/cosb)
             +rSigmaggH_b8*cosa/cosb*(cosa/cosb-sina/sinb);
        rHH_QdQd=cosa/cosb*cosa/cosb;
        rHH_ll=cosa/cosb*cosa/cosb;
        I_HH_F=sina/sinb*I_HH_Ux+cosa/cosb*(I_HH_Dx+I_HH_Lx);
        A_HH_F = (sina/sinb*A_HH_Ux+cosa/cosb*(A_HH_Dx+A_HH_Lx))/sqrt(sW2*cW2);
    }
    else if( modelflag == "typeX" ) {
        rHH_gg=sina/sinb*sina/sinb;
        rHH_QdQd=sina/sinb*sina/sinb;
        rHH_ll=cosa/cosb*cosa/cosb;
        I_HH_F=sina/sinb*(I_HH_Ux+I_HH_Dx)+cosa/cosb*I_HH_Lx;
        A_HH_F = (sina/sinb*(A_HH_Ux+A_HH_Dx)+cosa/cosb*A_HH_Lx)/sqrt(sW2*cW2);
    }
    else if( modelflag == "typeY" ) {
        rHH_gg=sina/sinb*cosa/cosb+rSigmaggH_t8*sina/sinb*(sina/sinb-cosa/cosb)
             +rSigmaggH_b8*cosa/cosb*(cosa/cosb-sina/sinb);
        rHH_QdQd=cosa/cosb*cosa/cosb;
        rHH_ll=sina/sinb*sina/sinb;
        I_HH_F=sina/sinb*(I_HH_Ux+I_HH_Lx)+cosa/cosb*I_HH_Dx;
        A_HH_F = (sina/sinb*(A_HH_Ux+A_HH_Lx)+cosa/cosb*A_HH_Dx)/sqrt(sW2*cW2);
    }
    else {
        throw std::runtime_error("modelflag can be only any of \"type1\", \"type2\", \"typeX\" or \"typeY\"");
    }

    /*Gamma_Hgaga expression can be found in arXiv:0902.4665v3, Appendix A, A.8*/
    double Gamma_Hgaga=GF*Ale*Ale*mHh*mHh*mHh/(sqrt(2.0)*128.0*M_PI*M_PI*M_PI)
                *(I_HH_F+I_HH_W+I_HH_Hp).abs2();

    /*Gamma_HZga expression can be found in arXiv:0902.4665v3, Appendix A, A.9*/
    double Gamma_HZga=HSTheta(mHh-MZ)*GF*Ale*Ale*mHh*mHh*mHh/(sqrt(2.0)*64.0*M_PI*M_PI*M_PI)
               *(1.0-MZ*MZ/(mHh*mHh))*(1.0-MZ*MZ/(mHh*mHh))*(1.0-MZ*MZ/(mHh*mHh))
               *(A_HH_F+A_HH_W+A_HH_Hp).abs2();

    /*Gamma_Hgg expression can be found in arXiv:0902.4665v3, Appendix A, A.10 or in the Higgs Hunter's Guide (2.30); relative coupling see above*/
    double Gamma_Hgg=rHH_gg*GF*Als*Als*mHh*mHh*mHh/(sqrt(2.0)*16.0*M_PI*M_PI*M_PI)
                     *(9.0/4.0)*(I_HH_Ux/4.0+I_HH_Dx).abs2();

    /*Calulation of rHH_QdQd, rHH_ll, rHH_gg, Gamma_Hgaga, Gamma_HZga, Gamma_Hgg: END*/

    SigmaggF_H8=ip_cs_ggtoH_8(mHh)*rHH_gg;
    SigmabbF_H8=ip_cs_pptobbH_8(mHh)*rHH_QdQd;
    SigmaVBF_H8=ip_cs_VBFtoH_8(mHh)*rHH_VV;
    double SigmattF_H8=ip_cs_pptottH_8(mHh)*rHH_QuQu;
    double SigmaVH_H8=(ip_cs_WtoWH_8(mHh)+ip_cs_ZtoZH_8(mHh))*rHH_VV;
    SigmaTotSM_H8 = 1.0e-15;
    if (mHh>=20. && mHh <=2000.) {
            SigmaTotSM_H8=ip_cs_ggtoH_8(mHh)+ip_cs_VBFtoH_8(mHh)+ip_cs_WtoWH_8(mHh)+ip_cs_ZtoZH_8(mHh)+ip_cs_pptottH_8(mHh)+ip_cs_pptobbH_8(mHh);
    }
    SigmaSumH8 = SigmaggF_H8 + SigmaVBF_H8 + SigmaVH_H8 + SigmattF_H8 + SigmabbF_H8;

    SigmaggF_H13=ip_cs_ggtoH_13(mHh)*rHH_gg;
    SigmabbF_H13=ip_cs_pptobbH_13(mHh)*rHH_QdQd;
    SigmaVBF_H13=ip_cs_VBFtoH_13(mHh)*rHH_VV;
    SigmattF_H13=ip_cs_pptottH_13(mHh)*rHH_QuQu;
    SigmaVH_H13=(ip_cs_WtoWH_13(mHh)+ip_cs_ZtoZH_13(mHh))*rHH_VV;
//    double SigmaTotSM_H13 = 1.0e-15;
//    if (mHh>=20. && mHh <=2000.) {
//            SigmaTotSM_H13=ip_cs_ggtoH_13(mHh)+ip_cs_VBFtoH_13(mHh)+ip_cs_WtoWH_13(mHh)+ip_cs_ZtoZH_13(mHh)+ip_cs_pptottH_13(mHh)+ip_cs_pptobbH_13(mHh);
//    }
    SigmaSumH13 = SigmaggF_H13 + SigmaVBF_H13 + SigmaVH_H13 + SigmattF_H13 + SigmabbF_H13;

    double BrSM_Htott=ip_Br_HPtott(mHh);
    double BrSM_Htocc=ip_Br_HPtocc(mHh);
    double BrSM_Htobb=ip_Br_HPtobb(mHh);
    double BrSM_Htotautau=ip_Br_HPtotautau(mHh);
    double BrSM_Htomumu=ip_Br_HPtomumu(mHh);
    double BrSM_HtoWW =ip_Br_HPtoWW(mHh);
    double BrSM_HtoZZ =ip_Br_HPtoZZ(mHh);

    double GammaHtotSM=ip_GammaHPtotSM(mHh);

    double GammaHhh=HSTheta(mHh - 2.0*sqrt(mHl2))*sqrt(std::fabs(1.0 - (4.0*mHl2)/mHh2))
                    *std::fabs((cos_ba*cos_ba/(4.0*sinb*cosb*sinb*cosb)
                    *pow(m12_2 + mHh2*cosa*sina + (2.0*mHl2 - 3.0*m12_2/(sinb*cosb))
                    *sina*cosa,2))/(vev*vev))/(8.0*mHh*M_PI);

    double GammaHHpHm=HSTheta(mHh - 2.0*sqrt(mHp2))*sqrt(std::fabs(1.0 - (4.0*mHp2)/mHh2))
                      *std::fabs(pow(cos_ba*(mHh2 + 2.0*mHp2 - 2.0*m12_2/sinb/cosb)
                                -cos_2b/(sinb*cosb)*(mHh2 - m12_2/sinb/cosb)*sin_ba,2)/(vev*vev))
                      /(16.0*mHh*M_PI);

    double GammaHAA=HSTheta(mHh-2.0*sqrt(mA2))*sqrt(std::fabs(1.0 - (4.0*mA2)/mHh2))
                    *std::fabs(pow(cos_ba*(2.0*mA2 + mHh2 - 2.0*m12_2/sinb/cosb)
                    - cos_2b/(sinb*cosb)*(mHh2 - m12_2/sinb/cosb)*sin_ba,2)/(vev*vev))
                    /(32.0*mHh*M_PI);

    double GammaHAZ=HSTheta(mHh-sqrt(mA2)-MZ)*pow(KaellenFunction(mHh2,MZ*MZ,mA2),3)
                    *sin_ba*sin_ba/(2.0*M_PI*vev*vev);

    double GammaHHpW=HSTheta(mHh-sqrt(mHp2)-MW)*pow(KaellenFunction(mHh2,MW*MW,mHp2),3)*sin_ba*sin_ba/(M_PI*vev*vev);

    GammaHtot= ((BrSM_Htott+BrSM_Htocc)*rHH_QuQu
                    +BrSM_Htobb*rHH_QdQd
                    +(BrSM_Htotautau+BrSM_Htomumu)*rHH_ll
                    +(BrSM_HtoWW+BrSM_HtoZZ)*rHH_VV)*GammaHtotSM
               +Gamma_Hgaga+Gamma_HZga+Gamma_Hgg
               +GammaHhh+GammaHHpHm+GammaHAA+GammaHAZ+GammaHHpW;
    
    Br_Htott=BrSM_Htott*rHH_QuQu*GammaHtotSM/GammaHtot;
    Br_Htobb=BrSM_Htobb*rHH_QdQd*GammaHtotSM/GammaHtot;
    Br_Htotautau=BrSM_Htotautau*rHH_ll*GammaHtotSM/GammaHtot;
    Br_HtoWW=BrSM_HtoWW*rHH_VV*GammaHtotSM/GammaHtot;
    Br_HtoZZ=BrSM_HtoZZ*rHH_VV*GammaHtotSM/GammaHtot;
    Br_Htogaga=Gamma_Hgaga/GammaHtot;
    Br_HtoZga=Gamma_HZga/GammaHtot;
    Br_Htohh=GammaHhh/GammaHtot;
    Br_HtoAA=GammaHAA/GammaHtot;
    Br_HtoHpHm=GammaHHpHm/GammaHtot;
    Br_HtoAZ=GammaHAZ/GammaHtot;
    Br_HtoHpW=GammaHHpW/GammaHtot;
}

void THDMcache::computeHHlimits()
{
    double mHh=sqrt(mHh2);
    double mA=sqrt(mA2);

    double Br_Ztoee=0.03363; //C. Patrignani et al.(Particle Data Group), Chin. Phys. C, 40, 100001 (2016)
    double Br_Ztomumu=0.03366; //C. Patrignani et al.(Particle Data Group), Chin. Phys. C, 40, 100001 (2016)
    double Br_Ztotautau=0.0337; //C. Patrignani et al.(Particle Data Group), Chin. Phys. C, 40, 100001 (2016)
    double Br_Ztoinv=0.2; //C. Patrignani et al.(Particle Data Group), Chin. Phys. C, 40, 100001 (2016)
    double Br_Wtoenu=0.1071; //C. Patrignani et al.(Particle Data Group), Chin. Phys. C, 40, 100001 (2016)
    double Br_Wtomunu=0.1063; //C. Patrignani et al.(Particle Data Group), Chin. Phys. C, 40, 100001 (2016)
    double Br_Wtotaunu=0.1138; //C. Patrignani et al.(Particle Data Group), Chin. Phys. C, 40, 100001 (2016)

    //Theoretical expressions for the Heavy Higgs cross sections times branching ratios at 8 TeV

    ggF_H_tautau_TH8=SigmaggF_H8*Br_Htotautau;
    bbF_H_tautau_TH8=SigmabbF_H8*Br_Htotautau;
    pp_H_gaga_TH8=SigmaSumH8*Br_Htogaga;
    ggF_H_gaga_TH8=SigmaggF_H8*Br_Htogaga;
    pp_H_Zga_llga_TH8=SigmaSumH8*Br_HtoZga*(Br_Ztoee+Br_Ztomumu);
    mu_pp_H_VV_TH8=SigmaSumH8/SigmaTotSM_H8*rHH_VV*GammaHtotSM/GammaHtot;
    ggF_H_ZZ_TH8=SigmaggF_H8*Br_HtoZZ;
    VBF_H_ZZ_TH8=SigmaVBF_H8*Br_HtoZZ;
    ggF_H_WW_TH8=SigmaggF_H8*Br_HtoWW;
    VBF_H_WW_TH8=SigmaVBF_H8*Br_HtoWW;
    ggF_H_hh_TH8=SigmaggF_H8*Br_Htohh;
    pp_H_hh_TH8=SigmaSumH8*Br_Htohh;
    ggF_H_hh_bbtautau_TH8=SigmaggF_H8*Br_Htohh*THDM_BR_h_bb*THDM_BR_h_tautau*2.0;
    pp_H_hh_bbbb_TH8=SigmaSumH8*Br_Htohh*THDM_BR_h_bb*THDM_BR_h_bb;
    pp_H_hh_gagabb_TH8=SigmaSumH8*Br_Htohh*THDM_BR_h_gaga*THDM_BR_h_bb*2.0;
    ggF_H_tt_TH8=SigmaggF_H8*Br_Htott;
    bbF_H_bb_TH8=SigmabbF_H8*Br_Htobb;
    pp_H_AZ_bbll_TH8=SigmaSumH8*Br_HtoAZ*Br_Atobb*(Br_Ztoee+Br_Ztomumu);
    pp_H_AZ_tautaull_TH8=SigmaSumH8*Br_HtoAZ*Br_Atotautau*(Br_Ztoee+Br_Ztomumu);
    
    //Ratios of theoretical Heavy Higgs cross sections and experimental upper limits at 8 TeV

    THoEX_ggF_H_tautau_ATLAS8=0.0;
    R_ggF_H_tautau_ATLAS8=0.0;
    THoEX_ggF_H_tautau_CMS8=0.0;
    R_ggF_H_tautau_CMS8=0.0;
    THoEX_bbF_H_tautau_ATLAS8=0.0;
    R_bbF_H_tautau_ATLAS8=0.0;
    THoEX_bbF_H_tautau_CMS8=0.0;
    R_bbF_H_tautau_CMS8=0.0;
    THoEX_pp_H_gaga_ATLAS8=0.0;
    R_pp_H_gaga_ATLAS8=0.0;
    THoEX_ggF_H_gaga_CMS8=0.0;
    R_ggF_H_gaga_CMS8=0.0;
//    LIMIT_ggF_H_gaga_CMS8=0.0;
//    LIMEST_ggF_H_gaga_CMS8=0.0;
//    DEVIATION_ggF_H_gaga_CMS8=0.0;
//    BANDSIZE_ggF_H_gaga_CMS8=0.0;
    THoEX_pp_H_Zga_llga_CMS8=0.0;
    R_pp_H_Zga_llga_CMS8=0.0;
    THoEX_pp_H_Zga_llga_ATLAS8=0.0;
    R_pp_H_Zga_llga_ATLAS8=0.0;
    THoEX_mu_pp_H_VV_CMS8=0.0;
    R_mu_pp_H_VV_CMS8=0.0;
    THoEX_ggF_H_WW_ATLAS8=0.0;
    R_ggF_H_WW_ATLAS8=0.0;
    THoEX_VBF_H_WW_ATLAS8=0.0;
    R_VBF_H_WW_ATLAS8=0.0;
    THoEX_ggF_H_ZZ_ATLAS8=0.0;
    R_ggF_H_ZZ_ATLAS8=0.0;
    THoEX_VBF_H_ZZ_ATLAS8=0.0;
    R_VBF_H_ZZ_ATLAS8=0.0;
    THoEX_ggF_H_hh_ATLAS8=0.0;
    R_ggF_H_hh_ATLAS8=0.0;
    THoEX_pp_H_hh_CMS8=0.0;
    R_pp_H_hh_CMS8=0.0;
    THoEX_ggF_H_hh_bbtautau_CMS8=0.0;
    R_ggF_H_hh_bbtautau_CMS8=0.0;
    THoEX_pp_H_hh_bbbb_CMS8=0.0;
    R_pp_H_hh_bbbb_CMS8=0.0;
    THoEX_pp_H_hh_gagabb_CMS8=0.0;
    R_pp_H_hh_gagabb_CMS8=0.0;
    THoEX_ggF_H_tt_ATLAS8=0.0;
    R_ggF_H_tt_ATLAS8=0.0;
    THoEX_bbF_H_bb_CMS8=0.0;
    R_bbF_H_bb_CMS8=0.0;
    THoEX_pp_H_AZ_bbll_CMS8=0.0;
    R_pp_H_AZ_bbll_CMS8=0.0;
    THoEX_pp_H_AZ_tautaull_CMS8=0.0;
    R_pp_H_AZ_tautaull_CMS8=0.0;

    //Theoretical expressions for the Heavy Higgs cross sections times branching ratios at 13 TeV

    ggF_H_tautau_TH13=SigmaggF_H13*Br_Htotautau;
    bbF_H_tautau_TH13=SigmabbF_H13*Br_Htotautau;
    pp_H_gaga_TH13=SigmaSumH13*Br_Htogaga;
    ggF_H_gaga_TH13=SigmaggF_H13*Br_Htogaga;
    pp_H_Zga_TH13=SigmaSumH13*Br_HtoZga;
    ggF_H_Zga_TH13=SigmaggF_H13*Br_HtoZga;
    pp_H_ZZ_TH13=SigmaSumH13*Br_HtoZZ;
    ggF_H_ZZ_TH13=SigmaggF_H13*Br_HtoZZ;
    VBF_H_ZZ_TH13=SigmaVBF_H13*Br_HtoZZ;
    ggF_H_ZZ_llll_TH13=SigmaggF_H13*Br_HtoZZ*(Br_Ztoee+Br_Ztomumu)*(Br_Ztoee+Br_Ztomumu);
    VBF_H_ZZ_llll_TH13=SigmaVBF_H13*Br_HtoZZ*(Br_Ztoee+Br_Ztomumu)*(Br_Ztoee+Br_Ztomumu);
    pp_H_ZZ_llll_TH13=SigmaSumH13*Br_HtoZZ*(Br_Ztoee+Br_Ztomumu)*(Br_Ztoee+Br_Ztomumu);
    VBF_VH_H_ZZ_llll_TH13=(SigmaVBF_H13+SigmaVH_H13)*Br_HtoZZ*(Br_Ztoee+Br_Ztomumu)*(Br_Ztoee+Br_Ztomumu);
    ggF_H_WW_TH13=SigmaggF_H13*Br_HtoWW;
    VBF_H_WW_TH13=SigmaVBF_H13*Br_HtoWW;
    ggF_VBF_H_WW_lnulnu_TH13=(SigmaggF_H13+SigmaVBF_H13)*Br_HtoWW*(Br_Wtoenu+Br_Wtomunu)*(Br_Wtoenu+Br_Wtomunu);
    pp_H_VV_TH13=SigmaSumH13*(Br_HtoZZ+Br_HtoWW);
    ggF_H_hh_TH13=SigmaggF_H13*Br_Htohh;
    pp_H_hh_TH13=SigmaSumH13*Br_Htohh;
    pp_H_hh_bbbb_TH13=SigmaSumH13*Br_Htohh*THDM_BR_h_bb*THDM_BR_h_bb;
    ggF_H_hh_bbbb_TH13=SigmaggF_H13*Br_Htohh*THDM_BR_h_bb*THDM_BR_h_bb;
    pp_H_hh_gagabb_TH13=SigmaSumH13*Br_Htohh*THDM_BR_h_gaga*THDM_BR_h_bb*2.0;
    pp_H_hh_bbtautau_TH13=SigmaSumH13*Br_Htohh*THDM_BR_h_bb*THDM_BR_h_tautau*2.0;
    pp_H_hh_bblnulnu_TH13=SigmaSumH13*Br_Htohh*5.77e-1*2.15e-1*(Br_Wtoenu+Br_Wtomunu)*(Br_Wtoenu+Br_Wtomunu)*2.0;/*SM BR assumed in the CMS analysis!*/
    pp_H_hh_bbVV_TH13=SigmaSumH13*Br_Htohh*2.0*THDM_BR_h_bb*(THDM_BR_h_WW*pow(Br_Wtoenu+Br_Wtomunu+Br_Wtotaunu*0.352,2)
                                                            +THDM_BR_h_ZZ*2.0*Br_Ztoinv*(Br_Ztoee+Br_Ztomumu+Br_Ztotautau*0.124));
    ttF_H_tt_TH13=SigmattF_H13*Br_Htott;
    bbF_H_tt_TH13=SigmabbF_H13*Br_Htott;
    pp_H_bb_TH13=SigmaSumH13*Br_Htobb;

    //Ratios of theoretical Heavy Higgs cross sections and experimental upper limits at 13 TeV

    THoEX_ttF_H_tt_ATLAS13=0.0;
    R_ttF_H_tt_ATLAS13=0.0;
    THoEX_bbF_H_tt_ATLAS13=0.0;
    R_bbF_H_tt_ATLAS13=0.0;
    THoEX_ggF_H_tautau_ATLAS13=0.0;
    R_ggF_H_tautau_ATLAS13=0.0;
    THoEX_bbF_H_tautau_ATLAS13=0.0;
    R_bbF_H_tautau_ATLAS13=0.0;
    THoEX_ggF_H_tautau_CMS13=0.0;
    R_ggF_H_tautau_CMS13=0.0;
    THoEX_bbF_H_tautau_CMS13=0.0;
    R_bbF_H_tautau_CMS13=0.0;
    THoEX_pp_H_gaga_ATLAS13=0.0;
    R_pp_H_gaga_ATLAS13=0.0;
    THoEX_ggF_H_gaga_CMS13=0.0;
    R_ggF_H_gaga_CMS13=0.0;
    THoEX_pp_H_Zga_ATLAS13=0.0;
    R_pp_H_Zga_ATLAS13=0.0;
    THoEX_ggF_H_Zga_llga_ATLAS13=0.0;
    R_ggF_H_Zga_llga_ATLAS13=0.0;
    THoEX_pp_H_Zga_llga_CMS13=0.0;
    R_pp_H_Zga_llga_CMS13=0.0;
    THoEX_pp_H_Zga_qqga_CMS13=0.0;
    R_pp_H_Zga_qqga_CMS13=0.0;
    THoEX_ggF_H_Zga_CMS13=0.0;
    R_ggF_H_Zga_CMS13=0.0;
    THoEX_ggF_H_ZZ_llllnunu_ATLAS13=0.0;
    R_ggF_H_ZZ_llllnunu_ATLAS13=0.0;
    THoEX_VBF_H_ZZ_llllnunu_ATLAS13=0.0;
    R_VBF_H_ZZ_llllnunu_ATLAS13=0.0;
    THoEX_ggF_H_ZZ_llnunu_ATLAS13=0.0;
    R_ggF_H_ZZ_llnunu_ATLAS13=0.0;
    THoEX_pp_H_ZZ_llnunu_CMS13=0.0;
    R_pp_H_ZZ_llnunu_CMS13=0.0;
    THoEX_ggF_H_ZZ_llnunu_CMS13=0.0;
    R_ggF_H_ZZ_llnunu_CMS13=0.0;
    THoEX_VBF_H_ZZ_llnunu_CMS13=0.0;
    R_VBF_H_ZZ_llnunu_CMS13=0.0;
    THoEX_ggF_H_ZZ_llll_ATLAS13=0.0;
    R_ggF_H_ZZ_llll_ATLAS13=0.0;
    THoEX_VBF_H_ZZ_llll_ATLAS13=0.0;
    R_VBF_H_ZZ_llll_ATLAS13=0.0;
    THoEX_pp_H_ZZ_llll_CMS13=0.0;
    R_pp_H_ZZ_llll_CMS13=0.0;
    THoEX_VBF_VH_H_ZZ_llll_CMS13=0.0;
    R_VBF_VH_H_ZZ_llll_CMS13=0.0;
    THoEX_ggF_H_ZZ_qqllnunu_ATLAS13=0.0;
    R_ggF_H_ZZ_qqllnunu_ATLAS13=0.0;
    THoEX_VBF_H_ZZ_qqllnunu_ATLAS13=0.0;
    R_VBF_H_ZZ_qqllnunu_ATLAS13=0.0;
    THoEX_ggF_H_ZZ_llqq_ATLAS13=0.0;
    R_ggF_H_ZZ_llqq_ATLAS13=0.0;
    THoEX_VBF_H_ZZ_llqq_ATLAS13=0.0;
    R_VBF_H_ZZ_llqq_ATLAS13=0.0;
    THoEX_ggF_H_ZZ_nunuqq_ATLAS13=0.0;
    R_ggF_H_ZZ_nunuqq_ATLAS13=0.0;
    THoEX_pp_H_ZZ_llqq_CMS13=0.0;
    R_pp_H_ZZ_llqq_CMS13=0.0;
    THoEX_ggF_H_WW_lnuqq_ATLAS13=0.0;
    R_ggF_H_WW_lnuqq_ATLAS13=0.0;
    THoEX_VBF_H_WW_lnuqq_ATLAS13=0.0;
    R_VBF_H_WW_lnuqq_ATLAS13=0.0;
    THoEX_ggF_H_WW_enumunu_ATLAS13=0.0;
    R_ggF_H_WW_enumunu_ATLAS13=0.0;
    THoEX_VBF_H_WW_enumunu_ATLAS13=0.0;
    R_VBF_H_WW_enumunu_ATLAS13=0.0;
    THoEX_ggF_VBF_H_WW_lnulnu_CMS13=0.0;
    R_ggF_VBF_H_WW_lnulnu_CMS13=0.0;
    THoEX_pp_H_VV_qqqq_ATLAS13=0.0;
    R_pp_H_VV_qqqq_ATLAS13=0.0;
    THoEX_pp_H_hh_bbgaga_ATLAS13=0.0;
    R_pp_H_hh_bbgaga_ATLAS13=0.0;
    THoEX_pp_H_hh_bbgaga_CMS13=0.0;
    R_pp_H_hh_bbgaga_CMS13=0.0;
    THoEX_pp_H_hh_bbbb_ATLAS13=0.0;
    R_pp_H_hh_bbbb_ATLAS13=0.0;
    THoEX_pp_H_hh_bbbb_CMS13=0.0;
    R_pp_H_hh_bbbb_CMS13=0.0;
    THoEX_ggF_H_hh_bbbb_CMS13=0.0;
    R_ggF_H_hh_bbbb_CMS13=0.0;
    THoEX_ggF_H_hh_gagaWW_ATLAS13=0.0;
    R_ggF_H_hh_gagaWW_ATLAS13=0.0;
    THoEX_pp_H_hh_bbtautau_CMS13=0.0;
    R_pp_H_hh_bbtautau_CMS13=0.0;
    THoEX_pp_H_hh_bbtautau1_CMS13=0.0;
    R_pp_H_hh_bbtautau1_CMS13=0.0;
    THoEX_pp_H_hh_bblnulnu_CMS13=0.0;
    R_pp_H_hh_bblnulnu_CMS13=0.0;
    THoEX_pp_H_hh_bbVV_CMS13=0.0;
    R_pp_H_hh_bbVV_CMS13=0.0;
    THoEX_pp_H_bb_CMS13=0.0;
    R_pp_H_bb_CMS13=0.0;

    //95% to 1 sigma conversion factor, roughly sqrt(3.84)
    double nftos=1.95996398454;

    if(mHh>=65.0 && mHh<90.0)
    {
        THoEX_pp_H_gaga_ATLAS8=pp_H_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mHh);
        R_pp_H_gaga_ATLAS8=(1+(pp_H_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mHh))/ip_ex_pp_phi_gaga_ATLAS8_e(mHh) ) * nftos;
        if(mA>=50.0 && mA<=1000.0)
        {
            THoEX_pp_H_AZ_tautaull_CMS8=pp_H_AZ_tautaull_TH8/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh);
            R_pp_H_AZ_tautaull_CMS8=(1+(pp_H_AZ_tautaull_TH8-ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }
    }
    else if(mHh>=90.0 && mHh<100.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mHh);
        R_ggF_H_tautau_CMS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mHh))/ip_ex_gg_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mHh);
        R_bbF_H_tautau_CMS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mHh))/ip_ex_bb_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS8=pp_H_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mHh);
        R_pp_H_gaga_ATLAS8=(1+(pp_H_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mHh))/ip_ex_pp_phi_gaga_ATLAS8_e(mHh) ) * nftos;
        if(mA>=50.0 && mA<=1000.0)
        {
            THoEX_pp_H_AZ_tautaull_CMS8=pp_H_AZ_tautaull_TH8/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh);
            R_pp_H_AZ_tautaull_CMS8=(1+(pp_H_AZ_tautaull_TH8-ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=100.0 && mHh<130.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mHh);
        R_ggF_H_tautau_CMS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mHh))/ip_ex_gg_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mHh);
        R_bbF_H_tautau_CMS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mHh))/ip_ex_bb_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS8=pp_H_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mHh);
        R_pp_H_gaga_ATLAS8=(1+(pp_H_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mHh))/ip_ex_pp_phi_gaga_ATLAS8_e(mHh) ) * nftos;
        if(mA>=50.0 && mA<=1000.0)
        {
            THoEX_pp_H_AZ_tautaull_CMS8=pp_H_AZ_tautaull_TH8/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh);
            R_pp_H_AZ_tautaull_CMS8=(1+(pp_H_AZ_tautaull_TH8-ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=130.0 && mHh<140.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mHh);
        R_ggF_H_tautau_CMS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mHh))/ip_ex_gg_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mHh);
        R_bbF_H_tautau_CMS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mHh))/ip_ex_bb_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS8=pp_H_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mHh);
        R_pp_H_gaga_ATLAS8=(1+(pp_H_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mHh))/ip_ex_pp_phi_gaga_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_bb_CMS8=bbF_H_bb_TH8/ip_ex_bb_phi_bb_CMS8(mHh);
        R_bbF_H_bb_CMS8=(1+(bbF_H_bb_TH8-ip_ex_bb_phi_bb_CMS8(mHh))/ip_ex_bb_phi_bb_CMS8_e(mHh) ) * nftos;
        if(mA>=50.0 && mA<=1000.0)
        {
            THoEX_pp_H_AZ_tautaull_CMS8=pp_H_AZ_tautaull_TH8/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh);
            R_pp_H_AZ_tautaull_CMS8=(1+(pp_H_AZ_tautaull_TH8-ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=140.0 && mHh<145.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mHh);
        R_ggF_H_tautau_CMS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mHh))/ip_ex_gg_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mHh);
        R_bbF_H_tautau_CMS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mHh))/ip_ex_bb_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS8=pp_H_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mHh);
        R_pp_H_gaga_ATLAS8=(1+(pp_H_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mHh))/ip_ex_pp_phi_gaga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_ATLAS8=ggF_H_ZZ_TH8/ip_ex_gg_H_ZZ_ATLAS8(mHh);
        R_ggF_H_ZZ_ATLAS8=(1+(ggF_H_ZZ_TH8-ip_ex_gg_H_ZZ_ATLAS8(mHh))/ip_ex_gg_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_ATLAS8=VBF_H_ZZ_TH8/ip_ex_VBF_H_ZZ_ATLAS8(mHh);
        R_VBF_H_ZZ_ATLAS8=(1+(VBF_H_ZZ_TH8-ip_ex_VBF_H_ZZ_ATLAS8(mHh))/ip_ex_VBF_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_bb_CMS8=bbF_H_bb_TH8/ip_ex_bb_phi_bb_CMS8(mHh);
        R_bbF_H_bb_CMS8=(1+(bbF_H_bb_TH8-ip_ex_bb_phi_bb_CMS8(mHh))/ip_ex_bb_phi_bb_CMS8_e(mHh) ) * nftos;
        if(mA>=50.0 && mA<=1000.0)
        {
            THoEX_pp_H_AZ_tautaull_CMS8=pp_H_AZ_tautaull_TH8/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh);
            R_pp_H_AZ_tautaull_CMS8=(1+(pp_H_AZ_tautaull_TH8-ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=145.0 && mHh<150.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mHh);
        R_ggF_H_tautau_CMS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mHh))/ip_ex_gg_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mHh);
        R_bbF_H_tautau_CMS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mHh))/ip_ex_bb_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS8=pp_H_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mHh);
        R_pp_H_gaga_ATLAS8=(1+(pp_H_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mHh))/ip_ex_pp_phi_gaga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_ATLAS8=ggF_H_ZZ_TH8/ip_ex_gg_H_ZZ_ATLAS8(mHh);
        R_ggF_H_ZZ_ATLAS8=(1+(ggF_H_ZZ_TH8-ip_ex_gg_H_ZZ_ATLAS8(mHh))/ip_ex_gg_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_ATLAS8=VBF_H_ZZ_TH8/ip_ex_VBF_H_ZZ_ATLAS8(mHh);
        R_VBF_H_ZZ_ATLAS8=(1+(VBF_H_ZZ_TH8-ip_ex_VBF_H_ZZ_ATLAS8(mHh))/ip_ex_VBF_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_mu_pp_H_VV_CMS8=mu_pp_H_VV_TH8/ip_ex_mu_pp_H_VV_CMS8(mHh);
        R_mu_pp_H_VV_CMS8=(1+(mu_pp_H_VV_TH8-ip_ex_mu_pp_H_VV_CMS8(mHh))/ip_ex_mu_pp_H_VV_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_bb_CMS8=bbF_H_bb_TH8/ip_ex_bb_phi_bb_CMS8(mHh);
        R_bbF_H_bb_CMS8=(1+(bbF_H_bb_TH8-ip_ex_bb_phi_bb_CMS8(mHh))/ip_ex_bb_phi_bb_CMS8_e(mHh) ) * nftos;
        if(mA>=50.0 && mA<=1000.0)
        {
            THoEX_pp_H_AZ_tautaull_CMS8=pp_H_AZ_tautaull_TH8/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh);
            R_pp_H_AZ_tautaull_CMS8=(1+(pp_H_AZ_tautaull_TH8-ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=150.0 && mHh<175.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mHh);
        R_ggF_H_tautau_CMS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mHh))/ip_ex_gg_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mHh);
        R_bbF_H_tautau_CMS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mHh))/ip_ex_bb_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS8=pp_H_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mHh);
        R_pp_H_gaga_ATLAS8=(1+(pp_H_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mHh))/ip_ex_pp_phi_gaga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_ATLAS8=ggF_H_ZZ_TH8/ip_ex_gg_H_ZZ_ATLAS8(mHh);
        R_ggF_H_ZZ_ATLAS8=(1+(ggF_H_ZZ_TH8-ip_ex_gg_H_ZZ_ATLAS8(mHh))/ip_ex_gg_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_ATLAS8=VBF_H_ZZ_TH8/ip_ex_VBF_H_ZZ_ATLAS8(mHh);
        R_VBF_H_ZZ_ATLAS8=(1+(VBF_H_ZZ_TH8-ip_ex_VBF_H_ZZ_ATLAS8(mHh))/ip_ex_VBF_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_mu_pp_H_VV_CMS8=mu_pp_H_VV_TH8/ip_ex_mu_pp_H_VV_CMS8(mHh);
        R_mu_pp_H_VV_CMS8=(1+(mu_pp_H_VV_TH8-ip_ex_mu_pp_H_VV_CMS8(mHh))/ip_ex_mu_pp_H_VV_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_bb_CMS8=bbF_H_bb_TH8/ip_ex_bb_phi_bb_CMS8(mHh);
        R_bbF_H_bb_CMS8=(1+(bbF_H_bb_TH8-ip_ex_bb_phi_bb_CMS8(mHh))/ip_ex_bb_phi_bb_CMS8_e(mHh) ) * nftos;
        if(mA>=50.0 && mA<=1000.0)
        {
            THoEX_pp_H_AZ_tautaull_CMS8=pp_H_AZ_tautaull_TH8/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh);
            R_pp_H_AZ_tautaull_CMS8=(1+(pp_H_AZ_tautaull_TH8-ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=175.0 && mHh<200.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mHh);
        R_ggF_H_tautau_CMS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mHh))/ip_ex_gg_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mHh);
        R_bbF_H_tautau_CMS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mHh))/ip_ex_bb_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS8=pp_H_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mHh);
        R_pp_H_gaga_ATLAS8=(1+(pp_H_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mHh))/ip_ex_pp_phi_gaga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_ATLAS8=ggF_H_ZZ_TH8/ip_ex_gg_H_ZZ_ATLAS8(mHh);
        R_ggF_H_ZZ_ATLAS8=(1+(ggF_H_ZZ_TH8-ip_ex_gg_H_ZZ_ATLAS8(mHh))/ip_ex_gg_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_ATLAS8=VBF_H_ZZ_TH8/ip_ex_VBF_H_ZZ_ATLAS8(mHh);
        R_VBF_H_ZZ_ATLAS8=(1+(VBF_H_ZZ_TH8-ip_ex_VBF_H_ZZ_ATLAS8(mHh))/ip_ex_VBF_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_mu_pp_H_VV_CMS8=mu_pp_H_VV_TH8/ip_ex_mu_pp_H_VV_CMS8(mHh);
        R_mu_pp_H_VV_CMS8=(1+(mu_pp_H_VV_TH8-ip_ex_mu_pp_H_VV_CMS8(mHh))/ip_ex_mu_pp_H_VV_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_bb_CMS8=bbF_H_bb_TH8/ip_ex_bb_phi_bb_CMS8(mHh);
        R_bbF_H_bb_CMS8=(1+(bbF_H_bb_TH8-ip_ex_bb_phi_bb_CMS8(mHh))/ip_ex_bb_phi_bb_CMS8_e(mHh) ) * nftos;
        if(mA>=40.0 && mA<=910.0)
        {
            THoEX_pp_H_AZ_bbll_CMS8=pp_H_AZ_bbll_TH8/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh);
            R_pp_H_AZ_bbll_CMS8=(1+(pp_H_AZ_bbll_TH8-ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh))/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mA>=50.0 && mA<=1000.0)
        {
            THoEX_pp_H_AZ_tautaull_CMS8=pp_H_AZ_tautaull_TH8/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh);
            R_pp_H_AZ_tautaull_CMS8=(1+(pp_H_AZ_tautaull_TH8-ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=200.0 && mHh<220.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mHh);
        R_ggF_H_tautau_CMS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mHh))/ip_ex_gg_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mHh);
        R_bbF_H_tautau_CMS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mHh))/ip_ex_bb_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS8=pp_H_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mHh);
        R_pp_H_gaga_ATLAS8=(1+(pp_H_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mHh))/ip_ex_pp_phi_gaga_ATLAS8_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_CMS8=pp_H_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mHh);
        R_pp_H_Zga_llga_CMS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mHh))/ip_ex_pp_A_Zga_llga_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_ATLAS8=pp_H_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mHh);
        R_pp_H_Zga_llga_ATLAS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_ATLAS8=ggF_H_ZZ_TH8/ip_ex_gg_H_ZZ_ATLAS8(mHh);
        R_ggF_H_ZZ_ATLAS8=(1+(ggF_H_ZZ_TH8-ip_ex_gg_H_ZZ_ATLAS8(mHh))/ip_ex_gg_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_ATLAS8=VBF_H_ZZ_TH8/ip_ex_VBF_H_ZZ_ATLAS8(mHh);
        R_VBF_H_ZZ_ATLAS8=(1+(VBF_H_ZZ_TH8-ip_ex_VBF_H_ZZ_ATLAS8(mHh))/ip_ex_VBF_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_mu_pp_H_VV_CMS8=mu_pp_H_VV_TH8/ip_ex_mu_pp_H_VV_CMS8(mHh);
        R_mu_pp_H_VV_CMS8=(1+(mu_pp_H_VV_TH8-ip_ex_mu_pp_H_VV_CMS8(mHh))/ip_ex_mu_pp_H_VV_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_bb_CMS8=bbF_H_bb_TH8/ip_ex_bb_phi_bb_CMS8(mHh);
        R_bbF_H_bb_CMS8=(1+(bbF_H_bb_TH8-ip_ex_bb_phi_bb_CMS8(mHh))/ip_ex_bb_phi_bb_CMS8_e(mHh) ) * nftos;
        if(mA>=40.0 && mA<=910.0)
        {
            THoEX_pp_H_AZ_bbll_CMS8=pp_H_AZ_bbll_TH8/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh);
            R_pp_H_AZ_bbll_CMS8=(1+(pp_H_AZ_bbll_TH8-ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh))/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mA>=50.0 && mA<=1000.0)
        {
            THoEX_pp_H_AZ_tautaull_CMS8=pp_H_AZ_tautaull_TH8/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh);
            R_pp_H_AZ_tautaull_CMS8=(1+(pp_H_AZ_tautaull_TH8-ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ggF_H_tautau_ATLAS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mHh);
        R_ggF_H_tautau_ATLAS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mHh))/ip_ex_gg_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mHh);
        R_bbF_H_tautau_ATLAS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mHh))/ip_ex_bb_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS13=pp_H_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mHh);
        R_pp_H_gaga_ATLAS13=(1+(pp_H_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mHh))/ip_ex_pp_phi_gaga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_llllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llnunu_CMS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llnunu_CMS13(mHh);
        R_ggF_H_ZZ_llnunu_CMS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llnunu_CMS13(mHh))/ip_ex_gg_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llnunu_CMS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llnunu_CMS13(mHh);
        R_VBF_H_ZZ_llnunu_CMS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llnunu_CMS13(mHh))/ip_ex_VBF_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llll_ATLAS13=ggF_H_ZZ_llll_TH13/ip_ex_gg_H_ZZ_llll_ATLAS13(mHh);
        R_ggF_H_ZZ_llll_ATLAS13=(1+(ggF_H_ZZ_llll_TH13-ip_ex_gg_H_ZZ_llll_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llll_ATLAS13=VBF_H_ZZ_llll_TH13/ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh);
        R_VBF_H_ZZ_llll_ATLAS13=(1+(VBF_H_ZZ_llll_TH13-ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_VBF_H_WW_lnulnu_CMS13=ggF_VBF_H_WW_lnulnu_TH13/ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh);
        R_ggF_VBF_H_WW_lnulnu_CMS13=(1+(ggF_VBF_H_WW_lnulnu_TH13-ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh))/ip_ex_ggVV_H_WW_lnulnu_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=220.0 && mHh<250.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mHh);
        R_ggF_H_tautau_CMS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mHh))/ip_ex_gg_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mHh);
        R_bbF_H_tautau_CMS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mHh))/ip_ex_bb_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS8=pp_H_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mHh);
        R_pp_H_gaga_ATLAS8=(1+(pp_H_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mHh))/ip_ex_pp_phi_gaga_ATLAS8_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_CMS8=pp_H_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mHh);
        R_pp_H_Zga_llga_CMS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mHh))/ip_ex_pp_A_Zga_llga_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_ATLAS8=pp_H_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mHh);
        R_pp_H_Zga_llga_ATLAS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_ATLAS8=ggF_H_ZZ_TH8/ip_ex_gg_H_ZZ_ATLAS8(mHh);
        R_ggF_H_ZZ_ATLAS8=(1+(ggF_H_ZZ_TH8-ip_ex_gg_H_ZZ_ATLAS8(mHh))/ip_ex_gg_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_ATLAS8=VBF_H_ZZ_TH8/ip_ex_VBF_H_ZZ_ATLAS8(mHh);
        R_VBF_H_ZZ_ATLAS8=(1+(VBF_H_ZZ_TH8-ip_ex_VBF_H_ZZ_ATLAS8(mHh))/ip_ex_VBF_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_mu_pp_H_VV_CMS8=mu_pp_H_VV_TH8/ip_ex_mu_pp_H_VV_CMS8(mHh);
        R_mu_pp_H_VV_CMS8=(1+(mu_pp_H_VV_TH8-ip_ex_mu_pp_H_VV_CMS8(mHh))/ip_ex_mu_pp_H_VV_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_bb_CMS8=bbF_H_bb_TH8/ip_ex_bb_phi_bb_CMS8(mHh);
        R_bbF_H_bb_CMS8=(1+(bbF_H_bb_TH8-ip_ex_bb_phi_bb_CMS8(mHh))/ip_ex_bb_phi_bb_CMS8_e(mHh) ) * nftos;
        if(mA>=40.0 && mA<=910.0)
        {
            THoEX_pp_H_AZ_bbll_CMS8=pp_H_AZ_bbll_TH8/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh);
            R_pp_H_AZ_bbll_CMS8=(1+(pp_H_AZ_bbll_TH8-ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh))/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mA>=50.0 && mA<=1000.0)
        {
            THoEX_pp_H_AZ_tautaull_CMS8=pp_H_AZ_tautaull_TH8/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh);
            R_pp_H_AZ_tautaull_CMS8=(1+(pp_H_AZ_tautaull_TH8-ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ggF_H_tautau_ATLAS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mHh);
        R_ggF_H_tautau_ATLAS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mHh))/ip_ex_gg_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mHh);
        R_bbF_H_tautau_ATLAS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mHh))/ip_ex_bb_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS13=pp_H_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mHh);
        R_pp_H_gaga_ATLAS13=(1+(pp_H_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mHh))/ip_ex_pp_phi_gaga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_llllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llnunu_CMS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llnunu_CMS13(mHh);
        R_ggF_H_ZZ_llnunu_CMS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llnunu_CMS13(mHh))/ip_ex_gg_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llnunu_CMS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llnunu_CMS13(mHh);
        R_VBF_H_ZZ_llnunu_CMS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llnunu_CMS13(mHh))/ip_ex_VBF_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llll_ATLAS13=ggF_H_ZZ_llll_TH13/ip_ex_gg_H_ZZ_llll_ATLAS13(mHh);
        R_ggF_H_ZZ_llll_ATLAS13=(1+(ggF_H_ZZ_llll_TH13-ip_ex_gg_H_ZZ_llll_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llll_ATLAS13=VBF_H_ZZ_llll_TH13/ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh);
        R_VBF_H_ZZ_llll_ATLAS13=(1+(VBF_H_ZZ_llll_TH13-ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_VBF_H_WW_lnulnu_CMS13=ggF_VBF_H_WW_lnulnu_TH13/ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh);
        R_ggF_VBF_H_WW_lnulnu_CMS13=(1+(ggF_VBF_H_WW_lnulnu_TH13-ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh))/ip_ex_ggVV_H_WW_lnulnu_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=250.0 && mHh<260.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mHh);
        R_ggF_H_tautau_CMS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mHh))/ip_ex_gg_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mHh);
        R_bbF_H_tautau_CMS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mHh))/ip_ex_bb_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS8=pp_H_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mHh);
        R_pp_H_gaga_ATLAS8=(1+(pp_H_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mHh))/ip_ex_pp_phi_gaga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS8=ggF_H_gaga_TH8/ip_ex_gg_phi_gaga_CMS8(mHh);
        R_ggF_H_gaga_CMS8=(1+(ggF_H_gaga_TH8-ip_ex_gg_phi_gaga_CMS8(mHh))/ip_ex_gg_phi_gaga_CMS8_e(mHh) ) * nftos;
//    LIMIT_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8(mHh);
//    LIMEST_ggF_H_gaga_CMS8=(ggF_H_gaga_TH8+ip_ex_gg_phi_gaga_CMS8_e(mHh)-ip_ex_gg_phi_gaga_CMS8(mHh))/ip_ex_gg_phi_gaga_CMS8_e(mHh) ;
//    DEVIATION_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8(mHh)-ip_ex_gg_phi_gaga_CMS8_e(mHh);
//    BANDSIZE_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8_ep2(mHh)-ip_ex_gg_phi_gaga_CMS8_em2(mHh);
        THoEX_pp_H_Zga_llga_CMS8=pp_H_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mHh);
        R_pp_H_Zga_llga_CMS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mHh))/ip_ex_pp_A_Zga_llga_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_ATLAS8=pp_H_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mHh);
        R_pp_H_Zga_llga_ATLAS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_ATLAS8=ggF_H_ZZ_TH8/ip_ex_gg_H_ZZ_ATLAS8(mHh);
        R_ggF_H_ZZ_ATLAS8=(1+(ggF_H_ZZ_TH8-ip_ex_gg_H_ZZ_ATLAS8(mHh))/ip_ex_gg_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_ATLAS8=VBF_H_ZZ_TH8/ip_ex_VBF_H_ZZ_ATLAS8(mHh);
        R_VBF_H_ZZ_ATLAS8=(1+(VBF_H_ZZ_TH8-ip_ex_VBF_H_ZZ_ATLAS8(mHh))/ip_ex_VBF_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_mu_pp_H_VV_CMS8=mu_pp_H_VV_TH8/ip_ex_mu_pp_H_VV_CMS8(mHh);
        R_mu_pp_H_VV_CMS8=(1+(mu_pp_H_VV_TH8-ip_ex_mu_pp_H_VV_CMS8(mHh))/ip_ex_mu_pp_H_VV_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_bb_CMS8=bbF_H_bb_TH8/ip_ex_bb_phi_bb_CMS8(mHh);
        R_bbF_H_bb_CMS8=(1+(bbF_H_bb_TH8-ip_ex_bb_phi_bb_CMS8(mHh))/ip_ex_bb_phi_bb_CMS8_e(mHh) ) * nftos;
        if(mA>=40.0 && mA<=910.0)
        {
            THoEX_pp_H_AZ_bbll_CMS8=pp_H_AZ_bbll_TH8/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh);
            R_pp_H_AZ_bbll_CMS8=(1+(pp_H_AZ_bbll_TH8-ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh))/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mA>=50.0 && mA<=1000.0)
        {
            THoEX_pp_H_AZ_tautaull_CMS8=pp_H_AZ_tautaull_TH8/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh);
            R_pp_H_AZ_tautaull_CMS8=(1+(pp_H_AZ_tautaull_TH8-ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ggF_H_tautau_ATLAS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mHh);
        R_ggF_H_tautau_ATLAS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mHh))/ip_ex_gg_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mHh);
        R_bbF_H_tautau_ATLAS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mHh))/ip_ex_bb_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS13=pp_H_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mHh);
        R_pp_H_gaga_ATLAS13=(1+(pp_H_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mHh))/ip_ex_pp_phi_gaga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_llga_ATLAS13=ggF_H_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mHh);
        R_ggF_H_Zga_llga_ATLAS13=(1+(ggF_H_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mHh))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_ATLAS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mHh);
        R_pp_H_Zga_ATLAS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_llllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llnunu_CMS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llnunu_CMS13(mHh);
        R_ggF_H_ZZ_llnunu_CMS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llnunu_CMS13(mHh))/ip_ex_gg_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llnunu_CMS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llnunu_CMS13(mHh);
        R_VBF_H_ZZ_llnunu_CMS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llnunu_CMS13(mHh))/ip_ex_VBF_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llll_ATLAS13=ggF_H_ZZ_llll_TH13/ip_ex_gg_H_ZZ_llll_ATLAS13(mHh);
        R_ggF_H_ZZ_llll_ATLAS13=(1+(ggF_H_ZZ_llll_TH13-ip_ex_gg_H_ZZ_llll_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llll_ATLAS13=VBF_H_ZZ_llll_TH13/ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh);
        R_VBF_H_ZZ_llll_ATLAS13=(1+(VBF_H_ZZ_llll_TH13-ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_VBF_H_WW_lnulnu_CMS13=ggF_VBF_H_WW_lnulnu_TH13/ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh);
        R_ggF_VBF_H_WW_lnulnu_CMS13=(1+(ggF_VBF_H_WW_lnulnu_TH13-ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh))/ip_ex_ggVV_H_WW_lnulnu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbgaga_CMS13=pp_H_hh_gagabb_TH13/ip_ex_pp_H_hh_gagabb_CMS13(mHh);
        R_pp_H_hh_bbgaga_CMS13=(1+(pp_H_hh_gagabb_TH13-ip_ex_pp_H_hh_gagabb_CMS13(mHh))/ip_ex_pp_H_hh_gagabb_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbtautau_CMS13=pp_H_hh_bbtautau_TH13/ip_ex_pp_H_hh_bbtautau_CMS13(mHh);
        R_pp_H_hh_bbtautau_CMS13=(1+(pp_H_hh_bbtautau_TH13-ip_ex_pp_H_hh_bbtautau_CMS13(mHh))/ip_ex_pp_H_hh_bbtautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbtautau1_CMS13=pp_H_hh_bbtautau_TH13/ip_ex_pp_H_hh_bbtautau1_CMS13(mHh);
        R_pp_H_hh_bbtautau1_CMS13=(1+(pp_H_hh_bbtautau_TH13-ip_ex_pp_H_hh_bbtautau1_CMS13(mHh))/ip_ex_pp_H_hh_bbtautau1_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=260.0 && mHh<270.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mHh);
        R_ggF_H_tautau_CMS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mHh))/ip_ex_gg_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mHh);
        R_bbF_H_tautau_CMS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mHh))/ip_ex_bb_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS8=pp_H_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mHh);
        R_pp_H_gaga_ATLAS8=(1+(pp_H_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mHh))/ip_ex_pp_phi_gaga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS8=ggF_H_gaga_TH8/ip_ex_gg_phi_gaga_CMS8(mHh);
        R_ggF_H_gaga_CMS8=(1+(ggF_H_gaga_TH8-ip_ex_gg_phi_gaga_CMS8(mHh))/ip_ex_gg_phi_gaga_CMS8_e(mHh) ) * nftos;
//    LIMIT_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8(mHh);
//    LIMEST_ggF_H_gaga_CMS8=0.0;
//    DEVIATION_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8(mHh)-ip_ex_gg_phi_gaga_CMS8_e(mHh);
//    BANDSIZE_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8_ep2(mHh)-ip_ex_gg_phi_gaga_CMS8_em2(mHh);
        THoEX_pp_H_Zga_llga_CMS8=pp_H_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mHh);
        R_pp_H_Zga_llga_CMS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mHh))/ip_ex_pp_A_Zga_llga_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_ATLAS8=pp_H_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mHh);
        R_pp_H_Zga_llga_ATLAS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_ATLAS8=ggF_H_ZZ_TH8/ip_ex_gg_H_ZZ_ATLAS8(mHh);
        R_ggF_H_ZZ_ATLAS8=(1+(ggF_H_ZZ_TH8-ip_ex_gg_H_ZZ_ATLAS8(mHh))/ip_ex_gg_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_ATLAS8=VBF_H_ZZ_TH8/ip_ex_VBF_H_ZZ_ATLAS8(mHh);
        R_VBF_H_ZZ_ATLAS8=(1+(VBF_H_ZZ_TH8-ip_ex_VBF_H_ZZ_ATLAS8(mHh))/ip_ex_VBF_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_mu_pp_H_VV_CMS8=mu_pp_H_VV_TH8/ip_ex_mu_pp_H_VV_CMS8(mHh);
        R_mu_pp_H_VV_CMS8=(1+(mu_pp_H_VV_TH8-ip_ex_mu_pp_H_VV_CMS8(mHh))/ip_ex_mu_pp_H_VV_CMS8_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_ATLAS8=ggF_H_hh_TH8/ip_ex_gg_H_hh_ATLAS8(mHh);
        R_ggF_H_hh_ATLAS8=(1+(ggF_H_hh_TH8-ip_ex_gg_H_hh_ATLAS8(mHh))/ip_ex_gg_H_hh_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_bbtautau_CMS8=ggF_H_hh_bbtautau_TH8/ip_ex_gg_H_hh_bbtautau_CMS8(mHh);
        R_ggF_H_hh_bbtautau_CMS8=(1+(ggF_H_hh_bbtautau_TH8-ip_ex_gg_H_hh_bbtautau_CMS8(mHh))/ip_ex_gg_H_hh_bbtautau_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_gagabb_CMS8=pp_H_hh_gagabb_TH8/ip_ex_pp_phi_hh_gagabb_CMS8(mHh);
        R_pp_H_hh_gagabb_CMS8=(1+(pp_H_hh_gagabb_TH8-ip_ex_pp_phi_hh_gagabb_CMS8(mHh))/ip_ex_pp_phi_hh_gagabb_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_bb_CMS8=bbF_H_bb_TH8/ip_ex_bb_phi_bb_CMS8(mHh);
        R_bbF_H_bb_CMS8=(1+(bbF_H_bb_TH8-ip_ex_bb_phi_bb_CMS8(mHh))/ip_ex_bb_phi_bb_CMS8_e(mHh) ) * nftos;
        if(mA>=40.0 && mA<=910.0)
        {
            THoEX_pp_H_AZ_bbll_CMS8=pp_H_AZ_bbll_TH8/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh);
            R_pp_H_AZ_bbll_CMS8=(1+(pp_H_AZ_bbll_TH8-ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh))/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mA>=50.0 && mA<=1000.0)
        {
            THoEX_pp_H_AZ_tautaull_CMS8=pp_H_AZ_tautaull_TH8/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh);
            R_pp_H_AZ_tautaull_CMS8=(1+(pp_H_AZ_tautaull_TH8-ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ggF_H_tautau_ATLAS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mHh);
        R_ggF_H_tautau_ATLAS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mHh))/ip_ex_gg_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mHh);
        R_bbF_H_tautau_ATLAS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mHh))/ip_ex_bb_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS13=pp_H_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mHh);
        R_pp_H_gaga_ATLAS13=(1+(pp_H_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mHh))/ip_ex_pp_phi_gaga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_llga_ATLAS13=ggF_H_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mHh);
        R_ggF_H_Zga_llga_ATLAS13=(1+(ggF_H_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mHh))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_ATLAS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mHh);
        R_pp_H_Zga_ATLAS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_llllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llnunu_CMS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llnunu_CMS13(mHh);
        R_ggF_H_ZZ_llnunu_CMS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llnunu_CMS13(mHh))/ip_ex_gg_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llnunu_CMS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llnunu_CMS13(mHh);
        R_VBF_H_ZZ_llnunu_CMS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llnunu_CMS13(mHh))/ip_ex_VBF_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llll_ATLAS13=ggF_H_ZZ_llll_TH13/ip_ex_gg_H_ZZ_llll_ATLAS13(mHh);
        R_ggF_H_ZZ_llll_ATLAS13=(1+(ggF_H_ZZ_llll_TH13-ip_ex_gg_H_ZZ_llll_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llll_ATLAS13=VBF_H_ZZ_llll_TH13/ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh);
        R_VBF_H_ZZ_llll_ATLAS13=(1+(VBF_H_ZZ_llll_TH13-ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_VBF_H_WW_lnulnu_CMS13=ggF_VBF_H_WW_lnulnu_TH13/ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh);
        R_ggF_VBF_H_WW_lnulnu_CMS13=(1+(ggF_VBF_H_WW_lnulnu_TH13-ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh))/ip_ex_ggVV_H_WW_lnulnu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbgaga_CMS13=pp_H_hh_gagabb_TH13/ip_ex_pp_H_hh_gagabb_CMS13(mHh);
        R_pp_H_hh_bbgaga_CMS13=(1+(pp_H_hh_gagabb_TH13-ip_ex_pp_H_hh_gagabb_CMS13(mHh))/ip_ex_pp_H_hh_gagabb_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_gagaWW_ATLAS13=ggF_H_hh_TH13/ip_ex_pp_H_hh_gagaWW_ATLAS13(mHh);
        R_ggF_H_hh_gagaWW_ATLAS13=(1+(ggF_H_hh_TH13-ip_ex_pp_H_hh_gagaWW_ATLAS13(mHh))/ip_ex_pp_H_hh_gagaWW_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbtautau_CMS13=pp_H_hh_bbtautau_TH13/ip_ex_pp_H_hh_bbtautau_CMS13(mHh);
        R_pp_H_hh_bbtautau_CMS13=(1+(pp_H_hh_bbtautau_TH13-ip_ex_pp_H_hh_bbtautau_CMS13(mHh))/ip_ex_pp_H_hh_bbtautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbtautau1_CMS13=pp_H_hh_bbtautau_TH13/ip_ex_pp_H_hh_bbtautau1_CMS13(mHh);
        R_pp_H_hh_bbtautau1_CMS13=(1+(pp_H_hh_bbtautau_TH13-ip_ex_pp_H_hh_bbtautau1_CMS13(mHh))/ip_ex_pp_H_hh_bbtautau1_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bblnulnu_CMS13=pp_H_hh_bblnulnu_TH13/ip_ex_pp_H_hh_bblnulnu_CMS13(mHh);
        R_pp_H_hh_bblnulnu_CMS13=(1+(pp_H_hh_bblnulnu_TH13-ip_ex_pp_H_hh_bblnulnu_CMS13(mHh))/ip_ex_pp_H_hh_bblnulnu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbVV_CMS13=pp_H_hh_bbVV_TH13/ip_ex_pp_H_hh_bbVV_CMS13(mHh);
        R_pp_H_hh_bbVV_CMS13=(1+(pp_H_hh_bbVV_TH13-ip_ex_pp_H_hh_bbVV_CMS13(mHh))/ip_ex_pp_H_hh_bbVV_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=270.0 && mHh<275.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mHh);
        R_ggF_H_tautau_CMS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mHh))/ip_ex_gg_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mHh);
        R_bbF_H_tautau_CMS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mHh))/ip_ex_bb_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS8=pp_H_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mHh);
        R_pp_H_gaga_ATLAS8=(1+(pp_H_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mHh))/ip_ex_pp_phi_gaga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS8=ggF_H_gaga_TH8/ip_ex_gg_phi_gaga_CMS8(mHh);
        R_ggF_H_gaga_CMS8=(1+(ggF_H_gaga_TH8-ip_ex_gg_phi_gaga_CMS8(mHh))/ip_ex_gg_phi_gaga_CMS8_e(mHh) ) * nftos;
//    LIMIT_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8(mHh);
//    LIMEST_ggF_H_gaga_CMS8=0.0;
//    DEVIATION_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8(mHh)-ip_ex_gg_phi_gaga_CMS8_e(mHh);
//    BANDSIZE_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8_ep2(mHh)-ip_ex_gg_phi_gaga_CMS8_em2(mHh);
        THoEX_pp_H_Zga_llga_CMS8=pp_H_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mHh);
        R_pp_H_Zga_llga_CMS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mHh))/ip_ex_pp_A_Zga_llga_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_ATLAS8=pp_H_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mHh);
        R_pp_H_Zga_llga_ATLAS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_ATLAS8=ggF_H_ZZ_TH8/ip_ex_gg_H_ZZ_ATLAS8(mHh);
        R_ggF_H_ZZ_ATLAS8=(1+(ggF_H_ZZ_TH8-ip_ex_gg_H_ZZ_ATLAS8(mHh))/ip_ex_gg_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_ATLAS8=VBF_H_ZZ_TH8/ip_ex_VBF_H_ZZ_ATLAS8(mHh);
        R_VBF_H_ZZ_ATLAS8=(1+(VBF_H_ZZ_TH8-ip_ex_VBF_H_ZZ_ATLAS8(mHh))/ip_ex_VBF_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_mu_pp_H_VV_CMS8=mu_pp_H_VV_TH8/ip_ex_mu_pp_H_VV_CMS8(mHh);
        R_mu_pp_H_VV_CMS8=(1+(mu_pp_H_VV_TH8-ip_ex_mu_pp_H_VV_CMS8(mHh))/ip_ex_mu_pp_H_VV_CMS8_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_ATLAS8=ggF_H_hh_TH8/ip_ex_gg_H_hh_ATLAS8(mHh);
        R_ggF_H_hh_ATLAS8=(1+(ggF_H_hh_TH8-ip_ex_gg_H_hh_ATLAS8(mHh))/ip_ex_gg_H_hh_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_bbtautau_CMS8=ggF_H_hh_bbtautau_TH8/ip_ex_gg_H_hh_bbtautau_CMS8(mHh);
        R_ggF_H_hh_bbtautau_CMS8=(1+(ggF_H_hh_bbtautau_TH8-ip_ex_gg_H_hh_bbtautau_CMS8(mHh))/ip_ex_gg_H_hh_bbtautau_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_gagabb_CMS8=pp_H_hh_gagabb_TH8/ip_ex_pp_phi_hh_gagabb_CMS8(mHh);
        R_pp_H_hh_gagabb_CMS8=(1+(pp_H_hh_gagabb_TH8-ip_ex_pp_phi_hh_gagabb_CMS8(mHh))/ip_ex_pp_phi_hh_gagabb_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_bb_CMS8=bbF_H_bb_TH8/ip_ex_bb_phi_bb_CMS8(mHh);
        R_bbF_H_bb_CMS8=(1+(bbF_H_bb_TH8-ip_ex_bb_phi_bb_CMS8(mHh))/ip_ex_bb_phi_bb_CMS8_e(mHh) ) * nftos;
        if(mA>=40.0 && mA<=910.0)
        {
            THoEX_pp_H_AZ_bbll_CMS8=pp_H_AZ_bbll_TH8/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh);
            R_pp_H_AZ_bbll_CMS8=(1+(pp_H_AZ_bbll_TH8-ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh))/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mA>=50.0 && mA<=1000.0)
        {
            THoEX_pp_H_AZ_tautaull_CMS8=pp_H_AZ_tautaull_TH8/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh);
            R_pp_H_AZ_tautaull_CMS8=(1+(pp_H_AZ_tautaull_TH8-ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ggF_H_tautau_ATLAS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mHh);
        R_ggF_H_tautau_ATLAS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mHh))/ip_ex_gg_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mHh);
        R_bbF_H_tautau_ATLAS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mHh))/ip_ex_bb_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS13=pp_H_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mHh);
        R_pp_H_gaga_ATLAS13=(1+(pp_H_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mHh))/ip_ex_pp_phi_gaga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_llga_ATLAS13=ggF_H_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mHh);
        R_ggF_H_Zga_llga_ATLAS13=(1+(ggF_H_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mHh))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_ATLAS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mHh);
        R_pp_H_Zga_ATLAS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_llllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llnunu_CMS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llnunu_CMS13(mHh);
        R_ggF_H_ZZ_llnunu_CMS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llnunu_CMS13(mHh))/ip_ex_gg_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llnunu_CMS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llnunu_CMS13(mHh);
        R_VBF_H_ZZ_llnunu_CMS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llnunu_CMS13(mHh))/ip_ex_VBF_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llll_ATLAS13=ggF_H_ZZ_llll_TH13/ip_ex_gg_H_ZZ_llll_ATLAS13(mHh);
        R_ggF_H_ZZ_llll_ATLAS13=(1+(ggF_H_ZZ_llll_TH13-ip_ex_gg_H_ZZ_llll_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llll_ATLAS13=VBF_H_ZZ_llll_TH13/ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh);
        R_VBF_H_ZZ_llll_ATLAS13=(1+(VBF_H_ZZ_llll_TH13-ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_VBF_H_WW_lnulnu_CMS13=ggF_VBF_H_WW_lnulnu_TH13/ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh);
        R_ggF_VBF_H_WW_lnulnu_CMS13=(1+(ggF_VBF_H_WW_lnulnu_TH13-ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh))/ip_ex_ggVV_H_WW_lnulnu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbgaga_CMS13=pp_H_hh_gagabb_TH13/ip_ex_pp_H_hh_gagabb_CMS13(mHh);
        R_pp_H_hh_bbgaga_CMS13=(1+(pp_H_hh_gagabb_TH13-ip_ex_pp_H_hh_gagabb_CMS13(mHh))/ip_ex_pp_H_hh_gagabb_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_gagaWW_ATLAS13=ggF_H_hh_TH13/ip_ex_pp_H_hh_gagaWW_ATLAS13(mHh);
        R_ggF_H_hh_gagaWW_ATLAS13=(1+(ggF_H_hh_TH13-ip_ex_pp_H_hh_gagaWW_ATLAS13(mHh))/ip_ex_pp_H_hh_gagaWW_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbtautau_CMS13=pp_H_hh_bbtautau_TH13/ip_ex_pp_H_hh_bbtautau_CMS13(mHh);
        R_pp_H_hh_bbtautau_CMS13=(1+(pp_H_hh_bbtautau_TH13-ip_ex_pp_H_hh_bbtautau_CMS13(mHh))/ip_ex_pp_H_hh_bbtautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbtautau1_CMS13=pp_H_hh_bbtautau_TH13/ip_ex_pp_H_hh_bbtautau1_CMS13(mHh);
        R_pp_H_hh_bbtautau1_CMS13=(1+(pp_H_hh_bbtautau_TH13-ip_ex_pp_H_hh_bbtautau1_CMS13(mHh))/ip_ex_pp_H_hh_bbtautau1_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bblnulnu_CMS13=pp_H_hh_bblnulnu_TH13/ip_ex_pp_H_hh_bblnulnu_CMS13(mHh);
        R_pp_H_hh_bblnulnu_CMS13=(1+(pp_H_hh_bblnulnu_TH13-ip_ex_pp_H_hh_bblnulnu_CMS13(mHh))/ip_ex_pp_H_hh_bblnulnu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbVV_CMS13=pp_H_hh_bbVV_TH13/ip_ex_pp_H_hh_bbVV_CMS13(mHh);
        R_pp_H_hh_bbVV_CMS13=(1+(pp_H_hh_bbVV_TH13-ip_ex_pp_H_hh_bbVV_CMS13(mHh))/ip_ex_pp_H_hh_bbVV_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=275.0 && mHh<300.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mHh);
        R_ggF_H_tautau_CMS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mHh))/ip_ex_gg_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mHh);
        R_bbF_H_tautau_CMS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mHh))/ip_ex_bb_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS8=pp_H_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mHh);
        R_pp_H_gaga_ATLAS8=(1+(pp_H_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mHh))/ip_ex_pp_phi_gaga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS8=ggF_H_gaga_TH8/ip_ex_gg_phi_gaga_CMS8(mHh);
        R_ggF_H_gaga_CMS8=(1+(ggF_H_gaga_TH8-ip_ex_gg_phi_gaga_CMS8(mHh))/ip_ex_gg_phi_gaga_CMS8_e(mHh) ) * nftos;
//    LIMIT_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8(mHh);
//    LIMEST_ggF_H_gaga_CMS8=0.0;
//    DEVIATION_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8(mHh)-ip_ex_gg_phi_gaga_CMS8_e(mHh);
//    BANDSIZE_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8_ep2(mHh)-ip_ex_gg_phi_gaga_CMS8_em2(mHh);
        THoEX_pp_H_Zga_llga_CMS8=pp_H_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mHh);
        R_pp_H_Zga_llga_CMS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mHh))/ip_ex_pp_A_Zga_llga_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_ATLAS8=pp_H_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mHh);
        R_pp_H_Zga_llga_ATLAS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_ATLAS8=ggF_H_ZZ_TH8/ip_ex_gg_H_ZZ_ATLAS8(mHh);
        R_ggF_H_ZZ_ATLAS8=(1+(ggF_H_ZZ_TH8-ip_ex_gg_H_ZZ_ATLAS8(mHh))/ip_ex_gg_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_ATLAS8=VBF_H_ZZ_TH8/ip_ex_VBF_H_ZZ_ATLAS8(mHh);
        R_VBF_H_ZZ_ATLAS8=(1+(VBF_H_ZZ_TH8-ip_ex_VBF_H_ZZ_ATLAS8(mHh))/ip_ex_VBF_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_mu_pp_H_VV_CMS8=mu_pp_H_VV_TH8/ip_ex_mu_pp_H_VV_CMS8(mHh);
        R_mu_pp_H_VV_CMS8=(1+(mu_pp_H_VV_TH8-ip_ex_mu_pp_H_VV_CMS8(mHh))/ip_ex_mu_pp_H_VV_CMS8_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_ATLAS8=ggF_H_hh_TH8/ip_ex_gg_H_hh_ATLAS8(mHh);
        R_ggF_H_hh_ATLAS8=(1+(ggF_H_hh_TH8-ip_ex_gg_H_hh_ATLAS8(mHh))/ip_ex_gg_H_hh_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_bbtautau_CMS8=ggF_H_hh_bbtautau_TH8/ip_ex_gg_H_hh_bbtautau_CMS8(mHh);
        R_ggF_H_hh_bbtautau_CMS8=(1+(ggF_H_hh_bbtautau_TH8-ip_ex_gg_H_hh_bbtautau_CMS8(mHh))/ip_ex_gg_H_hh_bbtautau_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS8=pp_H_hh_bbbb_TH8/ip_ex_pp_phi_hh_bbbb_CMS8(mHh);
        R_pp_H_hh_bbbb_CMS8=(1+(pp_H_hh_bbbb_TH8-ip_ex_pp_phi_hh_bbbb_CMS8(mHh))/ip_ex_pp_phi_hh_bbbb_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_gagabb_CMS8=pp_H_hh_gagabb_TH8/ip_ex_pp_phi_hh_gagabb_CMS8(mHh);
        R_pp_H_hh_gagabb_CMS8=(1+(pp_H_hh_gagabb_TH8-ip_ex_pp_phi_hh_gagabb_CMS8(mHh))/ip_ex_pp_phi_hh_gagabb_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_bb_CMS8=bbF_H_bb_TH8/ip_ex_bb_phi_bb_CMS8(mHh);
        R_bbF_H_bb_CMS8=(1+(bbF_H_bb_TH8-ip_ex_bb_phi_bb_CMS8(mHh))/ip_ex_bb_phi_bb_CMS8_e(mHh) ) * nftos;
        if(mA>=40.0 && mA<=910.0)
        {
            THoEX_pp_H_AZ_bbll_CMS8=pp_H_AZ_bbll_TH8/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh);
            R_pp_H_AZ_bbll_CMS8=(1+(pp_H_AZ_bbll_TH8-ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh))/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mA>=50.0 && mA<=1000.0)
        {
            THoEX_pp_H_AZ_tautaull_CMS8=pp_H_AZ_tautaull_TH8/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh);
            R_pp_H_AZ_tautaull_CMS8=(1+(pp_H_AZ_tautaull_TH8-ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ggF_H_tautau_ATLAS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mHh);
        R_ggF_H_tautau_ATLAS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mHh))/ip_ex_gg_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mHh);
        R_bbF_H_tautau_ATLAS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mHh))/ip_ex_bb_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS13=pp_H_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mHh);
        R_pp_H_gaga_ATLAS13=(1+(pp_H_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mHh))/ip_ex_pp_phi_gaga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_llga_ATLAS13=ggF_H_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mHh);
        R_ggF_H_Zga_llga_ATLAS13=(1+(ggF_H_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mHh))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_ATLAS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mHh);
        R_pp_H_Zga_ATLAS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_llllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llnunu_CMS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llnunu_CMS13(mHh);
        R_ggF_H_ZZ_llnunu_CMS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llnunu_CMS13(mHh))/ip_ex_gg_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llnunu_CMS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llnunu_CMS13(mHh);
        R_VBF_H_ZZ_llnunu_CMS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llnunu_CMS13(mHh))/ip_ex_VBF_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llll_ATLAS13=ggF_H_ZZ_llll_TH13/ip_ex_gg_H_ZZ_llll_ATLAS13(mHh);
        R_ggF_H_ZZ_llll_ATLAS13=(1+(ggF_H_ZZ_llll_TH13-ip_ex_gg_H_ZZ_llll_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llll_ATLAS13=VBF_H_ZZ_llll_TH13/ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh);
        R_VBF_H_ZZ_llll_ATLAS13=(1+(VBF_H_ZZ_llll_TH13-ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_VBF_H_WW_lnulnu_CMS13=ggF_VBF_H_WW_lnulnu_TH13/ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh);
        R_ggF_VBF_H_WW_lnulnu_CMS13=(1+(ggF_VBF_H_WW_lnulnu_TH13-ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh))/ip_ex_ggVV_H_WW_lnulnu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbgaga_ATLAS13=pp_H_hh_TH13/ip_ex_pp_H_hh_gagabb_ATLAS13(mHh);
        R_pp_H_hh_bbgaga_ATLAS13=(1+(pp_H_hh_TH13-ip_ex_pp_H_hh_gagabb_ATLAS13(mHh))/ip_ex_pp_H_hh_gagabb_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbgaga_CMS13=pp_H_hh_gagabb_TH13/ip_ex_pp_H_hh_gagabb_CMS13(mHh);
        R_pp_H_hh_bbgaga_CMS13=(1+(pp_H_hh_gagabb_TH13-ip_ex_pp_H_hh_gagabb_CMS13(mHh))/ip_ex_pp_H_hh_gagabb_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_gagaWW_ATLAS13=ggF_H_hh_TH13/ip_ex_pp_H_hh_gagaWW_ATLAS13(mHh);
        R_ggF_H_hh_gagaWW_ATLAS13=(1+(ggF_H_hh_TH13-ip_ex_pp_H_hh_gagaWW_ATLAS13(mHh))/ip_ex_pp_H_hh_gagaWW_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbtautau_CMS13=pp_H_hh_bbtautau_TH13/ip_ex_pp_H_hh_bbtautau_CMS13(mHh);
        R_pp_H_hh_bbtautau_CMS13=(1+(pp_H_hh_bbtautau_TH13-ip_ex_pp_H_hh_bbtautau_CMS13(mHh))/ip_ex_pp_H_hh_bbtautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbtautau1_CMS13=pp_H_hh_bbtautau_TH13/ip_ex_pp_H_hh_bbtautau1_CMS13(mHh);
        R_pp_H_hh_bbtautau1_CMS13=(1+(pp_H_hh_bbtautau_TH13-ip_ex_pp_H_hh_bbtautau1_CMS13(mHh))/ip_ex_pp_H_hh_bbtautau1_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bblnulnu_CMS13=pp_H_hh_bblnulnu_TH13/ip_ex_pp_H_hh_bblnulnu_CMS13(mHh);
        R_pp_H_hh_bblnulnu_CMS13=(1+(pp_H_hh_bblnulnu_TH13-ip_ex_pp_H_hh_bblnulnu_CMS13(mHh))/ip_ex_pp_H_hh_bblnulnu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbVV_CMS13=pp_H_hh_bbVV_TH13/ip_ex_pp_H_hh_bbVV_CMS13(mHh);
        R_pp_H_hh_bbVV_CMS13=(1+(pp_H_hh_bbVV_TH13-ip_ex_pp_H_hh_bbVV_CMS13(mHh))/ip_ex_pp_H_hh_bbVV_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=300.0 && mHh<350.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mHh);
        R_ggF_H_tautau_CMS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mHh))/ip_ex_gg_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mHh);
        R_bbF_H_tautau_CMS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mHh))/ip_ex_bb_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS8=pp_H_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mHh);
        R_pp_H_gaga_ATLAS8=(1+(pp_H_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mHh))/ip_ex_pp_phi_gaga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS8=ggF_H_gaga_TH8/ip_ex_gg_phi_gaga_CMS8(mHh);
        R_ggF_H_gaga_CMS8=(1+(ggF_H_gaga_TH8-ip_ex_gg_phi_gaga_CMS8(mHh))/ip_ex_gg_phi_gaga_CMS8_e(mHh) ) * nftos;
//    LIMIT_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8(mHh);
//    LIMEST_ggF_H_gaga_CMS8=0.0;
//    DEVIATION_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8(mHh)-ip_ex_gg_phi_gaga_CMS8_e(mHh);
//    BANDSIZE_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8_ep2(mHh)-ip_ex_gg_phi_gaga_CMS8_em2(mHh);
        THoEX_pp_H_Zga_llga_CMS8=pp_H_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mHh);
        R_pp_H_Zga_llga_CMS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mHh))/ip_ex_pp_A_Zga_llga_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_ATLAS8=pp_H_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mHh);
        R_pp_H_Zga_llga_ATLAS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_ATLAS8=ggF_H_ZZ_TH8/ip_ex_gg_H_ZZ_ATLAS8(mHh);
        R_ggF_H_ZZ_ATLAS8=(1+(ggF_H_ZZ_TH8-ip_ex_gg_H_ZZ_ATLAS8(mHh))/ip_ex_gg_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_ATLAS8=VBF_H_ZZ_TH8/ip_ex_VBF_H_ZZ_ATLAS8(mHh);
        R_VBF_H_ZZ_ATLAS8=(1+(VBF_H_ZZ_TH8-ip_ex_VBF_H_ZZ_ATLAS8(mHh))/ip_ex_VBF_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_mu_pp_H_VV_CMS8=mu_pp_H_VV_TH8/ip_ex_mu_pp_H_VV_CMS8(mHh);
        R_mu_pp_H_VV_CMS8=(1+(mu_pp_H_VV_TH8-ip_ex_mu_pp_H_VV_CMS8(mHh))/ip_ex_mu_pp_H_VV_CMS8_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_ATLAS8=ggF_H_WW_TH8/ip_ex_gg_H_WW_ATLAS8(mHh);
        R_ggF_H_WW_ATLAS8=(1+(ggF_H_WW_TH8-ip_ex_gg_H_WW_ATLAS8(mHh))/ip_ex_gg_H_WW_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_ATLAS8=VBF_H_WW_TH8/ip_ex_VBF_H_WW_ATLAS8(mHh);
        R_VBF_H_WW_ATLAS8=(1+(VBF_H_WW_TH8-ip_ex_VBF_H_WW_ATLAS8(mHh))/ip_ex_VBF_H_WW_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_ATLAS8=ggF_H_hh_TH8/ip_ex_gg_H_hh_ATLAS8(mHh);
        R_ggF_H_hh_ATLAS8=(1+(ggF_H_hh_TH8-ip_ex_gg_H_hh_ATLAS8(mHh))/ip_ex_gg_H_hh_ATLAS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_CMS8=pp_H_hh_TH8/ip_ex_pp_H_hh_CMS8(mHh);
        R_pp_H_hh_CMS8=(1+(pp_H_hh_TH8-ip_ex_pp_H_hh_CMS8(mHh))/ip_ex_pp_H_hh_CMS8_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_bbtautau_CMS8=ggF_H_hh_bbtautau_TH8/ip_ex_gg_H_hh_bbtautau_CMS8(mHh);
        R_ggF_H_hh_bbtautau_CMS8=(1+(ggF_H_hh_bbtautau_TH8-ip_ex_gg_H_hh_bbtautau_CMS8(mHh))/ip_ex_gg_H_hh_bbtautau_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS8=pp_H_hh_bbbb_TH8/ip_ex_pp_phi_hh_bbbb_CMS8(mHh);
        R_pp_H_hh_bbbb_CMS8=(1+(pp_H_hh_bbbb_TH8-ip_ex_pp_phi_hh_bbbb_CMS8(mHh))/ip_ex_pp_phi_hh_bbbb_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_gagabb_CMS8=pp_H_hh_gagabb_TH8/ip_ex_pp_phi_hh_gagabb_CMS8(mHh);
        R_pp_H_hh_gagabb_CMS8=(1+(pp_H_hh_gagabb_TH8-ip_ex_pp_phi_hh_gagabb_CMS8(mHh))/ip_ex_pp_phi_hh_gagabb_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_bb_CMS8=bbF_H_bb_TH8/ip_ex_bb_phi_bb_CMS8(mHh);
        R_bbF_H_bb_CMS8=(1+(bbF_H_bb_TH8-ip_ex_bb_phi_bb_CMS8(mHh))/ip_ex_bb_phi_bb_CMS8_e(mHh) ) * nftos;
        if(mA>=40.0 && mA<=910.0)
        {
            THoEX_pp_H_AZ_bbll_CMS8=pp_H_AZ_bbll_TH8/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh);
            R_pp_H_AZ_bbll_CMS8=(1+(pp_H_AZ_bbll_TH8-ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh))/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mA>=50.0 && mA<=1000.0)
        {
            THoEX_pp_H_AZ_tautaull_CMS8=pp_H_AZ_tautaull_TH8/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh);
            R_pp_H_AZ_tautaull_CMS8=(1+(pp_H_AZ_tautaull_TH8-ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ggF_H_tautau_ATLAS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mHh);
        R_ggF_H_tautau_ATLAS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mHh))/ip_ex_gg_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mHh);
        R_bbF_H_tautau_ATLAS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mHh))/ip_ex_bb_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS13=pp_H_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mHh);
        R_pp_H_gaga_ATLAS13=(1+(pp_H_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mHh))/ip_ex_pp_phi_gaga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_llga_ATLAS13=ggF_H_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mHh);
        R_ggF_H_Zga_llga_ATLAS13=(1+(ggF_H_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mHh))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_ATLAS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mHh);
        R_pp_H_Zga_ATLAS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_CMS13(mHh);
        R_pp_H_Zga_llga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_CMS13(mHh))/ip_ex_pp_phi_Zga_llga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_llllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llnunu_CMS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llnunu_CMS13(mHh);
        R_ggF_H_ZZ_llnunu_CMS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llnunu_CMS13(mHh))/ip_ex_gg_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llnunu_CMS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llnunu_CMS13(mHh);
        R_VBF_H_ZZ_llnunu_CMS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llnunu_CMS13(mHh))/ip_ex_VBF_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llll_ATLAS13=ggF_H_ZZ_llll_TH13/ip_ex_gg_H_ZZ_llll_ATLAS13(mHh);
        R_ggF_H_ZZ_llll_ATLAS13=(1+(ggF_H_ZZ_llll_TH13-ip_ex_gg_H_ZZ_llll_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llll_ATLAS13=VBF_H_ZZ_llll_TH13/ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh);
        R_VBF_H_ZZ_llll_ATLAS13=(1+(VBF_H_ZZ_llll_TH13-ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_qqllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_qqllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_qqllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_qqllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh);
        R_ggF_H_ZZ_llqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llqq_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh);
        R_VBF_H_ZZ_llqq_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_lnuqq_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh);
        R_ggF_H_WW_lnuqq_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_gg_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_lnuqq_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh);
        R_VBF_H_WW_lnuqq_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_VBF_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_enumunu_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_enumunu_ATLAS13(mHh);
        R_ggF_H_WW_enumunu_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_enumunu_ATLAS13(mHh))/ip_ex_gg_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_enumunu_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh);
        R_VBF_H_WW_enumunu_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh))/ip_ex_VBF_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_VBF_H_WW_lnulnu_CMS13=ggF_VBF_H_WW_lnulnu_TH13/ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh);
        R_ggF_VBF_H_WW_lnulnu_CMS13=(1+(ggF_VBF_H_WW_lnulnu_TH13-ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh))/ip_ex_ggVV_H_WW_lnulnu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbgaga_ATLAS13=pp_H_hh_TH13/ip_ex_pp_H_hh_gagabb_ATLAS13(mHh);
        R_pp_H_hh_bbgaga_ATLAS13=(1+(pp_H_hh_TH13-ip_ex_pp_H_hh_gagabb_ATLAS13(mHh))/ip_ex_pp_H_hh_gagabb_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbgaga_CMS13=pp_H_hh_gagabb_TH13/ip_ex_pp_H_hh_gagabb_CMS13(mHh);
        R_pp_H_hh_bbgaga_CMS13=(1+(pp_H_hh_gagabb_TH13-ip_ex_pp_H_hh_gagabb_CMS13(mHh))/ip_ex_pp_H_hh_gagabb_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_ATLAS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_ATLAS13(mHh);
        R_pp_H_hh_bbbb_ATLAS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_ATLAS13(mHh))/ip_ex_pp_H_hh_bbbb_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_gagaWW_ATLAS13=ggF_H_hh_TH13/ip_ex_pp_H_hh_gagaWW_ATLAS13(mHh);
        R_ggF_H_hh_gagaWW_ATLAS13=(1+(ggF_H_hh_TH13-ip_ex_pp_H_hh_gagaWW_ATLAS13(mHh))/ip_ex_pp_H_hh_gagaWW_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbtautau_CMS13=pp_H_hh_bbtautau_TH13/ip_ex_pp_H_hh_bbtautau_CMS13(mHh);
        R_pp_H_hh_bbtautau_CMS13=(1+(pp_H_hh_bbtautau_TH13-ip_ex_pp_H_hh_bbtautau_CMS13(mHh))/ip_ex_pp_H_hh_bbtautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbtautau1_CMS13=pp_H_hh_bbtautau_TH13/ip_ex_pp_H_hh_bbtautau1_CMS13(mHh);
        R_pp_H_hh_bbtautau1_CMS13=(1+(pp_H_hh_bbtautau_TH13-ip_ex_pp_H_hh_bbtautau1_CMS13(mHh))/ip_ex_pp_H_hh_bbtautau1_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bblnulnu_CMS13=pp_H_hh_bblnulnu_TH13/ip_ex_pp_H_hh_bblnulnu_CMS13(mHh);
        R_pp_H_hh_bblnulnu_CMS13=(1+(pp_H_hh_bblnulnu_TH13-ip_ex_pp_H_hh_bblnulnu_CMS13(mHh))/ip_ex_pp_H_hh_bblnulnu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbVV_CMS13=pp_H_hh_bbVV_TH13/ip_ex_pp_H_hh_bbVV_CMS13(mHh);
        R_pp_H_hh_bbVV_CMS13=(1+(pp_H_hh_bbVV_TH13-ip_ex_pp_H_hh_bbVV_CMS13(mHh))/ip_ex_pp_H_hh_bbVV_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=350.0 && mHh<400.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mHh);
        R_ggF_H_tautau_CMS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mHh))/ip_ex_gg_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mHh);
        R_bbF_H_tautau_CMS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mHh))/ip_ex_bb_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS8=pp_H_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mHh);
        R_pp_H_gaga_ATLAS8=(1+(pp_H_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mHh))/ip_ex_pp_phi_gaga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS8=ggF_H_gaga_TH8/ip_ex_gg_phi_gaga_CMS8(mHh);
        R_ggF_H_gaga_CMS8=(1+(ggF_H_gaga_TH8-ip_ex_gg_phi_gaga_CMS8(mHh))/ip_ex_gg_phi_gaga_CMS8_e(mHh) ) * nftos;
//    LIMIT_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8(mHh);
//    LIMEST_ggF_H_gaga_CMS8=0.0;
//    DEVIATION_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8(mHh)-ip_ex_gg_phi_gaga_CMS8_e(mHh);
//    BANDSIZE_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8_ep2(mHh)-ip_ex_gg_phi_gaga_CMS8_em2(mHh);
        THoEX_pp_H_Zga_llga_CMS8=pp_H_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mHh);
        R_pp_H_Zga_llga_CMS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mHh))/ip_ex_pp_A_Zga_llga_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_ATLAS8=pp_H_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mHh);
        R_pp_H_Zga_llga_ATLAS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_ATLAS8=ggF_H_ZZ_TH8/ip_ex_gg_H_ZZ_ATLAS8(mHh);
        R_ggF_H_ZZ_ATLAS8=(1+(ggF_H_ZZ_TH8-ip_ex_gg_H_ZZ_ATLAS8(mHh))/ip_ex_gg_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_ATLAS8=VBF_H_ZZ_TH8/ip_ex_VBF_H_ZZ_ATLAS8(mHh);
        R_VBF_H_ZZ_ATLAS8=(1+(VBF_H_ZZ_TH8-ip_ex_VBF_H_ZZ_ATLAS8(mHh))/ip_ex_VBF_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_mu_pp_H_VV_CMS8=mu_pp_H_VV_TH8/ip_ex_mu_pp_H_VV_CMS8(mHh);
        R_mu_pp_H_VV_CMS8=(1+(mu_pp_H_VV_TH8-ip_ex_mu_pp_H_VV_CMS8(mHh))/ip_ex_mu_pp_H_VV_CMS8_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_ATLAS8=ggF_H_WW_TH8/ip_ex_gg_H_WW_ATLAS8(mHh);
        R_ggF_H_WW_ATLAS8=(1+(ggF_H_WW_TH8-ip_ex_gg_H_WW_ATLAS8(mHh))/ip_ex_gg_H_WW_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_ATLAS8=VBF_H_WW_TH8/ip_ex_VBF_H_WW_ATLAS8(mHh);
        R_VBF_H_WW_ATLAS8=(1+(VBF_H_WW_TH8-ip_ex_VBF_H_WW_ATLAS8(mHh))/ip_ex_VBF_H_WW_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_ATLAS8=ggF_H_hh_TH8/ip_ex_gg_H_hh_ATLAS8(mHh);
        R_ggF_H_hh_ATLAS8=(1+(ggF_H_hh_TH8-ip_ex_gg_H_hh_ATLAS8(mHh))/ip_ex_gg_H_hh_ATLAS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_CMS8=pp_H_hh_TH8/ip_ex_pp_H_hh_CMS8(mHh);
        R_pp_H_hh_CMS8=(1+(pp_H_hh_TH8-ip_ex_pp_H_hh_CMS8(mHh))/ip_ex_pp_H_hh_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS8=pp_H_hh_bbbb_TH8/ip_ex_pp_phi_hh_bbbb_CMS8(mHh);
        R_pp_H_hh_bbbb_CMS8=(1+(pp_H_hh_bbbb_TH8-ip_ex_pp_phi_hh_bbbb_CMS8(mHh))/ip_ex_pp_phi_hh_bbbb_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_gagabb_CMS8=pp_H_hh_gagabb_TH8/ip_ex_pp_phi_hh_gagabb_CMS8(mHh);
        R_pp_H_hh_gagabb_CMS8=(1+(pp_H_hh_gagabb_TH8-ip_ex_pp_phi_hh_gagabb_CMS8(mHh))/ip_ex_pp_phi_hh_gagabb_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_bb_CMS8=bbF_H_bb_TH8/ip_ex_bb_phi_bb_CMS8(mHh);
        R_bbF_H_bb_CMS8=(1+(bbF_H_bb_TH8-ip_ex_bb_phi_bb_CMS8(mHh))/ip_ex_bb_phi_bb_CMS8_e(mHh) ) * nftos;
        if(mA>=40.0 && mA<=910.0)
        {
            THoEX_pp_H_AZ_bbll_CMS8=pp_H_AZ_bbll_TH8/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh);
            R_pp_H_AZ_bbll_CMS8=(1+(pp_H_AZ_bbll_TH8-ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh))/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mA>=50.0 && mA<=1000.0)
        {
            THoEX_pp_H_AZ_tautaull_CMS8=pp_H_AZ_tautaull_TH8/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh);
            R_pp_H_AZ_tautaull_CMS8=(1+(pp_H_AZ_tautaull_TH8-ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ggF_H_tautau_ATLAS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mHh);
        R_ggF_H_tautau_ATLAS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mHh))/ip_ex_gg_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mHh);
        R_bbF_H_tautau_ATLAS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mHh))/ip_ex_bb_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS13=pp_H_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mHh);
        R_pp_H_gaga_ATLAS13=(1+(pp_H_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mHh))/ip_ex_pp_phi_gaga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_llga_ATLAS13=ggF_H_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mHh);
        R_ggF_H_Zga_llga_ATLAS13=(1+(ggF_H_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mHh))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_ATLAS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mHh);
        R_pp_H_Zga_ATLAS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_CMS13(mHh);
        R_pp_H_Zga_llga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_CMS13(mHh))/ip_ex_pp_phi_Zga_llga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_CMS13=ggF_H_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mHh);
        R_ggF_H_Zga_CMS13=(1+(ggF_H_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mHh))/ip_ex_ggF_phi_Zga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_llllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llnunu_CMS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llnunu_CMS13(mHh);
        R_ggF_H_ZZ_llnunu_CMS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llnunu_CMS13(mHh))/ip_ex_gg_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llnunu_CMS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llnunu_CMS13(mHh);
        R_VBF_H_ZZ_llnunu_CMS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llnunu_CMS13(mHh))/ip_ex_VBF_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llll_ATLAS13=ggF_H_ZZ_llll_TH13/ip_ex_gg_H_ZZ_llll_ATLAS13(mHh);
        R_ggF_H_ZZ_llll_ATLAS13=(1+(ggF_H_ZZ_llll_TH13-ip_ex_gg_H_ZZ_llll_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llll_ATLAS13=VBF_H_ZZ_llll_TH13/ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh);
        R_VBF_H_ZZ_llll_ATLAS13=(1+(VBF_H_ZZ_llll_TH13-ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_qqllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_qqllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_qqllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_qqllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh);
        R_ggF_H_ZZ_llqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llqq_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh);
        R_VBF_H_ZZ_llqq_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_lnuqq_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh);
        R_ggF_H_WW_lnuqq_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_gg_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_lnuqq_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh);
        R_VBF_H_WW_lnuqq_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_VBF_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_enumunu_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_enumunu_ATLAS13(mHh);
        R_ggF_H_WW_enumunu_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_enumunu_ATLAS13(mHh))/ip_ex_gg_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_enumunu_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh);
        R_VBF_H_WW_enumunu_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh))/ip_ex_VBF_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_VBF_H_WW_lnulnu_CMS13=ggF_VBF_H_WW_lnulnu_TH13/ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh);
        R_ggF_VBF_H_WW_lnulnu_CMS13=(1+(ggF_VBF_H_WW_lnulnu_TH13-ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh))/ip_ex_ggVV_H_WW_lnulnu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbgaga_ATLAS13=pp_H_hh_TH13/ip_ex_pp_H_hh_gagabb_ATLAS13(mHh);
        R_pp_H_hh_bbgaga_ATLAS13=(1+(pp_H_hh_TH13-ip_ex_pp_H_hh_gagabb_ATLAS13(mHh))/ip_ex_pp_H_hh_gagabb_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbgaga_CMS13=pp_H_hh_gagabb_TH13/ip_ex_pp_H_hh_gagabb_CMS13(mHh);
        R_pp_H_hh_bbgaga_CMS13=(1+(pp_H_hh_gagabb_TH13-ip_ex_pp_H_hh_gagabb_CMS13(mHh))/ip_ex_pp_H_hh_gagabb_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_ATLAS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_ATLAS13(mHh);
        R_pp_H_hh_bbbb_ATLAS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_ATLAS13(mHh))/ip_ex_pp_H_hh_bbbb_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_gagaWW_ATLAS13=ggF_H_hh_TH13/ip_ex_pp_H_hh_gagaWW_ATLAS13(mHh);
        R_ggF_H_hh_gagaWW_ATLAS13=(1+(ggF_H_hh_TH13-ip_ex_pp_H_hh_gagaWW_ATLAS13(mHh))/ip_ex_pp_H_hh_gagaWW_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbtautau_CMS13=pp_H_hh_bbtautau_TH13/ip_ex_pp_H_hh_bbtautau_CMS13(mHh);
        R_pp_H_hh_bbtautau_CMS13=(1+(pp_H_hh_bbtautau_TH13-ip_ex_pp_H_hh_bbtautau_CMS13(mHh))/ip_ex_pp_H_hh_bbtautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbtautau1_CMS13=pp_H_hh_bbtautau_TH13/ip_ex_pp_H_hh_bbtautau1_CMS13(mHh);
        R_pp_H_hh_bbtautau1_CMS13=(1+(pp_H_hh_bbtautau_TH13-ip_ex_pp_H_hh_bbtautau1_CMS13(mHh))/ip_ex_pp_H_hh_bbtautau1_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bblnulnu_CMS13=pp_H_hh_bblnulnu_TH13/ip_ex_pp_H_hh_bblnulnu_CMS13(mHh);
        R_pp_H_hh_bblnulnu_CMS13=(1+(pp_H_hh_bblnulnu_TH13-ip_ex_pp_H_hh_bblnulnu_CMS13(mHh))/ip_ex_pp_H_hh_bblnulnu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbVV_CMS13=pp_H_hh_bbVV_TH13/ip_ex_pp_H_hh_bbVV_CMS13(mHh);
        R_pp_H_hh_bbVV_CMS13=(1+(pp_H_hh_bbVV_TH13-ip_ex_pp_H_hh_bbVV_CMS13(mHh))/ip_ex_pp_H_hh_bbVV_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=400.0 && mHh<500.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mHh);
        R_ggF_H_tautau_CMS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mHh))/ip_ex_gg_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mHh);
        R_bbF_H_tautau_CMS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mHh))/ip_ex_bb_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS8=pp_H_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mHh);
        R_pp_H_gaga_ATLAS8=(1+(pp_H_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mHh))/ip_ex_pp_phi_gaga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS8=ggF_H_gaga_TH8/ip_ex_gg_phi_gaga_CMS8(mHh);
        R_ggF_H_gaga_CMS8=(1+(ggF_H_gaga_TH8-ip_ex_gg_phi_gaga_CMS8(mHh))/ip_ex_gg_phi_gaga_CMS8_e(mHh) ) * nftos;
//    LIMIT_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8(mHh);
//    LIMEST_ggF_H_gaga_CMS8=0.0;
//    DEVIATION_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8(mHh)-ip_ex_gg_phi_gaga_CMS8_e(mHh);
//    BANDSIZE_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8_ep2(mHh)-ip_ex_gg_phi_gaga_CMS8_em2(mHh);
        THoEX_pp_H_Zga_llga_CMS8=pp_H_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mHh);
        R_pp_H_Zga_llga_CMS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mHh))/ip_ex_pp_A_Zga_llga_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_ATLAS8=pp_H_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mHh);
        R_pp_H_Zga_llga_ATLAS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_ATLAS8=ggF_H_ZZ_TH8/ip_ex_gg_H_ZZ_ATLAS8(mHh);
        R_ggF_H_ZZ_ATLAS8=(1+(ggF_H_ZZ_TH8-ip_ex_gg_H_ZZ_ATLAS8(mHh))/ip_ex_gg_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_ATLAS8=VBF_H_ZZ_TH8/ip_ex_VBF_H_ZZ_ATLAS8(mHh);
        R_VBF_H_ZZ_ATLAS8=(1+(VBF_H_ZZ_TH8-ip_ex_VBF_H_ZZ_ATLAS8(mHh))/ip_ex_VBF_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_mu_pp_H_VV_CMS8=mu_pp_H_VV_TH8/ip_ex_mu_pp_H_VV_CMS8(mHh);
        R_mu_pp_H_VV_CMS8=(1+(mu_pp_H_VV_TH8-ip_ex_mu_pp_H_VV_CMS8(mHh))/ip_ex_mu_pp_H_VV_CMS8_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_ATLAS8=ggF_H_WW_TH8/ip_ex_gg_H_WW_ATLAS8(mHh);
        R_ggF_H_WW_ATLAS8=(1+(ggF_H_WW_TH8-ip_ex_gg_H_WW_ATLAS8(mHh))/ip_ex_gg_H_WW_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_ATLAS8=VBF_H_WW_TH8/ip_ex_VBF_H_WW_ATLAS8(mHh);
        R_VBF_H_WW_ATLAS8=(1+(VBF_H_WW_TH8-ip_ex_VBF_H_WW_ATLAS8(mHh))/ip_ex_VBF_H_WW_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_ATLAS8=ggF_H_hh_TH8/ip_ex_gg_H_hh_ATLAS8(mHh);
        R_ggF_H_hh_ATLAS8=(1+(ggF_H_hh_TH8-ip_ex_gg_H_hh_ATLAS8(mHh))/ip_ex_gg_H_hh_ATLAS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_CMS8=pp_H_hh_TH8/ip_ex_pp_H_hh_CMS8(mHh);
        R_pp_H_hh_CMS8=(1+(pp_H_hh_TH8-ip_ex_pp_H_hh_CMS8(mHh))/ip_ex_pp_H_hh_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS8=pp_H_hh_bbbb_TH8/ip_ex_pp_phi_hh_bbbb_CMS8(mHh);
        R_pp_H_hh_bbbb_CMS8=(1+(pp_H_hh_bbbb_TH8-ip_ex_pp_phi_hh_bbbb_CMS8(mHh))/ip_ex_pp_phi_hh_bbbb_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_gagabb_CMS8=pp_H_hh_gagabb_TH8/ip_ex_pp_phi_hh_gagabb_CMS8(mHh);
        R_pp_H_hh_gagabb_CMS8=(1+(pp_H_hh_gagabb_TH8-ip_ex_pp_phi_hh_gagabb_CMS8(mHh))/ip_ex_pp_phi_hh_gagabb_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_bb_CMS8=bbF_H_bb_TH8/ip_ex_bb_phi_bb_CMS8(mHh);
        R_bbF_H_bb_CMS8=(1+(bbF_H_bb_TH8-ip_ex_bb_phi_bb_CMS8(mHh))/ip_ex_bb_phi_bb_CMS8_e(mHh) ) * nftos;
        if(mA>=40.0 && mA<=910.0)
        {
            THoEX_pp_H_AZ_bbll_CMS8=pp_H_AZ_bbll_TH8/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh);
            R_pp_H_AZ_bbll_CMS8=(1+(pp_H_AZ_bbll_TH8-ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh))/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mA>=50.0 && mA<=1000.0)
        {
            THoEX_pp_H_AZ_tautaull_CMS8=pp_H_AZ_tautaull_TH8/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh);
            R_pp_H_AZ_tautaull_CMS8=(1+(pp_H_AZ_tautaull_TH8-ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ttF_H_tt_ATLAS13=ttF_H_tt_TH13/ip_ex_tt_phi_tt_ATLAS13(mHh);
        R_ttF_H_tt_ATLAS13=(1+(ttF_H_tt_TH13-ip_ex_tt_phi_tt_ATLAS13(mHh))/ip_ex_tt_phi_tt_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tt_ATLAS13=bbF_H_tt_TH13/ip_ex_bb_phi_tt_ATLAS13(mHh);
        R_bbF_H_tt_ATLAS13=(1+(bbF_H_tt_TH13-ip_ex_bb_phi_tt_ATLAS13(mHh))/ip_ex_bb_phi_tt_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_ATLAS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mHh);
        R_ggF_H_tautau_ATLAS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mHh))/ip_ex_gg_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mHh);
        R_bbF_H_tautau_ATLAS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mHh))/ip_ex_bb_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS13=pp_H_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mHh);
        R_pp_H_gaga_ATLAS13=(1+(pp_H_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mHh))/ip_ex_pp_phi_gaga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_llga_ATLAS13=ggF_H_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mHh);
        R_ggF_H_Zga_llga_ATLAS13=(1+(ggF_H_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mHh))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_ATLAS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mHh);
        R_pp_H_Zga_ATLAS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_CMS13(mHh);
        R_pp_H_Zga_llga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_CMS13(mHh))/ip_ex_pp_phi_Zga_llga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_CMS13=ggF_H_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mHh);
        R_ggF_H_Zga_CMS13=(1+(ggF_H_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mHh))/ip_ex_ggF_phi_Zga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_llllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llnunu_CMS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llnunu_CMS13(mHh);
        R_ggF_H_ZZ_llnunu_CMS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llnunu_CMS13(mHh))/ip_ex_gg_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llnunu_CMS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llnunu_CMS13(mHh);
        R_VBF_H_ZZ_llnunu_CMS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llnunu_CMS13(mHh))/ip_ex_VBF_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llll_ATLAS13=ggF_H_ZZ_llll_TH13/ip_ex_gg_H_ZZ_llll_ATLAS13(mHh);
        R_ggF_H_ZZ_llll_ATLAS13=(1+(ggF_H_ZZ_llll_TH13-ip_ex_gg_H_ZZ_llll_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llll_ATLAS13=VBF_H_ZZ_llll_TH13/ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh);
        R_VBF_H_ZZ_llll_ATLAS13=(1+(VBF_H_ZZ_llll_TH13-ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_qqllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_qqllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_qqllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_qqllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh);
        R_ggF_H_ZZ_llqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llqq_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh);
        R_VBF_H_ZZ_llqq_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_lnuqq_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh);
        R_ggF_H_WW_lnuqq_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_gg_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_lnuqq_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh);
        R_VBF_H_WW_lnuqq_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_VBF_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_enumunu_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_enumunu_ATLAS13(mHh);
        R_ggF_H_WW_enumunu_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_enumunu_ATLAS13(mHh))/ip_ex_gg_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_enumunu_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh);
        R_VBF_H_WW_enumunu_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh))/ip_ex_VBF_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_VBF_H_WW_lnulnu_CMS13=ggF_VBF_H_WW_lnulnu_TH13/ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh);
        R_ggF_VBF_H_WW_lnulnu_CMS13=(1+(ggF_VBF_H_WW_lnulnu_TH13-ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh))/ip_ex_ggVV_H_WW_lnulnu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbgaga_CMS13=pp_H_hh_gagabb_TH13/ip_ex_pp_H_hh_gagabb_CMS13(mHh);
        R_pp_H_hh_bbgaga_CMS13=(1+(pp_H_hh_gagabb_TH13-ip_ex_pp_H_hh_gagabb_CMS13(mHh))/ip_ex_pp_H_hh_gagabb_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_ATLAS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_ATLAS13(mHh);
        R_pp_H_hh_bbbb_ATLAS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_ATLAS13(mHh))/ip_ex_pp_H_hh_bbbb_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_gagaWW_ATLAS13=ggF_H_hh_TH13/ip_ex_pp_H_hh_gagaWW_ATLAS13(mHh);
        R_ggF_H_hh_gagaWW_ATLAS13=(1+(ggF_H_hh_TH13-ip_ex_pp_H_hh_gagaWW_ATLAS13(mHh))/ip_ex_pp_H_hh_gagaWW_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbtautau_CMS13=pp_H_hh_bbtautau_TH13/ip_ex_pp_H_hh_bbtautau_CMS13(mHh);
        R_pp_H_hh_bbtautau_CMS13=(1+(pp_H_hh_bbtautau_TH13-ip_ex_pp_H_hh_bbtautau_CMS13(mHh))/ip_ex_pp_H_hh_bbtautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbtautau1_CMS13=pp_H_hh_bbtautau_TH13/ip_ex_pp_H_hh_bbtautau1_CMS13(mHh);
        R_pp_H_hh_bbtautau1_CMS13=(1+(pp_H_hh_bbtautau_TH13-ip_ex_pp_H_hh_bbtautau1_CMS13(mHh))/ip_ex_pp_H_hh_bbtautau1_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bblnulnu_CMS13=pp_H_hh_bblnulnu_TH13/ip_ex_pp_H_hh_bblnulnu_CMS13(mHh);
        R_pp_H_hh_bblnulnu_CMS13=(1+(pp_H_hh_bblnulnu_TH13-ip_ex_pp_H_hh_bblnulnu_CMS13(mHh))/ip_ex_pp_H_hh_bblnulnu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbVV_CMS13=pp_H_hh_bbVV_TH13/ip_ex_pp_H_hh_bbVV_CMS13(mHh);
        R_pp_H_hh_bbVV_CMS13=(1+(pp_H_hh_bbVV_TH13-ip_ex_pp_H_hh_bbVV_CMS13(mHh))/ip_ex_pp_H_hh_bbVV_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=500.0 && mHh<550.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mHh);
        R_ggF_H_tautau_CMS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mHh))/ip_ex_gg_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mHh);
        R_bbF_H_tautau_CMS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mHh))/ip_ex_bb_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS8=pp_H_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mHh);
        R_pp_H_gaga_ATLAS8=(1+(pp_H_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mHh))/ip_ex_pp_phi_gaga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS8=ggF_H_gaga_TH8/ip_ex_gg_phi_gaga_CMS8(mHh);
        R_ggF_H_gaga_CMS8=(1+(ggF_H_gaga_TH8-ip_ex_gg_phi_gaga_CMS8(mHh))/ip_ex_gg_phi_gaga_CMS8_e(mHh) ) * nftos;
//    LIMIT_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8(mHh);
//    LIMEST_ggF_H_gaga_CMS8=0.0;
//    DEVIATION_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8(mHh)-ip_ex_gg_phi_gaga_CMS8_e(mHh);
//    BANDSIZE_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8_ep2(mHh)-ip_ex_gg_phi_gaga_CMS8_em2(mHh);
        THoEX_pp_H_Zga_llga_CMS8=pp_H_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mHh);
        R_pp_H_Zga_llga_CMS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mHh))/ip_ex_pp_A_Zga_llga_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_ATLAS8=pp_H_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mHh);
        R_pp_H_Zga_llga_ATLAS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_ATLAS8=ggF_H_ZZ_TH8/ip_ex_gg_H_ZZ_ATLAS8(mHh);
        R_ggF_H_ZZ_ATLAS8=(1+(ggF_H_ZZ_TH8-ip_ex_gg_H_ZZ_ATLAS8(mHh))/ip_ex_gg_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_ATLAS8=VBF_H_ZZ_TH8/ip_ex_VBF_H_ZZ_ATLAS8(mHh);
        R_VBF_H_ZZ_ATLAS8=(1+(VBF_H_ZZ_TH8-ip_ex_VBF_H_ZZ_ATLAS8(mHh))/ip_ex_VBF_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_mu_pp_H_VV_CMS8=mu_pp_H_VV_TH8/ip_ex_mu_pp_H_VV_CMS8(mHh);
        R_mu_pp_H_VV_CMS8=(1+(mu_pp_H_VV_TH8-ip_ex_mu_pp_H_VV_CMS8(mHh))/ip_ex_mu_pp_H_VV_CMS8_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_ATLAS8=ggF_H_WW_TH8/ip_ex_gg_H_WW_ATLAS8(mHh);
        R_ggF_H_WW_ATLAS8=(1+(ggF_H_WW_TH8-ip_ex_gg_H_WW_ATLAS8(mHh))/ip_ex_gg_H_WW_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_ATLAS8=VBF_H_WW_TH8/ip_ex_VBF_H_WW_ATLAS8(mHh);
        R_VBF_H_WW_ATLAS8=(1+(VBF_H_WW_TH8-ip_ex_VBF_H_WW_ATLAS8(mHh))/ip_ex_VBF_H_WW_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_ATLAS8=ggF_H_hh_TH8/ip_ex_gg_H_hh_ATLAS8(mHh);
        R_ggF_H_hh_ATLAS8=(1+(ggF_H_hh_TH8-ip_ex_gg_H_hh_ATLAS8(mHh))/ip_ex_gg_H_hh_ATLAS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_CMS8=pp_H_hh_TH8/ip_ex_pp_H_hh_CMS8(mHh);
        R_pp_H_hh_CMS8=(1+(pp_H_hh_TH8-ip_ex_pp_H_hh_CMS8(mHh))/ip_ex_pp_H_hh_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS8=pp_H_hh_bbbb_TH8/ip_ex_pp_phi_hh_bbbb_CMS8(mHh);
        R_pp_H_hh_bbbb_CMS8=(1+(pp_H_hh_bbbb_TH8-ip_ex_pp_phi_hh_bbbb_CMS8(mHh))/ip_ex_pp_phi_hh_bbbb_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_gagabb_CMS8=pp_H_hh_gagabb_TH8/ip_ex_pp_phi_hh_gagabb_CMS8(mHh);
        R_pp_H_hh_gagabb_CMS8=(1+(pp_H_hh_gagabb_TH8-ip_ex_pp_phi_hh_gagabb_CMS8(mHh))/ip_ex_pp_phi_hh_gagabb_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_bb_CMS8=bbF_H_bb_TH8/ip_ex_bb_phi_bb_CMS8(mHh);
        R_bbF_H_bb_CMS8=(1+(bbF_H_bb_TH8-ip_ex_bb_phi_bb_CMS8(mHh))/ip_ex_bb_phi_bb_CMS8_e(mHh) ) * nftos;
        if(mA>=40.0 && mA<=910.0)
        {
            THoEX_pp_H_AZ_bbll_CMS8=pp_H_AZ_bbll_TH8/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh);
            R_pp_H_AZ_bbll_CMS8=(1+(pp_H_AZ_bbll_TH8-ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh))/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mA>=50.0 && mA<=1000.0)
        {
            THoEX_pp_H_AZ_tautaull_CMS8=pp_H_AZ_tautaull_TH8/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh);
            R_pp_H_AZ_tautaull_CMS8=(1+(pp_H_AZ_tautaull_TH8-ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ttF_H_tt_ATLAS13=ttF_H_tt_TH13/ip_ex_tt_phi_tt_ATLAS13(mHh);
        R_ttF_H_tt_ATLAS13=(1+(ttF_H_tt_TH13-ip_ex_tt_phi_tt_ATLAS13(mHh))/ip_ex_tt_phi_tt_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tt_ATLAS13=bbF_H_tt_TH13/ip_ex_bb_phi_tt_ATLAS13(mHh);
        R_bbF_H_tt_ATLAS13=(1+(bbF_H_tt_TH13-ip_ex_bb_phi_tt_ATLAS13(mHh))/ip_ex_bb_phi_tt_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_ATLAS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mHh);
        R_ggF_H_tautau_ATLAS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mHh))/ip_ex_gg_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mHh);
        R_bbF_H_tautau_ATLAS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mHh))/ip_ex_bb_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS13=pp_H_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mHh);
        R_pp_H_gaga_ATLAS13=(1+(pp_H_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mHh))/ip_ex_pp_phi_gaga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS13=ggF_H_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mHh);
        R_ggF_H_gaga_CMS13=(1+(ggF_H_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mHh))/ip_ex_gg_phi_gaga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_llga_ATLAS13=ggF_H_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mHh);
        R_ggF_H_Zga_llga_ATLAS13=(1+(ggF_H_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mHh))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_ATLAS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mHh);
        R_pp_H_Zga_ATLAS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_CMS13(mHh);
        R_pp_H_Zga_llga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_CMS13(mHh))/ip_ex_pp_phi_Zga_llga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_CMS13=ggF_H_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mHh);
        R_ggF_H_Zga_CMS13=(1+(ggF_H_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mHh))/ip_ex_ggF_phi_Zga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_llllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llnunu_CMS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llnunu_CMS13(mHh);
        R_ggF_H_ZZ_llnunu_CMS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llnunu_CMS13(mHh))/ip_ex_gg_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llnunu_CMS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llnunu_CMS13(mHh);
        R_VBF_H_ZZ_llnunu_CMS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llnunu_CMS13(mHh))/ip_ex_VBF_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llll_ATLAS13=ggF_H_ZZ_llll_TH13/ip_ex_gg_H_ZZ_llll_ATLAS13(mHh);
        R_ggF_H_ZZ_llll_ATLAS13=(1+(ggF_H_ZZ_llll_TH13-ip_ex_gg_H_ZZ_llll_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llll_ATLAS13=VBF_H_ZZ_llll_TH13/ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh);
        R_VBF_H_ZZ_llll_ATLAS13=(1+(VBF_H_ZZ_llll_TH13-ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_qqllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_qqllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_qqllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_qqllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh);
        R_ggF_H_ZZ_llqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llqq_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh);
        R_VBF_H_ZZ_llqq_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_nunuqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh);
        R_ggF_H_ZZ_nunuqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_nunuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llqq_CMS13=pp_H_ZZ_TH13/ip_ex_pp_H_ZZ_llqq_CMS13(mHh);
        R_pp_H_ZZ_llqq_CMS13=(1+(pp_H_ZZ_TH13-ip_ex_pp_H_ZZ_llqq_CMS13(mHh))/ip_ex_pp_H_ZZ_llqq_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_lnuqq_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh);
        R_ggF_H_WW_lnuqq_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_gg_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_lnuqq_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh);
        R_VBF_H_WW_lnuqq_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_VBF_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_enumunu_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_enumunu_ATLAS13(mHh);
        R_ggF_H_WW_enumunu_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_enumunu_ATLAS13(mHh))/ip_ex_gg_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_enumunu_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh);
        R_VBF_H_WW_enumunu_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh))/ip_ex_VBF_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_VBF_H_WW_lnulnu_CMS13=ggF_VBF_H_WW_lnulnu_TH13/ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh);
        R_ggF_VBF_H_WW_lnulnu_CMS13=(1+(ggF_VBF_H_WW_lnulnu_TH13-ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh))/ip_ex_ggVV_H_WW_lnulnu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbgaga_CMS13=pp_H_hh_gagabb_TH13/ip_ex_pp_H_hh_gagabb_CMS13(mHh);
        R_pp_H_hh_bbgaga_CMS13=(1+(pp_H_hh_gagabb_TH13-ip_ex_pp_H_hh_gagabb_CMS13(mHh))/ip_ex_pp_H_hh_gagabb_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_ATLAS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_ATLAS13(mHh);
        R_pp_H_hh_bbbb_ATLAS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_ATLAS13(mHh))/ip_ex_pp_H_hh_bbbb_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbtautau_CMS13=pp_H_hh_bbtautau_TH13/ip_ex_pp_H_hh_bbtautau_CMS13(mHh);
        R_pp_H_hh_bbtautau_CMS13=(1+(pp_H_hh_bbtautau_TH13-ip_ex_pp_H_hh_bbtautau_CMS13(mHh))/ip_ex_pp_H_hh_bbtautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbtautau1_CMS13=pp_H_hh_bbtautau_TH13/ip_ex_pp_H_hh_bbtautau1_CMS13(mHh);
        R_pp_H_hh_bbtautau1_CMS13=(1+(pp_H_hh_bbtautau_TH13-ip_ex_pp_H_hh_bbtautau1_CMS13(mHh))/ip_ex_pp_H_hh_bbtautau1_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bblnulnu_CMS13=pp_H_hh_bblnulnu_TH13/ip_ex_pp_H_hh_bblnulnu_CMS13(mHh);
        R_pp_H_hh_bblnulnu_CMS13=(1+(pp_H_hh_bblnulnu_TH13-ip_ex_pp_H_hh_bblnulnu_CMS13(mHh))/ip_ex_pp_H_hh_bblnulnu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbVV_CMS13=pp_H_hh_bbVV_TH13/ip_ex_pp_H_hh_bbVV_CMS13(mHh);
        R_pp_H_hh_bbVV_CMS13=(1+(pp_H_hh_bbVV_TH13-ip_ex_pp_H_hh_bbVV_CMS13(mHh))/ip_ex_pp_H_hh_bbVV_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=550.0 && mHh<600.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mHh);
        R_ggF_H_tautau_CMS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mHh))/ip_ex_gg_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mHh);
        R_bbF_H_tautau_CMS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mHh))/ip_ex_bb_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS8=pp_H_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mHh);
        R_pp_H_gaga_ATLAS8=(1+(pp_H_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mHh))/ip_ex_pp_phi_gaga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS8=ggF_H_gaga_TH8/ip_ex_gg_phi_gaga_CMS8(mHh);
        R_ggF_H_gaga_CMS8=(1+(ggF_H_gaga_TH8-ip_ex_gg_phi_gaga_CMS8(mHh))/ip_ex_gg_phi_gaga_CMS8_e(mHh) ) * nftos;
//    LIMIT_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8(mHh);
//    LIMEST_ggF_H_gaga_CMS8=0.0;
//    DEVIATION_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8(mHh)-ip_ex_gg_phi_gaga_CMS8_e(mHh);
//    BANDSIZE_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8_ep2(mHh)-ip_ex_gg_phi_gaga_CMS8_em2(mHh);
        THoEX_pp_H_Zga_llga_CMS8=pp_H_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mHh);
        R_pp_H_Zga_llga_CMS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mHh))/ip_ex_pp_A_Zga_llga_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_ATLAS8=pp_H_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mHh);
        R_pp_H_Zga_llga_ATLAS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_ATLAS8=ggF_H_ZZ_TH8/ip_ex_gg_H_ZZ_ATLAS8(mHh);
        R_ggF_H_ZZ_ATLAS8=(1+(ggF_H_ZZ_TH8-ip_ex_gg_H_ZZ_ATLAS8(mHh))/ip_ex_gg_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_ATLAS8=VBF_H_ZZ_TH8/ip_ex_VBF_H_ZZ_ATLAS8(mHh);
        R_VBF_H_ZZ_ATLAS8=(1+(VBF_H_ZZ_TH8-ip_ex_VBF_H_ZZ_ATLAS8(mHh))/ip_ex_VBF_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_mu_pp_H_VV_CMS8=mu_pp_H_VV_TH8/ip_ex_mu_pp_H_VV_CMS8(mHh);
        R_mu_pp_H_VV_CMS8=(1+(mu_pp_H_VV_TH8-ip_ex_mu_pp_H_VV_CMS8(mHh))/ip_ex_mu_pp_H_VV_CMS8_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_ATLAS8=ggF_H_WW_TH8/ip_ex_gg_H_WW_ATLAS8(mHh);
        R_ggF_H_WW_ATLAS8=(1+(ggF_H_WW_TH8-ip_ex_gg_H_WW_ATLAS8(mHh))/ip_ex_gg_H_WW_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_ATLAS8=VBF_H_WW_TH8/ip_ex_VBF_H_WW_ATLAS8(mHh);
        R_VBF_H_WW_ATLAS8=(1+(VBF_H_WW_TH8-ip_ex_VBF_H_WW_ATLAS8(mHh))/ip_ex_VBF_H_WW_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_ATLAS8=ggF_H_hh_TH8/ip_ex_gg_H_hh_ATLAS8(mHh);
        R_ggF_H_hh_ATLAS8=(1+(ggF_H_hh_TH8-ip_ex_gg_H_hh_ATLAS8(mHh))/ip_ex_gg_H_hh_ATLAS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_CMS8=pp_H_hh_TH8/ip_ex_pp_H_hh_CMS8(mHh);
        R_pp_H_hh_CMS8=(1+(pp_H_hh_TH8-ip_ex_pp_H_hh_CMS8(mHh))/ip_ex_pp_H_hh_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS8=pp_H_hh_bbbb_TH8/ip_ex_pp_phi_hh_bbbb_CMS8(mHh);
        R_pp_H_hh_bbbb_CMS8=(1+(pp_H_hh_bbbb_TH8-ip_ex_pp_phi_hh_bbbb_CMS8(mHh))/ip_ex_pp_phi_hh_bbbb_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_gagabb_CMS8=pp_H_hh_gagabb_TH8/ip_ex_pp_phi_hh_gagabb_CMS8(mHh);
        R_pp_H_hh_gagabb_CMS8=(1+(pp_H_hh_gagabb_TH8-ip_ex_pp_phi_hh_gagabb_CMS8(mHh))/ip_ex_pp_phi_hh_gagabb_CMS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tt_ATLAS8=ggF_H_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mHh);
        R_ggF_H_tt_ATLAS8=(1+(ggF_H_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mHh))/ip_ex_gg_phi_tt_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_bb_CMS8=bbF_H_bb_TH8/ip_ex_bb_phi_bb_CMS8(mHh);
        R_bbF_H_bb_CMS8=(1+(bbF_H_bb_TH8-ip_ex_bb_phi_bb_CMS8(mHh))/ip_ex_bb_phi_bb_CMS8_e(mHh) ) * nftos;
        if(mA>=40.0 && mA<=910.0)
        {
            THoEX_pp_H_AZ_bbll_CMS8=pp_H_AZ_bbll_TH8/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh);
            R_pp_H_AZ_bbll_CMS8=(1+(pp_H_AZ_bbll_TH8-ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh))/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mA>=50.0 && mA<=1000.0)
        {
            THoEX_pp_H_AZ_tautaull_CMS8=pp_H_AZ_tautaull_TH8/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh);
            R_pp_H_AZ_tautaull_CMS8=(1+(pp_H_AZ_tautaull_TH8-ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ttF_H_tt_ATLAS13=ttF_H_tt_TH13/ip_ex_tt_phi_tt_ATLAS13(mHh);
        R_ttF_H_tt_ATLAS13=(1+(ttF_H_tt_TH13-ip_ex_tt_phi_tt_ATLAS13(mHh))/ip_ex_tt_phi_tt_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tt_ATLAS13=bbF_H_tt_TH13/ip_ex_bb_phi_tt_ATLAS13(mHh);
        R_bbF_H_tt_ATLAS13=(1+(bbF_H_tt_TH13-ip_ex_bb_phi_tt_ATLAS13(mHh))/ip_ex_bb_phi_tt_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_bb_CMS13=pp_H_bb_TH13/ip_ex_pp_phi_bb_CMS13(mHh);
        R_pp_H_bb_CMS13=(1+(pp_H_bb_TH13-ip_ex_pp_phi_bb_CMS13(mHh))/ip_ex_pp_phi_bb_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_ATLAS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mHh);
        R_ggF_H_tautau_ATLAS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mHh))/ip_ex_gg_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mHh);
        R_bbF_H_tautau_ATLAS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mHh))/ip_ex_bb_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS13=pp_H_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mHh);
        R_pp_H_gaga_ATLAS13=(1+(pp_H_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mHh))/ip_ex_pp_phi_gaga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS13=ggF_H_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mHh);
        R_ggF_H_gaga_CMS13=(1+(ggF_H_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mHh))/ip_ex_gg_phi_gaga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_llga_ATLAS13=ggF_H_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mHh);
        R_ggF_H_Zga_llga_ATLAS13=(1+(ggF_H_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mHh))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_ATLAS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mHh);
        R_pp_H_Zga_ATLAS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_CMS13(mHh);
        R_pp_H_Zga_llga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_CMS13(mHh))/ip_ex_pp_phi_Zga_llga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_CMS13=ggF_H_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mHh);
        R_ggF_H_Zga_CMS13=(1+(ggF_H_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mHh))/ip_ex_ggF_phi_Zga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_llllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llnunu_CMS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llnunu_CMS13(mHh);
        R_ggF_H_ZZ_llnunu_CMS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llnunu_CMS13(mHh))/ip_ex_gg_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llnunu_CMS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llnunu_CMS13(mHh);
        R_VBF_H_ZZ_llnunu_CMS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llnunu_CMS13(mHh))/ip_ex_VBF_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llll_ATLAS13=ggF_H_ZZ_llll_TH13/ip_ex_gg_H_ZZ_llll_ATLAS13(mHh);
        R_ggF_H_ZZ_llll_ATLAS13=(1+(ggF_H_ZZ_llll_TH13-ip_ex_gg_H_ZZ_llll_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llll_ATLAS13=VBF_H_ZZ_llll_TH13/ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh);
        R_VBF_H_ZZ_llll_ATLAS13=(1+(VBF_H_ZZ_llll_TH13-ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_qqllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_qqllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_qqllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_qqllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh);
        R_ggF_H_ZZ_llqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llqq_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh);
        R_VBF_H_ZZ_llqq_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_nunuqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh);
        R_ggF_H_ZZ_nunuqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_nunuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llqq_CMS13=pp_H_ZZ_TH13/ip_ex_pp_H_ZZ_llqq_CMS13(mHh);
        R_pp_H_ZZ_llqq_CMS13=(1+(pp_H_ZZ_TH13-ip_ex_pp_H_ZZ_llqq_CMS13(mHh))/ip_ex_pp_H_ZZ_llqq_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_lnuqq_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh);
        R_ggF_H_WW_lnuqq_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_gg_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_lnuqq_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh);
        R_VBF_H_WW_lnuqq_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_VBF_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_enumunu_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_enumunu_ATLAS13(mHh);
        R_ggF_H_WW_enumunu_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_enumunu_ATLAS13(mHh))/ip_ex_gg_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_enumunu_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh);
        R_VBF_H_WW_enumunu_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh))/ip_ex_VBF_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_VBF_H_WW_lnulnu_CMS13=ggF_VBF_H_WW_lnulnu_TH13/ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh);
        R_ggF_VBF_H_WW_lnulnu_CMS13=(1+(ggF_VBF_H_WW_lnulnu_TH13-ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh))/ip_ex_ggVV_H_WW_lnulnu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbgaga_CMS13=pp_H_hh_gagabb_TH13/ip_ex_pp_H_hh_gagabb_CMS13(mHh);
        R_pp_H_hh_bbgaga_CMS13=(1+(pp_H_hh_gagabb_TH13-ip_ex_pp_H_hh_gagabb_CMS13(mHh))/ip_ex_pp_H_hh_gagabb_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_ATLAS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_ATLAS13(mHh);
        R_pp_H_hh_bbbb_ATLAS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_ATLAS13(mHh))/ip_ex_pp_H_hh_bbbb_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbtautau_CMS13=pp_H_hh_bbtautau_TH13/ip_ex_pp_H_hh_bbtautau_CMS13(mHh);
        R_pp_H_hh_bbtautau_CMS13=(1+(pp_H_hh_bbtautau_TH13-ip_ex_pp_H_hh_bbtautau_CMS13(mHh))/ip_ex_pp_H_hh_bbtautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbtautau1_CMS13=pp_H_hh_bbtautau_TH13/ip_ex_pp_H_hh_bbtautau1_CMS13(mHh);
        R_pp_H_hh_bbtautau1_CMS13=(1+(pp_H_hh_bbtautau_TH13-ip_ex_pp_H_hh_bbtautau1_CMS13(mHh))/ip_ex_pp_H_hh_bbtautau1_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bblnulnu_CMS13=pp_H_hh_bblnulnu_TH13/ip_ex_pp_H_hh_bblnulnu_CMS13(mHh);
        R_pp_H_hh_bblnulnu_CMS13=(1+(pp_H_hh_bblnulnu_TH13-ip_ex_pp_H_hh_bblnulnu_CMS13(mHh))/ip_ex_pp_H_hh_bblnulnu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbVV_CMS13=pp_H_hh_bbVV_TH13/ip_ex_pp_H_hh_bbVV_CMS13(mHh);
        R_pp_H_hh_bbVV_CMS13=(1+(pp_H_hh_bbVV_TH13-ip_ex_pp_H_hh_bbVV_CMS13(mHh))/ip_ex_pp_H_hh_bbVV_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=600.0 && mHh<650.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mHh);
        R_ggF_H_tautau_CMS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mHh))/ip_ex_gg_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mHh);
        R_bbF_H_tautau_CMS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mHh))/ip_ex_bb_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS8=ggF_H_gaga_TH8/ip_ex_gg_phi_gaga_CMS8(mHh);
        R_ggF_H_gaga_CMS8=(1+(ggF_H_gaga_TH8-ip_ex_gg_phi_gaga_CMS8(mHh))/ip_ex_gg_phi_gaga_CMS8_e(mHh) ) * nftos;
//    LIMIT_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8(mHh);
//    LIMEST_ggF_H_gaga_CMS8=0.0;
//    DEVIATION_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8(mHh)-ip_ex_gg_phi_gaga_CMS8_e(mHh);
//    BANDSIZE_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8_ep2(mHh)-ip_ex_gg_phi_gaga_CMS8_em2(mHh);
        THoEX_pp_H_Zga_llga_CMS8=pp_H_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mHh);
        R_pp_H_Zga_llga_CMS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mHh))/ip_ex_pp_A_Zga_llga_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_ATLAS8=pp_H_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mHh);
        R_pp_H_Zga_llga_ATLAS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_ATLAS8=ggF_H_ZZ_TH8/ip_ex_gg_H_ZZ_ATLAS8(mHh);
        R_ggF_H_ZZ_ATLAS8=(1+(ggF_H_ZZ_TH8-ip_ex_gg_H_ZZ_ATLAS8(mHh))/ip_ex_gg_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_ATLAS8=VBF_H_ZZ_TH8/ip_ex_VBF_H_ZZ_ATLAS8(mHh);
        R_VBF_H_ZZ_ATLAS8=(1+(VBF_H_ZZ_TH8-ip_ex_VBF_H_ZZ_ATLAS8(mHh))/ip_ex_VBF_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_mu_pp_H_VV_CMS8=mu_pp_H_VV_TH8/ip_ex_mu_pp_H_VV_CMS8(mHh);
        R_mu_pp_H_VV_CMS8=(1+(mu_pp_H_VV_TH8-ip_ex_mu_pp_H_VV_CMS8(mHh))/ip_ex_mu_pp_H_VV_CMS8_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_ATLAS8=ggF_H_WW_TH8/ip_ex_gg_H_WW_ATLAS8(mHh);
        R_ggF_H_WW_ATLAS8=(1+(ggF_H_WW_TH8-ip_ex_gg_H_WW_ATLAS8(mHh))/ip_ex_gg_H_WW_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_ATLAS8=VBF_H_WW_TH8/ip_ex_VBF_H_WW_ATLAS8(mHh);
        R_VBF_H_WW_ATLAS8=(1+(VBF_H_WW_TH8-ip_ex_VBF_H_WW_ATLAS8(mHh))/ip_ex_VBF_H_WW_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_ATLAS8=ggF_H_hh_TH8/ip_ex_gg_H_hh_ATLAS8(mHh);
        R_ggF_H_hh_ATLAS8=(1+(ggF_H_hh_TH8-ip_ex_gg_H_hh_ATLAS8(mHh))/ip_ex_gg_H_hh_ATLAS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_CMS8=pp_H_hh_TH8/ip_ex_pp_H_hh_CMS8(mHh);
        R_pp_H_hh_CMS8=(1+(pp_H_hh_TH8-ip_ex_pp_H_hh_CMS8(mHh))/ip_ex_pp_H_hh_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS8=pp_H_hh_bbbb_TH8/ip_ex_pp_phi_hh_bbbb_CMS8(mHh);
        R_pp_H_hh_bbbb_CMS8=(1+(pp_H_hh_bbbb_TH8-ip_ex_pp_phi_hh_bbbb_CMS8(mHh))/ip_ex_pp_phi_hh_bbbb_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_gagabb_CMS8=pp_H_hh_gagabb_TH8/ip_ex_pp_phi_hh_gagabb_CMS8(mHh);
        R_pp_H_hh_gagabb_CMS8=(1+(pp_H_hh_gagabb_TH8-ip_ex_pp_phi_hh_gagabb_CMS8(mHh))/ip_ex_pp_phi_hh_gagabb_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_bb_CMS8=bbF_H_bb_TH8/ip_ex_bb_phi_bb_CMS8(mHh);
        R_bbF_H_bb_CMS8=(1+(bbF_H_bb_TH8-ip_ex_bb_phi_bb_CMS8(mHh))/ip_ex_bb_phi_bb_CMS8_e(mHh) ) * nftos;
        if(mA>=40.0 && mA<=910.0)
        {
            THoEX_pp_H_AZ_bbll_CMS8=pp_H_AZ_bbll_TH8/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh);
            R_pp_H_AZ_bbll_CMS8=(1+(pp_H_AZ_bbll_TH8-ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh))/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mA>=50.0 && mA<=1000.0)
        {
            THoEX_pp_H_AZ_tautaull_CMS8=pp_H_AZ_tautaull_TH8/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh);
            R_pp_H_AZ_tautaull_CMS8=(1+(pp_H_AZ_tautaull_TH8-ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ttF_H_tt_ATLAS13=ttF_H_tt_TH13/ip_ex_tt_phi_tt_ATLAS13(mHh);
        R_ttF_H_tt_ATLAS13=(1+(ttF_H_tt_TH13-ip_ex_tt_phi_tt_ATLAS13(mHh))/ip_ex_tt_phi_tt_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tt_ATLAS13=bbF_H_tt_TH13/ip_ex_bb_phi_tt_ATLAS13(mHh);
        R_bbF_H_tt_ATLAS13=(1+(bbF_H_tt_TH13-ip_ex_bb_phi_tt_ATLAS13(mHh))/ip_ex_bb_phi_tt_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_bb_CMS13=pp_H_bb_TH13/ip_ex_pp_phi_bb_CMS13(mHh);
        R_pp_H_bb_CMS13=(1+(pp_H_bb_TH13-ip_ex_pp_phi_bb_CMS13(mHh))/ip_ex_pp_phi_bb_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_ATLAS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mHh);
        R_ggF_H_tautau_ATLAS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mHh))/ip_ex_gg_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mHh);
        R_bbF_H_tautau_ATLAS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mHh))/ip_ex_bb_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS13=pp_H_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mHh);
        R_pp_H_gaga_ATLAS13=(1+(pp_H_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mHh))/ip_ex_pp_phi_gaga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS13=ggF_H_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mHh);
        R_ggF_H_gaga_CMS13=(1+(ggF_H_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mHh))/ip_ex_gg_phi_gaga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_llga_ATLAS13=ggF_H_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mHh);
        R_ggF_H_Zga_llga_ATLAS13=(1+(ggF_H_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mHh))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_ATLAS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mHh);
        R_pp_H_Zga_ATLAS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_CMS13(mHh);
        R_pp_H_Zga_llga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_CMS13(mHh))/ip_ex_pp_phi_Zga_llga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_CMS13=ggF_H_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mHh);
        R_ggF_H_Zga_CMS13=(1+(ggF_H_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mHh))/ip_ex_ggF_phi_Zga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_llllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llnunu_CMS13=pp_H_ZZ_TH13/ip_ex_pp_H_ZZ_llnunu_CMS13(mHh);
        R_pp_H_ZZ_llnunu_CMS13=(1+(pp_H_ZZ_TH13-ip_ex_pp_H_ZZ_llnunu_CMS13(mHh))/ip_ex_pp_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llll_ATLAS13=ggF_H_ZZ_llll_TH13/ip_ex_gg_H_ZZ_llll_ATLAS13(mHh);
        R_ggF_H_ZZ_llll_ATLAS13=(1+(ggF_H_ZZ_llll_TH13-ip_ex_gg_H_ZZ_llll_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llll_ATLAS13=VBF_H_ZZ_llll_TH13/ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh);
        R_VBF_H_ZZ_llll_ATLAS13=(1+(VBF_H_ZZ_llll_TH13-ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_qqllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_qqllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_qqllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_qqllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh);
        R_ggF_H_ZZ_llqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llqq_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh);
        R_VBF_H_ZZ_llqq_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_nunuqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh);
        R_ggF_H_ZZ_nunuqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_nunuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llqq_CMS13=pp_H_ZZ_TH13/ip_ex_pp_H_ZZ_llqq_CMS13(mHh);
        R_pp_H_ZZ_llqq_CMS13=(1+(pp_H_ZZ_TH13-ip_ex_pp_H_ZZ_llqq_CMS13(mHh))/ip_ex_pp_H_ZZ_llqq_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_lnuqq_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh);
        R_ggF_H_WW_lnuqq_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_gg_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_lnuqq_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh);
        R_VBF_H_WW_lnuqq_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_VBF_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_enumunu_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_enumunu_ATLAS13(mHh);
        R_ggF_H_WW_enumunu_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_enumunu_ATLAS13(mHh))/ip_ex_gg_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_enumunu_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh);
        R_VBF_H_WW_enumunu_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh))/ip_ex_VBF_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_VBF_H_WW_lnulnu_CMS13=ggF_VBF_H_WW_lnulnu_TH13/ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh);
        R_ggF_VBF_H_WW_lnulnu_CMS13=(1+(ggF_VBF_H_WW_lnulnu_TH13-ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh))/ip_ex_ggVV_H_WW_lnulnu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbgaga_CMS13=pp_H_hh_gagabb_TH13/ip_ex_pp_H_hh_gagabb_CMS13(mHh);
        R_pp_H_hh_bbgaga_CMS13=(1+(pp_H_hh_gagabb_TH13-ip_ex_pp_H_hh_gagabb_CMS13(mHh))/ip_ex_pp_H_hh_gagabb_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_ATLAS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_ATLAS13(mHh);
        R_pp_H_hh_bbbb_ATLAS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_ATLAS13(mHh))/ip_ex_pp_H_hh_bbbb_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbtautau_CMS13=pp_H_hh_bbtautau_TH13/ip_ex_pp_H_hh_bbtautau_CMS13(mHh);
        R_pp_H_hh_bbtautau_CMS13=(1+(pp_H_hh_bbtautau_TH13-ip_ex_pp_H_hh_bbtautau_CMS13(mHh))/ip_ex_pp_H_hh_bbtautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbtautau1_CMS13=pp_H_hh_bbtautau_TH13/ip_ex_pp_H_hh_bbtautau1_CMS13(mHh);
        R_pp_H_hh_bbtautau1_CMS13=(1+(pp_H_hh_bbtautau_TH13-ip_ex_pp_H_hh_bbtautau1_CMS13(mHh))/ip_ex_pp_H_hh_bbtautau1_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bblnulnu_CMS13=pp_H_hh_bblnulnu_TH13/ip_ex_pp_H_hh_bblnulnu_CMS13(mHh);
        R_pp_H_hh_bblnulnu_CMS13=(1+(pp_H_hh_bblnulnu_TH13-ip_ex_pp_H_hh_bblnulnu_CMS13(mHh))/ip_ex_pp_H_hh_bblnulnu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbVV_CMS13=pp_H_hh_bbVV_TH13/ip_ex_pp_H_hh_bbVV_CMS13(mHh);
        R_pp_H_hh_bbVV_CMS13=(1+(pp_H_hh_bbVV_TH13-ip_ex_pp_H_hh_bbVV_CMS13(mHh))/ip_ex_pp_H_hh_bbVV_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=650.0 && mHh<760.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mHh);
        R_ggF_H_tautau_CMS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mHh))/ip_ex_gg_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mHh);
        R_bbF_H_tautau_CMS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mHh))/ip_ex_bb_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS8=ggF_H_gaga_TH8/ip_ex_gg_phi_gaga_CMS8(mHh);
        R_ggF_H_gaga_CMS8=(1+(ggF_H_gaga_TH8-ip_ex_gg_phi_gaga_CMS8(mHh))/ip_ex_gg_phi_gaga_CMS8_e(mHh) ) * nftos;
//    LIMIT_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8(mHh);
//    LIMEST_ggF_H_gaga_CMS8=0.0;
//    DEVIATION_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8(mHh)-ip_ex_gg_phi_gaga_CMS8_e(mHh);
//    BANDSIZE_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8_ep2(mHh)-ip_ex_gg_phi_gaga_CMS8_em2(mHh);
        THoEX_pp_H_Zga_llga_CMS8=pp_H_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mHh);
        R_pp_H_Zga_llga_CMS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mHh))/ip_ex_pp_A_Zga_llga_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_ATLAS8=pp_H_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mHh);
        R_pp_H_Zga_llga_ATLAS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_ATLAS8=ggF_H_ZZ_TH8/ip_ex_gg_H_ZZ_ATLAS8(mHh);
        R_ggF_H_ZZ_ATLAS8=(1+(ggF_H_ZZ_TH8-ip_ex_gg_H_ZZ_ATLAS8(mHh))/ip_ex_gg_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_ATLAS8=VBF_H_ZZ_TH8/ip_ex_VBF_H_ZZ_ATLAS8(mHh);
        R_VBF_H_ZZ_ATLAS8=(1+(VBF_H_ZZ_TH8-ip_ex_VBF_H_ZZ_ATLAS8(mHh))/ip_ex_VBF_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_mu_pp_H_VV_CMS8=mu_pp_H_VV_TH8/ip_ex_mu_pp_H_VV_CMS8(mHh);
        R_mu_pp_H_VV_CMS8=(1+(mu_pp_H_VV_TH8-ip_ex_mu_pp_H_VV_CMS8(mHh))/ip_ex_mu_pp_H_VV_CMS8_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_ATLAS8=ggF_H_WW_TH8/ip_ex_gg_H_WW_ATLAS8(mHh);
        R_ggF_H_WW_ATLAS8=(1+(ggF_H_WW_TH8-ip_ex_gg_H_WW_ATLAS8(mHh))/ip_ex_gg_H_WW_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_ATLAS8=VBF_H_WW_TH8/ip_ex_VBF_H_WW_ATLAS8(mHh);
        R_VBF_H_WW_ATLAS8=(1+(VBF_H_WW_TH8-ip_ex_VBF_H_WW_ATLAS8(mHh))/ip_ex_VBF_H_WW_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_ATLAS8=ggF_H_hh_TH8/ip_ex_gg_H_hh_ATLAS8(mHh);
        R_ggF_H_hh_ATLAS8=(1+(ggF_H_hh_TH8-ip_ex_gg_H_hh_ATLAS8(mHh))/ip_ex_gg_H_hh_ATLAS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_CMS8=pp_H_hh_TH8/ip_ex_pp_H_hh_CMS8(mHh);
        R_pp_H_hh_CMS8=(1+(pp_H_hh_TH8-ip_ex_pp_H_hh_CMS8(mHh))/ip_ex_pp_H_hh_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS8=pp_H_hh_bbbb_TH8/ip_ex_pp_phi_hh_bbbb_CMS8(mHh);
        R_pp_H_hh_bbbb_CMS8=(1+(pp_H_hh_bbbb_TH8-ip_ex_pp_phi_hh_bbbb_CMS8(mHh))/ip_ex_pp_phi_hh_bbbb_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_gagabb_CMS8=pp_H_hh_gagabb_TH8/ip_ex_pp_phi_hh_gagabb_CMS8(mHh);
        R_pp_H_hh_gagabb_CMS8=(1+(pp_H_hh_gagabb_TH8-ip_ex_pp_phi_hh_gagabb_CMS8(mHh))/ip_ex_pp_phi_hh_gagabb_CMS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tt_ATLAS8=ggF_H_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mHh);
        R_ggF_H_tt_ATLAS8=(1+(ggF_H_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mHh))/ip_ex_gg_phi_tt_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_bb_CMS8=bbF_H_bb_TH8/ip_ex_bb_phi_bb_CMS8(mHh);
        R_bbF_H_bb_CMS8=(1+(bbF_H_bb_TH8-ip_ex_bb_phi_bb_CMS8(mHh))/ip_ex_bb_phi_bb_CMS8_e(mHh) ) * nftos;
        if(mA>=40.0 && mA<=910.0)
        {
            THoEX_pp_H_AZ_bbll_CMS8=pp_H_AZ_bbll_TH8/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh);
            R_pp_H_AZ_bbll_CMS8=(1+(pp_H_AZ_bbll_TH8-ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh))/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mA>=50.0 && mA<=1000.0)
        {
            THoEX_pp_H_AZ_tautaull_CMS8=pp_H_AZ_tautaull_TH8/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh);
            R_pp_H_AZ_tautaull_CMS8=(1+(pp_H_AZ_tautaull_TH8-ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ttF_H_tt_ATLAS13=ttF_H_tt_TH13/ip_ex_tt_phi_tt_ATLAS13(mHh);
        R_ttF_H_tt_ATLAS13=(1+(ttF_H_tt_TH13-ip_ex_tt_phi_tt_ATLAS13(mHh))/ip_ex_tt_phi_tt_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tt_ATLAS13=bbF_H_tt_TH13/ip_ex_bb_phi_tt_ATLAS13(mHh);
        R_bbF_H_tt_ATLAS13=(1+(bbF_H_tt_TH13-ip_ex_bb_phi_tt_ATLAS13(mHh))/ip_ex_bb_phi_tt_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_bb_CMS13=pp_H_bb_TH13/ip_ex_pp_phi_bb_CMS13(mHh);
        R_pp_H_bb_CMS13=(1+(pp_H_bb_TH13-ip_ex_pp_phi_bb_CMS13(mHh))/ip_ex_pp_phi_bb_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_ATLAS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mHh);
        R_ggF_H_tautau_ATLAS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mHh))/ip_ex_gg_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mHh);
        R_bbF_H_tautau_ATLAS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mHh))/ip_ex_bb_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS13=pp_H_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mHh);
        R_pp_H_gaga_ATLAS13=(1+(pp_H_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mHh))/ip_ex_pp_phi_gaga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS13=ggF_H_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mHh);
        R_ggF_H_gaga_CMS13=(1+(ggF_H_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mHh))/ip_ex_gg_phi_gaga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_llga_ATLAS13=ggF_H_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mHh);
        R_ggF_H_Zga_llga_ATLAS13=(1+(ggF_H_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mHh))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_ATLAS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mHh);
        R_pp_H_Zga_ATLAS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_CMS13(mHh);
        R_pp_H_Zga_llga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_CMS13(mHh))/ip_ex_pp_phi_Zga_llga_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_qqga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_qqga_CMS13(mHh);
        R_pp_H_Zga_qqga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_qqga_CMS13(mHh))/ip_ex_pp_phi_Zga_qqga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_CMS13=ggF_H_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mHh);
        R_ggF_H_Zga_CMS13=(1+(ggF_H_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mHh))/ip_ex_ggF_phi_Zga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_llllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llnunu_CMS13=pp_H_ZZ_TH13/ip_ex_pp_H_ZZ_llnunu_CMS13(mHh);
        R_pp_H_ZZ_llnunu_CMS13=(1+(pp_H_ZZ_TH13-ip_ex_pp_H_ZZ_llnunu_CMS13(mHh))/ip_ex_pp_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llll_ATLAS13=ggF_H_ZZ_llll_TH13/ip_ex_gg_H_ZZ_llll_ATLAS13(mHh);
        R_ggF_H_ZZ_llll_ATLAS13=(1+(ggF_H_ZZ_llll_TH13-ip_ex_gg_H_ZZ_llll_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llll_ATLAS13=VBF_H_ZZ_llll_TH13/ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh);
        R_VBF_H_ZZ_llll_ATLAS13=(1+(VBF_H_ZZ_llll_TH13-ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_qqllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_qqllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_qqllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_qqllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh);
        R_ggF_H_ZZ_llqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llqq_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh);
        R_VBF_H_ZZ_llqq_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_nunuqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh);
        R_ggF_H_ZZ_nunuqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_nunuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llqq_CMS13=pp_H_ZZ_TH13/ip_ex_pp_H_ZZ_llqq_CMS13(mHh);
        R_pp_H_ZZ_llqq_CMS13=(1+(pp_H_ZZ_TH13-ip_ex_pp_H_ZZ_llqq_CMS13(mHh))/ip_ex_pp_H_ZZ_llqq_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_lnuqq_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh);
        R_ggF_H_WW_lnuqq_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_gg_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_lnuqq_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh);
        R_VBF_H_WW_lnuqq_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_VBF_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_enumunu_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_enumunu_ATLAS13(mHh);
        R_ggF_H_WW_enumunu_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_enumunu_ATLAS13(mHh))/ip_ex_gg_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_enumunu_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh);
        R_VBF_H_WW_enumunu_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh))/ip_ex_VBF_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_VBF_H_WW_lnulnu_CMS13=ggF_VBF_H_WW_lnulnu_TH13/ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh);
        R_ggF_VBF_H_WW_lnulnu_CMS13=(1+(ggF_VBF_H_WW_lnulnu_TH13-ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh))/ip_ex_ggVV_H_WW_lnulnu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbgaga_CMS13=pp_H_hh_gagabb_TH13/ip_ex_pp_H_hh_gagabb_CMS13(mHh);
        R_pp_H_hh_bbgaga_CMS13=(1+(pp_H_hh_gagabb_TH13-ip_ex_pp_H_hh_gagabb_CMS13(mHh))/ip_ex_pp_H_hh_gagabb_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_ATLAS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_ATLAS13(mHh);
        R_pp_H_hh_bbbb_ATLAS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_ATLAS13(mHh))/ip_ex_pp_H_hh_bbbb_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbtautau_CMS13=pp_H_hh_bbtautau_TH13/ip_ex_pp_H_hh_bbtautau_CMS13(mHh);
        R_pp_H_hh_bbtautau_CMS13=(1+(pp_H_hh_bbtautau_TH13-ip_ex_pp_H_hh_bbtautau_CMS13(mHh))/ip_ex_pp_H_hh_bbtautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbtautau1_CMS13=pp_H_hh_bbtautau_TH13/ip_ex_pp_H_hh_bbtautau1_CMS13(mHh);
        R_pp_H_hh_bbtautau1_CMS13=(1+(pp_H_hh_bbtautau_TH13-ip_ex_pp_H_hh_bbtautau1_CMS13(mHh))/ip_ex_pp_H_hh_bbtautau1_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bblnulnu_CMS13=pp_H_hh_bblnulnu_TH13/ip_ex_pp_H_hh_bblnulnu_CMS13(mHh);
        R_pp_H_hh_bblnulnu_CMS13=(1+(pp_H_hh_bblnulnu_TH13-ip_ex_pp_H_hh_bblnulnu_CMS13(mHh))/ip_ex_pp_H_hh_bblnulnu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbVV_CMS13=pp_H_hh_bbVV_TH13/ip_ex_pp_H_hh_bbVV_CMS13(mHh);
        R_pp_H_hh_bbVV_CMS13=(1+(pp_H_hh_bbVV_TH13-ip_ex_pp_H_hh_bbVV_CMS13(mHh))/ip_ex_pp_H_hh_bbVV_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=760.0 && mHh<850.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mHh);
        R_ggF_H_tautau_CMS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mHh))/ip_ex_gg_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mHh);
        R_bbF_H_tautau_CMS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mHh))/ip_ex_bb_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS8=ggF_H_gaga_TH8/ip_ex_gg_phi_gaga_CMS8(mHh);
        R_ggF_H_gaga_CMS8=(1+(ggF_H_gaga_TH8-ip_ex_gg_phi_gaga_CMS8(mHh))/ip_ex_gg_phi_gaga_CMS8_e(mHh) ) * nftos;
//    LIMIT_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8(mHh);
//    LIMEST_ggF_H_gaga_CMS8=0.0;
//    DEVIATION_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8(mHh)-ip_ex_gg_phi_gaga_CMS8_e(mHh);
//    BANDSIZE_ggF_H_gaga_CMS8=ip_ex_gg_phi_gaga_CMS8_ep2(mHh)-ip_ex_gg_phi_gaga_CMS8_em2(mHh);
        THoEX_pp_H_Zga_llga_CMS8=pp_H_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mHh);
        R_pp_H_Zga_llga_CMS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mHh))/ip_ex_pp_A_Zga_llga_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_ATLAS8=pp_H_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mHh);
        R_pp_H_Zga_llga_ATLAS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_ATLAS8=ggF_H_ZZ_TH8/ip_ex_gg_H_ZZ_ATLAS8(mHh);
        R_ggF_H_ZZ_ATLAS8=(1+(ggF_H_ZZ_TH8-ip_ex_gg_H_ZZ_ATLAS8(mHh))/ip_ex_gg_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_ATLAS8=VBF_H_ZZ_TH8/ip_ex_VBF_H_ZZ_ATLAS8(mHh);
        R_VBF_H_ZZ_ATLAS8=(1+(VBF_H_ZZ_TH8-ip_ex_VBF_H_ZZ_ATLAS8(mHh))/ip_ex_VBF_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_mu_pp_H_VV_CMS8=mu_pp_H_VV_TH8/ip_ex_mu_pp_H_VV_CMS8(mHh);
        R_mu_pp_H_VV_CMS8=(1+(mu_pp_H_VV_TH8-ip_ex_mu_pp_H_VV_CMS8(mHh))/ip_ex_mu_pp_H_VV_CMS8_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_ATLAS8=ggF_H_WW_TH8/ip_ex_gg_H_WW_ATLAS8(mHh);
        R_ggF_H_WW_ATLAS8=(1+(ggF_H_WW_TH8-ip_ex_gg_H_WW_ATLAS8(mHh))/ip_ex_gg_H_WW_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_ATLAS8=VBF_H_WW_TH8/ip_ex_VBF_H_WW_ATLAS8(mHh);
        R_VBF_H_WW_ATLAS8=(1+(VBF_H_WW_TH8-ip_ex_VBF_H_WW_ATLAS8(mHh))/ip_ex_VBF_H_WW_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_ATLAS8=ggF_H_hh_TH8/ip_ex_gg_H_hh_ATLAS8(mHh);
        R_ggF_H_hh_ATLAS8=(1+(ggF_H_hh_TH8-ip_ex_gg_H_hh_ATLAS8(mHh))/ip_ex_gg_H_hh_ATLAS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_CMS8=pp_H_hh_TH8/ip_ex_pp_H_hh_CMS8(mHh);
        R_pp_H_hh_CMS8=(1+(pp_H_hh_TH8-ip_ex_pp_H_hh_CMS8(mHh))/ip_ex_pp_H_hh_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS8=pp_H_hh_bbbb_TH8/ip_ex_pp_phi_hh_bbbb_CMS8(mHh);
        R_pp_H_hh_bbbb_CMS8=(1+(pp_H_hh_bbbb_TH8-ip_ex_pp_phi_hh_bbbb_CMS8(mHh))/ip_ex_pp_phi_hh_bbbb_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_gagabb_CMS8=pp_H_hh_gagabb_TH8/ip_ex_pp_phi_hh_gagabb_CMS8(mHh);
        R_pp_H_hh_gagabb_CMS8=(1+(pp_H_hh_gagabb_TH8-ip_ex_pp_phi_hh_gagabb_CMS8(mHh))/ip_ex_pp_phi_hh_gagabb_CMS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tt_ATLAS8=ggF_H_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mHh);
        R_ggF_H_tt_ATLAS8=(1+(ggF_H_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mHh))/ip_ex_gg_phi_tt_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_bb_CMS8=bbF_H_bb_TH8/ip_ex_bb_phi_bb_CMS8(mHh);
        R_bbF_H_bb_CMS8=(1+(bbF_H_bb_TH8-ip_ex_bb_phi_bb_CMS8(mHh))/ip_ex_bb_phi_bb_CMS8_e(mHh) ) * nftos;
        if(mA>=40.0 && mA<=910.0)
        {
            THoEX_pp_H_AZ_bbll_CMS8=pp_H_AZ_bbll_TH8/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh);
            R_pp_H_AZ_bbll_CMS8=(1+(pp_H_AZ_bbll_TH8-ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh))/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mA>=50.0 && mA<=1000.0)
        {
            THoEX_pp_H_AZ_tautaull_CMS8=pp_H_AZ_tautaull_TH8/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh);
            R_pp_H_AZ_tautaull_CMS8=(1+(pp_H_AZ_tautaull_TH8-ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ttF_H_tt_ATLAS13=ttF_H_tt_TH13/ip_ex_tt_phi_tt_ATLAS13(mHh);
        R_ttF_H_tt_ATLAS13=(1+(ttF_H_tt_TH13-ip_ex_tt_phi_tt_ATLAS13(mHh))/ip_ex_tt_phi_tt_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tt_ATLAS13=bbF_H_tt_TH13/ip_ex_bb_phi_tt_ATLAS13(mHh);
        R_bbF_H_tt_ATLAS13=(1+(bbF_H_tt_TH13-ip_ex_bb_phi_tt_ATLAS13(mHh))/ip_ex_bb_phi_tt_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_bb_CMS13=pp_H_bb_TH13/ip_ex_pp_phi_bb_CMS13(mHh);
        R_pp_H_bb_CMS13=(1+(pp_H_bb_TH13-ip_ex_pp_phi_bb_CMS13(mHh))/ip_ex_pp_phi_bb_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_ATLAS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mHh);
        R_ggF_H_tautau_ATLAS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mHh))/ip_ex_gg_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mHh);
        R_bbF_H_tautau_ATLAS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mHh))/ip_ex_bb_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS13=pp_H_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mHh);
        R_pp_H_gaga_ATLAS13=(1+(pp_H_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mHh))/ip_ex_pp_phi_gaga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS13=ggF_H_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mHh);
        R_ggF_H_gaga_CMS13=(1+(ggF_H_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mHh))/ip_ex_gg_phi_gaga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_llga_ATLAS13=ggF_H_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mHh);
        R_ggF_H_Zga_llga_ATLAS13=(1+(ggF_H_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mHh))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_ATLAS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mHh);
        R_pp_H_Zga_ATLAS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_CMS13(mHh);
        R_pp_H_Zga_llga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_CMS13(mHh))/ip_ex_pp_phi_Zga_llga_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_qqga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_qqga_CMS13(mHh);
        R_pp_H_Zga_qqga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_qqga_CMS13(mHh))/ip_ex_pp_phi_Zga_qqga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_CMS13=ggF_H_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mHh);
        R_ggF_H_Zga_CMS13=(1+(ggF_H_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mHh))/ip_ex_ggF_phi_Zga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_llllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llnunu_CMS13=pp_H_ZZ_TH13/ip_ex_pp_H_ZZ_llnunu_CMS13(mHh);
        R_pp_H_ZZ_llnunu_CMS13=(1+(pp_H_ZZ_TH13-ip_ex_pp_H_ZZ_llnunu_CMS13(mHh))/ip_ex_pp_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llll_ATLAS13=ggF_H_ZZ_llll_TH13/ip_ex_gg_H_ZZ_llll_ATLAS13(mHh);
        R_ggF_H_ZZ_llll_ATLAS13=(1+(ggF_H_ZZ_llll_TH13-ip_ex_gg_H_ZZ_llll_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llll_ATLAS13=VBF_H_ZZ_llll_TH13/ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh);
        R_VBF_H_ZZ_llll_ATLAS13=(1+(VBF_H_ZZ_llll_TH13-ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_qqllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_qqllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_qqllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_qqllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh);
        R_ggF_H_ZZ_llqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llqq_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh);
        R_VBF_H_ZZ_llqq_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_nunuqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh);
        R_ggF_H_ZZ_nunuqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_nunuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llqq_CMS13=pp_H_ZZ_TH13/ip_ex_pp_H_ZZ_llqq_CMS13(mHh);
        R_pp_H_ZZ_llqq_CMS13=(1+(pp_H_ZZ_TH13-ip_ex_pp_H_ZZ_llqq_CMS13(mHh))/ip_ex_pp_H_ZZ_llqq_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_lnuqq_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh);
        R_ggF_H_WW_lnuqq_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_gg_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_lnuqq_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh);
        R_VBF_H_WW_lnuqq_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_VBF_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_enumunu_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_enumunu_ATLAS13(mHh);
        R_ggF_H_WW_enumunu_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_enumunu_ATLAS13(mHh))/ip_ex_gg_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_enumunu_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh);
        R_VBF_H_WW_enumunu_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh))/ip_ex_VBF_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_VBF_H_WW_lnulnu_CMS13=ggF_VBF_H_WW_lnulnu_TH13/ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh);
        R_ggF_VBF_H_WW_lnulnu_CMS13=(1+(ggF_VBF_H_WW_lnulnu_TH13-ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh))/ip_ex_ggVV_H_WW_lnulnu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbgaga_CMS13=pp_H_hh_gagabb_TH13/ip_ex_pp_H_hh_gagabb_CMS13(mHh);
        R_pp_H_hh_bbgaga_CMS13=(1+(pp_H_hh_gagabb_TH13-ip_ex_pp_H_hh_gagabb_CMS13(mHh))/ip_ex_pp_H_hh_gagabb_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_ATLAS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_ATLAS13(mHh);
        R_pp_H_hh_bbbb_ATLAS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_ATLAS13(mHh))/ip_ex_pp_H_hh_bbbb_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbtautau_CMS13=pp_H_hh_bbtautau_TH13/ip_ex_pp_H_hh_bbtautau_CMS13(mHh);
        R_pp_H_hh_bbtautau_CMS13=(1+(pp_H_hh_bbtautau_TH13-ip_ex_pp_H_hh_bbtautau_CMS13(mHh))/ip_ex_pp_H_hh_bbtautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbtautau1_CMS13=pp_H_hh_bbtautau_TH13/ip_ex_pp_H_hh_bbtautau1_CMS13(mHh);
        R_pp_H_hh_bbtautau1_CMS13=(1+(pp_H_hh_bbtautau_TH13-ip_ex_pp_H_hh_bbtautau1_CMS13(mHh))/ip_ex_pp_H_hh_bbtautau1_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bblnulnu_CMS13=pp_H_hh_bblnulnu_TH13/ip_ex_pp_H_hh_bblnulnu_CMS13(mHh);
        R_pp_H_hh_bblnulnu_CMS13=(1+(pp_H_hh_bblnulnu_TH13-ip_ex_pp_H_hh_bblnulnu_CMS13(mHh))/ip_ex_pp_H_hh_bblnulnu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbVV_CMS13=pp_H_hh_bbVV_TH13/ip_ex_pp_H_hh_bbVV_CMS13(mHh);
        R_pp_H_hh_bbVV_CMS13=(1+(pp_H_hh_bbVV_TH13-ip_ex_pp_H_hh_bbVV_CMS13(mHh))/ip_ex_pp_H_hh_bbVV_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=850.0 && mHh<900.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mHh);
        R_ggF_H_tautau_CMS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mHh))/ip_ex_gg_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mHh);
        R_bbF_H_tautau_CMS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mHh))/ip_ex_bb_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_CMS8=pp_H_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mHh);
        R_pp_H_Zga_llga_CMS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mHh))/ip_ex_pp_A_Zga_llga_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_ATLAS8=pp_H_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mHh);
        R_pp_H_Zga_llga_ATLAS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_ATLAS8=ggF_H_ZZ_TH8/ip_ex_gg_H_ZZ_ATLAS8(mHh);
        R_ggF_H_ZZ_ATLAS8=(1+(ggF_H_ZZ_TH8-ip_ex_gg_H_ZZ_ATLAS8(mHh))/ip_ex_gg_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_ATLAS8=VBF_H_ZZ_TH8/ip_ex_VBF_H_ZZ_ATLAS8(mHh);
        R_VBF_H_ZZ_ATLAS8=(1+(VBF_H_ZZ_TH8-ip_ex_VBF_H_ZZ_ATLAS8(mHh))/ip_ex_VBF_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_mu_pp_H_VV_CMS8=mu_pp_H_VV_TH8/ip_ex_mu_pp_H_VV_CMS8(mHh);
        R_mu_pp_H_VV_CMS8=(1+(mu_pp_H_VV_TH8-ip_ex_mu_pp_H_VV_CMS8(mHh))/ip_ex_mu_pp_H_VV_CMS8_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_ATLAS8=ggF_H_WW_TH8/ip_ex_gg_H_WW_ATLAS8(mHh);
        R_ggF_H_WW_ATLAS8=(1+(ggF_H_WW_TH8-ip_ex_gg_H_WW_ATLAS8(mHh))/ip_ex_gg_H_WW_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_ATLAS8=VBF_H_WW_TH8/ip_ex_VBF_H_WW_ATLAS8(mHh);
        R_VBF_H_WW_ATLAS8=(1+(VBF_H_WW_TH8-ip_ex_VBF_H_WW_ATLAS8(mHh))/ip_ex_VBF_H_WW_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_ATLAS8=ggF_H_hh_TH8/ip_ex_gg_H_hh_ATLAS8(mHh);
        R_ggF_H_hh_ATLAS8=(1+(ggF_H_hh_TH8-ip_ex_gg_H_hh_ATLAS8(mHh))/ip_ex_gg_H_hh_ATLAS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_CMS8=pp_H_hh_TH8/ip_ex_pp_H_hh_CMS8(mHh);
        R_pp_H_hh_CMS8=(1+(pp_H_hh_TH8-ip_ex_pp_H_hh_CMS8(mHh))/ip_ex_pp_H_hh_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS8=pp_H_hh_bbbb_TH8/ip_ex_pp_phi_hh_bbbb_CMS8(mHh);
        R_pp_H_hh_bbbb_CMS8=(1+(pp_H_hh_bbbb_TH8-ip_ex_pp_phi_hh_bbbb_CMS8(mHh))/ip_ex_pp_phi_hh_bbbb_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_gagabb_CMS8=pp_H_hh_gagabb_TH8/ip_ex_pp_phi_hh_gagabb_CMS8(mHh);
        R_pp_H_hh_gagabb_CMS8=(1+(pp_H_hh_gagabb_TH8-ip_ex_pp_phi_hh_gagabb_CMS8(mHh))/ip_ex_pp_phi_hh_gagabb_CMS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tt_ATLAS8=ggF_H_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mHh);
        R_ggF_H_tt_ATLAS8=(1+(ggF_H_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mHh))/ip_ex_gg_phi_tt_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_bb_CMS8=bbF_H_bb_TH8/ip_ex_bb_phi_bb_CMS8(mHh);
        R_bbF_H_bb_CMS8=(1+(bbF_H_bb_TH8-ip_ex_bb_phi_bb_CMS8(mHh))/ip_ex_bb_phi_bb_CMS8_e(mHh) ) * nftos;
        if(mA>=40.0 && mA<=910.0)
        {
            THoEX_pp_H_AZ_bbll_CMS8=pp_H_AZ_bbll_TH8/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh);
            R_pp_H_AZ_bbll_CMS8=(1+(pp_H_AZ_bbll_TH8-ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh))/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mA>=50.0 && mA<=1000.0)
        {
            THoEX_pp_H_AZ_tautaull_CMS8=pp_H_AZ_tautaull_TH8/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh);
            R_pp_H_AZ_tautaull_CMS8=(1+(pp_H_AZ_tautaull_TH8-ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ttF_H_tt_ATLAS13=ttF_H_tt_TH13/ip_ex_tt_phi_tt_ATLAS13(mHh);
        R_ttF_H_tt_ATLAS13=(1+(ttF_H_tt_TH13-ip_ex_tt_phi_tt_ATLAS13(mHh))/ip_ex_tt_phi_tt_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tt_ATLAS13=bbF_H_tt_TH13/ip_ex_bb_phi_tt_ATLAS13(mHh);
        R_bbF_H_tt_ATLAS13=(1+(bbF_H_tt_TH13-ip_ex_bb_phi_tt_ATLAS13(mHh))/ip_ex_bb_phi_tt_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_bb_CMS13=pp_H_bb_TH13/ip_ex_pp_phi_bb_CMS13(mHh);
        R_pp_H_bb_CMS13=(1+(pp_H_bb_TH13-ip_ex_pp_phi_bb_CMS13(mHh))/ip_ex_pp_phi_bb_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_ATLAS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mHh);
        R_ggF_H_tautau_ATLAS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mHh))/ip_ex_gg_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mHh);
        R_bbF_H_tautau_ATLAS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mHh))/ip_ex_bb_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS13=pp_H_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mHh);
        R_pp_H_gaga_ATLAS13=(1+(pp_H_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mHh))/ip_ex_pp_phi_gaga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS13=ggF_H_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mHh);
        R_ggF_H_gaga_CMS13=(1+(ggF_H_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mHh))/ip_ex_gg_phi_gaga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_llga_ATLAS13=ggF_H_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mHh);
        R_ggF_H_Zga_llga_ATLAS13=(1+(ggF_H_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mHh))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_ATLAS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mHh);
        R_pp_H_Zga_ATLAS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_CMS13(mHh);
        R_pp_H_Zga_llga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_CMS13(mHh))/ip_ex_pp_phi_Zga_llga_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_qqga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_qqga_CMS13(mHh);
        R_pp_H_Zga_qqga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_qqga_CMS13(mHh))/ip_ex_pp_phi_Zga_qqga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_CMS13=ggF_H_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mHh);
        R_ggF_H_Zga_CMS13=(1+(ggF_H_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mHh))/ip_ex_ggF_phi_Zga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_llllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llnunu_CMS13=pp_H_ZZ_TH13/ip_ex_pp_H_ZZ_llnunu_CMS13(mHh);
        R_pp_H_ZZ_llnunu_CMS13=(1+(pp_H_ZZ_TH13-ip_ex_pp_H_ZZ_llnunu_CMS13(mHh))/ip_ex_pp_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llll_ATLAS13=ggF_H_ZZ_llll_TH13/ip_ex_gg_H_ZZ_llll_ATLAS13(mHh);
        R_ggF_H_ZZ_llll_ATLAS13=(1+(ggF_H_ZZ_llll_TH13-ip_ex_gg_H_ZZ_llll_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llll_ATLAS13=VBF_H_ZZ_llll_TH13/ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh);
        R_VBF_H_ZZ_llll_ATLAS13=(1+(VBF_H_ZZ_llll_TH13-ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_qqllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_qqllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_qqllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_qqllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh);
        R_ggF_H_ZZ_llqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llqq_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh);
        R_VBF_H_ZZ_llqq_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_nunuqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh);
        R_ggF_H_ZZ_nunuqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_nunuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llqq_CMS13=pp_H_ZZ_TH13/ip_ex_pp_H_ZZ_llqq_CMS13(mHh);
        R_pp_H_ZZ_llqq_CMS13=(1+(pp_H_ZZ_TH13-ip_ex_pp_H_ZZ_llqq_CMS13(mHh))/ip_ex_pp_H_ZZ_llqq_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_lnuqq_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh);
        R_ggF_H_WW_lnuqq_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_gg_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_lnuqq_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh);
        R_VBF_H_WW_lnuqq_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_VBF_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_enumunu_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_enumunu_ATLAS13(mHh);
        R_ggF_H_WW_enumunu_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_enumunu_ATLAS13(mHh))/ip_ex_gg_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_enumunu_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh);
        R_VBF_H_WW_enumunu_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh))/ip_ex_VBF_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_VBF_H_WW_lnulnu_CMS13=ggF_VBF_H_WW_lnulnu_TH13/ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh);
        R_ggF_VBF_H_WW_lnulnu_CMS13=(1+(ggF_VBF_H_WW_lnulnu_TH13-ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh))/ip_ex_ggVV_H_WW_lnulnu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbgaga_CMS13=pp_H_hh_gagabb_TH13/ip_ex_pp_H_hh_gagabb_CMS13(mHh);
        R_pp_H_hh_bbgaga_CMS13=(1+(pp_H_hh_gagabb_TH13-ip_ex_pp_H_hh_gagabb_CMS13(mHh))/ip_ex_pp_H_hh_gagabb_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_ATLAS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_ATLAS13(mHh);
        R_pp_H_hh_bbbb_ATLAS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_ATLAS13(mHh))/ip_ex_pp_H_hh_bbbb_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbtautau_CMS13=pp_H_hh_bbtautau_TH13/ip_ex_pp_H_hh_bbtautau_CMS13(mHh);
        R_pp_H_hh_bbtautau_CMS13=(1+(pp_H_hh_bbtautau_TH13-ip_ex_pp_H_hh_bbtautau_CMS13(mHh))/ip_ex_pp_H_hh_bbtautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbtautau1_CMS13=pp_H_hh_bbtautau_TH13/ip_ex_pp_H_hh_bbtautau1_CMS13(mHh);
        R_pp_H_hh_bbtautau1_CMS13=(1+(pp_H_hh_bbtautau_TH13-ip_ex_pp_H_hh_bbtautau1_CMS13(mHh))/ip_ex_pp_H_hh_bbtautau1_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bblnulnu_CMS13=pp_H_hh_bblnulnu_TH13/ip_ex_pp_H_hh_bblnulnu_CMS13(mHh);
        R_pp_H_hh_bblnulnu_CMS13=(1+(pp_H_hh_bblnulnu_TH13-ip_ex_pp_H_hh_bblnulnu_CMS13(mHh))/ip_ex_pp_H_hh_bblnulnu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbVV_CMS13=pp_H_hh_bbVV_TH13/ip_ex_pp_H_hh_bbVV_CMS13(mHh);
        R_pp_H_hh_bbVV_CMS13=(1+(pp_H_hh_bbVV_TH13-ip_ex_pp_H_hh_bbVV_CMS13(mHh))/ip_ex_pp_H_hh_bbVV_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=900.0 && mHh<1000.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mHh);
        R_ggF_H_tautau_CMS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mHh))/ip_ex_gg_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mHh);
        R_bbF_H_tautau_CMS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mHh))/ip_ex_bb_phi_tautau_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_CMS8=pp_H_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mHh);
        R_pp_H_Zga_llga_CMS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mHh))/ip_ex_pp_A_Zga_llga_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_ATLAS8=pp_H_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mHh);
        R_pp_H_Zga_llga_ATLAS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_ATLAS8=ggF_H_ZZ_TH8/ip_ex_gg_H_ZZ_ATLAS8(mHh);
        R_ggF_H_ZZ_ATLAS8=(1+(ggF_H_ZZ_TH8-ip_ex_gg_H_ZZ_ATLAS8(mHh))/ip_ex_gg_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_ATLAS8=VBF_H_ZZ_TH8/ip_ex_VBF_H_ZZ_ATLAS8(mHh);
        R_VBF_H_ZZ_ATLAS8=(1+(VBF_H_ZZ_TH8-ip_ex_VBF_H_ZZ_ATLAS8(mHh))/ip_ex_VBF_H_ZZ_ATLAS8_e(mHh) ) * nftos;
        THoEX_mu_pp_H_VV_CMS8=mu_pp_H_VV_TH8/ip_ex_mu_pp_H_VV_CMS8(mHh);
        R_mu_pp_H_VV_CMS8=(1+(mu_pp_H_VV_TH8-ip_ex_mu_pp_H_VV_CMS8(mHh))/ip_ex_mu_pp_H_VV_CMS8_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_ATLAS8=ggF_H_WW_TH8/ip_ex_gg_H_WW_ATLAS8(mHh);
        R_ggF_H_WW_ATLAS8=(1+(ggF_H_WW_TH8-ip_ex_gg_H_WW_ATLAS8(mHh))/ip_ex_gg_H_WW_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_ATLAS8=VBF_H_WW_TH8/ip_ex_VBF_H_WW_ATLAS8(mHh);
        R_VBF_H_WW_ATLAS8=(1+(VBF_H_WW_TH8-ip_ex_VBF_H_WW_ATLAS8(mHh))/ip_ex_VBF_H_WW_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_ATLAS8=ggF_H_hh_TH8/ip_ex_gg_H_hh_ATLAS8(mHh);
        R_ggF_H_hh_ATLAS8=(1+(ggF_H_hh_TH8-ip_ex_gg_H_hh_ATLAS8(mHh))/ip_ex_gg_H_hh_ATLAS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_CMS8=pp_H_hh_TH8/ip_ex_pp_H_hh_CMS8(mHh);
        R_pp_H_hh_CMS8=(1+(pp_H_hh_TH8-ip_ex_pp_H_hh_CMS8(mHh))/ip_ex_pp_H_hh_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS8=pp_H_hh_bbbb_TH8/ip_ex_pp_phi_hh_bbbb_CMS8(mHh);
        R_pp_H_hh_bbbb_CMS8=(1+(pp_H_hh_bbbb_TH8-ip_ex_pp_phi_hh_bbbb_CMS8(mHh))/ip_ex_pp_phi_hh_bbbb_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_gagabb_CMS8=pp_H_hh_gagabb_TH8/ip_ex_pp_phi_hh_gagabb_CMS8(mHh);
        R_pp_H_hh_gagabb_CMS8=(1+(pp_H_hh_gagabb_TH8-ip_ex_pp_phi_hh_gagabb_CMS8(mHh))/ip_ex_pp_phi_hh_gagabb_CMS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tt_ATLAS8=ggF_H_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mHh);
        R_ggF_H_tt_ATLAS8=(1+(ggF_H_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mHh))/ip_ex_gg_phi_tt_ATLAS8_e(mHh) ) * nftos;
        if(mA>=40.0 && mA<=910.0)
        {
            THoEX_pp_H_AZ_bbll_CMS8=pp_H_AZ_bbll_TH8/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh);
            R_pp_H_AZ_bbll_CMS8=(1+(pp_H_AZ_bbll_TH8-ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh))/ip_ex_pp_H_AZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mA>=50.0 && mA<=1000.0)
        {
            THoEX_pp_H_AZ_tautaull_CMS8=pp_H_AZ_tautaull_TH8/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh);
            R_pp_H_AZ_tautaull_CMS8=(1+(pp_H_AZ_tautaull_TH8-ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_H_AZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ttF_H_tt_ATLAS13=ttF_H_tt_TH13/ip_ex_tt_phi_tt_ATLAS13(mHh);
        R_ttF_H_tt_ATLAS13=(1+(ttF_H_tt_TH13-ip_ex_tt_phi_tt_ATLAS13(mHh))/ip_ex_tt_phi_tt_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tt_ATLAS13=bbF_H_tt_TH13/ip_ex_bb_phi_tt_ATLAS13(mHh);
        R_bbF_H_tt_ATLAS13=(1+(bbF_H_tt_TH13-ip_ex_bb_phi_tt_ATLAS13(mHh))/ip_ex_bb_phi_tt_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_bb_CMS13=pp_H_bb_TH13/ip_ex_pp_phi_bb_CMS13(mHh);
        R_pp_H_bb_CMS13=(1+(pp_H_bb_TH13-ip_ex_pp_phi_bb_CMS13(mHh))/ip_ex_pp_phi_bb_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_ATLAS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mHh);
        R_ggF_H_tautau_ATLAS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mHh))/ip_ex_gg_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mHh);
        R_bbF_H_tautau_ATLAS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mHh))/ip_ex_bb_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS13=pp_H_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mHh);
        R_pp_H_gaga_ATLAS13=(1+(pp_H_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mHh))/ip_ex_pp_phi_gaga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS13=ggF_H_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mHh);
        R_ggF_H_gaga_CMS13=(1+(ggF_H_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mHh))/ip_ex_gg_phi_gaga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_llga_ATLAS13=ggF_H_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mHh);
        R_ggF_H_Zga_llga_ATLAS13=(1+(ggF_H_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mHh))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_ATLAS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mHh);
        R_pp_H_Zga_ATLAS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_CMS13(mHh);
        R_pp_H_Zga_llga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_CMS13(mHh))/ip_ex_pp_phi_Zga_llga_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_qqga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_qqga_CMS13(mHh);
        R_pp_H_Zga_qqga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_qqga_CMS13(mHh))/ip_ex_pp_phi_Zga_qqga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_CMS13=ggF_H_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mHh);
        R_ggF_H_Zga_CMS13=(1+(ggF_H_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mHh))/ip_ex_ggF_phi_Zga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_llllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llnunu_CMS13=pp_H_ZZ_TH13/ip_ex_pp_H_ZZ_llnunu_CMS13(mHh);
        R_pp_H_ZZ_llnunu_CMS13=(1+(pp_H_ZZ_TH13-ip_ex_pp_H_ZZ_llnunu_CMS13(mHh))/ip_ex_pp_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llll_ATLAS13=ggF_H_ZZ_llll_TH13/ip_ex_gg_H_ZZ_llll_ATLAS13(mHh);
        R_ggF_H_ZZ_llll_ATLAS13=(1+(ggF_H_ZZ_llll_TH13-ip_ex_gg_H_ZZ_llll_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llll_ATLAS13=VBF_H_ZZ_llll_TH13/ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh);
        R_VBF_H_ZZ_llll_ATLAS13=(1+(VBF_H_ZZ_llll_TH13-ip_ex_VBF_H_ZZ_llll_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llll_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_qqllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_qqllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_qqllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_qqllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh);
        R_ggF_H_ZZ_llqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llqq_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh);
        R_VBF_H_ZZ_llqq_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_nunuqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh);
        R_ggF_H_ZZ_nunuqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_nunuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llqq_CMS13=pp_H_ZZ_TH13/ip_ex_pp_H_ZZ_llqq_CMS13(mHh);
        R_pp_H_ZZ_llqq_CMS13=(1+(pp_H_ZZ_TH13-ip_ex_pp_H_ZZ_llqq_CMS13(mHh))/ip_ex_pp_H_ZZ_llqq_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_lnuqq_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh);
        R_ggF_H_WW_lnuqq_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_gg_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_lnuqq_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh);
        R_VBF_H_WW_lnuqq_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_VBF_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_enumunu_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_enumunu_ATLAS13(mHh);
        R_ggF_H_WW_enumunu_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_enumunu_ATLAS13(mHh))/ip_ex_gg_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_enumunu_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh);
        R_VBF_H_WW_enumunu_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh))/ip_ex_VBF_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_VBF_H_WW_lnulnu_CMS13=ggF_VBF_H_WW_lnulnu_TH13/ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh);
        R_ggF_VBF_H_WW_lnulnu_CMS13=(1+(ggF_VBF_H_WW_lnulnu_TH13-ip_ex_ggVV_H_WW_lnulnu_CMS13(mHh))/ip_ex_ggVV_H_WW_lnulnu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_ATLAS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_ATLAS13(mHh);
        R_pp_H_hh_bbbb_ATLAS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_ATLAS13(mHh))/ip_ex_pp_H_hh_bbbb_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=1000.0 && mHh<1100.0)
    {
        THoEX_pp_H_Zga_llga_CMS8=pp_H_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mHh);
        R_pp_H_Zga_llga_CMS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mHh))/ip_ex_pp_A_Zga_llga_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_ATLAS8=pp_H_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mHh);
        R_pp_H_Zga_llga_ATLAS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_ATLAS8=ggF_H_WW_TH8/ip_ex_gg_H_WW_ATLAS8(mHh);
        R_ggF_H_WW_ATLAS8=(1+(ggF_H_WW_TH8-ip_ex_gg_H_WW_ATLAS8(mHh))/ip_ex_gg_H_WW_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_ATLAS8=VBF_H_WW_TH8/ip_ex_VBF_H_WW_ATLAS8(mHh);
        R_VBF_H_WW_ATLAS8=(1+(VBF_H_WW_TH8-ip_ex_VBF_H_WW_ATLAS8(mHh))/ip_ex_VBF_H_WW_ATLAS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS8=pp_H_hh_bbbb_TH8/ip_ex_pp_phi_hh_bbbb_CMS8(mHh);
        R_pp_H_hh_bbbb_CMS8=(1+(pp_H_hh_bbbb_TH8-ip_ex_pp_phi_hh_bbbb_CMS8(mHh))/ip_ex_pp_phi_hh_bbbb_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_gagabb_CMS8=pp_H_hh_gagabb_TH8/ip_ex_pp_phi_hh_gagabb_CMS8(mHh);
        R_pp_H_hh_gagabb_CMS8=(1+(pp_H_hh_gagabb_TH8-ip_ex_pp_phi_hh_gagabb_CMS8(mHh))/ip_ex_pp_phi_hh_gagabb_CMS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tt_ATLAS8=ggF_H_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mHh);
        R_ggF_H_tt_ATLAS8=(1+(ggF_H_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mHh))/ip_ex_gg_phi_tt_ATLAS8_e(mHh) ) * nftos;

        THoEX_pp_H_bb_CMS13=pp_H_bb_TH13/ip_ex_pp_phi_bb_CMS13(mHh);
        R_pp_H_bb_CMS13=(1+(pp_H_bb_TH13-ip_ex_pp_phi_bb_CMS13(mHh))/ip_ex_pp_phi_bb_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_ATLAS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mHh);
        R_ggF_H_tautau_ATLAS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mHh))/ip_ex_gg_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mHh);
        R_bbF_H_tautau_ATLAS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mHh))/ip_ex_bb_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS13=pp_H_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mHh);
        R_pp_H_gaga_ATLAS13=(1+(pp_H_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mHh))/ip_ex_pp_phi_gaga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS13=ggF_H_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mHh);
        R_ggF_H_gaga_CMS13=(1+(ggF_H_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mHh))/ip_ex_gg_phi_gaga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_llga_ATLAS13=ggF_H_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mHh);
        R_ggF_H_Zga_llga_ATLAS13=(1+(ggF_H_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mHh))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_ATLAS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mHh);
        R_pp_H_Zga_ATLAS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_CMS13(mHh);
        R_pp_H_Zga_llga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_CMS13(mHh))/ip_ex_pp_phi_Zga_llga_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_qqga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_qqga_CMS13(mHh);
        R_pp_H_Zga_qqga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_qqga_CMS13(mHh))/ip_ex_pp_phi_Zga_qqga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_CMS13=ggF_H_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mHh);
        R_ggF_H_Zga_CMS13=(1+(ggF_H_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mHh))/ip_ex_ggF_phi_Zga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_llllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llnunu_CMS13=pp_H_ZZ_TH13/ip_ex_pp_H_ZZ_llnunu_CMS13(mHh);
        R_pp_H_ZZ_llnunu_CMS13=(1+(pp_H_ZZ_TH13-ip_ex_pp_H_ZZ_llnunu_CMS13(mHh))/ip_ex_pp_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_qqllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_qqllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_qqllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_qqllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh);
        R_ggF_H_ZZ_llqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llqq_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh);
        R_VBF_H_ZZ_llqq_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_nunuqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh);
        R_ggF_H_ZZ_nunuqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_nunuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llqq_CMS13=pp_H_ZZ_TH13/ip_ex_pp_H_ZZ_llqq_CMS13(mHh);
        R_pp_H_ZZ_llqq_CMS13=(1+(pp_H_ZZ_TH13-ip_ex_pp_H_ZZ_llqq_CMS13(mHh))/ip_ex_pp_H_ZZ_llqq_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_lnuqq_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh);
        R_ggF_H_WW_lnuqq_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_gg_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_lnuqq_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh);
        R_VBF_H_WW_lnuqq_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_VBF_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_enumunu_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_enumunu_ATLAS13(mHh);
        R_ggF_H_WW_enumunu_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_enumunu_ATLAS13(mHh))/ip_ex_gg_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_enumunu_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh);
        R_VBF_H_WW_enumunu_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh))/ip_ex_VBF_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_ATLAS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_ATLAS13(mHh);
        R_pp_H_hh_bbbb_ATLAS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_ATLAS13(mHh))/ip_ex_pp_H_hh_bbbb_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=1100.0 && mHh<1200.0)
    {
        THoEX_pp_H_Zga_llga_CMS8=pp_H_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mHh);
        R_pp_H_Zga_llga_CMS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mHh))/ip_ex_pp_A_Zga_llga_CMS8_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_ATLAS8=pp_H_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mHh);
        R_pp_H_Zga_llga_ATLAS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_ATLAS8=ggF_H_WW_TH8/ip_ex_gg_H_WW_ATLAS8(mHh);
        R_ggF_H_WW_ATLAS8=(1+(ggF_H_WW_TH8-ip_ex_gg_H_WW_ATLAS8(mHh))/ip_ex_gg_H_WW_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_ATLAS8=VBF_H_WW_TH8/ip_ex_VBF_H_WW_ATLAS8(mHh);
        R_VBF_H_WW_ATLAS8=(1+(VBF_H_WW_TH8-ip_ex_VBF_H_WW_ATLAS8(mHh))/ip_ex_VBF_H_WW_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tt_ATLAS8=ggF_H_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mHh);
        R_ggF_H_tt_ATLAS8=(1+(ggF_H_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mHh))/ip_ex_gg_phi_tt_ATLAS8_e(mHh) ) * nftos;

        THoEX_pp_H_bb_CMS13=pp_H_bb_TH13/ip_ex_pp_phi_bb_CMS13(mHh);
        R_pp_H_bb_CMS13=(1+(pp_H_bb_TH13-ip_ex_pp_phi_bb_CMS13(mHh))/ip_ex_pp_phi_bb_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_ATLAS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mHh);
        R_ggF_H_tautau_ATLAS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mHh))/ip_ex_gg_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mHh);
        R_bbF_H_tautau_ATLAS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mHh))/ip_ex_bb_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS13=pp_H_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mHh);
        R_pp_H_gaga_ATLAS13=(1+(pp_H_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mHh))/ip_ex_pp_phi_gaga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS13=ggF_H_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mHh);
        R_ggF_H_gaga_CMS13=(1+(ggF_H_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mHh))/ip_ex_gg_phi_gaga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_llga_ATLAS13=ggF_H_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mHh);
        R_ggF_H_Zga_llga_ATLAS13=(1+(ggF_H_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mHh))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_ATLAS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mHh);
        R_pp_H_Zga_ATLAS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_CMS13(mHh);
        R_pp_H_Zga_llga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_CMS13(mHh))/ip_ex_pp_phi_Zga_llga_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_qqga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_qqga_CMS13(mHh);
        R_pp_H_Zga_qqga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_qqga_CMS13(mHh))/ip_ex_pp_phi_Zga_qqga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_CMS13=ggF_H_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mHh);
        R_ggF_H_Zga_CMS13=(1+(ggF_H_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mHh))/ip_ex_ggF_phi_Zga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_llllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_llllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llnunu_CMS13=pp_H_ZZ_TH13/ip_ex_pp_H_ZZ_llnunu_CMS13(mHh);
        R_pp_H_ZZ_llnunu_CMS13=(1+(pp_H_ZZ_TH13-ip_ex_pp_H_ZZ_llnunu_CMS13(mHh))/ip_ex_pp_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_qqllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_qqllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_qqllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_qqllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh);
        R_ggF_H_ZZ_llqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llqq_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh);
        R_VBF_H_ZZ_llqq_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_nunuqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh);
        R_ggF_H_ZZ_nunuqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_nunuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llqq_CMS13=pp_H_ZZ_TH13/ip_ex_pp_H_ZZ_llqq_CMS13(mHh);
        R_pp_H_ZZ_llqq_CMS13=(1+(pp_H_ZZ_TH13-ip_ex_pp_H_ZZ_llqq_CMS13(mHh))/ip_ex_pp_H_ZZ_llqq_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_lnuqq_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh);
        R_ggF_H_WW_lnuqq_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_gg_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_lnuqq_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh);
        R_VBF_H_WW_lnuqq_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_VBF_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_enumunu_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_enumunu_ATLAS13(mHh);
        R_ggF_H_WW_enumunu_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_enumunu_ATLAS13(mHh))/ip_ex_gg_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_enumunu_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh);
        R_VBF_H_WW_enumunu_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh))/ip_ex_VBF_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_ATLAS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_ATLAS13(mHh);
        R_pp_H_hh_bbbb_ATLAS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_ATLAS13(mHh))/ip_ex_pp_H_hh_bbbb_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=1200.0 && mHh<1500.0)
    {
        THoEX_pp_H_Zga_llga_ATLAS8=pp_H_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mHh);
        R_pp_H_Zga_llga_ATLAS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_ATLAS8=ggF_H_WW_TH8/ip_ex_gg_H_WW_ATLAS8(mHh);
        R_ggF_H_WW_ATLAS8=(1+(ggF_H_WW_TH8-ip_ex_gg_H_WW_ATLAS8(mHh))/ip_ex_gg_H_WW_ATLAS8_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_ATLAS8=VBF_H_WW_TH8/ip_ex_VBF_H_WW_ATLAS8(mHh);
        R_VBF_H_WW_ATLAS8=(1+(VBF_H_WW_TH8-ip_ex_VBF_H_WW_ATLAS8(mHh))/ip_ex_VBF_H_WW_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tt_ATLAS8=ggF_H_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mHh);
        R_ggF_H_tt_ATLAS8=(1+(ggF_H_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mHh))/ip_ex_gg_phi_tt_ATLAS8_e(mHh) ) * nftos;

        THoEX_ggF_H_tautau_ATLAS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mHh);
        R_ggF_H_tautau_ATLAS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mHh))/ip_ex_gg_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mHh);
        R_bbF_H_tautau_ATLAS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mHh))/ip_ex_bb_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS13=pp_H_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mHh);
        R_pp_H_gaga_ATLAS13=(1+(pp_H_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mHh))/ip_ex_pp_phi_gaga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS13=ggF_H_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mHh);
        R_ggF_H_gaga_CMS13=(1+(ggF_H_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mHh))/ip_ex_gg_phi_gaga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_llga_ATLAS13=ggF_H_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mHh);
        R_ggF_H_Zga_llga_ATLAS13=(1+(ggF_H_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mHh))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_ATLAS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mHh);
        R_pp_H_Zga_ATLAS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_CMS13(mHh);
        R_pp_H_Zga_llga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_CMS13(mHh))/ip_ex_pp_phi_Zga_llga_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_qqga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_qqga_CMS13(mHh);
        R_pp_H_Zga_qqga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_qqga_CMS13(mHh))/ip_ex_pp_phi_Zga_qqga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_CMS13=ggF_H_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mHh);
        R_ggF_H_Zga_CMS13=(1+(ggF_H_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mHh))/ip_ex_ggF_phi_Zga_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llnunu_CMS13=pp_H_ZZ_TH13/ip_ex_pp_H_ZZ_llnunu_CMS13(mHh);
        R_pp_H_ZZ_llnunu_CMS13=(1+(pp_H_ZZ_TH13-ip_ex_pp_H_ZZ_llnunu_CMS13(mHh))/ip_ex_pp_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_qqllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_qqllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_qqllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_qqllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh);
        R_ggF_H_ZZ_llqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llqq_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh);
        R_VBF_H_ZZ_llqq_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_nunuqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh);
        R_ggF_H_ZZ_nunuqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_nunuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llqq_CMS13=pp_H_ZZ_TH13/ip_ex_pp_H_ZZ_llqq_CMS13(mHh);
        R_pp_H_ZZ_llqq_CMS13=(1+(pp_H_ZZ_TH13-ip_ex_pp_H_ZZ_llqq_CMS13(mHh))/ip_ex_pp_H_ZZ_llqq_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_lnuqq_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh);
        R_ggF_H_WW_lnuqq_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_gg_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_lnuqq_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh);
        R_VBF_H_WW_lnuqq_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_VBF_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_enumunu_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_enumunu_ATLAS13(mHh);
        R_ggF_H_WW_enumunu_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_enumunu_ATLAS13(mHh))/ip_ex_gg_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_enumunu_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh);
        R_VBF_H_WW_enumunu_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh))/ip_ex_VBF_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_VV_qqqq_ATLAS13=pp_H_VV_TH13/ip_ex_pp_H_VV_qqqq_ATLAS13(mHh);
        R_pp_H_VV_qqqq_ATLAS13=(1+(pp_H_VV_TH13-ip_ex_pp_H_VV_qqqq_ATLAS13(mHh))/ip_ex_pp_H_VV_qqqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_ATLAS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_ATLAS13(mHh);
        R_pp_H_hh_bbbb_ATLAS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_ATLAS13(mHh))/ip_ex_pp_H_hh_bbbb_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_bbbb_CMS13=ggF_H_hh_bbbb_TH13/ip_ex_ggF_H_hh_bbbb_CMS13(mHh);
        R_ggF_H_hh_bbbb_CMS13=(1+(ggF_H_hh_bbbb_TH13-ip_ex_ggF_H_hh_bbbb_CMS13(mHh))/ip_ex_ggF_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=1500.0 && mHh<1600.0)
    {
        THoEX_pp_H_Zga_llga_ATLAS8=pp_H_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mHh);
        R_pp_H_Zga_llga_ATLAS8=(1+(pp_H_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mHh) ) * nftos;
        THoEX_ggF_H_tt_ATLAS8=ggF_H_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mHh);
        R_ggF_H_tt_ATLAS8=(1+(ggF_H_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mHh))/ip_ex_gg_phi_tt_ATLAS8_e(mHh) ) * nftos;

        THoEX_ggF_H_tautau_ATLAS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mHh);
        R_ggF_H_tautau_ATLAS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mHh))/ip_ex_gg_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mHh);
        R_bbF_H_tautau_ATLAS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mHh))/ip_ex_bb_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS13=pp_H_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mHh);
        R_pp_H_gaga_ATLAS13=(1+(pp_H_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mHh))/ip_ex_pp_phi_gaga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS13=ggF_H_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mHh);
        R_ggF_H_gaga_CMS13=(1+(ggF_H_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mHh))/ip_ex_gg_phi_gaga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_llga_ATLAS13=ggF_H_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mHh);
        R_ggF_H_Zga_llga_ATLAS13=(1+(ggF_H_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mHh))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_ATLAS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mHh);
        R_pp_H_Zga_ATLAS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_CMS13(mHh);
        R_pp_H_Zga_llga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_CMS13(mHh))/ip_ex_pp_phi_Zga_llga_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_qqga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_qqga_CMS13(mHh);
        R_pp_H_Zga_qqga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_qqga_CMS13(mHh))/ip_ex_pp_phi_Zga_qqga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_CMS13=ggF_H_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mHh);
        R_ggF_H_Zga_CMS13=(1+(ggF_H_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mHh))/ip_ex_ggF_phi_Zga_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llnunu_CMS13=pp_H_ZZ_TH13/ip_ex_pp_H_ZZ_llnunu_CMS13(mHh);
        R_pp_H_ZZ_llnunu_CMS13=(1+(pp_H_ZZ_TH13-ip_ex_pp_H_ZZ_llnunu_CMS13(mHh))/ip_ex_pp_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_qqllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_qqllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_qqllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_qqllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh);
        R_ggF_H_ZZ_llqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llqq_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh);
        R_VBF_H_ZZ_llqq_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_nunuqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh);
        R_ggF_H_ZZ_nunuqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_nunuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llqq_CMS13=pp_H_ZZ_TH13/ip_ex_pp_H_ZZ_llqq_CMS13(mHh);
        R_pp_H_ZZ_llqq_CMS13=(1+(pp_H_ZZ_TH13-ip_ex_pp_H_ZZ_llqq_CMS13(mHh))/ip_ex_pp_H_ZZ_llqq_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_lnuqq_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh);
        R_ggF_H_WW_lnuqq_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_gg_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_lnuqq_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh);
        R_VBF_H_WW_lnuqq_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_VBF_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_enumunu_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_enumunu_ATLAS13(mHh);
        R_ggF_H_WW_enumunu_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_enumunu_ATLAS13(mHh))/ip_ex_gg_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_enumunu_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh);
        R_VBF_H_WW_enumunu_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh))/ip_ex_VBF_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_VV_qqqq_ATLAS13=pp_H_VV_TH13/ip_ex_pp_H_VV_qqqq_ATLAS13(mHh);
        R_pp_H_VV_qqqq_ATLAS13=(1+(pp_H_VV_TH13-ip_ex_pp_H_VV_qqqq_ATLAS13(mHh))/ip_ex_pp_H_VV_qqqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_ATLAS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_ATLAS13(mHh);
        R_pp_H_hh_bbbb_ATLAS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_ATLAS13(mHh))/ip_ex_pp_H_hh_bbbb_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_bbbb_CMS13=ggF_H_hh_bbbb_TH13/ip_ex_ggF_H_hh_bbbb_CMS13(mHh);
        R_ggF_H_hh_bbbb_CMS13=(1+(ggF_H_hh_bbbb_TH13-ip_ex_ggF_H_hh_bbbb_CMS13(mHh))/ip_ex_ggF_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=1600.0 && mHh<2000.0)
    {
        THoEX_ggF_H_tt_ATLAS8=ggF_H_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mHh);
        R_ggF_H_tt_ATLAS8=(1+(ggF_H_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mHh))/ip_ex_gg_phi_tt_ATLAS8_e(mHh) ) * nftos;

        THoEX_ggF_H_tautau_ATLAS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mHh);
        R_ggF_H_tautau_ATLAS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mHh))/ip_ex_gg_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mHh);
        R_bbF_H_tautau_ATLAS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mHh))/ip_ex_bb_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS13=pp_H_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mHh);
        R_pp_H_gaga_ATLAS13=(1+(pp_H_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mHh))/ip_ex_pp_phi_gaga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS13=ggF_H_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mHh);
        R_ggF_H_gaga_CMS13=(1+(ggF_H_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mHh))/ip_ex_gg_phi_gaga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_llga_ATLAS13=ggF_H_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mHh);
        R_ggF_H_Zga_llga_ATLAS13=(1+(ggF_H_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mHh))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_ATLAS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mHh);
        R_pp_H_Zga_ATLAS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_llga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_CMS13(mHh);
        R_pp_H_Zga_llga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_CMS13(mHh))/ip_ex_pp_phi_Zga_llga_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_qqga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_qqga_CMS13(mHh);
        R_pp_H_Zga_qqga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_qqga_CMS13(mHh))/ip_ex_pp_phi_Zga_qqga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_CMS13=ggF_H_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mHh);
        R_ggF_H_Zga_CMS13=(1+(ggF_H_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mHh))/ip_ex_ggF_phi_Zga_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llnunu_CMS13=pp_H_ZZ_TH13/ip_ex_pp_H_ZZ_llnunu_CMS13(mHh);
        R_pp_H_ZZ_llnunu_CMS13=(1+(pp_H_ZZ_TH13-ip_ex_pp_H_ZZ_llnunu_CMS13(mHh))/ip_ex_pp_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_qqllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_qqllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_qqllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_qqllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh);
        R_ggF_H_ZZ_llqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llqq_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh);
        R_VBF_H_ZZ_llqq_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_nunuqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh);
        R_ggF_H_ZZ_nunuqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_nunuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llqq_CMS13=pp_H_ZZ_TH13/ip_ex_pp_H_ZZ_llqq_CMS13(mHh);
        R_pp_H_ZZ_llqq_CMS13=(1+(pp_H_ZZ_TH13-ip_ex_pp_H_ZZ_llqq_CMS13(mHh))/ip_ex_pp_H_ZZ_llqq_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_lnuqq_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh);
        R_ggF_H_WW_lnuqq_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_gg_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_lnuqq_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh);
        R_VBF_H_WW_lnuqq_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_VBF_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_enumunu_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_enumunu_ATLAS13(mHh);
        R_ggF_H_WW_enumunu_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_enumunu_ATLAS13(mHh))/ip_ex_gg_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_enumunu_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh);
        R_VBF_H_WW_enumunu_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh))/ip_ex_VBF_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_VV_qqqq_ATLAS13=pp_H_VV_TH13/ip_ex_pp_H_VV_qqqq_ATLAS13(mHh);
        R_pp_H_VV_qqqq_ATLAS13=(1+(pp_H_VV_TH13-ip_ex_pp_H_VV_qqqq_ATLAS13(mHh))/ip_ex_pp_H_VV_qqqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_ATLAS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_ATLAS13(mHh);
        R_pp_H_hh_bbbb_ATLAS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_ATLAS13(mHh))/ip_ex_pp_H_hh_bbbb_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_bbbb_CMS13=ggF_H_hh_bbbb_TH13/ip_ex_ggF_H_hh_bbbb_CMS13(mHh);
        R_ggF_H_hh_bbbb_CMS13=(1+(ggF_H_hh_bbbb_TH13-ip_ex_ggF_H_hh_bbbb_CMS13(mHh))/ip_ex_ggF_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=2000.0 && mHh<2250.0)
    {
        THoEX_ggF_H_tt_ATLAS8=ggF_H_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mHh);
        R_ggF_H_tt_ATLAS8=(1+(ggF_H_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mHh))/ip_ex_gg_phi_tt_ATLAS8_e(mHh) ) * nftos;

        THoEX_ggF_H_tautau_ATLAS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mHh);
        R_ggF_H_tautau_ATLAS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mHh))/ip_ex_gg_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mHh);
        R_bbF_H_tautau_ATLAS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mHh))/ip_ex_bb_phi_tautau_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS13=pp_H_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mHh);
        R_pp_H_gaga_ATLAS13=(1+(pp_H_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mHh))/ip_ex_pp_phi_gaga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS13=ggF_H_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mHh);
        R_ggF_H_gaga_CMS13=(1+(ggF_H_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mHh))/ip_ex_gg_phi_gaga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_llga_ATLAS13=ggF_H_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mHh);
        R_ggF_H_Zga_llga_ATLAS13=(1+(ggF_H_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mHh))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_ATLAS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mHh);
        R_pp_H_Zga_ATLAS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_qqga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_qqga_CMS13(mHh);
        R_pp_H_Zga_qqga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_qqga_CMS13(mHh))/ip_ex_pp_phi_Zga_qqga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_CMS13=ggF_H_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mHh);
        R_ggF_H_Zga_CMS13=(1+(ggF_H_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mHh))/ip_ex_ggF_phi_Zga_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llnunu_CMS13=pp_H_ZZ_TH13/ip_ex_pp_H_ZZ_llnunu_CMS13(mHh);
        R_pp_H_ZZ_llnunu_CMS13=(1+(pp_H_ZZ_TH13-ip_ex_pp_H_ZZ_llnunu_CMS13(mHh))/ip_ex_pp_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_qqllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_qqllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_qqllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_qqllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh);
        R_ggF_H_ZZ_llqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llqq_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh);
        R_VBF_H_ZZ_llqq_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_nunuqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh);
        R_ggF_H_ZZ_nunuqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_nunuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_lnuqq_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh);
        R_ggF_H_WW_lnuqq_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_gg_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_lnuqq_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh);
        R_VBF_H_WW_lnuqq_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_VBF_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_enumunu_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_enumunu_ATLAS13(mHh);
        R_ggF_H_WW_enumunu_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_enumunu_ATLAS13(mHh))/ip_ex_gg_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_enumunu_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh);
        R_VBF_H_WW_enumunu_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh))/ip_ex_VBF_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_VV_qqqq_ATLAS13=pp_H_VV_TH13/ip_ex_pp_H_VV_qqqq_ATLAS13(mHh);
        R_pp_H_VV_qqqq_ATLAS13=(1+(pp_H_VV_TH13-ip_ex_pp_H_VV_qqqq_ATLAS13(mHh))/ip_ex_pp_H_VV_qqqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_ATLAS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_ATLAS13(mHh);
        R_pp_H_hh_bbbb_ATLAS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_ATLAS13(mHh))/ip_ex_pp_H_hh_bbbb_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_bbbb_CMS13=ggF_H_hh_bbbb_TH13/ip_ex_ggF_H_hh_bbbb_CMS13(mHh);
        R_ggF_H_hh_bbbb_CMS13=(1+(ggF_H_hh_bbbb_TH13-ip_ex_ggF_H_hh_bbbb_CMS13(mHh))/ip_ex_ggF_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=2250.0 && mHh<2400.0)
    {
        THoEX_ggF_H_tt_ATLAS8=ggF_H_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mHh);
        R_ggF_H_tt_ATLAS8=(1+(ggF_H_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mHh))/ip_ex_gg_phi_tt_ATLAS8_e(mHh) ) * nftos;

        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS13=pp_H_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mHh);
        R_pp_H_gaga_ATLAS13=(1+(pp_H_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mHh))/ip_ex_pp_phi_gaga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS13=ggF_H_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mHh);
        R_ggF_H_gaga_CMS13=(1+(ggF_H_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mHh))/ip_ex_gg_phi_gaga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_llga_ATLAS13=ggF_H_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mHh);
        R_ggF_H_Zga_llga_ATLAS13=(1+(ggF_H_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mHh))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_ATLAS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mHh);
        R_pp_H_Zga_ATLAS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mHh))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_qqga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_qqga_CMS13(mHh);
        R_pp_H_Zga_qqga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_qqga_CMS13(mHh))/ip_ex_pp_phi_Zga_qqga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_CMS13=ggF_H_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mHh);
        R_ggF_H_Zga_CMS13=(1+(ggF_H_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mHh))/ip_ex_ggF_phi_Zga_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llnunu_CMS13=pp_H_ZZ_TH13/ip_ex_pp_H_ZZ_llnunu_CMS13(mHh);
        R_pp_H_ZZ_llnunu_CMS13=(1+(pp_H_ZZ_TH13-ip_ex_pp_H_ZZ_llnunu_CMS13(mHh))/ip_ex_pp_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_qqllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_qqllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_qqllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_qqllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh);
        R_ggF_H_ZZ_llqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llqq_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh);
        R_VBF_H_ZZ_llqq_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_nunuqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh);
        R_ggF_H_ZZ_nunuqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_nunuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_lnuqq_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh);
        R_ggF_H_WW_lnuqq_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_gg_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_lnuqq_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh);
        R_VBF_H_WW_lnuqq_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_VBF_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_enumunu_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_enumunu_ATLAS13(mHh);
        R_ggF_H_WW_enumunu_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_enumunu_ATLAS13(mHh))/ip_ex_gg_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_enumunu_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh);
        R_VBF_H_WW_enumunu_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh))/ip_ex_VBF_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_VV_qqqq_ATLAS13=pp_H_VV_TH13/ip_ex_pp_H_VV_qqqq_ATLAS13(mHh);
        R_pp_H_VV_qqqq_ATLAS13=(1+(pp_H_VV_TH13-ip_ex_pp_H_VV_qqqq_ATLAS13(mHh))/ip_ex_pp_H_VV_qqqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_ATLAS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_ATLAS13(mHh);
        R_pp_H_hh_bbbb_ATLAS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_ATLAS13(mHh))/ip_ex_pp_H_hh_bbbb_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_bbbb_CMS13=ggF_H_hh_bbbb_TH13/ip_ex_ggF_H_hh_bbbb_CMS13(mHh);
        R_ggF_H_hh_bbbb_CMS13=(1+(ggF_H_hh_bbbb_TH13-ip_ex_ggF_H_hh_bbbb_CMS13(mHh))/ip_ex_ggF_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=2400.0 && mHh<2500.0)
    {
        THoEX_ggF_H_tt_ATLAS8=ggF_H_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mHh);
        R_ggF_H_tt_ATLAS8=(1+(ggF_H_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mHh))/ip_ex_gg_phi_tt_ATLAS8_e(mHh) ) * nftos;

        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS13=pp_H_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mHh);
        R_pp_H_gaga_ATLAS13=(1+(pp_H_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mHh))/ip_ex_pp_phi_gaga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS13=ggF_H_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mHh);
        R_ggF_H_gaga_CMS13=(1+(ggF_H_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mHh))/ip_ex_gg_phi_gaga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_llga_ATLAS13=ggF_H_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mHh);
        R_ggF_H_Zga_llga_ATLAS13=(1+(ggF_H_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mHh))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_qqga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_qqga_CMS13(mHh);
        R_pp_H_Zga_qqga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_qqga_CMS13(mHh))/ip_ex_pp_phi_Zga_qqga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_CMS13=ggF_H_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mHh);
        R_ggF_H_Zga_CMS13=(1+(ggF_H_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mHh))/ip_ex_ggF_phi_Zga_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llnunu_CMS13=pp_H_ZZ_TH13/ip_ex_pp_H_ZZ_llnunu_CMS13(mHh);
        R_pp_H_ZZ_llnunu_CMS13=(1+(pp_H_ZZ_TH13-ip_ex_pp_H_ZZ_llnunu_CMS13(mHh))/ip_ex_pp_H_ZZ_llnunu_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_qqllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_qqllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_qqllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_qqllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh);
        R_ggF_H_ZZ_llqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llqq_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh);
        R_VBF_H_ZZ_llqq_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_nunuqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh);
        R_ggF_H_ZZ_nunuqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_nunuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_lnuqq_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh);
        R_ggF_H_WW_lnuqq_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_gg_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_lnuqq_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh);
        R_VBF_H_WW_lnuqq_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_VBF_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_enumunu_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_enumunu_ATLAS13(mHh);
        R_ggF_H_WW_enumunu_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_enumunu_ATLAS13(mHh))/ip_ex_gg_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_enumunu_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh);
        R_VBF_H_WW_enumunu_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh))/ip_ex_VBF_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_VV_qqqq_ATLAS13=pp_H_VV_TH13/ip_ex_pp_H_VV_qqqq_ATLAS13(mHh);
        R_pp_H_VV_qqqq_ATLAS13=(1+(pp_H_VV_TH13-ip_ex_pp_H_VV_qqqq_ATLAS13(mHh))/ip_ex_pp_H_VV_qqqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_ATLAS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_ATLAS13(mHh);
        R_pp_H_hh_bbbb_ATLAS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_ATLAS13(mHh))/ip_ex_pp_H_hh_bbbb_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_bbbb_CMS13=ggF_H_hh_bbbb_TH13/ip_ex_ggF_H_hh_bbbb_CMS13(mHh);
        R_ggF_H_hh_bbbb_CMS13=(1+(ggF_H_hh_bbbb_TH13-ip_ex_ggF_H_hh_bbbb_CMS13(mHh))/ip_ex_ggF_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=2500.0 && mHh<2530.0)
    {
        THoEX_ggF_H_tt_ATLAS8=ggF_H_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mHh);
        R_ggF_H_tt_ATLAS8=(1+(ggF_H_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mHh))/ip_ex_gg_phi_tt_ATLAS8_e(mHh) ) * nftos;

        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS13=pp_H_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mHh);
        R_pp_H_gaga_ATLAS13=(1+(pp_H_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mHh))/ip_ex_pp_phi_gaga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS13=ggF_H_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mHh);
        R_ggF_H_gaga_CMS13=(1+(ggF_H_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mHh))/ip_ex_gg_phi_gaga_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_qqga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_qqga_CMS13(mHh);
        R_pp_H_Zga_qqga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_qqga_CMS13(mHh))/ip_ex_pp_phi_Zga_qqga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_CMS13=ggF_H_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mHh);
        R_ggF_H_Zga_CMS13=(1+(ggF_H_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mHh))/ip_ex_ggF_phi_Zga_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_ZZ_llll_CMS13=pp_H_ZZ_llll_TH13/ip_ex_pp_H_ZZ_llll_CMS13(mHh);
        R_pp_H_ZZ_llll_CMS13=(1+(pp_H_ZZ_llll_TH13-ip_ex_pp_H_ZZ_llll_CMS13(mHh))/ip_ex_pp_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_VBF_VH_H_ZZ_llll_CMS13=VBF_VH_H_ZZ_llll_TH13/ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh);
        R_VBF_VH_H_ZZ_llll_CMS13=(1+(VBF_VH_H_ZZ_llll_TH13-ip_ex_VBF_VH_H_ZZ_llll_CMS13(mHh))/ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_qqllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_qqllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_qqllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_qqllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh);
        R_ggF_H_ZZ_llqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llqq_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh);
        R_VBF_H_ZZ_llqq_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_nunuqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh);
        R_ggF_H_ZZ_nunuqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_nunuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_lnuqq_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh);
        R_ggF_H_WW_lnuqq_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_gg_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_lnuqq_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh);
        R_VBF_H_WW_lnuqq_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_VBF_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_enumunu_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_enumunu_ATLAS13(mHh);
        R_ggF_H_WW_enumunu_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_enumunu_ATLAS13(mHh))/ip_ex_gg_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_enumunu_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh);
        R_VBF_H_WW_enumunu_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh))/ip_ex_VBF_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_VV_qqqq_ATLAS13=pp_H_VV_TH13/ip_ex_pp_H_VV_qqqq_ATLAS13(mHh);
        R_pp_H_VV_qqqq_ATLAS13=(1+(pp_H_VV_TH13-ip_ex_pp_H_VV_qqqq_ATLAS13(mHh))/ip_ex_pp_H_VV_qqqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_ATLAS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_ATLAS13(mHh);
        R_pp_H_hh_bbbb_ATLAS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_ATLAS13(mHh))/ip_ex_pp_H_hh_bbbb_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_bbbb_CMS13=ggF_H_hh_bbbb_TH13/ip_ex_ggF_H_hh_bbbb_CMS13(mHh);
        R_ggF_H_hh_bbbb_CMS13=(1+(ggF_H_hh_bbbb_TH13-ip_ex_ggF_H_hh_bbbb_CMS13(mHh))/ip_ex_ggF_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=2530.0 && mHh<2700.0)
    {
        THoEX_ggF_H_tt_ATLAS8=ggF_H_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mHh);
        R_ggF_H_tt_ATLAS8=(1+(ggF_H_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mHh))/ip_ex_gg_phi_tt_ATLAS8_e(mHh) ) * nftos;

        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_gaga_ATLAS13=pp_H_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mHh);
        R_pp_H_gaga_ATLAS13=(1+(pp_H_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mHh))/ip_ex_pp_phi_gaga_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS13=ggF_H_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mHh);
        R_ggF_H_gaga_CMS13=(1+(ggF_H_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mHh))/ip_ex_gg_phi_gaga_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_qqga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_qqga_CMS13(mHh);
        R_pp_H_Zga_qqga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_qqga_CMS13(mHh))/ip_ex_pp_phi_Zga_qqga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_CMS13=ggF_H_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mHh);
        R_ggF_H_Zga_CMS13=(1+(ggF_H_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mHh))/ip_ex_ggF_phi_Zga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_qqllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_qqllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_qqllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_qqllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh);
        R_ggF_H_ZZ_llqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llqq_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh);
        R_VBF_H_ZZ_llqq_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_nunuqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh);
        R_ggF_H_ZZ_nunuqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_nunuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_lnuqq_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh);
        R_ggF_H_WW_lnuqq_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_gg_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_lnuqq_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh);
        R_VBF_H_WW_lnuqq_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_VBF_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_enumunu_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_enumunu_ATLAS13(mHh);
        R_ggF_H_WW_enumunu_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_enumunu_ATLAS13(mHh))/ip_ex_gg_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_enumunu_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh);
        R_VBF_H_WW_enumunu_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh))/ip_ex_VBF_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_VV_qqqq_ATLAS13=pp_H_VV_TH13/ip_ex_pp_H_VV_qqqq_ATLAS13(mHh);
        R_pp_H_VV_qqqq_ATLAS13=(1+(pp_H_VV_TH13-ip_ex_pp_H_VV_qqqq_ATLAS13(mHh))/ip_ex_pp_H_VV_qqqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_ATLAS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_ATLAS13(mHh);
        R_pp_H_hh_bbbb_ATLAS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_ATLAS13(mHh))/ip_ex_pp_H_hh_bbbb_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_bbbb_CMS13=ggF_H_hh_bbbb_TH13/ip_ex_ggF_H_hh_bbbb_CMS13(mHh);
        R_ggF_H_hh_bbbb_CMS13=(1+(ggF_H_hh_bbbb_TH13-ip_ex_ggF_H_hh_bbbb_CMS13(mHh))/ip_ex_ggF_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=2700.0 && mHh<3000.0)
    {
        THoEX_ggF_H_tt_ATLAS8=ggF_H_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mHh);
        R_ggF_H_tt_ATLAS8=(1+(ggF_H_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mHh))/ip_ex_gg_phi_tt_ATLAS8_e(mHh) ) * nftos;

        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS13=ggF_H_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mHh);
        R_ggF_H_gaga_CMS13=(1+(ggF_H_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mHh))/ip_ex_gg_phi_gaga_CMS13_e(mHh) ) * nftos;
        THoEX_pp_H_Zga_qqga_CMS13=pp_H_Zga_TH13/ip_ex_pp_phi_Zga_qqga_CMS13(mHh);
        R_pp_H_Zga_qqga_CMS13=(1+(pp_H_Zga_TH13-ip_ex_pp_phi_Zga_qqga_CMS13(mHh))/ip_ex_pp_phi_Zga_qqga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_CMS13=ggF_H_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mHh);
        R_ggF_H_Zga_CMS13=(1+(ggF_H_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mHh))/ip_ex_ggF_phi_Zga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_qqllnunu_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_ggF_H_ZZ_qqllnunu_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_gg_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_qqllnunu_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh);
        R_VBF_H_ZZ_qqllnunu_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_qqllnunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_llqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh);
        R_ggF_H_ZZ_llqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_ZZ_llqq_ATLAS13=VBF_H_ZZ_TH13/ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh);
        R_VBF_H_ZZ_llqq_ATLAS13=(1+(VBF_H_ZZ_TH13-ip_ex_VBF_H_ZZ_llqq_ATLAS13(mHh))/ip_ex_VBF_H_ZZ_llqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_ZZ_nunuqq_ATLAS13=ggF_H_ZZ_TH13/ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh);
        R_ggF_H_ZZ_nunuqq_ATLAS13=(1+(ggF_H_ZZ_TH13-ip_ex_gg_H_ZZ_nunuqq_ATLAS13(mHh))/ip_ex_gg_H_ZZ_nunuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_lnuqq_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh);
        R_ggF_H_WW_lnuqq_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_gg_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_lnuqq_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh);
        R_VBF_H_WW_lnuqq_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_lnuqq_ATLAS13(mHh))/ip_ex_VBF_H_WW_lnuqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_enumunu_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_enumunu_ATLAS13(mHh);
        R_ggF_H_WW_enumunu_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_enumunu_ATLAS13(mHh))/ip_ex_gg_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_VBF_H_WW_enumunu_ATLAS13=VBF_H_WW_TH13/ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh);
        R_VBF_H_WW_enumunu_ATLAS13=(1+(VBF_H_WW_TH13-ip_ex_VBF_H_WW_enumunu_ATLAS13(mHh))/ip_ex_VBF_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_VV_qqqq_ATLAS13=pp_H_VV_TH13/ip_ex_pp_H_VV_qqqq_ATLAS13(mHh);
        R_pp_H_VV_qqqq_ATLAS13=(1+(pp_H_VV_TH13-ip_ex_pp_H_VV_qqqq_ATLAS13(mHh))/ip_ex_pp_H_VV_qqqq_ATLAS13_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_ATLAS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_ATLAS13(mHh);
        R_pp_H_hh_bbbb_ATLAS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_ATLAS13(mHh))/ip_ex_pp_H_hh_bbbb_ATLAS13_e(mHh) ) * nftos;
        THoEX_ggF_H_hh_bbbb_CMS13=ggF_H_hh_bbbb_TH13/ip_ex_ggF_H_hh_bbbb_CMS13(mHh);
        R_ggF_H_hh_bbbb_CMS13=(1+(ggF_H_hh_bbbb_TH13-ip_ex_ggF_H_hh_bbbb_CMS13(mHh))/ip_ex_ggF_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=3000.0 && mHh<3200.0)
    {
        THoEX_ggF_H_tautau_CMS13=ggF_H_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mHh);
        R_ggF_H_tautau_CMS13=(1+(ggF_H_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mHh))/ip_ex_gg_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_CMS13=bbF_H_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mHh);
        R_bbF_H_tautau_CMS13=(1+(bbF_H_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mHh))/ip_ex_bb_phi_tautau_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_gaga_CMS13=ggF_H_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mHh);
        R_ggF_H_gaga_CMS13=(1+(ggF_H_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mHh))/ip_ex_gg_phi_gaga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_CMS13=ggF_H_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mHh);
        R_ggF_H_Zga_CMS13=(1+(ggF_H_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mHh))/ip_ex_ggF_phi_Zga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_enumunu_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_enumunu_ATLAS13(mHh);
        R_ggF_H_WW_enumunu_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_enumunu_ATLAS13(mHh))/ip_ex_gg_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
    }
    else if(mHh>=3200.0 && mHh<4000.0)
    {
        THoEX_ggF_H_gaga_CMS13=ggF_H_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mHh);
        R_ggF_H_gaga_CMS13=(1+(ggF_H_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mHh))/ip_ex_gg_phi_gaga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_Zga_CMS13=ggF_H_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mHh);
        R_ggF_H_Zga_CMS13=(1+(ggF_H_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mHh))/ip_ex_ggF_phi_Zga_CMS13_e(mHh) ) * nftos;
        THoEX_ggF_H_WW_enumunu_ATLAS13=ggF_H_WW_TH13/ip_ex_gg_H_WW_enumunu_ATLAS13(mHh);
        R_ggF_H_WW_enumunu_ATLAS13=(1+(ggF_H_WW_TH13-ip_ex_gg_H_WW_enumunu_ATLAS13(mHh))/ip_ex_gg_H_WW_enumunu_ATLAS13_e(mHh) ) * nftos;
    }
}

void THDMcache::computeAquantities()
{
    double GF=1/(sqrt(2.0)*vev*vev);
    double mHh=sqrt(mHh2);
    double mA=sqrt(mA2);
    double sinb=tanb/sqrt(1.0+tanb*tanb);
    double cosb=1.0/sqrt(1.0+tanb*tanb);
    double mHp=sqrt(mHp2);
    double sin_ba=sin(bma);
    double cos_ba=cos(bma);

    //These cross sections ratios are necessary for rA_gg
    //gg -> A production cross section ratio at 8 TeV, top loop only over total
    double rSigmaggA_t8 = ip_csr_ggA_t_8(mA);
    //gg -> A production cross section ratio at 8 TeV, bottom loop only over total
    double rSigmaggA_b8 = ip_csr_ggA_b_8(mA);

    /* r_ii is the ratio between the squared 2HDM vertex coupling of the CP-odd
     * Higgs to the particle i and the corresponding coupling of the SM Higgs boson.*/
    double rA_QuQu=(cosb*cosb)/(sinb*sinb);
    double rA_QdQd=0.0;//It depends on the modelType
    double rA_ll=0.0;//It depends on the modelType
    rA_gg=0.0;//It depends on the modelType

    /*Calulation of rA_QdQd, rA_ll, rA_gg, Gamma_Agaga, Gamma_AZga, Gamma_Agg
     * (they depend on the model type): START*/

    /*Gamma_Agaga and Gamma_AZga expressions can be found in
     "The Higgs Hunter's Guide", Appendix C and in arXiv:0902.4665v3, Appendix A;
     *Gamma_Agg expression can be found in arXiv:0902.4665v3, Appendix A*/

    /*I_A_F is needed for Gamma_Agaga;
     * The expression can be found in "The Higgs Hunter's Guide", Appendix C, C.4*/
    gslpp::complex I_A_F=0.0;//It depends on the modelType
    gslpp::complex I_A_Ux=I_A_U(mA2,Mc,Mt);
    gslpp::complex I_A_Dx=I_A_D(mA2,Ms,Mb);
    gslpp::complex I_A_Lx=I_A_L(mA2,Mmu,Mtau);

    /*A_A_F is needed for Gamma_AZga*/
    /*The expression can be found in "The Higgs Hunter's Guide", Appendix C, C.12*/
    gslpp::complex A_A_F = 0.0;//It depends on the modelType
    gslpp::complex A_A_Ux = A_A_U(mA2,cW2,Mc,Mt,MZ);
    gslpp::complex A_A_Dx = A_A_D(mA2,cW2,Ms,Mb,MZ);
    gslpp::complex A_A_Lx = A_A_L(mA2,cW2,Mmu,Mtau,MZ);

    if( modelflag == "type1" ) {
        rA_gg=-cosb/sinb*cosb/sinb+2.0*cosb/sinb*cosb/sinb*(rSigmaggA_t8+rSigmaggA_b8);
        rA_QdQd=cosb/sinb*cosb/sinb;
        rA_ll=cosb/sinb*cosb/sinb;
        I_A_F=cosb/sinb*(I_A_Ux-I_A_Dx-I_A_Lx);
        A_A_F=cosb/sinb*(A_A_Ux-A_A_Dx-A_A_Lx);
    }
    else if( modelflag == "type2" ) {
        rA_gg= 1.0+(cosb/sinb-sinb/cosb)*(rSigmaggA_t8*cosb/sinb-rSigmaggA_b8*sinb/cosb);
        rA_QdQd=sinb/cosb*sinb/cosb;
        rA_ll=sinb/cosb*sinb/cosb;
        I_A_F=cosb/sinb*I_A_Ux+sinb/cosb*(I_A_Dx+I_A_Lx);
        A_A_F=cosb/sinb*A_A_Ux+sinb/cosb*(A_A_Dx+A_A_Lx);
    }
    else if( modelflag == "typeX" ) {
        rA_gg=-cosb/sinb*cosb/sinb+2.0*cosb/sinb*cosb/sinb*(rSigmaggA_t8+rSigmaggA_b8);
        rA_QdQd=cosb/sinb*cosb/sinb;
        rA_ll=sinb/cosb*sinb/cosb;
        I_A_F=cosb/sinb*(I_A_Ux-I_A_Dx)+sinb/cosb*I_A_Lx;
        A_A_F=cosb/sinb*(A_A_Ux-A_A_Dx)+sinb/cosb*A_A_Lx;
    }
    else if( modelflag == "typeY" ) {
        rA_gg=1.0+(cosb/sinb-sinb/cosb)*(rSigmaggA_t8*cosb/sinb-rSigmaggA_b8*sinb/cosb);
        rA_QdQd=sinb/cosb*sinb/cosb;
        rA_ll=cosb/sinb*cosb/sinb;
        I_A_F=cosb/sinb*(I_A_Ux-I_A_Lx)+sinb/cosb*I_A_Dx;
        A_A_F=cosb/sinb*(A_A_Ux-A_A_Lx)+sinb/cosb*A_A_Dx;
    }
    else {
        throw std::runtime_error("modelflag can be only any of \"type1\", \"type2\", \"typeX\" or \"typeY\"");
    }

    /*Gamma_Agaga expression can be found in in arXiv:0902.4665v3, Appendix A, A.8*/
    double Gamma_Agaga=GF*Ale*Ale*mA*mA*mA/(sqrt(2.0)*128.0*M_PI*M_PI*M_PI)
                *(I_A_F).abs2();
    /*Gamma_AZga expression can be found in in arXiv:0902.4665v3, Appendix A, A.9*/
    Gamma_AZga=HSTheta(mA-MZ)*GF*Ale*Ale*mA*mA*mA/(sqrt(2.0)*64.0*M_PI*M_PI*M_PI)
               *(1.0-MZ*MZ/(mA*mA))*(1.0-MZ*MZ/(mA*mA))*(1.0-MZ*MZ/(mA*mA))
               *(A_A_F).abs2();
    /*Gamma_Agg expression can be found in in arXiv:0902.4665v3, Appendix A, A.10*/
    double Gamma_Agg=rA_gg*GF*Als*Als*mA*mA*mA/(sqrt(2.0)*16.0*M_PI*M_PI*M_PI)
                     *(9.0/4.0)*(I_A_Ux/4.0+I_A_Dx).abs2();

    /*Calulation of rA_QdQd, rA_ll, rA_gg, Gamma_Agaga, Gamma_AZga, Gamma_Agg: END*/

    SigmaggF_A8=ip_cs_ggtoA_8(mA)*rA_gg;
    double SigmattF_A8=ip_cs_pptottA_8(mA)*rA_QuQu;
    SigmabbF_A8=ip_cs_pptobbA_8(mA)*rA_QdQd;
    SigmaSumA8 = SigmaggF_A8 + SigmattF_A8 + SigmabbF_A8;

    SigmaggF_A13=ip_cs_ggtoA_13(mA)*rA_gg;
    SigmattF_A13=ip_cs_pptottA_13(mA)*rA_QuQu;
    SigmabbF_A13=ip_cs_pptobbA_13(mA)*rA_QdQd;
    SigmaSumA13 = SigmaggF_A13 + SigmattF_A13 + SigmabbF_A13;

    double BrSM_Atocc=ip_Br_HPtocc(mA);
    double BrSM_Atobb=ip_Br_HPtobb(mA);
    double BrSM_Atott=ip_Br_HPtott(mA);
    double BrSM_Atomumu=ip_Br_HPtomumu(mA);
    double BrSM_Atotautau=ip_Br_HPtotautau(mA);

    double GammaAtotSM=ip_GammaHPtotSM(mA);

    double GammaAHZ=HSTheta(mA-MZ-mHh)*pow(KaellenFunction(mA2,MZ*MZ,mHh*mHh),3)
                    *sin_ba*sin_ba/(2.0*M_PI*vev*vev);

    double GammaAhZ=HSTheta(mA-MZ-sqrt(mHl2))*pow(KaellenFunction(mA2,MZ*MZ,mHl2),3)
                    *cos_ba*cos_ba/(2.0*M_PI*vev*vev);

    double GammaAHpW=2.*HSTheta(mA-MW-mHp)*pow(KaellenFunction(mA2,MW*MW,mHp*mHp),3)
                     /(2.0*M_PI*vev*vev);

    GammaAtot= ((BrSM_Atott+BrSM_Atocc)*rA_QuQu
                    +BrSM_Atobb*rA_QdQd
                    +(BrSM_Atotautau+BrSM_Atomumu)*rA_ll)*GammaAtotSM
               +Gamma_Agaga+Gamma_AZga+Gamma_Agg+GammaAHZ+GammaAhZ+GammaAHpW;

    Br_Atott=BrSM_Atott*rA_QuQu*GammaAtotSM/GammaAtot;
    Br_Atobb=BrSM_Atobb*rA_QdQd*GammaAtotSM/GammaAtot;
    Br_Atotautau=BrSM_Atotautau*rA_ll*GammaAtotSM/GammaAtot;
    Br_Atogaga=Gamma_Agaga/GammaAtot;
    Br_AtoHZ=GammaAHZ/GammaAtot;
    Br_AtohZ=GammaAhZ/GammaAtot;
    Br_AtoHpW=GammaAHpW/GammaAtot;
}

void THDMcache::computeAlimits()
{
    double mHh=sqrt(mHh2);
    double mA=sqrt(mA2);

    double Br_Ztoee=0.03363; //K.A. Olive et al. (Particle Data Group), Chin. Phys. C38, 090001 (2014)
    double Br_Ztomumu=0.03366; //K.A. Olive et al. (Particle Data Group), Chin. Phys. C38, 090001 (2014)

    //Theoretical expressions for the CP-odd Higgs cross sections times branching ratios at 8 TeV

    ggF_A_tautau_TH8=SigmaggF_A8*Br_Atotautau;
    bbF_A_tautau_TH8=SigmabbF_A8*Br_Atotautau;
    pp_A_gaga_TH8=SigmaSumA8*Br_Atogaga;
    ggF_A_gaga_TH8=SigmaggF_A8*Br_Atogaga;
    pp_A_Zga_llga_TH8=SigmaSumA8*Gamma_AZga/GammaAtot*(Br_Ztoee+Br_Ztomumu);
    ggF_A_hZ_bbll_TH8=SigmaggF_A8*Br_AtohZ*THDM_BR_h_bb*(Br_Ztoee+Br_Ztomumu);
    ggF_A_hZ_bbZ_TH8=SigmaggF_A8*Br_AtohZ*THDM_BR_h_bb;
    ggF_A_hZ_tautaull_TH8=SigmaggF_A8*Br_AtohZ*THDM_BR_h_tautau*(Br_Ztoee+Br_Ztomumu);
    ggF_A_hZ_tautauZ_TH8=SigmaggF_A8*Br_AtohZ*THDM_BR_h_tautau;
    ggF_A_tt_TH8=SigmaggF_A8*Br_Atott;
    bbF_A_bb_TH8=SigmabbF_A8*Br_Atobb;
    pp_A_HZ_bbll_TH8=SigmaSumA8*Br_AtoHZ*Br_Htobb*(Br_Ztoee+Br_Ztomumu);
    pp_A_HZ_tautaull_TH8=SigmaSumA8*Br_AtoHZ*Br_Htotautau*(Br_Ztoee+Br_Ztomumu);

    //Ratios of theoretical CP-odd Higgs cross sections and experimental upper limits at 8 TeV

    THoEX_ggF_A_tautau_ATLAS8=0.0;
    R_ggF_A_tautau_ATLAS8=0.0;
    THoEX_ggF_A_tautau_CMS8=0.0;
    R_ggF_A_tautau_CMS8=0.0;
    THoEX_bbF_A_tautau_ATLAS8=0.0;
    R_bbF_A_tautau_ATLAS8=0.0;
    THoEX_bbF_A_tautau_CMS8=0.0;
    R_bbF_A_tautau_CMS8=0.0;
    THoEX_pp_A_gaga_ATLAS8=0.0;
    R_pp_A_gaga_ATLAS8=0.0;
    THoEX_ggF_A_gaga_CMS8=0.0;
    R_ggF_A_gaga_CMS8=0.0;
    THoEX_pp_A_Zga_llga_ATLAS8=0.0;
    R_pp_A_Zga_llga_ATLAS8=0.0;
    THoEX_ggF_A_hZ_bbll_CMS8=0.0;
    R_ggF_A_hZ_bbll_CMS8=0.0;
    THoEX_ggF_A_hZ_bbZ_ATLAS8=0.0;
    R_ggF_A_hZ_bbZ_ATLAS8=0.0;
    THoEX_ggF_A_hZ_tautaull_CMS8=0.0;
    R_ggF_A_hZ_tautaull_CMS8=0.0;
    THoEX_ggF_A_hZ_tautauZ_ATLAS8=0.0;
    R_ggF_A_hZ_tautauZ_ATLAS8=0.0;
    THoEX_ggF_A_tt_ATLAS8=0.0;
    R_ggF_A_tt_ATLAS8=0.0;
    THoEX_bbF_A_bb_CMS8=0.0;
    R_bbF_A_bb_CMS8=0.0;
    THoEX_pp_A_HZ_bbll_CMS8=0.0;
    R_pp_A_HZ_bbll_CMS8=0.0;
    THoEX_pp_A_HZ_tautaull_CMS8=0.0;
    R_pp_A_HZ_tautaull_CMS8=0.0;

    //Theoretical expressions for the CP-odd Higgs cross sections times branching ratios at 13 TeV

    ggF_A_tautau_TH13=SigmaggF_A13*Br_Atotautau;
    bbF_A_tautau_TH13=SigmabbF_A13*Br_Atotautau;
    pp_A_gaga_TH13=SigmaSumA13*Br_Atogaga;
    ggF_A_gaga_TH13=SigmaggF_A13*Br_Atogaga;
    pp_A_Zga_TH13=SigmaSumA13*Gamma_AZga/GammaAtot;
    ggF_A_Zga_TH13=SigmaggF_A13*Gamma_AZga/GammaAtot;
    ggF_A_hZ_bbZ_TH13=SigmaggF_A13*Br_AtohZ*THDM_BR_h_bb;
    bbF_A_hZ_bbZ_TH13=SigmabbF_A13*Br_AtohZ*THDM_BR_h_bb;
    ttF_A_tt_TH13=SigmattF_A13*Br_Atott;
    bbF_A_tt_TH13=SigmabbF_A13*Br_Atott;
    pp_A_bb_TH13=SigmaSumA13*Br_Atobb;

    //Ratios of theoretical CP-odd Higgs cross sections and experimental upper limits at 13 TeV

    THoEX_ttF_A_tt_ATLAS13=0.0;
    R_ttF_A_tt_ATLAS13=0.0;
    THoEX_bbF_A_tt_ATLAS13=0.0;
    R_bbF_A_tt_ATLAS13=0.0;
    THoEX_ggF_A_tautau_ATLAS13=0.0;
    R_ggF_A_tautau_ATLAS13=0.0;
    THoEX_bbF_A_tautau_ATLAS13=0.0;
    R_bbF_A_tautau_ATLAS13=0.0;
    THoEX_ggF_A_tautau_CMS13=0.0;
    R_ggF_A_tautau_CMS13=0.0;
    THoEX_bbF_A_tautau_CMS13=0.0;
    R_bbF_A_tautau_CMS13=0.0;
    THoEX_pp_A_gaga_ATLAS13=0.0;
    R_pp_A_gaga_ATLAS13=0.0;
    THoEX_ggF_A_gaga_CMS13=0.0;
    R_ggF_A_gaga_CMS13=0.0;
    THoEX_pp_A_Zga_llga_ATLAS13=0.0;
    R_pp_A_Zga_llga_ATLAS13=0.0;
    THoEX_ggF_A_Zga_llga_ATLAS13=0.0;
    R_ggF_A_Zga_llga_ATLAS13=0.0;
    THoEX_pp_A_Zga_llga_CMS13=0.0;
    R_pp_A_Zga_llga_CMS13=0.0;
    THoEX_pp_A_Zga_qqga_CMS13=0.0;
    R_pp_A_Zga_qqga_CMS13=0.0;
    THoEX_ggF_A_Zga_CMS13=0.0;
    R_ggF_A_Zga_CMS13=0.0;
    THoEX_ggF_A_hZ_bbZ_ATLAS13=0.0;
    R_ggF_A_hZ_bbZ_ATLAS13=0.0;
    THoEX_bbF_A_hZ_bbZ_ATLAS13=0.0;
    R_bbF_A_hZ_bbZ_ATLAS13=0.0;
    THoEX_pp_A_bb_CMS13=0.0;
    R_pp_A_bb_CMS13=0.0;

    //95% to 1 sigma conversion factor, roughly sqrt(3.84)
    double nftos=1.95996398454;

    if(mA>=65.0 && mA<90.0)
    {
        THoEX_pp_A_gaga_ATLAS8=pp_A_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mA);
        R_pp_A_gaga_ATLAS8=(1+(pp_A_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mA))/ip_ex_pp_phi_gaga_ATLAS8_e(mA) ) * nftos;
        if(mHh>=50.0 && mHh<=1000.0)
        {
            THoEX_pp_A_HZ_tautaull_CMS8=pp_A_HZ_tautaull_TH8/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh);
            R_pp_A_HZ_tautaull_CMS8=(1+(pp_A_HZ_tautaull_TH8-ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }
    }
    else if(mA>=90.0 && mA<100.0)
    {
        THoEX_ggF_A_tautau_ATLAS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mA);
        R_ggF_A_tautau_ATLAS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mA))/ip_ex_gg_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mA);
        R_ggF_A_tautau_CMS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mA))/ip_ex_gg_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mA);
        R_bbF_A_tautau_ATLAS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mA))/ip_ex_bb_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mA);
        R_bbF_A_tautau_CMS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mA))/ip_ex_bb_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS8=pp_A_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mA);
        R_pp_A_gaga_ATLAS8=(1+(pp_A_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mA))/ip_ex_pp_phi_gaga_ATLAS8_e(mA) ) * nftos;
        if(mHh>=50.0 && mHh<=1000.0)
        {
            THoEX_pp_A_HZ_tautaull_CMS8=pp_A_HZ_tautaull_TH8/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh);
            R_pp_A_HZ_tautaull_CMS8=(1+(pp_A_HZ_tautaull_TH8-ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ggF_A_tautau_CMS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mA);
        R_ggF_A_tautau_CMS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mA))/ip_ex_gg_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mA);
        R_bbF_A_tautau_CMS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mA))/ip_ex_bb_phi_tautau_CMS13_e(mA) ) * nftos;
    }
    else if(mA>=100.0 && mA<150.0)
    {
        THoEX_ggF_A_tautau_ATLAS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mA);
        R_ggF_A_tautau_ATLAS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mA))/ip_ex_gg_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mA);
        R_ggF_A_tautau_CMS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mA))/ip_ex_gg_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mA);
        R_bbF_A_tautau_ATLAS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mA))/ip_ex_bb_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mA);
        R_bbF_A_tautau_CMS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mA))/ip_ex_bb_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS8=pp_A_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mA);
        R_pp_A_gaga_ATLAS8=(1+(pp_A_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mA))/ip_ex_pp_phi_gaga_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_bb_CMS8=bbF_A_bb_TH8/ip_ex_bb_phi_bb_CMS8(mA);
        R_bbF_A_bb_CMS8=(1+(bbF_A_bb_TH8-ip_ex_bb_phi_bb_CMS8(mA))/ip_ex_bb_phi_bb_CMS8_e(mA) ) * nftos;
        if(mHh>=50.0 && mHh<=1000.0)
        {
            THoEX_pp_A_HZ_tautaull_CMS8=pp_A_HZ_tautaull_TH8/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh);
            R_pp_A_HZ_tautaull_CMS8=(1+(pp_A_HZ_tautaull_TH8-ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ggF_A_tautau_CMS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mA);
        R_ggF_A_tautau_CMS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mA))/ip_ex_gg_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mA);
        R_bbF_A_tautau_CMS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mA))/ip_ex_bb_phi_tautau_CMS13_e(mA) ) * nftos;
    }
    else if(mA>=150.0 && mA<175.0)
    {
        THoEX_ggF_A_tautau_ATLAS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mA);
        R_ggF_A_tautau_ATLAS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mA))/ip_ex_gg_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mA);
        R_ggF_A_tautau_CMS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mA))/ip_ex_gg_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mA);
        R_bbF_A_tautau_ATLAS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mA))/ip_ex_bb_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mA);
        R_bbF_A_tautau_CMS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mA))/ip_ex_bb_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS8=pp_A_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mA);
        R_pp_A_gaga_ATLAS8=(1+(pp_A_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mA))/ip_ex_pp_phi_gaga_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS8=ggF_A_gaga_TH8/ip_ex_gg_phi_gaga_CMS8(mA);
        R_ggF_A_gaga_CMS8=(1+(ggF_A_gaga_TH8-ip_ex_gg_phi_gaga_CMS8(mA))/ip_ex_gg_phi_gaga_CMS8_e(mA) ) * nftos;
        THoEX_bbF_A_bb_CMS8=bbF_A_bb_TH8/ip_ex_bb_phi_bb_CMS8(mA);
        R_bbF_A_bb_CMS8=(1+(bbF_A_bb_TH8-ip_ex_bb_phi_bb_CMS8(mA))/ip_ex_bb_phi_bb_CMS8_e(mA) ) * nftos;
        if(mHh>=50.0 && mHh<=1000.0)
        {
            THoEX_pp_A_HZ_tautaull_CMS8=pp_A_HZ_tautaull_TH8/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh);
            R_pp_A_HZ_tautaull_CMS8=(1+(pp_A_HZ_tautaull_TH8-ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ggF_A_tautau_CMS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mA);
        R_ggF_A_tautau_CMS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mA))/ip_ex_gg_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mA);
        R_bbF_A_tautau_CMS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mA))/ip_ex_bb_phi_tautau_CMS13_e(mA) ) * nftos;
    }
    else if(mA>=175.0 && mA<200.0)
    {
        THoEX_ggF_A_tautau_ATLAS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mA);
        R_ggF_A_tautau_ATLAS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mA))/ip_ex_gg_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mA);
        R_ggF_A_tautau_CMS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mA))/ip_ex_gg_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mA);
        R_bbF_A_tautau_ATLAS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mA))/ip_ex_bb_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mA);
        R_bbF_A_tautau_CMS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mA))/ip_ex_bb_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS8=pp_A_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mA);
        R_pp_A_gaga_ATLAS8=(1+(pp_A_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mA))/ip_ex_pp_phi_gaga_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS8=ggF_A_gaga_TH8/ip_ex_gg_phi_gaga_CMS8(mA);
        R_ggF_A_gaga_CMS8=(1+(ggF_A_gaga_TH8-ip_ex_gg_phi_gaga_CMS8(mA))/ip_ex_gg_phi_gaga_CMS8_e(mA) ) * nftos;
        THoEX_bbF_A_bb_CMS8=bbF_A_bb_TH8/ip_ex_bb_phi_bb_CMS8(mA);
        R_bbF_A_bb_CMS8=(1+(bbF_A_bb_TH8-ip_ex_bb_phi_bb_CMS8(mA))/ip_ex_bb_phi_bb_CMS8_e(mA) ) * nftos;
        if(mHh>=50.0 && mHh<=910.0)
        {
            THoEX_pp_A_HZ_bbll_CMS8=pp_A_HZ_bbll_TH8/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh);
            R_pp_A_HZ_bbll_CMS8=(1+(pp_A_HZ_bbll_TH8-ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh))/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mHh>=50.0 && mHh<=1000.0)
        {
            THoEX_pp_A_HZ_tautaull_CMS8=pp_A_HZ_tautaull_TH8/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh);
            R_pp_A_HZ_tautaull_CMS8=(1+(pp_A_HZ_tautaull_TH8-ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ggF_A_tautau_CMS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mA);
        R_ggF_A_tautau_CMS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mA))/ip_ex_gg_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mA);
        R_bbF_A_tautau_CMS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mA))/ip_ex_bb_phi_tautau_CMS13_e(mA) ) * nftos;
    }
    else if(mA>=200.0 && mA<220.0)
    {
        THoEX_ggF_A_tautau_ATLAS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mA);
        R_ggF_A_tautau_ATLAS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mA))/ip_ex_gg_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mA);
        R_ggF_A_tautau_CMS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mA))/ip_ex_gg_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mA);
        R_bbF_A_tautau_ATLAS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mA))/ip_ex_bb_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mA);
        R_bbF_A_tautau_CMS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mA))/ip_ex_bb_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS8=pp_A_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mA);
        R_pp_A_gaga_ATLAS8=(1+(pp_A_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mA))/ip_ex_pp_phi_gaga_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS8=ggF_A_gaga_TH8/ip_ex_gg_phi_gaga_CMS8(mA);
        R_ggF_A_gaga_CMS8=(1+(ggF_A_gaga_TH8-ip_ex_gg_phi_gaga_CMS8(mA))/ip_ex_gg_phi_gaga_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_CMS8=pp_A_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mA);
        R_pp_A_Zga_llga_CMS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mA))/ip_ex_pp_A_Zga_llga_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS8=pp_A_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mA);
        R_pp_A_Zga_llga_ATLAS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mA))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_bb_CMS8=bbF_A_bb_TH8/ip_ex_bb_phi_bb_CMS8(mA);
        R_bbF_A_bb_CMS8=(1+(bbF_A_bb_TH8-ip_ex_bb_phi_bb_CMS8(mA))/ip_ex_bb_phi_bb_CMS8_e(mA) ) * nftos;
        if(mHh>=50.0 && mHh<=910.0)
        {
            THoEX_pp_A_HZ_bbll_CMS8=pp_A_HZ_bbll_TH8/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh);
            R_pp_A_HZ_bbll_CMS8=(1+(pp_A_HZ_bbll_TH8-ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh))/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mHh>=50.0 && mHh<=1000.0)
        {
            THoEX_pp_A_HZ_tautaull_CMS8=pp_A_HZ_tautaull_TH8/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh);
            R_pp_A_HZ_tautaull_CMS8=(1+(pp_A_HZ_tautaull_TH8-ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ggF_A_tautau_ATLAS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mA);
        R_ggF_A_tautau_ATLAS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mA))/ip_ex_gg_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mA);
        R_bbF_A_tautau_ATLAS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mA))/ip_ex_bb_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mA);
        R_ggF_A_tautau_CMS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mA))/ip_ex_gg_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mA);
        R_bbF_A_tautau_CMS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mA))/ip_ex_bb_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS13=pp_A_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mA);
        R_pp_A_gaga_ATLAS13=(1+(pp_A_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mA))/ip_ex_pp_phi_gaga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS13=ggF_A_hZ_bbZ_TH13/ip_ex_gg_A_Zh_Zbb_ATLAS13(mA);
        R_ggF_A_hZ_bbZ_ATLAS13=(1+(ggF_A_hZ_bbZ_TH13-ip_ex_gg_A_Zh_Zbb_ATLAS13(mA))/ip_ex_gg_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_hZ_bbZ_ATLAS13=bbF_A_hZ_bbZ_TH13/ip_ex_bb_A_Zh_Zbb_ATLAS13(mA);
        R_bbF_A_hZ_bbZ_ATLAS13=(1+(bbF_A_hZ_bbZ_TH13-ip_ex_bb_A_Zh_Zbb_ATLAS13(mA))/ip_ex_bb_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
    }
    else if(mA>=220.0 && mA<225.0)
    {
        THoEX_ggF_A_tautau_ATLAS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mA);
        R_ggF_A_tautau_ATLAS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mA))/ip_ex_gg_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mA);
        R_ggF_A_tautau_CMS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mA))/ip_ex_gg_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mA);
        R_bbF_A_tautau_ATLAS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mA))/ip_ex_bb_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mA);
        R_bbF_A_tautau_CMS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mA))/ip_ex_bb_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS8=pp_A_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mA);
        R_pp_A_gaga_ATLAS8=(1+(pp_A_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mA))/ip_ex_pp_phi_gaga_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS8=ggF_A_gaga_TH8/ip_ex_gg_phi_gaga_CMS8(mA);
        R_ggF_A_gaga_CMS8=(1+(ggF_A_gaga_TH8-ip_ex_gg_phi_gaga_CMS8(mA))/ip_ex_gg_phi_gaga_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_CMS8=pp_A_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mA);
        R_pp_A_Zga_llga_CMS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mA))/ip_ex_pp_A_Zga_llga_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS8=pp_A_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mA);
        R_pp_A_Zga_llga_ATLAS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mA))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS8=ggF_A_hZ_bbZ_TH8/ip_ex_gg_A_hZ_bbZ_ATLAS8(mA);
        R_ggF_A_hZ_bbZ_ATLAS8=(1+(ggF_A_hZ_bbZ_TH8-ip_ex_gg_A_hZ_bbZ_ATLAS8(mA))/ip_ex_gg_A_hZ_bbZ_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_tautaull_CMS8=ggF_A_hZ_tautaull_TH8/ip_ex_gg_A_hZ_tautaull_CMS8(mA);
        R_ggF_A_hZ_tautaull_CMS8=(1+(ggF_A_hZ_tautaull_TH8-ip_ex_gg_A_hZ_tautaull_CMS8(mA))/ip_ex_gg_A_hZ_tautaull_CMS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_tautauZ_ATLAS8=ggF_A_hZ_tautauZ_TH8/ip_ex_gg_A_hZ_tautauZ_ATLAS8(mA);
        R_ggF_A_hZ_tautauZ_ATLAS8=(1+(ggF_A_hZ_tautauZ_TH8-ip_ex_gg_A_hZ_tautauZ_ATLAS8(mA))/ip_ex_gg_A_hZ_tautauZ_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_bb_CMS8=bbF_A_bb_TH8/ip_ex_bb_phi_bb_CMS8(mA);
        R_bbF_A_bb_CMS8=(1+(bbF_A_bb_TH8-ip_ex_bb_phi_bb_CMS8(mA))/ip_ex_bb_phi_bb_CMS8_e(mA) ) * nftos;
        if(mHh>=50.0 && mHh<=910.0)
        {
            THoEX_pp_A_HZ_bbll_CMS8=pp_A_HZ_bbll_TH8/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh);
            R_pp_A_HZ_bbll_CMS8=(1+(pp_A_HZ_bbll_TH8-ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh))/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mHh>=50.0 && mHh<=1000.0)
        {
            THoEX_pp_A_HZ_tautaull_CMS8=pp_A_HZ_tautaull_TH8/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh);
            R_pp_A_HZ_tautaull_CMS8=(1+(pp_A_HZ_tautaull_TH8-ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ggF_A_tautau_ATLAS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mA);
        R_ggF_A_tautau_ATLAS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mA))/ip_ex_gg_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mA);
        R_bbF_A_tautau_ATLAS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mA))/ip_ex_bb_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mA);
        R_ggF_A_tautau_CMS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mA))/ip_ex_gg_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mA);
        R_bbF_A_tautau_CMS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mA))/ip_ex_bb_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS13=pp_A_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mA);
        R_pp_A_gaga_ATLAS13=(1+(pp_A_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mA))/ip_ex_pp_phi_gaga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS13=ggF_A_hZ_bbZ_TH13/ip_ex_gg_A_Zh_Zbb_ATLAS13(mA);
        R_ggF_A_hZ_bbZ_ATLAS13=(1+(ggF_A_hZ_bbZ_TH13-ip_ex_gg_A_Zh_Zbb_ATLAS13(mA))/ip_ex_gg_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_hZ_bbZ_ATLAS13=bbF_A_hZ_bbZ_TH13/ip_ex_bb_A_Zh_Zbb_ATLAS13(mA);
        R_bbF_A_hZ_bbZ_ATLAS13=(1+(bbF_A_hZ_bbZ_TH13-ip_ex_bb_A_Zh_Zbb_ATLAS13(mA))/ip_ex_bb_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
    }
    else if(mA>=225.0 && mA<250.0)
    {
        THoEX_ggF_A_tautau_ATLAS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mA);
        R_ggF_A_tautau_ATLAS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mA))/ip_ex_gg_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mA);
        R_ggF_A_tautau_CMS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mA))/ip_ex_gg_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mA);
        R_bbF_A_tautau_ATLAS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mA))/ip_ex_bb_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mA);
        R_bbF_A_tautau_CMS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mA))/ip_ex_bb_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS8=pp_A_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mA);
        R_pp_A_gaga_ATLAS8=(1+(pp_A_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mA))/ip_ex_pp_phi_gaga_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS8=ggF_A_gaga_TH8/ip_ex_gg_phi_gaga_CMS8(mA);
        R_ggF_A_gaga_CMS8=(1+(ggF_A_gaga_TH8-ip_ex_gg_phi_gaga_CMS8(mA))/ip_ex_gg_phi_gaga_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_CMS8=pp_A_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mA);
        R_pp_A_Zga_llga_CMS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mA))/ip_ex_pp_A_Zga_llga_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS8=pp_A_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mA);
        R_pp_A_Zga_llga_ATLAS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mA))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbll_CMS8=ggF_A_hZ_bbll_TH8/ip_ex_gg_A_hZ_bbll_CMS8(mA);
        R_ggF_A_hZ_bbll_CMS8=(1+(ggF_A_hZ_bbll_TH8-ip_ex_gg_A_hZ_bbll_CMS8(mA))/ip_ex_gg_A_hZ_bbll_CMS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS8=ggF_A_hZ_bbZ_TH8/ip_ex_gg_A_hZ_bbZ_ATLAS8(mA);
        R_ggF_A_hZ_bbZ_ATLAS8=(1+(ggF_A_hZ_bbZ_TH8-ip_ex_gg_A_hZ_bbZ_ATLAS8(mA))/ip_ex_gg_A_hZ_bbZ_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_tautaull_CMS8=ggF_A_hZ_tautaull_TH8/ip_ex_gg_A_hZ_tautaull_CMS8(mA);
        R_ggF_A_hZ_tautaull_CMS8=(1+(ggF_A_hZ_tautaull_TH8-ip_ex_gg_A_hZ_tautaull_CMS8(mA))/ip_ex_gg_A_hZ_tautaull_CMS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_tautauZ_ATLAS8=ggF_A_hZ_tautauZ_TH8/ip_ex_gg_A_hZ_tautauZ_ATLAS8(mA);
        R_ggF_A_hZ_tautauZ_ATLAS8=(1+(ggF_A_hZ_tautauZ_TH8-ip_ex_gg_A_hZ_tautauZ_ATLAS8(mA))/ip_ex_gg_A_hZ_tautauZ_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_bb_CMS8=bbF_A_bb_TH8/ip_ex_bb_phi_bb_CMS8(mA);
        R_bbF_A_bb_CMS8=(1+(bbF_A_bb_TH8-ip_ex_bb_phi_bb_CMS8(mA))/ip_ex_bb_phi_bb_CMS8_e(mA) ) * nftos;
        if(mHh>=50.0 && mHh<=910.0)
        {
            THoEX_pp_A_HZ_bbll_CMS8=pp_A_HZ_bbll_TH8/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh);
            R_pp_A_HZ_bbll_CMS8=(1+(pp_A_HZ_bbll_TH8-ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh))/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mHh>=50.0 && mHh<=1000.0)
        {
            THoEX_pp_A_HZ_tautaull_CMS8=pp_A_HZ_tautaull_TH8/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh);
            R_pp_A_HZ_tautaull_CMS8=(1+(pp_A_HZ_tautaull_TH8-ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ggF_A_tautau_ATLAS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mA);
        R_ggF_A_tautau_ATLAS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mA))/ip_ex_gg_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mA);
        R_bbF_A_tautau_ATLAS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mA))/ip_ex_bb_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mA);
        R_ggF_A_tautau_CMS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mA))/ip_ex_gg_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mA);
        R_bbF_A_tautau_CMS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mA))/ip_ex_bb_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS13=pp_A_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mA);
        R_pp_A_gaga_ATLAS13=(1+(pp_A_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mA))/ip_ex_pp_phi_gaga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS13=ggF_A_hZ_bbZ_TH13/ip_ex_gg_A_Zh_Zbb_ATLAS13(mA);
        R_ggF_A_hZ_bbZ_ATLAS13=(1+(ggF_A_hZ_bbZ_TH13-ip_ex_gg_A_Zh_Zbb_ATLAS13(mA))/ip_ex_gg_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_hZ_bbZ_ATLAS13=bbF_A_hZ_bbZ_TH13/ip_ex_bb_A_Zh_Zbb_ATLAS13(mA);
        R_bbF_A_hZ_bbZ_ATLAS13=(1+(bbF_A_hZ_bbZ_TH13-ip_ex_bb_A_Zh_Zbb_ATLAS13(mA))/ip_ex_bb_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
    }
    else if(mA>=250.0 && mA<300.0)
    {
        THoEX_ggF_A_tautau_ATLAS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mA);
        R_ggF_A_tautau_ATLAS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mA))/ip_ex_gg_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mA);
        R_ggF_A_tautau_CMS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mA))/ip_ex_gg_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mA);
        R_bbF_A_tautau_ATLAS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mA))/ip_ex_bb_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mA);
        R_bbF_A_tautau_CMS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mA))/ip_ex_bb_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS8=pp_A_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mA);
        R_pp_A_gaga_ATLAS8=(1+(pp_A_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mA))/ip_ex_pp_phi_gaga_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS8=ggF_A_gaga_TH8/ip_ex_gg_phi_gaga_CMS8(mA);
        R_ggF_A_gaga_CMS8=(1+(ggF_A_gaga_TH8-ip_ex_gg_phi_gaga_CMS8(mA))/ip_ex_gg_phi_gaga_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_CMS8=pp_A_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mA);
        R_pp_A_Zga_llga_CMS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mA))/ip_ex_pp_A_Zga_llga_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS8=pp_A_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mA);
        R_pp_A_Zga_llga_ATLAS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mA))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbll_CMS8=ggF_A_hZ_bbll_TH8/ip_ex_gg_A_hZ_bbll_CMS8(mA);
        R_ggF_A_hZ_bbll_CMS8=(1+(ggF_A_hZ_bbll_TH8-ip_ex_gg_A_hZ_bbll_CMS8(mA))/ip_ex_gg_A_hZ_bbll_CMS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS8=ggF_A_hZ_bbZ_TH8/ip_ex_gg_A_hZ_bbZ_ATLAS8(mA);
        R_ggF_A_hZ_bbZ_ATLAS8=(1+(ggF_A_hZ_bbZ_TH8-ip_ex_gg_A_hZ_bbZ_ATLAS8(mA))/ip_ex_gg_A_hZ_bbZ_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_tautaull_CMS8=ggF_A_hZ_tautaull_TH8/ip_ex_gg_A_hZ_tautaull_CMS8(mA);
        R_ggF_A_hZ_tautaull_CMS8=(1+(ggF_A_hZ_tautaull_TH8-ip_ex_gg_A_hZ_tautaull_CMS8(mA))/ip_ex_gg_A_hZ_tautaull_CMS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_tautauZ_ATLAS8=ggF_A_hZ_tautauZ_TH8/ip_ex_gg_A_hZ_tautauZ_ATLAS8(mA);
        R_ggF_A_hZ_tautauZ_ATLAS8=(1+(ggF_A_hZ_tautauZ_TH8-ip_ex_gg_A_hZ_tautauZ_ATLAS8(mA))/ip_ex_gg_A_hZ_tautauZ_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_bb_CMS8=bbF_A_bb_TH8/ip_ex_bb_phi_bb_CMS8(mA);
        R_bbF_A_bb_CMS8=(1+(bbF_A_bb_TH8-ip_ex_bb_phi_bb_CMS8(mA))/ip_ex_bb_phi_bb_CMS8_e(mA) ) * nftos;
        if(mHh>=50.0 && mHh<=910.0)
        {
            THoEX_pp_A_HZ_bbll_CMS8=pp_A_HZ_bbll_TH8/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh);
            R_pp_A_HZ_bbll_CMS8=(1+(pp_A_HZ_bbll_TH8-ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh))/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mHh>=50.0 && mHh<=1000.0)
        {
            THoEX_pp_A_HZ_tautaull_CMS8=pp_A_HZ_tautaull_TH8/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh);
            R_pp_A_HZ_tautaull_CMS8=(1+(pp_A_HZ_tautaull_TH8-ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ggF_A_tautau_ATLAS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mA);
        R_ggF_A_tautau_ATLAS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mA))/ip_ex_gg_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mA);
        R_bbF_A_tautau_ATLAS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mA))/ip_ex_bb_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mA);
        R_ggF_A_tautau_CMS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mA))/ip_ex_gg_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mA);
        R_bbF_A_tautau_CMS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mA))/ip_ex_bb_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS13=pp_A_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mA);
        R_pp_A_gaga_ATLAS13=(1+(pp_A_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mA))/ip_ex_pp_phi_gaga_ATLAS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mA);
        R_pp_A_Zga_llga_ATLAS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mA))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_llga_ATLAS13=ggF_A_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mA);
        R_ggF_A_Zga_llga_ATLAS13=(1+(ggF_A_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mA))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS13=ggF_A_hZ_bbZ_TH13/ip_ex_gg_A_Zh_Zbb_ATLAS13(mA);
        R_ggF_A_hZ_bbZ_ATLAS13=(1+(ggF_A_hZ_bbZ_TH13-ip_ex_gg_A_Zh_Zbb_ATLAS13(mA))/ip_ex_gg_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_hZ_bbZ_ATLAS13=bbF_A_hZ_bbZ_TH13/ip_ex_bb_A_Zh_Zbb_ATLAS13(mA);
        R_bbF_A_hZ_bbZ_ATLAS13=(1+(bbF_A_hZ_bbZ_TH13-ip_ex_bb_A_Zh_Zbb_ATLAS13(mA))/ip_ex_bb_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
    }
    else if(mA>=300.0 && mA<350.0)
    {
        THoEX_ggF_A_tautau_ATLAS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mA);
        R_ggF_A_tautau_ATLAS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mA))/ip_ex_gg_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mA);
        R_ggF_A_tautau_CMS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mA))/ip_ex_gg_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mA);
        R_bbF_A_tautau_ATLAS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mA))/ip_ex_bb_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mA);
        R_bbF_A_tautau_CMS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mA))/ip_ex_bb_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS8=pp_A_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mA);
        R_pp_A_gaga_ATLAS8=(1+(pp_A_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mA))/ip_ex_pp_phi_gaga_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS8=ggF_A_gaga_TH8/ip_ex_gg_phi_gaga_CMS8(mA);
        R_ggF_A_gaga_CMS8=(1+(ggF_A_gaga_TH8-ip_ex_gg_phi_gaga_CMS8(mA))/ip_ex_gg_phi_gaga_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_CMS8=pp_A_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mA);
        R_pp_A_Zga_llga_CMS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mA))/ip_ex_pp_A_Zga_llga_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS8=pp_A_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mA);
        R_pp_A_Zga_llga_ATLAS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mA))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbll_CMS8=ggF_A_hZ_bbll_TH8/ip_ex_gg_A_hZ_bbll_CMS8(mA);
        R_ggF_A_hZ_bbll_CMS8=(1+(ggF_A_hZ_bbll_TH8-ip_ex_gg_A_hZ_bbll_CMS8(mA))/ip_ex_gg_A_hZ_bbll_CMS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS8=ggF_A_hZ_bbZ_TH8/ip_ex_gg_A_hZ_bbZ_ATLAS8(mA);
        R_ggF_A_hZ_bbZ_ATLAS8=(1+(ggF_A_hZ_bbZ_TH8-ip_ex_gg_A_hZ_bbZ_ATLAS8(mA))/ip_ex_gg_A_hZ_bbZ_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_tautaull_CMS8=ggF_A_hZ_tautaull_TH8/ip_ex_gg_A_hZ_tautaull_CMS8(mA);
        R_ggF_A_hZ_tautaull_CMS8=(1+(ggF_A_hZ_tautaull_TH8-ip_ex_gg_A_hZ_tautaull_CMS8(mA))/ip_ex_gg_A_hZ_tautaull_CMS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_tautauZ_ATLAS8=ggF_A_hZ_tautauZ_TH8/ip_ex_gg_A_hZ_tautauZ_ATLAS8(mA);
        R_ggF_A_hZ_tautauZ_ATLAS8=(1+(ggF_A_hZ_tautauZ_TH8-ip_ex_gg_A_hZ_tautauZ_ATLAS8(mA))/ip_ex_gg_A_hZ_tautauZ_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_bb_CMS8=bbF_A_bb_TH8/ip_ex_bb_phi_bb_CMS8(mA);
        R_bbF_A_bb_CMS8=(1+(bbF_A_bb_TH8-ip_ex_bb_phi_bb_CMS8(mA))/ip_ex_bb_phi_bb_CMS8_e(mA) ) * nftos;
        if(mHh>=50.0 && mHh<=910.0)
        {
            THoEX_pp_A_HZ_bbll_CMS8=pp_A_HZ_bbll_TH8/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh);
            R_pp_A_HZ_bbll_CMS8=(1+(pp_A_HZ_bbll_TH8-ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh))/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mHh>=50.0 && mHh<=1000.0)
        {
            THoEX_pp_A_HZ_tautaull_CMS8=pp_A_HZ_tautaull_TH8/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh);
            R_pp_A_HZ_tautaull_CMS8=(1+(pp_A_HZ_tautaull_TH8-ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ggF_A_tautau_ATLAS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mA);
        R_ggF_A_tautau_ATLAS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mA))/ip_ex_gg_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mA);
        R_bbF_A_tautau_ATLAS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mA))/ip_ex_bb_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mA);
        R_ggF_A_tautau_CMS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mA))/ip_ex_gg_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mA);
        R_bbF_A_tautau_CMS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mA))/ip_ex_bb_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS13=pp_A_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mA);
        R_pp_A_gaga_ATLAS13=(1+(pp_A_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mA))/ip_ex_pp_phi_gaga_ATLAS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mA);
        R_pp_A_Zga_llga_ATLAS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mA))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_llga_ATLAS13=ggF_A_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mA);
        R_ggF_A_Zga_llga_ATLAS13=(1+(ggF_A_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mA))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_CMS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_llga_CMS13(mA);
        R_pp_A_Zga_llga_CMS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_llga_CMS13(mA))/ip_ex_pp_phi_Zga_llga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS13=ggF_A_hZ_bbZ_TH13/ip_ex_gg_A_Zh_Zbb_ATLAS13(mA);
        R_ggF_A_hZ_bbZ_ATLAS13=(1+(ggF_A_hZ_bbZ_TH13-ip_ex_gg_A_Zh_Zbb_ATLAS13(mA))/ip_ex_gg_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_hZ_bbZ_ATLAS13=bbF_A_hZ_bbZ_TH13/ip_ex_bb_A_Zh_Zbb_ATLAS13(mA);
        R_bbF_A_hZ_bbZ_ATLAS13=(1+(bbF_A_hZ_bbZ_TH13-ip_ex_bb_A_Zh_Zbb_ATLAS13(mA))/ip_ex_bb_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
    }
    else if(mA>=350.0 && mA<400.0)
    {
        THoEX_ggF_A_tautau_ATLAS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mA);
        R_ggF_A_tautau_ATLAS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mA))/ip_ex_gg_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mA);
        R_ggF_A_tautau_CMS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mA))/ip_ex_gg_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mA);
        R_bbF_A_tautau_ATLAS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mA))/ip_ex_bb_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mA);
        R_bbF_A_tautau_CMS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mA))/ip_ex_bb_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS8=pp_A_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mA);
        R_pp_A_gaga_ATLAS8=(1+(pp_A_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mA))/ip_ex_pp_phi_gaga_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS8=ggF_A_gaga_TH8/ip_ex_gg_phi_gaga_CMS8(mA);
        R_ggF_A_gaga_CMS8=(1+(ggF_A_gaga_TH8-ip_ex_gg_phi_gaga_CMS8(mA))/ip_ex_gg_phi_gaga_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_CMS8=pp_A_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mA);
        R_pp_A_Zga_llga_CMS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mA))/ip_ex_pp_A_Zga_llga_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS8=pp_A_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mA);
        R_pp_A_Zga_llga_ATLAS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mA))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbll_CMS8=ggF_A_hZ_bbll_TH8/ip_ex_gg_A_hZ_bbll_CMS8(mA);
        R_ggF_A_hZ_bbll_CMS8=(1+(ggF_A_hZ_bbll_TH8-ip_ex_gg_A_hZ_bbll_CMS8(mA))/ip_ex_gg_A_hZ_bbll_CMS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS8=ggF_A_hZ_bbZ_TH8/ip_ex_gg_A_hZ_bbZ_ATLAS8(mA);
        R_ggF_A_hZ_bbZ_ATLAS8=(1+(ggF_A_hZ_bbZ_TH8-ip_ex_gg_A_hZ_bbZ_ATLAS8(mA))/ip_ex_gg_A_hZ_bbZ_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_tautauZ_ATLAS8=ggF_A_hZ_tautauZ_TH8/ip_ex_gg_A_hZ_tautauZ_ATLAS8(mA);
        R_ggF_A_hZ_tautauZ_ATLAS8=(1+(ggF_A_hZ_tautauZ_TH8-ip_ex_gg_A_hZ_tautauZ_ATLAS8(mA))/ip_ex_gg_A_hZ_tautauZ_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_bb_CMS8=bbF_A_bb_TH8/ip_ex_bb_phi_bb_CMS8(mA);
        R_bbF_A_bb_CMS8=(1+(bbF_A_bb_TH8-ip_ex_bb_phi_bb_CMS8(mA))/ip_ex_bb_phi_bb_CMS8_e(mA) ) * nftos;
        if(mHh>=50.0 && mHh<=910.0)
        {
            THoEX_pp_A_HZ_bbll_CMS8=pp_A_HZ_bbll_TH8/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh);
            R_pp_A_HZ_bbll_CMS8=(1+(pp_A_HZ_bbll_TH8-ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh))/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mHh>=50.0 && mHh<=1000.0)
        {
            THoEX_pp_A_HZ_tautaull_CMS8=pp_A_HZ_tautaull_TH8/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh);
            R_pp_A_HZ_tautaull_CMS8=(1+(pp_A_HZ_tautaull_TH8-ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ggF_A_tautau_ATLAS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mA);
        R_ggF_A_tautau_ATLAS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mA))/ip_ex_gg_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mA);
        R_bbF_A_tautau_ATLAS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mA))/ip_ex_bb_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mA);
        R_ggF_A_tautau_CMS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mA))/ip_ex_gg_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mA);
        R_bbF_A_tautau_CMS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mA))/ip_ex_bb_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS13=pp_A_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mA);
        R_pp_A_gaga_ATLAS13=(1+(pp_A_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mA))/ip_ex_pp_phi_gaga_ATLAS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mA);
        R_pp_A_Zga_llga_ATLAS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mA))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_llga_ATLAS13=ggF_A_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mA);
        R_ggF_A_Zga_llga_ATLAS13=(1+(ggF_A_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mA))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_CMS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_llga_CMS13(mA);
        R_pp_A_Zga_llga_CMS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_llga_CMS13(mA))/ip_ex_pp_phi_Zga_llga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_CMS13=ggF_A_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mA);
        R_ggF_A_Zga_CMS13=(1+(ggF_A_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mA))/ip_ex_ggF_phi_Zga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS13=ggF_A_hZ_bbZ_TH13/ip_ex_gg_A_Zh_Zbb_ATLAS13(mA);
        R_ggF_A_hZ_bbZ_ATLAS13=(1+(ggF_A_hZ_bbZ_TH13-ip_ex_gg_A_Zh_Zbb_ATLAS13(mA))/ip_ex_gg_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_hZ_bbZ_ATLAS13=bbF_A_hZ_bbZ_TH13/ip_ex_bb_A_Zh_Zbb_ATLAS13(mA);
        R_bbF_A_hZ_bbZ_ATLAS13=(1+(bbF_A_hZ_bbZ_TH13-ip_ex_bb_A_Zh_Zbb_ATLAS13(mA))/ip_ex_bb_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
    }
    else if(mA>=400.0 && mA<500.0)
    {
        THoEX_ggF_A_tautau_ATLAS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mA);
        R_ggF_A_tautau_ATLAS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mA))/ip_ex_gg_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mA);
        R_ggF_A_tautau_CMS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mA))/ip_ex_gg_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mA);
        R_bbF_A_tautau_ATLAS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mA))/ip_ex_bb_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mA);
        R_bbF_A_tautau_CMS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mA))/ip_ex_bb_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS8=pp_A_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mA);
        R_pp_A_gaga_ATLAS8=(1+(pp_A_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mA))/ip_ex_pp_phi_gaga_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS8=ggF_A_gaga_TH8/ip_ex_gg_phi_gaga_CMS8(mA);
        R_ggF_A_gaga_CMS8=(1+(ggF_A_gaga_TH8-ip_ex_gg_phi_gaga_CMS8(mA))/ip_ex_gg_phi_gaga_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_CMS8=pp_A_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mA);
        R_pp_A_Zga_llga_CMS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mA))/ip_ex_pp_A_Zga_llga_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS8=pp_A_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mA);
        R_pp_A_Zga_llga_ATLAS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mA))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbll_CMS8=ggF_A_hZ_bbll_TH8/ip_ex_gg_A_hZ_bbll_CMS8(mA);
        R_ggF_A_hZ_bbll_CMS8=(1+(ggF_A_hZ_bbll_TH8-ip_ex_gg_A_hZ_bbll_CMS8(mA))/ip_ex_gg_A_hZ_bbll_CMS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS8=ggF_A_hZ_bbZ_TH8/ip_ex_gg_A_hZ_bbZ_ATLAS8(mA);
        R_ggF_A_hZ_bbZ_ATLAS8=(1+(ggF_A_hZ_bbZ_TH8-ip_ex_gg_A_hZ_bbZ_ATLAS8(mA))/ip_ex_gg_A_hZ_bbZ_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_tautauZ_ATLAS8=ggF_A_hZ_tautauZ_TH8/ip_ex_gg_A_hZ_tautauZ_ATLAS8(mA);
        R_ggF_A_hZ_tautauZ_ATLAS8=(1+(ggF_A_hZ_tautauZ_TH8-ip_ex_gg_A_hZ_tautauZ_ATLAS8(mA))/ip_ex_gg_A_hZ_tautauZ_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_tt_ATLAS8=ggF_A_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mA);
        R_ggF_A_tt_ATLAS8=(1+(ggF_A_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mA))/ip_ex_gg_phi_tt_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_bb_CMS8=bbF_A_bb_TH8/ip_ex_bb_phi_bb_CMS8(mA);
        R_bbF_A_bb_CMS8=(1+(bbF_A_bb_TH8-ip_ex_bb_phi_bb_CMS8(mA))/ip_ex_bb_phi_bb_CMS8_e(mA) ) * nftos;
        if(mHh>=50.0 && mHh<=910.0)
        {
            THoEX_pp_A_HZ_bbll_CMS8=pp_A_HZ_bbll_TH8/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh);
            R_pp_A_HZ_bbll_CMS8=(1+(pp_A_HZ_bbll_TH8-ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh))/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mHh>=50.0 && mHh<=1000.0)
        {
            THoEX_pp_A_HZ_tautaull_CMS8=pp_A_HZ_tautaull_TH8/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh);
            R_pp_A_HZ_tautaull_CMS8=(1+(pp_A_HZ_tautaull_TH8-ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ttF_A_tt_ATLAS13=ttF_A_tt_TH13/ip_ex_tt_phi_tt_ATLAS13(mA);
        R_ttF_A_tt_ATLAS13=(1+(ttF_A_tt_TH13-ip_ex_tt_phi_tt_ATLAS13(mA))/ip_ex_tt_phi_tt_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_tt_ATLAS13=bbF_A_tt_TH13/ip_ex_bb_phi_tt_ATLAS13(mA);
        R_bbF_A_tt_ATLAS13=(1+(bbF_A_tt_TH13-ip_ex_bb_phi_tt_ATLAS13(mA))/ip_ex_bb_phi_tt_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_ATLAS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mA);
        R_ggF_A_tautau_ATLAS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mA))/ip_ex_gg_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mA);
        R_bbF_A_tautau_ATLAS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mA))/ip_ex_bb_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mA);
        R_ggF_A_tautau_CMS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mA))/ip_ex_gg_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mA);
        R_bbF_A_tautau_CMS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mA))/ip_ex_bb_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS13=pp_A_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mA);
        R_pp_A_gaga_ATLAS13=(1+(pp_A_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mA))/ip_ex_pp_phi_gaga_ATLAS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mA);
        R_pp_A_Zga_llga_ATLAS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mA))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_llga_ATLAS13=ggF_A_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mA);
        R_ggF_A_Zga_llga_ATLAS13=(1+(ggF_A_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mA))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_CMS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_llga_CMS13(mA);
        R_pp_A_Zga_llga_CMS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_llga_CMS13(mA))/ip_ex_pp_phi_Zga_llga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_CMS13=ggF_A_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mA);
        R_ggF_A_Zga_CMS13=(1+(ggF_A_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mA))/ip_ex_ggF_phi_Zga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS13=ggF_A_hZ_bbZ_TH13/ip_ex_gg_A_Zh_Zbb_ATLAS13(mA);
        R_ggF_A_hZ_bbZ_ATLAS13=(1+(ggF_A_hZ_bbZ_TH13-ip_ex_gg_A_Zh_Zbb_ATLAS13(mA))/ip_ex_gg_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_hZ_bbZ_ATLAS13=bbF_A_hZ_bbZ_TH13/ip_ex_bb_A_Zh_Zbb_ATLAS13(mA);
        R_bbF_A_hZ_bbZ_ATLAS13=(1+(bbF_A_hZ_bbZ_TH13-ip_ex_bb_A_Zh_Zbb_ATLAS13(mA))/ip_ex_bb_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
    }
    else if(mA>=500.0 && mA<550.0)
    {
        THoEX_ggF_A_tautau_ATLAS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mA);
        R_ggF_A_tautau_ATLAS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mA))/ip_ex_gg_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mA);
        R_ggF_A_tautau_CMS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mA))/ip_ex_gg_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mA);
        R_bbF_A_tautau_ATLAS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mA))/ip_ex_bb_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mA);
        R_bbF_A_tautau_CMS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mA))/ip_ex_bb_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS8=pp_A_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mA);
        R_pp_A_gaga_ATLAS8=(1+(pp_A_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mA))/ip_ex_pp_phi_gaga_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS8=ggF_A_gaga_TH8/ip_ex_gg_phi_gaga_CMS8(mA);
        R_ggF_A_gaga_CMS8=(1+(ggF_A_gaga_TH8-ip_ex_gg_phi_gaga_CMS8(mA))/ip_ex_gg_phi_gaga_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_CMS8=pp_A_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mA);
        R_pp_A_Zga_llga_CMS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mA))/ip_ex_pp_A_Zga_llga_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS8=pp_A_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mA);
        R_pp_A_Zga_llga_ATLAS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mA))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbll_CMS8=ggF_A_hZ_bbll_TH8/ip_ex_gg_A_hZ_bbll_CMS8(mA);
        R_ggF_A_hZ_bbll_CMS8=(1+(ggF_A_hZ_bbll_TH8-ip_ex_gg_A_hZ_bbll_CMS8(mA))/ip_ex_gg_A_hZ_bbll_CMS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS8=ggF_A_hZ_bbZ_TH8/ip_ex_gg_A_hZ_bbZ_ATLAS8(mA);
        R_ggF_A_hZ_bbZ_ATLAS8=(1+(ggF_A_hZ_bbZ_TH8-ip_ex_gg_A_hZ_bbZ_ATLAS8(mA))/ip_ex_gg_A_hZ_bbZ_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_tautauZ_ATLAS8=ggF_A_hZ_tautauZ_TH8/ip_ex_gg_A_hZ_tautauZ_ATLAS8(mA);
        R_ggF_A_hZ_tautauZ_ATLAS8=(1+(ggF_A_hZ_tautauZ_TH8-ip_ex_gg_A_hZ_tautauZ_ATLAS8(mA))/ip_ex_gg_A_hZ_tautauZ_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_tt_ATLAS8=ggF_A_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mA);
        R_ggF_A_tt_ATLAS8=(1+(ggF_A_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mA))/ip_ex_gg_phi_tt_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_bb_CMS8=bbF_A_bb_TH8/ip_ex_bb_phi_bb_CMS8(mA);
        R_bbF_A_bb_CMS8=(1+(bbF_A_bb_TH8-ip_ex_bb_phi_bb_CMS8(mA))/ip_ex_bb_phi_bb_CMS8_e(mA) ) * nftos;
        if(mHh>=50.0 && mHh<=910.0)
        {
            THoEX_pp_A_HZ_bbll_CMS8=pp_A_HZ_bbll_TH8/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh);
            R_pp_A_HZ_bbll_CMS8=(1+(pp_A_HZ_bbll_TH8-ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh))/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mHh>=50.0 && mHh<=1000.0)
        {
            THoEX_pp_A_HZ_tautaull_CMS8=pp_A_HZ_tautaull_TH8/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh);
            R_pp_A_HZ_tautaull_CMS8=(1+(pp_A_HZ_tautaull_TH8-ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ttF_A_tt_ATLAS13=ttF_A_tt_TH13/ip_ex_tt_phi_tt_ATLAS13(mA);
        R_ttF_A_tt_ATLAS13=(1+(ttF_A_tt_TH13-ip_ex_tt_phi_tt_ATLAS13(mA))/ip_ex_tt_phi_tt_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_tt_ATLAS13=bbF_A_tt_TH13/ip_ex_bb_phi_tt_ATLAS13(mA);
        R_bbF_A_tt_ATLAS13=(1+(bbF_A_tt_TH13-ip_ex_bb_phi_tt_ATLAS13(mA))/ip_ex_bb_phi_tt_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_ATLAS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mA);
        R_ggF_A_tautau_ATLAS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mA))/ip_ex_gg_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mA);
        R_bbF_A_tautau_ATLAS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mA))/ip_ex_bb_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mA);
        R_ggF_A_tautau_CMS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mA))/ip_ex_gg_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mA);
        R_bbF_A_tautau_CMS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mA))/ip_ex_bb_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS13=pp_A_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mA);
        R_pp_A_gaga_ATLAS13=(1+(pp_A_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mA))/ip_ex_pp_phi_gaga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS13=ggF_A_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mA);
        R_ggF_A_gaga_CMS13=(1+(ggF_A_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mA))/ip_ex_gg_phi_gaga_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mA);
        R_pp_A_Zga_llga_ATLAS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mA))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_llga_ATLAS13=ggF_A_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mA);
        R_ggF_A_Zga_llga_ATLAS13=(1+(ggF_A_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mA))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_CMS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_llga_CMS13(mA);
        R_pp_A_Zga_llga_CMS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_llga_CMS13(mA))/ip_ex_pp_phi_Zga_llga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_CMS13=ggF_A_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mA);
        R_ggF_A_Zga_CMS13=(1+(ggF_A_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mA))/ip_ex_ggF_phi_Zga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS13=ggF_A_hZ_bbZ_TH13/ip_ex_gg_A_Zh_Zbb_ATLAS13(mA);
        R_ggF_A_hZ_bbZ_ATLAS13=(1+(ggF_A_hZ_bbZ_TH13-ip_ex_gg_A_Zh_Zbb_ATLAS13(mA))/ip_ex_gg_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_hZ_bbZ_ATLAS13=bbF_A_hZ_bbZ_TH13/ip_ex_bb_A_Zh_Zbb_ATLAS13(mA);
        R_bbF_A_hZ_bbZ_ATLAS13=(1+(bbF_A_hZ_bbZ_TH13-ip_ex_bb_A_Zh_Zbb_ATLAS13(mA))/ip_ex_bb_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
    }
    else if(mA>=550.0 && mA<600.0)
    {
        THoEX_ggF_A_tautau_ATLAS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mA);
        R_ggF_A_tautau_ATLAS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mA))/ip_ex_gg_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mA);
        R_ggF_A_tautau_CMS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mA))/ip_ex_gg_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mA);
        R_bbF_A_tautau_ATLAS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mA))/ip_ex_bb_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mA);
        R_bbF_A_tautau_CMS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mA))/ip_ex_bb_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS8=pp_A_gaga_TH8/ip_ex_pp_phi_gaga_ATLAS8(mA);
        R_pp_A_gaga_ATLAS8=(1+(pp_A_gaga_TH8-ip_ex_pp_phi_gaga_ATLAS8(mA))/ip_ex_pp_phi_gaga_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS8=ggF_A_gaga_TH8/ip_ex_gg_phi_gaga_CMS8(mA);
        R_ggF_A_gaga_CMS8=(1+(ggF_A_gaga_TH8-ip_ex_gg_phi_gaga_CMS8(mA))/ip_ex_gg_phi_gaga_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_CMS8=pp_A_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mA);
        R_pp_A_Zga_llga_CMS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mA))/ip_ex_pp_A_Zga_llga_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS8=pp_A_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mA);
        R_pp_A_Zga_llga_ATLAS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mA))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbll_CMS8=ggF_A_hZ_bbll_TH8/ip_ex_gg_A_hZ_bbll_CMS8(mA);
        R_ggF_A_hZ_bbll_CMS8=(1+(ggF_A_hZ_bbll_TH8-ip_ex_gg_A_hZ_bbll_CMS8(mA))/ip_ex_gg_A_hZ_bbll_CMS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS8=ggF_A_hZ_bbZ_TH8/ip_ex_gg_A_hZ_bbZ_ATLAS8(mA);
        R_ggF_A_hZ_bbZ_ATLAS8=(1+(ggF_A_hZ_bbZ_TH8-ip_ex_gg_A_hZ_bbZ_ATLAS8(mA))/ip_ex_gg_A_hZ_bbZ_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_tautauZ_ATLAS8=ggF_A_hZ_tautauZ_TH8/ip_ex_gg_A_hZ_tautauZ_ATLAS8(mA);
        R_ggF_A_hZ_tautauZ_ATLAS8=(1+(ggF_A_hZ_tautauZ_TH8-ip_ex_gg_A_hZ_tautauZ_ATLAS8(mA))/ip_ex_gg_A_hZ_tautauZ_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_tt_ATLAS8=ggF_A_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mA);
        R_ggF_A_tt_ATLAS8=(1+(ggF_A_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mA))/ip_ex_gg_phi_tt_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_bb_CMS8=bbF_A_bb_TH8/ip_ex_bb_phi_bb_CMS8(mA);
        R_bbF_A_bb_CMS8=(1+(bbF_A_bb_TH8-ip_ex_bb_phi_bb_CMS8(mA))/ip_ex_bb_phi_bb_CMS8_e(mA) ) * nftos;
        if(mHh>=50.0 && mHh<=910.0)
        {
            THoEX_pp_A_HZ_bbll_CMS8=pp_A_HZ_bbll_TH8/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh);
            R_pp_A_HZ_bbll_CMS8=(1+(pp_A_HZ_bbll_TH8-ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh))/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mHh>=50.0 && mHh<=1000.0)
        {
            THoEX_pp_A_HZ_tautaull_CMS8=pp_A_HZ_tautaull_TH8/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh);
            R_pp_A_HZ_tautaull_CMS8=(1+(pp_A_HZ_tautaull_TH8-ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ttF_A_tt_ATLAS13=ttF_A_tt_TH13/ip_ex_tt_phi_tt_ATLAS13(mA);
        R_ttF_A_tt_ATLAS13=(1+(ttF_A_tt_TH13-ip_ex_tt_phi_tt_ATLAS13(mA))/ip_ex_tt_phi_tt_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_tt_ATLAS13=bbF_A_tt_TH13/ip_ex_bb_phi_tt_ATLAS13(mA);
        R_bbF_A_tt_ATLAS13=(1+(bbF_A_tt_TH13-ip_ex_bb_phi_tt_ATLAS13(mA))/ip_ex_bb_phi_tt_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_ATLAS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mA);
        R_ggF_A_tautau_ATLAS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mA))/ip_ex_gg_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mA);
        R_bbF_A_tautau_ATLAS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mA))/ip_ex_bb_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mA);
        R_ggF_A_tautau_CMS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mA))/ip_ex_gg_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mA);
        R_bbF_A_tautau_CMS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mA))/ip_ex_bb_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS13=pp_A_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mA);
        R_pp_A_gaga_ATLAS13=(1+(pp_A_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mA))/ip_ex_pp_phi_gaga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS13=ggF_A_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mA);
        R_ggF_A_gaga_CMS13=(1+(ggF_A_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mA))/ip_ex_gg_phi_gaga_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mA);
        R_pp_A_Zga_llga_ATLAS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mA))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_llga_ATLAS13=ggF_A_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mA);
        R_ggF_A_Zga_llga_ATLAS13=(1+(ggF_A_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mA))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_CMS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_llga_CMS13(mA);
        R_pp_A_Zga_llga_CMS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_llga_CMS13(mA))/ip_ex_pp_phi_Zga_llga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_CMS13=ggF_A_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mA);
        R_ggF_A_Zga_CMS13=(1+(ggF_A_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mA))/ip_ex_ggF_phi_Zga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS13=ggF_A_hZ_bbZ_TH13/ip_ex_gg_A_Zh_Zbb_ATLAS13(mA);
        R_ggF_A_hZ_bbZ_ATLAS13=(1+(ggF_A_hZ_bbZ_TH13-ip_ex_gg_A_Zh_Zbb_ATLAS13(mA))/ip_ex_gg_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_hZ_bbZ_ATLAS13=bbF_A_hZ_bbZ_TH13/ip_ex_bb_A_Zh_Zbb_ATLAS13(mA);
        R_bbF_A_hZ_bbZ_ATLAS13=(1+(bbF_A_hZ_bbZ_TH13-ip_ex_bb_A_Zh_Zbb_ATLAS13(mA))/ip_ex_bb_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
        THoEX_pp_A_bb_CMS13=pp_A_bb_TH13/ip_ex_pp_phi_bb_CMS13(mA);
        R_pp_A_bb_CMS13=(1+(pp_A_bb_TH13-ip_ex_pp_phi_bb_CMS13(mA))/ip_ex_pp_phi_bb_CMS13_e(mA) ) * nftos;
    }
    else if(mA>=600.0 && mA<650.0)
    {
        THoEX_ggF_A_tautau_ATLAS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mA);
        R_ggF_A_tautau_ATLAS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mA))/ip_ex_gg_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mA);
        R_ggF_A_tautau_CMS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mA))/ip_ex_gg_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mA);
        R_bbF_A_tautau_ATLAS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mA))/ip_ex_bb_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mA);
        R_bbF_A_tautau_CMS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mA))/ip_ex_bb_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS8=ggF_A_gaga_TH8/ip_ex_gg_phi_gaga_CMS8(mA);
        R_ggF_A_gaga_CMS8=(1+(ggF_A_gaga_TH8-ip_ex_gg_phi_gaga_CMS8(mA))/ip_ex_gg_phi_gaga_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_CMS8=pp_A_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mA);
        R_pp_A_Zga_llga_CMS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mA))/ip_ex_pp_A_Zga_llga_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS8=pp_A_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mA);
        R_pp_A_Zga_llga_ATLAS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mA))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS8=ggF_A_hZ_bbZ_TH8/ip_ex_gg_A_hZ_bbZ_ATLAS8(mA);
        R_ggF_A_hZ_bbZ_ATLAS8=(1+(ggF_A_hZ_bbZ_TH8-ip_ex_gg_A_hZ_bbZ_ATLAS8(mA))/ip_ex_gg_A_hZ_bbZ_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_tautauZ_ATLAS8=ggF_A_hZ_tautauZ_TH8/ip_ex_gg_A_hZ_tautauZ_ATLAS8(mA);
        R_ggF_A_hZ_tautauZ_ATLAS8=(1+(ggF_A_hZ_tautauZ_TH8-ip_ex_gg_A_hZ_tautauZ_ATLAS8(mA))/ip_ex_gg_A_hZ_tautauZ_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_tt_ATLAS8=ggF_A_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mA);
        R_ggF_A_tt_ATLAS8=(1+(ggF_A_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mA))/ip_ex_gg_phi_tt_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_bb_CMS8=bbF_A_bb_TH8/ip_ex_bb_phi_bb_CMS8(mA);
        R_bbF_A_bb_CMS8=(1+(bbF_A_bb_TH8-ip_ex_bb_phi_bb_CMS8(mA))/ip_ex_bb_phi_bb_CMS8_e(mA) ) * nftos;
        if(mHh>=50.0 && mHh<=910.0)
        {
            THoEX_pp_A_HZ_bbll_CMS8=pp_A_HZ_bbll_TH8/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh);
            R_pp_A_HZ_bbll_CMS8=(1+(pp_A_HZ_bbll_TH8-ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh))/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mHh>=50.0 && mHh<=1000.0)
        {
            THoEX_pp_A_HZ_tautaull_CMS8=pp_A_HZ_tautaull_TH8/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh);
            R_pp_A_HZ_tautaull_CMS8=(1+(pp_A_HZ_tautaull_TH8-ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ttF_A_tt_ATLAS13=ttF_A_tt_TH13/ip_ex_tt_phi_tt_ATLAS13(mA);
        R_ttF_A_tt_ATLAS13=(1+(ttF_A_tt_TH13-ip_ex_tt_phi_tt_ATLAS13(mA))/ip_ex_tt_phi_tt_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_tt_ATLAS13=bbF_A_tt_TH13/ip_ex_bb_phi_tt_ATLAS13(mA);
        R_bbF_A_tt_ATLAS13=(1+(bbF_A_tt_TH13-ip_ex_bb_phi_tt_ATLAS13(mA))/ip_ex_bb_phi_tt_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_ATLAS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mA);
        R_ggF_A_tautau_ATLAS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mA))/ip_ex_gg_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mA);
        R_bbF_A_tautau_ATLAS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mA))/ip_ex_bb_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mA);
        R_ggF_A_tautau_CMS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mA))/ip_ex_gg_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mA);
        R_bbF_A_tautau_CMS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mA))/ip_ex_bb_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS13=pp_A_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mA);
        R_pp_A_gaga_ATLAS13=(1+(pp_A_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mA))/ip_ex_pp_phi_gaga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS13=ggF_A_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mA);
        R_ggF_A_gaga_CMS13=(1+(ggF_A_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mA))/ip_ex_gg_phi_gaga_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mA);
        R_pp_A_Zga_llga_ATLAS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mA))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_llga_ATLAS13=ggF_A_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mA);
        R_ggF_A_Zga_llga_ATLAS13=(1+(ggF_A_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mA))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_CMS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_llga_CMS13(mA);
        R_pp_A_Zga_llga_CMS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_llga_CMS13(mA))/ip_ex_pp_phi_Zga_llga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_CMS13=ggF_A_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mA);
        R_ggF_A_Zga_CMS13=(1+(ggF_A_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mA))/ip_ex_ggF_phi_Zga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS13=ggF_A_hZ_bbZ_TH13/ip_ex_gg_A_Zh_Zbb_ATLAS13(mA);
        R_ggF_A_hZ_bbZ_ATLAS13=(1+(ggF_A_hZ_bbZ_TH13-ip_ex_gg_A_Zh_Zbb_ATLAS13(mA))/ip_ex_gg_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_hZ_bbZ_ATLAS13=bbF_A_hZ_bbZ_TH13/ip_ex_bb_A_Zh_Zbb_ATLAS13(mA);
        R_bbF_A_hZ_bbZ_ATLAS13=(1+(bbF_A_hZ_bbZ_TH13-ip_ex_bb_A_Zh_Zbb_ATLAS13(mA))/ip_ex_bb_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
        THoEX_pp_A_bb_CMS13=pp_A_bb_TH13/ip_ex_pp_phi_bb_CMS13(mA);
        R_pp_A_bb_CMS13=(1+(pp_A_bb_TH13-ip_ex_pp_phi_bb_CMS13(mA))/ip_ex_pp_phi_bb_CMS13_e(mA) ) * nftos;
    }
    else if(mA>=650.0 && mA<850.0)
    {
        THoEX_ggF_A_tautau_ATLAS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mA);
        R_ggF_A_tautau_ATLAS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mA))/ip_ex_gg_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mA);
        R_ggF_A_tautau_CMS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mA))/ip_ex_gg_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mA);
        R_bbF_A_tautau_ATLAS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mA))/ip_ex_bb_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mA);
        R_bbF_A_tautau_CMS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mA))/ip_ex_bb_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS8=ggF_A_gaga_TH8/ip_ex_gg_phi_gaga_CMS8(mA);
        R_ggF_A_gaga_CMS8=(1+(ggF_A_gaga_TH8-ip_ex_gg_phi_gaga_CMS8(mA))/ip_ex_gg_phi_gaga_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_CMS8=pp_A_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mA);
        R_pp_A_Zga_llga_CMS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mA))/ip_ex_pp_A_Zga_llga_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS8=pp_A_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mA);
        R_pp_A_Zga_llga_ATLAS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mA))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS8=ggF_A_hZ_bbZ_TH8/ip_ex_gg_A_hZ_bbZ_ATLAS8(mA);
        R_ggF_A_hZ_bbZ_ATLAS8=(1+(ggF_A_hZ_bbZ_TH8-ip_ex_gg_A_hZ_bbZ_ATLAS8(mA))/ip_ex_gg_A_hZ_bbZ_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_tautauZ_ATLAS8=ggF_A_hZ_tautauZ_TH8/ip_ex_gg_A_hZ_tautauZ_ATLAS8(mA);
        R_ggF_A_hZ_tautauZ_ATLAS8=(1+(ggF_A_hZ_tautauZ_TH8-ip_ex_gg_A_hZ_tautauZ_ATLAS8(mA))/ip_ex_gg_A_hZ_tautauZ_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_tt_ATLAS8=ggF_A_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mA);
        R_ggF_A_tt_ATLAS8=(1+(ggF_A_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mA))/ip_ex_gg_phi_tt_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_bb_CMS8=bbF_A_bb_TH8/ip_ex_bb_phi_bb_CMS8(mA);
        R_bbF_A_bb_CMS8=(1+(bbF_A_bb_TH8-ip_ex_bb_phi_bb_CMS8(mA))/ip_ex_bb_phi_bb_CMS8_e(mA) ) * nftos;
        if(mHh>=50.0 && mHh<=910.0)
        {
            THoEX_pp_A_HZ_bbll_CMS8=pp_A_HZ_bbll_TH8/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh);
            R_pp_A_HZ_bbll_CMS8=(1+(pp_A_HZ_bbll_TH8-ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh))/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mHh>=50.0 && mHh<=1000.0)
        {
            THoEX_pp_A_HZ_tautaull_CMS8=pp_A_HZ_tautaull_TH8/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh);
            R_pp_A_HZ_tautaull_CMS8=(1+(pp_A_HZ_tautaull_TH8-ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ttF_A_tt_ATLAS13=ttF_A_tt_TH13/ip_ex_tt_phi_tt_ATLAS13(mA);
        R_ttF_A_tt_ATLAS13=(1+(ttF_A_tt_TH13-ip_ex_tt_phi_tt_ATLAS13(mA))/ip_ex_tt_phi_tt_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_tt_ATLAS13=bbF_A_tt_TH13/ip_ex_bb_phi_tt_ATLAS13(mA);
        R_bbF_A_tt_ATLAS13=(1+(bbF_A_tt_TH13-ip_ex_bb_phi_tt_ATLAS13(mA))/ip_ex_bb_phi_tt_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_ATLAS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mA);
        R_ggF_A_tautau_ATLAS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mA))/ip_ex_gg_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mA);
        R_bbF_A_tautau_ATLAS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mA))/ip_ex_bb_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mA);
        R_ggF_A_tautau_CMS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mA))/ip_ex_gg_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mA);
        R_bbF_A_tautau_CMS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mA))/ip_ex_bb_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS13=pp_A_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mA);
        R_pp_A_gaga_ATLAS13=(1+(pp_A_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mA))/ip_ex_pp_phi_gaga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS13=ggF_A_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mA);
        R_ggF_A_gaga_CMS13=(1+(ggF_A_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mA))/ip_ex_gg_phi_gaga_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mA);
        R_pp_A_Zga_llga_ATLAS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mA))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_llga_ATLAS13=ggF_A_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mA);
        R_ggF_A_Zga_llga_ATLAS13=(1+(ggF_A_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mA))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_CMS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_llga_CMS13(mA);
        R_pp_A_Zga_llga_CMS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_llga_CMS13(mA))/ip_ex_pp_phi_Zga_llga_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_qqga_CMS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_qqga_CMS13(mA);
        R_pp_A_Zga_qqga_CMS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_qqga_CMS13(mA))/ip_ex_pp_phi_Zga_qqga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_CMS13=ggF_A_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mA);
        R_ggF_A_Zga_CMS13=(1+(ggF_A_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mA))/ip_ex_ggF_phi_Zga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS13=ggF_A_hZ_bbZ_TH13/ip_ex_gg_A_Zh_Zbb_ATLAS13(mA);
        R_ggF_A_hZ_bbZ_ATLAS13=(1+(ggF_A_hZ_bbZ_TH13-ip_ex_gg_A_Zh_Zbb_ATLAS13(mA))/ip_ex_gg_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_hZ_bbZ_ATLAS13=bbF_A_hZ_bbZ_TH13/ip_ex_bb_A_Zh_Zbb_ATLAS13(mA);
        R_bbF_A_hZ_bbZ_ATLAS13=(1+(bbF_A_hZ_bbZ_TH13-ip_ex_bb_A_Zh_Zbb_ATLAS13(mA))/ip_ex_bb_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
        THoEX_pp_A_bb_CMS13=pp_A_bb_TH13/ip_ex_pp_phi_bb_CMS13(mA);
        R_pp_A_bb_CMS13=(1+(pp_A_bb_TH13-ip_ex_pp_phi_bb_CMS13(mA))/ip_ex_pp_phi_bb_CMS13_e(mA) ) * nftos;
    }
    else if(mA>=850.0 && mA<900.0)
    {
        THoEX_ggF_A_tautau_ATLAS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mA);
        R_ggF_A_tautau_ATLAS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mA))/ip_ex_gg_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mA);
        R_ggF_A_tautau_CMS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mA))/ip_ex_gg_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mA);
        R_bbF_A_tautau_ATLAS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mA))/ip_ex_bb_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mA);
        R_bbF_A_tautau_CMS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mA))/ip_ex_bb_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_CMS8=pp_A_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mA);
        R_pp_A_Zga_llga_CMS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mA))/ip_ex_pp_A_Zga_llga_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS8=pp_A_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mA);
        R_pp_A_Zga_llga_ATLAS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mA))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS8=ggF_A_hZ_bbZ_TH8/ip_ex_gg_A_hZ_bbZ_ATLAS8(mA);
        R_ggF_A_hZ_bbZ_ATLAS8=(1+(ggF_A_hZ_bbZ_TH8-ip_ex_gg_A_hZ_bbZ_ATLAS8(mA))/ip_ex_gg_A_hZ_bbZ_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_tautauZ_ATLAS8=ggF_A_hZ_tautauZ_TH8/ip_ex_gg_A_hZ_tautauZ_ATLAS8(mA);
        R_ggF_A_hZ_tautauZ_ATLAS8=(1+(ggF_A_hZ_tautauZ_TH8-ip_ex_gg_A_hZ_tautauZ_ATLAS8(mA))/ip_ex_gg_A_hZ_tautauZ_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_tt_ATLAS8=ggF_A_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mA);
        R_ggF_A_tt_ATLAS8=(1+(ggF_A_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mA))/ip_ex_gg_phi_tt_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_bb_CMS8=bbF_A_bb_TH8/ip_ex_bb_phi_bb_CMS8(mA);
        R_bbF_A_bb_CMS8=(1+(bbF_A_bb_TH8-ip_ex_bb_phi_bb_CMS8(mA))/ip_ex_bb_phi_bb_CMS8_e(mA) ) * nftos;
        if(mHh>=50.0 && mHh<=910.0)
        {
            THoEX_pp_A_HZ_bbll_CMS8=pp_A_HZ_bbll_TH8/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh);
            R_pp_A_HZ_bbll_CMS8=(1+(pp_A_HZ_bbll_TH8-ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh))/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mHh>=50.0 && mHh<=1000.0)
        {
            THoEX_pp_A_HZ_tautaull_CMS8=pp_A_HZ_tautaull_TH8/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh);
            R_pp_A_HZ_tautaull_CMS8=(1+(pp_A_HZ_tautaull_TH8-ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ttF_A_tt_ATLAS13=ttF_A_tt_TH13/ip_ex_tt_phi_tt_ATLAS13(mA);
        R_ttF_A_tt_ATLAS13=(1+(ttF_A_tt_TH13-ip_ex_tt_phi_tt_ATLAS13(mA))/ip_ex_tt_phi_tt_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_tt_ATLAS13=bbF_A_tt_TH13/ip_ex_bb_phi_tt_ATLAS13(mA);
        R_bbF_A_tt_ATLAS13=(1+(bbF_A_tt_TH13-ip_ex_bb_phi_tt_ATLAS13(mA))/ip_ex_bb_phi_tt_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_ATLAS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mA);
        R_ggF_A_tautau_ATLAS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mA))/ip_ex_gg_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mA);
        R_bbF_A_tautau_ATLAS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mA))/ip_ex_bb_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mA);
        R_ggF_A_tautau_CMS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mA))/ip_ex_gg_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mA);
        R_bbF_A_tautau_CMS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mA))/ip_ex_bb_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS13=pp_A_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mA);
        R_pp_A_gaga_ATLAS13=(1+(pp_A_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mA))/ip_ex_pp_phi_gaga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS13=ggF_A_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mA);
        R_ggF_A_gaga_CMS13=(1+(ggF_A_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mA))/ip_ex_gg_phi_gaga_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mA);
        R_pp_A_Zga_llga_ATLAS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mA))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_llga_ATLAS13=ggF_A_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mA);
        R_ggF_A_Zga_llga_ATLAS13=(1+(ggF_A_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mA))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_CMS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_llga_CMS13(mA);
        R_pp_A_Zga_llga_CMS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_llga_CMS13(mA))/ip_ex_pp_phi_Zga_llga_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_qqga_CMS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_qqga_CMS13(mA);
        R_pp_A_Zga_qqga_CMS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_qqga_CMS13(mA))/ip_ex_pp_phi_Zga_qqga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_CMS13=ggF_A_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mA);
        R_ggF_A_Zga_CMS13=(1+(ggF_A_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mA))/ip_ex_ggF_phi_Zga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS13=ggF_A_hZ_bbZ_TH13/ip_ex_gg_A_Zh_Zbb_ATLAS13(mA);
        R_ggF_A_hZ_bbZ_ATLAS13=(1+(ggF_A_hZ_bbZ_TH13-ip_ex_gg_A_Zh_Zbb_ATLAS13(mA))/ip_ex_gg_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_hZ_bbZ_ATLAS13=bbF_A_hZ_bbZ_TH13/ip_ex_bb_A_Zh_Zbb_ATLAS13(mA);
        R_bbF_A_hZ_bbZ_ATLAS13=(1+(bbF_A_hZ_bbZ_TH13-ip_ex_bb_A_Zh_Zbb_ATLAS13(mA))/ip_ex_bb_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
        THoEX_pp_A_bb_CMS13=pp_A_bb_TH13/ip_ex_pp_phi_bb_CMS13(mA);
        R_pp_A_bb_CMS13=(1+(pp_A_bb_TH13-ip_ex_pp_phi_bb_CMS13(mA))/ip_ex_pp_phi_bb_CMS13_e(mA) ) * nftos;
    }
    else if(mA>=900.0 && mA<1000.0)
    {
        THoEX_ggF_A_tautau_ATLAS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mA);
        R_ggF_A_tautau_ATLAS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mA))/ip_ex_gg_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS8=ggF_A_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(mA);
        R_ggF_A_tautau_CMS8=(1+(ggF_A_tautau_TH8-ip_ex_gg_phi_tautau_CMS8(mA))/ip_ex_gg_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mA);
        R_bbF_A_tautau_ATLAS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mA))/ip_ex_bb_phi_tautau_ATLAS8_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS8=bbF_A_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(mA);
        R_bbF_A_tautau_CMS8=(1+(bbF_A_tautau_TH8-ip_ex_bb_phi_tautau_CMS8(mA))/ip_ex_bb_phi_tautau_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_CMS8=pp_A_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mA);
        R_pp_A_Zga_llga_CMS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mA))/ip_ex_pp_A_Zga_llga_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS8=pp_A_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mA);
        R_pp_A_Zga_llga_ATLAS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mA))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS8=ggF_A_hZ_bbZ_TH8/ip_ex_gg_A_hZ_bbZ_ATLAS8(mA);
        R_ggF_A_hZ_bbZ_ATLAS8=(1+(ggF_A_hZ_bbZ_TH8-ip_ex_gg_A_hZ_bbZ_ATLAS8(mA))/ip_ex_gg_A_hZ_bbZ_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_tautauZ_ATLAS8=ggF_A_hZ_tautauZ_TH8/ip_ex_gg_A_hZ_tautauZ_ATLAS8(mA);
        R_ggF_A_hZ_tautauZ_ATLAS8=(1+(ggF_A_hZ_tautauZ_TH8-ip_ex_gg_A_hZ_tautauZ_ATLAS8(mA))/ip_ex_gg_A_hZ_tautauZ_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_tt_ATLAS8=ggF_A_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mA);
        R_ggF_A_tt_ATLAS8=(1+(ggF_A_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mA))/ip_ex_gg_phi_tt_ATLAS8_e(mA) ) * nftos;
        if(mHh>=50.0 && mHh<=910.0)
        {
            THoEX_pp_A_HZ_bbll_CMS8=pp_A_HZ_bbll_TH8/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh);
            R_pp_A_HZ_bbll_CMS8=(1+(pp_A_HZ_bbll_TH8-ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh))/ip_ex_pp_A_HZ_bbll_CMS8(mA,mHh) ) * nftos;
        }
        if(mHh>=50.0 && mHh<=1000.0)
        {
            THoEX_pp_A_HZ_tautaull_CMS8=pp_A_HZ_tautaull_TH8/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh);
            R_pp_A_HZ_tautaull_CMS8=(1+(pp_A_HZ_tautaull_TH8-ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh))/ip_ex_pp_A_HZ_tautaull_CMS8(mA,mHh) ) * nftos;
        }

        THoEX_ttF_A_tt_ATLAS13=ttF_A_tt_TH13/ip_ex_tt_phi_tt_ATLAS13(mA);
        R_ttF_A_tt_ATLAS13=(1+(ttF_A_tt_TH13-ip_ex_tt_phi_tt_ATLAS13(mA))/ip_ex_tt_phi_tt_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_tt_ATLAS13=bbF_A_tt_TH13/ip_ex_bb_phi_tt_ATLAS13(mA);
        R_bbF_A_tt_ATLAS13=(1+(bbF_A_tt_TH13-ip_ex_bb_phi_tt_ATLAS13(mA))/ip_ex_bb_phi_tt_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_ATLAS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mA);
        R_ggF_A_tautau_ATLAS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mA))/ip_ex_gg_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mA);
        R_bbF_A_tautau_ATLAS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mA))/ip_ex_bb_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mA);
        R_ggF_A_tautau_CMS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mA))/ip_ex_gg_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mA);
        R_bbF_A_tautau_CMS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mA))/ip_ex_bb_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS13=pp_A_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mA);
        R_pp_A_gaga_ATLAS13=(1+(pp_A_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mA))/ip_ex_pp_phi_gaga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS13=ggF_A_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mA);
        R_ggF_A_gaga_CMS13=(1+(ggF_A_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mA))/ip_ex_gg_phi_gaga_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mA);
        R_pp_A_Zga_llga_ATLAS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mA))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_llga_ATLAS13=ggF_A_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mA);
        R_ggF_A_Zga_llga_ATLAS13=(1+(ggF_A_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mA))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_CMS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_llga_CMS13(mA);
        R_pp_A_Zga_llga_CMS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_llga_CMS13(mA))/ip_ex_pp_phi_Zga_llga_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_qqga_CMS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_qqga_CMS13(mA);
        R_pp_A_Zga_qqga_CMS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_qqga_CMS13(mA))/ip_ex_pp_phi_Zga_qqga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_CMS13=ggF_A_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mA);
        R_ggF_A_Zga_CMS13=(1+(ggF_A_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mA))/ip_ex_ggF_phi_Zga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS13=ggF_A_hZ_bbZ_TH13/ip_ex_gg_A_Zh_Zbb_ATLAS13(mA);
        R_ggF_A_hZ_bbZ_ATLAS13=(1+(ggF_A_hZ_bbZ_TH13-ip_ex_gg_A_Zh_Zbb_ATLAS13(mA))/ip_ex_gg_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_hZ_bbZ_ATLAS13=bbF_A_hZ_bbZ_TH13/ip_ex_bb_A_Zh_Zbb_ATLAS13(mA);
        R_bbF_A_hZ_bbZ_ATLAS13=(1+(bbF_A_hZ_bbZ_TH13-ip_ex_bb_A_Zh_Zbb_ATLAS13(mA))/ip_ex_bb_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
        THoEX_pp_A_bb_CMS13=pp_A_bb_TH13/ip_ex_pp_phi_bb_CMS13(mA);
        R_pp_A_bb_CMS13=(1+(pp_A_bb_TH13-ip_ex_pp_phi_bb_CMS13(mA))/ip_ex_pp_phi_bb_CMS13_e(mA) ) * nftos;
    }
    else if(mA>=1000.0 && mA<1200.0)
    {
        THoEX_pp_A_Zga_llga_CMS8=pp_A_Zga_llga_TH8/ip_ex_pp_A_Zga_llga_CMS8(mA);
        R_pp_A_Zga_llga_CMS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_A_Zga_llga_CMS8(mA))/ip_ex_pp_A_Zga_llga_CMS8_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS8=pp_A_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mA);
        R_pp_A_Zga_llga_ATLAS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mA))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_tt_ATLAS8=ggF_A_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mA);
        R_ggF_A_tt_ATLAS8=(1+(ggF_A_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mA))/ip_ex_gg_phi_tt_ATLAS8_e(mA) ) * nftos;

        THoEX_ggF_A_tautau_ATLAS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mA);
        R_ggF_A_tautau_ATLAS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mA))/ip_ex_gg_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mA);
        R_bbF_A_tautau_ATLAS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mA))/ip_ex_bb_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mA);
        R_ggF_A_tautau_CMS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mA))/ip_ex_gg_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mA);
        R_bbF_A_tautau_CMS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mA))/ip_ex_bb_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS13=pp_A_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mA);
        R_pp_A_gaga_ATLAS13=(1+(pp_A_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mA))/ip_ex_pp_phi_gaga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS13=ggF_A_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mA);
        R_ggF_A_gaga_CMS13=(1+(ggF_A_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mA))/ip_ex_gg_phi_gaga_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mA);
        R_pp_A_Zga_llga_ATLAS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mA))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_llga_ATLAS13=ggF_A_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mA);
        R_ggF_A_Zga_llga_ATLAS13=(1+(ggF_A_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mA))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_CMS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_llga_CMS13(mA);
        R_pp_A_Zga_llga_CMS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_llga_CMS13(mA))/ip_ex_pp_phi_Zga_llga_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_qqga_CMS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_qqga_CMS13(mA);
        R_pp_A_Zga_qqga_CMS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_qqga_CMS13(mA))/ip_ex_pp_phi_Zga_qqga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_CMS13=ggF_A_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mA);
        R_ggF_A_Zga_CMS13=(1+(ggF_A_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mA))/ip_ex_ggF_phi_Zga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS13=ggF_A_hZ_bbZ_TH13/ip_ex_gg_A_Zh_Zbb_ATLAS13(mA);
        R_ggF_A_hZ_bbZ_ATLAS13=(1+(ggF_A_hZ_bbZ_TH13-ip_ex_gg_A_Zh_Zbb_ATLAS13(mA))/ip_ex_gg_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_hZ_bbZ_ATLAS13=bbF_A_hZ_bbZ_TH13/ip_ex_bb_A_Zh_Zbb_ATLAS13(mA);
        R_bbF_A_hZ_bbZ_ATLAS13=(1+(bbF_A_hZ_bbZ_TH13-ip_ex_bb_A_Zh_Zbb_ATLAS13(mA))/ip_ex_bb_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
        THoEX_pp_A_bb_CMS13=pp_A_bb_TH13/ip_ex_pp_phi_bb_CMS13(mA);
        R_pp_A_bb_CMS13=(1+(pp_A_bb_TH13-ip_ex_pp_phi_bb_CMS13(mA))/ip_ex_pp_phi_bb_CMS13_e(mA) ) * nftos;
    }
    else if(mA>=1200.0 && mA<1600.0)
    {
        THoEX_pp_A_Zga_llga_ATLAS8=pp_A_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(mA);
        R_pp_A_Zga_llga_ATLAS8=(1+(pp_A_Zga_llga_TH8-ip_ex_pp_phi_Zga_llga_ATLAS8(mA))/ip_ex_pp_phi_Zga_llga_ATLAS8_e(mA) ) * nftos;
        THoEX_ggF_A_tt_ATLAS8=ggF_A_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mA);
        R_ggF_A_tt_ATLAS8=(1+(ggF_A_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mA))/ip_ex_gg_phi_tt_ATLAS8_e(mA) ) * nftos;

        THoEX_ggF_A_tautau_ATLAS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mA);
        R_ggF_A_tautau_ATLAS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mA))/ip_ex_gg_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mA);
        R_bbF_A_tautau_ATLAS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mA))/ip_ex_bb_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mA);
        R_ggF_A_tautau_CMS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mA))/ip_ex_gg_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mA);
        R_bbF_A_tautau_CMS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mA))/ip_ex_bb_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS13=pp_A_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mA);
        R_pp_A_gaga_ATLAS13=(1+(pp_A_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mA))/ip_ex_pp_phi_gaga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS13=ggF_A_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mA);
        R_ggF_A_gaga_CMS13=(1+(ggF_A_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mA))/ip_ex_gg_phi_gaga_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mA);
        R_pp_A_Zga_llga_ATLAS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mA))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_llga_ATLAS13=ggF_A_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mA);
        R_ggF_A_Zga_llga_ATLAS13=(1+(ggF_A_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mA))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_CMS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_llga_CMS13(mA);
        R_pp_A_Zga_llga_CMS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_llga_CMS13(mA))/ip_ex_pp_phi_Zga_llga_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_qqga_CMS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_qqga_CMS13(mA);
        R_pp_A_Zga_qqga_CMS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_qqga_CMS13(mA))/ip_ex_pp_phi_Zga_qqga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_CMS13=ggF_A_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mA);
        R_ggF_A_Zga_CMS13=(1+(ggF_A_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mA))/ip_ex_ggF_phi_Zga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS13=ggF_A_hZ_bbZ_TH13/ip_ex_gg_A_Zh_Zbb_ATLAS13(mA);
        R_ggF_A_hZ_bbZ_ATLAS13=(1+(ggF_A_hZ_bbZ_TH13-ip_ex_gg_A_Zh_Zbb_ATLAS13(mA))/ip_ex_gg_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_hZ_bbZ_ATLAS13=bbF_A_hZ_bbZ_TH13/ip_ex_bb_A_Zh_Zbb_ATLAS13(mA);
        R_bbF_A_hZ_bbZ_ATLAS13=(1+(bbF_A_hZ_bbZ_TH13-ip_ex_bb_A_Zh_Zbb_ATLAS13(mA))/ip_ex_bb_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
    }
    else if(mA>=1600.0 && mA<2000.0)
    {
        THoEX_ggF_A_tt_ATLAS8=ggF_A_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mA);
        R_ggF_A_tt_ATLAS8=(1+(ggF_A_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mA))/ip_ex_gg_phi_tt_ATLAS8_e(mA) ) * nftos;

        THoEX_ggF_A_tautau_ATLAS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mA);
        R_ggF_A_tautau_ATLAS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mA))/ip_ex_gg_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mA);
        R_bbF_A_tautau_ATLAS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mA))/ip_ex_bb_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mA);
        R_ggF_A_tautau_CMS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mA))/ip_ex_gg_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mA);
        R_bbF_A_tautau_CMS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mA))/ip_ex_bb_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS13=pp_A_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mA);
        R_pp_A_gaga_ATLAS13=(1+(pp_A_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mA))/ip_ex_pp_phi_gaga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS13=ggF_A_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mA);
        R_ggF_A_gaga_CMS13=(1+(ggF_A_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mA))/ip_ex_gg_phi_gaga_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mA);
        R_pp_A_Zga_llga_ATLAS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mA))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_llga_ATLAS13=ggF_A_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mA);
        R_ggF_A_Zga_llga_ATLAS13=(1+(ggF_A_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mA))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_CMS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_llga_CMS13(mA);
        R_pp_A_Zga_llga_CMS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_llga_CMS13(mA))/ip_ex_pp_phi_Zga_llga_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_qqga_CMS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_qqga_CMS13(mA);
        R_pp_A_Zga_qqga_CMS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_qqga_CMS13(mA))/ip_ex_pp_phi_Zga_qqga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_CMS13=ggF_A_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mA);
        R_ggF_A_Zga_CMS13=(1+(ggF_A_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mA))/ip_ex_ggF_phi_Zga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_hZ_bbZ_ATLAS13=ggF_A_hZ_bbZ_TH13/ip_ex_gg_A_Zh_Zbb_ATLAS13(mA);
        R_ggF_A_hZ_bbZ_ATLAS13=(1+(ggF_A_hZ_bbZ_TH13-ip_ex_gg_A_Zh_Zbb_ATLAS13(mA))/ip_ex_gg_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_hZ_bbZ_ATLAS13=bbF_A_hZ_bbZ_TH13/ip_ex_bb_A_Zh_Zbb_ATLAS13(mA);
        R_bbF_A_hZ_bbZ_ATLAS13=(1+(bbF_A_hZ_bbZ_TH13-ip_ex_bb_A_Zh_Zbb_ATLAS13(mA))/ip_ex_bb_A_Zh_Zbb_ATLAS13_e(mA) ) * nftos;
    }
    else if(mA>=2000.0 && mA<2250.0)
    {
        THoEX_ggF_A_tt_ATLAS8=ggF_A_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mA);
        R_ggF_A_tt_ATLAS8=(1+(ggF_A_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mA))/ip_ex_gg_phi_tt_ATLAS8_e(mA) ) * nftos;

        THoEX_ggF_A_tautau_ATLAS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(mA);
        R_ggF_A_tautau_ATLAS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_ATLAS13(mA))/ip_ex_gg_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_ATLAS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(mA);
        R_bbF_A_tautau_ATLAS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_ATLAS13(mA))/ip_ex_bb_phi_tautau_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_tautau_CMS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mA);
        R_ggF_A_tautau_CMS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mA))/ip_ex_gg_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mA);
        R_bbF_A_tautau_CMS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mA))/ip_ex_bb_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS13=pp_A_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mA);
        R_pp_A_gaga_ATLAS13=(1+(pp_A_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mA))/ip_ex_pp_phi_gaga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS13=ggF_A_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mA);
        R_ggF_A_gaga_CMS13=(1+(ggF_A_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mA))/ip_ex_gg_phi_gaga_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mA);
        R_pp_A_Zga_llga_ATLAS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mA))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_llga_ATLAS13=ggF_A_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mA);
        R_ggF_A_Zga_llga_ATLAS13=(1+(ggF_A_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mA))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_qqga_CMS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_qqga_CMS13(mA);
        R_pp_A_Zga_qqga_CMS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_qqga_CMS13(mA))/ip_ex_pp_phi_Zga_qqga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_CMS13=ggF_A_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mA);
        R_ggF_A_Zga_CMS13=(1+(ggF_A_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mA))/ip_ex_ggF_phi_Zga_CMS13_e(mA) ) * nftos;
    }
    else if(mA>=2250.0 && mA<2400.0)
    {
        THoEX_ggF_A_tt_ATLAS8=ggF_A_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mA);
        R_ggF_A_tt_ATLAS8=(1+(ggF_A_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mA))/ip_ex_gg_phi_tt_ATLAS8_e(mA) ) * nftos;

        THoEX_ggF_A_tautau_CMS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mA);
        R_ggF_A_tautau_CMS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mA))/ip_ex_gg_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mA);
        R_bbF_A_tautau_CMS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mA))/ip_ex_bb_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS13=pp_A_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mA);
        R_pp_A_gaga_ATLAS13=(1+(pp_A_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mA))/ip_ex_pp_phi_gaga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS13=ggF_A_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mA);
        R_ggF_A_gaga_CMS13=(1+(ggF_A_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mA))/ip_ex_gg_phi_gaga_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_llga_ATLAS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_llga_ATLAS13(mA);
        R_pp_A_Zga_llga_ATLAS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_llga_ATLAS13(mA))/ip_ex_pp_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_llga_ATLAS13=ggF_A_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mA);
        R_ggF_A_Zga_llga_ATLAS13=(1+(ggF_A_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mA))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_qqga_CMS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_qqga_CMS13(mA);
        R_pp_A_Zga_qqga_CMS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_qqga_CMS13(mA))/ip_ex_pp_phi_Zga_qqga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_CMS13=ggF_A_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mA);
        R_ggF_A_Zga_CMS13=(1+(ggF_A_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mA))/ip_ex_ggF_phi_Zga_CMS13_e(mA) ) * nftos;
    }
    else if(mA>=2400.0 && mA<2500.0)
    {
        THoEX_ggF_A_tt_ATLAS8=ggF_A_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mA);
        R_ggF_A_tt_ATLAS8=(1+(ggF_A_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mA))/ip_ex_gg_phi_tt_ATLAS8_e(mA) ) * nftos;

        THoEX_ggF_A_tautau_CMS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mA);
        R_ggF_A_tautau_CMS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mA))/ip_ex_gg_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mA);
        R_bbF_A_tautau_CMS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mA))/ip_ex_bb_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS13=pp_A_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mA);
        R_pp_A_gaga_ATLAS13=(1+(pp_A_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mA))/ip_ex_pp_phi_gaga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS13=ggF_A_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mA);
        R_ggF_A_gaga_CMS13=(1+(ggF_A_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mA))/ip_ex_gg_phi_gaga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_llga_ATLAS13=ggF_A_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(mA);
        R_ggF_A_Zga_llga_ATLAS13=(1+(ggF_A_Zga_TH13-ip_ex_gg_phi_Zga_llga_ATLAS13(mA))/ip_ex_gg_phi_Zga_llga_ATLAS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_qqga_CMS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_qqga_CMS13(mA);
        R_pp_A_Zga_qqga_CMS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_qqga_CMS13(mA))/ip_ex_pp_phi_Zga_qqga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_CMS13=ggF_A_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mA);
        R_ggF_A_Zga_CMS13=(1+(ggF_A_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mA))/ip_ex_ggF_phi_Zga_CMS13_e(mA) ) * nftos;
    }
    else if(mA>=2500.0 && mA<2700.0)
    {
        THoEX_ggF_A_tt_ATLAS8=ggF_A_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mA);
        R_ggF_A_tt_ATLAS8=(1+(ggF_A_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mA))/ip_ex_gg_phi_tt_ATLAS8_e(mA) ) * nftos;

        THoEX_ggF_A_tautau_CMS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mA);
        R_ggF_A_tautau_CMS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mA))/ip_ex_gg_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mA);
        R_bbF_A_tautau_CMS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mA))/ip_ex_bb_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_gaga_ATLAS13=pp_A_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(mA);
        R_pp_A_gaga_ATLAS13=(1+(pp_A_gaga_TH13-ip_ex_pp_phi_gaga_ATLAS13(mA))/ip_ex_pp_phi_gaga_ATLAS13_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS13=ggF_A_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mA);
        R_ggF_A_gaga_CMS13=(1+(ggF_A_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mA))/ip_ex_gg_phi_gaga_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_qqga_CMS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_qqga_CMS13(mA);
        R_pp_A_Zga_qqga_CMS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_qqga_CMS13(mA))/ip_ex_pp_phi_Zga_qqga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_CMS13=ggF_A_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mA);
        R_ggF_A_Zga_CMS13=(1+(ggF_A_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mA))/ip_ex_ggF_phi_Zga_CMS13_e(mA) ) * nftos;
    }
    else if(mA>=2700.0 && mA<3000.0)
    {
        THoEX_ggF_A_tt_ATLAS8=ggF_A_tt_TH8/ip_ex_gg_phi_tt_ATLAS8(mA);
        R_ggF_A_tt_ATLAS8=(1+(ggF_A_tt_TH8-ip_ex_gg_phi_tt_ATLAS8(mA))/ip_ex_gg_phi_tt_ATLAS8_e(mA) ) * nftos;

        THoEX_ggF_A_tautau_CMS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mA);
        R_ggF_A_tautau_CMS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mA))/ip_ex_gg_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mA);
        R_bbF_A_tautau_CMS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mA))/ip_ex_bb_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS13=ggF_A_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mA);
        R_ggF_A_gaga_CMS13=(1+(ggF_A_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mA))/ip_ex_gg_phi_gaga_CMS13_e(mA) ) * nftos;
        THoEX_pp_A_Zga_qqga_CMS13=pp_A_Zga_TH13/ip_ex_pp_phi_Zga_qqga_CMS13(mA);
        R_pp_A_Zga_qqga_CMS13=(1+(pp_A_Zga_TH13-ip_ex_pp_phi_Zga_qqga_CMS13(mA))/ip_ex_pp_phi_Zga_qqga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_CMS13=ggF_A_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mA);
        R_ggF_A_Zga_CMS13=(1+(ggF_A_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mA))/ip_ex_ggF_phi_Zga_CMS13_e(mA) ) * nftos;
    }
    else if(mA>=3000.0 && mA<3200.0)
    {
        THoEX_ggF_A_tautau_CMS13=ggF_A_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(mA);
        R_ggF_A_tautau_CMS13=(1+(ggF_A_tautau_TH13-ip_ex_gg_phi_tautau_CMS13(mA))/ip_ex_gg_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_bbF_A_tautau_CMS13=bbF_A_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(mA);
        R_bbF_A_tautau_CMS13=(1+(bbF_A_tautau_TH13-ip_ex_bb_phi_tautau_CMS13(mA))/ip_ex_bb_phi_tautau_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_gaga_CMS13=ggF_A_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mA);
        R_ggF_A_gaga_CMS13=(1+(ggF_A_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mA))/ip_ex_gg_phi_gaga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_CMS13=ggF_A_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mA);
        R_ggF_A_Zga_CMS13=(1+(ggF_A_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mA))/ip_ex_ggF_phi_Zga_CMS13_e(mA) ) * nftos;
    }
    else if(mA>=3200.0 && mA<4000.0)
    {
        THoEX_ggF_A_gaga_CMS13=ggF_A_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(mA);
        R_ggF_A_gaga_CMS13=(1+(ggF_A_gaga_TH13-ip_ex_gg_phi_gaga_CMS13(mA))/ip_ex_gg_phi_gaga_CMS13_e(mA) ) * nftos;
        THoEX_ggF_A_Zga_CMS13=ggF_A_Zga_TH13/ip_ex_ggF_phi_Zga_CMS13(mA);
        R_ggF_A_Zga_CMS13=(1+(ggF_A_Zga_TH13-ip_ex_ggF_phi_Zga_CMS13(mA))/ip_ex_ggF_phi_Zga_CMS13_e(mA) ) * nftos;
    }
}

void THDMcache::computeHpquantities()
{
    double MW2=MW*MW;
    double mHh=sqrt(mHh2);
    double mA=sqrt(mA2);
    double mHp=sqrt(mHp2);
    double Mt2=Mt*Mt;
    double Mb2=Mb*Mb;
    double Mtau2=Mtau*Mtau;
    double Vtb=myTHDM->getCKM().getV_tb().abs();
    double gHpt=-1.0/tanb;
    double gHpb=0.0;
    double gHptau=0.0;
    double vev2=vev*vev;

    double SigmaHp8=0.0;
    double SigmaHp13=0.0;
    
    if( modelflag == "type1" ) {
        SigmaHp8=ip_cs_ggtoHp_8(mHp,0.0)/(tanb*tanb);
        SigmaHp13=ip_cs_ggtoHp_13(mHp,0.0)/(tanb*tanb);
        gHpb=1.0/tanb;
        gHptau=1.0/tanb;
    }
    else if( modelflag == "type2" ) {
        if(logtb>=-1.0 && logtb<=1.75)
        {
            SigmaHp8=ip_cs_ggtoHp_8(mHp,logtb);
            SigmaHp13=ip_cs_ggtoHp_13(mHp,logtb);
        }
        gHpb=-tanb;
        gHptau=-tanb;
    }
    else if( modelflag == "typeX" ) {
        SigmaHp8=ip_cs_ggtoHp_8(mHp,0.0)/(tanb*tanb);
        SigmaHp13=ip_cs_ggtoHp_13(mHp,0.0)/(tanb*tanb);
        gHpb=1.0/tanb;
        gHptau=-tanb;
    }
    else if( modelflag == "typeY" ) {
        if(logtb>=-1.0 && logtb<=1.75)
        {
            SigmaHp8=ip_cs_ggtoHp_8(mHp,logtb);
            SigmaHp13=ip_cs_ggtoHp_13(mHp,logtb);
        }
        gHpb=-tanb;
        gHptau=1.0/tanb;
    }
    else {
        throw std::runtime_error("modelflag can be only any of \"type1\", \"type2\", \"typeX\" or \"typeY\"");
    }

    double GammaHptaunu=HSTheta(mHp-Mtau)*(Mtau2*(mHp2-Mtau2)*(mHp2-Mtau2)*gHptau*gHptau)/(8.0*mHp*mHp2*M_PI*vev2);
    double GammaHptb=HSTheta(mHp-Mt-Mb)*(Vtb*Vtb/(8.0*mHp*M_PI*vev2))*3.0*(-4.0*gHpb*gHpt*Mb2*Mt2
                        -gHpb*gHpb*Mb2*(Mb2-mHp2+Mt2)-gHpt*gHpt*Mt2*(Mb2-mHp2+Mt2))
                      *sqrt((Mb2*Mb2+(mHp2-Mt2)*(mHp2-Mt2)-2.0*Mb2*(mHp2+Mt2))/(mHp2*mHp2));
    double GammaHpHlW=KaellenFunction(1.0,mHl/mHp,MW/mHp)*KaellenFunction(1.0,mHp/MW,mHl/MW)*KaellenFunction(1.0,mHp/MW,mHl/MW)
                      *MW2*MW2/mHp*cos(bma)*cos(bma)/(2.0*M_PI*vev2);
    double GammaHpHhW=KaellenFunction(1.0,mHh/mHp,MW/mHp)*KaellenFunction(1.0,mHp/MW,mHh/MW)*KaellenFunction(1.0,mHp/MW,mHh/MW)
                      *MW2*MW2/mHp*sin(bma)*sin(bma)/(2.0*M_PI*vev2);
    double GammaHpAW=KaellenFunction(1.0,mA/mHp,MW/mHp)*KaellenFunction(1.0,mHp/MW,mA/MW)*KaellenFunction(1.0,mHp/MW,mA/MW)
                      *MW2*MW2/mHp/(2.0*M_PI*vev2);
    GammaHptot= GammaHptaunu + GammaHptb + GammaHpHlW + GammaHpHhW + GammaHpAW;

    double Br_Hptotaunu=GammaHptaunu/GammaHptot;
    double Br_Hptotb=GammaHptb/GammaHptot;

    //Theoretical expressions for the charged Higgs cross sections times branching ratios

    pp_Hpm_taunu_TH8=2.0*SigmaHp8*Br_Hptotaunu;
    pp_Hp_taunu_TH8=SigmaHp8*Br_Hptotaunu;
    pp_Hpm_tb_TH8=2.0*SigmaHp8*Br_Hptotb;
    pp_Hp_tb_TH8=SigmaHp8*Br_Hptotb;
    pp_Hpm_taunu_TH13=2.0*SigmaHp13*Br_Hptotaunu;
    pp_Hp_tb_TH13=SigmaHp13*Br_Hptotb;

    //Ratios of theoretical Heavy Higgs cross sections and experimental upper limits

    THoEX_pp_Hpm_taunu_ATLAS8=0.0;
    R_pp_Hpm_taunu_ATLAS8=0.0;
    THoEX_pp_Hp_taunu_CMS8=0.0;
    R_pp_Hp_taunu_CMS8=0.0;
    THoEX_pp_Hpm_tb_ATLAS8=0.0;
    R_pp_Hpm_tb_ATLAS8=0.0;
    THoEX_pp_Hp_tb_CMS8=0.0;
    R_pp_Hp_tb_CMS8=0.0;
    THoEX_pp_Hpm_taunu_ATLAS13=0.0;
    R_pp_Hpm_taunu_ATLAS13=0.0;
    THoEX_pp_Hpm_taunu_CMS13=0.0;
    R_pp_Hpm_taunu_CMS13=0.0;
    THoEX_pp_Hp_tb_ATLAS13_1=0.0;
    R_pp_Hp_tb_ATLAS13_1=0.0;
    THoEX_pp_Hp_tb_ATLAS13_2=0.0;
    R_pp_Hp_tb_ATLAS13_2=0.0;
    THoEX_pp_Hp_tb_ATLAS13=0.0;
    R_pp_Hp_tb_ATLAS13=0.0;

    //95% to 1 sigma conversion factor, roughly sqrt(3.84)
    double nftos=1.95996398454;

    if(mHp>=180.0 && mHp<200.0)
    {
        THoEX_pp_Hpm_taunu_ATLAS8=pp_Hpm_taunu_TH8/ip_ex_pp_Hpm_taunu_ATLAS8(mHp);
        R_pp_Hpm_taunu_ATLAS8=(1+(pp_Hpm_taunu_TH8-ip_ex_pp_Hpm_taunu_ATLAS8(mHp))/ip_ex_pp_Hpm_taunu_ATLAS8_e(mHp) ) * nftos;
        THoEX_pp_Hp_taunu_CMS8=pp_Hp_taunu_TH8/ip_ex_pp_Hp_taunu_CMS8(mHp);
        R_pp_Hp_taunu_CMS8=(1+(pp_Hp_taunu_TH8-ip_ex_pp_Hp_taunu_CMS8(mHp))/ip_ex_pp_Hp_taunu_CMS8_e(mHp) ) * nftos;
        THoEX_pp_Hp_tb_CMS8=pp_Hp_tb_TH8/ip_ex_pp_Hp_tb_CMS8(mHp);
        R_pp_Hp_tb_CMS8=(1+(pp_Hp_tb_TH8-ip_ex_pp_Hp_tb_CMS8(mHp))/ip_ex_pp_Hp_tb_CMS8_e(mHp) ) * nftos;

        THoEX_pp_Hpm_taunu_CMS13=pp_Hpm_taunu_TH13/ip_ex_pp_Hpm_taunu_CMS13(mHp);
        R_pp_Hpm_taunu_CMS13=(1+(pp_Hpm_taunu_TH13-ip_ex_pp_Hpm_taunu_CMS13(mHp))/ip_ex_pp_Hpm_taunu_CMS13_e(mHp) ) * nftos;
    }
    else if(mHp>=200.0 && mHp<300.0)
    {
        THoEX_pp_Hpm_taunu_ATLAS8=pp_Hpm_taunu_TH8/ip_ex_pp_Hpm_taunu_ATLAS8(mHp);
        R_pp_Hpm_taunu_ATLAS8=(1+(pp_Hpm_taunu_TH8-ip_ex_pp_Hpm_taunu_ATLAS8(mHp))/ip_ex_pp_Hpm_taunu_ATLAS8_e(mHp) ) * nftos;
        THoEX_pp_Hp_taunu_CMS8=pp_Hp_taunu_TH8/ip_ex_pp_Hp_taunu_CMS8(mHp);
        R_pp_Hp_taunu_CMS8=(1+(pp_Hp_taunu_TH8-ip_ex_pp_Hp_taunu_CMS8(mHp))/ip_ex_pp_Hp_taunu_CMS8_e(mHp) ) * nftos;
        THoEX_pp_Hpm_tb_ATLAS8=pp_Hpm_tb_TH8/ip_ex_pp_Hpm_tb_ATLAS8(mHp);
        R_pp_Hpm_tb_ATLAS8=(1+(pp_Hpm_tb_TH8-ip_ex_pp_Hpm_tb_ATLAS8(mHp))/ip_ex_pp_Hpm_tb_ATLAS8_e(mHp) ) * nftos;
        THoEX_pp_Hp_tb_CMS8=pp_Hp_tb_TH8/ip_ex_pp_Hp_tb_CMS8(mHp);
        R_pp_Hp_tb_CMS8=(1+(pp_Hp_tb_TH8-ip_ex_pp_Hp_tb_CMS8(mHp))/ip_ex_pp_Hp_tb_CMS8_e(mHp) ) * nftos;

        THoEX_pp_Hpm_taunu_ATLAS13=pp_Hpm_taunu_TH13/ip_ex_pp_Hpm_taunu_ATLAS13(mHp);
        R_pp_Hpm_taunu_ATLAS13=(1+(pp_Hpm_taunu_TH13-ip_ex_pp_Hpm_taunu_ATLAS13(mHp))/ip_ex_pp_Hpm_taunu_ATLAS13_e(mHp) ) * nftos;
        THoEX_pp_Hpm_taunu_CMS13=pp_Hpm_taunu_TH13/ip_ex_pp_Hpm_taunu_CMS13(mHp);
        R_pp_Hpm_taunu_CMS13=(1+(pp_Hpm_taunu_TH13-ip_ex_pp_Hpm_taunu_CMS13(mHp))/ip_ex_pp_Hpm_taunu_CMS13_e(mHp) ) * nftos;
        THoEX_pp_Hp_tb_ATLAS13_2=pp_Hp_tb_TH13/ip_ex_pp_Hp_tb_ATLAS13_2(mHp);
        R_pp_Hp_tb_ATLAS13_2=(1+(pp_Hp_tb_TH13-ip_ex_pp_Hp_tb_ATLAS13_2(mHp))/ip_ex_pp_Hp_tb_ATLAS13_2_e(mHp) ) * nftos;
        THoEX_pp_Hp_tb_ATLAS13=pp_Hp_tb_TH13/ip_ex_pp_Hp_tb_ATLAS13_2(mHp);
        R_pp_Hp_tb_ATLAS13=(1+(pp_Hp_tb_TH13-ip_ex_pp_Hp_tb_ATLAS13_2(mHp))/ip_ex_pp_Hp_tb_ATLAS13_2_e(mHp) ) * nftos;
    }
    else if(mHp>=300.0 && mHp<600.0)
    {
        THoEX_pp_Hpm_taunu_ATLAS8=pp_Hpm_taunu_TH8/ip_ex_pp_Hpm_taunu_ATLAS8(mHp);
        R_pp_Hpm_taunu_ATLAS8=(1+(pp_Hpm_taunu_TH8-ip_ex_pp_Hpm_taunu_ATLAS8(mHp))/ip_ex_pp_Hpm_taunu_ATLAS8_e(mHp) ) * nftos;
        THoEX_pp_Hp_taunu_CMS8=pp_Hp_taunu_TH8/ip_ex_pp_Hp_taunu_CMS8(mHp);
        R_pp_Hp_taunu_CMS8=(1+(pp_Hp_taunu_TH8-ip_ex_pp_Hp_taunu_CMS8(mHp))/ip_ex_pp_Hp_taunu_CMS8_e(mHp) ) * nftos;
        THoEX_pp_Hpm_tb_ATLAS8=pp_Hpm_tb_TH8/ip_ex_pp_Hpm_tb_ATLAS8(mHp);
        R_pp_Hpm_tb_ATLAS8=(1+(pp_Hpm_tb_TH8-ip_ex_pp_Hpm_tb_ATLAS8(mHp))/ip_ex_pp_Hpm_tb_ATLAS8_e(mHp) ) * nftos;
        THoEX_pp_Hp_tb_CMS8=pp_Hp_tb_TH8/ip_ex_pp_Hp_tb_CMS8(mHp);
        R_pp_Hp_tb_CMS8=(1+(pp_Hp_tb_TH8-ip_ex_pp_Hp_tb_CMS8(mHp))/ip_ex_pp_Hp_tb_CMS8_e(mHp) ) * nftos;

        THoEX_pp_Hpm_taunu_ATLAS13=pp_Hpm_taunu_TH13/ip_ex_pp_Hpm_taunu_ATLAS13(mHp);
        R_pp_Hpm_taunu_ATLAS13=(1+(pp_Hpm_taunu_TH13-ip_ex_pp_Hpm_taunu_ATLAS13(mHp))/ip_ex_pp_Hpm_taunu_ATLAS13_e(mHp) ) * nftos;
        THoEX_pp_Hpm_taunu_CMS13=pp_Hpm_taunu_TH13/ip_ex_pp_Hpm_taunu_CMS13(mHp);
        R_pp_Hpm_taunu_CMS13=(1+(pp_Hpm_taunu_TH13-ip_ex_pp_Hpm_taunu_CMS13(mHp))/ip_ex_pp_Hpm_taunu_CMS13_e(mHp) ) * nftos;
        THoEX_pp_Hp_tb_ATLAS13_1=pp_Hp_tb_TH13/ip_ex_pp_Hp_tb_ATLAS13_1(mHp);
        R_pp_Hp_tb_ATLAS13_1=(1+(pp_Hp_tb_TH13-ip_ex_pp_Hp_tb_ATLAS13_1(mHp))/ip_ex_pp_Hp_tb_ATLAS13_1_e(mHp) ) * nftos;
        THoEX_pp_Hp_tb_ATLAS13_2=pp_Hp_tb_TH13/ip_ex_pp_Hp_tb_ATLAS13_2(mHp);
        R_pp_Hp_tb_ATLAS13_2=(1+(pp_Hp_tb_TH13-ip_ex_pp_Hp_tb_ATLAS13_2(mHp))/ip_ex_pp_Hp_tb_ATLAS13_2_e(mHp) ) * nftos;
        THoEX_pp_Hp_tb_ATLAS13=pp_Hp_tb_TH13/ip_ex_pp_Hp_tb_ATLAS13_1(mHp);
        R_pp_Hp_tb_ATLAS13=(1+(pp_Hp_tb_TH13-ip_ex_pp_Hp_tb_ATLAS13_1(mHp))/ip_ex_pp_Hp_tb_ATLAS13_1_e(mHp) ) * nftos;
    }
    else if(mHp>=600.0 && mHp<1000.0)
    {
        THoEX_pp_Hpm_taunu_ATLAS8=pp_Hpm_taunu_TH8/ip_ex_pp_Hpm_taunu_ATLAS8(mHp);
        R_pp_Hpm_taunu_ATLAS8=(1+(pp_Hpm_taunu_TH8-ip_ex_pp_Hpm_taunu_ATLAS8(mHp))/ip_ex_pp_Hpm_taunu_ATLAS8_e(mHp) ) * nftos;

        THoEX_pp_Hpm_taunu_ATLAS13=pp_Hpm_taunu_TH13/ip_ex_pp_Hpm_taunu_ATLAS13(mHp);
        R_pp_Hpm_taunu_ATLAS13=(1+(pp_Hpm_taunu_TH13-ip_ex_pp_Hpm_taunu_ATLAS13(mHp))/ip_ex_pp_Hpm_taunu_ATLAS13_e(mHp) ) * nftos;
        THoEX_pp_Hpm_taunu_CMS13=pp_Hpm_taunu_TH13/ip_ex_pp_Hpm_taunu_CMS13(mHp);
        R_pp_Hpm_taunu_CMS13=(1+(pp_Hpm_taunu_TH13-ip_ex_pp_Hpm_taunu_CMS13(mHp))/ip_ex_pp_Hpm_taunu_CMS13_e(mHp) ) * nftos;
        THoEX_pp_Hp_tb_ATLAS13_1=pp_Hp_tb_TH13/ip_ex_pp_Hp_tb_ATLAS13_1(mHp);
        R_pp_Hp_tb_ATLAS13_1=(1+(pp_Hp_tb_TH13-ip_ex_pp_Hp_tb_ATLAS13_1(mHp))/ip_ex_pp_Hp_tb_ATLAS13_1_e(mHp) ) * nftos;
        THoEX_pp_Hp_tb_ATLAS13_2=pp_Hp_tb_TH13/ip_ex_pp_Hp_tb_ATLAS13_2(mHp);
        R_pp_Hp_tb_ATLAS13_2=(1+(pp_Hp_tb_TH13-ip_ex_pp_Hp_tb_ATLAS13_2(mHp))/ip_ex_pp_Hp_tb_ATLAS13_2_e(mHp) ) * nftos;
        THoEX_pp_Hp_tb_ATLAS13=pp_Hp_tb_TH13/ip_ex_pp_Hp_tb_ATLAS13_1(mHp);
        R_pp_Hp_tb_ATLAS13=(1+(pp_Hp_tb_TH13-ip_ex_pp_Hp_tb_ATLAS13_1(mHp))/ip_ex_pp_Hp_tb_ATLAS13_1_e(mHp) ) * nftos;
    }
    else if(mHp>=1000.0 && mHp<1400.0)
    {
        THoEX_pp_Hpm_taunu_ATLAS13=pp_Hpm_taunu_TH13/ip_ex_pp_Hpm_taunu_ATLAS13(mHp);
        R_pp_Hpm_taunu_ATLAS13=(1+(pp_Hpm_taunu_TH13-ip_ex_pp_Hpm_taunu_ATLAS13(mHp))/ip_ex_pp_Hpm_taunu_ATLAS13_e(mHp) ) * nftos;
        THoEX_pp_Hpm_taunu_CMS13=pp_Hpm_taunu_TH13/ip_ex_pp_Hpm_taunu_CMS13(mHp);
        R_pp_Hpm_taunu_CMS13=(1+(pp_Hpm_taunu_TH13-ip_ex_pp_Hpm_taunu_CMS13(mHp))/ip_ex_pp_Hpm_taunu_CMS13_e(mHp) ) * nftos;
        THoEX_pp_Hp_tb_ATLAS13_2=pp_Hp_tb_TH13/ip_ex_pp_Hp_tb_ATLAS13_2(mHp);
        R_pp_Hp_tb_ATLAS13_2=(1+(pp_Hp_tb_TH13-ip_ex_pp_Hp_tb_ATLAS13_2(mHp))/ip_ex_pp_Hp_tb_ATLAS13_2_e(mHp) ) * nftos;
        THoEX_pp_Hp_tb_ATLAS13=pp_Hp_tb_TH13/ip_ex_pp_Hp_tb_ATLAS13_2(mHp);
        R_pp_Hp_tb_ATLAS13=(1+(pp_Hp_tb_TH13-ip_ex_pp_Hp_tb_ATLAS13_2(mHp))/ip_ex_pp_Hp_tb_ATLAS13_2_e(mHp) ) * nftos;
    }
    else if(mHp>=1400.0 && mHp<2000.0)
    {
        THoEX_pp_Hpm_taunu_ATLAS13=pp_Hpm_taunu_TH13/ip_ex_pp_Hpm_taunu_ATLAS13(mHp);
        R_pp_Hpm_taunu_ATLAS13=(1+(pp_Hpm_taunu_TH13-ip_ex_pp_Hpm_taunu_ATLAS13(mHp))/ip_ex_pp_Hpm_taunu_ATLAS13_e(mHp) ) * nftos;
        THoEX_pp_Hpm_taunu_CMS13=pp_Hpm_taunu_TH13/ip_ex_pp_Hpm_taunu_CMS13(mHp);
        R_pp_Hpm_taunu_CMS13=(1+(pp_Hpm_taunu_TH13-ip_ex_pp_Hpm_taunu_CMS13(mHp))/ip_ex_pp_Hpm_taunu_CMS13_e(mHp) ) * nftos;
        THoEX_pp_Hp_tb_ATLAS13_2=pp_Hp_tb_TH13/ip_ex_pp_Hp_tb_ATLAS13_2(mHp);
        R_pp_Hp_tb_ATLAS13_2=(1+(pp_Hp_tb_TH13-ip_ex_pp_Hp_tb_ATLAS13_2(mHp))/ip_ex_pp_Hp_tb_ATLAS13_2_e(mHp) ) * nftos;
        THoEX_pp_Hp_tb_ATLAS13=pp_Hp_tb_TH13/ip_ex_pp_Hp_tb_ATLAS13_2(mHp);
        R_pp_Hp_tb_ATLAS13=(1+(pp_Hp_tb_TH13-ip_ex_pp_Hp_tb_ATLAS13_2(mHp))/ip_ex_pp_Hp_tb_ATLAS13_2_e(mHp) ) * nftos;
    }
}

void THDMcache::runTHDMparameters()
{
    vev=myTHDM->v();
    double cosb=myTHDM->getcosb();
    double sinb=myTHDM->getsinb();
    modelflag=myTHDM->getModelTypeflag();
    std::string RGEorder=myTHDM->getRGEorderflag();
    //flag will be used to transport information about model and RGEorder to the Runner:
    //flag=3*(0 for type I, 1 for type II, 2 for type X and 3 for type Y) + (0 for LO, 1 for approxNLO and 2 for NLO)
    int flag;
    if( RGEorder == "LO" ) flag=0;
    else if( RGEorder == "approxNLO" ) flag=1;
    else if( RGEorder == "NLO" ) flag=2;
    else {
        throw std::runtime_error("RGEorder can be only any of \"LO\", \"approxNLO\" or \"NLO\"");
    }

    double g1_at_MZ=sqrt(4.0*M_PI*Ale/cW2);
    double g2_at_MZ=sqrt(4.0*M_PI*Ale/(1-cW2));
    double g3_at_MZ=sqrt(4.0*M_PI*Als);

    double Ytop_at_MZ=(sqrt(2.0)*myTHDM->getQuarks(QCD::TOP).getMass())/(vev*sinb);
    double Ybottom1_at_MZ=0.0;
    double Ybottom2_at_MZ=0.0;
    double Ytau1_at_MZ=0.0;
    double Ytau2_at_MZ=0.0;

    /*link these to the SM values*/
    double Mb_at_MZ=2.96;//GeV
    double Mtau_at_MZ=1.75;//GeV

    if( modelflag == "type1" ) {
        Ybottom2_at_MZ=(sqrt(2.0)*Mb_at_MZ)/(vev*sinb);
        Ytau2_at_MZ=(sqrt(2.0)*Mtau_at_MZ)/(vev*sinb);
    }
    else if( modelflag == "type2" ) {
        Ybottom1_at_MZ=(sqrt(2.0)*Mb_at_MZ)/(vev*cosb);
        Ytau1_at_MZ=(sqrt(2.0)*Mtau_at_MZ)/(vev*cosb);
        flag +=3;
    }
    else if( modelflag == "typeX" ) {
        Ybottom2_at_MZ=(sqrt(2.0)*Mb_at_MZ)/(vev*sinb);
        Ytau1_at_MZ=(sqrt(2.0)*Mtau_at_MZ)/(vev*cosb);
        flag +=6;
    }
    else if( modelflag == "typeY" ) {
        Ybottom1_at_MZ=(sqrt(2.0)*Mb_at_MZ)/(vev*cosb);
        Ytau2_at_MZ=(sqrt(2.0)*Mtau_at_MZ)/(vev*sinb);
        flag +=9;
    }
    else {
        throw std::runtime_error("modelflag can be only any of \"type1\", \"type2\", \"typeX\" or \"typeY\"");
    }

    double m11_2_at_MZ=mym11_2->computeThValue();
    double m22_2_at_MZ=mym22_2->computeThValue();
    double m12_2_at_MZ=myTHDM->getm12_2();
    double lambda1_at_MZ=mylambda1->computeThValue();
    double lambda2_at_MZ=mylambda2->computeThValue();
    double lambda3_at_MZ=mylambda3->computeThValue();
    double lambda4_at_MZ=mylambda4->computeThValue();
    double lambda5_at_MZ=mylambda5->computeThValue();

    double Rpeps=myTHDM->getRpeps();
    double NLOuniscale=myTHDM->getNLOuniscale();

    if(fabs(Q_THDM-log10(MZ))<0.005)   //at MZ scale
    {
        Q_cutoff=log10(MZ);

        g1_at_Q = g1_at_MZ;
        g2_at_Q = g2_at_MZ;
        g3_at_Q = g3_at_MZ;
        Ytop_at_Q = Ytop_at_MZ;
        Ybottom1_at_Q = Ybottom1_at_MZ;
        Ybottom2_at_Q = Ybottom2_at_MZ;
        Ytau1_at_Q = Ytau1_at_MZ;
        Ytau2_at_Q = Ytau2_at_MZ;
        m11_2_at_Q = m11_2_at_MZ;
        m22_2_at_Q = m22_2_at_MZ;
        m12_2_at_Q = m12_2_at_MZ;
        lambda1_at_Q = lambda1_at_MZ;
        lambda2_at_Q = lambda2_at_MZ;
        lambda3_at_Q = lambda3_at_MZ;
        lambda4_at_Q = lambda4_at_MZ;
        lambda5_at_Q = lambda5_at_MZ;
    }
    else   //at some other scale
    {
        double InitVals[14];
        InitVals[0]=g1_at_MZ;
        InitVals[1]=g2_at_MZ;
        InitVals[2]=g3_at_MZ;
        InitVals[3]=Ytop_at_MZ;
        InitVals[4]=Ybottom1_at_MZ+Ybottom2_at_MZ;
        InitVals[5]=Ytau1_at_MZ+Ytau2_at_MZ;
        InitVals[6]=m11_2_at_MZ;
        InitVals[7]=m22_2_at_MZ;
        InitVals[8]=m12_2_at_MZ;
        InitVals[9]=lambda1_at_MZ;
        InitVals[10]=lambda2_at_MZ;
        InitVals[11]=lambda3_at_MZ;
        InitVals[12]=lambda4_at_MZ;
        InitVals[13]=lambda5_at_MZ;

        Q_cutoff=myRunner->RGERunner(InitVals, 14, log10(MZ), Q_THDM, flag, Rpeps, NLOuniscale);  //Running up to Q_cutoff<=Q_THDM

        g1_at_Q = InitVals[0];
        g2_at_Q = InitVals[1];
        g3_at_Q = InitVals[2];
        Ytop_at_Q = InitVals[3];
        Ybottom1_at_Q = 0.0;
        Ybottom2_at_Q = 0.0;
        Ytau1_at_Q = 0.0;
        Ytau2_at_Q = 0.0;
        if( modelflag == "type1" ) {
            Ybottom2_at_Q=InitVals[4];
            Ytau2_at_Q=InitVals[5];
        }
        else if( modelflag == "type2" ) {
            Ybottom1_at_Q=InitVals[4];
            Ytau1_at_Q=InitVals[5];
        }
        else if( modelflag == "typeX" ) {
            Ybottom2_at_Q=InitVals[4];
            Ytau1_at_Q=InitVals[5];
        }
        else if( modelflag == "typeY" ) {
            Ybottom1_at_Q=InitVals[4];
            Ytau2_at_Q=InitVals[5];
        }
        else {
            throw std::runtime_error("modelflag can be only any of \"type1\", \"type2\", \"typeX\" or \"typeY\"");
        }
        m11_2_at_Q = InitVals[6];
        m22_2_at_Q = InitVals[7];
        m12_2_at_Q = InitVals[8];
        lambda1_at_Q = InitVals[9];
        lambda2_at_Q = InitVals[10];
        lambda3_at_Q = InitVals[11];
        lambda4_at_Q = InitVals[12];
        lambda5_at_Q = InitVals[13];
    }

}

void THDMcache::computeWFRcombinations()
{
    double WFRcomb1a = 0.0;
    double WFRcomb1b = 0.0;
    double WFRcomb2a = 0.0;
    double WFRcomb3a = 0.0;
    double WFRcomb3b = 0.0;
    double WFRcomb4a = 0.0;

    if(WFRflag)
    {
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double beta=atan(tanb);
    double alpha=beta-bma;
    double MZ2=MZ*MZ;

    double B000mh = B0_MZ2_0_0_mHl2(MZ2,mHl2).real();
    double B000mH = B0_MZ2_0_0_mHh2(MZ2,mHh2).real();
    double B00mHpmh = B0_MZ2_0_mHp2_mHl2(MZ2,mHp2,mHl2).real();
    double B00mHpmH = B0_MZ2_0_mHp2_mHh2(MZ2,mHp2,mHh2).real();
    double B00mAmh = B0_MZ2_0_mA2_mHl2(MZ2,mA2,mHl2).real();
    double B00mAmH = B0_MZ2_0_mA2_mHh2(MZ2,mA2,mHh2).real();
    double B0mh00 = B0_MZ2_mHl2_0_0(MZ2,mHl2).real();
    double B0mh0mHp = B0_MZ2_mHl2_0_mHp2(MZ2,mHl2,mHp2).real();
    double B0mh0mA = B0_MZ2_mHl2_0_mA2(MZ2,mHl2,mA2).real();
    double B0mhmhmh = B0_MZ2_mHl2_mHl2_mHl2(MZ2,mHl2).real();
    double B0mhmHmh = B0_MZ2_mHl2_mHh2_mHl2(MZ2,mHl2,mHh2).real();
    double B0mhmHmH = B0_MZ2_mHl2_mHh2_mHh2(MZ2,mHl2,mHh2).real();
    double B0mhmHpmHp = B0_MZ2_mHl2_mHp2_mHp2(MZ2,mHl2,mHp2).real();
    double B0mhmAmA = B0_MZ2_mHl2_mA2_mA2(MZ2,mHl2,mA2).real();
    double B0mH00 = B0_MZ2_mHh2_0_0(MZ2,mHh2).real();
    double B0mH0mHp = B0_MZ2_mHh2_0_mHp2(MZ2,mHh2,mHp2).real();
    double B0mH0mA = B0_MZ2_mHh2_0_mA2(MZ2,mHh2,mA2).real();
    double B0mHmhmh = B0_MZ2_mHh2_mHl2_mHl2(MZ2,mHh2,mHl2).real();
    double B0mHmHmh = B0_MZ2_mHh2_mHh2_mHl2(MZ2,mHh2,mHl2).real();
    double B0mHmHmH = B0_MZ2_mHh2_mHh2_mHh2(MZ2,mHh2).real();
    double B0mHmHpmHp = B0_MZ2_mHh2_mHp2_mHp2(MZ2,mHh2,mHp2).real();
    double B0mHmAmA = B0_MZ2_mHh2_mA2_mA2(MZ2,mHh2,mA2).real();
    double B0mHp0mh = B0_MZ2_mHp2_0_mHl2(MZ2,mHp2,mHl2).real();
    double B0mHp0mH = B0_MZ2_mHp2_0_mHh2(MZ2,mHp2,mHh2).real();
    double B0mHpmHpmh = B0_MZ2_mHp2_mHp2_mHl2(MZ2,mHp2,mHl2).real();
    double B0mHpmHpmH = B0_MZ2_mHp2_mHp2_mHh2(MZ2,mHp2,mHh2).real();
    double B0mA0mh = B0_MZ2_mA2_0_mHl2(MZ2,mA2,mHl2).real();
    double B0mA0mH = B0_MZ2_mA2_0_mHh2(MZ2,mA2,mHh2).real();
    double B0mAmAmh = B0_MZ2_mA2_mA2_mHl2(MZ2,mA2,mHl2).real();
    double B0mAmAmH = B0_MZ2_mA2_mA2_mHh2(MZ2,mA2,mHh2).real();

    double ddpB000mh = B0p_MZ2_0_0_mHl2(MZ2,mHl2).real();
    double ddpB000mH = B0p_MZ2_0_0_mHh2(MZ2,mHh2).real();
    double ddpB00mHpmh = B0p_MZ2_0_mHp2_mHl2(MZ2,mHp2,mHl2).real();
    double ddpB00mHpmH = B0p_MZ2_0_mHp2_mHh2(MZ2,mHp2,mHh2).real();
    double ddpB00mHpmA = B0p_MZ2_0_mHp2_mA2(MZ2,mHp2,mA2).real();
    double ddpB00mAmh = B0p_MZ2_0_mA2_mHl2(MZ2,mA2,mHl2).real();
    double ddpB00mAmH = B0p_MZ2_0_mA2_mHh2(MZ2,mA2,mHh2).real();
    double ddpB0mh00 = B0p_MZ2_mHl2_0_0(MZ2,mHl2).real();
    double ddpB0mh0mHp = B0p_MZ2_mHl2_0_mHp2(MZ2,mHl2,mHp2).real();
    double ddpB0mh0mA = B0p_MZ2_mHl2_0_mA2(MZ2,mHl2,mA2).real();
    double ddpB0mhmhmh = B0p_MZ2_mHl2_mHl2_mHl2(MZ2,mHl2).real();
    double ddpB0mhmHmh = B0p_MZ2_mHl2_mHh2_mHl2(MZ2,mHl2,mHh2).real();
    double ddpB0mhmHmH = B0p_MZ2_mHl2_mHh2_mHh2(MZ2,mHl2,mHh2).real();
    double ddpB0mhmHpmHp = B0p_MZ2_mHl2_mHp2_mHp2(MZ2,mHl2,mHp2).real();
    double ddpB0mhmAmA = B0p_MZ2_mHl2_mA2_mA2(MZ2,mHl2,mA2).real();
    double ddpB0mH00 = B0p_MZ2_mHh2_0_0(MZ2,mHh2).real();
    double ddpB0mH0mHp = B0p_MZ2_mHh2_0_mHp2(MZ2,mHh2,mHp2).real();
    double ddpB0mH0mA = B0p_MZ2_mHh2_0_mA2(MZ2,mHh2,mA2).real();
    double ddpB0mHmhmh = B0p_MZ2_mHh2_mHl2_mHl2(MZ2,mHh2,mHl2).real();
    double ddpB0mHmHmh = B0p_MZ2_mHh2_mHh2_mHl2(MZ2,mHh2,mHl2).real();
    double ddpB0mHmHmH = B0p_MZ2_mHh2_mHh2_mHh2(MZ2,mHh2).real();
    double ddpB0mHmHpmHp = B0p_MZ2_mHh2_mHp2_mHp2(MZ2,mHh2,mHp2).real();
    double ddpB0mHmAmA = B0p_MZ2_mHh2_mA2_mA2(MZ2,mHh2,mA2).real();
    double ddpB0mHp0mh = B0p_MZ2_mHp2_0_mHl2(MZ2,mHp2,mHl2).real();
    double ddpB0mHp0mH = B0p_MZ2_mHp2_0_mHh2(MZ2,mHp2,mHh2).real();
    double ddpB0mHp0mA = B0p_MZ2_mHp2_0_mA2(MZ2,mHp2,mA2).real();
    double ddpB0mHpmHpmh = B0p_MZ2_mHp2_mHp2_mHl2(MZ2,mHp2,mHl2).real();
    double ddpB0mHpmHpmH = B0p_MZ2_mHp2_mHp2_mHh2(MZ2,mHp2,mHh2).real();
    double ddpB0mA0mh = B0p_MZ2_mA2_0_mHl2(MZ2,mA2,mHl2).real();
    double ddpB0mA0mH = B0p_MZ2_mA2_0_mHh2(MZ2,mA2,mHh2).real();
    double ddpB0mA0mHp = B0p_MZ2_mA2_0_mHp2(MZ2,mA2,mHp2).real();
    double ddpB0mAmAmh = B0p_MZ2_mA2_mA2_mHl2(MZ2,mA2,mHl2).real();
    double ddpB0mAmAmH = B0p_MZ2_mA2_mA2_mHh2(MZ2,mA2,mHh2).real();

    WFRcomb1a = 3.0*mHl2*mHl2*cosb*cosb*sin(bma)*sin(bma) * ddpB000mh
        + 3.0*mHh2*mHh2*cos(bma)*cos(bma)*cosb*cosb * ddpB000mH
        + 2.0*(mHl2-mHp2)*(mHl2-mHp2)*cos(bma)*cos(bma)*cosb*cosb * ddpB00mHpmh
        + 2.0*(mHh2-mHp2)*(mHh2-mHp2)*sin(bma)*sin(bma)*cosb*cosb * ddpB00mHpmH
        + 2.0*(mA2-mHp2)*(mA2-mHp2)*cosb*cosb * ddpB00mHpmA
        + (mA2-mHl2)*(mA2-mHl2)*cos(bma)*cos(bma)*cosb*cosb * ddpB00mAmh
        + (mA2-mHh2)*(mA2-mHh2)*cosb*cosb*sin(bma)*sin(bma) * ddpB00mAmH
        + 1.5*mHl2*mHl2*sin(alpha)*sin(alpha)*sin(bma)*sin(bma) * ddpB0mh00
        + 2.0*(mHl2-mHp2)*(mHl2-mHp2)*cos(bma)*cos(bma)*sin(alpha)*sin(alpha) * ddpB0mh0mHp
        + (mA2-mHl2)*(mA2-mHl2)*cos(bma)*cos(bma)*sin(alpha)*sin(alpha) * ddpB0mh0mA
        + 9.0*sin(alpha)*sin(alpha)*pow(-mHl2*(3.0*sin(bma)+sin(3.0*bma)+sin(3.0*alpha+beta)+3.0*sin(alpha+3.0*beta))
                                        +16.0*m12_2*cos(bma)*cos(bma)*cos(alpha+beta),2)/(512.0*pow(cosb*sinb,4)) * ddpB0mhmhmh
        + sin(alpha)*sin(alpha)*pow((cos(alpha)/sinb+sin(alpha)/cosb)*(m12_2+cos(alpha)*sin(alpha)*(mHh2+2.0*mHl2-(3.0*m12_2)/(cosb*sinb))),2) * ddpB0mhmHmh
        + sin(alpha)*sin(alpha)*sin(bma)*sin(bma)*pow(-2.0*m12_2+(2.0*mHh2+mHl2-(3.0*m12_2)/(cosb*sinb))*sin(2.0*alpha),2)/(8.0*cosb*cosb*sinb*sinb) * ddpB0mhmHmH
        + (sin(alpha)*sin(alpha)*pow((mHl2-2.0*mHp2)*cos(alpha-3.0*beta)*2.0*sinb*cosb
                                     +cos(alpha+beta)*(-8.0*m12_2+(3.0*mHl2+2.0*mHp2)*2.0*sinb*cosb),2))/(64.0*pow(cosb*sinb,4)) * ddpB0mhmHpmHp
        + (sin(alpha)*sin(alpha)*pow((2.0*mA2-mHl2)*cos(alpha-3.0*beta)*2.0*sinb*cosb+cos(alpha+beta)*(8.0*m12_2-(2.0*mA2+3.0*mHl2)*2.0*sinb*cosb),2))/(128.0*pow(cosb*sinb,4)) * ddpB0mhmAmA
        + 1.5*mHh2*mHh2*cos(alpha)*cos(alpha)*cos(bma)*cos(bma) * ddpB0mH00
        + 2.0*(mHh2-mHp2)*(mHh2-mHp2)*cos(alpha)*cos(alpha)*sin(bma)*sin(bma) * ddpB0mH0mHp
        + (mA2-mHh2)*(mA2-mHh2)*cos(alpha)*cos(alpha)*sin(bma)*sin(bma) * ddpB0mH0mA
        + cos(alpha)*cos(alpha)*cos(bma)*cos(bma)*pow(m12_2+cos(alpha)*sin(alpha)*(mHh2-3.0*m12_2/(cosb*sinb))+mHl2*sin(2.0*alpha),2)/(2.0*cosb*cosb*sinb*sinb) * ddpB0mHmhmh
        + cos(alpha)*cos(alpha)*sin(bma)*sin(bma)*pow(m12_2*cosb*sinb+0.5*sin(2.0*alpha)*(3.0*m12_2-(2.0*mHh2+mHl2)*cosb*sinb),2)/pow(cosb*sinb,4) * ddpB0mHmHmh
        + 9.0*cos(alpha)*cos(alpha)*pow(mHh2*(-3.0*cos(bma)+cos(3.0*bma)-cos(3.0*alpha+beta)+3.0*cos(alpha+3.0*beta))
                                        +16.0*m12_2*sin(bma)*sin(bma)*sin(alpha+beta),2)/(512.0*pow(cosb*sinb,4)) * ddpB0mHmHmH
        + (cos(alpha)*cos(alpha)*pow((mHh2-2.0*mHp2)*cos(alpha-5.0*beta)
                                     +2.0*(mHh2+2.0*mHp2)*cos(bma)-(3.0*mHh2+2.0*mHp2)*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta),2))/(256.0*pow(cosb*sinb,4)) * ddpB0mHmHpmHp
        + (cos(alpha)*cos(alpha)*pow((2.0*mA2-mHh2)*cos(alpha-5.0*beta)
                                     -2.0*(2.0*mA2+mHh2)*cos(bma)+(2.0*mA2+3.0*mHh2)*cos(alpha+3.0*beta)+16.0*m12_2*sin(alpha+beta),2))/(512.0*pow(cosb*sinb,4)) * ddpB0mHmAmA
        + 2.0*(mHl2-mHp2)*(mHl2-mHp2)*cos(bma)*cos(bma)*sinb*sinb * ddpB0mHp0mh
        + 2.0*(mHh2-mHp2)*(mHh2-mHp2)*sin(bma)*sin(bma)*sinb*sinb * ddpB0mHp0mH
        + 2.0*(mA2-mHp2)*(mA2-mHp2)*sinb*sinb * ddpB0mHp0mA
        + 2.0*pow((m12_2*cos(alpha+beta))/(sinb*sinb*cosb*cosb)-(mHl2*cos(bma)*cos(2.0*beta))/(cosb*sinb)-(mHl2+2.0*mHp2)*sin(bma),2)*sinb*sinb * ddpB0mHpmHpmh
        + 2.0*pow(sinb*((mHh2+2.0*mHp2)*cos(bma)-mHh2*cos(2.0*beta)*sin(bma)/(cosb*sinb)-m12_2*sin(alpha+beta)/(sinb*sinb*cosb*cosb)),2) * ddpB0mHpmHpmH
        + (mA2-mHl2)*(mA2-mHl2)*cos(bma)*cos(bma)*sinb*sinb * ddpB0mA0mh
        + (mA2-mHh2)*(mA2-mHh2)*sin(bma)*sin(bma)*sinb*sinb * ddpB0mA0mH
        + 2.0*(mA2-mHp2)*(mA2-mHp2)*sinb*sinb * ddpB0mA0mHp
        + pow((2.0*mA2-mHl2)*cos(alpha-3.0*beta)*2.0*sinb*cosb+cos(alpha+beta)*(8.0*m12_2-(2.0*mA2+3.0*mHl2)*2.0*sinb*cosb),2)/(64.0*pow(cosb,4)*sinb*sinb) * ddpB0mAmAmh
        + pow((2.0*mA2-mHh2)*cos(alpha-5.0*beta)-2.0*(2.0*mA2+mHh2)*cos(bma)+2.0*mA2*cos(alpha+3.0*beta) 
               + 3.0*mHh2*cos(alpha+3.0*beta)+16.0*m12_2*sin(alpha+beta),2)/(256.0*pow(cosb,4)*sinb*sinb) * ddpB0mAmAmH;

    WFRcomb1b = (mHl2*(mA2*(2.0*mHl2-3.0*mHp2)+mHl2*mHp2)*sin(2.0*bma)*2.0*sinb*cosb)/(2.0*mA2*mHp2) * B000mh
        - (mHh2*(mA2*(2.0*mHh2-3.0*mHp2)+mHh2*mHp2)*sin(2.0*bma)*2.0*sinb*cosb)/(2.0*mA2*mHp2) * B000mH
        + ((mHl2-mHp2)*cos(bma)*((mHl2-2.0*mHp2)*cos(alpha-3.0*beta)*2.0*sinb*cosb
                                 +cos(alpha+beta)*(-8.0*m12_2+(3.0*mHl2+2.0*mHp2)*2.0*sinb*cosb)))/(2.0*mHp2*sinb*cosb) * B00mHpmh
        + ((mHp2-mHh2)*sin(bma)*((mHh2-2.0*mHp2)*cos(alpha-5.0*beta)+2.0*(mHh2+2.0*mHp2)*cos(bma)-3.0*mHh2*cos(alpha+3.0*beta)
                                 -2.0*mHp2*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta)))/(4.0*mHp2*sinb*cosb) * B00mHpmH
        + ((mHl2-mA2)*cos(bma)*((mHl2-2.0*mA2)*cos(alpha-3.0*beta)*2.0*sinb*cosb
                                +cos(alpha+beta)*(-8.0*m12_2+(2.0*mA2+3.0*mHl2)*2.0*sinb*cosb)))/(4.0*mA2*sinb*cosb) * B00mAmh
        + ((mHh2-mA2)*sin(bma)*((2.0*mA2-mHh2)*cos(alpha-5.0*beta)-2.0*(2.0*mA2+mHh2)*cos(bma)+2.0*mA2*cos(alpha+3.0*beta)
                                +3.0*mHh2*cos(alpha+3.0*beta)+16.0*m12_2*sin(alpha+beta)))/(8.0*mA2*sinb*cosb) * B00mAmH
        + (3.0*mHh2*mHl2*sin(2.0*alpha)*sin(2.0*bma))/(4.0*(mHh2-mHl2)) * B0mh00
        + ((mHl2-mHp2)*(mHp2-mHh2)*sin(2.0*alpha)*sin(2.0*bma))/(mHh2-mHl2) * B0mh0mHp
        - ((mA2-mHh2)*(mA2-mHl2)*sin(2.0*alpha)*sin(2.0*bma))/(2.0*(mHh2-mHl2)) * B0mh0mA
        + 3.0*cos(bma)*sin(2.0*alpha)*(m12_2+cos(alpha)*sin(alpha)*(mHh2-(3.0*m12_2)/(sinb*cosb))+mHl2*sin(2.0*alpha))
          *(-mHl2*(3.0*sin(bma)+sin(3.0*bma)+sin(3.0*alpha+beta)+3.0*sin(alpha+3.0*beta))
            +16.0*m12_2*cos(bma)*cos(bma)*cos(alpha+beta))/(32.0*(mHl2-mHh2)*pow(sinb*cosb,3)) * B0mhmhmh
        + (sin(2.0*bma)*sin(2.0*alpha)
           *(4.0*m12_2*m12_2+2.0*m12_2*(mHl2-mHh2)*sin(2.0*alpha)
             -(2.0*mHh2*mHh2+5.0*mHh2*mHl2+2.0*mHl2*mHl2-(9.0*m12_2*(mHh2+mHl2))/(sinb*cosb)
               +(9.0*m12_2*m12_2)/(sinb*sinb*cosb*cosb))*sin(2.0*alpha)*sin(2.0*alpha)))/(8.0*(mHh2-mHl2)*sinb*sinb*cosb*cosb) * B0mhmHmh
        - 3.0*sin(2.0*alpha)*sin(bma)*(m12_2-cos(alpha)*sin(alpha)*(2.0*mHh2+mHl2-(3.0*m12_2)/(sinb*cosb)))
          *(mHh2*(-3.0*cos(bma)+cos(3.0*bma)-cos(3.0*alpha+beta)+3.0*cos(alpha+3.0*beta))
            +16.0*m12_2*sin(bma)*sin(bma)*sin(alpha+beta))/(32.0*(mHh2-mHl2)*pow(sinb*cosb,3)) * B0mhmHmH
        + sin(2.0*alpha)*((mHl2-2.0*mHp2)*cos(alpha-3.0*beta)*2.0*sinb*cosb+cos(alpha+beta)*(-8.0*m12_2+(3.0*mHl2+2.0*mHp2)*2.0*sinb*cosb))
          *((mHh2-2.0*mHp2)*cos(alpha-5.0*beta)+2.0*(mHh2+2.0*mHp2)*cos(bma)-3.0*mHh2*cos(alpha+3.0*beta)
             -2.0*mHp2*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta))/(128.0*(mHh2-mHl2)*pow(cosb*sinb,4)) * B0mhmHpmHp
        + sin(2.0*alpha)*((mHl2-2.0*mA2)*cos(alpha-3.0*beta)*2.0*sinb*cosb+cos(alpha+beta)*(-8.0*m12_2+(2.0*mA2+3.0*mHl2)*2.0*sinb*cosb))
          *((mHh2-2.0*mA2)*cos(alpha-5.0*beta)+2.0*(2.0*mA2+mHh2)*cos(bma)-2.0*mA2*cos(alpha+3.0*beta)
             -3.0*mHh2*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta))/(256.0*(mHh2-mHl2)*pow(cosb*sinb,4)) * B0mhmAmA
        - (3.0*mHh2*mHl2*sin(2.0*alpha)*sin(2.0*bma))/(4.0*(mHh2-mHl2)) * B0mH00
        + ((mHh2-mHp2)*(mHl2-mHp2)*sin(2.0*alpha)*sin(2.0*bma))/(mHh2-mHl2) * B0mH0mHp
        + ((mA2-mHh2)*(mA2-mHl2)*sin(2.0*alpha)*sin(2.0*bma))/(2.0*(mHh2-mHl2)) * B0mH0mA
        + 3.0*cos(bma)*sin(2.0*alpha)*(m12_2+cos(alpha)*sin(alpha)*(mHh2-(3.0*m12_2)/(sinb*cosb))+mHl2*sin(2.0*alpha))
          *(-mHl2*(3.0*sin(bma)+sin(3.0*bma)+sin(3.0*alpha+beta)+3.0*sin(alpha+3.0*beta))
            +16.0*m12_2*cos(bma)*cos(bma)*cos(alpha+beta))/(32.0*(mHh2-mHl2)*pow(sinb*cosb,3)) * B0mHmhmh
        + (sin(2.0*bma)*sin(2.0*alpha)*(-4.0*m12_2*m12_2+2.0*m12_2*(mHh2-mHl2)*sin(2.0*alpha)
                                        +(2.0*mHh2*mHh2+5.0*mHh2*mHl2+2.0*mHl2*mHl2
                                          -(9.0*m12_2*(mHh2+mHl2))/(sinb*cosb)
                                          +(9.0*m12_2*m12_2)/(sinb*sinb*cosb*cosb))*sin(2.0*alpha)*sin(2.0*alpha)))
          /(8.0*(mHh2-mHl2)*sinb*sinb*cosb*cosb) * B0mHmHmh
        + 3.0*sin(bma)*sin(2.0*alpha)*(m12_2-cos(alpha)*sin(alpha)*(2.0*mHh2+mHl2-(3.0*m12_2)/(sinb*cosb)))
          *(mHh2*(-3.0*cos(bma)+cos(3.0*bma)-cos(3.0*alpha+beta)+3.0*cos(alpha+3.0*beta))
            +16.0*m12_2*sin(bma)*sin(bma)*sin(alpha+beta))/(32.0*(mHh2-mHl2)*pow(sinb*cosb,3)) * B0mHmHmH
        - sin(2.0*alpha)*((mHl2-2.0*mHp2)*cos(alpha-3.0*beta)*2.0*sinb*cosb+cos(alpha+beta)*(-8.0*m12_2+(3.0*mHl2+2.0*mHp2)*2.0*sinb*cosb))
          *((mHh2-2.0*mHp2)*cos(alpha-5.0*beta)+2.0*(mHh2+2.0*mHp2)*cos(bma)-3.0*mHh2*cos(alpha+3.0*beta)
            -2.0*mHp2*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta))/(128.0*(mHh2-mHl2)*pow(cosb*sinb,4)) * B0mHmHpmHp
        - sin(2.0*alpha)*((mHl2-2.0*mA2)*cos(alpha-3.0*beta)*2.0*sinb*cosb+cos(alpha+beta)*(-8.0*m12_2+(3.0*mHl2+2.0*mA2)*2.0*sinb*cosb))
          *((mHh2-2.0*mA2)*cos(alpha-5.0*beta)+2.0*(mHh2+2.0*mA2)*cos(bma)-3.0*mHh2*cos(alpha+3.0*beta)
            -2.0*mA2*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta))/(256.0*(mHh2-mHl2)*pow(cosb*sinb,4)) * B0mHmAmA
        - (mHl2*(mHl2-mHp2)*sin(2.0*bma)*2.0*sinb*cosb)/mHp2 * B0mHp0mh
        + (mHh2*(mHh2-mHp2)*sin(2.0*bma)*2.0*sinb*cosb)/mHp2 * B0mHp0mH
        + ((mHp2-mHl2)*cos(bma)*((mHl2-2.0*mHp2)*cos(alpha-3.0*beta)*2.0*sinb*cosb
                                 +cos(alpha+beta)*(-8.0*m12_2+(3.0*mHl2+2.0*mHp2)*2.0*sinb*cosb)))/(2.0*mHp2*sinb*cosb) * B0mHpmHpmh
        + ((mHh2-mHp2)*sin(bma)*((mHh2-2.0*mHp2)*cos(alpha-5.0*beta)+2.0*(mHh2+2.0*mHp2)*cos(bma)-3.0*mHh2*cos(alpha+3.0*beta)
                                 -2.0*mHp2*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta)))/(4.0*mHp2*sinb*cosb) * B0mHpmHpmH
        + (mHl2*(mA2-mHl2)*sin(bma)*cos(bma)*2.0*sinb*cosb)/mA2 * B0mA0mh
        + ((mHh2-mA2)*mHh2*sin(bma)*cos(bma)*2.0*sinb*cosb)/mA2 * B0mA0mH
        + ((mA2-mHl2)*cos(bma)*((-2.0*mA2+mHl2)*cos(alpha-3.0*beta)*2.0*sinb*cosb
                                +cos(alpha+beta)*(-8.0*m12_2+(2.0*mA2+3.0*mHl2)*2.0*sinb*cosb)))/(4.0*mA2*sinb*cosb) * B0mAmAmh
        +((mA2-mHh2)*sin(bma)*((2.0*mA2-mHh2)*cos(alpha-5.0*beta)-2.0*(2.0*mA2+mHh2)*cos(bma)+2.0*mA2*cos(alpha+3.0*beta)
                               +3.0*mHh2*cos(alpha+3.0*beta)+16.0*m12_2*sin(alpha+beta)))/(8.0*mA2*sinb*cosb) * B0mAmAmH;

    WFRcomb2a = 1.5*mHl2*mHl2*sin(bma)*sin(bma) * ddpB000mh
        + 1.5*mHh2*mHh2*cos(bma)*cos(bma) * ddpB000mH
        + (mHl2-mHp2)*(mHl2-mHp2)*cos(bma)*cos(bma) * ddpB00mHpmh
        + (mHh2-mHp2)*(mHh2-mHp2)*sin(bma)*sin(bma) * ddpB00mHpmH
        + (mA2-mHp2)*(mA2-mHp2) * ddpB00mHpmA
        + 0.5*(mA2-mHl2)*(mA2-mHl2)*cos(bma)*cos(bma) * ddpB00mAmh
        + 0.5*(mA2-mHh2)*(mA2-mHh2)*sin(bma)*sin(bma) * ddpB00mAmH
        + 0.75*mHl2*mHl2*sin(bma)*sin(bma) * ddpB0mh00
        + (mHl2-mHp2)*(mHl2-mHp2)*cos(bma)*cos(bma) * ddpB0mh0mHp
        + 0.5*(mA2-mHl2)*(mA2-mHl2)*cos(bma)*cos(bma) * ddpB0mh0mA
        + 9.0*pow(16.0*m12_2*cos(bma)*cos(bma)*cos(alpha+beta)
                  -mHl2*(3.0*sin(bma)+sin(3.0*bma)+sin(3.0*alpha+beta)+3.0*sin(alpha+3.0*beta)),2)/(1024.0*pow(cosb*sinb,4)) * ddpB0mhmhmh
        + 0.5*pow(cos(alpha)/sinb + sin(alpha)/cosb,2)
          *pow(m12_2+cos(alpha)*sin(alpha)*(mHh2+2.0*mHl2-3.0*m12_2/(cosb*sinb)),2) * ddpB0mhmHmh
        + (pow(-2.0*m12_2+(2.0*mHh2+mHl2-3.0*m12_2/(cosb*sinb))*sin(2.0*alpha),2)
           *sin(bma)*sin(bma))/(16.0*cosb*cosb*sinb*sinb) * ddpB0mhmHmH
        + pow((mHl2-2.0*mHp2)*cos(alpha-3.0*beta)*cosb*sinb
              +cos(alpha+beta)*(-4.0*m12_2+(3.0*mHl2+2.0*mHp2)*cosb*sinb),2)/(32.0*pow(cosb*sinb,4)) * ddpB0mhmHpmHp
        + pow((2.0*mA2-mHl2)*cos(alpha-3.0*beta)*cosb*sinb
              +cos(alpha+beta)*(4.0*m12_2-(2.0*mA2+3.0*mHl2)*cosb*sinb),2)/(64.0*pow(cosb*sinb,4)) * ddpB0mhmAmA
        + 0.75*mHh2*mHh2*cos(bma)*cos(bma) * ddpB0mH00
        + (mHh2 - mHp2)*(mHh2 - mHp2)*sin(bma)*sin(bma) * ddpB0mH0mHp
        + 0.5*(mA2-mHh2)*(mA2-mHh2)*sin(bma)*sin(bma) * ddpB0mH0mA
        + 0.25*pow(cos(alpha)/sinb + sin(alpha)/cosb,2)
          *pow(m12_2+cos(alpha)*sin(alpha)*(mHh2+2.0*mHl2-3.0*m12_2/(cosb*sinb)),2) * ddpB0mHmhmh
        + (pow(-2.0*m12_2+(2.0*mHh2+mHl2-3.0*m12_2/(cosb*sinb))*sin(2.0*alpha),2)
           *sin(bma)*sin(bma))/(8.0*cosb*cosb*sinb*sinb) * ddpB0mHmHmh
        + 9.0*pow(mHh2*(-3.0*cos(bma)+cos(3.0*bma)-cos(3.0*alpha+beta)+3.0*cos(alpha+3.0*beta))
                  +16.0*m12_2*sin(bma)*sin(bma)*sin(alpha+beta),2)/(1024.0*pow(cosb*sinb,4)) * ddpB0mHmHmH
        + pow((mHh2-2.0*mHp2)*cos(alpha-5.0*beta)+2.0*(mHh2+2.0*mHp2)*cos(bma)-(3.0*mHh2+2.0*mHp2)*cos(alpha+3.0*beta)
           -16.0*m12_2*sin(alpha+beta),2)/(512*pow(cosb*sinb,4)) * ddpB0mHmHpmHp
        + pow((2.0*mA2-mHh2)*cos(alpha-5.0*beta)-2.0*(2.0*mA2+mHh2)*cos(bma)+(2.0*mA2+3.0*mHh2)*cos(alpha+3.0*beta)
           +16.0*m12_2*sin(alpha+beta),2)/(1024*pow(cosb*sinb,4)) * ddpB0mHmAmA
        + (mHl2-mHp2)*(mHl2-mHp2)*cos(bma)*cos(bma) * ddpB0mHp0mh
        + (mHh2-mHp2)*(mHh2-mHp2)*sin(bma)*sin(bma) * ddpB0mHp0mH
        + (mA2-mHp2)*(mA2-mHp2) * ddpB0mHp0mA
        + pow((m12_2*cos(alpha+beta))/(cosb*cosb*sinb*sinb)-(mHl2*cos(bma)*cos(2.0*beta))/(cosb*sinb)
              -(mHl2+2.0*mHp2)*sin(bma),2) * ddpB0mHpmHpmh
        + pow((mHh2+2.0*mHp2)*cos(bma)-(mHh2*cos(2.0*beta)*sin(bma))/(cosb*sinb)
              -(m12_2*sin(alpha+beta))/(cosb*cosb*sinb*sinb),2) * ddpB0mHpmHpmH
        + 0.5*(mA2-mHl2)*(mA2-mHl2)*cos(bma)*cos(bma) * ddpB0mA0mh
        + 0.5*(mA2-mHh2)*(mA2-mHh2)*sin(bma)*sin(bma) * ddpB0mA0mH
        + (mA2-mHp2)*(mA2-mHp2) * ddpB0mA0mHp
        + pow((2.0*mA2-mHl2)*cos(alpha-3.0*beta)*2.0*cosb*sinb
              +cos(alpha+beta)*(8.0*m12_2-(2.0*mA2+3.0*mHl2)*2.0*cosb*sinb),2)/(128.0*pow(cosb*sinb,4)) * ddpB0mAmAmh
        + pow((2.0*mA2-mHh2)*cos(alpha-5.0*beta)-2.0*(2.0*mA2+mHh2)*cos(bma)+(2.0*mA2+3.0*mHh2)*cos(alpha+3.0*beta)
              +16.0*m12_2*sin(alpha+beta),2)/(512.0*pow(cosb*sinb,4)) * ddpB0mAmAmH;

    WFRcomb3a = 0.5*mHl2*mHl2*(3.0-cos(2.0*beta))*sin(bma)*sin(bma) * ddpB000mh
        + 0.5*mHh2*mHh2*(3.0-cos(2.0*beta))*cos(bma)*cos(bma) * ddpB000mH
        + 2.0*(mHl2-mHp2)*(mHl2-mHp2)*cos(bma)*cos(bma)*sinb*sinb * ddpB00mHpmh
        + 2.0*(mHh2-mHp2)*(mHh2-mHp2)*sin(bma)*sin(bma)*sinb*sinb * ddpB00mHpmH
        + 2.0*(mA2-mHp2)*(mA2-mHp2)*sinb*sinb * ddpB00mHpmA
        + (mA2-mHl2)*(mA2-mHl2)*cos(bma)*cos(bma)*cosb*cosb * ddpB00mAmh
        + (mA2-mHh2)*(mA2-mHh2)*cosb*cosb*sin(bma)*sin(bma) * ddpB00mAmH
        + 1.5*mHl2*mHl2*sin(alpha)*sin(alpha)*sin(bma)*sin(bma) * ddpB0mh00
        + 2.0*(mHl2-mHp2)*(mHl2-mHp2)*cos(bma)*cos(bma)*sin(alpha)*sin(alpha) * ddpB0mh0mHp
        + (mA2-mHl2)*(mA2-mHl2)*cos(bma)*cos(bma)*sin(alpha)*sin(alpha) * ddpB0mh0mA
        + 9.0*sin(alpha)*sin(alpha)*pow(-mHl2*(3.0*sin(bma)+sin(3.0*bma)+sin(3.0*alpha+beta)+3.0*sin(alpha+3.0*beta))
                                        +16.0*m12_2*cos(bma)*cos(bma)*cos(alpha+beta),2)/(512.0*pow(cosb*sinb,4)) * ddpB0mhmhmh
        + sin(alpha)*sin(alpha)*pow((cos(alpha)/sinb+sin(alpha)/cosb)*(m12_2+cos(alpha)*sin(alpha)*(mHh2+2.0*mHl2-(3.0*m12_2)/(cosb*sinb))),2) * ddpB0mhmHmh
        + sin(alpha)*sin(alpha)*sin(bma)*sin(bma)*pow(-2.0*m12_2+(2.0*mHh2+mHl2-(3.0*m12_2)/(cosb*sinb))*sin(2.0*alpha),2)/(8.0*cosb*cosb*sinb*sinb) * ddpB0mhmHmH
        + (sin(alpha)*sin(alpha)*pow((mHl2-2.0*mHp2)*cos(alpha-3.0*beta)*2.0*sinb*cosb
                                     +cos(alpha+beta)*(-8.0*m12_2+(3.0*mHl2+2.0*mHp2)*2.0*sinb*cosb),2))/(64.0*pow(cosb*sinb,4)) * ddpB0mhmHpmHp
        + (sin(alpha)*sin(alpha)*pow((2.0*mA2-mHl2)*cos(alpha-3.0*beta)*2.0*sinb*cosb+cos(alpha+beta)*(8.0*m12_2-(2.0*mA2+3.0*mHl2)*2.0*sinb*cosb),2))/(128.0*pow(cosb*sinb,4)) * ddpB0mhmAmA
        + 1.5*mHh2*mHh2*cos(alpha)*cos(alpha)*cos(bma)*cos(bma) * ddpB0mH00
        + 2.0*(mHh2-mHp2)*(mHh2-mHp2)*cos(alpha)*cos(alpha)*sin(bma)*sin(bma) * ddpB0mH0mHp
        + (mA2-mHh2)*(mA2-mHh2)*cos(alpha)*cos(alpha)*sin(bma)*sin(bma) * ddpB0mH0mA
        + cos(alpha)*cos(alpha)*cos(bma)*cos(bma)*pow(m12_2+cos(alpha)*sin(alpha)*(mHh2-3.0*m12_2/(cosb*sinb))+mHl2*sin(2.0*alpha),2)/(2.0*cosb*cosb*sinb*sinb) * ddpB0mHmhmh
        + cos(alpha)*cos(alpha)*sin(bma)*sin(bma)*pow(m12_2*cosb*sinb+0.5*sin(2.0*alpha)*(3.0*m12_2-(2.0*mHh2+mHl2)*cosb*sinb),2)/pow(cosb*sinb,4) * ddpB0mHmHmh
        + 9.0*cos(alpha)*cos(alpha)*pow(mHh2*(-3.0*cos(bma)+cos(3.0*bma)-cos(3.0*alpha+beta)+3.0*cos(alpha+3.0*beta))
                                        +16.0*m12_2*sin(bma)*sin(bma)*sin(alpha+beta),2)/(512.0*pow(cosb*sinb,4)) * ddpB0mHmHmH
        + (cos(alpha)*cos(alpha)*pow((mHh2-2.0*mHp2)*cos(alpha-5.0*beta)
                                     +2.0*(mHh2+2.0*mHp2)*cos(bma)-(3.0*mHh2+2.0*mHp2)*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta),2))/(256.0*pow(cosb*sinb,4)) * ddpB0mHmHpmHp
        + (cos(alpha)*cos(alpha)*pow((2.0*mA2-mHh2)*cos(alpha-5.0*beta)
                                     -2.0*(2.0*mA2+mHh2)*cos(bma)+(2.0*mA2+3.0*mHh2)*cos(alpha+3.0*beta)+16.0*m12_2*sin(alpha+beta),2))/(512.0*pow(cosb*sinb,4)) * ddpB0mHmAmA
        + 2.0*(mHl2-mHp2)*(mHl2-mHp2)*cos(bma)*cos(bma)*cosb*cosb * ddpB0mHp0mh
        + 2.0*(mHh2-mHp2)*(mHh2-mHp2)*sin(bma)*sin(bma)*cosb*cosb * ddpB0mHp0mH
        + 2.0*(mA2-mHp2)*(mA2-mHp2)*cosb*cosb * ddpB0mHp0mA
        + 2.0*cosb*cosb*pow((m12_2*cos(alpha+beta))/(cosb*cosb*sinb*sinb)-(mHl2*cos(bma)*cos(2.0*beta))/(cosb*sinb)-(mHl2+2.0*mHp2)*sin(bma),2) * ddpB0mHpmHpmh
        + 2.0*cosb*cosb*pow((mHh2+2.0*mHp2)*cos(bma)-(mHh2*cos(2.0*beta)*sin(bma))/(cosb*sinb)-(m12_2*sin(alpha+beta))/(cosb*cosb*sinb*sinb),2) * ddpB0mHpmHpmH
        + (mA2-mHl2)*(mA2-mHl2)*cos(bma)*cos(bma)*sinb*sinb * ddpB0mA0mh
        + (mA2-mHh2)*(mA2-mHh2)*sin(bma)*sin(bma)*sinb*sinb * ddpB0mA0mH
        + 2.0*(mA2-mHp2)*(mA2-mHp2)*sinb*sinb * ddpB0mA0mHp
        + pow((2.0*mA2-mHl2)*cos(alpha-3.0*beta)*2.0*sinb*cosb+cos(alpha+beta)*(8.0*m12_2-(2.0*mA2+3.0*mHl2)*2.0*sinb*cosb),2)/(64.0*pow(cosb,4)*sinb*sinb) * ddpB0mAmAmh
        + pow((2.0*mA2-mHh2)*cos(alpha-5.0*beta)-2.0*(2.0*mA2+mHh2)*cos(bma)+2.0*mA2*cos(alpha+3.0*beta) 
               + 3.0*mHh2*cos(alpha+3.0*beta)+16.0*m12_2*sin(alpha+beta),2)/(256.0*pow(cosb,4)*sinb*sinb) * ddpB0mAmAmH;

    WFRcomb3b = ((mHl2*mHl2*mHp2+mA2*(-2.0*mHl2*mHl2+mHl2*mHp2))*sin(2.0*bma)*sinb*cosb)/(mA2*mHp2) * B000mh
        + ((-mHh2*mHh2*mHp2+mA2*(2.0*mHh2*mHh2-mHh2*mHp2))*sin(2.0*bma)*sinb*cosb)/(mA2*mHp2) * B000mH
        + ((mHp2-mHl2)*cos(bma)*((mHl2-2.0*mHp2)*cos(alpha-3.0*beta)*2.0*sinb*cosb
                                 +cos(alpha+beta)*(-8.0*m12_2+(3.0*mHl2+2.0*mHp2)*2.0*sinb*cosb)))/(2.0*mHp2*sinb*cosb) * B00mHpmh
        + ((mHh2-mHp2)*sin(bma)*((mHh2-2.0*mHp2)*cos(alpha-5.0*beta)+2.0*(mHh2+2.0*mHp2)*cos(bma)-3.0*mHh2*cos(alpha+3.0*beta)
                                 -2.0*mHp2*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta)))/(4.0*mHp2*sinb*cosb) * B00mHpmH
        + ((mHl2-mA2)*cos(bma)*((mHl2-2.0*mA2)*cos(alpha-3.0*beta)*2.0*sinb*cosb
                                +cos(alpha+beta)*(-8.0*m12_2+(2.0*mA2+3.0*mHl2)*2.0*sinb*cosb)))/(4.0*mA2*sinb*cosb) * B00mAmh
        + ((mHh2-mA2)*sin(bma)*((2.0*mA2-mHh2)*cos(alpha-5.0*beta)-2.0*(2.0*mA2+mHh2)*cos(bma)+2.0*mA2*cos(alpha+3.0*beta)
                                +3.0*mHh2*cos(alpha+3.0*beta)+16.0*m12_2*sin(alpha+beta)))/(8.0*mA2*sinb*cosb) * B00mAmH
        + (3.0*mHh2*mHl2*sin(2.0*alpha)*sin(2.0*bma))/(4.0*(mHh2-mHl2)) * B0mh00
        + ((mHl2-mHp2)*(mHp2-mHh2)*sin(2.0*alpha)*sin(2.0*bma))/(mHh2-mHl2) * B0mh0mHp
        - ((mA2-mHh2)*(mA2-mHl2)*sin(2.0*alpha)*sin(2.0*bma))/(2.0*(mHh2-mHl2)) * B0mh0mA
        + 3.0*cos(bma)*sin(2.0*alpha)*(m12_2+cos(alpha)*sin(alpha)*(mHh2-(3.0*m12_2)/(sinb*cosb))+mHl2*sin(2.0*alpha))
          *(-mHl2*(3.0*sin(bma)+sin(3.0*bma)+sin(3.0*alpha+beta)+3.0*sin(alpha+3.0*beta))
            +16.0*m12_2*cos(bma)*cos(bma)*cos(alpha+beta))/(32.0*(mHl2-mHh2)*pow(sinb*cosb,3)) * B0mhmhmh
        + (sin(2.0*bma)*sin(2.0*alpha)
           *(4.0*m12_2*m12_2+2.0*m12_2*(mHl2-mHh2)*sin(2.0*alpha)
             -(2.0*mHh2*mHh2+5.0*mHh2*mHl2+2.0*mHl2*mHl2-(9.0*m12_2*(mHh2+mHl2))/(sinb*cosb)
               +(9.0*m12_2*m12_2)/(sinb*sinb*cosb*cosb))*sin(2.0*alpha)*sin(2.0*alpha)))/(8.0*(mHh2-mHl2)*sinb*sinb*cosb*cosb) * B0mhmHmh
        - 3.0*sin(2.0*alpha)*sin(bma)*(m12_2-cos(alpha)*sin(alpha)*(2.0*mHh2+mHl2-(3.0*m12_2)/(sinb*cosb)))
          *(mHh2*(-3.0*cos(bma)+cos(3.0*bma)-cos(3.0*alpha+beta)+3.0*cos(alpha+3.0*beta))
            +16.0*m12_2*sin(bma)*sin(bma)*sin(alpha+beta))/(32.0*(mHh2-mHl2)*pow(sinb*cosb,3)) * B0mhmHmH
        + sin(2.0*alpha)*((mHl2-2.0*mHp2)*cos(alpha-3.0*beta)*2.0*sinb*cosb+cos(alpha+beta)*(-8.0*m12_2+(3.0*mHl2+2.0*mHp2)*2.0*sinb*cosb))
          *((mHh2-2.0*mHp2)*cos(alpha-5.0*beta)+2.0*(mHh2+2.0*mHp2)*cos(bma)-3.0*mHh2*cos(alpha+3.0*beta)
             -2.0*mHp2*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta))/(128.0*(mHh2-mHl2)*pow(cosb*sinb,4)) * B0mhmHpmHp
        + sin(2.0*alpha)*((mHl2-2.0*mA2)*cos(alpha-3.0*beta)*2.0*sinb*cosb+cos(alpha+beta)*(-8.0*m12_2+(2.0*mA2+3.0*mHl2)*2.0*sinb*cosb))
          *((mHh2-2.0*mA2)*cos(alpha-5.0*beta)+2.0*(2.0*mA2+mHh2)*cos(bma)-2.0*mA2*cos(alpha+3.0*beta)
             -3.0*mHh2*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta))/(256.0*(mHh2-mHl2)*pow(cosb*sinb,4)) * B0mhmAmA
        - (3.0*mHh2*mHl2*sin(2.0*alpha)*sin(2.0*bma))/(4.0*(mHh2-mHl2)) * B0mH00
        + ((mHh2-mHp2)*(mHl2-mHp2)*sin(2.0*alpha)*sin(2.0*bma))/(mHh2-mHl2) * B0mH0mHp
        + ((mA2-mHh2)*(mA2-mHl2)*sin(2.0*alpha)*sin(2.0*bma))/(2.0*(mHh2-mHl2)) * B0mH0mA
        + 3.0*cos(bma)*sin(2.0*alpha)*(m12_2+cos(alpha)*sin(alpha)*(mHh2-(3.0*m12_2)/(sinb*cosb))+mHl2*sin(2.0*alpha))
          *(-mHl2*(3.0*sin(bma)+sin(3.0*bma)+sin(3.0*alpha+beta)+3.0*sin(alpha+3.0*beta))
            +16.0*m12_2*cos(bma)*cos(bma)*cos(alpha+beta))/(32.0*(mHh2-mHl2)*pow(sinb*cosb,3)) * B0mHmhmh
        + (sin(2.0*bma)*sin(2.0*alpha)*(-4.0*m12_2*m12_2+2.0*m12_2*(mHh2-mHl2)*sin(2.0*alpha)
                                        +(2.0*mHh2*mHh2+5.0*mHh2*mHl2+2.0*mHl2*mHl2
                                          -(9.0*m12_2*(mHh2+mHl2))/(sinb*cosb)
                                          +(9.0*m12_2*m12_2)/(sinb*sinb*cosb*cosb))*sin(2.0*alpha)*sin(2.0*alpha)))
          /(8.0*(mHh2-mHl2)*sinb*sinb*cosb*cosb) * B0mHmHmh
        + 3.0*sin(bma)*sin(2.0*alpha)*(m12_2-cos(alpha)*sin(alpha)*(2.0*mHh2+mHl2-(3.0*m12_2)/(sinb*cosb)))
          *(mHh2*(-3.0*cos(bma)+cos(3.0*bma)-cos(3.0*alpha+beta)+3.0*cos(alpha+3.0*beta))
            +16.0*m12_2*sin(bma)*sin(bma)*sin(alpha+beta))/(32.0*(mHh2-mHl2)*pow(sinb*cosb,3)) * B0mHmHmH
        - sin(2.0*alpha)*((mHl2-2.0*mHp2)*cos(alpha-3.0*beta)*2.0*sinb*cosb+cos(alpha+beta)*(-8.0*m12_2+(3.0*mHl2+2.0*mHp2)*2.0*sinb*cosb))
          *((mHh2-2.0*mHp2)*cos(alpha-5.0*beta)+2.0*(mHh2+2.0*mHp2)*cos(bma)-3.0*mHh2*cos(alpha+3.0*beta)
            -2.0*mHp2*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta))/(128.0*(mHh2-mHl2)*pow(cosb*sinb,4)) * B0mHmHpmHp
        - sin(2.0*alpha)*((mHl2-2.0*mA2)*cos(alpha-3.0*beta)*2.0*sinb*cosb+cos(alpha+beta)*(-8.0*m12_2+(3.0*mHl2+2.0*mA2)*2.0*sinb*cosb))
          *((mHh2-2.0*mA2)*cos(alpha-5.0*beta)+2.0*(mHh2+2.0*mA2)*cos(bma)-3.0*mHh2*cos(alpha+3.0*beta)
            -2.0*mA2*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta))/(256.0*(mHh2-mHl2)*pow(cosb*sinb,4)) * B0mHmAmA
        + (mHl2*(mHl2-mHp2)*sin(2.0*bma)*2.0*sinb*cosb)/mHp2 * B0mHp0mh
        - (mHh2*(mHh2-mHp2)*sin(2.0*bma)*2.0*sinb*cosb)/mHp2 * B0mHp0mH
        + (mHl2-mHp2)*cos(bma)*((mHl2-2.0*mHp2)*cos(alpha-3.0*beta)*2.0*sinb*cosb
                                +cos(alpha+beta)*(-8.0*m12_2+(3.0*mHl2+2.0*mHp2)*2.0*sinb*cosb))/(2.0*mHp2*sinb*cosb) * B0mHpmHpmh
        - (mHh2-mHp2)*sin(bma)*((mHh2-2.0*mHp2)*cos(alpha-5.0*beta)+2.0*(mHh2+2.0*mHp2)*cos(bma)
                                -3.0*mHh2*cos(alpha+3.0*beta)-2.0*mHp2*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta))/(4.0*mHp2*sinb*cosb) * B0mHpmHpmH
        + (mHl2*(mA2-mHl2)*sin(bma)*cos(bma)*2.0*sinb*cosb)/mA2 * B0mA0mh
        + ((mHh2-mA2)*mHh2*sin(bma)*cos(bma)*2.0*sinb*cosb)/mA2 * B0mA0mH
        + ((mA2-mHl2)*cos(bma)*((-2.0*mA2+mHl2)*cos(alpha-3.0*beta)*2.0*sinb*cosb
                                +cos(alpha+beta)*(-8.0*m12_2+(2.0*mA2+3.0*mHl2)*2.0*sinb*cosb)))/(4.0*mA2*sinb*cosb) * B0mAmAmh
        +((mA2-mHh2)*sin(bma)*((2.0*mA2-mHh2)*cos(alpha-5.0*beta)-2.0*(2.0*mA2+mHh2)*cos(bma)+2.0*mA2*cos(alpha+3.0*beta)
                               +3.0*mHh2*cos(alpha+3.0*beta)+16.0*m12_2*sin(alpha+beta)))/(8.0*mA2*sinb*cosb) * B0mAmAmH;

    WFRcomb4a = 0.5*mHl2*mHl2*sin(bma)*sin(bma) * ddpB000mh
        + 0.5*mHh2*mHh2*cos(bma)*cos(bma) * ddpB000mH
        + 0.5*(mA2-mHl2)*(mA2-mHl2)*cos(bma)*cos(bma) * ddpB00mAmh
        + 0.5*(mA2-mHh2)*(mA2-mHh2)*sin(bma)*sin(bma) * ddpB00mAmH
        + 0.75*mHl2*mHl2*sin(bma)*sin(bma) * ddpB0mh00
        + (mHl2-mHp2)*(mHl2-mHp2)*cos(bma)*cos(bma) * ddpB0mh0mHp
        + 0.5*(mA2-mHl2)*(mA2-mHl2)*cos(bma)*cos(bma) * ddpB0mh0mA
        + 9.0*pow(16.0*m12_2*cos(bma)*cos(bma)*cos(alpha+beta)
                  -mHl2*(3.0*sin(bma)+sin(3.0*bma)+sin(3.0*alpha+beta)+3.0*sin(alpha+3.0*beta)),2)/(1024.0*pow(cosb*sinb,4)) * ddpB0mhmhmh
        + 0.5*pow(cos(alpha)/sinb + sin(alpha)/cosb,2)
          *pow(m12_2+cos(alpha)*sin(alpha)*(mHh2+2.0*mHl2-3.0*m12_2/(cosb*sinb)),2) * ddpB0mhmHmh
        + (pow(-2.0*m12_2+(2.0*mHh2+mHl2-3.0*m12_2/(cosb*sinb))*sin(2.0*alpha),2)
           *sin(bma)*sin(bma))/(16.0*cosb*cosb*sinb*sinb) * ddpB0mhmHmH
        + pow((mHl2-2.0*mHp2)*cos(alpha-3.0*beta)*cosb*sinb
              +cos(alpha+beta)*(-4.0*m12_2+(3.0*mHl2+2.0*mHp2)*cosb*sinb),2)/(32.0*pow(cosb*sinb,4)) * ddpB0mhmHpmHp
        + pow((2.0*mA2-mHl2)*cos(alpha-3.0*beta)*cosb*sinb
              +cos(alpha+beta)*(4.0*m12_2-(2.0*mA2+3.0*mHl2)*cosb*sinb),2)/(64.0*pow(cosb*sinb,4)) * ddpB0mhmAmA
        + 0.75*mHh2*mHh2*cos(bma)*cos(bma) * ddpB0mH00
        + (mHh2 - mHp2)*(mHh2 - mHp2)*sin(bma)*sin(bma) * ddpB0mH0mHp
        + 0.5*(mA2-mHh2)*(mA2-mHh2)*sin(bma)*sin(bma) * ddpB0mH0mA
        + 0.25*pow(cos(alpha)/sinb + sin(alpha)/cosb,2)
          *pow(m12_2+cos(alpha)*sin(alpha)*(mHh2+2.0*mHl2-3.0*m12_2/(cosb*sinb)),2) * ddpB0mHmhmh
        + (pow(-2.0*m12_2+(2.0*mHh2+mHl2-3.0*m12_2/(cosb*sinb))*sin(2.0*alpha),2)
           *sin(bma)*sin(bma))/(8.0*cosb*cosb*sinb*sinb) * ddpB0mHmHmh
        + 9.0*pow(mHh2*(-3.0*cos(bma)+cos(3.0*bma)-cos(3.0*alpha+beta)+3.0*cos(alpha+3.0*beta))
                  +16.0*m12_2*sin(bma)*sin(bma)*sin(alpha+beta),2)/(1024.0*pow(cosb*sinb,4)) * ddpB0mHmHmH
        + pow((mHh2-2.0*mHp2)*cos(alpha-5.0*beta)+2.0*(mHh2+2.0*mHp2)*cos(bma)-(3.0*mHh2+2.0*mHp2)*cos(alpha+3.0*beta)
           -16.0*m12_2*sin(alpha+beta),2)/(512*pow(cosb*sinb,4)) * ddpB0mHmHpmHp
        + pow((2.0*mA2-mHh2)*cos(alpha-5.0*beta)-2.0*(2.0*mA2+mHh2)*cos(bma)+(2.0*mA2+3.0*mHh2)*cos(alpha+3.0*beta)
           +16.0*m12_2*sin(alpha+beta),2)/(1024*pow(cosb*sinb,4)) * ddpB0mHmAmA
        + 0.5*(mA2-mHl2)*(mA2-mHl2)*cos(bma)*cos(bma) * ddpB0mA0mh
        + 0.5*(mA2-mHh2)*(mA2-mHh2)*sin(bma)*sin(bma) * ddpB0mA0mH
        + (mA2-mHp2)*(mA2-mHp2) * ddpB0mA0mHp
        + pow((2.0*mA2-mHl2)*cos(alpha-3.0*beta)*2.0*cosb*sinb
              +cos(alpha+beta)*(8.0*m12_2-(2.0*mA2+3.0*mHl2)*2.0*cosb*sinb),2)/(128.0*pow(cosb*sinb,4)) * ddpB0mAmAmh
        + pow((2.0*mA2-mHh2)*cos(alpha-5.0*beta)-2.0*(2.0*mA2+mHh2)*cos(bma)+(2.0*mA2+3.0*mHh2)*cos(alpha+3.0*beta)
              +16.0*m12_2*sin(alpha+beta),2)/(512.0*pow(cosb*sinb,4)) * ddpB0mAmAmH;
    }
    WFRcomb1=-(WFRcomb1a+WFRcomb1b)/(vev*vev);
    WFRcomb2=-WFRcomb2a/(vev*vev);
    WFRcomb3=-(WFRcomb3a+WFRcomb3b)/(vev*vev);
    WFRcomb4=-WFRcomb4a/(vev*vev);
}

double THDMcache::updateCache()
{
    Q_THDM=myTHDM->getQ_THDM();
    bma=myTHDM->getbma();
    logtb=myTHDM->getlogtb();
    tanb=myTHDM->gettanb();
    m12_2=myTHDM->getm12_2();
    mA2=myTHDM->getmA2();
    mHp2=myTHDM->getmHp2();
    Rpeps=myTHDM->getRpeps();
    MW=MWTHDM(myTHDM->Mw_tree());
    cW2=cW2THDM(myTHDM->c02());
    /*This should be the only reference to the SM Higgs mass!*/
    mHl=myTHDM->getMHl();
    mHl2=mHl*mHl;
    mHh2=myTHDM->getmHh2();
    BDtaunu_SM=myTHDM->getBDtaunu_SM();
    BDtaunu_A=myTHDM->getBDtaunu_A();
    BDtaunu_B=myTHDM->getBDtaunu_B();
    BDstartaunu_SM=myTHDM->getBDstartaunu_SM();
    BDstartaunu_A=myTHDM->getBDstartaunu_A();
    BDstartaunu_B=myTHDM->getBDstartaunu_B();
    vev=myTHDM->v();
    Ale=myTHDM->getAle();
    Als=myTHDM->getAlsMz();
    Mt=myTHDM->getQuarks(QCD::TOP).getMass();
    Mb=myTHDM->getQuarks(QCD::BOTTOM).getMass();
    Mtau=myTHDM->getLeptons(StandardModel::TAU).getMass();
    Mc=myTHDM->getQuarks(QCD::CHARM).getMass();
    Ms=myTHDM->getQuarks(QCD::STRANGE).getMass();
    Mmu=myTHDM->getLeptons(StandardModel::MU).getMass();
    Mu=myTHDM->getQuarks(QCD::UP).getMass();
    Md=myTHDM->getQuarks(QCD::DOWN).getMass();
    Me=myTHDM->getLeptons(StandardModel::ELECTRON).getMass();
    MZ=myTHDM->getMz();
    modelflag=myTHDM->getModelTypeflag();
    WFRflag=myTHDM->getWFRflag();

    computeSignalStrengthQuantities();
    computeHHquantities();
    computeAquantities();
    computeHHlimits();
    computeAlimits();
    computeHpquantities();
    runTHDMparameters();
    computeWFRcombinations();
    
    return mHl2;
}
