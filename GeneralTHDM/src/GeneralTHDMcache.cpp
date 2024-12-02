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


//#include <boost/stacktrace.hpp>

//#include "log_cs_ggH_13.h"



GeneralTHDMcache::GeneralTHDMcache(const StandardModel& SM_i)
:       br_tt(19981, 2, 0.),
        br_bb(19981, 2, 0.),
        br_tautau(19981, 2, 0.),
        br_cc(19981, 2, 0.),
        br_mumu(19981, 2, 0.),
        br_ZZ(19981, 2, 0.),
        br_WW(19981, 2, 0.),
        GammaHtot_SM(19981, 2, 0.),
        log_cs_ggH_8(200, 2, 0.),
        log_cs_VBF_8(200, 2, 0.),
        log_cs_WH_8(200, 2, 0.),
        log_cs_ZH_8(200, 2, 0.),
        log_cs_ggH_13(200, 2, 0.),
        log_cs_VBF_13(200, 2, 0.),
        log_cs_WH_13(200, 2, 0.),
        log_cs_ZH_13(200, 2, 0.),
        log_cs_ttH_8(200, 2, 0.),
        log_cs_ttH_13(200, 2, 0.),
        log_cs_bbH_8(200, 2, 0.),
        log_cs_bbH_13(200, 2, 0.),
        log_cs_ggA_8(200, 2, 0.),
        log_cs_ttA_8(200, 2, 0.),
        log_cs_bbA_8(200, 2, 0.),
        log_cs_ggA_13(200, 2, 0.),
        log_cs_ttA_13(200, 2, 0.),
        log_cs_bbA_13(200, 2, 0.),
        log_cs_ggHp_8(744, 3, 0.),
        log_cs_ggHp_13(1104, 3, 0.),
        csrH_top_charm_8(200, 2, 0.),
        csrH_bottom_8(200, 2, 0.),
        csrA_top_charm_8(200, 2, 0.),
        csrA_bottom_8(200, 2, 0.),
        csrH_top_charm_13(200, 2, 0.),
        csrH_bottom_13(200, 2, 0.),
        csrA_top_charm_13(200, 2, 0.),
        csrA_bottom_13(200, 2, 0.),
        //
        CMS8_gg_phi_mumu(78, 2, 0.),
        CMS8_bb_phi_mumu(78, 2, 0.),
        CMS13_gg_phi_mumu(175, 2, 0.),
        CMS13_bb_phi_mumu(175, 2, 0.),
        ATLAS13_gg_phi_mumu(81, 2, 0.),
        ATLAS13_bb_phi_mumu(81, 2, 0.),
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
        CMS8_pp_phi_VV(172, 2, 0.),
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
        ATLAS13_bb_phi_bb(96, 2, 0.),               //Included in mid 2022
        CMS13_gg_phi_WW_heavy(71, 2, 0.),           //Included in mid 2022
        CMS13_VV_phi_WW_heavy(71, 2, 0.),           //Included in mid 2022
        CMS13_gg_phi_WW(281, 2, 0.),                //Included in mid 2022
        CMS13_VV_phi_WW(561, 2, 0.),                //Included in mid 2022
        ATLAS13_gg_phi_VV_llqq(95, 2, 0.),
        ATLAS13_VV_phi_VV_llqq(95, 2, 0.),
        //ATLAS13_tt_phi_tt(61, 2, 0.),             //OLD before mid 2022
        ATLAS13_tt_phi_tt(13, 2, 0.),               //Updated in mid 2022
        ATLAS13_bb_phi_tt(61, 2, 0.),
        //ATLAS13_gg_phi_tautau(206, 2, 0.),        //OLD before mid 2022
        ATLAS13_gg_phi_tautau(47, 2, 0.),           //Updated in mid 2022
        //ATLAS13_bb_phi_tautau(206, 2, 0.),        //OLD before mid 2022
        ATLAS13_bb_phi_tautau(47, 2, 0.),           //Updated in mid 2022
        //ATLAS13_pp_phi_gaga(251, 2, 0.),          //OLD before mid 2022
        ATLAS13_pp_phi_gaga(285, 2, 0.),            //Updated in mid 2022
        ATLAS13_gg_phi_Zga_llga(1590, 2, 0.),
        ATLAS13_gg_phi_Zga_qqga(291, 2, 0.),
        //ATLAS13_gg_phi_ZZ_llllnunu(101, 2, 0.),   //OLD before mid 2022
        //ATLAS13_VV_phi_ZZ_llllnunu(101, 2, 0.),   //OLD before mid 2022
        ATLAS13_gg_phi_ZZ_llllnunu(359, 2, 0.),     //Updated in mid 2022
        ATLAS13_VV_phi_ZZ_llllnunu(359, 2, 0.),     //Updated in mid 2022
        ATLAS13_gg_phi_ZZ_qqllnunu(271, 2, 0.),
        ATLAS13_VV_phi_ZZ_qqllnunu(271, 2, 0.),
        //ATLAS13_gg_phi_WW_enumunu(381, 2, 0.),    //OLD before mid 2024
        //ATLAS13_VV_phi_WW_enumunu(281, 2, 0.),    //OLD before mid 2024
        ATLAS13_gg_phi_WW_enumunu(72, 2, 0.),       //Updated in mid 2024
        ATLAS13_VV_phi_WW_enumunu(72, 2, 0.),       //Updated in mid 2024
        ATLAS13_gg_phi_WW_lnuqq(271, 2, 0.),
        ATLAS13_VV_phi_WW_lnuqq(271, 2, 0.),
        ATLAS13_pp_phi_VV_qqqq(181, 2, 0.),
        //ATLAS13_pp_phi_phi1phi1_bbbb(275, 2, 0.),     //OLD before mid 2022
        ATLAS13_pp_phi_phi1phi1_bbbb(476, 2, 0.),       //Updated in mid 2022
        //ATLAS13_pp_phi_phi1phi1_bbgaga(75, 2, 0.),    //OLD before mid 2022
        ATLAS13_pp_phi_phi1phi1_bbgaga(76, 2, 0.),      //Updated in mid 2022
        //ATLAS13_pp_phi_phi1phi1_bbtautau(75, 2, 0.),  //OLD before mid 2022
        ATLAS13_pp_phi_phi1phi1_bbtautau_1(136, 2, 0.), //Updated in mid 2022
        ATLAS13_pp_phi_phi1phi1_bbtautau_2(41, 2, 0.),  //Updated in mid 2022
        ATLAS13_pp_phi_phi1phi1_bbWW(51, 2, 0.),
        ATLAS13_gg_phi_phi1phi1_gagaWW(25, 2, 0.),
        ATLAS13_gg_phi_phi1Z_bbZ(45, 2, 0.),
        ATLAS13_bb_phi_phi1Z_bbZ(45, 2, 0.),
        //ATLAS13_gg_phii_phijZ_bbZ(3364, 3, 0.),   //OLD before mid 2022
        //ATLAS13_bb_phii_phijZ_bbZ(3364, 3, 0.),   //OLD before mid 2022
        ATLAS13_gg_phii_phijZ_bbZ(1711, 3, 0.),     //Updated in mid 2022
        ATLAS13_bb_phii_phijZ_bbZ(1711, 3, 0.),     //Updated in mid 2022
        ATLAS13_gg_phii_phijZ_WWZ(1326, 3, 0.),     //Included in mid 2024
        CMS13_gg_phi_phi1Z_tautaull(19, 2, 0.),
        CMS13_pp_phi2_bb_light(61, 2, 0.),          //Included in mid 2022
        CMS13_pp_phi3_bb_light(61, 2, 0.),          //Included in mid 2022
        CMS13_tt_phi2_tt(31, 2, 0.),               //Included in mid 2022
        CMS13_tt_phi3_tt(31, 2, 0.),               //Included in mid 2022
        CMS13_pp_phi_bb(66, 2, 0.),
        CMS13_bb_phi_bb(101, 2, 0.),
        //CMS13_gg_phi_tautau(312, 2, 0.),          //OLD before mid 2022
        //CMS13_bb_phi_tautau(312, 2, 0.),          //OLD before mid 2022
        //CMS13_gg_phi_tautau(689, 2, 0.),            //Updated in mid 2022
        //CMS13_bb_phi_tautau(689, 2, 0.),            //Updated in mid 2022
        CMS13_gg_phi_tautau(172, 2, 0.),            //Updated in mid 2024
        CMS13_bb_phi_tautau(172, 2, 0.),            //Updated in mid 2024
        CMS13_gg_phi_gaga(176, 2, 0.),
        CMS13_gg_phi_Zga(366, 2, 0.),
        CMS13_pp_phi_ZZ_llqqnunull(288, 2, 0.),
        CMS13_pp_phi_ZZ_qqnunu(301, 2, 0.),
        CMS13_ggVV_phi_WW_lnulnu(81, 2, 0.),
        CMS13_pp_phi_WW_lnuqq(341, 2, 0.),
        CMS13_pp_phi_phi1phi1_bbbb_1(95, 2, 0.),
        //CMS13_pp_phi_phi1phi1_bbbb_2(181, 2, 0.),     //OLD before mid 2022
        CMS13_pp_phi_phi1phi1_bbbb_2(41, 2, 0.),        //Updated in mid 2022 
        //CMS13_pp_phi_phi1phi1_bbgaga(66, 2, 0.),      //OLD before mid 2024
        CMS13_pp_phi_phi1phi1_bbgaga(74, 2, 0.),        //Updated in mid 2024 
        CMS13_pp_phi_phi1phi1_bbtautau_1(66, 2, 0.),
        CMS13_pp_phi_phi1phi1_bbtautau_2(311, 2, 0.),
        CMS13_pp_phi_phi1phi1_bbVV(65, 2, 0.),
        CMS13_gg_phi_phi1Z_bbZ_1(79, 2, 0.),
        CMS13_gg_phi_phi1Z_bbZ_2(121, 2, 0.),
        CMS13_bb_phi_phi1Z_bbZ_1(79, 2, 0.),
        CMS13_bb_phi_phi1Z_bbZ_2(121, 2, 0.),
        CMS13_pp_phi_phi1phi1_4WOr2W2tauOr4tau(76, 2, 0.), //Included in mid 2022 
        CMS13_pp_phi_phi1phi1_bbWW_qqlnu(55, 2, 0.), //Included in mid 2022 
        CMS13_pp_phi_phi1phi1_bbZZ_lljj(149, 2, 0.), //Included in mid 2022 
        CMS13_pp_phi_phi1phi1_bbZZ_llnunu(151, 2, 0.), //Included in mid 2022 
        CMS13_pp_phi_phi1phi1_bbWWorbbtautau(75, 2, 0.), //Included in mid 2022
        CMS13_pp_phi_phi1phi1_bbWWorbbtautau_low_masses(65, 2, 0.), //Included in mid 2024
        ATLAS8_pp_Hpm_taunu(83, 2, 0.),
        ATLAS8_pp_Hpm_tb(41, 2, 0.),
        CMS8_pp_Hp_taunu(43, 2, 0.),
        CMS8_pp_Hp_tb(43, 2, 0.),
        ATLAS13_pp_Hpm_taunu(192, 2, 0.),
        //ATLAS13_pp_Hpm_tb(181, 2, 0.),    //OLD before mid 2022
        ATLAS13_pp_Hpm_tb(181, 2, 0.),      //Updated in mid 2022 (in this case both have the same size)
        //CMS13_pp_Hpm_taunu(283, 2, 0.),   //OLD before mid 2022
        CMS13_pp_Hpm_taunu(585, 2, 0.),     //Updated in mid 2022
        CMS13_pp_Hpm_tb(281, 2, 0.),       //Included in mid 2022
        //
        CMS13_pp_h_phi3phi3_mumutautau(48, 2, 0.), //Added in late 2023
        //
    //  CMS13_pp_h_phi3phi3_bbtautau(10, 2, 0.), //from 1805.10191
        CMS13_pp_h_phi3phi3_bbtautau(11, 2, 0.), //Added in 2024
        //
    //  CMS13_pp_h_phi3phi3_bbmumu(86, 2, 0.), //from 1812.06359
        CMS13_pp_h_phi3phi3_bbmumu(95, 2, 0.), //Added in 2024
        //
        CMS13_pp_h_phi23Z_mumull(352, 2, 0.), //Added in 2024
        //
        CMS13_pp_h_phi23phi23_mumumumu(462, 2, 0.), //Added in 2024
        //
        CMS13_pp_h_phi3phi3_gagagaga(95, 2, 0.), //Added in 2024
        //
        CMS13_pp_h_phi3phi3_tautautautau(12, 2, 0.), //Added in 2024
        //
        CMS13_pp_phi2_gaga(81, 2, 0.), //Added in 2024
        //
        CMS13_pp_bbphi3_bbtautau(10, 2, 0.), //Added in 2024
        //
        CMS13_pp_ttphi3_ttmumu(220, 2, 0.), //Added in 2024
        //
        CMS13_pp_ttphi3_tttautau(18, 2, 0.), //Added in 2024
        //
        CMS13_Vh_h_phi3phi3_bbbb(10, 2, 0.), //Added in 2024
        //
        ATLAS13_pp_h_phi3phi3_bbmumu(47, 2, 0.), //Added in 2024
        //
        ATLAS13_gg_h_phi23phi23_mumumumu(110, 2, 0.), //Added in 2024
        //
        ATLAS13_gg_h_phi23Z_mumull(16, 2, 0.), //Added in 2024
        //
        ATLAS13_Vh_h_phi23phi23_bbbb(41, 2, 0.), //Added in 2024
        //
        ATLAS13_Zh_h_phi23phi23_bbbb(7, 2, 0.), //Added in 2024
        //
        ATLAS13_pp_h_phi23phi23_bbmumu_old(101, 2, 0.), //Added in 2024
        //
        ATLAS13_pp_h_phi23phi23_gagagg(13, 2, 0.), //Added in 2024
        //
        ATLAS13_pp_phi2_gaga_low(89, 2, 0.), //Added in 2024
        //
        ATLAS13_pp_ttphi3_ttmumu(37, 2, 0.), //Added in 2024
        //
        ATLAS13_pp_h_phi3phi3_gagagaga(53, 2, 0.), //Added in 2024
        // 
        ATLAS13_pp_h_phi3phi3_bbtautau(8, 2, 0.), //Added in 2024
        //
        ATLAS13_gg_phi3_tautau_low(12, 2, 0.),     //Added in 2024
        //
        ATLAS8_pp_h_phi3phi3_gagagaga(105, 2, 0.), //Added in 2024
        //
        ATLAS8_gg_h_phi3phi3_tautautautau(101, 2, 0.), //Added in 2024
        //
        CMS8_pp_h_phi3phi3_tautautautau(6, 2, 0.), //Added in 2024
        //
        CMS8_pp_h_phi3phi3_bbmumu(76, 2, 0.), //Added in 2024
        //
        CMS8_pp_h_phi3phi3_mumutautau(96, 2, 0.), //Added in 2024
        //
        CMS8_pp_phi2_gaga(300, 2, 0.), //Added in 2024
        //
        CMS8_pp_bbphi3_bbtautau(12, 2, 0.), //Added in 2024
        //
        CMS8_pp_bbphi3_bbmumu(36, 2, 0.), //Added in 2024
        //
        LEP209_phi2Z_gagaZ(49, 2, 0.),            //Added in 2024
        LEP209_phi2Z_bbZ(217, 2, 0.),             //Added in 2024
        LEP209_phi2Z_tautauZ(233, 2, 0.),         //Added in 2024
        LEP209_phiEphi3_bbbb(324, 3, 0.),         //Added in 2024
        LEP209_phiEphi3_tautautautau(365, 3, 0.), //Added in 2024
        //
        CMS8_t_Hpb_csb(7, 2, 0.),                 //Added in 2024
        CMS8_t_Hpb_taunub(8, 2, 0.),              //Added in 2024
        CMS8_t_Hpb_cbb(7, 2, 0.),                 //Added in 2024
        CMS13_t_Hpb_WAb_Wmumub(95, 2, 0.),        //Added in 2024
        CMS13_t_Hpb_csb(8, 2, 0.),                //Added in 2024
        ATLAS8_t_Hpb_taunub(9, 2, 0.),            //Added in 2024
        ATLAS13_t_Hpb_cbb(11, 2, 0.),             //Added in 2024
        ATLAS13_t_Hpb_WAb_Wmumub(31, 2, 0.),      //Added in 2024
        ATLAS13_t_Hpb_csb(12, 2, 0.),             //Added in 2024
        //
        //
        LEP209_HpHm_taunutaunu(105, 2, 0.),          //Added in 2024
        LEP209_HpHm_qqqq(105, 2, 0.),                //Added in 2024
        OPAL209_HpHm_qqtaunu(87, 2, 0.),             //Added in 2024
        OPAL172_HpHm_qqtaunu(41, 2, 0.),             //Added in 2024
        OPAL209_HpHm_AWstAWst_bbbb(2835, 3, 0.),     //Added in 2024
        OPAL209_HpHm_taunuAWst_taunubb(2835, 3, 0.), //Added in 2024
        //
        ATLAS13_pp_HpHm_taunutaunu(19, 2, 0.),    //Added in 2024
        ATLAS13_pp_HpHm_munumunu(11, 2, 0.),      //Added in 2024
        CMS13_pp_HpHm_taunutaunu(16, 2, 0.),      //Added in 2024
        CMS13_pp_HpHm_munumunu(22, 2, 0.),        //Added in 2024
        //
        csrHpHm_pp_13_LO(12, 2, 0.),
        csrHpHm_pp_13_NLO(31, 2, 0.),
        //Tables of integrals for g-2
        integral_x2_1mx_G_log(62500, 3, 0.),
        integral_x2_1px_G_log(62500, 3, 0.),
        integral_x2_G_log(62500, 3, 0.),
        integral_x_1mx2_G_log(62500, 3, 0.),
        integral_x_1mx_1px_G_log(62500, 3, 0.),
        
        integral_x2_1mx_G_variable_set_1_log(2000, 2, 0.),
        integral_x2_G_variable_set_1_log(2000, 2, 0.),
        integral_x_1mx2_G_variable_set_0_log(2000, 2, 0.),
        integral_perturb_ATHDM_values_log(2000, 2, 0.),
        //
        integral_Hp_A_Wstar(2835, 3, 0.),
        //
        //
        arraybsgamma(1111, 3, 0.),
        Rij_GTHDM(3, 3, 0.),
        //The below matrices are not used anywhere, the first ones are the mass matrices and the second ones not sure
        //In principle they should be the Yukawas in the Higgs basis but the up-type NP coupling (Nu) should be conjugated!!!
        Mu_GTHDM(3,3,0.), Md_GTHDM(3,3,0.), Ml_GTHDM(3,3,0.),
        Nu_GTHDM(3,3,0.), Nd_GTHDM(3,3,0.), Nl_GTHDM(3,3,0.),
        //The below matrices doesn't make sense for this model (they weren't defined in the higgs basis)
        //Yu1_GTHDM(3,3,0.), Yu2_GTHDM(3,3,0.), Yd1_GTHDM(3,3,0.), Yd2_GTHDM(3,3,0.),
        //Yl1_GTHDM(3,3,0.), Yl2_GTHDM(3,3,0.),
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
       
    if(x<xmin){
        std::cout<<"warning: your table parameter value is smaller than the minimum allowed value for m = "<<x<<"\n while xmin = "<<xmin<<std::endl;
        //std::cout<<"\033[1;36m   warning: your table parameter value is smaller than the minimum allowed value for m =  \033[0m "<<  x  <<std::endl;
        //std::cout<<"\033[1;36m   while min =  \033[0m "<<  xmin  <<std::endl;
        //std::cout << boost::stacktrace::stacktrace() << std::endl;
        
        
        return 0.;
    }
    else if(x>xmax){
        std::cout<<"warning: your table parameter value is greater than the maximum allowed value for m = "<<x<<"\n while xmax = "<<xmax<<std::endl;
        //std::cout<<"\033[1;36m   warning: your table parameter value is greater than the maximum allowed value for m =  \033[0m "<<  x  <<std::endl;
        //std::cout<<"\033[1;36m   while max =  \033[0m "<<  xmax  <<std::endl;
        //std::cout << boost::stacktrace::stacktrace() << std::endl;
        return 0.;
    }
    else if(x==xmax){
        return arrayTab(rowN-1,1);
    }
    else{
        return ((arrayTab(Nintervals+1,1)-arrayTab(Nintervals,1))/(arrayTab(Nintervals+1,0)
                -arrayTab(Nintervals,0))*(x-arrayTab(Nintervals,0))+arrayTab(Nintervals,1));
    }
}

//1D non-uniform interval interpolation

double GeneralTHDMcache::interpolateNU(gslpp::matrix<double> arrayTab, double x){

    int rowN=arrayTab.size_i();

    double xmin = arrayTab(0,0);
    double xmax = arrayTab(rowN-1,0);
    double test = x - (xmax + xmin) / 2;
    int row;

    if(test <= 0)
    {
        for(row = 1; row <= rowN-1; row++)
        {
            if(arrayTab(row,0) - x > 0)
            break;
        }
    }
    else
    {
        for(row = rowN-2; row >= 0; row--)
        {
            if(x - arrayTab(row,0) > 0)
            break;
        }

        row = row + 1;
    }

    if(x<xmin)
        return 0.;
    else if(x>xmax)
        return 0.;
    else if(x==xmax)
        return arrayTab(rowN-1,1);
    else
        return ((arrayTab(row,1)-arrayTab(row-1,1))/(arrayTab(row,0)
                -arrayTab(row-1,0))*(x-arrayTab(row-1,0))+arrayTab(row-1,1));
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
        std::cout<<"warning: the parameter point lies outside the table on interpolate2D"<<std::endl;
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

//2D interpolation with "y" changing first

double GeneralTHDMcache::interpolate2Dv2(gslpp::matrix<double> arrayTab, double x, double y){
    int rowN = arrayTab.size_i();

    double xmin = arrayTab(0,0);
    double xmax = arrayTab(rowN-1,0);
    double ymin = arrayTab(0,1);
    double ymax = arrayTab(rowN-1,1);
    double intervaly = arrayTab(1,1) - arrayTab(0,1);
    int i = 1;
    do i++;
    while (arrayTab(i,0) - arrayTab(i-1,0) == 0 && i < 30000);
    double intervalx = arrayTab(i,0) - arrayTab(i-1,0);
    
    int Nintervalsy = (y - ymin) / intervaly;
    int Nintervalsx = (x - xmin) / intervalx;
    if (x < xmin || x > xmax || y < ymin || y > ymax) {
        std::cout << "Warning: the parameter point lies outside the table on interpolate2Dv2" << std::endl;
        return 0.;
    } else {
        
        double x1 = arrayTab(Nintervalsy + Nintervalsx*i, 0);
        double x2 = arrayTab(Nintervalsy + (Nintervalsx+1)*i, 0);
        double y1 = arrayTab(Nintervalsy + Nintervalsx*i, 1);
        double y2 = arrayTab(Nintervalsy + 1 + Nintervalsx*i, 1);
        
        return ((arrayTab(Nintervalsy + Nintervalsx*i, 2) * (x2 - x) * (y2 - y)
                + arrayTab(Nintervalsy + (Nintervalsx + 1)*i, 2) * (x - x1) * (y2 - y)
                + arrayTab(Nintervalsy + 1 + Nintervalsx*i, 2) * (x2 - x) * (y - y1)
                + arrayTab((Nintervalsy + 1) + (Nintervalsx + 1)*i, 2) * (x - x1) * (y - y1))
               / ((x2 - x1) * (y2 - y1)));
        
        
    }
}

//2D interpolation change starts in y axis and triangular data

double GeneralTHDMcache::interpolate2DtriangularData(gslpp::matrix<double> arrayTab, double x, double y){

    int rowN=arrayTab.size_i();

    double xmin = arrayTab(0,0);
    double xmax = arrayTab(rowN-1,0);
    double ymin = arrayTab(0,1);
    double ymax = arrayTab(rowN-1,1);
    double intervaly = arrayTab(1,1)-arrayTab(0,1);
    int i=1;
    do i++;
    while(arrayTab(i,0)-arrayTab(i-1,0)==0&&i<30000);
    double intervalx = arrayTab(i,0)-arrayTab(i-1,0);
    int Nintervalsx = (x-xmin)/intervalx;
    int Nintervalsy = (y-ymin)/intervaly;
    if(x<xmin||x>xmax||y<ymin||y>ymax){
        std::cout<<"warning: the parameter point lies outside the table on interpolate2DtriangularData"<<std::endl;
        return 0.;
    }
    else{
    
        int Nx1y1=i*Nintervalsx-(Nintervalsx*(Nintervalsx+1)/2) +Nintervalsy;
        int Nx1y2=i*Nintervalsx-(Nintervalsx*(Nintervalsx+1)/2) +Nintervalsy+1;
        int Nx2y1=i*(Nintervalsx+1)-((Nintervalsx+1)*((Nintervalsx+1)+1)/2) +Nintervalsy;
        int Nx2y2=i*(Nintervalsx+1)-((Nintervalsx+1)*((Nintervalsx+1)+1)/2) +Nintervalsy+1;
        
        
    //std::cout<<" "<<std::endl;
    //std::cout<<" intervalx= "<<intervalx<<std::endl;
    //std::cout<<" intervaly= "<<intervaly<<std::endl;
    //std::cout<<" Nintervalsx= "<<Nintervalsx<<std::endl;
    //std::cout<<" Nintervalsy= "<<Nintervalsy<<std::endl;

    //std::cout<<" xmin= "<<xmin<<std::endl;
    //std::cout<<" xmax= "<<xmax<<std::endl;
    //std::cout<<" ymin= "<<ymin<<std::endl;
    //std::cout<<" ymax= "<<ymax<<std::endl;
    //std::cout<<" x= "<<x<<std::endl;
    //std::cout<<" y= "<<y<<std::endl;
    //std::cout<<"i*Nintervalsx+Nintervalsy ="<< i*Nintervalsx+Nintervalsy<<std::endl;
    //std::cout<<"i*Nintervalsx+Nintervalsy+1 ="<< i*Nintervalsx+Nintervalsy+1<<std::endl;
    //std::cout<<"i*Nintervalsx+Nintervalsy ="<< i*Nintervalsx+Nintervalsy<<std::endl;
    //std::cout<<"i*(Nintervalsx+1)+Nintervalsy ="<< i*(Nintervalsx+1)+Nintervalsy<<std::endl;
         
    //std::cout<<"y1(Nx1y1) = "<< arrayTab(Nx1y1,1)<<std::endl;
    //std::cout<<"y1(Nx2y1) = "<< arrayTab(Nx2y1,1)<<std::endl;
    //std::cout<<"y2(Nx1y2) = "<< arrayTab(Nx1y2,1)<<std::endl;
    //std::cout<<"y2(Nx2y2) = "<< arrayTab(Nx2y2,1)<<std::endl;
    
    
    //std::cout<<"x1(Nx1y1) = "<< arrayTab(Nx1y1,0)<<std::endl;
    //std::cout<<"x1(Nx1y2) = "<< arrayTab(Nx1y2,0)<<std::endl;
    //std::cout<<"x2(Nx2y1) = "<< arrayTab(Nx2y1,0)<<std::endl;
    //std::cout<<"x2(Nx2y2) = "<< arrayTab(Nx2y2,0)<<std::endl;
    //std::cout<<"Nx1y1 = "<< Nx1y1 <<std::endl;
    //std::cout<<"Nx1y2 = "<< Nx1y2 <<std::endl;
    //std::cout<<"Nx2y1 = "<< Nx2y1 <<std::endl;
    //std::cout<<"Nx2y2 = "<< Nx2y2 <<std::endl;
      
    double y1 = arrayTab(Nx1y1,1);
    double y2 = arrayTab(Nx1y2,1);
    double x1 = arrayTab(Nx1y1,0);  
    double x2 = arrayTab(Nx2y1,0);
    

    return (arrayTab(Nx2y2,2) * (y2-y) * (x2-x)
            +arrayTab(Nx2y1,2) * (y-y1) * (x2-x)
            +arrayTab(Nx1y2,2) * (y2-y) * (x-x1)
            +arrayTab(Nx1y1,2) * (y-y1) * (x-x1))
           /((x2-x1)*(y2-y1));
    }
}

//2D interpolation with fixed steps

double GeneralTHDMcache::interpolate2Dsteps(gslpp::matrix<double> arrayTab, double x, double y, double dx, double dy)
{
    unsigned int interpol, Npoints, Nbreak;
    std::vector<double> row, aux;
    double xdx = x / dx;
    double ydy = y / dy;

    if(xdx - ((int) xdx) == 0. && ydy - ((int) ydy) == 0.){
        interpol = 0;
        Npoints  = 1;
    }
    else if(xdx - ((int) xdx) == 0.){
        interpol = 1;
        Npoints  = 2;
    }
    else if(ydy - ((int) ydy) == 0.){
        interpol = 2;
        Npoints  = 2;
    }
    else {
        interpol = 3;
        Npoints  = 4;
    }

    Nbreak = 0;

    for(unsigned int i = 0; Nbreak < Npoints; i++)
    {
        if(i == arrayTab.size_i())
            break;

        if(std::fabs(x-arrayTab(i,0)) < dx && std::fabs(y-arrayTab(i,1)) < dy)
        {
            row.push_back(i);
            Nbreak = Nbreak + 1;
        }
    }

    if(row.size() == 3)
        interpol = 4;

    switch(interpol)
    {
        case 0:
            return (arrayTab(row[0],2));
        case 1:
            return (arrayTab(row[0],2) + (arrayTab(row[1],2)-arrayTab(row[0],2)) * (y-arrayTab(row[0],1)) / dy);
        case 2:
            return (arrayTab(row[0],2) + (arrayTab(row[1],2)-arrayTab(row[0],2)) * (x-arrayTab(row[0],0)) / dx);
        case 3:
            aux.push_back(arrayTab(row[0],2) + (arrayTab(row[1],2)-arrayTab(row[0],2)) * (y-arrayTab(row[0],1)) / dy);
            aux.push_back(arrayTab(row[2],2) + (arrayTab(row[3],2)-arrayTab(row[2],2)) * (y-arrayTab(row[2],1)) / dy);
            return (aux[0] + (aux[1]-aux[0]) * (x-arrayTab(row[0],0)) / dx);
        case 4:
            aux.push_back(arrayTab(row[0],2) + (arrayTab(row[1],2)-arrayTab(row[0],2)) * (y-arrayTab(row[0],1)) / dy);
            aux.push_back(arrayTab(row[1],2) + (arrayTab(row[2],2)-arrayTab(row[1],2)) * (y-arrayTab(row[0],1)) / dy);
            return (aux[0] + (aux[1]-aux[0]) * (x-arrayTab(row[0],0)) / dx);
        default:
            throw std::runtime_error("Review applicability of the mass ranges provided to GeneralTHDMcache::interpolate2Dsteps");
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
    std::stringstream ex1m6,ex1m5,ex1m4,ex1m3,ex1m2,ex1m1;
    std::stringstream ex1,ex2,ex3,ex4,ex5,ex6,ex7,ex8,ex9,ex10,ex11,ex12,ex13,ex14,ex15,ex16,ex17,ex18,ex19,ex20,ex21,ex22,ex23;
    std::stringstream ex24,ex25,ex26,ex27,ex28m1,ex28,ex29,ex30,ex31,ex32,ex33,ex34,ex35,ex36,ex37,ex38,ex39m4,ex39m3,ex39m2,ex39m1,ex39,ex40,ex41,ex42,ex43,ex43p1,ex43p2,ex44,\
            ex45,ex46n1,ex46n2,ex46a,ex47,ex48,ex49,ex49p2,ex50a,ex50b,ex51,ex52m2,ex52m1,ex52,ex52p1,ex52p2,ex53,ex54,ex55,ex56;
    std::stringstream ex57,ex58,ex59,ex60,ex61,ex62,ex63,ex64,ex65,ex66,ex67,ex67p1,ex67p2,ex67p3,ex67p4,ex67p5,ex67p6,ex68,ex69,ex70,ex71,ex72,ex73,ex74,ex75,ex76,ex77,\
            ex78,ex79;//,ex80,ex81,ex82,ex83,ex84,ex85,ex86,ex87,ex88,ex89,ex90,ex91,ex92,ex93,ex94,ex95,ex96,ex97,ex98
    std::stringstream lowC01,lowC02,lowC03,lowC04,lowC05,lowC06,lowC07,lowC08, lowC09, lowC10, lowC11, lowC12;
    std::stringstream lowA01, lowA02, lowA03, lowA04, lowA05, lowA06, lowA07, lowA08, lowA09, lowA10, lowA11, lowA12;
    std::stringstream lowA801, lowA802;
    std::stringstream lowC801, lowC802, lowC803, lowC804, lowC805, lowC806;
    std::stringstream low209a, low209b, low209c, low209d, low209e;
    std::stringstream lowHpC801, lowHpC802, lowHpC803, lowHpC1301, lowHpC1302;
    std::stringstream lowHpA801, lowHpA1301, lowHpA1302, lowHpA1303;
    std::stringstream lowHp209a, lowHp209b, lowHp209c, lowHp172a, lowHp209d, lowHp209e;
    std::stringstream thint01, thint02, thint03, thint04, thint05, thint06, thint07, thint08, thint09, thint10;
    std::stringstream csrslepLO, csrslepNLO, susyHpA01, susyHpA02, susyHpC01, susyHpC02;
    std::stringstream bsg1;

       std::cout<<"reading tables"<<std::endl;

//    std::cout << "HEPFITTABS = " << getenv("HEPFITPATH") << std::endl;
    std::stringstream path;
    path << getenv("HEPFITTABS") << "/THDM/tabs/";
    std::string tablepath=path.str();

    br1 << tablepath << "br1.dat";
    br_tt = readTable(br1.str(),19981,2);
    br2 << tablepath << "br2.dat";
    br_bb = readTable(br2.str(),19981,2);
    br3 << tablepath << "br3.dat";
    br_tautau = readTable(br3.str(),19981,2); 
    br4 << tablepath << "br4.dat";
    br_cc = readTable(br4.str(),19981,2);
    br5 << tablepath << "br5.dat";
    br_mumu = readTable(br5.str(),19981,2);
    br6 << tablepath << "br6.dat";
    br_ZZ = readTable(br6.str(),19981,2);
    br7 << tablepath << "br7.dat";
    br_WW = readTable(br7.str(),19981,2);
    dw1 << tablepath << "dw1.dat";
    GammaHtot_SM = readTable(dw1.str(),19981,2);
    cs1 << tablepath << "log_cs_ggH_8.dat";
    log_cs_ggH_8 = readTable(cs1.str(),200,2);
    cs11 << tablepath << "log_cs_ggH_13.dat";
    log_cs_ggH_13 = readTable(cs11.str(),200,2);
    cs2 << tablepath << "log_cs_VBF_8.dat";
    log_cs_VBF_8 = readTable(cs2.str(),200,2);
    cs12 << tablepath << "log_cs_VBF_13.dat";
    log_cs_VBF_13 = readTable(cs12.str(),200,2);
    cs3 << tablepath << "log_cs_WH_8.dat";
    log_cs_WH_8 = readTable(cs3.str(),200,2);
    cs13 << tablepath << "log_cs_WH_13.dat";
    log_cs_WH_13 = readTable(cs13.str(),200,2);
    cs4 << tablepath << "log_cs_ZH_8.dat";
    log_cs_ZH_8 = readTable(cs4.str(),200,2);
    cs14 << tablepath << "log_cs_ZH_13.dat";
    log_cs_ZH_13 = readTable(cs14.str(),200,2);
    cs5 << tablepath << "log_cs_ttH_8.dat";
    log_cs_ttH_8 = readTable(cs5.str(),200,2);
    cs15 << tablepath << "log_cs_ttH_13.dat";
    log_cs_ttH_13 = readTable(cs15.str(),200,2);
    cs6 << tablepath << "log_cs_bbH_8.dat";
    log_cs_bbH_8 = readTable(cs6.str(),200,2);
    cs16 << tablepath << "log_cs_bbH_13.dat";
    log_cs_bbH_13 = readTable(cs16.str(),200,2);
    cs7 << tablepath << "log_cs_ggA_8.dat";
    log_cs_ggA_8 = readTable(cs7.str(),200,2);
    cs17 << tablepath << "log_cs_ggA_13.dat";
    log_cs_ggA_13 = readTable(cs17.str(),200,2);
    cs8 << tablepath << "log_cs_ttA_8.dat";
    log_cs_ttA_8 = readTable(cs8.str(),200,2);
    cs18 << tablepath << "log_cs_ttA_13.dat";
    log_cs_ttA_13 = readTable(cs18.str(),200,2);
    cs9 << tablepath << "log_cs_bbA_8.dat";
    log_cs_bbA_8 = readTable(cs9.str(),200,2);
    cs19 << tablepath << "log_cs_bbA_13.dat";
    log_cs_bbA_13 = readTable(cs19.str(),200,2);
    cs20 << tablepath << "log_cs_ggHp_8.dat";
    log_cs_ggHp_8 = readTable(cs20.str(),744,3);
    cs21 << tablepath << "log_cs_ggHp_13.dat";
    log_cs_ggHp_13 = readTable(cs21.str(),1104,3);
    csr1 << tablepath << "csrH_top_charm_8.dat";
    csrH_top_charm_8 = readTable(csr1.str(),200,2);
    csr11 << tablepath << "csrH_top_charm_13.dat";
    csrH_top_charm_13 = readTable(csr11.str(),200,2);
    csr2 << tablepath << "csrH_bottom_8.dat";
    csrH_bottom_8 = readTable(csr2.str(),200,2);
    csr12 << tablepath << "csrH_bottom_13.dat";
    csrH_bottom_13 = readTable(csr12.str(),200,2);
    csr3 << tablepath << "csrA_top_charm_8.dat";
    csrA_top_charm_8 = readTable(csr3.str(),200,2);
    csr13 << tablepath << "csrA_top_charm_13.dat";
    csrA_top_charm_13 = readTable(csr13.str(),200,2);
    csr4 << tablepath << "csrA_bottom_8.dat";
    csrA_bottom_8 = readTable(csr4.str(),200,2);
    csr14 << tablepath << "csrA_bottom_13.dat";
    csrA_bottom_13 = readTable(csr14.str(),200,2);
    csrslepLO << tablepath << "csrHpHm_pp_slepton_13_LO.dat";
    csrHpHm_pp_13_LO = readTable(csrslepLO.str(),12,2);
    csrslepNLO << tablepath << "csrHpHm_pp_slepton_13_NLONLL.dat";
    csrHpHm_pp_13_NLO = readTable(csrslepNLO.str(),31,2);

    
    ex1m6<< tablepath << "150801437_9b.dat";                //Included in mid 2022
    CMS8_gg_phi_mumu = readTable(ex1m6.str(),78,2);         //Included in mid 2022
    ex1m5<< tablepath << "150801437_9a.dat";                //Included in mid 2022
    CMS8_bb_phi_mumu = readTable(ex1m5.str(),78,2);         //Included in mid 2022
    ex1m4<< tablepath << "190703152_6b.dat";                //Included in mid 2022
    CMS13_gg_phi_mumu = readTable(ex1m4.str(),175,2);       //Included in mid 2022
    ex1m3<< tablepath << "190703152_6a.dat";                //Included in mid 2022
    CMS13_bb_phi_mumu = readTable(ex1m3.str(),175,2);       //Included in mid 2022
    ex1m2<< tablepath << "190108144_4b.dat";                //Included in mid 2022
    ATLAS13_gg_phi_mumu = readTable(ex1m2.str(),81,2);      //Included in mid 2022
    ex1m1<< tablepath << "190108144_4a.dat";                //Included in mid 2022
    ATLAS13_bb_phi_mumu = readTable(ex1m1.str(),81,2);      //Included in mid 2022
    
    
    
    
    
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
    CMS8_pp_phi_VV = readTable(ex12.str(),172,2);
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

    
    
    ex28m1 << tablepath << "190702749.dat";                           //Included in mid 2022
    ATLAS13_bb_phi_bb = readTable(ex28m1.str(),96,2);                 //Included in mid 2022
    //ex28 << tablepath << "180711883.dat";                           //OLD previous to mid 2022
    //ATLAS13_tt_phi_tt = readTable(ex28.str(),61,2);                 //OLD previous to mid 2022
    ex28 << tablepath << "ATLAS_CONF_2022_008.dat";                   //Updated in mid 2022
    ATLAS13_tt_phi_tt = readTable(ex28.str(),13,2);                   //Updated in mid 2022
    ex29 << tablepath << "ATLAS-CONF-2016-104_b.dat";
    ATLAS13_bb_phi_tt = readTable(ex29.str(),61,2);
    
    
    
    //ex30 << tablepath << "170907242_a.dat";                       //OLD previous to mid 2022
    //ATLAS13_gg_phi_tautau = readTable(ex30.str(),206,2);          //OLD previous to mid 2022
    //ex31 << tablepath << "170907242_b.dat";                       //OLD previous to mid 2022
    //ATLAS13_bb_phi_tautau = readTable(ex31.str(),206,2);          //OLD previous to mid 2022
    
    
    
    ex30 << tablepath << "200212223_2a.dat";                       //Updated in mid 2022
    ATLAS13_gg_phi_tautau = readTable(ex30.str(),47,2);            //Updated in mid 2022
    ex31 << tablepath << "200212223_2b.dat";                       //Updated in mid 2022
    ATLAS13_bb_phi_tautau = readTable(ex31.str(),47,2);            //Updated in mid 2022
    
    
    
    //ex32 << tablepath << "170704147.dat";                         //OLD previous to mid 2022
    //ATLAS13_pp_phi_gaga = readTable(ex32.str(),251,2);            //OLD previous to mid 2022
    ex32 << tablepath << "210213405.dat";                           //Updated in mid 2022
    ATLAS13_pp_phi_gaga = readTable(ex32.str(),285,2);              //Updated in mid 2022
    
    
    
    //ex33 << tablepath << "170800212.dat";                           //OLD previous to mid 2024
    //ATLAS13_gg_phi_Zga_llga = readTable(ex33.str(),216,2);          //OLD previous to mid 2024
    ex33 << tablepath << "230904364.dat";                           //Updated in mid 2024
    ATLAS13_gg_phi_Zga_llga = readTable(ex33.str(),1590,2);          //Updated in mid 2024
    
    
    
    
    //ex34 << tablepath << "180501908.dat";                           //OLD previous to mid 2024
    //ATLAS13_gg_phi_Zga_qqga = readTable(ex34.str(),581,2);          //OLD previous to mid 2024
    ex34 << tablepath << "230411962.dat";                           //Updated in mid 2024
    ATLAS13_gg_phi_Zga_qqga = readTable(ex34.str(),291,2);          //Updated in mid 2024
    
    //ex35 << tablepath << "171206386_a.dat";                           //OLD previous to mid 2022
    //ATLAS13_gg_phi_ZZ_llllnunu = readTable(ex35.str(),101,2);         //OLD previous to mid 2022
    //ex36 << tablepath << "171206386_b.dat";                           //OLD previous to mid 2022
    //ATLAS13_VV_phi_ZZ_llllnunu = readTable(ex36.str(),101,2);         //OLD previous to mid 2022
    
    
    ex35 << tablepath << "200914791_4a.dat";                            //Updated in mid 2022
    ATLAS13_gg_phi_ZZ_llllnunu = readTable(ex35.str(),359,2);           //Updated in mid 2022
    ex36 << tablepath << "200914791_4b.dat";                            //Updated in mid 2022
    ATLAS13_VV_phi_ZZ_llllnunu = readTable(ex36.str(),359,2);           //Updated in mid 2022
    
    
    
    ex37 << tablepath << "170809638_a.dat";
    ATLAS13_gg_phi_ZZ_qqllnunu = readTable(ex37.str(),271,2);
    ex38 << tablepath << "170809638_b.dat";
    ATLAS13_VV_phi_ZZ_qqllnunu = readTable(ex38.str(),271,2);
    
    
    ex39m4 << tablepath << "210906055_7a.dat";                  //Included in mid 2022
    CMS13_gg_phi_WW_heavy = readTable(ex39m4.str(),71,2);       //Included in mid 2022
    ex39m3 << tablepath << "210906055_7b.dat";                  //Included in mid 2022
    CMS13_VV_phi_WW_heavy = readTable(ex39m3.str(),71,2);       //Included in mid 2022
    
    
    ex39m2 << tablepath << "191201594_6c.dat";                  //Included in mid 2022
    CMS13_gg_phi_WW = readTable(ex39m2.str(),281,2);            //Included in mid 2022
    ex39m1 << tablepath << "191201594_6d.dat";                  //Included in mid 2022
    CMS13_VV_phi_WW = readTable(ex39m1.str(),561,2);            //Included in mid 2022
    
    
    //ex39 << tablepath << "171001123_a.dat";                     //OLD previous to mid 2024
    //ATLAS13_gg_phi_WW_enumunu = readTable(ex39.str(),381,2);    //OLD previous to mid 2024
    //ex40 << tablepath << "171001123_b.dat";                     //OLD previous to mid 2024
    //ATLAS13_VV_phi_WW_enumunu = readTable(ex40.str(),281,2);    //OLD previous to mid 2024
    
    ex39 << tablepath << "ATLAS_CONF_2022_066_ggF.dat";                     //Updated in mid 2024
    ATLAS13_gg_phi_WW_enumunu = readTable(ex39.str(),72,2);    //Updated in mid 2024
    ex40 << tablepath << "ATLAS_CONF_2022_066_VBF.dat";                     //Updated in mid 2024
    ATLAS13_VV_phi_WW_enumunu = readTable(ex40.str(),72,2);    //Updated in mid 2024
    
    
    
    ex41 << tablepath << "171007235_a.dat";
    ATLAS13_gg_phi_WW_lnuqq = readTable(ex41.str(),271,2);
    ex42 << tablepath << "171007235_b.dat";
    ATLAS13_VV_phi_WW_lnuqq = readTable(ex42.str(),271,2);
    ex43 << tablepath << "170804445.dat";
    ATLAS13_pp_phi_VV_qqqq = readTable(ex43.str(),181,2);
    
    
    ex43p1 << tablepath << "200414636_12a.dat";
    ATLAS13_gg_phi_VV_llqq = readTable(ex43p1.str(),95,2);
    ex43p2 << tablepath << "200414636_12b.dat";
    ATLAS13_VV_phi_VV_llqq = readTable(ex43p2.str(),95,2);
    
    //ex44 << tablepath << "180406174.dat";                         //OLD previous mid 2022
    //ATLAS13_pp_phi_phi1phi1_bbbb = readTable(ex44.str(),275,2);   //OLD previous mid 2022
    ex44 << tablepath << "220207288.dat";                           //Updated in mid 2022
    ATLAS13_pp_phi_phi1phi1_bbbb = readTable(ex44.str(),476,2);     //Updated in mid 2022
    //ex45 << tablepath << "180704873.dat";                         //OLD previous mid 2022
    //ATLAS13_pp_phi_phi1phi1_bbgaga = readTable(ex45.str(),75,2);  //OLD previous mid 2022
    ex45 << tablepath << "211211876.dat";                           //Updated in mid 2022
    ATLAS13_pp_phi_phi1phi1_bbgaga = readTable(ex45.str(),76,2);    //Updated in mid 2022
    //ex46 << tablepath << "180800336.dat";                         //OLD previous mid 2022
    //ATLAS13_pp_phi_phi1phi1_bbtautau = readTable(ex46.str(),75,2);//OLD previous mid 2022
    ex46n1 << tablepath << "ATLAS-CONF-2021-030.dat";                 //Updated in mid 2022
    ATLAS13_pp_phi_phi1phi1_bbtautau_1 = readTable(ex46n1.str(),136,2);  //Updated in mid 2022
    ex46n2 << tablepath << "200714811.dat";                           //Updated in mid 2022
    ATLAS13_pp_phi_phi1phi1_bbtautau_2 = readTable(ex46n2.str(),41,2);  //Updated in mid 2022
    ex46a << tablepath << "181104671.dat";
    ATLAS13_pp_phi_phi1phi1_bbWW = readTable(ex46a.str(),51,2);
    ex47 << tablepath << "180708567.dat";
    ATLAS13_gg_phi_phi1phi1_gagaWW = readTable(ex47.str(),25,2);
    //ex48 << tablepath << "171206518_a.dat";                       //Superseded
    //ATLAS13_gg_phi_phi1Z_bbZ = readTable(ex48.str(),181,2);       //Superseded
    //ex49 << tablepath << "171206518_b.dat";                       //Superseded
    //ATLAS13_bb_phi_phi1Z_bbZ = readTable(ex49.str(),181,2);       //Superseded
    ex48 << tablepath << "220700230_ggF.dat";                       //Updated in mid 2024
    ATLAS13_gg_phi_phi1Z_bbZ = readTable(ex48.str(),45,2);          //Updated in mid 2024
    ex49 << tablepath << "220700230_bb.dat";                        //Updated in mid 2024
    ATLAS13_bb_phi_phi1Z_bbZ = readTable(ex49.str(),45,2);          //Updated in mid 2024
    
    ex49p2 << tablepath << "191011634.dat";
    CMS13_gg_phi_phi1Z_tautaull = readTable(ex49p2.str(),19,2);
    
    
    
    
    //ex50 << tablepath << "180401126_a.dat";                       //OLD previous mid 2022
    //ATLAS13_gg_phii_phijZ_bbZ = readTable(ex50.str(),3364,3);     //OLD previous mid 2022
    ex50a << tablepath << "201105639_ggF_decay_bb.dat";                       //Updated in mid 2022
    ATLAS13_gg_phii_phijZ_bbZ = readTable(ex50a.str(),1711,3);       //Updated in mid 2022
    //ex51 << tablepath << "180401126_b.dat";                       //OLD previous mid 2022
    //ATLAS13_bb_phii_phijZ_bbZ = readTable(ex51.str(),3364,3);     //OLD previous mid 2022
    ex51 << tablepath << "201105639_bb_decay_bb.dat";               //Updated in mid 2022
    ATLAS13_bb_phii_phijZ_bbZ = readTable(ex51.str(),1711,3);       //Updated in mid 2022
    
    ex50b << tablepath << "201105639_ggF_decay_WW.dat";              //Included in mid 2024
    ATLAS13_gg_phii_phijZ_WWZ = readTable(ex50b.str(),1326,3);       //Included in mid 2024

    
    ex52m2 << tablepath << "190806463_7a.dat";                      //Included in mid 2022
    CMS13_tt_phi2_tt = readTable(ex52m2.str(),31,2);                //Included in mid 2022
    ex52m1 << tablepath << "190806463_7b.dat";                      //Included in mid 2022
    CMS13_tt_phi3_tt = readTable(ex52m1.str(),31,2);                //Included in mid 2022
    
    
    ex52 << tablepath << "CMS-PAS-HIG-16-025.dat";
    CMS13_pp_phi_bb = readTable(ex52.str(),66,2);
    
    
    ex52p1 << tablepath << "181011822_7a.dat";                         //Included in mid 2022
    CMS13_pp_phi2_bb_light = readTable(ex52p1.str(),61,2);             //Included in mid 2022
    ex52p2 << tablepath << "181011822_8a.dat";                         //Included in mid 2022
    CMS13_pp_phi3_bb_light = readTable(ex52p2.str(),61,2);             //Included in mid 2022
    
    
    
    ex53 << tablepath << "180512191.dat";
    CMS13_bb_phi_bb = readTable(ex53.str(),101,2);
    
    
    //ex54 << tablepath << "180306553_a.dat";                       //OLD previous to mid 2022
    //CMS13_gg_phi_tautau = readTable(ex54.str(),312,2);            //OLD previous to mid 2022
    //ex55 << tablepath << "180306553_b.dat";                       //OLD previous to mid 2022
    //CMS13_bb_phi_tautau = readTable(ex55.str(),312,2);            //OLD previous to mid 2022
    
    
    //CMS_PAS_HIG_21_001_9b
    //ex54 << tablepath << "CMS_PAS_HIG_21_001_9a.dat";               //Updated in mid 2022
    //CMS13_gg_phi_tautau = readTable(ex54.str(),689,2);              //Updated in mid 2022
    //ex55 << tablepath << "CMS_PAS_HIG_21_001_9b.dat";               //Updated in mid 2022
    //CMS13_bb_phi_tautau = readTable(ex55.str(),689,2);              //Updated in mid 2022
    
    
    ex54 << tablepath << "220802717_ggF.dat";              //Updated in mid 2024
    CMS13_gg_phi_tautau = readTable(ex54.str(),172,2);     //Updated in mid 2024
    ex55 << tablepath << "220802717_bb.dat";               //Updated in mid 2024
    CMS13_bb_phi_tautau = readTable(ex55.str(),172,2);     //Updated in mid 2024
    
    
    
    //ex56 << tablepath << "160902507.dat";                         //OLD previous to mid 2022
    //CMS13_gg_phi_gaga = readTable(ex56.str(),351,2);              //OLD previous to mid 2022
    
    
    //ex56 << tablepath << "180900327.dat";                           //Updated in mid 2022
    //CMS13_gg_phi_gaga = readTable(ex56.str(),901,2);                //Updated in mid 2022
    
    ex56 << tablepath << "240509320.dat";                           //Updated in mid 2024
    CMS13_gg_phi_gaga = readTable(ex56.str(),176,2);                //Updated in mid 2024
    
    
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
    
        
    //ex63 << tablepath << "180801473.dat";                            //OLD previous to mid 2022
    //CMS13_pp_phi_phi1phi1_bbbb_2 = readTable(ex63.str(),181,2);      //OLD previous to mid 2022
    ex63 << tablepath << "CMS-PAS-B2G-20-004.dat";                     //Updated in mid 2022
    CMS13_pp_phi_phi1phi1_bbbb_2 = readTable(ex63.str(),41,2);         //Updated in mid 2022
    
    //ex64 << tablepath << "180600408.dat";                            //OLD previous to mid 2024
    //CMS13_pp_phi_phi1phi1_bbgaga = readTable(ex64.str(),66,2);       //OLD previous to mid 2024
    ex64 << tablepath << "231001643.dat";                              //Updated in mid 2024
    CMS13_pp_phi_phi1phi1_bbgaga = readTable(ex64.str(),74,2);         //Updated in mid 2024
    
    
    ex65 << tablepath << "170702909.dat";
    CMS13_pp_phi_phi1phi1_bbtautau_1 = readTable(ex65.str(),66,2);
    ex66 << tablepath << "180801365.dat";
    CMS13_pp_phi_phi1phi1_bbtautau_2 = readTable(ex66.str(),311,2);
    ex67 << tablepath << "170804188.dat";
    CMS13_pp_phi_phi1phi1_bbVV = readTable(ex67.str(),65,2);
    ex67p1 << tablepath << "220610268.dat";                                         //Included in mid 2022
    CMS13_pp_phi_phi1phi1_4WOr2W2tauOr4tau = readTable(ex67p1.str(),76,2);          //Included in mid 2022
    ex67p2 << tablepath << "190404193.dat";                                         //Included in mid 2022
    CMS13_pp_phi_phi1phi1_bbWW_qqlnu = readTable(ex67p2.str(),55,2);                //Included in mid 2022
    
    
    ex67p3 << tablepath << "200606391_bblljj.dat";                                  //Included in mid 2022
    CMS13_pp_phi_phi1phi1_bbZZ_lljj = readTable(ex67p3.str(),149,2);                //Included in mid 2022
    
    ex67p4 << tablepath << "200606391_bbllnunu.dat";                                //Included in mid 2022
    CMS13_pp_phi_phi1phi1_bbZZ_llnunu = readTable(ex67p4.str(),151,2);              //Included in mid 2022
    
    
    ex67p5 << tablepath << "211203161.dat";                                         //Included in mid 2022
    CMS13_pp_phi_phi1phi1_bbWWorbbtautau = readTable(ex67p5.str(),75,2);            //Included in mid 2022
    
    
    ex67p6 << tablepath << "240309430.dat";                                         //Included in mid 2024
    CMS13_pp_phi_phi1phi1_bbWWorbbtautau_low_masses = readTable(ex67p6.str(),65,2); //Included in mid 2024
    
    
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
    //ex77 << tablepath << "180803599.dat";             //OLD Previous to mid 2022
    //ATLAS13_pp_Hpm_tb = readTable(ex77.str(),181,2);  //OLD Previous to mid 2022
    ex77 << tablepath << "210210076.dat";               //Updated in mid 2022
    ATLAS13_pp_Hpm_tb = readTable(ex77.str(),181,2);    //Updated in mid 2022
    ex78 << tablepath << "200107763.dat";               //Included in mid 2022
    CMS13_pp_Hpm_tb = readTable(ex78.str(),281,2);    //Included in mid 2022
    //ex78 << tablepath << "CMS-PAS-HIG-16-031.dat";    //OLD Previous to mid 2022
    //CMS13_pp_Hpm_taunu = readTable(ex78.str(),283,2); //OLD Previous to mid 2022
    ex79 << tablepath << "190304560.dat";               //Updated in mid 2022
    CMS13_pp_Hpm_taunu = readTable(ex79.str(),585,2);   //Updated in mid 2022

    lowC01 << tablepath << "CMS-HIG-17-029_4e.dat";               //Added in late 2023
    CMS13_pp_h_phi3phi3_mumutautau = readTable(lowC01.str(),48,2);

//  lowC02 << tablepath << "CMS-HIG-17-024_7d.dat";               //Added in 2024
//  CMS13_pp_h_phi3phi3_bbtautau = readTable(lowC02.str(),10,2);

    lowC02 << tablepath << "CMS-HIG-22-007_11d.dat";               //Added in 2024
    CMS13_pp_h_phi3phi3_bbtautau = readTable(lowC02.str(),11,2);

//  lowC03 << tablepath << "CMS-HIG-18-011_5b.dat";               //Added in 2024
//  CMS13_pp_h_phi3phi3_bbmumu = readTable(lowC03.str(),86,2);

    lowC03 << tablepath << "CMS-HIG-22-007_10.dat";               //Added in 2024
    CMS13_pp_h_phi3phi3_bbmumu = readTable(lowC03.str(),95,2);

    lowC04 << tablepath << "CMS-HIG-19-007_4a.dat";               //Added in 2024
    CMS13_pp_h_phi23Z_mumull = readTable(lowC04.str(),352,2);

    lowC05 << tablepath << "CMS-HIG-19-007_5a.dat";               //Added in 2024
    CMS13_pp_h_phi23phi23_mumumumu = readTable(lowC05.str(),462,2);

    lowC06 << tablepath << "CMS-HIG-21-003_6.dat";               //Added in 2024
    CMS13_pp_h_phi3phi3_gagagaga = readTable(lowC06.str(),95,2);

    lowC07 << tablepath << "CMS-HIG-18-006_10.dat";               //Added in 2024
    CMS13_pp_h_phi3phi3_tautautautau = readTable(lowC07.str(),12,2);

    lowC08 << tablepath << "CMS-PAS-HIG-20-002_5b.dat";               //Added in 2024
    CMS13_pp_phi2_gaga = readTable(lowC08.str(),81,2);

    lowC09 << tablepath << "CMS-HIG-17-014_4.dat";               //Added in 2024
    CMS13_pp_bbphi3_bbtautau = readTable(lowC09.str(),10,2);

    lowC10 << tablepath << "CMS-EXO-21-018_14d.dat";               //Added in 2024
    CMS13_pp_ttphi3_ttmumu = readTable(lowC10.str(),220,2);

    lowC11 << tablepath << "CMS-EXO-21-018_14f.dat";               //Added in 2024
    CMS13_pp_ttphi3_tttautau = readTable(lowC11.str(),18,2);

    lowC12 << tablepath << "CMS-HIG-18-026_3c.dat";               //Added in 2024
    CMS13_Vh_h_phi3phi3_bbbb = readTable(lowC12.str(),10,2);

    lowA01 << tablepath << "ATLAS_CERN-EP-2021-157_9.dat";               //Added in 2024
    ATLAS13_pp_h_phi3phi3_bbmumu = readTable(lowA01.str(),47,2);

    lowA02 << tablepath << "ATLAS_CERN-EP-2021-193_14b.dat";               //Added in 2024
    ATLAS13_gg_h_phi23phi23_mumumumu = readTable(lowA02.str(),110,2);

    lowA03 << tablepath << "ATLAS_CERN-EP-2021-193_17b.dat";               //Added in 2024
    ATLAS13_gg_h_phi23Z_mumull = readTable(lowA03.str(),16,2);

    lowA04 << tablepath << "ATLAS_CERN-EP-2018-128_9c.dat";               //Added in 2024
    ATLAS13_Vh_h_phi23phi23_bbbb = readTable(lowA04.str(),41,2);

    lowA05 << tablepath << "ATLAS_CERN-EP-2020-067_10.dat";               //Added in 2024
    ATLAS13_Zh_h_phi23phi23_bbbb = readTable(lowA05.str(),7,2);

    lowA06 << tablepath << "ATLAS_CERN-EP-2018-153_6a.dat";               //Added in 2024
    ATLAS13_pp_h_phi23phi23_bbmumu_old = readTable(lowA06.str(),101,2);

    lowA07 << tablepath << "ATLAS_CERN-EP-2017-295_2.dat";               //Added in 2024
    ATLAS13_pp_h_phi23phi23_gagagg = readTable(lowA07.str(),13,2);

    lowA08 << tablepath << "ATLAS-CONF-2023-035_7b.dat";               //Added in 2024
    ATLAS13_pp_phi2_gaga_low = readTable(lowA08.str(),89,2);

    lowA09 << tablepath << "ATLAS_CERN-EP-2023-070_5a.dat";               //Added in 2024
    ATLAS13_pp_ttphi3_ttmumu = readTable(lowA09.str(),37,2);
    
    lowA10 << tablepath << "ATLAS_CERN-EP-2024-235_7.dat";               //Added in 2024
    ATLAS13_gg_phi3_tautau_low = readTable(lowA10.str(),12,2);

    lowA11 << tablepath << "ATLAS_CERN-EP-2023-202_5d.dat";               //Added in 2024
    ATLAS13_pp_h_phi3phi3_gagagaga = readTable(lowA11.str(),53,2);

    lowA12 << tablepath << "ATLAS_CERN-EP-2024-164_10a.dat";              //Added in 2024
    ATLAS13_pp_h_phi3phi3_bbtautau = readTable(lowA12.str(),8,2);

    lowA801 << tablepath << "ATLAS_CERN-PH-EP-2015-187_4b.dat";               //Added in 2024
    ATLAS8_pp_h_phi3phi3_gagagaga = readTable(lowA801.str(),105,2);

    lowA802 << tablepath << "ATLAS_CERN-PH-EP-2015-057_6a.dat";               //Added in 2024
    ATLAS8_gg_h_phi3phi3_tautautautau = readTable(lowA802.str(),101,2);

    lowC801 << tablepath << "CMS-HIG-16-015_4.dat";               //Added in 2024
    CMS8_pp_h_phi3phi3_tautautautau = readTable(lowC801.str(),6,2);

    lowC802 << tablepath << "CMS-HIG-16-015_5.dat";               //Added in 2024
    CMS8_pp_h_phi3phi3_bbmumu = readTable(lowC802.str(),76,2);

    lowC803 << tablepath << "CMS-HIG-16-015_6.dat";               //Added in 2024
    CMS8_pp_h_phi3phi3_mumutautau = readTable(lowC803.str(),96,2);

    lowC804 << tablepath << "CMS-HIG-17-013_4a.dat";               //Added in 2024
    CMS8_pp_phi2_gaga = readTable(lowC804.str(),300,2);

    lowC805 << tablepath << "CMS-HIG-14-033_3.dat";               //Added in 2024
    CMS8_pp_bbphi3_bbtautau = readTable(lowC805.str(),12,2);

    lowC806 << tablepath << "CMS-HIG-15-009_5b.dat";               //Added in 2024
    CMS8_pp_bbphi3_bbmumu = readTable(lowC806.str(),36,2);

    low209a << tablepath << "LEP_LHWG-NOTE-2001-08_4.dat";         //Added in 2024
    LEP209_phi2Z_gagaZ = readTable(low209a.str(),49,2);

    low209b << tablepath << "LEP_CERN-PH-EP-2006-001_t14b.dat";    //Added in 2024
    LEP209_phi2Z_bbZ = readTable(low209b.str(),217,2);

    low209c << tablepath << "LEP_CERN-PH-EP-2006-001_t14c.dat";    //Added in 2024
    LEP209_phi2Z_tautauZ = readTable(low209c.str(),233,2);

    low209d << tablepath << "LEP_CERN-PH-EP-2006-001_t18.dat";     //Added in 2024
    LEP209_phiEphi3_bbbb = readTable(low209d.str(),324,3);

    low209e << tablepath << "LEP_CERN-PH-EP-2006-001_t19.dat";     //Added in 2024
    LEP209_phiEphi3_tautautautau = readTable(low209e.str(),365,3);

    lowHpC801 << tablepath << "CMS-HIG-13-035_t4.dat";             //Added in 2024
    CMS8_t_Hpb_csb = readTable(lowHpC801.str(),7,2);

    lowHpC802 << tablepath << "CMS-HIG-14-023_t10.dat";            //Added in 2024
    CMS8_t_Hpb_taunub = readTable(lowHpC802.str(),8,2);

    lowHpC803 << tablepath << "CMS-HIG-16-030_3.dat";              //Added in 2024
    CMS8_t_Hpb_cbb = readTable(lowHpC803.str(),7,2);

    lowHpC1301 << tablepath << "CMS-HIG-18-020_2a.dat";            //Added in 2024
    CMS13_t_Hpb_WAb_Wmumub = readTable(lowHpC1301.str(),95,2);

    lowHpC1302 << tablepath << "CMS-HIG-18-021_4c.dat";            //Added in 2024
    CMS13_t_Hpb_csb = readTable(lowHpC1302.str(),8,2);

    lowHpA801 << tablepath << "ATLAS_CERN-PH-EP-2014-274_7a.dat";  //Added in 2024
    ATLAS8_t_Hpb_taunub = readTable(lowHpA801.str(),9,2);

    lowHpA1301 << tablepath << "ATLAS_CERN-EP-2022-207_8.dat";     //Added in 2024
    ATLAS13_t_Hpb_cbb = readTable(lowHpA1301.str(),11,2);

    lowHpA1302 << tablepath << "ATLAS_CERN-EP-2023-070_5b.dat";    //Added in 2024
    ATLAS13_t_Hpb_WAb_Wmumub = readTable(lowHpA1302.str(),31,2);

    lowHpA1303 << tablepath << "ATLAS_CERN-EP-2024-185_6.dat";    //Added in 2024
    ATLAS13_t_Hpb_csb = readTable(lowHpA1303.str(),12,2);

    lowHp209a << tablepath << "LEP_CERN-PH-EP-2012-369_4a.dat";    //Added in 2024
    LEP209_HpHm_taunutaunu = readTable(lowHp209a.str(),105,2);

    lowHp209b << tablepath << "LEP_CERN-PH-EP-2012-369_4d.dat";    //Added in 2024
    LEP209_HpHm_qqqq = readTable(lowHp209b.str(),105,2);

    lowHp209c << tablepath << "LEP_CERN-PH-EP-2008-016_7b.dat";    //Added in 2024
    OPAL209_HpHm_qqtaunu = readTable(lowHp209c.str(),87,2);

    lowHp172a << tablepath << "LEP_CERN-PPE-97-168_3b_tnqq.dat";   //Added in 2024
    OPAL172_HpHm_qqtaunu = readTable(lowHp172a.str(),41,2);

    lowHp209d << tablepath << "LEP_CERN-PH-EP-2008-016_10a.dat";   //Added in 2024
    OPAL209_HpHm_AWstAWst_bbbb = readTable(lowHp209d.str(),2835,3);

    lowHp209e << tablepath << "LEP_CERN-PH-EP-2008-016_10b.dat";   //Added in 2024
    OPAL209_HpHm_taunuAWst_taunubb = readTable(lowHp209e.str(),2835,3);

    susyHpA01 << tablepath << "ATLAS_CERN-EP-2019-191_7b.dat";     //Added in 2024
    ATLAS13_pp_HpHm_taunutaunu = readTable(susyHpA01.str(),19,2);

    susyHpA02 << tablepath << "ATLAS_CERN-EP-2022-132_6a.dat";     //Added in 2024
    ATLAS13_pp_HpHm_munumunu = readTable(susyHpA02.str(),11,2);

    susyHpC01 << tablepath << "CMS-SUS-21-001_5a.dat";             //Added in 2024
    CMS13_pp_HpHm_taunutaunu = readTable(susyHpC01.str(),16,2);

    susyHpC02 << tablepath << "CMS-SUS-21-008_16a.dat";            //Added in 2024
    CMS13_pp_HpHm_munumunu = readTable(susyHpC02.str(),22,2);

    thint01 << tablepath << "integral_x2_1mx_G_values_log.dat";
    integral_x2_1mx_G_log = readTable(thint01.str(),62500,3);

    thint02 << tablepath << "integral_x2_1px_G_values_log.dat";
    integral_x2_1px_G_log = readTable(thint02.str(),62500,3);

    thint03 << tablepath << "integral_x2_G_values_log.dat";
    integral_x2_G_log = readTable(thint03.str(),62500,3);

    thint04 << tablepath << "integral_x_1mx2_G_values_log.dat";
    integral_x_1mx2_G_log = readTable(thint04.str(),62500,3);

    thint05 << tablepath << "integral_x_1mx_1px_G_values_log.dat";
    integral_x_1mx_1px_G_log = readTable(thint05.str(),62500,3);

    thint06 << tablepath << "integral_x2_1mx_G_variable_set_1_values_log.dat";
    integral_x2_1mx_G_variable_set_1_log = readTable(thint06.str(),2000,2);

    thint07 << tablepath << "integral_x2_G_variable_set_1_values_log.dat";
    integral_x2_G_variable_set_1_log = readTable(thint07.str(),2000,2);

    thint08 << tablepath << "integral_x_1mx2_G_variable_set_0_values_log.dat";
    integral_x_1mx2_G_variable_set_0_log = readTable(thint08.str(),2000,2);

    thint09 << tablepath << "integral_perturb_ATHDM_values_log.dat";
    integral_perturb_ATHDM_values_log = readTable(thint09.str(),2000,2);

    thint10 << tablepath << "integral_Hp_A_Wstar.dat";
    integral_Hp_A_Wstar = readTable(thint10.str(),2835,3);

    //std::cout<< CMS13_pp_Hpm_taunu<<std::endl;

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
        if (mass>=10. && mass <=2000.) {
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
        if (mass>=10. && mass <=2000.) {
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
        if (mass>=10. && mass <=2000.) {
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
        if (mass>=10. && mass <=2000.) {
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
        if (mass>=10. && mass <=2000.) {
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
        if (mass>=10. && mass <=2000.) {
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
        if (mass>=10. && mass <=2000.) {
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
        if (mass>=10. && mass <=2000.) {
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
        if (mass>=10. && mass <=2000.) {
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
        if (mass>=10. && mass <=2000.) {
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
        if (mass>=10. && mass <=2000.) {
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
        if (mass>=10. && mass <=2000.) {
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



double GeneralTHDMcache::ip_csr_ggH_tc_8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_csr_ggH_tc_8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_csr_ggH_tc_8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (csrH_top_charm_8,mass);
        CacheShiftReal(ip_csr_ggH_tc_8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_csr_ggH_tc_13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_csr_ggH_tc_13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_csr_ggH_tc_13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (csrH_top_charm_13,mass);
        CacheShiftReal(ip_csr_ggH_tc_13_cache, NumPar, params, newResult);
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



double GeneralTHDMcache::ip_csr_ggA_tc_8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_csr_ggA_tc_8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_csr_ggA_tc_8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (csrA_top_charm_8,mass);
        CacheShiftReal(ip_csr_ggA_tc_8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_csr_ggA_tc_13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_csr_ggA_tc_13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_csr_ggA_tc_13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (csrA_top_charm_13,mass);
        CacheShiftReal(ip_csr_ggA_tc_13_cache, NumPar, params, newResult);
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


double GeneralTHDMcache::ip_csr_HpHm_pp_13_LO(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_csr_HpHm_pp_13_LO_cache, NumPar, params);
    if (i>=0) {
        return ( ip_csr_HpHm_pp_13_LO_cache[NumPar][i] );
    } else {
        double newResult = interpolate (csrHpHm_pp_13_LO,mass);
        CacheShiftReal(ip_csr_HpHm_pp_13_LO_cache, NumPar, params, newResult);
        return newResult;
    }
}


double GeneralTHDMcache::ip_csr_HpHm_pp_13_NLO(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_csr_HpHm_pp_13_NLO_cache, NumPar, params);
    if (i>=0) {
        return ( ip_csr_HpHm_pp_13_NLO_cache[NumPar][i] );
    } else {
        double newResult = interpolateNU (csrHpHm_pp_13_NLO,mass);
        CacheShiftReal(ip_csr_HpHm_pp_13_NLO_cache, NumPar, params, newResult);
        return newResult;
    }
}


//ATLAS13_bb_phi_bb
double GeneralTHDMcache::ip_ex_bb_phi_bb_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_bb_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phi_bb_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS13_bb_phi_bb,mass);
        CacheShiftReal(ip_ex_bb_phi_bb_ATLAS13_cache, NumPar, params, newResult);
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

//CMS13_tt_phi3_tt

double GeneralTHDMcache::ip_ex_tt_phi2_tt_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_tt_phi2_tt_CMS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_tt_phi2_tt_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS13_tt_phi2_tt,mass);
        CacheShiftReal(ip_ex_tt_phi2_tt_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}


double GeneralTHDMcache::ip_ex_tt_phi3_tt_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_tt_phi3_tt_CMS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_tt_phi3_tt_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS13_tt_phi3_tt,mass);
        CacheShiftReal(ip_ex_tt_phi3_tt_CMS13_cache, NumPar, params, newResult);
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








double GeneralTHDMcache::ip_ex_pp_phi2_bb_light_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi2_bb_light_CMS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_phi2_bb_light_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS13_pp_phi2_bb_light,mass);
        CacheShiftReal(ip_ex_pp_phi2_bb_light_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}




double GeneralTHDMcache::ip_ex_pp_phi3_bb_light_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi3_bb_light_CMS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_phi3_bb_light_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS13_pp_phi3_bb_light,mass);
        CacheShiftReal(ip_ex_pp_phi3_bb_light_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}


















double GeneralTHDMcache::ip_ex_gg_phi_mumu_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_mumu_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_mumu_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS8_gg_phi_mumu,mass);
        CacheShiftReal(ip_ex_gg_phi_mumu_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}


double GeneralTHDMcache::ip_ex_bb_phi_mumu_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_mumu_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phi_mumu_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS8_bb_phi_mumu,mass);
        CacheShiftReal(ip_ex_bb_phi_mumu_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}


double GeneralTHDMcache::ip_ex_gg_phi_mumu_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_mumu_CMS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_mumu_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS13_gg_phi_mumu,mass);
        CacheShiftReal(ip_ex_gg_phi_mumu_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_ex_bb_phi_mumu_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_mumu_CMS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phi_mumu_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS13_bb_phi_mumu,mass);
        CacheShiftReal(ip_ex_bb_phi_mumu_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}


double GeneralTHDMcache::ip_ex_gg_phi_mumu_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_mumu_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_mumu_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS13_gg_phi_mumu,mass);
        CacheShiftReal(ip_ex_gg_phi_mumu_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_ex_bb_phi_mumu_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_mumu_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phi_mumu_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS13_bb_phi_mumu,mass);
        CacheShiftReal(ip_ex_bb_phi_mumu_ATLAS13_cache, NumPar, params, newResult);
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






double GeneralTHDMcache::ip_ex_gg_phi_WW_heavy_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_WW_heavy_CMS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_WW_heavy_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS13_gg_phi_WW_heavy,mass);
        CacheShiftReal(ip_ex_gg_phi_WW_heavy_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}





double GeneralTHDMcache::ip_ex_VV_phi_WW_heavy_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VV_phi_WW_heavy_CMS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_VV_phi_WW_heavy_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS13_VV_phi_WW_heavy,mass);
        CacheShiftReal(ip_ex_VV_phi_WW_heavy_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}











double GeneralTHDMcache::ip_ex_gg_phi_WW_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_WW_CMS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_WW_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS13_gg_phi_WW,mass);
        CacheShiftReal(ip_ex_gg_phi_WW_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}





double GeneralTHDMcache::ip_ex_VV_phi_WW_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VV_phi_WW_CMS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_VV_phi_WW_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS13_VV_phi_WW,mass);
        CacheShiftReal(ip_ex_VV_phi_WW_CMS13_cache, NumPar, params, newResult);
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

double GeneralTHDMcache::ip_ex_pp_phi_VV_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_VV_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_phi_VV_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS8_pp_phi_VV,mass);
        CacheShiftReal(ip_ex_pp_phi_VV_CMS8_cache, NumPar, params, newResult);
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







double GeneralTHDMcache::ip_ex_gg_phi_VV_llqq_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_VV_llqq_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_VV_llqq_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS13_gg_phi_VV_llqq,mass);
        CacheShiftReal(ip_ex_gg_phi_VV_llqq_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}


double GeneralTHDMcache::ip_ex_VV_phi_VV_llqq_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_VV_phi_VV_llqq_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_VV_phi_VV_llqq_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS13_VV_phi_VV_llqq,mass);
        CacheShiftReal(ip_ex_VV_phi_VV_llqq_ATLAS13_cache, NumPar, params, newResult);
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
        double newResult = interpolate(ATLAS13_pp_phi_phi1phi1_bbgaga,mass);
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
        double newResult = interpolate(CMS13_pp_phi_phi1phi1_bbgaga,mass);
        CacheShiftReal(ip_ex_pp_phi_phi1phi1_bbgaga_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}


/*
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
*/



double GeneralTHDMcache::ip_ex_pp_phi_phi1phi1_bbtautau_1_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_phi1phi1_bbtautau_1_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_phi1phi1_bbtautau_1_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_phi_phi1phi1_bbtautau_1,mass);
        CacheShiftReal(ip_ex_pp_phi_phi1phi1_bbtautau_1_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}


double GeneralTHDMcache::ip_ex_pp_phi_phi1phi1_bbtautau_2_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_phi1phi1_bbtautau_2_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_phi1phi1_bbtautau_2_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_phi_phi1phi1_bbtautau_2,mass);
        CacheShiftReal(ip_ex_pp_phi_phi1phi1_bbtautau_2_ATLAS13_cache, NumPar, params, newResult);
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



double GeneralTHDMcache::ip_ex_pp_phi_phi1phi1_4WOr2W2tauOr4tau_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_phi1phi1_4WOr2W2tauOr4tau_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_phi1phi1_4WOr2W2tauOr4tau_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_phi_phi1phi1_4WOr2W2tauOr4tau,mass);
        CacheShiftReal(ip_ex_pp_phi_phi1phi1_4WOr2W2tauOr4tau_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}




//CMS13_pp_phi_phi1phi1_bbWW_qqlnu
double GeneralTHDMcache::ip_ex_pp_phi_phi1phi1_bbWW_qqlnu_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_phi1phi1_bbWW_qqlnu_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_phi1phi1_bbWW_qqlnu_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_phi_phi1phi1_bbWW_qqlnu,mass);
        CacheShiftReal(ip_ex_pp_phi_phi1phi1_bbWW_qqlnu_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}
















//CMS13_pp_phi_phi1phi1_bbZZ_lljj
double GeneralTHDMcache::ip_ex_pp_phi_phi1phi1_bbZZ_lljj_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_phi1phi1_bbZZ_lljj_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_phi1phi1_bbZZ_lljj_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_phi_phi1phi1_bbZZ_lljj,mass);
        CacheShiftReal(ip_ex_pp_phi_phi1phi1_bbZZ_lljj_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}




//CMS13_pp_phi_phi1phi1_bbZZ_llnunu
double GeneralTHDMcache::ip_ex_pp_phi_phi1phi1_bbZZ_llnunu_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_phi1phi1_bbZZ_llnunu_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_phi1phi1_bbZZ_llnunu_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_phi_phi1phi1_bbZZ_llnunu,mass);
        CacheShiftReal(ip_ex_pp_phi_phi1phi1_bbZZ_llnunu_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}




//CMS13_pp_phi_phi1phi1_bbWWorbbtautau
double GeneralTHDMcache::ip_ex_pp_phi_phi1phi1_bbWWorbbtautau_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_phi1phi1_bbWWorbbtautau_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_phi1phi1_bbWWorbbtautau_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_phi_phi1phi1_bbWWorbbtautau,mass);
        CacheShiftReal(ip_ex_pp_phi_phi1phi1_bbWWorbbtautau_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}


//CMS13_pp_phi_phi1phi1_bbWWorbbtautau_low_masses
double GeneralTHDMcache::ip_ex_pp_phi_phi1phi1_bbWWorbbtautau_low_masses_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_phi1phi1_bbWWorbbtautau_low_masses_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_phi1phi1_bbWWorbbtautau_low_masses_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_phi_phi1phi1_bbWWorbbtautau_low_masses,mass);
        CacheShiftReal(ip_ex_pp_phi_phi1phi1_bbWWorbbtautau_low_masses_CMS13_cache, NumPar, params, newResult);
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




double GeneralTHDMcache::ip_ex_gg_phi_phi1Z_tautaull_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_phi1Z_tautaull_CMS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_phi1Z_tautaull_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS13_gg_phi_phi1Z_tautaull,mass);
        CacheShiftReal(ip_ex_gg_phi_phi1Z_tautaull_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}




//ATLAS13_bb_phi_phi1Z_tautaull


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
        double newResult = interpolate2DtriangularData(ATLAS13_gg_phii_phijZ_bbZ, mj, mi);
        
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
        double newResult = interpolate2DtriangularData(ATLAS13_bb_phii_phijZ_bbZ, mj, mi);
        CacheShiftReal(ip_ex_bb_phii_phijZ_bbZ_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}




double GeneralTHDMcache::ip_ex_gg_phii_phijZ_WWZ_ATLAS13(double mj, double mi){
    int NumPar = 2;
    double params[] = {mj, mi};

    int i = CacheCheckReal(ip_ex_gg_phii_phijZ_WWZ_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phii_phijZ_WWZ_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate2DtriangularData(ATLAS13_gg_phii_phijZ_WWZ, mj, mi);
        CacheShiftReal(ip_ex_gg_phii_phijZ_WWZ_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}




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



double GeneralTHDMcache::ip_ex_pp_Hpm_tb_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_Hpm_tb_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_Hpm_tb_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_Hpm_tb,mass);
        CacheShiftReal(ip_ex_pp_Hpm_tb_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}


double GeneralTHDMcache::ip_low_pp_h_phi3phi3_mumutautau_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_low_pp_h_phi3phi3_mumutautau_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_pp_h_phi3phi3_mumutautau_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_h_phi3phi3_mumutautau,mass);
        CacheShiftReal(ip_low_pp_h_phi3phi3_mumutautau_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_pp_h_phi3phi3_bbtautau_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_low_pp_h_phi3phi3_bbtautau_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_pp_h_phi3phi3_bbtautau_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolateNU (CMS13_pp_h_phi3phi3_bbtautau,mass);
        CacheShiftReal(ip_low_pp_h_phi3phi3_bbtautau_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_pp_h_phi3phi3_bbmumu_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_low_pp_h_phi3phi3_bbmumu_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_pp_h_phi3phi3_bbmumu_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_h_phi3phi3_bbmumu,mass);
        CacheShiftReal(ip_low_pp_h_phi3phi3_bbmumu_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_pp_h_phi23Z_mumull_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_low_pp_h_phi23Z_mumull_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_pp_h_phi23Z_mumull_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolateNU (CMS13_pp_h_phi23Z_mumull,mass);
        CacheShiftReal(ip_low_pp_h_phi23Z_mumull_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_pp_h_phi23phi23_mumumumu_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_low_pp_h_phi23phi23_mumumumu_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_pp_h_phi23phi23_mumumumu_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolateNU (CMS13_pp_h_phi23phi23_mumumumu,mass);
        CacheShiftReal(ip_low_pp_h_phi23phi23_mumumumu_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_pp_h_phi3phi3_gagagaga_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_low_pp_h_phi3phi3_gagagaga_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_pp_h_phi3phi3_gagagaga_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_h_phi3phi3_gagagaga,mass);
        CacheShiftReal(ip_low_pp_h_phi3phi3_gagagaga_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_pp_h_phi3phi3_tautautautau_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_low_pp_h_phi3phi3_tautautautau_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_pp_h_phi3phi3_tautautautau_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_h_phi3phi3_tautautautau,mass);
        CacheShiftReal(ip_low_pp_h_phi3phi3_tautautautau_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_pp_ttphi3_ttmumu_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_low_pp_ttphi3_ttmumu_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_pp_ttphi3_ttmumu_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolateNU (CMS13_pp_ttphi3_ttmumu,mass);
        CacheShiftReal(ip_low_pp_ttphi3_ttmumu_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_pp_ttphi3_tttautau_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_low_pp_ttphi3_tttautau_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_pp_ttphi3_tttautau_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolateNU (CMS13_pp_ttphi3_tttautau,mass);
        CacheShiftReal(ip_low_pp_ttphi3_tttautau_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_Vh_h_phi3phi3_bbbb_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_low_Vh_h_phi3phi3_bbbb_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_Vh_h_phi3phi3_bbbb_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_Vh_h_phi3phi3_bbbb,mass);
        CacheShiftReal(ip_low_Vh_h_phi3phi3_bbbb_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_pp_h_phi3phi3_bbmumu_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_pp_h_phi3phi3_bbmumu_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_pp_h_phi3phi3_bbmumu_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_h_phi3phi3_bbmumu,mass);
        CacheShiftReal(ip_low_pp_h_phi3phi3_bbmumu_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_gg_h_phi23phi23_mumumumu_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_gg_h_phi23phi23_mumumumu_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_gg_h_phi23phi23_mumumumu_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolateNU (ATLAS13_gg_h_phi23phi23_mumumumu,mass);
        CacheShiftReal(ip_low_gg_h_phi23phi23_mumumumu_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_gg_h_phi23Z_mumull_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_gg_h_phi23Z_mumull_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_gg_h_phi23Z_mumull_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_gg_h_phi23Z_mumull,mass);
        CacheShiftReal(ip_low_gg_h_phi23Z_mumull_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_Vh_h_phi23phi23_bbbb_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_Vh_h_phi23phi23_bbbb_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_Vh_h_phi23phi23_bbbb_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_Vh_h_phi23phi23_bbbb,mass);
        CacheShiftReal(ip_low_Vh_h_phi23phi23_bbbb_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_Zh_h_phi23phi23_bbbb_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_Zh_h_phi23phi23_bbbb_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_Zh_h_phi23phi23_bbbb_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_Zh_h_phi23phi23_bbbb,mass);
        CacheShiftReal(ip_low_Zh_h_phi23phi23_bbbb_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_pp_h_phi23phi23_bbmumu_ATLAS13_old(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_pp_h_phi23phi23_bbmumu_ATLAS13_old_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_pp_h_phi23phi23_bbmumu_ATLAS13_old_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_h_phi23phi23_bbmumu_old,mass);
        CacheShiftReal(ip_low_pp_h_phi23phi23_bbmumu_ATLAS13_old_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_pp_h_phi23phi23_gagagg_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_pp_h_phi23phi23_gagagg_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_pp_h_phi23phi23_gagagg_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolateNU (ATLAS13_pp_h_phi23phi23_gagagg,mass);
        CacheShiftReal(ip_low_pp_h_phi23phi23_gagagg_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_gg_phi3_tautau_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_gg_phi3_tautau_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_gg_phi3_tautau_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolateNU (ATLAS13_gg_phi3_tautau_low,mass);
        CacheShiftReal(ip_low_gg_phi3_tautau_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_pp_h_phi3phi3_gagagaga_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_pp_h_phi3phi3_gagagaga_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_pp_h_phi3phi3_gagagaga_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_h_phi3phi3_gagagaga,mass);
        CacheShiftReal(ip_low_pp_h_phi3phi3_gagagaga_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_pp_h_phi3phi3_bbtautau_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_pp_h_phi3phi3_bbtautau_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_pp_h_phi3phi3_bbtautau_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolateNU (ATLAS13_pp_h_phi3phi3_bbtautau,mass);
        CacheShiftReal(ip_low_pp_h_phi3phi3_bbtautau_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_pp_h_phi3phi3_gagagaga_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_pp_h_phi3phi3_gagagaga_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_pp_h_phi3phi3_gagagaga_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS8_pp_h_phi3phi3_gagagaga,mass);
        CacheShiftReal(ip_low_pp_h_phi3phi3_gagagaga_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_gg_h_phi3phi3_tautautautau_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_gg_h_phi3phi3_tautautautau_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_gg_h_phi3phi3_tautautautau_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolateNU (ATLAS8_gg_h_phi3phi3_tautautautau,mass);
        CacheShiftReal(ip_low_gg_h_phi3phi3_tautautautau_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_pp_h_phi3phi3_tautautautau_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_pp_h_phi3phi3_tautautautau_CMS8_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_pp_h_phi3phi3_tautautautau_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_pp_h_phi3phi3_tautautautau,mass);
        CacheShiftReal(ip_low_pp_h_phi3phi3_tautautautau_CMS8_cache, NumPar, params, newResult);
                return newResult;
    }
}

double GeneralTHDMcache::ip_low_pp_h_phi3phi3_bbmumu_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_pp_h_phi3phi3_bbmumu_CMS8_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_pp_h_phi3phi3_bbmumu_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_pp_h_phi3phi3_bbmumu,mass);
        CacheShiftReal(ip_low_pp_h_phi3phi3_bbmumu_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_pp_h_phi3phi3_mumutautau_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_pp_h_phi3phi3_mumutautau_CMS8_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_pp_h_phi3phi3_mumutautau_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_pp_h_phi3phi3_mumutautau,mass);
        CacheShiftReal(ip_low_pp_h_phi3phi3_mumutautau_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}


double GeneralTHDMcache::ip_low_pp_phi2_gaga_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_pp_phi2_gaga_CMS8_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_pp_phi2_gaga_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS8_pp_phi2_gaga,mass);
        CacheShiftReal(ip_low_pp_phi2_gaga_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_pp_phi2_gaga_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_pp_phi2_gaga_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_pp_phi2_gaga_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_phi2_gaga,mass);
        CacheShiftReal(ip_low_pp_phi2_gaga_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_pp_phi2_gaga_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_pp_phi2_gaga_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_pp_phi2_gaga_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_phi2_gaga_low,mass);
        CacheShiftReal(ip_low_pp_phi2_gaga_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_pp_bbphi3_bbtautau_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_pp_bbphi3_bbtautau_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_pp_bbphi3_bbtautau_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_bbphi3_bbtautau,mass);
        CacheShiftReal(ip_low_pp_bbphi3_bbtautau_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_pp_bbphi3_bbtautau_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_pp_bbphi3_bbtautau_CMS8_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_pp_bbphi3_bbtautau_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_pp_bbphi3_bbtautau,mass);
        CacheShiftReal(ip_low_pp_bbphi3_bbtautau_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_pp_ttphi3_ttmumu_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_pp_ttphi3_ttmumu_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_pp_ttphi3_ttmumu_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolateNU (ATLAS13_pp_ttphi3_ttmumu,mass);
        CacheShiftReal(ip_low_pp_ttphi3_ttmumu_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_pp_bbphi3_bbmumu_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_pp_bbphi3_bbmumu_CMS8_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_pp_bbphi3_bbmumu_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_pp_bbphi3_bbmumu,mass);
        CacheShiftReal(ip_low_pp_bbphi3_bbmumu_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_phi2Z_gagaZ_LEP209(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_phi2Z_gagaZ_LEP209_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_phi2Z_gagaZ_LEP209_cache[NumPar][i] );
    } else {
        double newResult = interpolate (LEP209_phi2Z_gagaZ,mass);
        CacheShiftReal(ip_low_phi2Z_gagaZ_LEP209_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_phi2Z_bbZ_LEP209(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_phi2Z_bbZ_LEP209_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_phi2Z_bbZ_LEP209_cache[NumPar][i] );
    } else {
        double newResult = interpolate (LEP209_phi2Z_bbZ,mass);
        CacheShiftReal(ip_low_phi2Z_bbZ_LEP209_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_phi2Z_tautauZ_LEP209(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_phi2Z_tautauZ_LEP209_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_phi2Z_tautauZ_LEP209_cache[NumPar][i] );
    } else {
        double newResult = interpolate (LEP209_phi2Z_tautauZ,mass);
        CacheShiftReal(ip_low_phi2Z_tautauZ_LEP209_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_phiEphi3_bbbb_LEP209(double mi, double mj, double di, double dj){
    int NumPar = 2;
    double params[] = {mi, mj};
    int i = CacheCheckReal(ip_low_phiEphi3_bbbb_LEP209_cache, NumPar, params);
    if (i>=0) {
        return ( ip_low_phiEphi3_bbbb_LEP209_cache[NumPar][i] );
    } else {
        double newResult = interpolate2Dsteps(LEP209_phiEphi3_bbbb, mi, mj, di, dj);
        CacheShiftReal(ip_low_phiEphi3_bbbb_LEP209_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_phiEphi3_tautautautau_LEP209(double mi, double mj, double di, double dj){
    int NumPar = 2;
    double params[] = {mi, mj};
    int i = CacheCheckReal(ip_low_phiEphi3_tautautautau_LEP209_cache, NumPar, params);
    if (i>=0) {
        return ( ip_low_phiEphi3_tautautautau_LEP209_cache[NumPar][i] );
    } else {
        double newResult = interpolate2Dsteps (LEP209_phiEphi3_tautautautau, mi, mj, di, dj);
        CacheShiftReal(ip_low_phiEphi3_tautautautau_LEP209_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_t_Hpb_csb_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_t_Hpb_csb_CMS8_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_t_Hpb_csb_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolateNU (CMS8_t_Hpb_csb,mass);
        CacheShiftReal(ip_low_t_Hpb_csb_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_t_Hpb_taunub_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_t_Hpb_taunub_CMS8_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_t_Hpb_taunub_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolateNU (CMS8_t_Hpb_taunub,mass);
        CacheShiftReal(ip_low_t_Hpb_taunub_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_t_Hpb_cbb_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_t_Hpb_cbb_CMS8_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_t_Hpb_cbb_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_t_Hpb_cbb,mass);
        CacheShiftReal(ip_low_t_Hpb_cbb_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_t_Hpb_WAb_Wmumub_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_t_Hpb_WAb_Wmumub_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_t_Hpb_WAb_Wmumub_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolateNU (CMS13_t_Hpb_WAb_Wmumub,mass);
        CacheShiftReal(ip_low_t_Hpb_WAb_Wmumub_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_t_Hpb_csb_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_t_Hpb_csb_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_t_Hpb_csb_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolateNU (CMS13_t_Hpb_csb,mass);
        CacheShiftReal(ip_low_t_Hpb_csb_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_t_Hpb_taunub_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_t_Hpb_taunub_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_t_Hpb_taunub_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS8_t_Hpb_taunub,mass);
        CacheShiftReal(ip_low_t_Hpb_taunub_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_t_Hpb_cbb_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_t_Hpb_cbb_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_t_Hpb_cbb_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_t_Hpb_cbb,mass);
        CacheShiftReal(ip_low_t_Hpb_cbb_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_t_Hpb_WAb_Wmumub_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_t_Hpb_WAb_Wmumub_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_t_Hpb_WAb_Wmumub_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolateNU (ATLAS13_t_Hpb_WAb_Wmumub,mass);
        CacheShiftReal(ip_low_t_Hpb_WAb_Wmumub_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_t_Hpb_csb_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_t_Hpb_csb_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_t_Hpb_csb_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolateNU (ATLAS13_t_Hpb_csb,mass);
        CacheShiftReal(ip_low_t_Hpb_csb_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_HpHm_taunutaunu_LEP209(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_HpHm_taunutaunu_LEP209_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_HpHm_taunutaunu_LEP209_cache[NumPar][i] );
    } else {
        double newResult = interpolate (LEP209_HpHm_taunutaunu,mass);
        CacheShiftReal(ip_low_HpHm_taunutaunu_LEP209_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_HpHm_qqqq_LEP209(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_HpHm_qqqq_LEP209_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_HpHm_qqqq_LEP209_cache[NumPar][i] );
    } else {
        double newResult = interpolate (LEP209_HpHm_qqqq,mass);
        CacheShiftReal(ip_low_HpHm_qqqq_LEP209_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_HpHm_qqtaunu_OPAL209(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_HpHm_qqtaunu_OPAL209_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_HpHm_qqtaunu_OPAL209_cache[NumPar][i] );
    } else {
        double newResult = interpolate (OPAL209_HpHm_qqtaunu,mass);
        CacheShiftReal(ip_low_HpHm_qqtaunu_OPAL209_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_HpHm_qqtaunu_OPAL172(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_low_HpHm_qqtaunu_OPAL172_cache, NumPar, params);
    if (i>=0) {
        return(ip_low_HpHm_qqtaunu_OPAL172_cache[NumPar][i] );
    } else {
        double newResult = interpolate (OPAL172_HpHm_qqtaunu,mass);
        CacheShiftReal(ip_low_HpHm_qqtaunu_OPAL172_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_HpHm_AWstAWst_bbbb_OPAL209(double mi, double mj, double di, double dj){
    int NumPar = 2;
    double params[] = {mi, mj};
    int i = CacheCheckReal(ip_low_HpHm_AWstAWst_bbbb_OPAL209_cache, NumPar, params);
    if (i>=0) {
        return ( ip_low_HpHm_AWstAWst_bbbb_OPAL209_cache[NumPar][i] );
    } else {
        double newResult = interpolate2Dsteps (OPAL209_HpHm_AWstAWst_bbbb, mi, mj, di, dj);
        CacheShiftReal(ip_low_HpHm_AWstAWst_bbbb_OPAL209_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_low_HpHm_taunuAWst_taunubb_OPAL209(double mi, double mj, double di, double dj){
    int NumPar = 2;
    double params[] = {mi, mj};
    int i = CacheCheckReal(ip_low_HpHm_taunuAWst_taunubb_OPAL209_cache, NumPar, params);
    if (i>=0) {
        return ( ip_low_HpHm_taunuAWst_taunubb_OPAL209_cache[NumPar][i] );
    } else {
        double newResult = interpolate2Dsteps (OPAL209_HpHm_taunuAWst_taunubb, mi, mj, di, dj);
        CacheShiftReal(ip_low_HpHm_taunuAWst_taunubb_OPAL209_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_susy_pp_HpHm_taunutaunu_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_susy_pp_HpHm_taunutaunu_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_susy_pp_HpHm_taunutaunu_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_HpHm_taunutaunu,mass);
        CacheShiftReal(ip_susy_pp_HpHm_taunutaunu_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_susy_pp_HpHm_taunutaunu_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_susy_pp_HpHm_taunutaunu_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_susy_pp_HpHm_taunutaunu_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolateNU (CMS13_pp_HpHm_taunutaunu,mass);
        CacheShiftReal(ip_susy_pp_HpHm_taunutaunu_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_susy_pp_HpHm_munumunu_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_susy_pp_HpHm_munumunu_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_susy_pp_HpHm_munumunu_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolateNU (ATLAS13_pp_HpHm_munumunu,mass);
        CacheShiftReal(ip_susy_pp_HpHm_munumunu_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_susy_pp_HpHm_munumunu_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};
    int i = CacheCheckReal(ip_susy_pp_HpHm_munumunu_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_susy_pp_HpHm_munumunu_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolateNU (CMS13_pp_HpHm_munumunu,mass);
        CacheShiftReal(ip_susy_pp_HpHm_munumunu_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_susy_HpHm_munumunu_LEP208(double mass){
    // data from lepsusy.web.cern.ch/lepsusy/www/sleptons_summer04/slep_final.html
    // taking mid-value from each coloured range

    double sigxBRxBR;

    if(mass >= 45. && mass <= 65.)
        sigxBRxBR = 0.025;
    else if(mass > 65.  && mass <= 76.)
        sigxBRxBR = 0.035;
    else if(mass > 76.  && mass <= 86.)
        sigxBRxBR = 0.025;
    else if(mass > 86.  && mass <= 89.)
        sigxBRxBR = 0.035;
    else if(mass > 89.  && mass <= 96.)
        sigxBRxBR = 0.025;
    else if(mass > 96.  && mass <= 97.)
        sigxBRxBR = 0.015;
    else if(mass > 97.  && mass <= 98.)
        sigxBRxBR = 0.025;
    else if(mass > 98.  && mass <= 99.)
        sigxBRxBR = 0.015;
    else if(mass > 99.  && mass <= 101.)
        sigxBRxBR = 0.025;
    else if(mass > 101. && mass <= 102.)
        sigxBRxBR = 0.035;
    else
        sigxBRxBR = 1.e-10;

    return sigxBRxBR;
}


//Interpolation function for the integrals

double GeneralTHDMcache::ip_integral_x2_1mx_G(double wa, double wb){
    int NumPar = 2;
    double params[] = {wa, wb};
    int i = CacheCheckReal(ip_integral_x2_1mx_G_cache, NumPar, params);
    if (i>=0) {
        return(ip_integral_x2_1mx_G_cache[NumPar][i] );
    } else {
        double log10wa = log10(wa);
        double log10wb = log10(wb);
        double newResult = interpolate2Dv2(integral_x2_1mx_G_log,log10wa,log10wb);
        CacheShiftReal(ip_integral_x2_1mx_G_cache, NumPar, params, newResult);
        return newResult;
    }
}


double GeneralTHDMcache::ip_integral_x2_1px_G(double wa, double wb){
    int NumPar = 2;
    double params[] = {wa, wb};
    int i = CacheCheckReal(ip_integral_x2_1px_G_cache, NumPar, params);
    if (i>=0) {
        return(ip_integral_x2_1px_G_cache[NumPar][i] );
    } else {
        double log10wa = log10(wa);
        double log10wb = log10(wb);
        double newResult = interpolate2Dv2(integral_x2_1px_G_log,log10wa,log10wb);
        CacheShiftReal(ip_integral_x2_1px_G_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GeneralTHDMcache::ip_integral_x2_G(double wa, double wb){
    int NumPar = 2;
    double params[] = {wa, wb};
    int i = CacheCheckReal(ip_integral_x2_G_cache, NumPar, params);
    if (i>=0) {
        return(ip_integral_x2_G_cache[NumPar][i] );
    } else {
        double log10wa = log10(wa);
        double log10wb = log10(wb);
        double newResult = interpolate2Dv2(integral_x2_G_log,log10wa,log10wb);
        CacheShiftReal(ip_integral_x2_G_cache, NumPar, params, newResult);
        return newResult;
    }
}

double GeneralTHDMcache::ip_integral_x_1mx2_G(double wa, double wb){
    int NumPar = 2;
    double params[] = {wa, wb};
    int i = CacheCheckReal(ip_integral_x_1mx2_G_cache, NumPar, params);
    if (i>=0) {
        return(ip_integral_x_1mx2_G_cache[NumPar][i] );
    } else {
        double log10wa = log10(wa);
        double log10wb = log10(wb);
        double newResult = interpolate2Dv2(integral_x_1mx2_G_log,log10wa,log10wb);
        CacheShiftReal(ip_integral_x_1mx2_G_cache, NumPar, params, newResult);
        return newResult;
    }
}


double GeneralTHDMcache::ip_integral_x_1mx_1px_G(double wa, double wb){
    int NumPar = 2;
    double params[] = {wa, wb};
    int i = CacheCheckReal(ip_integral_x_1mx_1px_G_cache, NumPar, params);
    if (i>=0) {
        return(ip_integral_x_1mx_1px_G_cache[NumPar][i] );
    } else {
        double log10wa = log10(wa);
        double log10wb = log10(wb);
        double newResult = interpolate2Dv2(integral_x_1mx_1px_G_log,log10wa,log10wb);
        CacheShiftReal(ip_integral_x_1mx_1px_G_cache, NumPar, params, newResult);
        return newResult;
    }
}


double GeneralTHDMcache::ip_integral_x2_1mx_G_variable_set_1(double wb){
    int NumPar = 1;
    double params[] = {wb};
    int i = CacheCheckReal(ip_integral_x2_1mx_G_variable_set_1_cache, NumPar, params);
    if (i>=0) {
        return(ip_integral_x2_1mx_G_variable_set_1_cache[NumPar][i] );
    } else {
        double log10wb = log10(wb);
        double newResult = interpolate(integral_x2_1mx_G_variable_set_1_log,log10wb);
        CacheShiftReal(ip_integral_x2_1mx_G_variable_set_1_cache, NumPar, params, newResult);
        return newResult;
    }
}


double GeneralTHDMcache::ip_integral_x2_G_variable_set_1(double wb){
    int NumPar = 1;
    double params[] = {wb};
    int i = CacheCheckReal(ip_integral_x2_G_variable_set_1_cache, NumPar, params);
    if (i>=0) {
        return(ip_integral_x2_G_variable_set_1_cache[NumPar][i] );
    } else {
        double log10wb = log10(wb);
        double newResult = interpolate(integral_x2_G_variable_set_1_log,log10wb);
        CacheShiftReal(ip_integral_x2_G_variable_set_1_cache, NumPar, params, newResult);
        return newResult;
    }
}


double GeneralTHDMcache::ip_integral_x_1mx2_G_variable_set_0(double wb){
    int NumPar = 1;
    double params[] = {wb};
    int i = CacheCheckReal(ip_integral_x_1mx2_G_variable_set_0_cache, NumPar, params);
    if (i>=0) {
        return(ip_integral_x_1mx2_G_variable_set_0_cache[NumPar][i] );
    } else {
        double log10wb = log10(wb);
        double newResult = interpolate(integral_x_1mx2_G_variable_set_0_log,log10wb);
        CacheShiftReal(ip_integral_x_1mx2_G_variable_set_0_cache, NumPar, params, newResult);
        return newResult;
    }
}


double GeneralTHDMcache::ip_integral_perturb_ATHDM(double ratio){
    int NumPar = 1;
    double params[] = {ratio};
    int i = CacheCheckReal(ip_integral_perturb_ATHDM_cache, NumPar, params);
    if (i>=0) {
        return(ip_integral_perturb_ATHDM_cache[NumPar][i] );
    } else {
        double log10ratio = log10(ratio);
        double newResult = interpolate(integral_perturb_ATHDM_values_log,log10ratio);
        CacheShiftReal(ip_integral_perturb_ATHDM_cache, NumPar, params, newResult);
        return newResult;
    }
}


double GeneralTHDMcache::ip_integral_Hp_A_Wstar(double mi, double mj, double di, double dj){
    int NumPar = 2;
    double params[] = {mi, mj};
    int i = CacheCheckReal(ip_integral_Hp_A_Wstar_cache, NumPar, params);
    if (i>=0) {
        return ( ip_integral_Hp_A_Wstar_cache[NumPar][i] );
    } else {
        double newResult = interpolate2Dsteps (integral_Hp_A_Wstar, mi, mj, di, dj);
        CacheShiftReal(ip_integral_Hp_A_Wstar_cache, NumPar, params, newResult);
        return newResult;
    }
}


//This seems not to be used by the code. Check if this was only something included of the THDM
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
        


        //std::cout<<"\033[1;36m   sW2= \033[0m "<<  sW2  <<std::endl;
        
        //std::cout<<"\033[1;36m   cW2= \033[0m "<<  cW2  <<std::endl;
        
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


double GeneralTHDMcache::Sigma_HpHm_LEP(const double mass2, const double sqrts) const
{
    double s_LEP     = sqrts * sqrts;
    double one_MZ2_s = 1. - MZ*MZ / s_LEP;

    double cWsW   = sqrt(cW2) * sqrt(sW2);
    double ahat_e = -1. / 4. / cWsW;
    double vhat_e = (-1. + 4. * sW2) / 4. / cWsW;
    double vhat_H = (-1. + 2. * sW2) / 2. / cWsW;
    double beta_H = sqrt(std::fabs(1. - 4. * mass2 / s_LEP));

    return (2. * GF*GF * MW*MW*MW*MW * sW2*sW2 * (1. + 2. * vhat_e * vhat_H / one_MZ2_s +
           (ahat_e*ahat_e + vhat_e*vhat_e) * vhat_H*vhat_H / one_MZ2_s / one_MZ2_s) *
            beta_H*beta_H*beta_H / (3. * M_PI * s_LEP));
};


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
        //std::cout<<"\033[1;32m   c02= \033[0m "<<  c02  <<std::endl;
    return c02;
    }



    double GeneralTHDMcache::MWGTHDM(const double MW) const{
    return MW;
    }



    double GeneralTHDMcache::beta(const double mf, const double m_2) const
    {
        // We take the absolute value; this beta should always be multiplied by Heaviside theta
        return sqrt(std::fabs(1.0 - 4.0*mf*mf / m_2));
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
        return vev*(lambda3*Ri1 + Relambda7*Ri2 - Imlambda7*Ri3);
    }
       
void GeneralTHDMcache::computeSignalStrengths()
{
    //FLAG to select only the model in which all the couplings are the same (by families)
    if (!myGTHDM->getATHDMflag())
    {
        throw std::runtime_error("Signal strengths are only available in the A2HDM.");
    }

    /*complex i */
    gslpp::complex i = gslpp::complex::i();

    //fermionic couplings for phi1
    yu1 = myGTHDM->getyu1();
    yd1 = myGTHDM->getyd1();
    yl1 = myGTHDM->getyl1();

    yu1R = myGTHDM->getyu1R();
    yd1R = myGTHDM->getyd1R();
    yl1R = myGTHDM->getyl1R();
    
    
    //std::cout<<"\033[1;30m yu1 = \033[0m "<< yu1 <<std::endl;
    //std::cout<<"\033[1;30m yd1 = \033[0m "<< yd1 <<std::endl;
    //std::cout<<"\033[1;30m yl1 = \033[0m "<< yl1 <<std::endl;
    
    //std::cout<<"\033[1;30m yu1R = \033[0m "<< yu1R <<std::endl;
    //std::cout<<"\033[1;30m yd1R = \033[0m "<< yd1R <<std::endl;
    //std::cout<<"\033[1;30m yl1R = \033[0m "<< yl1R <<std::endl;

    //The ggH cross section in the SM at 8 TeV
    double SigmaggF8 = myGTHDM->computeSigmaggH(8.0);
    //The square of the top-quark contribution to the ggH cross section in the SM at 8 TeV
    double Sigmaggh_tt8 = myGTHDM->computeSigmaggH_tt(8.0);
    //The square of the top-quark contribution to the ggH cross section in the SM at 13 TeV
    //double Sigmaggh_tt13 = myTHDM->computeSigmaggH_tt(13.0);
    //The square of the bottom-quark contribution to the ggH cross section in the SM at 8 TeV
    double Sigmaggh_bb8 = myGTHDM->computeSigmaggH_bb(8.0);
    //The square of the bottom-quark contribution to the ggH cross section in the SM at 13 TeV
    //double Sigmaggh_bb13 = myTHDM->computeSigmaggH_bb(13.0);
    //The ttH production cross section in the SM at 8 TeV
    double Sigmatth8 = myGTHDM->computeSigmattH(8.0);
    //The Vh cross section in the SM at 8 TeV.
    SigmaVh8 = myGTHDM->computeSigmaWH(8.0)+myGTHDM->computeSigmaZH(8.0);
    //The VBF cross section in the SM at 8 TeV.
    SigmaVBF8 = myGTHDM->computeSigmaVBF(8.0);
    //The bbH production cross section in the SM at 8 TeV
    double Sigmabbh8 = ip_cs_pptobbH_8(mHl);

    //The ggH cross section in the SM at 13 TeV.
    SigmaggF13 = myGTHDM->computeSigmaggH(13.0);
    //The Vh cross section in the SM at 13 TeV.
    SigmaVh13 = myGTHDM->computeSigmaWH(13.0)+myGTHDM->computeSigmaZH(13.0);
    //The Zh cross section in the SM at 13 TeV.
    SigmaZh13 = myGTHDM->computeSigmaZH(13.0);
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
       the particle i and the respective squared SM coupling. Where separated
       E means the even part (coming from the CP-even scalar) and O
       the odd part (coming from the CP-odd scalar) */
                 
    rh_QuQuE= yu1.real()*yu1.real(); 
    rh_QuQuO= yu1.imag()*yu1.imag(); 
    rh_QdQdE= yd1.real()*yd1.real(); 
    rh_QdQdO= yd1.imag()*yd1.imag(); 
    rh_QlQlE= yl1.real()*yl1.real(); 
    rh_QlQlO= yl1.imag()*yl1.imag(); 
    //Isn't the interference missing here???
    rh_ggE = yu1.real()*yd1.real() + (yu1.real()*yu1.real() - yu1.real()*yd1.real())*(Sigmaggh_tt8/SigmaggF8)  + (yd1.real()*yd1.real() - yu1.real()*yd1.real())*(Sigmaggh_bb8/SigmaggF8);
    rh_ggO = yu1.imag()*yu1.imag() + (yu1.imag()*yu1.imag() - yu1.imag()*yd1.imag())*rSigmagghO_b8  + (yd1.imag()*yd1.imag() - yu1.imag()*yd1.imag())*rSigmagghO_b8;
    rh_gg = rh_ggE+rh_ggO*(SigmagghO_8/SigmaggF8);
    rh_VV=Rij_GTHDM(0,0)*Rij_GTHDM(0,0);        



    /*Loop functions needed to rh_gaga and rh_Zga ...*/
    
    
    gslpp::complex fermU=I_h_U(m1_2,Mu,Mc,Mt);
    gslpp::complex fermD=I_h_D(m1_2,Md,Ms,Mb);
    gslpp::complex fermL=I_h_L(m1_2,Me,Mmu,Mtau);
    gslpp::complex I_hSM_W=I_H_W(mHl,MW);
    gslpp::complex I_h_W=Rij_GTHDM(0,0)*I_hSM_W;
     
    
    gslpp::complex I_hSM_F= fermU+fermD+fermL;
    gslpp::complex I_hE_F= yu1.real()*fermU+ yd1.real()*fermD+yl1.real()*fermL;
    
                                                                               
    /*Coupling between h and two charged Higgs*/
    
       
    double lambdahHpHm = lambdaipm(Rij_GTHDM(0,0), Rij_GTHDM(0,1), Rij_GTHDM(0,2));

    gslpp::complex I_h_Hp=(vev)/(2.0*mHp2)*I_H_Hp(mHp2,m1)*lambdahHpHm;
       
    /*CP ODD */

    gslpp::complex I_h_Ux=I_A_U(m1_2,Mc,Mt);
    gslpp::complex I_h_Dx=I_A_D(m1_2,Ms,Mb);
    gslpp::complex I_h_Lx=I_A_L(m1_2,Mmu,Mtau);

    gslpp::complex I_hO_F = yu1.imag()*I_h_Ux + yd1.imag()*I_h_Dx + yl1.imag()*I_h_Lx;

   // double Gamma_hgaga=(GF*Ale*Ale*m1*m1*m1/(sqrt(2.0)*128.0*M_PI*M_PI*M_PI))*((I_hE_F+I_h_W+I_h_Hp).abs2()+ (I_hO_F).abs2());
    rh_gaga = ((I_hE_F+I_h_W+I_h_Hp).abs2()+ (I_hO_F).abs2())/(I_hSM_F +I_hSM_W).abs2();    
    
   /* std::cout << "rh_gaga = " << rh_gaga << std::endl;
    std::cout << "I_hE_F = " << I_hE_F << std::endl;
    std::cout << "I_h_W = " << I_h_W << std::endl;
    std::cout << "I_h_Hp = " << I_h_Hp << std::endl;
    std::cout << "I_hO_F = " << I_hO_F << std::endl;
    std::cout << "I_hSM_F = " << I_hSM_F << std::endl;
    std::cout << "I_hSM_W = " << I_hSM_W << std::endl;*/

   /* std::cout << "yu1 c = " << yu1 << std::endl;
    std::cout << "yd1 c = " << yd1 << std::endl;
    std::cout << "yl1 c = " << yl1 << std::endl;*/

    
    /*Decay to Z gamma
    CP-EVEN PART*/

    gslpp::complex A_hE_Ux = A_h_U(m1_2,cW2,Mu,Mc,Mt,MZ);
    gslpp::complex A_hE_Dx = A_h_D(m1_2,cW2,Md,Ms,Mb,MZ);
    gslpp::complex A_hE_Lx = A_h_L(m1_2,cW2,Me,Mmu,Mtau,MZ);
    gslpp::complex A_hSM_W = A_H_W(m1,cW2,MW,MZ);
    gslpp::complex A_h_W = Rij_GTHDM(0,0)*A_hSM_W;
    
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
    //std::cout<<"\033[1;32m   rh_Zga= \033[0m "<< rh_Zga <<std::endl;
    /*Decay to gluons*/
    //std::cout<<"\033[1;32m   Ale= \033[0m "<< Ale <<std::endl;
    
    //std::cout<<"\033[1;32m   sW2= \033[0m "<< sW2 <<std::endl;
    
//    double Gamma_hgg=rh_gg*GF*Als*Als*m1*m1*m1/(sqrt(2.0)*16.0*M_PI*M_PI*M_PI)*(9.0/4.0)*(fermU/4.0+fermD).abs2();

//    double lambda122 = (2.0)*(lambdaijk(R11, R12, R13, R21, R22, R23, R21, R22, R23, lambda1, lambda3, lambda4, Relambda5,
//                                        Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7) +
//                              lambdaijk(R21, R22, R23, R11, R12, R13, R21, R22, R23,   lambda1, lambda3, lambda4, Relambda5,
//                                        Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7) +
//                              lambdaijk(R21, R22, R23, R21, R22, R23, R11, R21, R13,   lambda1, lambda3, lambda4, Relambda5,
//                                        Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7 ));

    double lambda345 = (lambda3 + lambda4 + Relambda5);

    //Modify for CPV cases
    double lambda122 = vev*(3.0*R21*R21*(lambda1*R11 + Relambda6*R12) + 2.0*R21*R22*(3.0*Relambda6*R11 + lambda345*R12) + R22*R22*(lambda345*R11 + 3.0*Relambda7*R12));

    double lambda133 = (vev)*(Relambda7*R12 + (lambda3 + lambda4 - Relambda5)*R11);

    double lambda1pm = vev*(lambda3*R11 + Relambda7*R12 - Imlambda7*R13);

    double Gamma_hHH = HSTheta(m1 - 2.0*m2)*sqrt(std::fabs(1 - 4*(m2*m2)/(m1*m1)))*(lambda122*lambda122)/(32.0*M_PI*m1);

    // KaellenFunction as defined above corresponds to the usual sqrt(lambda) divided by 2*m1
    //double Gamma_hHH = HSTheta(m1 - 2.0*m2)*(KaellenFunction(m1_2,m2_2,m2_2)*lambda122*lambda122)/(8.0*m1_2*M_PI);

    // KaellenFunction as defined above corresponds to the usual sqrt(lambda) divided by 2*m1
    //double Gamma_hAA = HSTheta(m1 - 2.0*m3)*(KaellenFunction(m1_2,m3_2,m3_2)*lambda133*lambda133)/(8.0*m1_2*M_PI);

    double Gamma_hAA = HSTheta(m1 - 2.0*m3)*sqrt(std::fabs(1 - 4*(m3*m3)/(m1*m1)))*(lambda133*lambda133)/(32.0*M_PI*m1);

    double Gamma_hHpHm  = HSTheta(m1-2.0*mHp)*sqrt(std::fabs(1 - 4*(mHp*mHp)/(m1*m1)))*(lambda1pm*lambda1pm)/(16.0*M_PI*m1);

    //Check the sign of Rij in CPV case
    // KaellenFunction as defined above corresponds to the usual sqrt(lambda) divided by 2*m1
    double Gamma_hHZ  = HSTheta(m1-(m2+MZ))*pow(KaellenFunction(m1_2,MZ*MZ,m2_2),3)*(R23*R12 + R22*R13)*(R23*R12 + R22*R13)/(2.0*M_PI*vev*vev);//h -> phi2 Z
    double Gamma_hAZ  = HSTheta(m1-(m3+MZ))*pow(KaellenFunction(m1_2,MZ*MZ,m3_2),3)*(R33*R12 + R32*R13)*(R33*R12 + R32*R13)/(2.0*M_PI*vev*vev);//h -> phi3 Z    
    
  //  /* ggF_tth8 is the ratio of the THDM and SM cross sections for ggF or tth production at 8 TeV*/
  //  ggF_tth8 = (SigmaggF8*rh_ggE + SigmagghO_8*rh_ggO + Sigmatth8*(rh_QuQuE + rh_QuQuO/(beta(Mc, m1_2)*beta(Mc, m1_2))))/(SigmaggF8 + Sigmatth8);
  //  /* ggF_tth13 is the ratio of the THDM and SM cross sections for ggF or tth production at 13 TeV */
  //  ggF_tth13 = (SigmaggF13*rh_ggE + SigmagghO_13*rh_ggO + Sigmatth8*(rh_QuQuE + rh_QuQuO/(beta(Mc, m1_2)*beta(Mc, m1_2))))/(SigmaggF13 + Sigmatth13);
  //  /* pph13 is the ratio of the THDM and SM cross sections for an h production at 13 TeV */
  //  pph13 = (SigmaggF13*rh_ggE + SigmagghO_13*rh_ggO+ SigmaVBFVh13*rh_VV + Sigmatth13*(rh_QuQuE + rh_QuQuO/(beta(Mc, m1_2)*beta(Mc, m1_2))) + Sigmabbh13*(rh_QdQdE + rh_QdQdE/(beta(Mb, m1_2)*beta(Mb, m1_2))))/(SigmaggF13 + SigmaVBFVh13 + Sigmatth13 + Sigmabbh13);
  //  /* VBF_Vh is the ratio of the THDM and SM cross sections for VBF or Vh production */
  //  VBF_Vh = rh_VV;

    /* ggF_tth8 is the ratio of the THDM and SM cross sections for ggF and tth production at 8 TeV*/
    //ggF_tth8 = (SigmaggF8*rh_gg + Sigmatth8*(rh_QuQuE + rh_QuQuO/(beta(Mc, m1_2)*beta(Mc, m1_2))))/(SigmaggF8 + Sigmatth8);
    ggF_tth8 = (SigmaggF8*rh_gg + Sigmatth8*(rh_QuQuE + rh_QuQuO))/(SigmaggF8 + Sigmatth8);
    
    /*SM cross-section of higgs production at 8 TeV at LHC*/
    SigSM_pph8 = (SigmaggF8 + SigmaVBF8 + SigmaVh8 + Sigmatth8 + Sigmabbh8);
    /* pph8 is the ratio of the THDM and SM cross sections for h production at 8 TeV */
    //pph8 = (SigmaggF8*rh_gg+ (SigmaVBF8 + SigmaVh8)*rh_VV + Sigmatth8*(rh_QuQuE + rh_QuQuO/(beta(Mc, m1_2)*beta(Mc, m1_2))) + Sigmabbh8*(rh_QdQdE + rh_QdQdE/(beta(Mb, m1_2)*beta(Mb, m1_2))))/SigSM_pph8;
    pph8 = (SigmaggF8*rh_gg+ (SigmaVBF8 + SigmaVh8)*rh_VV + Sigmatth8*(rh_QuQuE + rh_QuQuO) + Sigmabbh8*(rh_QdQdE + rh_QdQdE))/SigSM_pph8;
    
    /* ggF_tth13 is the ratio of the THDM and SM cross sections for ggF and tth production at 13 TeV */
    //ggF_tth13 = (SigmaggF13*rh_gg + Sigmatth8*(rh_QuQuE + rh_QuQuO/(beta(Mc, m1_2)*beta(Mc, m1_2))))/(SigmaggF13 + Sigmatth13);
    ggF_tth13 = (SigmaggF13*rh_gg + Sigmatth8*(rh_QuQuE + rh_QuQuO))/(SigmaggF13 + Sigmatth13);
    
    /*SM cross-section of higgs production at 13 TeV at LHC*/
    SigSM_pph13 = (SigmaggF13 + SigmaVBFVh13 + Sigmatth13 + Sigmabbh13);
    /* pph13 is the ratio of the THDM and SM cross sections for h production at 13 TeV */
    //pph13 = (SigmaggF13*rh_gg+ SigmaVBFVh13*rh_VV + Sigmatth13*(rh_QuQuE + rh_QuQuO/(beta(Mc, m1_2)*beta(Mc, m1_2))) + Sigmabbh13*(rh_QdQdE + rh_QdQdE/(beta(Mb, m1_2)*beta(Mb, m1_2))))/SigSM_pph13;
    pph13 = (SigmaggF13*rh_gg+ SigmaVBFVh13*rh_VV + Sigmatth13*(rh_QuQuE + rh_QuQuO) + Sigmabbh13*(rh_QdQdE + rh_QdQdE))/SigSM_pph13;
    
    /* VBF_Vh is the ratio of the THDM and SM cross sections for VBF and Vh production */
    VBF_Vh = rh_VV;

    /* ggF_VBF_Vh13 is the ratio of the GTHDM and SM cross sections for ggF + VBF + Vh production at 13 TeV */
    ggF_VBF_Vh13 = (SigmaggF13 * rh_gg + SigmaVBFVh13 * rh_VV)/(SigmaggF13 + SigmaVBFVh13);

    double Gamma_h_exp = myGTHDM->computeGammaHTotal();

    //double Gamma_htobb     = BrSM_htobb * (rh_QdQdE + rh_QdQdO/(beta(Mb, m1_2)*beta(Mb, m1_2))) * Gamma_h_exp;
    double Gamma_htobb     = BrSM_htobb * (rh_QdQdE + rh_QdQdO) * Gamma_h_exp;
    double Gamma_htoWW     = BrSM_htoWW * rh_VV * Gamma_h_exp;
    double Gamma_htoZZ     = BrSM_htoZZ * rh_VV * Gamma_h_exp;
    //double Gamma_htotautau = BrSM_htotautau * (rh_QlQlE + rh_QlQlO/(beta(Mtau, m1_2)*beta(Mtau, m1_2))) * Gamma_h_exp;
    double Gamma_htotautau = BrSM_htotautau * (rh_QlQlE + rh_QlQlO) * Gamma_h_exp;
    double Gamma_htogaga   = BrSM_htogaga * rh_gaga * Gamma_h_exp;
    double Gamma_htogg     = BrSM_htogg * rh_gg * Gamma_h_exp;
    double Gamma_htoZga    = BrSM_htoZga * rh_Zga * Gamma_h_exp;
    //double Gamma_htocc     = BrSM_htocc * (rh_QuQuE + rh_QuQuO/(beta(Mc, m1_2)*beta(Mc, m1_2))) * Gamma_h_exp;
    double Gamma_htocc     = BrSM_htocc * (rh_QuQuE + rh_QuQuO) * Gamma_h_exp;

    Gamma_h_inv = Gamma_hHpHm + Gamma_hHH + Gamma_hAA + Gamma_hHZ + Gamma_hAZ;

    Gamma_h = (Gamma_htobb + Gamma_htoWW + Gamma_htoZZ + Gamma_htotautau + Gamma_htogaga +
               Gamma_htogg + Gamma_htoZga + Gamma_htocc + Gamma_h_inv);

    GTHDM_BR_h_bb = Gamma_htobb / Gamma_h;
    GTHDM_BR_h_WW = Gamma_htoWW / Gamma_h;
    GTHDM_BR_h_ZZ = Gamma_htoZZ / Gamma_h;
    GTHDM_BR_h_tautau = Gamma_htotautau / Gamma_h;
    GTHDM_BR_h_cc = Gamma_htocc / Gamma_h;
    GTHDM_BR_h_gaga = Gamma_htogaga / Gamma_h;
    GTHDM_BR_h_gg = Gamma_htogg / Gamma_h;
    GTHDM_BR_h_HpHm = Gamma_hHpHm / Gamma_h;
    GTHDM_BR_h_HH = Gamma_hHH / Gamma_h;
    GTHDM_BR_h_AA = Gamma_hAA / Gamma_h;
    GTHDM_BR_h_HZ = Gamma_hHZ / Gamma_h;
    GTHDM_BR_h_AZ = Gamma_hAZ / Gamma_h;

    sumModBRs = Gamma_h / Gamma_h_exp;
}

    

void GeneralTHDMcache::computephi2quantities()
{
    //FLAG to select only the model in which all the couplings are the same (by families)
    if (!myGTHDM->getATHDMflag())
    {
        throw std::runtime_error("Direct Searches are only aviable in the A2HDM.");
    }

    /*complex i */
     gslpp::complex i = gslpp::complex::i();

    //Here we define the couplings to the fermions for the mass states
    yu2 = Rij_GTHDM(1,0) + (Rij_GTHDM(1,1) - i*Rij_GTHDM(1,2))*su.conjugate();
    yd2 = Rij_GTHDM(1,0) + (Rij_GTHDM(1,1) + i*Rij_GTHDM(1,2))*sd;
    yl2 = Rij_GTHDM(1,0) + (Rij_GTHDM(1,1) + i*Rij_GTHDM(1,2))*sl;

    yu2R = Rij_GTHDM(1,0) + (Rij_GTHDM(1,1))*su.real();
    yd2R = Rij_GTHDM(1,0) + (Rij_GTHDM(1,1))*sd.real();
    yl2R = Rij_GTHDM(1,0) + (Rij_GTHDM(1,1))*sl.real();

    //These cross sections ratios are necessary for rphi2_gg
    //At 8 TeV
    
    //SM gg -> H (phi even) production cross section ratio at 8 TeV, charm and top loop only over total
    double rSigmaggphi2E_tc8 = ip_csr_ggH_tc_8(m2);
    //SM gg -> H (phi even) production cross section ratio at 8 TeV, bottom loop only over total
    double rSigmaggphi2E_b8 = ip_csr_ggH_b_8(m2);
    //gg -> H (phieven) production cross section at 8 TeV, total
   //  double Sigmaggphi2E_8 = ip_cs_ggtoH_8(m2);

    
    //gg -> A (phi odd) production cross section ratio at 8 TeV, charm and top loop only over total
    double rSigmaggphi2O_tc8 = ip_csr_ggA_tc_8(m2);
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
    rphi2_ggE = yu2.real()*yd2.real() + (yu2.real()*yu2.real() - yu2.real()*yd2.real())*rSigmaggphi2E_tc8  + (yd2.real()*yd2.real() - yu2.real()*yd2.real())*rSigmaggphi2E_b8;
    rphi2_ggO = yu2.imag()*yu2.imag() + (yu2.imag()*yu2.imag() - yu2.imag()*yd2.imag())*rSigmaggphi2O_tc8  + (yd2.imag()*yd2.imag() - yu2.imag()*yd2.imag())*rSigmaggphi2O_b8;

    rphi2_VV=Rij_GTHDM(1,0)*Rij_GTHDM(1,0);
    

    /*Gamma_phi2gaga and Gamma_phi2Zga expressions ...*/
    
    /*Decay to photons. The fermionic contribution has a CP-even part (HH) and a CP-odd (A)*/
    /*CP EVEN*/
   
                                                                                  
    gslpp::complex I_HH2_Ux=I_HH_U(m2_2,Mc,Mt);
    gslpp::complex I_HH2_Dx=I_HH_D(m2_2,Ms,Mb);
    gslpp::complex I_HH2_Lx=I_HH_L(m2_2,Mmu,Mtau);
    gslpp::complex I_phi2E_F= yu2.real()*I_HH2_Ux+ yd2.real()*I_HH2_Dx+yl2.real()*I_HH2_Lx;
                                                                               
    gslpp::complex I_phi2_W=Rij_GTHDM(1,0)*I_H_W(m2,MW);
    
    
    double lambdaphi2HpHm = lambdaipm(R21, R22, R32);
    
    gslpp::complex I_phi2_Hp=(vev*vev)/(2.0*mHp2)*I_H_Hp(mHp2,m2)*(lambdaphi2HpHm);
    
    
    /*CP ODD */
            
    gslpp::complex I_A2_Ux=I_A_U(m2_2,Mc,Mt);
    gslpp::complex I_A2_Dx=I_A_D(m2_2,Ms,Mb);
    gslpp::complex I_A2_Lx=I_A_L(m2_2,Mmu,Mtau);
                                                                               
    gslpp::complex I_phi2O_F = yu2.imag()*I_A2_Ux + yd2.imag()*I_A2_Dx + yl2.imag()*I_A2_Lx;
                                                                             
    double Gamma_phi2gaga=(GF*Ale*Ale*m2*m2*m2/(sqrt(2.0)*128.0*M_PI*M_PI*M_PI))*((I_phi2E_F+I_phi2_W+I_phi2_Hp).abs2()
                        + (I_phi2O_F).abs2());
    
    //std::cout<<"\033[1;33m  Ale= \033[0m "<< Ale <<std::endl;
    
    /*Decay to Z gamma
    CP-EVEN PART*/

    gslpp::complex A_HH2_Ux = A_HH_U(m2_2,cW2,Mc,Mt,MZ);
    gslpp::complex A_HH2_Dx = A_HH_D(m2_2,cW2,Ms,Mb,MZ);
    gslpp::complex A_HH2_Lx = A_HH_L(m2_2,cW2,Mmu,Mtau,MZ);
                                                                               
    gslpp::complex A_phi2E_F = (yu2.real()*A_HH2_Ux+ yd2.real()*A_HH2_Dx+ yl2.real()*A_HH2_Lx)/sqrt(sW2*cW2);
    
    
    
    gslpp::complex A_phi2_W = Rij_GTHDM(1,0)*A_H_W(m2,cW2,MW,MZ);
    
    gslpp::complex A_phi2_Hp = (vev*vev)/(2.0*mHp2)*A_H_Hp(mHp2,m2,cW2,MZ)*(lambdaphi2HpHm);

    /*CP-ODD PART*/
                                                                               
    gslpp::complex A_A2_Ux = A_A_U(m2_2,cW2,Mc,Mt,MZ);
    gslpp::complex A_A2_Dx = A_A_D(m2_2,cW2,Ms,Mb,MZ);
    gslpp::complex A_A2_Lx = A_A_L(m2_2,cW2,Mmu,Mtau,MZ);
                                                                               
    gslpp::complex A_phi2O_F=yu2.imag()*A_A2_Ux + yd2.imag()*A_A2_Dx + yl2.imag()*A_A2_Lx;
    
                                                                             
                                                                               
    double Gamma_phi2Zga=HSTheta(m2-MZ)*GF*Ale*Ale*m2*m2*m2/(sqrt(2.0)*64.0*M_PI*M_PI*M_PI)*(1.0-MZ*MZ/(m2*m2))*(1.0-MZ*MZ/(m2*m2))*(1.0-MZ*MZ/(m2*m2))*((A_phi2E_F+A_phi2_W+A_phi2_Hp).abs2()
                        + A_phi2O_F.abs2());
       
    
    
    //std::cout<<"\033[1;35m   Gamma_phi2Zga= \033[0m "<<  Gamma_phi2Zga  <<std::endl;
    
    
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


    SigmaggF_phi2_8=ip_cs_ggtoH_8(m2)*rphi2_ggE + ip_cs_ggtoA_8(m2)* rphi2_ggO;
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

//    std::cout<< SigmaggF_phi2_13<<std::endl;
//    std::cout<<SigmaVBF_phi2_13;
    
    /*std::cout << "THoEX_bb_phi2_tautau_ATLAS13 = " << THoEX_bb_phi2_tautau_ATLAS13 << std::endl;
    std::cout << "bb_phi2_tautau_TH13 = " << bb_phi2_tautau_TH13 << std::endl;
    std::cout << "ip_ex_bb_phi_tautau_ATLAS13(m2) = " << ip_ex_bb_phi_tautau_ATLAS13(m2) << std::endl;
    std::cout << "SigmabbF_phi2_13 = " << SigmabbF_phi2_13 << std::endl;
    std::cout << "ip_cs_pptobbH_13(m2) = " << ip_cs_pptobbH_13(m2) << std::endl;
    std::cout << "rphi2_QdQdE = " << rphi2_QdQdE << std::endl;
    std::cout << "yd2 = " << yd2 << std::endl;
    std::cout << "R21 = " << R21 << std::endl;
    std::cout << "R22 = " << R22 << std::endl;
    std::cout << "R23 = " << R23 << std::endl;
    std::cout << "sd " << sd << std::endl;
    std::cout << "ip_cs_pptobbA_13(m2) = " << ip_cs_pptobbA_13(m2) << std::endl;
    std::cout << "rphi2_QdQdO = " << rphi2_QdQdO << std::endl;*/
    
    
 
//    double SigmaTotSM_H13 = 1.0e-15;
//    if (mHh>=20. && mHh <=2000.) {
//            SigmaTotSM_H13=ip_cs_ggtoH_13(mHh)+ip_cs_VBFtoH_13(mHh)+ip_cs_WtoWH_13(mHh)+ip_cs_ZtoZH_13(mHh)+ip_cs_pptottH_13(mHh)+ip_cs_pptobbH_13(mHh);
//    }
    SigmaSumphi2_13 = SigmaggF_phi2_13 + SigmaVBF_phi2_13 + SigmaVH_phi2_13 + SigmattF_phi2_13 + SigmabbF_phi2_13;
     
    double BrSM_phi2tott=ip_Br_HPtott(m2);
    double BrSM_phi2tocc=ip_Br_HPtocc(m2);
    double BrSM_phi2tobb=ip_Br_HPtobb(m2);
    double BrSM_phi2totautau=ip_Br_HPtotautau(m2);
    double BrSM_phi2tomumu=ip_Br_HPtomumu(m2);
    double BrSM_phi2toWW =ip_Br_HPtoWW(m2);
    double BrSM_phi2toZZ =ip_Br_HPtoZZ(m2);

Gammaphi2totSM=ip_GammaHPtotSM(m2);
    
 /*Decay of phi3 to the others scalars*/
double lambda132 = lambdaijk(R11, R12, R13, R31, R32, R33, R21, R22, R23, lambda1, lambda3, lambda4, Relambda5,
                             Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7) +
                   lambdaijk(R21, R22, R23, R11, R12, R13, R31, R32, R33, lambda1, lambda3, lambda4, Relambda5,
                             Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7) +
                   lambdaijk(R21, R22, R23, R31, R32, R33, R11, R12, R13, lambda1, lambda3, lambda4, Relambda5,
                             Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7) +
                   lambdaijk(R11, R12, R13, R21, R22, R23, R32, R32, R32, lambda1, lambda3, lambda4, Relambda5,
                             Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7) +
                   lambdaijk(R31, R32, R33, R11, R12, R13, R21, R22, R23, lambda1, lambda3, lambda4, Relambda5,
                             Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7) +
                   lambdaijk(R31, R32, R33, R21, R22, R23, R11, R12, R13, lambda1, lambda3, lambda4, Relambda5,
                             Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7);

double lambda332 = lambdaijk(R31, R32, R33, R31, R32, R33, R21, R22, R23, lambda1, lambda3, lambda4, Relambda5,
                             Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7) +
                   lambdaijk(R21, R22, R23, R31, R32, R33, R31, R32, R33, lambda1, lambda3, lambda4, Relambda5,
                             Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7) +
                   lambdaijk(R31, R32, R33, R21, R22, R23, R31, R32, R33, lambda1, lambda3, lambda4, Relambda5,
                             Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7);

//phi2 -> phi1 phi3
double Gammaphi2_phi1phi3=HSTheta(m2 - (m1+m3))*KaellenFunction(m2_2,m1_2,m3_2)*lambda132*lambda132/(8.0*m3_2*M_PI);
double Gammaphi2_phi3phi3=HSTheta(m2 - 2.0*m3)*sqrt(std::fabs(1.0 - (4.0*m3_2)/m2_2))*lambda332*lambda332/(32.0*m2*M_PI);
    
/*Decay of phi2 to the others scalars*/
double lambda112 = (2.0)*(lambdaijk(R11, R12, R13, R21, R22, R23, R11, R12, R13, lambda1, lambda3, lambda4, Relambda5,
                                    Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7) +
                          lambdaijk(R11, R12, R13, R11, R12, R13, R21, R22, R23, lambda1, lambda3, lambda4, Relambda5,
                                    Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7) +
                          lambdaijk(R21, R22, R23, R11, R12, R13, R11, R12, R13, lambda1, lambda3, lambda4, Relambda5,
                                    Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7));

//phi2 -> phi1phi1
double Gammaphi2_phi1phi1=HSTheta(m2 - 2.0*m1)*sqrt(std::fabs(1.0 - (4.0*m1_2)/m2_2))*lambda112*lambda112/(32.0*m2*M_PI);

//std::cout<<"\033[1;35m   m2 =\033[0m "<< m2 <<std::endl;
//std::cout<<"\033[1;35m   m1 =\033[0m "<< m1 <<std::endl;
//std::cout<<"\033[1;35m   m1_2 =\033[0m "<< m1_2 <<std::endl;
//std::cout<<"\033[1;35m   m2_2 =\033[0m "<< m2_2 <<std::endl;
//std::cout<<"\033[1;35m   sqrt(std::fabs(1.0 - (4.0*m1_2)/m2_2)) =\033[0m "<< sqrt(std::fabs(1.0 - (4.0*m1_2)/m2_2)) <<std::endl;


//std::cout<<"\033[1;35m   lambda112*lambda112 =\033[0m "<< lambda112*lambda112 <<std::endl;
//std::cout<<"\033[1;35m   1/(32.0*m2*M_PI) =\033[0m "<< 1/(32.0*m2*M_PI) <<std::endl;


//std::cout<<"\033[1;35m   HSTheta(m2 - 2.0*m1) =\033[0m "<< HSTheta(m2 - 2.0*m1) <<std::endl;


 //phi2 ->H+H-
double Gammaphi2_HpHm=HSTheta(m2 - 2.0*sqrt(mHp2))*sqrt(std::fabs(1.0 - (4.0*mHp2)/m2_2))*lambdaphi2HpHm*lambdaphi2HpHm/(8.0*m2*M_PI);
//phi2 -> phi1 Z
double Gammaphi2_phi1Z=HSTheta(m2-(m1+MZ))*pow(KaellenFunction(m2_2,MZ*MZ,m1_2),3)*(R23*R12 + R22*R13)*(R23*R12 + R22*R13)/(2.0*M_PI*vev*vev);

//phi2 -> phi3 Z
double Gammaphi2_phi3Z=HSTheta(m2-(m3+MZ))*pow(KaellenFunction(m2_2,MZ*MZ,m3_2),3)*(R33*R22 + R32*R23)*(R33*R22 + R32*R23)/(2.0*M_PI*vev*vev);


/* phi2 -> H+W- */
double Gammaphi2_HpW=HSTheta(m2-sqrt(mHp2)-MW)*pow(KaellenFunction(m2_2,MW*MW,mHp2),3)*(R23-i*R22).abs2()/(M_PI*vev*vev);

double Gammaphi2_tt = HSTheta(m2 - 2*Mt)*BrSM_phi2tott*(rphi2_QuQuE + rphi2_QuQuO/(beta(Mt, m2_2)*beta(Mt, m2_2)))*Gammaphi2totSM;

double Gammaphi2_cc = HSTheta(m2 - 2*Mc)*BrSM_phi2tocc*(rphi2_QuQuE + rphi2_QuQuO/(beta(Mc, m2_2)*beta(Mc, m2_2)))*Gammaphi2totSM;

double Gammaphi2_bb = HSTheta(m2 - 2*Mb)*BrSM_phi2tobb*(rphi2_QdQdE + rphi2_QdQdO/(beta(Mb, m2_2)*beta(Mb, m2_2)))*Gammaphi2totSM;

double Gammaphi2_mumu = HSTheta(m2 - 2*Mmu)*BrSM_phi2tomumu*(rphi2_QlQlE + rphi2_QlQlO/(beta(Mmu, m2_2)*beta(Mmu, m2_2)))*Gammaphi2totSM;

double Gammaphi2_tautau = HSTheta(m2 - 2*Mtau)*BrSM_phi2totautau*(rphi2_QlQlE + rphi2_QlQlO/(beta(Mtau, m2_2)*beta(Mtau, m2_2)))*Gammaphi2totSM;

double Gammaphi2_WW = HSTheta(m2 - 2*MW)*BrSM_phi2toWW*(rphi2_VV)*Gammaphi2totSM;

double Gammaphi2_ZZ = HSTheta(m2 - 2*MZ)*BrSM_phi2toZZ*(rphi2_VV)*Gammaphi2totSM;




Gammaphi2tot = 1.e-10;

Gammaphi2tot= Gammaphi2tot + Gammaphi2_tt+Gammaphi2_cc
            +Gammaphi2_bb+Gammaphi2_tautau+Gammaphi2_mumu
            +Gammaphi2_WW+Gammaphi2_ZZ+Gamma_phi2gaga
            +Gamma_phi2Zga+Gamma_phi2gg +Gammaphi2_phi1phi3
            +Gammaphi2_phi1phi1+Gammaphi2_phi3phi3+Gammaphi2_HpHm
            +Gammaphi2_phi1Z+Gammaphi2_phi3Z+Gammaphi2_HpW;

    
    //std::cout<<"\033[1;34m   R11 =\033[0m "<< R11 <<std::endl;
    
    //std::cout<<"\033[1;34m   R23 =\033[0m "<< R23 <<std::endl;
    //std::cout<<"\033[1;34m   R22 =\033[0m "<< R22 <<std::endl;

    //std::cout<<"\033[1;34m   lambdaphi2HpHm =\033[0m "<< lambdaphi2HpHm <<std::endl;
    //std::cout<<"\033[1;34m   lambdaphi2HpHm =\033[0m "<< lambdaphi2HpHm <<std::endl;
    //std::cout<<"\033[1;34m   Gammaphi2_HpHm =\033[0m "<< Gammaphi2_HpHm <<std::endl;
    //std::cout<<"\033[1;34m   Gammaphi2_HpW =\033[0m "<< Gammaphi2_HpW <<std::endl;
  
    


    Br_phi2tott=Gammaphi2_tt/Gammaphi2tot;
    Br_phi2tobb=Gammaphi2_bb/Gammaphi2tot;
    Br_phi2totautau=Gammaphi2_tautau/Gammaphi2tot;
    Br_phi2tomumu=Gammaphi2_mumu/Gammaphi2tot;
    
    Br_phi2toWW=Gammaphi2_WW/Gammaphi2tot;
    Br_phi2toZZ=Gammaphi2_ZZ/Gammaphi2tot;
    Br_phi2togaga=Gamma_phi2gaga/Gammaphi2tot;
    Br_phi2togg=Gamma_phi2gg/Gammaphi2tot;
    Br_phi2toZga=Gamma_phi2Zga/Gammaphi2tot;
    Br_phi2tophi1phi1=Gammaphi2_phi1phi1/Gammaphi2tot;
    Br_phi2tophi3phi3=Gammaphi2_phi3phi3/Gammaphi2tot;
    Br_phi2tophi1phi3=Gammaphi2_phi1phi3/Gammaphi2tot;
    Br_phi2toHpHm=Gammaphi2_HpHm/Gammaphi2tot;
    Br_phi2tophi1Z=Gammaphi2_phi1Z/Gammaphi2tot;
    Br_phi2tophi3Z=Gammaphi2_phi3Z/Gammaphi2tot;
    Br_phi2toHpW=Gammaphi2_HpW/Gammaphi2tot;
    
    
    //std::cout<<"\033[1;30m   Gammaphi2_tt =\033[0m "<< Gammaphi2_tt <<std::endl;
    //std::cout<<"\033[1;30m   Gammaphi2_bb =\033[0m "<< Gammaphi2_bb <<std::endl;
    //std::cout<<"\033[1;30m   Gammaphi2_cc =\033[0m "<< Gammaphi2_cc <<std::endl;
    //std::cout<<"\033[1;30m   Gammaphi2_tautau =\033[0m "<< Gammaphi2_tautau <<std::endl;
    //std::cout<<"\033[1;30m   Gammaphi2_mumu =\033[0m "<< Gammaphi2_mumu <<std::endl;
    
    //std::cout<<"\033[1;30m   Gammaphi2_WW =\033[0m "<< Gammaphi2_WW <<std::endl;
    //std::cout<<"\033[1;30m   Gammaphi2_ZZ =\033[0m "<< Gammaphi2_ZZ <<std::endl;
    //std::cout<<"\033[1;30m   Gamma_phi2gaga =\033[0m "<< Gamma_phi2gaga <<std::endl;
    //std::cout<<"\033[1;30m   Gamma_phi2Zga =\033[0m "<< Gamma_phi2Zga <<std::endl;
    //std::cout<<"\033[1;30m   Gamma_phi2gg =\033[0m "<< Gamma_phi2gg <<std::endl;
    
    //std::cout<<"\033[1;30m   Gammaphi2_phi1phi3 =\033[0m "<< Gammaphi2_phi1phi3 <<std::endl;
    //std::cout<<"\033[1;30m   Gammaphi2_phi1phi1 =\033[0m "<< Gammaphi2_phi1phi1 <<std::endl;
    //std::cout<<"\033[1;30m   Gammaphi2_phi3phi3 =\033[0m "<< Gammaphi2_phi3phi3 <<std::endl;
    //std::cout<<"\033[1;30m   Gammaphi2_HpHm =\033[0m "<< Gammaphi2_HpHm <<std::endl;
    //std::cout<<"\033[1;30m   Gammaphi2_phi1Z =\033[0m "<< Gammaphi2_phi1Z <<std::endl;
    
    //std::cout<<"\033[1;30m   Gammaphi2_phi3Z =\033[0m "<< Gammaphi2_phi3Z <<std::endl;
    //std::cout<<"\033[1;30m   Gammaphi2_HpW =\033[0m "<< Gammaphi2_HpW <<std::endl;
}


    
void GeneralTHDMcache::computephi3quantities()
{
    //FLAG to select only the model in which all the couplings are the same (by families)
    if (!myGTHDM->getATHDMflag())
    {
        throw std::runtime_error("Direct Searches are only aviable in the A2HDM.");
    }

    /*complex i */
    gslpp::complex i = gslpp::complex::i();

    //fermionic couplings for phi3
    yu3 = Rij_GTHDM(2,0) + (Rij_GTHDM(2,1) - i*Rij_GTHDM(2,2))*su.conjugate();
    yd3 = Rij_GTHDM(2,0) + (Rij_GTHDM(2,1) + i*Rij_GTHDM(2,2))*sd;
    yl3 = Rij_GTHDM(2,0) + (Rij_GTHDM(2,1) + i*Rij_GTHDM(2,2))*sl;

    yu3R = Rij_GTHDM(2,0) + (Rij_GTHDM(2,1))*su.real();
    yd3R = Rij_GTHDM(2,0) + (Rij_GTHDM(2,1))*sd.real();
    yl3R = Rij_GTHDM(2,0) + (Rij_GTHDM(2,1))*sl.real();

    //These cross sections ratios are necessary for rphi3_gg
    //At 8 TeV
    
    //SM gg -> H (phi even) production cross section ratio at 8 TeV, charm and top loop only over total
    double rSigmaggphi3E_tc8 = ip_csr_ggH_tc_8(m3);
    //SM gg -> H (phi even) production cross section ratio at 8 TeV, bottom loop only over total
    double rSigmaggphi3E_b8 = ip_csr_ggH_b_8(m3);
    //gg -> H (phieven) production cross section at 8 TeV, total
   // double Sigmaggphi3E_8 = ip_cs_ggtoH_8(m3);
    
    
    //gg -> A (phi odd) production cross section ratio at 8 TeV, charm and top loop only over total
    double rSigmaggphi3O_tc8 = ip_csr_ggA_tc_8(m3);
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
    rphi3_ggE = yu3.real()*yd3.real() + (yu3.real()*yu3.real() - yu3.real()*yd3.real())*rSigmaggphi3E_tc8  + (yd3.real()*yd3.real() - yu3.real()*yd3.real())*rSigmaggphi3E_b8;
    rphi3_ggO = yu3.imag()*yd3.imag() + (yu3.imag()*yu3.imag() - yu3.imag()*yd3.imag())*rSigmaggphi3O_tc8  + (yd3.imag()*yd3.imag() - yu3.imag()*yd3.imag())*rSigmaggphi3O_b8;
    rphi3_VV=Rij_GTHDM(2,0)*Rij_GTHDM(2,0);
     
      
    /*Gamma_phi3gaga and Gamma_phi3Zga expressions ...*/
    
    /*Decay to photons. The fermionic contribution has a CP-even part (HH) and a CP-odd (A)*/
    /*CP EVEN*/
                                                                                   
    gslpp::complex I_HH_Ux=I_HH_U(m3_2,Mc,Mt);
    gslpp::complex I_HH_Dx=I_HH_D(m3_2,Ms,Mb);
    gslpp::complex I_HH_Lx=I_HH_L(m3_2,Mmu,Mtau);
    gslpp::complex I_phi3E_F= yu3.real()*I_HH_Ux+ yd3.real()*I_HH_Dx+yl3.real()*I_HH_Lx;
                                                                               
    gslpp::complex I_phi3_W=Rij_GTHDM(2,0)*I_H_W(m3,MW);


                          
    double lambdaphi3HpHm =  lambdaipm(R31, R32, R33);
    gslpp::complex I_phi3_Hp=(vev*vev)/(2.0*mHp2)*I_H_Hp(mHp2,m3)*(lambdaphi3HpHm);
    
    
    /*CP ODD */
            
    gslpp::complex I_A_Ux=I_A_U(m3_2,Mc,Mt);
    gslpp::complex I_A_Dx=I_A_D(m3_2,Ms,Mb);
    gslpp::complex I_A_Lx=I_A_L(m3_2,Mmu,Mtau);
                                                                               
    gslpp::complex I_phi3O_F = yu3.imag()*I_A_Ux + yd3.imag()*I_A_Dx + yl3.imag()*I_A_Lx;
                                                                             
    double Gamma_phi3gaga=(GF*Ale*Ale*m3*m3*m3/(sqrt(2.0)*128.0*M_PI*M_PI*M_PI))*((I_phi3E_F+I_phi3_W+I_phi3_Hp).abs2() + (I_phi3O_F).abs2());


    /*Decay to Z gamma
    CP-EVEN PART*/

    gslpp::complex A_HH_Ux = A_HH_U(m3_2,cW2,Mc,Mt,MZ);
    gslpp::complex A_HH_Dx = A_HH_D(m3_2,cW2,Ms,Mb,MZ);
    gslpp::complex A_HH_Lx = A_HH_L(m3_2,cW2,Mmu,Mtau,MZ);
                                                                               
    gslpp::complex A_phi3E_F = (yu3.real()*A_HH_Ux+ yd3.real()*A_HH_Dx+ yl3.real()*A_HH_Lx)/sqrt(sW2*cW2);
    gslpp::complex A_phi3_W = Rij_GTHDM(2,0)*A_H_W(m3,cW2,MW,MZ);
  
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
    //std::cout<<"\033[1;33m   ip_cs_ggtoH_8(m3) =\033[0m "<< ip_cs_ggtoH_8(m3) <<std::endl;
    //std::cout<<"\033[1;33m   rphi3_ggE =\033[0m "<< rphi3_ggE <<std::endl;
    //std::cout<<"\033[1;33m   ip_cs_ggtoA_8(m3) =\033[0m "<< ip_cs_ggtoA_8(m3) <<std::endl;
    //std::cout<<"\033[1;33m   rphi3_ggO =\033[0m "<< rphi3_ggO <<std::endl;
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
double BrSM_phi3tomumu=ip_Br_HPtomumu(m3);
double BrSM_phi3toWW =ip_Br_HPtoWW(m3);
double BrSM_phi3toZZ =ip_Br_HPtoZZ(m3);

Gammaphi3totSM=ip_GammaHPtotSM(m3);

 
 
/*Decay of phi3 to the others scalars*/
double lambda123 = lambdaijk(R11, R12, R13, R21, R22, R23, R31, R32, R33, lambda1, lambda3, lambda4, Relambda5,
                             Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7) +
                   lambdaijk(R21, R22, R23, R11, R12, R13, R31, R32, R33, lambda1, lambda3, lambda4, Relambda5,
                             Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7) +
                   lambdaijk(R21, R22, R23, R31, R32, R33, R11, R12, R13, lambda1, lambda3, lambda4, Relambda5,
                             Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7) +
                   lambdaijk(R11, R12, R13, R31, R32, R33, R21, R22, R23, lambda1, lambda3, lambda4, Relambda5,
                             Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7) +
                   lambdaijk(R31, R32, R33, R11, R12, R13, R21, R22, R23, lambda1, lambda3, lambda4, Relambda5,
                             Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7) +
                   lambdaijk(R31, R32, R33, R21, R22, R23, R11, R12, R13, lambda1, lambda3, lambda4, Relambda5,
                             Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7);
double lambda223 = lambdaijk(R21, R22, R23, R21, R22, R23, R31, R32, R33, lambda1, lambda3, lambda4, Relambda5,
                             Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7) +
                   lambdaijk(R21, R22, R23, R31, R32, R33, R21, R22, R23, lambda1, lambda3, lambda4, Relambda5,
                             Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7) +
                   lambdaijk(R31, R32, R33, R21, R22, R23, R21, R22, R23, lambda1, lambda3, lambda4, Relambda5,
                             Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7);
double lambda113 = lambdaijk(R11, R12, R13, R11, R13, R13, R31, R32, R33, lambda1, lambda3, lambda4, Relambda5,
                             Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7) +
                   lambdaijk(R11, R11, R13, R31, R32, R33, R11, R12, R13, lambda1, lambda3, lambda4, Relambda5,
                             Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7) +
                   lambdaijk(R31, R32, R33, R11, R12, R13, R11, R12, R13, lambda1, lambda3, lambda4, Relambda5,
                             Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7);

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

double Gammaphi3_tt = HSTheta(m3 - 2*Mt)*BrSM_phi3tott*(rphi3_QuQuE + rphi3_QuQuO/(beta(Mt, m3_2)*beta(Mt, m3_2)))*Gammaphi3totSM;

double Gammaphi3_cc = HSTheta(m3 - 2*Mc)*BrSM_phi3tocc*(rphi3_QuQuE + rphi3_QuQuO/(beta(Mc, m3_2)*beta(Mc, m3_2)))*Gammaphi3totSM;

double Gammaphi3_bb = HSTheta(m3 - 2*Mb)*BrSM_phi3tobb*(rphi3_QdQdE + rphi3_QdQdO/(beta(Mb, m3_2)*beta(Mb, m3_2)))*Gammaphi3totSM;

double Gammaphi3_tautau = HSTheta(m3 - 2*Mtau)*BrSM_phi3totautau*(rphi3_QlQlE + rphi3_QlQlO/(beta(Mtau, m3_2)*beta(Mtau, m3_2)))*Gammaphi3totSM;

double Gammaphi3_mumu = HSTheta(m3 - 2*Mmu)*BrSM_phi3tomumu*(rphi3_QlQlE + rphi3_QlQlO/(beta(Mmu, m3_2)*beta(Mmu, m3_2)))*Gammaphi3totSM;

double Gammaphi3_WW = HSTheta(m3 - 2*MW)*BrSM_phi3toWW*(rphi3_VV)*Gammaphi3totSM;

double Gammaphi3_ZZ = HSTheta(m3 - 2*MZ)*BrSM_phi3toZZ*(rphi3_VV)*Gammaphi3totSM;


Gammaphi3tot  = 1.e-15;

Gammaphi3tot = Gammaphi3tot + Gammaphi3_tt + Gammaphi3_cc + Gammaphi3_bb
             + Gammaphi3_tautau + Gammaphi3_mumu + Gammaphi3_WW
             + Gammaphi3_ZZ + Gamma_phi3gaga + Gamma_phi3Zga + Gamma_phi3gg
             + Gammaphi3_phi1phi1 + Gammaphi3_phi2phi2 + Gammaphi3_phi1phi2
             + Gammaphi3_HpHm + Gammaphi3_phi1Z + Gammaphi3_phi2Z
             + Gammaphi3_HpW;



                
    //std::cout<<"\033[1;34m   R31 =\033[0m "<< R31 <<std::endl;   
    //std::cout<<"\033[1;34m   R32 =\033[0m "<< R32 <<std::endl;   
    //std::cout<<"\033[1;34m   R33 =\033[0m "<< R33 <<std::endl;   
            
    //std::cout<<"\033[1;34m   lambdaphi3HpHm =\033[0m "<< lambdaphi3HpHm <<std::endl;
    //std::cout<<"\033[1;34m   lambdaphi3HpHm =\033[0m "<< lambdaphi3HpHm <<std::endl;
    //std::cout<<"\033[1;34m   Gammaphi3_HpHm =\033[0m "<< Gammaphi3_HpHm <<std::endl;
    //std::cout<<"\033[1;34m   Gammaphi3_HpW =\033[0m "<< Gammaphi3_HpW <<std::endl;
    
    





    Br_phi3tott=Gammaphi3_tt/Gammaphi3tot;
    Br_phi3tobb=Gammaphi3_bb/Gammaphi3tot;
    Br_phi3tomumu=Gammaphi3_mumu/Gammaphi3tot;
    Br_phi3totautau=Gammaphi3_tautau/Gammaphi3tot;
    Br_phi3toWW=Gammaphi3_WW/Gammaphi3tot;
    Br_phi3toZZ=Gammaphi3_ZZ/Gammaphi3tot;
    Br_phi3togaga=Gamma_phi3gaga/Gammaphi3tot;
    Br_phi3togg=Gamma_phi3gg/Gammaphi3tot;
    Br_phi3toZga=Gamma_phi3Zga/Gammaphi3tot;
    Br_phi3tophi1phi1=Gammaphi3_phi1phi1/Gammaphi3tot;
    Br_phi3tophi2phi2=Gammaphi3_phi2phi2/Gammaphi3tot;
    Br_phi3tophi1phi2=Gammaphi3_phi1phi2/Gammaphi3tot;
    Br_phi3toHpHm=Gammaphi3_HpHm/Gammaphi3tot;
    Br_phi3tophi1Z=Gammaphi3_phi1Z/Gammaphi3tot;
    Br_phi3tophi2Z=Gammaphi3_phi2Z/Gammaphi3tot;
    Br_phi3toHpW=Gammaphi3_HpW/Gammaphi3tot;
    
    
    //std::cout<<"\033[1;34m   BrSM_phi3totautau =\033[0m "<< BrSM_phi3totautau <<std::endl;
}


void GeneralTHDMcache::computeHpquantities()
{
    //FLAG to select only the model in which all the couplings are the same (by families)
    if (!myGTHDM->getATHDMflag())
    {
        throw std::runtime_error("Direct Searches are only available in the A2HDM.");
    }

    /*complex i */
    gslpp::complex i = gslpp::complex::i();

    // Pole masses and CKM matrix elements relevant for Hp quantities
    double alsmuc = myGTHDM->Als(myGTHDM->getMuc(), FULLNNNLO, true, true);
    // Mbar2Mp does not receive Mc
    double Mcp = Mc*(1.+4.*alsmuc/3./M_PI+alsmuc*alsmuc/M_PI/M_PI*(-1.0414*(1.-4.*Ms/3.*Mc)+13.4434));
    double Mbp = myGTHDM->Mbar2Mp(Mb, QCD::BOTTOM, FULLNNLO);
    double Mtp = myGTHDM->getMtpole();
    double Mmu2  = Mmu*Mmu;
    double Mtau2 = Mtau*Mtau;
    double Ms2 = Ms*Ms;
    double Mc2 = Mcp*Mcp;
    double Mb2 = Mbp*Mbp;
    double Mt2 = Mtp*Mtp;
    double r2hot = mHp*mHp/Mtp/Mtp;
    double r2bot = Mbp*Mbp/Mtp/Mtp;
    double Vtb = myGTHDM->getCKM().getV_tb().abs();
    double Vcb = myGTHDM->getCKM().getV_cb().abs();
    double Vcs = myGTHDM->getCKM().getV_cs().abs();

    //In order to compute the xsection we use the xsections generated in the table log_cs_ggtoHp_8
    //such a table is generated for the type II mode, basically we take the value for tanb=1 and
    //then we rescale with the coupling of the top-quark, there should be a residual dependence on 
    //the coupling to the bottom-quarks which we neglect since it should be proportional to the bottom
    //quark (and the one included to the top quark)
    SigmaHp8=0.0;
    SigmaHpm13=0.0;
        
    SigmaHp8=ip_cs_ggtoHp_8(mHp,0.0)*su.abs2();
    SigmaHpm13=ip_cs_ggtoHp_13(mHp,0.0)*su.abs2();


    double GammaHpmunu  = HSTheta(mHp-Mmu)  * (Mmu2  * (mHp2-Mmu2)  * (mHp2-Mmu2)  * sl.abs2()) / (8.0*mHp*mHp2*M_PI*vev*vev);
    double GammaHptaunu = HSTheta(mHp-Mtau) * (Mtau2 * (mHp2-Mtau2) * (mHp2-Mtau2) * sl.abs2()) / (8.0*mHp*mHp2*M_PI*vev*vev);

    double GammaHpcs = 0.;
    double GammaHpcb = 0.;
    double GammaHptb = 0.;
    double GammaHpHlW = 0.;
    double GammaHpphi2W = 0.;
    double GammaHpphi3W = 0.;
    double GammaHpphi3Wst = 0.;
    
    if(mHp>=Mtp+Mbp)
    {
        GammaHptb = HSTheta(mHp-Mtp-Mbp)*(Vtb*Vtb/(8.0*mHp*M_PI*vev*vev))*3.0*(-4.0*(su*sd).real()*Mb2*Mt2
                    -sd.abs2()*Mb2*(Mb2-mHp2+Mt2)-su.abs2()*Mt2*(Mb2-mHp2+Mt2))
                    *sqrt((Mb2*Mb2+(mHp2-Mt2)*(mHp2-Mt2)-2.0*Mb2*(mHp2+Mt2))/(mHp2*mHp2));
    }

    if(mHp>=Mcp+Mbp)
    {
        GammaHpcb = HSTheta(mHp-Mcp-Mbp)*(Vcb*Vcb/(8.0*mHp*M_PI*vev*vev))*3.0*(-4.0*(su*sd).real()*Mb2*Mc2
                    -sd.abs2()*Mb2*(Mb2-mHp2+Mc2)-su.abs2()*Mc2*(Mb2-mHp2+Mc2))
                    *sqrt((Mb2*Mb2+(mHp2-Mc2)*(mHp2-Mc2)-2.0*Mb2*(mHp2+Mc2))/(mHp2*mHp2));
    }

    if(mHp>=Mcp+Ms)
    {
        GammaHpcs = HSTheta(mHp-Mcp-Ms)*(Vcs*Vcs/(8.0*mHp*M_PI*vev*vev))*3.0*(-4.0*(su*sd).real()*Ms2*Mc2
                    -sd.abs2()*Ms2*(Ms2-mHp2+Mc2)-su.abs2()*Mc2*(Ms2-mHp2+Mc2))
                    *sqrt((Ms2*Ms2+(mHp2-Mc2)*(mHp2-Mc2)-2.0*Ms2*(mHp2+Mc2))/(mHp2*mHp2));
    }

    //decay to phi1 (SM Higgs)
    if(mHp>=mHl+MW)
    {
        GammaHpHlW = KaellenFunction(1.0, mHl*mHl/mHp/mHp, MW*MW/mHp/mHp) *
                     KaellenFunction(1.0, mHl*mHl/mHp/mHp, MW*MW/mHp/mHp) *
                     KaellenFunction(1.0, mHl*mHl/mHp/mHp, MW*MW/mHp/mHp) *
                     mHp*mHp*mHp * (R12*R12 + R13*R13) / (2.0*M_PI*vev*vev);
    }

    //decay to phi2
    if(mHp>=m2+MW)
    {
        GammaHpphi2W = KaellenFunction(1.0, m2*m2/mHp/mHp, MW*MW/mHp/mHp) *
                       KaellenFunction(1.0, m2*m2/mHp/mHp, MW*MW/mHp/mHp) *
                       KaellenFunction(1.0, m2*m2/mHp/mHp, MW*MW/mHp/mHp) *
                       mHp*mHp*mHp * (R22*R22 + R23*R23) / (2.0*M_PI*vev*vev);
    }

    //decay to phi3
    if(mHp>=m3+MW)
    {
        GammaHpphi3W = KaellenFunction(1.0, m3*m3/mHp/mHp, MW*MW/mHp/mHp) *
                       KaellenFunction(1.0, m3*m3/mHp/mHp, MW*MW/mHp/mHp) *
                       KaellenFunction(1.0, m3*m3/mHp/mHp, MW*MW/mHp/mHp) *
                       mHp*mHp*mHp * (R32*R32 + R33*R33) / (2.0*M_PI*vev*vev);
    }

    //decay to phi3 + Wstar
    //computed for LEP interval of Hp and phi3 masses
    if(mHp >= 40.0 && mHp <= 93.0 && mH3 >= 12.0 && mH3 <= (mHp - 3.))
        GammaHpphi3Wst = ip_integral_Hp_A_Wstar(mHp,mH3,1.,1.);

    GammaHptot = 1.e-10;

    GammaHptot = GammaHptot + GammaHpmunu  + GammaHptaunu +
                 GammaHpcs  + GammaHpcb    + GammaHptb    +
                 GammaHpHlW + GammaHpphi2W + GammaHpphi3W + GammaHpphi3Wst;


    Br_Hptomunu    = GammaHpmunu/GammaHptot;
    Br_Hptotaunu   = GammaHptaunu/GammaHptot;
    Br_Hptocs      = GammaHpcs/GammaHptot;
    Br_Hptocb      = GammaHpcb/GammaHptot;
    Br_Hptotb      = GammaHptb/GammaHptot;
    Br_Hptophi3W   = GammaHpphi3W/GammaHptot;
    Br_Hptophi3Wst = GammaHpphi3Wst/GammaHptot;
    
    //HSTheta(mHp-Mt-Mb)*(Vtb*Vtb/(8.0*mHp*M_PI*vev*vev))*3.0*(-4.0*(su*sd).real()*Mb2*Mt2
    //                    -sd.abs2()*Mb2*(Mb2-mHp2+Mt2)-su.abs2()*Mt2*(Mb2-mHp2+Mt2))
    //                  *sqrt((Mb2*Mb2+(mHp2-Mt2)*(mHp2-Mt2)-2.0*Mb2*(mHp2+Mt2))/(mHp2*mHp2))


    // Hp production from top decay, for light charged scalars
    // Formula according to Czarnecki and Davidson, hep-ph/9301237
    //double Gammatoptot = 1.42; // PDG 2023
    double GammatHpb = 0.;
    double pre_tHpb = GF * Mt2 * Mtp * Vtb * Vtb / 4. / sqrt(2.) / M_PI;
    gslpp::complex aCD = su - sd * Mbp / Mtp;
    gslpp::complex bCD = su + sd * Mbp / Mtp;
    double aCD2 = (aCD * aCD.conjugate()).real();
    double bCD2 = (bCD * bCD.conjugate()).real();

    if(Mtp>=mHp+Mbp)
    {
        GammatHpb = pre_tHpb * ((1. - r2hot + r2bot) * (aCD2 + bCD2) / 2. + r2bot * (aCD2 - bCD2)) *
                    sqrt(1. + r2hot*r2hot + r2bot*r2bot - 2. * (r2hot + r2bot + r2hot*r2bot)) / 2.;
    }

    double GammatWb = 1.35; // (in GeV) at NLO-QCD from 2409.04179

    Gammatoptot = GammatWb + GammatHpb;

    Br_ttoHpb = GammatHpb / Gammatoptot;
}


void GeneralTHDMcache::computeWZquantities()
{
    double MZ2 = MZ * MZ;
    double MW2 = MW * MW;
    double preZ = GF * MZ * MZ2 / 24. / sqrt(2.) / M_PI;
    double preW = GF * MW * MW2 / 24. / sqrt(2.) / M_PI;

    double gZAH_2   = R22 * R22; // CP conserving case; needs to be generalised for CPV analyses
    double gZHpHm_2 = (sW2 - cW2) * (sW2 - cW2);
    double gWHpH_2  = R22 * R22; // CP conserving case; needs to be generalised for CPV analyses
    double gWHpA_2  = R33 * R33; // CP conserving case; needs to be generalised for CPV analyses

    double func_ZAH   = 1. + (mH3sq - mH2sq) * (mH3sq - mH2sq) / MZ2 / MZ2 - 2. * (mH3sq + mH2sq) / MZ2;
    double func_ZHpHm = 1. - 4. * mHp2 / MZ2;
    double func_WHpH  = 1. + (mHp2  - mH2sq) * (mHp2  - mH2sq) / MW2 / MW2 - 2. * (mHp2  + mH2sq) / MW2;
    double func_WHpA  = 1. + (mHp2  - mH3sq) * (mHp2  - mH3sq) / MW2 / MW2 - 2. * (mHp2  + mH3sq) / MW2;

    double Gamma_ZAH   = HSTheta(MZ-mH3-mH2) * preZ * gZAH_2   * sqrt(std::fabs(func_ZAH))   * func_ZAH;
    double Gamma_ZHpHm = HSTheta(MZ-2.*mHp)  * preZ * gZHpHm_2 * sqrt(std::fabs(func_ZHpHm)) * func_ZHpHm;
    double Gamma_WHpH  = HSTheta(MW-mHp-mH2) * preW * gWHpH_2  * sqrt(std::fabs(func_WHpH))  * func_WHpH;
    double Gamma_WHpA  = HSTheta(MW-mHp-mH3) * preW * gWHpA_2  * sqrt(std::fabs(func_WHpA))  * func_WHpA;

    Gamma_Z_inv = myGTHDM->Gamma_inv() + Gamma_ZAH + Gamma_ZHpHm;

    Gamma_W_inv = Gamma_WHpH + Gamma_WHpA;
}


//Higgs direct searches

void GeneralTHDMcache::computeHeavyHiggs()
{
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
//    double Br_Ztohadrons=0.69911; //PDG2022
    double Br_Wtoenu=0.1071; //C. Patrignani et al.(Particle Data Group), Chin. Phys. C, 40, 100001 (2016)
    double Br_Wtomunu=0.1063; //C. Patrignani et al.(Particle Data Group), Chin. Phys. C, 40, 100001 (2016)
    double Br_Wtotaunu=0.1138; //C. Patrignani et al.(Particle Data Group), Chin. Phys. C, 40, 100001 (2016)
    double Br_Wtohadrons=0.6741; // PDG2022
    double Br_tautoleptons=0.3521; // PDG2022
    double Br_tautohadrons=1-Br_tautoleptons; // PDG2022
    
    
    THoEX_bb_phi2_bb_ATLAS13=0.0;
    THoEX_bb_phi3_bb_ATLAS13=0.0;
    THoEX_tt_phi2_tt_ATLAS13=0.0;
    THoEX_tt_phi3_tt_ATLAS13=0.0;
    THoEX_bb_phi2_tt_ATLAS13=0.0;
    THoEX_bb_phi3_tt_ATLAS13=0.0;
    THoEX_bb_phi2_bb_CMS8=0.0;
    THoEX_bb_phi3_bb_CMS8=0.0;
    THoEX_gg_phi2_bb_CMS8=0.0;
    THoEX_gg_phi3_bb_CMS8=0.0;
    THoEX_tt_phi2_tt_CMS13=0.0;
    THoEX_tt_phi3_tt_CMS13=0.0;
    THoEX_pp_phi3_bb_light_CMS13=0.0;
    THoEX_pp_phi2_bb_light_CMS13=0.0;
    THoEX_pp_phi2_bb_CMS13=0.0;
    THoEX_pp_phi3_bb_CMS13=0.0;
    THoEX_bb_phi2_bb_CMS13=0.0;
    THoEX_bb_phi3_bb_CMS13=0.0;
    
    
    THoEX_gg_phi2_mumu_CMS8=0.0;
    THoEX_gg_phi3_mumu_CMS8=0.0;
    THoEX_bb_phi2_mumu_CMS8=0.0;
    THoEX_bb_phi3_mumu_CMS8=0.0;
    
    THoEX_gg_phi2_mumu_CMS13=0.0;
    THoEX_gg_phi3_mumu_CMS13=0.0;
    THoEX_bb_phi2_mumu_CMS13=0.0;
    THoEX_bb_phi3_mumu_CMS13=0.0;
    
    THoEX_gg_phi2_mumu_ATLAS13=0.0;
    THoEX_gg_phi3_mumu_ATLAS13=0.0;
    THoEX_bb_phi2_mumu_ATLAS13=0.0;
    THoEX_bb_phi3_mumu_ATLAS13=0.0;
            
    
    
    
    
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
    
    
    THoEX_gg_phi2_WW_heavy_CMS13=0.0;
    THoEX_gg_phi3_WW_heavy_CMS13=0.0;
    THoEX_VV_phi2_WW_heavy_CMS13=0.0;
    THoEX_VV_phi3_WW_heavy_CMS13=0.0;
    
    
    THoEX_gg_phi2_WW_CMS13=0.0;
    THoEX_gg_phi3_WW_CMS13=0.0;
    THoEX_VV_phi2_WW_CMS13=0.0;
    THoEX_VV_phi3_WW_CMS13=0.0;
    
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
    THoEX_pp_phi2_VV_CMS8=0.0;
    THoEX_pp_phi3_VV_CMS8=0.0;
    THoEX_pp_phi2_VV_qqqq_ATLAS13=0.0;
    THoEX_pp_phi3_VV_qqqq_ATLAS13=0.0;
    
    THoEX_gg_phi2_VV_llqq_ATLAS13=0.0;
    THoEX_gg_phi3_VV_llqq_ATLAS13=0.0;
    THoEX_VV_phi2_VV_llqq_ATLAS13=0.0;
    THoEX_VV_phi3_VV_llqq_ATLAS13=0.0;
    
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
    //THoEX_pp_phi2_phi1phi1_bbtautau_ATLAS13=0.0;        //OLD this has been splitted in two
    THoEX_pp_phi2_phi1phi1_bbtautau_1_ATLAS13=0.0;
    THoEX_pp_phi2_phi1phi1_bbtautau_2_ATLAS13=0.0;
    //THoEX_pp_phi3_phi1phi1_bbtautau_ATLAS13=0.0;      //OLD this has been splitted in two
    THoEX_pp_phi3_phi1phi1_bbtautau_1_ATLAS13=0.0;
    THoEX_pp_phi3_phi1phi1_bbtautau_2_ATLAS13=0.0;
    THoEX_pp_phi2_phi1phi1_bbtautau_1_CMS13=0.0;
    THoEX_pp_phi3_phi1phi1_bbtautau_1_CMS13=0.0;
    THoEX_pp_phi2_phi1phi1_bbtautau_2_CMS13=0.0;
    THoEX_pp_phi3_phi1phi1_bbtautau_2_CMS13=0.0;
    THoEX_pp_phi2_phi1phi1_bbVV_CMS13=0.0;
    THoEX_pp_phi3_phi1phi1_bbVV_CMS13=0.0;
    
    THoEX_pp_phi2_phi1phi1_4WOr2W2tauOr4tau_CMS13=0.0;
    THoEX_pp_phi3_phi1phi1_4WOr2W2tauOr4tau_CMS13=0.0;
    
    
    THoEX_pp_phi2_phi1phi1_bbWW_qqlnu_CMS13=0.0;
    THoEX_pp_phi3_phi1phi1_bbWW_qqlnu_CMS13=0.0;
    
    
    
    THoEX_pp_phi2_phi1phi1_bbZZ_lljj_CMS13=0.0;
    THoEX_pp_phi3_phi1phi1_bbZZ_llnunu_CMS13=0.0;
    THoEX_pp_phi2_phi1phi1_bbZZ_lljj_CMS13=0.0;
    THoEX_pp_phi3_phi1phi1_bbZZ_llnunu_CMS13=0.0;
    
    
    THoEX_pp_phi2_phi1phi1_bbWWorbbtautau_CMS13=0.0;
    THoEX_pp_phi3_phi1phi1_bbWWorbbtautau_CMS13=0.0;
            
            
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
    
    THoEX_gg_phi2_phi1Z_tautaull_CMS13=0.0;
    THoEX_gg_phi3_phi1Z_tautaull_CMS13=0.0;
    
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
    THoEX_gg_phi3_phi2Z_WWZ_ATLAS13=0.0;
    THoEX_gg_phi2_phi3Z_WWZ_ATLAS13=0.0;  
    THoEX_bb_phi3_phi2Z_bbZ_ATLAS13=0.0;
    THoEX_bb_phi2_phi3Z_bbZ_ATLAS13=0.0;
    
    THoEX_pp_Hpm_taunu_ATLAS8=0.0;
    THoEX_pp_Hp_taunu_CMS8=0.0;
    THoEX_pp_Hpm_taunu_ATLAS13=0.0;
    THoEX_pp_Hpm_taunu_CMS13=0.0;
    THoEX_pp_Hpm_tb_ATLAS8=0.0;
    THoEX_pp_Hp_tb_CMS8=0.0;
    THoEX_pp_Hpm_tb_ATLAS13=0.0;
    THoEX_pp_Hpm_tb_CMS13=0.0;

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
    pp_phi3_bb_TH13=SigmaggF_phi3_13*Br_phi3tobb;
    bb_phi2_bb_TH13=SigmabbF_phi2_13*Br_phi2tobb;
    bb_phi3_bb_TH13=SigmabbF_phi3_13*Br_phi3tobb;
    
    
    
    gg_phi2_mumu_TH8=SigmaggF_phi2_8*Br_phi2tomumu;
    gg_phi3_mumu_TH8=SigmaggF_phi3_8*Br_phi3tomumu;
    bb_phi2_mumu_TH8=SigmabbF_phi2_8*Br_phi2tomumu;
    bb_phi3_mumu_TH8=SigmabbF_phi3_8*Br_phi3tomumu;
    gg_phi2_mumu_TH13=SigmaggF_phi2_13*Br_phi2tomumu;
    gg_phi3_mumu_TH13=SigmaggF_phi3_13*Br_phi3tomumu;
    bb_phi2_mumu_TH13=SigmabbF_phi2_13*Br_phi2tomumu;
    bb_phi3_mumu_TH13=SigmabbF_phi3_13*Br_phi3tomumu;
    
    
    
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
    
    

    


    pp_phi2_VV_TH8=SigmaSumphi2_8/SigmaTotSM_phi2_8*rphi2_VV*Gammaphi2totSM/Gammaphi2tot;
    pp_phi3_VV_TH8=SigmaSumphi3_8/SigmaTotSM_phi3_8*rphi3_VV*Gammaphi3totSM/Gammaphi3tot;
    pp_phi2_VV_TH13=SigmaSumphi2_13*(Br_phi2toZZ+Br_phi2toWW);
    pp_phi3_VV_TH13=SigmaSumphi3_13*(Br_phi3toZZ+Br_phi3toWW);
    
    gg_phi2_VV_TH13=SigmaggF_phi2_13*(Br_phi2toZZ+Br_phi2toWW);
    gg_phi3_VV_TH13=SigmaggF_phi3_13*(Br_phi3toZZ+Br_phi3toWW);

    VV_phi2_VV_TH13=SigmaVBF_phi2_13*(Br_phi2toZZ+Br_phi2toWW);
    VV_phi3_VV_TH13=SigmaVBF_phi3_13*(Br_phi3toZZ+Br_phi3toWW);
    
    


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
    pp_phi2_phi1phi1_with_channel_bbbb_TH13=SigmaSumphi2_13*Br_phi2tophi1phi1*GTHDM_BR_h_bb*GTHDM_BR_h_bb/(BrSM_htobb*BrSM_htobb);
    pp_phi3_phi1phi1_with_channel_bbbb_TH13=SigmaSumphi3_13*Br_phi3tophi1phi1*GTHDM_BR_h_bb*GTHDM_BR_h_bb/(BrSM_htobb*BrSM_htobb);
    
    pp_phi2_phi1phi1_TH13=SigmaSumphi2_13*Br_phi2tophi1phi1;
    pp_phi3_phi1phi1_TH13=SigmaSumphi3_13*Br_phi3tophi1phi1;
    
    pp_phi2_phi1phi1_bbgaga_TH13=SigmaSumphi2_13*Br_phi2tophi1phi1*GTHDM_BR_h_bb*GTHDM_BR_h_gaga*2.0;
    pp_phi3_phi1phi1_bbgaga_TH13=SigmaSumphi3_13*Br_phi3tophi1phi1*GTHDM_BR_h_bb*GTHDM_BR_h_gaga*2.0;
    pp_phi2_phi1phi1_with_channel_bbgaga_TH13=SigmaSumphi2_13*Br_phi2tophi1phi1*GTHDM_BR_h_bb*GTHDM_BR_h_gaga/(BrSM_htobb*BrSM_htogaga);
    pp_phi3_phi1phi1_with_channel_bbgaga_TH13=SigmaSumphi3_13*Br_phi3tophi1phi1*GTHDM_BR_h_bb*GTHDM_BR_h_gaga/(BrSM_htobb*BrSM_htogaga);
    
    pp_phi2_phi1phi1_bbtautau_TH13=SigmaSumphi2_13*Br_phi2tophi1phi1*GTHDM_BR_h_bb*GTHDM_BR_h_tautau*2.0;
    pp_phi3_phi1phi1_bbtautau_TH13=SigmaSumphi3_13*Br_phi3tophi1phi1*GTHDM_BR_h_bb*GTHDM_BR_h_tautau*2.0;
    
    pp_phi2_phi1phi1_with_channel_bbtautau_TH13=SigmaSumphi2_13*Br_phi2tophi1phi1*GTHDM_BR_h_bb*GTHDM_BR_h_tautau/(BrSM_htobb*BrSM_htotautau);
    pp_phi3_phi1phi1_with_channel_bbtautau_TH13=SigmaSumphi3_13*Br_phi3tophi1phi1*GTHDM_BR_h_bb*GTHDM_BR_h_tautau/(BrSM_htobb*BrSM_htotautau);
    
    
    pp_phi2_phi1phi1_bbVV_TH13=SigmaSumphi2_13*Br_phi2tophi1phi1*2.0*GTHDM_BR_h_bb*(GTHDM_BR_h_WW*pow(Br_Wtoenu+Br_Wtomunu+Br_Wtotaunu*Br_tautoleptons,2)
                                                            +GTHDM_BR_h_ZZ*2.0*Br_Ztoinv*(Br_Ztoee+Br_Ztomumu+Br_Ztotautau*Br_tautoleptons*Br_tautoleptons));
    pp_phi3_phi1phi1_bbVV_TH13=SigmaSumphi3_13*Br_phi3tophi1phi1*2.0*GTHDM_BR_h_bb*(GTHDM_BR_h_WW*pow(Br_Wtoenu+Br_Wtomunu+Br_Wtotaunu*Br_tautoleptons,2)
                                                            +GTHDM_BR_h_ZZ*2.0*Br_Ztoinv*(Br_Ztoee+Br_Ztomumu+Br_Ztotautau*Br_tautoleptons*Br_tautoleptons));  
    
    
            
// I'd say that in the paper give values for \phi->HH and assume SM for the other decays so we correct those final dicays but divide by the 
// BR of the SM (which is what experimentalist have assumed)
    pp_phi2_phi1phi1_with_channel_4WOr2W2tauOr4tau_TH13=SigmaSumphi2_13*Br_phi2tophi1phi1*(GTHDM_BR_h_WW*GTHDM_BR_h_WW*pow((Br_Wtoenu+Br_Wtomunu+Br_Wtotaunu*Br_tautohadrons),4)
            +2.0*GTHDM_BR_h_WW*GTHDM_BR_h_tautau*Br_tautohadrons*Br_tautohadrons*(Br_Wtoenu+Br_Wtomunu+Br_Wtotaunu*Br_tautohadrons)*(Br_Wtoenu+Br_Wtomunu+Br_Wtotaunu*Br_tautohadrons)
            +GTHDM_BR_h_tautau*GTHDM_BR_h_tautau*Br_tautohadrons*Br_tautohadrons*Br_tautohadrons*Br_tautohadrons)/
            (BrSM_htoWW*BrSM_htoWW*pow((Br_Wtoenu+Br_Wtomunu+Br_Wtotaunu*Br_tautohadrons),4)
            +2.0*BrSM_htoWW*BrSM_htotautau*Br_tautohadrons*Br_tautohadrons*(Br_Wtoenu+Br_Wtomunu+Br_Wtotaunu*Br_tautohadrons)*(Br_Wtoenu+Br_Wtomunu+Br_Wtotaunu*Br_tautohadrons)
            +BrSM_htotautau*BrSM_htotautau*Br_tautohadrons*Br_tautohadrons*Br_tautohadrons*Br_tautohadrons);
    
    pp_phi3_phi1phi1_with_channel_4WOr2W2tauOr4tau_TH13=SigmaSumphi3_13*Br_phi3tophi1phi1*(GTHDM_BR_h_WW*GTHDM_BR_h_WW*pow((Br_Wtoenu+Br_Wtomunu+Br_Wtotaunu*Br_tautohadrons),4)
            +2.0*GTHDM_BR_h_WW*GTHDM_BR_h_tautau*Br_tautohadrons*Br_tautohadrons*(Br_Wtoenu+Br_Wtomunu+Br_Wtotaunu*Br_tautohadrons)*(Br_Wtoenu+Br_Wtomunu+Br_Wtotaunu*Br_tautohadrons)
            +GTHDM_BR_h_tautau*GTHDM_BR_h_tautau*Br_tautohadrons*Br_tautohadrons*Br_tautohadrons*Br_tautohadrons)/
            (BrSM_htoWW*BrSM_htoWW*pow((Br_Wtoenu+Br_Wtomunu+Br_Wtotaunu*Br_tautohadrons),4)
            +2.0*BrSM_htoWW*BrSM_htotautau*Br_tautohadrons*Br_tautohadrons*(Br_Wtoenu+Br_Wtomunu+Br_Wtotaunu*Br_tautohadrons)*(Br_Wtoenu+Br_Wtomunu+Br_Wtotaunu*Br_tautohadrons)
            +BrSM_htotautau*BrSM_htotautau*Br_tautohadrons*Br_tautohadrons*Br_tautohadrons*Br_tautohadrons);

            
    
        
    //pp_phi2_phi1phi1_bbWW_qqlnu_TH13=SigmaSumphi2_13*Br_phi2tophi1phi1*GTHDM_BR_h_bb*2.0*GTHDM_BR_h_WW*(Br_Wtoenu+Br_Wtomunu)*Br_Wtohadrons;
    //pp_phi3_phi1phi1_bbWW_qqlnu_TH13=SigmaSumphi3_13*Br_phi3tophi1phi1*GTHDM_BR_h_bb*2.0*GTHDM_BR_h_WW*(Br_Wtoenu+Br_Wtomunu)*Br_Wtohadrons;
    
    pp_phi2_phi1phi1_with_channel_bbWW_qqlnu_TH13=SigmaSumphi2_13*Br_phi2tophi1phi1*(GTHDM_BR_h_bb*GTHDM_BR_h_WW)/(BrSM_htobb*BrSM_htoWW);
    pp_phi3_phi1phi1_with_channel_bbWW_qqlnu_TH13=SigmaSumphi3_13*Br_phi3tophi1phi1*(GTHDM_BR_h_bb*GTHDM_BR_h_WW)/(BrSM_htobb*BrSM_htoWW);
    
    
    
    pp_phi2_phi1phi1_bbZZ_TH13=SigmaSumphi2_13*Br_phi2tophi1phi1*GTHDM_BR_h_bb*2.0*GTHDM_BR_h_ZZ;
    pp_phi3_phi1phi1_bbZZ_TH13=SigmaSumphi3_13*Br_phi3tophi1phi1*GTHDM_BR_h_bb*2.0*GTHDM_BR_h_ZZ;

    
    //pp_phi2_phi1phi1_bbZZ_lljj_TH13=SigmaSumphi2_13*Br_phi2tophi1phi1*GTHDM_BR_h_bb*2.0*GTHDM_BR_h_ZZ*(Br_Ztoee+Br_Ztomumu)*Br_Ztohadrons;
    //pp_phi3_phi1phi1_bbZZ_lljj_TH13=SigmaSumphi3_13*Br_phi3tophi1phi1*GTHDM_BR_h_bb*2.0*GTHDM_BR_h_ZZ*(Br_Ztoee+Br_Ztomumu)*Br_Ztohadrons;

    //pp_phi2_phi1phi1_bbZZ_llnunu_TH13=SigmaSumphi2_13*Br_phi2tophi1phi1*GTHDM_BR_h_bb*2.0*GTHDM_BR_h_ZZ*(Br_Ztoee+Br_Ztomumu)*Br_Ztoinv;
    //pp_phi3_phi1phi1_bbZZ_llnunu_TH13=SigmaSumphi3_13*Br_phi3tophi1phi1*GTHDM_BR_h_bb*2.0*GTHDM_BR_h_ZZ*(Br_Ztoee+Br_Ztomumu)*Br_Ztoinv;

    
    
    
    //pp_phi2_phi1phi1_bbWWorbbtautau_TH13=SigmaSumphi2_13*Br_phi2tophi1phi1*
    //        2.0*GTHDM_BR_h_bb*
    //        (GTHDM_BR_h_WW*((Br_Wtoenu+Br_Wtomunu)*(Br_Wtoenu+Br_Wtomunu)
    //        +2.0*(Br_Wtoenu+Br_Wtomunu)*Br_Wtohadrons)
    //        +GTHDM_BR_h_tautau*Br_tautoleptons*Br_tautoleptons);
    //pp_phi3_phi1phi1_bbWWorbbtautau_TH13=SigmaSumphi3_13*Br_phi3tophi1phi1*
    //        2.0*GTHDM_BR_h_bb*
    //        (GTHDM_BR_h_WW*((Br_Wtoenu+Br_Wtomunu)*(Br_Wtoenu+Br_Wtomunu)
    //        +2.0*(Br_Wtoenu+Br_Wtomunu)*Br_Wtohadrons)
    //        +GTHDM_BR_h_tautau*Br_tautoleptons*Br_tautoleptons);
    
    
    pp_phi2_phi1phi1_with_channel_bbWWorbbtautau_TH13=SigmaSumphi2_13*Br_phi2tophi1phi1*
            2.0*GTHDM_BR_h_bb*
            (GTHDM_BR_h_WW*((Br_Wtoenu+Br_Wtomunu)*(Br_Wtoenu+Br_Wtomunu)
            +2.0*(Br_Wtoenu+Br_Wtomunu)*Br_Wtohadrons)
            +GTHDM_BR_h_tautau*Br_tautoleptons*Br_tautoleptons)/(
            2.0*BrSM_htobb*
            (BrSM_htoWW*((Br_Wtoenu+Br_Wtomunu)*(Br_Wtoenu+Br_Wtomunu)
            +2.0*(Br_Wtoenu+Br_Wtomunu)*Br_Wtohadrons)
            +BrSM_htotautau*Br_tautoleptons*Br_tautoleptons)
            );
    pp_phi3_phi1phi1_with_channel_bbWWorbbtautau_TH13=SigmaSumphi3_13*Br_phi3tophi1phi1*
            2.0*GTHDM_BR_h_bb*
            (GTHDM_BR_h_WW*((Br_Wtoenu+Br_Wtomunu)*(Br_Wtoenu+Br_Wtomunu)
            +2.0*(Br_Wtoenu+Br_Wtomunu)*Br_Wtohadrons)
            +GTHDM_BR_h_tautau*Br_tautoleptons*Br_tautoleptons)/(
            2.0*BrSM_htobb*
            (BrSM_htoWW*((Br_Wtoenu+Br_Wtomunu)*(Br_Wtoenu+Br_Wtomunu)
            +2.0*(Br_Wtoenu+Br_Wtomunu)*Br_Wtohadrons)
            +BrSM_htotautau*Br_tautoleptons*Br_tautoleptons)
            );
    

    
    
    
    //pp_phi2_phi1phi1_bbWW_TH13=SigmaSumphi2_13*Br_phi2tophi1phi1*2.0*5.77e-1*2.15e-1; //SM Br of hh assumed
    //pp_phi3_phi1phi1_bbWW_TH13=SigmaSumphi3_13*Br_phi3tophi1phi1*2.0*5.77e-1*2.15e-1; //SM Br of hh assumed
    
    pp_phi2_phi1phi1_with_channel_bbWW_TH13=SigmaSumphi2_13*Br_phi2tophi1phi1*GTHDM_BR_h_bb*GTHDM_BR_h_WW/(BrSM_htobb*BrSM_htoWW); //SM Br of hh assumed
    pp_phi3_phi1phi1_with_channel_bbWW_TH13=SigmaSumphi3_13*Br_phi3tophi1phi1*GTHDM_BR_h_bb*GTHDM_BR_h_WW/(BrSM_htobb*BrSM_htoWW); //SM Br of hh assumed
    
    
    
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
    gg_phi2_phi1Z_with_channel_bbZ_TH13=SigmaggF_phi2_13*Br_phi2tophi1Z*GTHDM_BR_h_bb/BrSM_htobb;
    gg_phi3_phi1Z_with_channel_bbZ_TH13=SigmaggF_phi3_13*Br_phi3tophi1Z*GTHDM_BR_h_bb/BrSM_htobb;
    
    bb_phi2_phi1Z_bbZ_TH13=SigmabbF_phi2_13*Br_phi2tophi1Z*GTHDM_BR_h_bb;
    bb_phi3_phi1Z_bbZ_TH13=SigmabbF_phi3_13*Br_phi3tophi1Z*GTHDM_BR_h_bb;
    bb_phi2_phi1Z_with_channel_bbZ_TH13=SigmabbF_phi2_13*Br_phi2tophi1Z*GTHDM_BR_h_bb/BrSM_htobb;
    bb_phi3_phi1Z_with_channel_bbZ_TH13=SigmabbF_phi3_13*Br_phi3tophi1Z*GTHDM_BR_h_bb/BrSM_htobb;
    
    //gg_phi2_phi1Z_tautaull_TH13=SigmaggF_phi2_13*Br_phi2tophi1Z*GTHDM_BR_h_tautau*(Br_Ztoee+Br_Ztomumu)*Br_tautohadrons*Br_tautohadrons;
    //gg_phi3_phi1Z_tautaull_TH13=SigmaggF_phi3_13*Br_phi3tophi1Z*GTHDM_BR_h_tautau*(Br_Ztoee+Br_Ztomumu)*Br_tautohadrons*Br_tautohadrons;
    gg_phi2_phi1Z_tautaull_TH13=SigmaggF_phi2_13*Br_phi2tophi1Z*GTHDM_BR_h_tautau*(Br_Ztoee+Br_Ztomumu);
    gg_phi3_phi1Z_tautaull_TH13=SigmaggF_phi3_13*Br_phi3tophi1Z*GTHDM_BR_h_tautau*(Br_Ztoee+Br_Ztomumu);

    
    //std::cout<<"\033[1;34m gg_phi2_phi1Z_bbZ_TH8 = \033[0m "<<gg_phi2_phi1Z_bbZ_TH8<<std::endl;
    //std::cout<<"\033[1;34m gg_phi3_phi1Z_bbZ_TH8 = \033[0m "<<gg_phi3_phi1Z_bbZ_TH8<<std::endl;
    //std::cout<<"\033[1;34m gg_phi2_phi1Z_bbll_TH8 = \033[0m "<<gg_phi2_phi1Z_bbll_TH8<<std::endl;
    //std::cout<<"\033[1;34m gg_phi3_phi1Z_bbll_TH8 = \033[0m "<<gg_phi3_phi1Z_bbll_TH8<<std::endl;
    //std::cout<<"\033[1;34m gg_phi2_phi1Z_tautauZ_TH8 = \033[0m "<<gg_phi2_phi1Z_tautauZ_TH8<<std::endl;
    //std::cout<<"\033[1;34m gg_phi3_phi1Z_tautauZ_TH8 = \033[0m "<<gg_phi3_phi1Z_tautauZ_TH8<<std::endl;
    //std::cout<<"\033[1;34m gg_phi2_phi1Z_tautaull_TH8 = \033[0m "<<gg_phi2_phi1Z_tautaull_TH8<<std::endl;
    //std::cout<<"\033[1;34m gg_phi3_phi1Z_tautaull_TH8 = \033[0m "<<gg_phi3_phi1Z_tautaull_TH8<<std::endl;
    //std::cout<<"\033[1;34m gg_phi2_phi1Z_bbZ_TH13 = \033[0m "<<gg_phi2_phi1Z_bbZ_TH13<<std::endl;
    //std::cout<<"\033[1;34m gg_phi3_phi1Z_bbZ_TH13 = \033[0m "<<gg_phi3_phi1Z_bbZ_TH13<<std::endl;
    //std::cout<<"\033[1;34m bb_phi2_phi1Z_bbZ_TH13 = \033[0m "<<bb_phi2_phi1Z_bbZ_TH13<<std::endl;
    //std::cout<<"\033[1;34m bb_phi3_phi1Z_bbZ_TH13 = \033[0m "<<bb_phi3_phi1Z_bbZ_TH13<<std::endl;
    //std::cout<<"\033[1;31m Br_phi2tophi1Z = \033[0m "<<Br_phi2tophi1Z<<"\033[1;31m   THIS MUST BE SOLVED \033[0m "<<std::endl;
    //std::cout<< "\n" <<std::endl;
    
    
    
    
    pp_phi3_phi2Z_bbll_TH8=SigmaSumphi3_8*Br_phi3tophi2Z*Br_phi2tobb*(Br_Ztoee+Br_Ztomumu);
    pp_phi2_phi3Z_bbll_TH8=SigmaSumphi2_8*Br_phi2tophi3Z*Br_phi3tobb*(Br_Ztoee+Br_Ztomumu);     
    pp_phi3_phi2Z_tautaull_TH8=SigmaSumphi3_8*Br_phi3tophi2Z*Br_phi2totautau*(Br_Ztoee+Br_Ztomumu);
    pp_phi2_phi3Z_tautaull_TH8=SigmaSumphi2_8*Br_phi2tophi3Z*Br_phi3totautau*(Br_Ztoee+Br_Ztomumu);
    gg_phi3_phi2Z_bbZ_TH13=SigmaggF_phi3_13*Br_phi3tophi2Z*Br_phi2tobb;
    gg_phi2_phi3Z_bbZ_TH13=SigmaggF_phi2_13*Br_phi2tophi3Z*Br_phi3tobb;
    bb_phi3_phi2Z_bbZ_TH13=SigmabbF_phi3_13*Br_phi3tophi2Z*Br_phi2tobb;
    bb_phi2_phi3Z_bbZ_TH13=SigmabbF_phi2_13*Br_phi2tophi3Z*Br_phi3tobb;
    
    gg_phi3_phi2Z_WWZ_TH13=SigmaggF_phi3_13*Br_phi3tophi2Z*Br_phi2toWW;
    gg_phi2_phi3Z_WWZ_TH13=SigmaggF_phi2_13*Br_phi2tophi3Z*Br_phi3toWW;
    
    //std::cout<<"\033[1;34m SigmaggF_phi2_13*Br_phi2tophi3Z = \033[0m "<<SigmaggF_phi2_13*Br_phi2tophi3Z<<std::endl;
    //std::cout<<"\033[1;34m Br_phi3toWW = \033[0m "<<Br_phi3toWW<<std::endl;
    
    //std::cout<<"\033[1;34m pp_phi3_phi2Z_bbll_TH8 = \033[0m "<<pp_phi3_phi2Z_bbll_TH8<<std::endl;
    //std::cout<<"\033[1;34m pp_phi2_phi3Z_bbll_TH8 = \033[0m "<<pp_phi2_phi3Z_bbll_TH8<<std::endl;
    //std::cout<<"\033[1;34m pp_phi3_phi2Z_tautaull_TH8 = \033[0m "<<pp_phi3_phi2Z_tautaull_TH8<<std::endl;
    //std::cout<<"\033[1;34m pp_phi2_phi3Z_tautaull_TH8 = \033[0m "<<pp_phi2_phi3Z_tautaull_TH8<<std::endl;
    //std::cout<<"\033[1;34m gg_phi3_phi2Z_bbZ_TH13 = \033[0m "<<gg_phi3_phi2Z_bbZ_TH13<<std::endl;
    //std::cout<<"\033[1;34m gg_phi2_phi3Z_bbZ_TH13 = \033[0m "<<gg_phi2_phi3Z_bbZ_TH13<<std::endl;
    //std::cout<<"\033[1;34m bb_phi3_phi2Z_bbZ_TH13 = \033[0m "<<bb_phi3_phi2Z_bbZ_TH13<<std::endl;
    //std::cout<<"\033[1;34m bb_phi2_phi3Z_bbZ_TH13 = \033[0m "<<bb_phi2_phi3Z_bbZ_TH13<<std::endl;
    //std::cout<<"\033[1;31m Br_phi3tophi2Z = \033[0m "<<Br_phi3tophi2Z<<"\033[1;31m   THIS MUST BE SOLVED \033[0m "<<std::endl;
    //std::cout<< "\n" <<std::endl;
    
    
    //std::cout<<"\033[1;34m SigmaHp8 = \033[0m "<<SigmaHp8<<std::endl;
    //std::cout<<"\033[1;34m Br_Hptotaunu = \033[0m "<<Br_Hptotaunu<<std::endl;
        
    pp_Hpm_taunu_TH8=2.0*SigmaHp8*Br_Hptotaunu;
    pp_Hp_taunu_TH8=SigmaHp8*Br_Hptotaunu;
    pp_Hpm_taunu_TH13=2.0*SigmaHpm13*Br_Hptotaunu;
    pp_Hpm_tb_TH8=2.0*SigmaHp8*Br_Hptotb;
    pp_Hp_tb_TH8=SigmaHp8*Br_Hptotb;
    pp_Hpm_tb_TH13=2.0*SigmaHpm13*Br_Hptotb;
    
    //std::cout<<"\033[1;34m SigmaHp8 = \033[0m "<<SigmaHp8<<std::endl;
    //std::cout<<"\033[1;34m pp_Hpm_tb_TH13 = \033[0m "<<pp_Hpm_tb_TH13<<std::endl;
    
    //std::cout<<"mHp = "<<mHp<<std::endl;
    //std::cout<<"ip_ex_pp_Hpm_taunu_CMS13(100) = "<<ip_ex_pp_Hpm_taunu_CMS13(100)<<std::endl;
    //std::cout<<"ip_ex_pp_Hpm_taunu_CMS13(200) = "<<ip_ex_pp_Hpm_taunu_CMS13(200)<<std::endl;
    //std::cout<<"ip_ex_pp_Hpm_taunu_CMS13(300) = "<<ip_ex_pp_Hpm_taunu_CMS13(300)<<std::endl;
    //std::cout<<"ip_ex_pp_Hpm_taunu_CMS13(400) = "<<ip_ex_pp_Hpm_taunu_CMS13(400)<<std::endl;
    //std::cout<<"pp_Hpm_taunu_CMS13 = "<<ip_ex_pp_Hpm_taunu_CMS13(500)<<std::endl;
    //std::cout<<"pp_Hpm_taunu_CMS13 = "<<ip_ex_pp_Hpm_taunu_CMS13(1065.95)<<std::endl;
    
    //std::cout<<"ip_ex_pp_Hpm_taunu_CMS13 = "<<ip_ex_pp_Hpm_taunu_CMS13(mHp)<<std::endl;
    
    //std::cout<<"\033[1;34m pp_Hpm_taunu_TH8 = \033[0m "<<pp_Hpm_taunu_TH8<<std::endl;
    //std::cout<<"\033[1;34m pp_Hpm_taunu_TH8 = \033[0m "<<pp_Hpm_taunu_TH8<<std::endl;
    //std::cout<<"\033[1;34m pp_Hp_taunu_TH8 = \033[0m "<<pp_Hp_taunu_TH8<<std::endl;
    //std::cout<<"\033[1;34m pp_Hpm_taunu_TH13 = \033[0m "<<pp_Hpm_taunu_TH13<<std::endl;
    //std::cout<<"\033[1;34m pp_Hpm_tb_TH8 = \033[0m "<<pp_Hpm_tb_TH8<<std::endl;
    //std::cout<<"\033[1;34m pp_Hp_tb_TH8 = \033[0m "<<pp_Hp_tb_TH8<<std::endl;
    //std::cout<<"\033[1;34m pp_Hpm_tb_TH13 = \033[0m "<<pp_Hpm_tb_TH13<<std::endl;

    //std::cout<<"\033[1;31m SigmaHp8 = \033[0m "<<SigmaHp8<<"\033[1;31m   THIS MUST BE SOLVED \033[0m "<<std::endl;


    //95% to 1 sigma conversion factor, roughly sqrt(3.84)
//    double nftos=1.95996398454;

    
    //std::cout<<"\033[1;33m   m2 = \033[0m "<<m2<<std::endl;
    //std::cout<<"\033[1;33m   m3 = \033[0m "<<m3<<std::endl;
    //std::cout<<"\033[1;33m   mHp = \033[0m "<<mHp<<std::endl;
    
    //std::cout<<"\033[1;31m   mHp = \033[0m "<< mHp <<std::endl;
        
    if(m2>= 450.0 && m2<1400.0) {
                
        //std::cout<<"\033[1;31m   stop1 \033[0m "<<std::endl;
        
        THoEX_bb_phi2_bb_ATLAS13=bb_phi2_bb_TH13/ip_ex_bb_phi_bb_ATLAS13(m2);
    //    if(THoEX_bb_phi2_bb_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
    }                
    if(m3>= 450.0 && m3<1400.0){ 
                        
        //std::cout<<"\033[1;31m   stop2 \033[0m "<<std::endl;
        
        THoEX_bb_phi3_bb_ATLAS13=bb_phi3_bb_TH13/ip_ex_bb_phi_bb_ATLAS13(m3);
    //    if(THoEX_bb_phi3_bb_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
        }
    
    
    
    if(m2>= 400.0 && m2<1000.0) {
                        
        //std::cout<<"\033[1;31m   stop3 \033[0m "<<std::endl;
        
        THoEX_tt_phi2_tt_ATLAS13=tt_phi2_tt_TH13/ip_ex_tt_phi_tt_ATLAS13(m2);
    //    if(THoEX_tt_phi2_tt_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
    }                
    if(m3>= 400.0 && m3<1000.0){ 
                        
        //std::cout<<"\033[1;31m   stop4 \033[0m "<<std::endl;
        
        THoEX_tt_phi3_tt_ATLAS13=tt_phi3_tt_TH13/ip_ex_tt_phi_tt_ATLAS13(m3);
    //    if(THoEX_tt_phi3_tt_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
        }
        
    if(m2>= 400.0 && m2<1000.0) {
                        
        //std::cout<<"\033[1;31m   stop5 \033[0m "<<std::endl;
        
        THoEX_bb_phi2_tt_ATLAS13=bb_phi2_tt_TH13/ip_ex_bb_phi_tt_ATLAS13(m2);  
    //    if(THoEX_bb_phi2_tt_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
        }
        
    if(m3>= 400.0 && m3<1000.0){
                        
        //std::cout<<"\033[1;31m   stop6 \033[0m "<<std::endl;
        
        THoEX_bb_phi3_tt_ATLAS13=bb_phi3_tt_TH13/ip_ex_bb_phi_tt_ATLAS13(m3);  
    //    if(THoEX_bb_phi3_tt_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();

    }
        
    if(m2>= 100.0 && m2< 900.0) {
                        
        //std::cout<<"\033[1;31m   stop7 \033[0m "<<std::endl;
        
        THoEX_bb_phi2_bb_CMS8=bb_phi2_bb_TH8/ip_ex_bb_phi_bb_CMS8(m2);   
    //    if(THoEX_bb_phi2_bb_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();

    }
    if(m3>= 100.0 && m3< 900.0) {
                        
        //std::cout<<"\033[1;31m   stop8 \033[0m "<<std::endl;
        
        THoEX_bb_phi3_bb_CMS8=bb_phi3_bb_TH8/ip_ex_bb_phi_bb_CMS8(m3);    
    //    if(THoEX_bb_phi3_bb_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();
    }
    if(m2>= 330.0 && m2<1200.0) {
                        
        //std::cout<<"\033[1;31m   stop9 \033[0m "<<std::endl;
        
        THoEX_gg_phi2_bb_CMS8=gg_phi2_bb_TH8/ip_ex_gg_phi_bb_CMS8(m2);
    //    if(THoEX_gg_phi2_bb_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();
    }
    if(m3>= 330.0 && m3<1200.0) {
                        
        //std::cout<<"\033[1;31m   stop10 \033[0m "<<std::endl;
        
        THoEX_gg_phi3_bb_CMS8=gg_phi3_bb_TH8/ip_ex_gg_phi_bb_CMS8(m3);
    //    if(THoEX_gg_phi3_bb_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();
    }
    
    
    
    
    
    
    
    if(m2>= 350.0 && m2<650.0) {
                                
        //std::cout<<"\033[1;31m   stop11 \033[0m "<<std::endl;
        
        THoEX_tt_phi2_tt_CMS13=tt_phi2_tt_TH13/ip_ex_tt_phi2_tt_CMS13(m2);

    }
    if(m3>= 350.0 && m3<650.0) {
                                        
        //std::cout<<"\033[1;31m   stop12 \033[0m "<<std::endl;
        
        THoEX_tt_phi3_tt_CMS13=tt_phi3_tt_TH13/ip_ex_tt_phi3_tt_CMS13(m3); 
    //    if(THoEX_tt_phi3_tt_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
    }
    
    
    
    
    
    
    
    
    if(m2>= 550.0 && m2<1200.0) {
                                        
        //std::cout<<"\033[1;31m   stop13 \033[0m "<<std::endl;
        
        THoEX_pp_phi2_bb_CMS13=pp_phi2_bb_TH13/ip_ex_pp_phi_bb_CMS13(m2);
    //    if(THoEX_pp_phi2_bb_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
    }
    if(m3>= 550.0 && m3<1200.0) {
                                        
        //std::cout<<"\033[1;31m   stop14 \033[0m "<<std::endl;
        
        THoEX_pp_phi3_bb_CMS13=pp_phi3_bb_TH13/ip_ex_pp_phi_bb_CMS13(m3); 
    //    if(THoEX_pp_phi3_bb_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
    }
    
    
    
    
    
    
    
    if(m2>= 50.0 && m2<350.0) {
                                        
        //std::cout<<"\033[1;31m   stop15 \033[0m "<<std::endl;
        
        THoEX_pp_phi2_bb_light_CMS13=pp_phi2_bb_TH13/ip_ex_pp_phi2_bb_light_CMS13(m2);
    //    if(THoEX_pp_phi2_bb_light_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
    }
    if(m3>= 50.0 && m3<350.0) {
                                        
        //std::cout<<"\033[1;31m   stop16 \033[0m "<<std::endl;
        
        THoEX_pp_phi3_bb_light_CMS13=pp_phi3_bb_TH13/ip_ex_pp_phi3_bb_light_CMS13(m3); 
    //    if(THoEX_pp_phi3_bb_light_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
    }
    
    
    
    
    
    
    
    
    
    if(m2>= 300.0 && m2<1300.0) {
                                        
        //std::cout<<"\033[1;31m   stop17 \033[0m "<<std::endl;
        
        THoEX_bb_phi2_bb_CMS13=bb_phi2_bb_TH13/ip_ex_bb_phi_bb_CMS13(m2);  
    //    if(THoEX_bb_phi2_bb_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();      
    }
    if(m3>= 300.0 && m3<1300.0) {
                                        
        //std::cout<<"\033[1;31m   stop18 \033[0m "<<std::endl;
        
       THoEX_bb_phi3_bb_CMS13=bb_phi3_bb_TH13/ip_ex_bb_phi_bb_CMS13(m3);  
    //    if(THoEX_bb_phi3_bb_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();      
    }
        
    
    
    
    
    
    
    
     
    //1508.01437
    if(m2>=  120.0 && m2<500.0) {
                                        
        //std::cout<<"\033[1;31m   stop19 \033[0m "<<std::endl;
        
        THoEX_gg_phi2_mumu_CMS8=gg_phi2_mumu_TH8/ip_ex_gg_phi_mumu_CMS8(m2);   
    //    if(THoEX_gg_phi2_mumu_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();     
    }
    if(m3>=  120.0 && m3<500.0) {
                                        
        //std::cout<<"\033[1;31m   stop20 \033[0m "<<std::endl;
        
        THoEX_gg_phi3_mumu_CMS8=gg_phi3_mumu_TH8/ip_ex_gg_phi_mumu_CMS8(m3);     
    //    if(THoEX_gg_phi3_mumu_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();  
    }
    

    if(m2>=  120.0 && m2<500.0) {
                                                
        //std::cout<<"\033[1;31m   stop21 \033[0m "<<std::endl;
        
        THoEX_bb_phi2_mumu_CMS8=bb_phi2_mumu_TH8/ip_ex_bb_phi_mumu_CMS8(m2);   
    //    if(THoEX_bb_phi2_mumu_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();     
    }
    if(m3>=  120.0 && m3<500.0) {
                                                
        //std::cout<<"\033[1;31m   stop22 \033[0m "<<std::endl;
        
        THoEX_bb_phi3_mumu_CMS8=bb_phi3_mumu_TH8/ip_ex_bb_phi_mumu_CMS8(m3);     
    //    if(THoEX_bb_phi3_mumu_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();  
    }
    
    
     
     
     
     
     
     
     
     
     
     
    //1907.03152
    if(m2>=  140.0 && m2<1000.0) {
                                                
        //std::cout<<"\033[1;31m   stop23 \033[0m "<<std::endl;
        
        
        //std::cout<<"\033[1;32m   gg_phi2_mumu_TH13= \033[0m "<<  gg_phi2_mumu_TH13  <<std::endl;
        //std::cout<<"\033[1;32m   ip_ex_gg_phi_mumu_CMS13= \033[0m "<<  ip_ex_gg_phi_mumu_CMS13(m2)  <<std::endl;
        
        
        
        THoEX_gg_phi2_mumu_CMS13=gg_phi2_mumu_TH13/ip_ex_gg_phi_mumu_CMS13(m2);   
    //    if(THoEX_gg_phi2_mumu_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();     
    }
    if(m3>=  140.0 && m3<1000.0) {
                                                
        //std::cout<<"\033[1;31m   stop24 \033[0m "<<std::endl;
        
        THoEX_gg_phi3_mumu_CMS13=gg_phi3_mumu_TH13/ip_ex_gg_phi_mumu_CMS13(m3);     
    //    if(THoEX_gg_phi3_mumu_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();  
    }
    

    if(m2>=  140.0 && m2<1000.0) {
                                                
        //std::cout<<"\033[1;31m   stop25 \033[0m "<<std::endl;
        
        
        
        
        
        
        
        
        THoEX_bb_phi2_mumu_CMS13=bb_phi2_mumu_TH13/ip_ex_bb_phi_mumu_CMS13(m2);   
    //    if(THoEX_bb_phi2_mumu_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();     
    }
    if(m3>=  140.0 && m3<1000.0) {
                                                
        //std::cout<<"\033[1;31m   stop26 \033[0m "<<std::endl;
        
        
        
        
        
        
        THoEX_bb_phi3_mumu_CMS13=bb_phi3_mumu_TH13/ip_ex_bb_phi_mumu_CMS13(m3);     
    //    if(THoEX_bb_phi3_mumu_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();  
    }
    
    
    
    
    
    
    
    //1901.08144
    if(m2>=  200.0 && m2<1000.0) {
                                                
        //std::cout<<"\033[1;31m   stop27 \033[0m "<<std::endl;
        
        //std::cout<<"\033[1;32m   gg_phi2_mumu_TH13= \033[0m "<<  gg_phi2_mumu_TH13  <<std::endl;
        //std::cout<<"\033[1;32m   ip_ex_gg_phi_mumu_ATLAS13= \033[0m "<<  ip_ex_gg_phi_mumu_ATLAS13(m2)  <<std::endl;
        
        THoEX_gg_phi2_mumu_ATLAS13=gg_phi2_mumu_TH13/ip_ex_gg_phi_mumu_ATLAS13(m2);   
    //    if(THoEX_gg_phi2_mumu_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();     
    }
    if(m3>=  200.0 && m3<1000.0) {
                                                
        //std::cout<<"\033[1;31m   stop28 \033[0m "<<std::endl;
        
        THoEX_gg_phi3_mumu_ATLAS13=gg_phi3_mumu_TH13/ip_ex_gg_phi_mumu_ATLAS13(m3);     
    //    if(THoEX_gg_phi3_mumu_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();  
    }
    

    if(m2>=  200.0 && m2<1000.0) {
                                                
        //std::cout<<"\033[1;31m   stop29 \033[0m "<<std::endl;
        
        THoEX_bb_phi2_mumu_ATLAS13=bb_phi2_mumu_TH13/ip_ex_bb_phi_mumu_ATLAS13(m2);   
    //    if(THoEX_bb_phi2_mumu_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();     
    }
    if(m3>=  200.0 && m3<1000.0) {
                                                
        //std::cout<<"\033[1;31m   stop30 \033[0m "<<std::endl;
        
        THoEX_bb_phi3_mumu_ATLAS13=bb_phi3_mumu_TH13/ip_ex_bb_phi_mumu_ATLAS13(m3);     
    //    if(THoEX_bb_phi3_mumu_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();  
    }
    
    
    
    
    
    
    
    if(m2>=  90.0 && m2<1000.0) { //
                                                        
        //std::cout<<"\033[1;31m   gg_phi2_tautau_TH8 =\033[0m "<< gg_phi2_tautau_TH8 <<std::endl;
        //std::cout<<"\033[1;31m   ip_ex_gg_phi_tautau_ATLAS8(m2) =\033[0m "<< ip_ex_gg_phi_tautau_ATLAS8(m2) <<std::endl;
        
        THoEX_gg_phi2_tautau_ATLAS8=gg_phi2_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(m2);   
    //    if(THoEX_gg_phi2_tautau_ATLAS8 >5) return std::numeric_limits<double>::quiet_NaN();    
    }
    if(m3>=  90.0 && m3<1000.0) { //
                                                        
        //std::cout<<"\033[1;31m   stop32 \033[0m "<<std::endl;
        
        //std::cout<<"\033[1;31m   gg_phi3_tautau_TH8 = \033[0m "<< gg_phi3_tautau_TH8 <<std::endl;
        //std::cout<<"\033[1;31m   ip_ex_gg_phi_tautau_ATLAS8(m3) = \033[0m "<< ip_ex_gg_phi_tautau_ATLAS8(m3) <<std::endl;

        
        THoEX_gg_phi3_tautau_ATLAS8=gg_phi3_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(m3);   
    //    if(THoEX_gg_phi3_tautau_ATLAS8 >5) return std::numeric_limits<double>::quiet_NaN();    
    }
    if(m2>=  90.0 && m2<1000.0) {
                                                        
        //std::cout<<"\033[1;31m   stop33 \033[0m "<<std::endl;
        
        THoEX_gg_phi2_tautau_CMS8=gg_phi2_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(m2);   
    //    if(THoEX_gg_phi2_tautau_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();     
    }
    if(m3>=  90.0 && m3<1000.0) {
                                                        
        //std::cout<<"\033[1;31m   stop34 \033[0m "<<std::endl;
        
        THoEX_gg_phi3_tautau_CMS8=gg_phi3_tautau_TH8/ip_ex_gg_phi_tautau_CMS8(m3);     
    //    if(THoEX_gg_phi3_tautau_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();  
    }
    if(m2>=  90.0 && m2<1000.0) {
                                                        
        //std::cout<<"\033[1;31m   stop35 \033[0m "<<std::endl;
        
        THoEX_bb_phi2_tautau_ATLAS8=bb_phi2_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(m2);    
    //    if(THoEX_bb_phi2_tautau_ATLAS8 >5) return std::numeric_limits<double>::quiet_NaN();     
    }
    if(m3>=  90.0 && m3<1000.0) {
                                                        
        //std::cout<<"\033[1;31m   stop36 \033[0m "<<std::endl;
        
        THoEX_bb_phi3_tautau_ATLAS8=bb_phi3_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(m3);    
    //    if(THoEX_bb_phi3_tautau_ATLAS8 >5) return std::numeric_limits<double>::quiet_NaN();   
    }
    if(m2>=  90.0 && m2<1000.0) {
                                                        
        //std::cout<<"\033[1;31m   stop37 \033[0m "<<std::endl;
        
        THoEX_bb_phi2_tautau_CMS8=bb_phi2_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(m2);     
    //    if(THoEX_bb_phi2_tautau_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();   
    }
    if(m3>=  90.0 && m3<1000.0) {
                                                        
        //std::cout<<"\033[1;31m   stop38 \033[0m "<<std::endl;
        
        THoEX_bb_phi3_tautau_CMS8=bb_phi3_tautau_TH8/ip_ex_bb_phi_tautau_CMS8(m3);    
    //    if(THoEX_bb_phi3_tautau_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();    
    }
    if(m2>= 200.0 && m2<2250.0) {
                                                        
        //std::cout<<"\033[1;31m   stop39 \033[0m "<<std::endl;
        
        
        
        THoEX_gg_phi2_tautau_ATLAS13=gg_phi2_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(m2);  
    //    if(THoEX_gg_phi2_tautau_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();     
    }
        
        
    if(m3>= 200.0 && m3<2250.0) {
                                                        
        //std::cout<<"\033[1;31m   stop40 \033[0m "<<std::endl;
        
        //std::cout<<"\033[1;31m   ip_ex_gg_phi_tautau_ATLAS13(m2) \033[0m "<< ip_ex_gg_phi_tautau_ATLAS13(m3) <<std::endl;
        //std::cout<<"\033[1;31m   gg_phi3_tautau_TH13 \033[0m "<< gg_phi3_tautau_TH13 <<std::endl;
        
        THoEX_gg_phi3_tautau_ATLAS13=gg_phi3_tautau_TH13/ip_ex_gg_phi_tautau_ATLAS13(m3);     
    //    if(THoEX_gg_phi3_tautau_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
    }
    if(m2>=  90.0 && m2<3200.0) {
                                                                
        //std::cout<<"\033[1;31m   stop41 \033[0m "<<std::endl;
        
        //std::cout<<"\033[1;31m   ip_ex_gg_phi_tautau_CMS13(m2) \033[0m "<< ip_ex_gg_phi_tautau_CMS13(m2) <<std::endl;
        //std::cout<<"\033[1;31m   gg_phi2_tautau_TH13 \033[0m "<< gg_phi2_tautau_TH13 <<std::endl;
        
        THoEX_gg_phi2_tautau_CMS13=gg_phi2_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(m2);    
    //    if(THoEX_gg_phi2_tautau_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();    
    }
    if(m3>=  90.0 && m3<3200.0) {
                                                                
        //std::cout<<"\033[1;31m   stop42 \033[0m "<<std::endl;
        
        THoEX_gg_phi3_tautau_CMS13=gg_phi3_tautau_TH13/ip_ex_gg_phi_tautau_CMS13(m3);     
    //    if(THoEX_gg_phi3_tautau_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();  
    }
    if(m2>= 200.0 && m2<2250.0) {
                                                                
        //std::cout<<"\033[1;31m   stop43 \033[0m "<<std::endl;
        
        THoEX_bb_phi2_tautau_ATLAS13=bb_phi2_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(m2);   
    //    if(THoEX_bb_phi2_tautau_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();      
    }
        
        
    if(m3>= 200.0 && m3<2250.0) {
                                                                
        //std::cout<<"\033[1;31m   stop44 \033[0m "<<std::endl;
        
        THoEX_bb_phi3_tautau_ATLAS13=bb_phi3_tautau_TH13/ip_ex_bb_phi_tautau_ATLAS13(m3);    
    //    if(THoEX_bb_phi3_tautau_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();   
    }
    if(m2>=  60.0 && m2<3500.0) {
                                                                
        //std::cout<<"\033[1;31m   stop45 \033[0m "<<std::endl;
        
        THoEX_bb_phi2_tautau_CMS13=bb_phi2_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(m2);        
    //    if(THoEX_bb_phi2_tautau_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
    }
    if(m3>=  60.0 && m3<3500.0) {
                                                                
        //std::cout<<"\033[1;31m   stop46 \033[0m "<<std::endl;
        
        THoEX_bb_phi3_tautau_CMS13=bb_phi3_tautau_TH13/ip_ex_bb_phi_tautau_CMS13(m3);        
    //    if(THoEX_bb_phi3_tautau_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
    }
        
        
    if(m2>=  65.0 && m2< 600.0) {
                                                                
        //std::cout<<"\033[1;31m   stop47 \033[0m "<<std::endl;
        
        THoEX_gg_phi2_gaga_ATLAS8=gg_phi2_gaga_TH8/ip_ex_gg_phi_gaga_ATLAS8(m2);        
    //    if(THoEX_gg_phi2_gaga_ATLAS8 >5) return std::numeric_limits<double>::quiet_NaN();
    }
    if(m3>=  65.0 && m3< 600.0) {
                                                                
        //std::cout<<"\033[1;31m   stop48 \033[0m "<<std::endl;
        
        THoEX_gg_phi3_gaga_ATLAS8=gg_phi3_gaga_TH8/ip_ex_gg_phi_gaga_ATLAS8(m3);        
    //    if(THoEX_gg_phi3_gaga_ATLAS8 >5) return std::numeric_limits<double>::quiet_NaN();
    }
    if(m2>= 160.0 && m2<3000.0) {
                                                                
        //std::cout<<"\033[1;31m   stop49 \033[0m "<<std::endl;
        
        
        
        
        THoEX_pp_phi2_gaga_ATLAS13=pp_phi2_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(m2);     
    //    if(THoEX_pp_phi2_gaga_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();  
    }
        

    if(m3>= 160.0 && m3<3000.0)  
        {

        THoEX_pp_phi3_gaga_ATLAS13=pp_phi3_gaga_TH13/ip_ex_pp_phi_gaga_ATLAS13(m3);   

    }
    if(m2>= 600.0 && m2<5000.0)     
        {

        THoEX_gg_phi2_gaga_CMS13=gg_phi2_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(m2);        

    }
    if(m3>= 600.0 && m3<5000.0) 
        {

        THoEX_gg_phi3_gaga_CMS13=gg_phi3_gaga_TH13/ip_ex_gg_phi_gaga_CMS13(m3);     
   
    }
    if(m2>= 200.0 && m2<1600.0) 
        {

        THoEX_pp_phi2_Zga_llga_ATLAS8=pp_phi2_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(m2);       
    
    }
        
    if(m3>= 200.0 && m3<1600.0) 
        {
                                                                                
        //std::cout<<"\033[1;31m   stop54 \033[0m "<<std::endl;
        
        THoEX_pp_phi3_Zga_llga_ATLAS8=pp_phi3_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_ATLAS8(m3);     
    //    if(THoEX_pp_phi3_Zga_llga_ATLAS8 >5) return std::numeric_limits<double>::quiet_NaN();
    }
    if(m2>= 200.0 && m2<1200.0) 
        {
                                                                                
        //std::cout<<"\033[1;31m   stop55 \033[0m "<<std::endl;
        
        THoEX_pp_phi2_Zga_llga_CMS8=pp_phi2_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_CMS8(m2);     
    //    if(THoEX_pp_phi2_Zga_llga_CMS8 >5) return std::numeric_limits<double>::quiet_NaN(); 
    }
    if(m3>= 200.0 && m3<1200.0) 
        {
                                                                                
        //std::cout<<"\033[1;31m   stop56 \033[0m "<<std::endl;
        
        THoEX_pp_phi3_Zga_llga_CMS8=pp_phi3_Zga_llga_TH8/ip_ex_pp_phi_Zga_llga_CMS8(m3);     
    //    if(THoEX_pp_phi3_Zga_llga_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();  
    }
    if(m2>= 250.0 && m2<2400.0) 
        {
                                                                                
        //std::cout<<"\033[1;31m   stop57 \033[0m "<<std::endl;
        
        THoEX_gg_phi2_Zga_llga_ATLAS13=gg_phi2_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(m2);   
    //    if(THoEX_gg_phi2_Zga_llga_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();      
    }
    if(m3>= 250.0 && m3<2400.0) 
        {
                                                                                
        //std::cout<<"\033[1;31m   stop58 \033[0m "<<std::endl;
        
        THoEX_gg_phi3_Zga_llga_ATLAS13=gg_phi3_Zga_TH13/ip_ex_gg_phi_Zga_llga_ATLAS13(m3);   
    //    if(THoEX_gg_phi3_Zga_llga_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();     
    }
    if(m2>=1000.0 && m2<6800.0) 
        {
                                                                                
        //std::cout<<"\033[1;31m   stop59 \033[0m "<<std::endl;
        
        THoEX_gg_phi2_Zga_qqga_ATLAS13=gg_phi2_Zga_TH13/ip_ex_gg_phi_Zga_qqga_ATLAS13(m2);     
    //    if(THoEX_gg_phi2_Zga_qqga_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();  
    }
    if(m3>=1000.0 && m3<6800.0) 
        {
                                                                                
        //std::cout<<"\033[1;31m   stop60 \033[0m "<<std::endl;
        
        THoEX_gg_phi3_Zga_qqga_ATLAS13=gg_phi3_Zga_TH13/ip_ex_gg_phi_Zga_qqga_ATLAS13(m3);    
    //    if(THoEX_gg_phi3_Zga_qqga_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();    
    }
        
    if(m2>= 350.0 && m2<4000.0) 
        {
                                                                                        
        //std::cout<<"\033[1;31m   stop61 \033[0m "<<std::endl;
        
        THoEX_gg_phi2_Zga_CMS13=gg_phi2_Zga_TH13/ip_ex_gg_phi_Zga_CMS13(m2);      
    //    if(THoEX_gg_phi2_Zga_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();  
    }
    if(m3>= 350.0 && m3<4000.0) 
        {
                                                                                        
        //std::cout<<"\033[1;31m   stop62 \033[0m "<<std::endl;
        
        THoEX_gg_phi3_Zga_CMS13=gg_phi3_Zga_TH13/ip_ex_gg_phi_Zga_CMS13(m3);      
    //    if(THoEX_gg_phi3_Zga_CMS13>5) return std::numeric_limits<double>::quiet_NaN();  
    }
    if(m2>= 140.0 && m2<1000.0) 
        {
                                                                                        
        //std::cout<<"\033[1;31m   stop63 \033[0m "<<std::endl;
        
        THoEX_gg_phi2_ZZ_ATLAS8=gg_phi2_ZZ_TH8/ip_ex_gg_phi_ZZ_ATLAS8(m2);     
    //    if(THoEX_gg_phi2_ZZ_ATLAS8 >5) return std::numeric_limits<double>::quiet_NaN();    
    }
    if(m3>= 140.0 && m3<1000.0) 
        {
                                                                                        
        //std::cout<<"\033[1;31m   stop64 \033[0m "<<std::endl;
        
        THoEX_gg_phi3_ZZ_ATLAS8=gg_phi3_ZZ_TH8/ip_ex_gg_phi_ZZ_ATLAS8(m3);     
    //    if(THoEX_gg_phi3_ZZ_ATLAS8 >5) return std::numeric_limits<double>::quiet_NaN();    
    }
    if(m2>= 140.0 && m2<1000.0) 
        {
                                                                                        
        //std::cout<<"\033[1;31m   stop65 \033[0m "<<std::endl;
        
        THoEX_VV_phi2_ZZ_ATLAS8=VV_phi2_ZZ_TH8/ip_ex_VV_phi_ZZ_ATLAS8(m2);     
    //    if(THoEX_VV_phi2_ZZ_ATLAS8 >5) return std::numeric_limits<double>::quiet_NaN();  
    }     
        
    if(m3>= 140.0 && m3<1000.0) 
        {
                                                                                        
        //std::cout<<"\033[1;31m   stop66 \033[0m "<<std::endl;
        
        THoEX_VV_phi3_ZZ_ATLAS8=VV_phi3_ZZ_TH8/ip_ex_VV_phi_ZZ_ATLAS8(m3);     
    //    if(THoEX_VV_phi3_ZZ_ATLAS8 >5) return std::numeric_limits<double>::quiet_NaN();  
    }
    if(m2>= 210.0 && m2<1200.0) 
        {
                                                                                        
        //std::cout<<"\033[1;31m   stop67 \033[0m "<<std::endl;
        
        THoEX_gg_phi2_ZZ_llllnunu_ATLAS13=gg_phi2_ZZ_TH13/ip_ex_gg_phi_ZZ_llllnunu_ATLAS13(m2);   
    //    if(THoEX_gg_phi2_ZZ_llllnunu_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();     
    }
    if(m3>= 210.0 && m3<1200.0) 
        {
                                                                                        
        //std::cout<<"\033[1;31m   stop68 \033[0m "<<std::endl;
        
        THoEX_gg_phi3_ZZ_llllnunu_ATLAS13=gg_phi3_ZZ_TH13/ip_ex_gg_phi_ZZ_llllnunu_ATLAS13(m3);   
    //    if(THoEX_gg_phi3_ZZ_llllnunu_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();     
    }
    if(m2>= 210.0 && m2<2000.0) 
        {
                                                                                        
        //std::cout<<"\033[1;31m   stop69 \033[0m "<<std::endl;
        
        THoEX_VV_phi2_ZZ_llllnunu_ATLAS13=VV_phi2_ZZ_TH13/ip_ex_VV_phi_ZZ_llllnunu_ATLAS13(m2);     
    //    if(THoEX_VV_phi2_ZZ_llllnunu_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();   
    }
    if(m3>= 210.0 && m3<2000.0) 
        {
                                                                                        
        //std::cout<<"\033[1;31m   stop70 \033[0m "<<std::endl;
        
        THoEX_VV_phi3_ZZ_llllnunu_ATLAS13=VV_phi3_ZZ_TH13/ip_ex_VV_phi_ZZ_llllnunu_ATLAS13(m3);    
    //    if(THoEX_VV_phi3_ZZ_llllnunu_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();    
    }
    if(m2>= 300.0 && m2<3000.0) 
        {
                                                                                                
        //std::cout<<"\033[1;31m   stop71 \033[0m "<<std::endl;
        
        THoEX_gg_phi2_ZZ_qqllnunu_ATLAS13=gg_phi2_ZZ_TH13/ip_ex_gg_phi_ZZ_qqllnunu_ATLAS13(m2);    
    //    if(THoEX_gg_phi2_ZZ_qqllnunu_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();     
    }
    if(m3>= 300.0 && m3<3000.0) 
        {
                                                                                                
        //std::cout<<"\033[1;31m   stop72 \033[0m "<<std::endl;
        
        THoEX_gg_phi3_ZZ_qqllnunu_ATLAS13=gg_phi3_ZZ_TH13/ip_ex_gg_phi_ZZ_qqllnunu_ATLAS13(m3);    
    //    if(THoEX_gg_phi3_ZZ_qqllnunu_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();    
    }
    if(m2>= 300.0 && m2<3000.0) 
        {
                                                                                                
        //std::cout<<"\033[1;31m   stop73 \033[0m "<<std::endl;
        
        THoEX_VV_phi2_ZZ_qqllnunu_ATLAS13=VV_phi2_ZZ_TH13/ip_ex_VV_phi_ZZ_qqllnunu_ATLAS13(m2); 
    //    if(THoEX_VV_phi2_ZZ_qqllnunu_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
    }
    if(m3>= 300.0 && m3<3000.0) 
        {
                                                                                                
        //std::cout<<"\033[1;31m   stop74 \033[0m "<<std::endl;
        
        THoEX_VV_phi3_ZZ_qqllnunu_ATLAS13=VV_phi3_ZZ_TH13/ip_ex_VV_phi_ZZ_qqllnunu_ATLAS13(m3); 
    //    if(THoEX_VV_phi3_ZZ_qqllnunu_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
    }
    
    if(m2>= 130.0 && m2<3000.0)
         {
                                                                                                
        //std::cout<<"\033[1;31m   stop75 \033[0m "<<std::endl;
        
           THoEX_pp_phi2_ZZ_llqqnunull_CMS13=pp_phi2_ZZ_TH13/ip_ex_pp_phi_ZZ_llqqnunull_CMS13(m2);
    //       if(THoEX_pp_phi2_ZZ_llqqnunull_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
         }
    if(m3>= 130.0 && m3<3000.0)
          {
                                                                                                
        //std::cout<<"\033[1;31m   stop76 \033[0m "<<std::endl;
        
             THoEX_pp_phi3_ZZ_llqqnunull_CMS13=pp_phi3_ZZ_TH13/ip_ex_pp_phi_ZZ_llqqnunull_CMS13(m3);
    //         if(THoEX_pp_phi3_ZZ_llqqnunull_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m2>=1000.0 && m2<4000.0)
            {
                                                                                                
        //std::cout<<"\033[1;31m   stop77 \033[0m "<<std::endl;
        
             THoEX_pp_phi2_ZZ_qqnunu_CMS13=pp_phi2_ZZ_TH13/ip_ex_pp_phi_ZZ_qqnunu_CMS13(m2);
    //         if(THoEX_pp_phi2_ZZ_qqnunu_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m3>=1000.0 && m3<4000.0)
            {
                                                                                                
        //std::cout<<"\033[1;31m   stop78 \033[0m "<<std::endl;
        
             THoEX_pp_phi3_ZZ_qqnunu_CMS13=pp_phi3_ZZ_TH13/ip_ex_pp_phi_ZZ_qqnunu_CMS13(m3);
    //         if(THoEX_pp_phi3_ZZ_qqnunu_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m2>= 300.0 && m2<1500.0)
            {
                                                                                                
        //std::cout<<"\033[1;31m   stop79 \033[0m "<<std::endl;
        
             THoEX_gg_phi2_WW_ATLAS8=gg_phi2_WW_TH8/ip_ex_gg_phi_WW_ATLAS8(m2);
    //         if(THoEX_gg_phi2_WW_ATLAS8 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m3>= 300.0 && m3<1500.0)
            {
                                                                                                
        //std::cout<<"\033[1;31m   stop80 \033[0m "<<std::endl;
        
             THoEX_gg_phi3_WW_ATLAS8=gg_phi3_WW_TH8/ip_ex_gg_phi_WW_ATLAS8(m3);
    //         if(THoEX_gg_phi3_WW_ATLAS8 >5) return std::numeric_limits<double>::quiet_NaN();
             }  
    if(m2>= 300.0 && m2<1500.0)
            {
                                                                                                        
        //std::cout<<"\033[1;31m   stop81 \033[0m "<<std::endl;
        
             THoEX_VV_phi2_WW_ATLAS8=VV_phi2_WW_TH8/ip_ex_VV_phi_WW_ATLAS8(m2);
    //         if(THoEX_VV_phi2_WW_ATLAS8 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m3>= 300.0 && m3<1500.0)
            {
                                                                                                        
        //std::cout<<"\033[1;31m   stop82 \033[0m "<<std::endl;
        
             THoEX_VV_phi3_WW_ATLAS8=VV_phi3_WW_TH8/ip_ex_VV_phi_WW_ATLAS8(m3);
    //         if(THoEX_VV_phi3_WW_ATLAS8 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    
    //ip_ex_VV_phi_WW_CMS13
    
    if(m2>= 200.0 && m2<3000.0)
            {
                                                                                                        
        //std::cout<<"\033[1;31m   stop83 \033[0m "<<std::endl;
        
             THoEX_VV_phi2_WW_CMS13=VV_phi2_WW_TH13/ip_ex_VV_phi_WW_CMS13(m2);
    //         if(THoEX_VV_phi2_WW_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m3>= 200.0 && m3<3000.0)
            {
                                                                                                        
        //std::cout<<"\033[1;31m   stop84 \033[0m "<<std::endl;
        
             THoEX_VV_phi3_WW_CMS13=VV_phi3_WW_TH13/ip_ex_VV_phi_WW_CMS13(m3);
    //         if(THoEX_VV_phi3_WW_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m2>= 200.0 && m2<3000.0)
            {
                                                                                                        
        //std::cout<<"\033[1;31m   stop85 \033[0m "<<std::endl;
        
             THoEX_gg_phi2_WW_CMS13=gg_phi2_WW_TH13/ip_ex_gg_phi_WW_CMS13(m2);
    //         if(THoEX_gg_phi2_WW_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m3>= 200.0 && m3<3000.0)
            {
                                                                                                        
        //std::cout<<"\033[1;31m   stop86 \033[0m "<<std::endl;
        
             THoEX_gg_phi3_WW_CMS13=gg_phi3_WW_TH13/ip_ex_gg_phi_WW_CMS13(m3);
    //         if(THoEX_gg_phi3_WW_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    
    
    
    
    
    
    
    
    
    
    if(m2>= 1000.0 && m2<4500.0)
            {
                                                                                                        
        //std::cout<<"\033[1;31m   stop87 \033[0m "<<std::endl;
        
             THoEX_VV_phi2_WW_heavy_CMS13=VV_phi2_WW_TH13/ip_ex_VV_phi_WW_heavy_CMS13(m2);
    //         if(THoEX_VV_phi2_WW_heavy_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m3>= 1000.0 && m3<4500.0)
            {
                                                                                                        
        //std::cout<<"\033[1;31m   stop88 \033[0m "<<std::endl;
        
             THoEX_VV_phi3_WW_heavy_CMS13=VV_phi3_WW_TH13/ip_ex_VV_phi_WW_heavy_CMS13(m3);
    //         if(THoEX_VV_phi3_WW_heavy_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m2>= 1000.0 && m2<4500.0)
            {
                                                                                                        
        //std::cout<<"\033[1;31m   stop89 \033[0m "<<std::endl;
        
             THoEX_gg_phi2_WW_heavy_CMS13=gg_phi2_WW_TH13/ip_ex_gg_phi_WW_heavy_CMS13(m2);
    //         if(THoEX_gg_phi2_WW_heavy_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m3>= 1000.0 && m3<4500.0)
            {
                                                                                                        
        //std::cout<<"\033[1;31m   stop90 \033[0m "<<std::endl;
        
             THoEX_gg_phi3_WW_heavy_CMS13=gg_phi3_WW_TH13/ip_ex_gg_phi_WW_heavy_CMS13(m3);
    //         if(THoEX_gg_phi3_WW_heavy_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    
    
    
    
    
    
    
    
    
    
    
    
    if(m2>= 250.0 && m2<4000.0)
            {
                                                                                                                
        //std::cout<<"\033[1;31m   stop91 \033[0m "<<std::endl;
        
             THoEX_gg_phi2_WW_enumunu_ATLAS13=gg_phi2_WW_TH13/ip_ex_gg_phi_WW_enumunu_ATLAS13(m2);
    //         if(THoEX_gg_phi2_WW_enumunu_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m3>= 250.0 && m3<4000.0)
            {
                                                                                                                
        //std::cout<<"\033[1;31m   stop92 \033[0m "<<std::endl;
        
             THoEX_gg_phi3_WW_enumunu_ATLAS13=gg_phi3_WW_TH13/ip_ex_gg_phi_WW_enumunu_ATLAS13(m3);
    //         if(THoEX_gg_phi3_WW_enumunu_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m2>= 250.0 && m2<3000.0)
            {
                                                                                                                
        //std::cout<<"\033[1;31m   stop93 \033[0m "<<std::endl;
        
             THoEX_VV_phi2_WW_enumunu_ATLAS13=VV_phi2_WW_TH13/ip_ex_VV_phi_WW_enumunu_ATLAS13(m2);
    //         if(THoEX_VV_phi2_WW_enumunu_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m3>= 250.0 && m3<3000.0)
            {
                                                                                                                
        //std::cout<<"\033[1;31m   stop94 \033[0m "<<std::endl;
        
             THoEX_VV_phi3_WW_enumunu_ATLAS13=VV_phi3_WW_TH13/ip_ex_VV_phi_WW_enumunu_ATLAS13(m3);
    //         if(THoEX_VV_phi3_WW_enumunu_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m2>= 200.0 && m2<1000.0)
            {
                                                                                                                
        //std::cout<<"\033[1;31m   stop95 \033[0m "<<std::endl;
        
             THoEX_ggVV_phi2_WW_lnulnu_CMS13=ggVV_phi2_WW_lnulnu_TH13/ip_ex_ggVV_phi_WW_lnulnu_CMS13(m2);
    //         if(THoEX_ggVV_phi2_WW_lnulnu_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m3>= 200.0 && m3<1000.0)
            {
                                                                                                                
        //std::cout<<"\033[1;31m   stop96 \033[0m "<<std::endl;
        
             THoEX_ggVV_phi3_WW_lnulnu_CMS13=ggVV_phi3_WW_lnulnu_TH13/ip_ex_ggVV_phi_WW_lnulnu_CMS13(m3);
    //         if(THoEX_ggVV_phi3_WW_lnulnu_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m2>= 300.0 && m2<3000.0)
            {
                                                                                                                
        //std::cout<<"\033[1;31m   stop97 \033[0m "<<std::endl;
        
             THoEX_gg_phi2_WW_lnuqq_ATLAS13=gg_phi2_WW_TH13/ip_ex_gg_phi_WW_lnuqq_ATLAS13(m2);
    //         if(THoEX_gg_phi2_WW_lnuqq_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m3>= 300.0 && m3<3000.0)
            {
                                                                                                                
        //std::cout<<"\033[1;31m   stop98 \033[0m "<<std::endl;
        
             THoEX_gg_phi3_WW_lnuqq_ATLAS13=gg_phi3_WW_TH13/ip_ex_gg_phi_WW_lnuqq_ATLAS13(m3);
    //         if(THoEX_gg_phi3_WW_lnuqq_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m2>= 300.0 && m2<3000.0)
            {
                                                                                                                
        //std::cout<<"\033[1;31m   stop99 \033[0m "<<std::endl;
        
             THoEX_VV_phi2_WW_lnuqq_ATLAS13=VV_phi2_WW_TH13/ip_ex_VV_phi_WW_lnuqq_ATLAS13(m2);
    //         if(THoEX_VV_phi2_WW_lnuqq_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m3>= 300.0 && m3<3000.0)
            {
                                                                                                                
        //std::cout<<"\033[1;31m   stop100 \033[0m "<<std::endl;
        
             THoEX_VV_phi3_WW_lnuqq_ATLAS13=VV_phi3_WW_TH13/ip_ex_VV_phi_WW_lnuqq_ATLAS13(m3);
    //         if(THoEX_VV_phi3_WW_lnuqq_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m2>=1000.0 && m2<4400.0)
            {
                                                                                                                        
        //std::cout<<"\033[1;31m   stop101 \033[0m "<<std::endl;
        
             THoEX_pp_phi2_WW_lnuqq_CMS13=pp_phi2_WW_TH13/ip_ex_pp_phi_WW_lnuqq_CMS13(m2);
    //         if(THoEX_pp_phi2_WW_lnuqq_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m3>=1000.0 && m3<4400.0)
            {
                                                                                                                                
        //std::cout<<"\033[1;31m   stop102 \033[0m "<<std::endl;
        
             THoEX_pp_phi3_WW_lnuqq_CMS13=pp_phi3_WW_TH13/ip_ex_pp_phi_WW_lnuqq_CMS13(m3);
    //         if(THoEX_pp_phi3_WW_lnuqq_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m2>= 145.0 && m2<1000.0)
            {
                                                                                                                                
             THoEX_pp_phi2_VV_CMS8=pp_phi2_VV_TH8/ip_ex_pp_phi_VV_CMS8(m2);
    //         if(THoEX_mu_pp_phi2_VV_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m3>= 145.0 && m3<1000.0)
            {
                                                                                                                                
             THoEX_pp_phi3_VV_CMS8=pp_phi3_VV_TH8/ip_ex_pp_phi_VV_CMS8(m3);

             }
    if(m2>=1200.0 && m2<3000.0)
            {
                                                                                                                                

             THoEX_pp_phi2_VV_qqqq_ATLAS13=pp_phi2_VV_TH13/ip_ex_pp_phi_VV_qqqq_ATLAS13(m2);
    //         if(THoEX_pp_phi2_VV_qqqq_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m3>=1200.0 && m3<3000.0)
            {
                                                                                                                                
        //std::cout<<"\033[1;31m   stop106 \033[0m "<<std::endl;
        
             THoEX_pp_phi3_VV_qqqq_ATLAS13=pp_phi3_VV_TH13/ip_ex_pp_phi_VV_qqqq_ATLAS13(m3);
    //         if(THoEX_pp_phi3_VV_qqqq_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    
    
    
    
    
    if(m2>=300.0 && m2<5000.0)
            {
                                                                                                                                
        //std::cout<<"\033[1;31m   stop107 \033[0m "<<std::endl;
        
             THoEX_gg_phi2_VV_llqq_ATLAS13=gg_phi2_VV_TH13/ip_ex_gg_phi_VV_llqq_ATLAS13(m2);
    //         if(THoEX_gg_phi2_VV_llqq_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m3>=300.0 && m3<5000.0)
            {
                                                                                                                                
        //std::cout<<"\033[1;31m   stop108 \033[0m "<<std::endl;
        
             THoEX_gg_phi3_VV_llqq_ATLAS13=gg_phi3_VV_TH13/ip_ex_gg_phi_VV_llqq_ATLAS13(m3);
    //         if(THoEX_gg_phi3_VV_llqq_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    
    
    if(m2>=300.0 && m2<5000.0)
            {
                                                                                                                                
        //std::cout<<"\033[1;31m   stop109 \033[0m "<<std::endl;
        
             THoEX_VV_phi2_VV_llqq_ATLAS13=VV_phi2_VV_TH13/ip_ex_VV_phi_VV_llqq_ATLAS13(m2);
    //         if(THoEX_VV_phi2_VV_llqq_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m3>=300.0 && m3<5000.0)
            {
                                                                                                                                
        //std::cout<<"\033[1;31m   stop110 \033[0m "<<std::endl;
        
             THoEX_VV_phi3_VV_llqq_ATLAS13=VV_phi3_VV_TH13/ip_ex_VV_phi_VV_llqq_ATLAS13(m3);
    //         if(THoEX_VV_phi3_VV_llqq_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
            
            
            
            
    
    if(m2>= 260.0 && m2<1000.0)
            {
                                                                                                                                        
        
        THoEX_gg_phi2_phi1phi1_ATLAS8=gg_phi2_phi1phi1_TH8/ip_ex_gg_phi_phi1phi1_ATLAS8(m2);
    
    }   
    if(m3>= 260.0 && m3<1000.0)
            {
                                                                                                                                        
        //std::cout<<"\033[1;31m   stop112 \033[0m "<<std::endl;
        
             THoEX_gg_phi3_phi1phi1_ATLAS8=gg_phi3_phi1phi1_TH8/ip_ex_gg_phi_phi1phi1_ATLAS8(m3);
    //         if(THoEX_gg_phi3_phi1phi1_ATLAS8 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m2>= 270.0 && m2<1100.0)
            {
                                                                                                                                        
        //std::cout<<"\033[1;31m   stop113 \033[0m "<<std::endl;
        
             THoEX_pp_phi2_phi1phi1_bbbb_CMS8=pp_phi2_phi1phi1_bbbb_TH8/ip_ex_pp_phi_phi1phi1_bbbb_CMS8(m2);
    //         if(THoEX_pp_phi2_phi1phi1_bbbb_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m3>= 270.0 && m3<1100.0)
            {
                                                                                                                                        
        //std::cout<<"\033[1;31m   stop114 \033[0m "<<std::endl;
        
             THoEX_pp_phi3_phi1phi1_bbbb_CMS8=pp_phi3_phi1phi1_bbbb_TH8/ip_ex_pp_phi_phi1phi1_bbbb_CMS8(m3);
    //         if(THoEX_pp_phi3_phi1phi1_bbbb_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();
             }            
    if(m2>= 260.0 && m2<1100.0)
            {
                                                                                                                                        
        //std::cout<<"\033[1;31m   stop115 \033[0m "<<std::endl;
        
             THoEX_pp_phi2_phi1phi1_bbgaga_CMS8=pp_phi2_phi1phi1_bbgaga_TH8/ip_ex_pp_phi_phi1phi1_bbgaga_CMS8(m2);
    //         if(THoEX_pp_phi2_phi1phi1_bbgaga_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m3>= 260.0 && m3<1100.0)
            {
                                                                                                                                        
        //std::cout<<"\033[1;31m   stop116 \033[0m "<<std::endl;
        
             THoEX_pp_phi3_phi1phi1_bbgaga_CMS8=pp_phi3_phi1phi1_bbgaga_TH8/ip_ex_pp_phi_phi1phi1_bbgaga_CMS8(m3);
    //         if(THoEX_pp_phi3_phi1phi1_bbgaga_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m2>= 260.0 && m2< 350.0)
            {
                                                                                                                                        
        //std::cout<<"\033[1;31m   stop117 \033[0m "<<std::endl;
        
             THoEX_gg_phi2_phi1phi1_bbtautau_CMS8=gg_phi2_phi1phi1_bbtautau_TH8/ip_ex_gg_phi_phi1phi1_bbtautau_CMS8(m2);
    //         if(THoEX_gg_phi2_phi1phi1_bbtautau_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m3>= 260.0 && m3< 350.0)
            {
                                                                                                                                        
        //std::cout<<"\033[1;31m   stop118 \033[0m "<<std::endl;
        
             THoEX_gg_phi3_phi1phi1_bbtautau_CMS8=gg_phi3_phi1phi1_bbtautau_TH8/ip_ex_gg_phi_phi1phi1_bbtautau_CMS8(m3);
    //         if(THoEX_gg_phi3_phi1phi1_bbtautau_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m2>= 350.0 && m2<1000.0)
            {
                                                                                                                                        
        //std::cout<<"\033[1;31m   stop119 \033[0m "<<std::endl;
        
             THoEX_pp_phi2_phi1phi1_bbtautau_CMS8=pp_phi2_phi1phi1_TH8/ip_ex_pp_phi_phi1phi1_bbtautau_CMS8(m2);
    //         if(THoEX_pp_phi2_phi1phi1_bbtautau_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m3>= 350.0 && m3<1000.0)
            {
                                                                                                                                        
        //std::cout<<"\033[1;31m   stop120 \033[0m "<<std::endl;
        
             THoEX_pp_phi3_phi1phi1_bbtautau_CMS8=pp_phi3_phi1phi1_TH8/ip_ex_pp_phi_phi1phi1_bbtautau_CMS8(m3);
    //         if(THoEX_pp_phi3_phi1phi1_bbtautau_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m2>= 250.0 && m2<3000.0)
            {
                                                                                                                                                
        //std::cout<<"\033[1;31m   stop121 \033[0m "<<std::endl;
        
             THoEX_pp_phi2_phi1phi1_bbbb_ATLAS13=pp_phi2_phi1phi1_with_channel_bbbb_TH13/ip_ex_pp_phi_phi1phi1_bbbb_ATLAS13(m2);
    //         if(THoEX_pp_phi2_phi1phi1_bbbb_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m3>= 250.0 && m3<3000.0)
            {
                                                                                                                                                
        //std::cout<<"\033[1;31m   stop122 \033[0m "<<std::endl;
        
             THoEX_pp_phi3_phi1phi1_bbbb_ATLAS13=pp_phi3_phi1phi1_with_channel_bbbb_TH13/ip_ex_pp_phi_phi1phi1_bbbb_ATLAS13(m3);
    //         if(THoEX_pp_phi3_phi1phi1_bbbb_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m2>= 260.0 && m2<1200.0)
            {
                                                                                                                                                
        //std::cout<<"\033[1;31m   stop123 \033[0m "<<std::endl;
        
             THoEX_pp_phi2_phi1phi1_bbbb_1_CMS13=pp_phi2_phi1phi1_bbbb_TH13/ip_ex_pp_phi_phi1phi1_bbbb_1_CMS13(m2);
    //         if(THoEX_pp_phi2_phi1phi1_bbbb_1_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m3>= 260.0 && m3<1200.0)
            {
                                                                                                                                                
        //std::cout<<"\033[1;31m   stop124 \033[0m "<<std::endl;
        
             THoEX_pp_phi3_phi1phi1_bbbb_1_CMS13=pp_phi3_phi1phi1_bbbb_TH13/ip_ex_pp_phi_phi1phi1_bbbb_1_CMS13(m3);
    //         if(THoEX_pp_phi3_phi1phi1_bbbb_1_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m2>=1000.0 && m2<3000.0)
            {
                                                                                                                                                
        //std::cout<<"\033[1;31m   stop125 \033[0m "<<std::endl;
        
             THoEX_pp_phi2_phi1phi1_bbbb_2_CMS13=pp_phi2_phi1phi1_with_channel_bbbb_TH13/ip_ex_pp_phi_phi1phi1_bbbb_2_CMS13(m2);
    //         if(THoEX_pp_phi2_phi1phi1_bbbb_2_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m3>=1000.0 && m3<3000.0)
            {
                                                                                                                                                
        //std::cout<<"\033[1;31m   stop126 \033[0m "<<std::endl;
        
             THoEX_pp_phi3_phi1phi1_bbbb_2_CMS13=pp_phi3_phi1phi1_with_channel_bbbb_TH13/ip_ex_pp_phi_phi1phi1_bbbb_2_CMS13(m3);
    //         if(THoEX_pp_phi3_phi1phi1_bbbb_2_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m2>= 251.0 && m2<1000.0)
            {
                                                                                                                                                
        //std::cout<<"\033[1;31m   stop127 \033[0m "<<std::endl;
        
             THoEX_pp_phi2_phi1phi1_bbgaga_ATLAS13=pp_phi2_phi1phi1_with_channel_bbgaga_TH13/ip_ex_pp_phi_phi1phi1_bbgaga_ATLAS13(m2);
    //         if(THoEX_pp_phi2_phi1phi1_bbgaga_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    if(m3>= 251.0 && m3<1000.0)
            {
        
             THoEX_pp_phi3_phi1phi1_bbgaga_ATLAS13=pp_phi3_phi1phi1_with_channel_bbgaga_TH13/ip_ex_pp_phi_phi1phi1_bbgaga_ATLAS13(m3);
    //         if(THoEX_pp_phi3_phi1phi1_bbgaga_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
             }
    //if(m2>= 250.0 && m2< 900.0) //Old limits
    if(m2>= 260.0 && m2< 1000.0) 
        {
                                                                                                                    
        
        THoEX_pp_phi2_phi1phi1_bbgaga_CMS13=pp_phi2_phi1phi1_bbgaga_TH13/ip_ex_pp_phi_phi1phi1_bbgaga_CMS13(m2);
    
        }
    //if(m3>= 250.0 && m3< 900.0) //Old limits
    if(m3>= 260.0 && m3< 1000.0) 
        {

        THoEX_pp_phi3_phi1phi1_bbgaga_CMS13=pp_phi3_phi1phi1_bbgaga_TH13/ip_ex_pp_phi_phi1phi1_bbgaga_CMS13(m3);
    }
    /*
    if(m2>= 260.0 && m2<1000.0) 
        {
        THoEX_pp_phi2_phi1phi1_bbtautau_ATLAS13=pp_phi2_phi1phi1_bbtautau_TH13/ip_ex_pp_phi_phi1phi1_bbtautau_ATLAS13(m2);
             if(THoEX_pp_phi2_phi1phi1_bbtautau_ATLAS13 >2) return std::numeric_limits<double>::quiet_NaN();
        }
    if(m3>= 260.0 && m3<1000.0) 
        {
        THoEX_pp_phi3_phi1phi1_bbtautau_ATLAS13=pp_phi3_phi1phi1_bbtautau_TH13/ip_ex_pp_phi_phi1phi1_bbtautau_ATLAS13(m3);
             if(THoEX_pp_phi3_phi1phi1_bbtautau_ATLAS13 >2) return std::numeric_limits<double>::quiet_NaN();
        }
    */ //OLD We have split this one in two parts
    /////////
    if(m2>= 251.0 && m2<1600.0) 
        {
                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop131 \033[0m "<<std::endl;
        
        THoEX_pp_phi2_phi1phi1_bbtautau_1_ATLAS13=pp_phi2_phi1phi1_with_channel_bbtautau_TH13/ip_ex_pp_phi_phi1phi1_bbtautau_1_ATLAS13(m2);
    //         if(THoEX_pp_phi2_phi1phi1_bbtautau_1_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
        }
    if(m3>= 251.0 && m3<1600.0) 
        {
                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop132 \033[0m "<<std::endl;
        
        THoEX_pp_phi3_phi1phi1_bbtautau_1_ATLAS13=pp_phi3_phi1phi1_with_channel_bbtautau_TH13/ip_ex_pp_phi_phi1phi1_bbtautau_1_ATLAS13(m3);
    //         if(THoEX_pp_phi3_phi1phi1_bbtautau_1_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
        }
    if(m2>= 1000.0 && m2<3000.0) 
        {
                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop133 \033[0m "<<std::endl;
        
        THoEX_pp_phi2_phi1phi1_bbtautau_2_ATLAS13=pp_phi2_phi1phi1_with_channel_bbtautau_TH13/ip_ex_pp_phi_phi1phi1_bbtautau_2_ATLAS13(m2);
    //         if(THoEX_pp_phi2_phi1phi1_bbtautau_2_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
        }
    if(m3>= 1000.0 && m3<3000.0) 
        {
                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop134 \033[0m "<<std::endl;
        
        THoEX_pp_phi3_phi1phi1_bbtautau_2_ATLAS13=pp_phi3_phi1phi1_with_channel_bbtautau_TH13/ip_ex_pp_phi_phi1phi1_bbtautau_2_ATLAS13(m3);
    //         if(THoEX_pp_phi3_phi1phi1_bbtautau_2_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
        }
    /////////
    if(m2>= 250.0 && m2< 900.0) 
        {
                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop135 \033[0m "<<std::endl;
        
        THoEX_pp_phi2_phi1phi1_bbtautau_1_CMS13=pp_phi2_phi1phi1_bbtautau_TH13/ip_ex_pp_phi_phi1phi1_bbtautau_1_CMS13(m2);
    //         if(THoEX_pp_phi2_phi1phi1_bbtautau_1_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
        }
    if(m3>= 250.0 && m3< 900.0) 
        {
                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop136 \033[0m "<<std::endl;
        
        THoEX_pp_phi3_phi1phi1_bbtautau_1_CMS13=pp_phi3_phi1phi1_bbtautau_TH13/ip_ex_pp_phi_phi1phi1_bbtautau_1_CMS13(m3);
    //         if(THoEX_pp_phi3_phi1phi1_bbtautau_1_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
        }
    if(m2>= 900.0 && m2<4000.0) 
        {
                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop137 \033[0m "<<std::endl;
        
        THoEX_pp_phi2_phi1phi1_bbtautau_2_CMS13=pp_phi2_phi1phi1_with_channel_bbtautau_TH13/ip_ex_pp_phi_phi1phi1_bbtautau_2_CMS13(m2);
    //         if(THoEX_pp_phi2_phi1phi1_bbtautau_2_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
        }
    if(m3>= 900.0 && m3<4000.0) 
        {
                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop138 \033[0m "<<std::endl;
        
        THoEX_pp_phi3_phi1phi1_bbtautau_2_CMS13=pp_phi3_phi1phi1_with_channel_bbtautau_TH13/ip_ex_pp_phi_phi1phi1_bbtautau_2_CMS13(m3);
    //         if(THoEX_pp_phi3_phi1phi1_bbtautau_2_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
        }
    if(m2>= 260.0 && m2< 900.0) 
        {
                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop139 \033[0m "<<std::endl;
        
        THoEX_pp_phi2_phi1phi1_bbVV_CMS13=pp_phi2_phi1phi1_bbVV_TH13/ip_ex_pp_phi_phi1phi1_bbVV_CMS13(m2);
    //         if(THoEX_pp_phi2_phi1phi1_bbVV_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
        }
    if(m3>= 260.0 && m3< 900.0) 
        {
                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop140 \033[0m "<<std::endl;
        
        THoEX_pp_phi3_phi1phi1_bbVV_CMS13=pp_phi3_phi1phi1_bbVV_TH13/ip_ex_pp_phi_phi1phi1_bbVV_CMS13(m3);
    //         if(THoEX_pp_phi3_phi1phi1_bbVV_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
        } 
    
    
    if(m2>= 250.0 && m2< 1000.0) 
        {
                                                                                                                                                                
        //std::cout<<"\033[1;31m   stop141 \033[0m "<<std::endl;
        
        THoEX_pp_phi2_phi1phi1_4WOr2W2tauOr4tau_CMS13=pp_phi2_phi1phi1_with_channel_4WOr2W2tauOr4tau_TH13/ip_ex_pp_phi_phi1phi1_4WOr2W2tauOr4tau_CMS13(m2);
    //         if(THoEX_pp_phi2_phi1phi1_4WOr2W2tauOr4tau_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
        }
    if(m3>= 250.0 && m3< 1000.0) 
        {
                                                                                                                                                                
        //std::cout<<"\033[1;31m   stop142 \033[0m "<<std::endl;
        
        THoEX_pp_phi3_phi1phi1_4WOr2W2tauOr4tau_CMS13=pp_phi3_phi1phi1_with_channel_4WOr2W2tauOr4tau_TH13/ip_ex_pp_phi_phi1phi1_4WOr2W2tauOr4tau_CMS13(m3);
    //         if(THoEX_pp_phi3_phi1phi1_4WOr2W2tauOr4tau_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
        } 
    
    
    
    if(m2>= 800.0 && m2< 3500.0) 
        {
                                                                                                                                                                
        //std::cout<<"\033[1;31m   stop143 \033[0m "<<std::endl;
        
        THoEX_pp_phi2_phi1phi1_bbWW_qqlnu_CMS13=pp_phi2_phi1phi1_with_channel_bbWW_qqlnu_TH13/ip_ex_pp_phi_phi1phi1_bbWW_qqlnu_CMS13(m2);
    //         if(THoEX_pp_phi2_phi1phi1_bbWW_qqlnu_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
        }
    if(m3>= 800.0 && m3< 3500.0) 
        {
                                                                                                                                                                
        //std::cout<<"\033[1;31m   stop144 \033[0m "<<std::endl;
        
        THoEX_pp_phi3_phi1phi1_bbWW_qqlnu_CMS13=pp_phi3_phi1phi1_with_channel_bbWW_qqlnu_TH13/ip_ex_pp_phi_phi1phi1_bbWW_qqlnu_CMS13(m3);
    //         if(THoEX_pp_phi3_phi1phi1_bbWW_qqlnu_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
        } 
    
    
    
    
    
    
    if(m2>= 260.0 && m2< 1000.0) 
        {
                                                                                                                                                                
        //std::cout<<"\033[1;31m   stop145 \033[0m "<<std::endl;
        
        THoEX_pp_phi2_phi1phi1_bbZZ_lljj_CMS13=pp_phi2_phi1phi1_bbZZ_TH13/ip_ex_pp_phi_phi1phi1_bbZZ_lljj_CMS13(m2);
    //         if(THoEX_pp_phi2_phi1phi1_bbZZ_lljj_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
        }
    if(m3>= 260.0 && m3< 1000.0) 
        {
                                                                                                                                                                
        //std::cout<<"\033[1;31m   stop146 \033[0m "<<std::endl;
        
        THoEX_pp_phi3_phi1phi1_bbZZ_lljj_CMS13=pp_phi3_phi1phi1_bbZZ_TH13/ip_ex_pp_phi_phi1phi1_bbZZ_lljj_CMS13(m3);
    //         if(THoEX_pp_phi3_phi1phi1_bbZZ_lljj_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
        } 
    if(m2>= 250.0 && m2< 1000.0) 
        {
                                                                                                                                                                
        //std::cout<<"\033[1;31m   stop147 \033[0m "<<std::endl;
        
        THoEX_pp_phi2_phi1phi1_bbZZ_llnunu_CMS13=pp_phi2_phi1phi1_bbZZ_TH13/ip_ex_pp_phi_phi1phi1_bbZZ_llnunu_CMS13(m2);
    //         if(THoEX_pp_phi2_phi1phi1_bbZZ_llnunu_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
        }
    if(m3>= 250.0 && m3< 1000.0) 
        {
                                                                                                                                                                
        //std::cout<<"\033[1;31m   stop148 \033[0m "<<std::endl;
        
        THoEX_pp_phi3_phi1phi1_bbZZ_llnunu_CMS13=pp_phi3_phi1phi1_bbZZ_TH13/ip_ex_pp_phi_phi1phi1_bbZZ_llnunu_CMS13(m3);
    //         if(THoEX_pp_phi3_phi1phi1_bbZZ_llnunu_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
        } 
    
    
    //ip_ex_pp_phi_phi1phi1_bbWWorbbtautau_CMS13
      
    if(m2>= 800.0 && m2< 4500.0) 
        {
        
        THoEX_pp_phi2_phi1phi1_bbWWorbbtautau_CMS13=pp_phi2_phi1phi1_with_channel_bbWWorbbtautau_TH13/ip_ex_pp_phi_phi1phi1_bbWWorbbtautau_CMS13(m2);
    
    }
    if(m3>= 800.0 && m3< 4500.0) 
        {
        
        THoEX_pp_phi3_phi1phi1_bbWWorbbtautau_CMS13=pp_phi3_phi1phi1_with_channel_bbWWorbbtautau_TH13/ip_ex_pp_phi_phi1phi1_bbWWorbbtautau_CMS13(m3);
        } 

    
    //This will overwrite the previous one when the mass is in this range.
    //This way we consider this search only once so it's fine
    if(m2>= 250.0 && m2 < 899.0) 
        {
        
        THoEX_pp_phi2_phi1phi1_bbWWorbbtautau_CMS13=pp_phi2_phi1phi1_with_channel_bbWWorbbtautau_TH13/ip_ex_pp_phi_phi1phi1_bbWWorbbtautau_low_masses_CMS13(m2);
    
    }
    if(m3>= 250.0 && m3< 899.0) 
        {
        
        THoEX_pp_phi3_phi1phi1_bbWWorbbtautau_CMS13=pp_phi3_phi1phi1_with_channel_bbWWorbbtautau_TH13/ip_ex_pp_phi_phi1phi1_bbWWorbbtautau_low_masses_CMS13(m3);
        } 
    
    
    
    if(m2>= 500.0 && m2< 3000.0) 
        {
                                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop151 \033[0m "<<std::endl;
        
        THoEX_pp_phi2_phi1phi1_bbWW_ATLAS13=pp_phi2_phi1phi1_with_channel_bbWW_TH13/ip_ex_pp_phi_phi1phi1_bbWW_ATLAS13(m2);
    //         if(THoEX_pp_phi2_phi1phi1_bbWW_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
        }
    if(m3>= 500.0 && m3< 3000.0) 
        {
                                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop152 \033[0m "<<std::endl;
        
        THoEX_pp_phi3_phi1phi1_bbWW_ATLAS13=pp_phi3_phi1phi1_with_channel_bbWW_TH13/ip_ex_pp_phi_phi1phi1_bbWW_ATLAS13(m3);
    //         if(THoEX_pp_phi3_phi1phi1_bbWW_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
        }  
    if(m2>= 260.0 && m2< 500.0)  
        {
                                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop153 \033[0m "<<std::endl;
        
        THoEX_gg_phi2_phi1phi1_gagaWW_ATLAS13=gg_phi2_phi1phi1_gagaWW_TH13/ip_ex_gg_phi_phi1phi1_gagaWW_ATLAS13(m2);
    //     if(THoEX_gg_phi2_phi1phi1_gagaWW_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
        }
    if(m3>= 260.0 && m3< 500.0)  
        {
                                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop154 \033[0m "<<std::endl;
        
        THoEX_gg_phi3_phi1phi1_gagaWW_ATLAS13=gg_phi3_phi1phi1_gagaWW_TH13/ip_ex_gg_phi_phi1phi1_gagaWW_ATLAS13(m3);
    //     if(THoEX_gg_phi3_phi1phi1_gagaWW_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
        }
    if(m2>= 220.0 && m2<1000.0)
        {
                                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop155 \033[0m "<<std::endl;
        
        THoEX_gg_phi2_phi1Z_bbZ_ATLAS8=gg_phi2_phi1Z_bbZ_TH8/ip_ex_gg_phi_phi1Z_bbZ_ATLAS8(m2);
    //    if(THoEX_gg_phi2_phi1Z_bbZ_ATLAS8 >5) return std::numeric_limits<double>::quiet_NaN();
         }
    if(m3>= 220.0 && m3<1000.0)
        {
                                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop156 \033[0m "<<std::endl;
        
         THoEX_gg_phi3_phi1Z_bbZ_ATLAS8=gg_phi3_phi1Z_bbZ_TH8/ip_ex_gg_phi_phi1Z_bbZ_ATLAS8(m3);
    //     if(THoEX_gg_phi3_phi1Z_bbZ_ATLAS8 >5) return std::numeric_limits<double>::quiet_NaN();
         }
    if(m2>= 225.0 && m2< 600.0)
        {
                                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop157 \033[0m "<<std::endl;
        
         THoEX_gg_phi2_phi1Z_bbll_CMS8=gg_phi2_phi1Z_bbll_TH8/ip_ex_gg_phi_phi1Z_bbll_CMS8(m2);
    //    if(THoEX_gg_phi2_phi1Z_bbll_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();
         }
    if(m3>= 225.0 && m3< 600.0)
        {
                                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop158 \033[0m "<<std::endl;
        
         THoEX_gg_phi3_phi1Z_bbll_CMS8=gg_phi3_phi1Z_bbll_TH8/ip_ex_gg_phi_phi1Z_bbll_CMS8(m3);
    //     if(THoEX_gg_phi3_phi1Z_bbll_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();
         }
    if(m2>= 220.0 && m2<1000.0)
        {
                                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop159 \033[0m "<<std::endl;
        
         THoEX_gg_phi2_phi1Z_tautauZ_ATLAS8=gg_phi2_phi1Z_tautauZ_TH8/ip_ex_gg_phi_phi1Z_tautauZ_ATLAS8(m2);
    //     if(THoEX_gg_phi2_phi1Z_tautauZ_ATLAS8 >5) return std::numeric_limits<double>::quiet_NaN();
         }
    if(m3>= 220.0 && m3<1000.0)
        {
                                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop160 \033[0m "<<std::endl;
        
         THoEX_gg_phi3_phi1Z_tautauZ_ATLAS8=gg_phi3_phi1Z_tautauZ_TH8/ip_ex_gg_phi_phi1Z_tautauZ_ATLAS8(m3);
    //     if(THoEX_gg_phi3_phi1Z_tautauZ_ATLAS8 >5) return std::numeric_limits<double>::quiet_NaN();
         }
    if(m2>= 220.0 && m2< 350.0)
        {
                                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop161 \033[0m "<<std::endl;
        
         THoEX_gg_phi2_phi1Z_tautaull_CMS8=gg_phi2_phi1Z_tautaull_TH8/ip_ex_gg_phi_phi1Z_tautaull_CMS8(m2);
    //     if(THoEX_gg_phi2_phi1Z_tautaull_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();
         }
    if(m3>= 220.0 && m3< 350.0)
        {
                                                                                                                                                                                 THoEX_gg_phi3_phi1Z_tautaull_CMS8=gg_phi3_phi1Z_tautaull_TH8/ip_ex_gg_phi_phi1Z_tautaull_CMS8(m3);
    
    }
    //if(m2>= 200.0 && m2<2000.0)
    if(m2>= 220.0 && m2<2000.0)
    {
                                                                                                                                                                        
               
        //THoEX_gg_phi2_phi1Z_bbZ_ATLAS13=gg_phi2_phi1Z_bbZ_TH13/ip_ex_gg_phi_phi1Z_bbZ_ATLAS13(m2);
        THoEX_gg_phi2_phi1Z_bbZ_ATLAS13=gg_phi2_phi1Z_with_channel_bbZ_TH13/ip_ex_gg_phi_phi1Z_bbZ_ATLAS13(m2);
    
        
    }
    //if(m3>= 200.0 && m3<2000.0)
    if(m3>= 220.0 && m3<2000.0)
        {
        
        //THoEX_gg_phi3_phi1Z_bbZ_ATLAS13=gg_phi3_phi1Z_bbZ_TH13/ip_ex_gg_phi_phi1Z_bbZ_ATLAS13(m3);
        THoEX_gg_phi3_phi1Z_bbZ_ATLAS13=gg_phi3_phi1Z_with_channel_bbZ_TH13/ip_ex_gg_phi_phi1Z_bbZ_ATLAS13(m3);

        
         }
    if(m2>= 220.0 && m2< 800.0)
        {
                                                                                                                                                                        
        
         THoEX_gg_phi2_phi1Z_bbZ_1_CMS13=gg_phi2_phi1Z_bbZ_TH13/ip_ex_gg_phi_phi1Z_bbZ_1_CMS13(m2);

         }
    if(m3>= 220.0 && m3< 800.0)
        {

         THoEX_gg_phi3_phi1Z_bbZ_1_CMS13=gg_phi3_phi1Z_bbZ_TH13/ip_ex_gg_phi_phi1Z_bbZ_1_CMS13(m3);

         }
    if(m2>= 800.0 && m2<2000.0)
        {
                                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop167 \033[0m "<<std::endl;
        
         THoEX_gg_phi2_phi1Z_bbZ_2_CMS13=gg_phi2_phi1Z_bbZ_TH13/ip_ex_gg_phi_phi1Z_bbZ_2_CMS13(m2);
    //     if(THoEX_gg_phi2_phi1Z_bbZ_2_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
         }
    if(m3>= 800.0 && m3<2000.0)
        {
                                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop168 \033[0m "<<std::endl;
        
         THoEX_gg_phi3_phi1Z_bbZ_2_CMS13=gg_phi3_phi1Z_bbZ_TH13/ip_ex_gg_phi_phi1Z_bbZ_2_CMS13(m3);
    //     if(THoEX_gg_phi3_phi1Z_bbZ_2_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
         }
    //if(m2>= 200.0 && m2<2000.0)
    if(m2>= 220.0 && m2<2000.0)
        {
                                                                                                                                                                        
        //THoEX_bb_phi2_phi1Z_bbZ_ATLAS13=bb_phi2_phi1Z_bbZ_TH13/ip_ex_bb_phi_phi1Z_bbZ_ATLAS13(m2);
        THoEX_bb_phi2_phi1Z_bbZ_ATLAS13=bb_phi2_phi1Z_with_channel_bbZ_TH13/ip_ex_bb_phi_phi1Z_bbZ_ATLAS13(m2);

         }
    //if(m3>= 220.0 && m3<2000.0)
    if(m3>= 220.0 && m3<2000.0)
        {

        //THoEX_bb_phi3_phi1Z_bbZ_ATLAS13=bb_phi3_phi1Z_bbZ_TH13/ip_ex_bb_phi_phi1Z_bbZ_ATLAS13(m3);
        THoEX_bb_phi3_phi1Z_bbZ_ATLAS13=bb_phi3_phi1Z_with_channel_bbZ_TH13/ip_ex_bb_phi_phi1Z_bbZ_ATLAS13(m3);

         }
    
    //ip_ex_bb_phi_phi1Z_tautaull_ATLAS13
    if(m2>= 220.0 && m2<400.0)
        {
                                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop171 \033[0m "<<std::endl;
        
         THoEX_gg_phi2_phi1Z_tautaull_CMS13=gg_phi2_phi1Z_tautaull_TH13/ip_ex_gg_phi_phi1Z_tautaull_CMS13(m2);
    //     if(THoEX_gg_phi2_phi1Z_tautaull_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
         }
    if(m3>= 220.0 && m3<400.0)
        {
                                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop172 \033[0m "<<std::endl;
        
         THoEX_gg_phi3_phi1Z_tautaull_CMS13=gg_phi3_phi1Z_tautaull_TH13/ip_ex_gg_phi_phi1Z_tautaull_CMS13(m3);
    //     if(THoEX_gg_phi3_phi1Z_tautaull_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
         }
    
    
    
    
    
    
    if(m2>= 220.0 && m2< 800.0)
        {
                                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop173 \033[0m "<<std::endl;
        
         THoEX_bb_phi2_phi1Z_bbZ_1_CMS13=bb_phi2_phi1Z_bbZ_TH13/ip_ex_bb_phi_phi1Z_bbZ_1_CMS13(m2);
    //     if(THoEX_bb_phi2_phi1Z_bbZ_1_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
         }
    if(m3>= 220.0 && m3< 800.0)
        {
                                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop174 \033[0m "<<std::endl;
        
         THoEX_bb_phi3_phi1Z_bbZ_1_CMS13=bb_phi3_phi1Z_bbZ_TH13/ip_ex_bb_phi_phi1Z_bbZ_1_CMS13(m3);
    //     if(THoEX_bb_phi3_phi1Z_bbZ_1_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
         }
    if(m2>= 800.0 && m2<2000.0)
        {
                                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop175 \033[0m "<<std::endl;
        
         THoEX_bb_phi2_phi1Z_bbZ_2_CMS13=bb_phi2_phi1Z_bbZ_TH13/ip_ex_bb_phi_phi1Z_bbZ_2_CMS13(m2);
    //     if(THoEX_bb_phi2_phi1Z_bbZ_2_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
         }
    if(m3>= 800.0 && m3<2000.0)
        {
                                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop176 \033[0m "<<std::endl;
        
         THoEX_bb_phi3_phi1Z_bbZ_2_CMS13=bb_phi3_phi1Z_bbZ_TH13/ip_ex_bb_phi_phi1Z_bbZ_2_CMS13(m3);
    //     if(THoEX_bb_phi3_phi1Z_bbZ_2_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
         }
    if(m3>= 175.0 && m3<1000.0 && m2 >=50.0 && m2 <910.0)
        {
                                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop177 \033[0m "<<std::endl;
        
         THoEX_pp_phi3_phi2Z_bbll_1_CMS8=pp_phi3_phi2Z_bbll_TH8/ip_ex_pp_phii_phijZ_bbll_1_CMS8(m3,m2);
    //     if(THoEX_pp_phi3_phi2Z_bbll_1_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();
         } //mA=m3, mH=m2
    if(m2>= 175.0 && m2<1000.0 && m3 >=50.0 && m3 <910.0)
        {
                                                                                                                                                                        
        //std::cout<<"\033[1;31m   stop178 \033[0m "<<std::endl;
        
         THoEX_pp_phi2_phi3Z_bbll_1_CMS8=pp_phi2_phi3Z_bbll_TH8/ip_ex_pp_phii_phijZ_bbll_1_CMS8(m2,m3);
    //     if(THoEX_pp_phi2_phi3Z_bbll_1_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();
         } //mA=m3, mH=m2
    if(m3>=  50.0 && m3<1000.0 && m2 >=50.0 && m2 <1000.0)
        {

        
         THoEX_pp_phi3_phi2Z_tautaull_1_CMS8=pp_phi3_phi2Z_tautaull_TH8/ip_ex_pp_phii_phijZ_tautaull_1_CMS8(m3,m2);

    } //mA=m3, mH=m2
    if(m2>=  50.0 && m2<1000.0 && m3 >=50.0 && m3 <1000.0)
        {

        
         THoEX_pp_phi2_phi3Z_tautaull_1_CMS8=pp_phi2_phi3Z_tautaull_TH8/ip_ex_pp_phii_phijZ_tautaull_1_CMS8(m2,m3);
    
    } //mA=m3, mH=m2
    if(m3>=  50.0 && m3<1000.0 && m2 >=50.0 && m2 <1000.0)
        {

        
         THoEX_pp_phi3_phi2Z_tautaull_2_CMS8=pp_phi3_phi2Z_tautaull_TH8/ip_ex_pp_phii_phijZ_tautaull_2_CMS8(m2,m3);
    
    } //mA=m2, mH=m3
    if(m2>=  50.0 && m2<1000.0 && m3 >=50.0 && m3 <1000.0)
        {
                                                                                                                                                                        
    
         THoEX_pp_phi2_phi3Z_tautaull_2_CMS8=pp_phi2_phi3Z_tautaull_TH8/ip_ex_pp_phii_phijZ_tautaull_2_CMS8(m3,m2);

    } 
    if(m2 >= 230.0 && m2 <800.0 && m3>=130.0 && m3<700.0 && m2-m3>=100)
        {
        
         THoEX_gg_phi2_phi3Z_bbZ_ATLAS13=gg_phi2_phi3Z_bbZ_TH13/ip_ex_gg_phii_phijZ_bbZ_ATLAS13(m3,m2);
         //std::cout<<"\033[1;31m   Enters = \033[0m " <<std::endl;
         //std::cout<<"\033[1;31m   gg_phi3_phi2Z_bbZ_TH13 = \033[0m "<< gg_phi3_phi2Z_bbZ_TH13 <<std::endl;
         //std::cout<<"\033[1;31m   ip_ex_gg_phii_phijZ_bbZ_ATLAS13(m3,m2) = \033[0m "<< ip_ex_gg_phii_phijZ_bbZ_ATLAS13(m2,m3) <<std::endl;

    }
    if(m3 >= 230.0 && m3 <800.0 && m2>=130.0 && m2<700.0 && m3-m2>=100)
        {
                                                                                                                                                                        
        
        THoEX_gg_phi3_phi2Z_bbZ_ATLAS13=gg_phi3_phi2Z_bbZ_TH13/ip_ex_gg_phii_phijZ_bbZ_ATLAS13(m2,m3);
    }
    
    if(m2 >= 230.0 && m2 <800.0 && m3>=130.0 && m3<700.0 && m2-m3>=100)
        {
        
        THoEX_bb_phi2_phi3Z_bbZ_ATLAS13=bb_phi2_phi3Z_bbZ_TH13/ip_ex_bb_phii_phijZ_bbZ_ATLAS13(m3,m2);
    
    }
    if(m3 >= 230.0 && m3 <800.0 && m2>=130.0 && m2<700.0 && m3-m2>=100)
        {
                                                                                                                                                                                
        THoEX_bb_phi3_phi2Z_bbZ_ATLAS13=bb_phi3_phi2Z_bbZ_TH13/ip_ex_bb_phii_phijZ_bbZ_ATLAS13(m2,m3);
         
         
    }
    
    
    if(m2 >= 300.0 && m2 <800.0 && m3>=200.0 && m3<700.0 && m2-m3>=100)
        {
        
        THoEX_gg_phi2_phi3Z_WWZ_ATLAS13=gg_phi2_phi3Z_WWZ_TH13/ip_ex_gg_phii_phijZ_WWZ_ATLAS13(m3,m2);

    }
    if(m3 >= 300.0 && m3 <800.0 && m2>=200.0 && m2<700.0 && m3-m2>=100)
        {
                                                                                                                                                                        
        
        THoEX_gg_phi3_phi2Z_WWZ_ATLAS13=gg_phi3_phi2Z_WWZ_TH13/ip_ex_gg_phii_phijZ_WWZ_ATLAS13(m2,m3);
    }
    
    
    
    if(mHp>= 180.0 && mHp<1000.0)
        {
              
        
        
        
        //std::cout<<"\033[1;31m   pp_Hpm_taunu_TH8 = \033[0m "<< pp_Hpm_taunu_TH8 <<std::endl;
        
        //std::cout<<"\033[1;31m   ip_ex_pp_Hpm_taunu_ATLAS8(mHp) = \033[0m "<< ip_ex_pp_Hpm_taunu_ATLAS8(mHp) <<std::endl;
        
        
        //std::cout<<"\033[1;31m   pp_Hpm_taunu_TH8/ip_ex_pp_Hpm_taunu_ATLAS8(mHp) = \033[0m "<< pp_Hpm_taunu_TH8/ip_ex_pp_Hpm_taunu_ATLAS8(mHp) <<std::endl;
        
        THoEX_pp_Hpm_taunu_ATLAS8=pp_Hpm_taunu_TH8/ip_ex_pp_Hpm_taunu_ATLAS8(mHp);
        
        
        
    //    if(THoEX_pp_Hpm_taunu_ATLAS8 >5) return std::numeric_limits<double>::quiet_NaN();
         }
        
    if(mHp>= 180.0 && mHp< 600.0)
        {
                                                                                                                                                                        
        //////std::cout<<"\033[1;31m   stop188 \033[0m "<<std::endl;
        
        
        //////if(THoEX_pp_Hp_taunu_CMS8 != 0) std::cout<<"\033[1;31m   THoEX_pp_Hp_taunu_CMS8 = \033[0m "<< THoEX_pp_Hp_taunu_CMS8 <<std::endl;
        //std::cout<<"\033[1;31m   THoEX_pp_Hp_taunu_CMS8 = \033[0m "<< THoEX_pp_Hp_taunu_CMS8 <<std::endl;
        
        
        
        THoEX_pp_Hp_taunu_CMS8=pp_Hp_taunu_TH8/ip_ex_pp_Hp_taunu_CMS8(mHp);
        
    //    if(THoEX_pp_Hp_taunu_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();
         }
    if(mHp>= 180.0 && mHp<2000.0)//exp can start in 150 GeV but theoretical no, we should compute more points there
        {
                                                                                                                                                                        
        //////std::cout<<"\033[1;31m   stop189 \033[0m "<<std::endl;
        
        //////if(THoEX_pp_Hpm_taunu_ATLAS13 != 0) std::cout<<"\033[1;31m   THoEX_pp_Hpm_taunu_ATLAS13 = \033[0m "<< THoEX_pp_Hpm_taunu_ATLAS13 <<std::endl;
        //std::cout<<"\033[1;31m   THoEX_pp_Hpm_taunu_ATLAS13 = \033[0m "<< THoEX_pp_Hpm_taunu_ATLAS13 <<std::endl;
        
        
        
        
        
        THoEX_pp_Hpm_taunu_ATLAS13=pp_Hpm_taunu_TH13/ip_ex_pp_Hpm_taunu_ATLAS13(mHp);
    //    if(THoEX_pp_Hpm_taunu_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
         }
    if(mHp>= 180.0 && mHp<3000.0)//exp can start in 80 GeV but theoretical no, we should compute more points there
        {
                                                                                                                                                                        
        //////std::cout<<"\033[1;31m   stop190 \033[0m "<<std::endl;
        
        //if(THoEX_pp_Hpm_taunu_CMS13 != 0) std::cout<<"\033[1;31m   THoEX_pp_Hpm_taunu_CMS13 = \033[0m "<< THoEX_pp_Hpm_taunu_CMS13 <<std::endl;
        
        //std::cout<<"\033[1;31m   THoEX_pp_Hpm_taunu_CMS13 = \033[0m "<< THoEX_pp_Hpm_taunu_CMS13 <<std::endl;

        
        
        THoEX_pp_Hpm_taunu_CMS13=pp_Hpm_taunu_TH13/ip_ex_pp_Hpm_taunu_CMS13(mHp);
    //    if(THoEX_pp_Hpm_taunu_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
         }
    if(mHp>= 200.0 && mHp< 600.0)
        {
                                                                                                                                                                        
        //////std::cout<<"\033[1;31m   stop191 \033[0m "<<std::endl;
        
        //////if(THoEX_pp_Hpm_tb_ATLAS8 != 0) std::cout<<"\033[1;31m   THoEX_pp_Hpm_tb_ATLAS8 = \033[0m "<< THoEX_pp_Hpm_tb_ATLAS8 <<std::endl;

        //std::cout<<"\033[1;31m   THoEX_pp_Hpm_tb_ATLAS8 = \033[0m "<< THoEX_pp_Hpm_tb_ATLAS8 <<std::endl;

        
        
        THoEX_pp_Hpm_tb_ATLAS8=pp_Hpm_tb_TH8/ip_ex_pp_Hpm_tb_ATLAS8(mHp);
    //    if(THoEX_pp_Hpm_tb_ATLAS8 >5) return std::numeric_limits<double>::quiet_NaN();
         }
    if(mHp>= 180.0 && mHp< 600.0)
        {
                                                                                                                                                                        
        ////////std::cout<<"\033[1;31m   stop192 \033[0m "<<std::endl;
        
        ////////if(THoEX_pp_Hp_tb_CMS8 != 0) std::cout<<"\033[1;31m   THoEX_pp_Hp_tb_CMS8 = \033[0m "<< THoEX_pp_Hp_tb_CMS8 <<std::endl;

        //std::cout<<"\033[1;31m   THoEX_pp_Hp_tb_CMS8 = \033[0m "<< THoEX_pp_Hp_tb_CMS8 <<std::endl;

        
        THoEX_pp_Hp_tb_CMS8=pp_Hp_tb_TH8/ip_ex_pp_Hp_tb_CMS8(mHp);
    //    if(THoEX_pp_Hp_tb_CMS8 >5) return std::numeric_limits<double>::quiet_NaN();
         }
    if(mHp>= 200.0 && mHp<2000.0)
        {
                                                                                                                                                                        
        ////////std::cout<<"\033[1;31m   stop193 \033[0m "<<std::endl;
        
        ////////if(THoEX_pp_Hpm_tb_ATLAS13 != 0) std::cout<<"\033[1;31m   THoEX_pp_Hpm_tb_ATLAS13 = \033[0m "<< THoEX_pp_Hpm_tb_ATLAS13 <<std::endl;

        //std::cout<<"\033[1;31m   THoEX_pp_Hpm_tb_ATLAS13 = \033[0m "<< THoEX_pp_Hpm_tb_ATLAS13 <<std::endl;
        //std::cout<<"\033[1;31m   pp_Hpm_tb_TH13 = \033[0m "<< pp_Hpm_tb_TH13 <<std::endl;
        //std::cout<<"\033[1;31m   ip_ex_pp_Hpm_tb_ATLAS13(mHp) = \033[0m "<< ip_ex_pp_Hpm_tb_ATLAS13(mHp) <<std::endl;

        
        THoEX_pp_Hpm_tb_ATLAS13=pp_Hpm_tb_TH13/ip_ex_pp_Hpm_tb_ATLAS13(mHp);
        //std::cout<<"\033[1;31m   THoEX_pp_Hpm_tb_ATLAS13 = \033[0m "<< THoEX_pp_Hpm_tb_ATLAS13 <<std::endl;
        
    //    if(THoEX_pp_Hpm_tb_ATLAS13 >5) return std::numeric_limits<double>::quiet_NaN();
         }
    if(mHp>= 200.0 && mHp<3000.0)
        {
                                                                                                                                                                        
        ////////std::cout<<"\033[1;31m   stop194 \033[0m "<<std::endl;
        
        ////////if(THoEX_pp_Hpm_tb_CMS13 != 0) std::cout<<"\033[1;31m   THoEX_pp_Hpm_tb_CMS13 = \033[0m "<< THoEX_pp_Hpm_tb_CMS13 <<std::endl;

        //std::cout<<"\033[1;31m   THoEX_pp_Hpm_tb_CMS13 = \033[0m "<< THoEX_pp_Hpm_tb_CMS13 <<std::endl;
        
        
        
        
        //std::cout<<"\033[1;31m   GammaHptot = \033[0m "<< GammaHptot <<std::endl;
        //std::cout<<"\033[1;31m   GammaHptb = \033[0m "<< GammaHptb <<std::endl;
        //std::cout<<"\033[1;32m   Br_Hptotb = \033[0m "<< Br_Hptotb <<std::endl;
        //std::cout<<"\033[1;32m   SigmaHpm13 = \033[0m "<< SigmaHpm13 <<std::endl;
        //std::cout<<"\033[1;31m   pp_Hpm_tb_TH13 = \033[0m "<< pp_Hpm_tb_TH13 <<std::endl;
        //std::cout<<"\033[1;31m   ip_ex_pp_Hpm_tb_CMS13(mHp) = \033[0m "<< ip_ex_pp_Hpm_tb_CMS13(mHp) <<std::endl;
        //std::cout<<"\033[1;31m   THoEX_pp_Hpm_tb_CMS13 = \033[0m "<< THoEX_pp_Hpm_tb_CMS13 <<std::endl;
        
        
        
        THoEX_pp_Hpm_tb_CMS13=pp_Hpm_tb_TH13/ip_ex_pp_Hpm_tb_CMS13(mHp);
    //    if(THoEX_pp_Hpm_tb_CMS13 >5) return std::numeric_limits<double>::quiet_NaN();
         }
}

void GeneralTHDMcache::computeLowMass()
{
    /*********************************/
    /* Observables with phi_3, i.e A */
    /*********************************/

    THoEX_pp_h_phi3phi3_mumutautau_CMS13 = 0.0;

    if(mH3 >= 15.0 && mH3 <= 62.0)
    {
        THoEX_pp_h_phi3phi3_mumutautau_CMS13 = (2.0 * pph13 * GTHDM_BR_h_AA * Br_phi3tomumu * Br_phi3totautau) / ip_low_pp_h_phi3phi3_mumutautau_CMS13(mH3);
    }

//  THoEX_pp_h_phi3phi3_bbtautau_CMS13 = 0.0;
//
//  if(mH3 >= 15.0 && mH3 <= 60.0)
//  {
//     THoEX_pp_h_phi3phi3_bbtautau_CMS13 = (2.0 * pph13 * GTHDM_BR_h_AA * Br_phi3tobb * Br_phi3totautau) / ip_low_pp_h_phi3phi3_bbtautau_CMS13(mH3);
//  }

    THoEX_pp_h_phi3phi3_bbtautau_CMS13 = 0.0;

    if(mH3 >= 12.0 && mH3 <= 60.0)
    {
        THoEX_pp_h_phi3phi3_bbtautau_CMS13 = (2.0 * GTHDM_BR_h_AA * Br_phi3tobb * Br_phi3totautau) / ip_low_pp_h_phi3phi3_bbtautau_CMS13(mH3);
    }

//  THoEX_pp_h_phi3phi3_bbmumu_CMS13 = 0.0;
//
//  if(mH3 >= 20.0 && mH3 <= 62.5)
//  {
//      THoEX_pp_h_phi3phi3_bbmumu_CMS13 = (2.0 * pph13 * GTHDM_BR_h_AA * Br_phi3tobb * Br_phi3tomumu) / ip_low_pp_h_phi3phi3_bbmumu_CMS13(mH3);
//  }

    THoEX_pp_h_phi3phi3_bbmumu_CMS13 = 0.0;

    if(mH3 >= 15.0 && mH3 <= 62.0)
    {
        THoEX_pp_h_phi3phi3_bbmumu_CMS13 = (2.0 * GTHDM_BR_h_AA * Br_phi3tobb * Br_phi3tomumu) / ip_low_pp_h_phi3phi3_bbmumu_CMS13(mH3);
    }

    THoEX_pp_h_phi3Z_mumull_CMS13 = 0.0;

    if(mH3 >= 4.2045 && mH3 <= 34.843)
    {
        THoEX_pp_h_phi3Z_mumull_CMS13 = (GTHDM_BR_h_AZ * Br_phi3tomumu) / ip_low_pp_h_phi23Z_mumull_CMS13(mH3);
    }

    THoEX_pp_h_phi3phi3_mumumumu_CMS13 = 0.0;

    if(mH3 >= 4.2 && mH3 <= 59.946)
    {
        THoEX_pp_h_phi3phi3_mumumumu_CMS13 = (GTHDM_BR_h_AA * Br_phi3tomumu * Br_phi3tomumu) / ip_low_pp_h_phi23phi23_mumumumu_CMS13(mH3);
    }

    THoEX_pp_h_phi3phi3_gagagaga_CMS13 = 0.0;

    if(mH3 >= 15.0 && mH3 <= 62.0)
    {
        //The experimental values are in fb whereas the theoretical value is calculated in pb.
        THoEX_pp_h_phi3phi3_gagagaga_CMS13 = (1.0e3) * ((SigSM_pph13 * pph13) * Br_phi3togaga * Br_phi3togaga) / ip_low_pp_h_phi3phi3_gagagaga_CMS13(mH3);
    }

    THoEX_pp_h_phi3phi3_tautautautau_CMS13 = 0.0;

    if(mH3 >= 4.0 && mH3 <= 15.0)
    {
        THoEX_pp_h_phi3phi3_tautautautau_CMS13 = (pph13 * GTHDM_BR_h_AA * Br_phi3totautau * Br_phi3totautau) / ip_low_pp_h_phi3phi3_tautautautau_CMS13(mH3);
    }

    THoEX_pp_ttphi3_ttmumu_CMS13 = 0.0;

    if(mH3 >= 15.0 && mH3 <= 350.0)
    {
        THoEX_pp_ttphi3_ttmumu_CMS13 = (SigmattF_phi3_13 * Br_phi3tomumu) / ip_low_pp_ttphi3_ttmumu_CMS13(mH3);
    }

    THoEX_pp_ttphi3_tttautau_CMS13 = 0.0;

    if(mH3 >= 15.0 && mH3 <= 350.0)
    {
        THoEX_pp_ttphi3_tttautau_CMS13 = (SigmattF_phi3_13 * Br_phi3totautau) / ip_low_pp_ttphi3_tttautau_CMS13(mH3);
    }

    THoEX_Vh_h_phi3phi3_bbbb_CMS13 = 0.0;

    if(mH3 >= 15.0 && mH3 <= 60.0)
    {
        THoEX_Vh_h_phi3phi3_bbbb_CMS13 = (rh_VV * GTHDM_BR_h_AA * Br_phi3tobb * Br_phi3tobb) / ip_low_Vh_h_phi3phi3_bbbb_CMS13(mH3);
    }

    THoEX_pp_h_phi3phi3_bbmumu_ATLAS13 = 0.0;

    if(mH3 >= 16.0 && mH3 <= 62.0)
    {
        THoEX_pp_h_phi3phi3_bbmumu_ATLAS13 = (2.0 * GTHDM_BR_h_AA * Br_phi3tobb * Br_phi3tomumu) / ip_low_pp_h_phi3phi3_bbmumu_ATLAS13(mH3);
    }

    THoEX_gg_h_phi3phi3_mumumumu_ATLAS13 = 0.0;

    if((mH3 >= 1.5 && mH3 <= 2.0)||(mH3 >= 4.5 && mH3 <= 8.5)||(mH3 >= 11.5 && mH3 <= 14.999)||(mH3 >= 15.0 && mH3 <= 60.0))
    {
        //The experimental values are in fb whereas the theoretical value is calculated in pb.
        THoEX_gg_h_phi3phi3_mumumumu_ATLAS13 = (1.0e3) * ((SigmaggF13 * rh_gg) * GTHDM_BR_h_AA * Br_phi3tomumu * Br_phi3tomumu) / ip_low_gg_h_phi23phi23_mumumumu_ATLAS13(mH3);
    }

    THoEX_gg_h_phi3Z_mumull_ATLAS13 = 0.0;

    if(mH3 >= 15.0 && mH3 <= 30.0)
    {
        //The experimental values are in fb whereas the theoretical value is calculated in pb.
        THoEX_gg_h_phi3Z_mumull_ATLAS13 = (1.0e3) * ((SigmaggF13 * rh_gg) * GTHDM_BR_h_AZ * Br_phi3tomumu * BrSM_Ztoll) / ip_low_gg_h_phi23Z_mumull_ATLAS13(mH3);
    }

    THoEX_Vh_h_phi3phi3_bbbb_ATLAS13 = 0.0;

    if(mH3 >= 20.0 && mH3 <= 60.0)
    {
        THoEX_Vh_h_phi3phi3_bbbb_ATLAS13 = ((SigmaVh13 * rh_VV) * GTHDM_BR_h_AA * Br_phi3tobb * Br_phi3tobb) / ip_low_Vh_h_phi23phi23_bbbb_ATLAS13(mH3);
    }

    THoEX_Zh_h_phi3phi3_bbbb_ATLAS13 = 0.0;

    if(mH3 >= 15.0 && mH3 <= 30.0)
    {
        THoEX_Zh_h_phi3phi3_bbbb_ATLAS13 = ((SigmaZh13 * rh_VV) * GTHDM_BR_h_AA * Br_phi3tobb * Br_phi3tobb) / ip_low_Zh_h_phi23phi23_bbbb_ATLAS13(mH3);
    }

    THoEX_pp_h_phi3phi3_bbmumu_ATLAS13_old = 0.0;

    if(mH3 >= 20.0 && mH3 <= 60.0)
    {
        THoEX_pp_h_phi3phi3_bbmumu_ATLAS13_old = (2.0 * ggF_VBF_Vh13 * GTHDM_BR_h_AA * Br_phi3tobb * Br_phi3tomumu) / ip_low_pp_h_phi23phi23_bbmumu_ATLAS13_old(mH3);
    }

    THoEX_pp_h_phi3phi3_gagagg_ATLAS13 = 0.0;

    if((mH3 >= 20.0 && mH3 <= 24.999)||(mH3 >= 25.0 && mH3 <= 34.999)||(mH3 >= 35.0 && mH3 <= 44.999)||(mH3 >= 45.0 && mH3 <= 54.999)||(mH3 >= 55.0 && mH3 <= 60.0))
    {
        THoEX_pp_h_phi3phi3_gagagg_ATLAS13 = (2.0 * pph13 * GTHDM_BR_h_AA * Br_phi3togaga * Br_phi3togg) / ip_low_pp_h_phi23phi23_gagagg_ATLAS13(mH3);
    }

    THoEX_gg_phi3_tautau_ATLAS13_low = 0.0;

    if(mH3 >= 20.0 && mH3 <= 90.0)
    {
        THoEX_gg_phi3_tautau_ATLAS13_low = (SigmaggF_phi3_13 * Br_phi3totautau) / ip_low_gg_phi3_tautau_ATLAS13(mH3);
    }

    THoEX_pp_h_phi3phi3_gagagaga_ATLAS8 = 0.0;

    if(mH3 >= 10.0 && mH3 <= 62.0)
    {
        THoEX_pp_h_phi3phi3_gagagaga_ATLAS8 = (pph8 * GTHDM_BR_h_AA * Br_phi3togaga * Br_phi3togaga) / ip_low_pp_h_phi3phi3_gagagaga_ATLAS8(mH3);
    }

    THoEX_gg_h_phi3phi3_tautautautau_ATLAS8 = 0.0;

    if(mH3 >= 3.7 && mH3 <= 50.0)
    {
        THoEX_gg_h_phi3phi3_tautautautau_ATLAS8 = (rh_gg * GTHDM_BR_h_AA * Br_phi3totautau * Br_phi3totautau) / ip_low_gg_h_phi3phi3_tautautautau_ATLAS8(mH3);
    }

    THoEX_pp_h_phi3phi3_tautautautau_CMS8 = 0.0;

    if(mH3 >= 5.0 && mH3 <= 15.0)
    {
        THoEX_pp_h_phi3phi3_tautautautau_CMS8 = (pph8 * GTHDM_BR_h_AA * Br_phi3totautau * Br_phi3totautau) / ip_low_pp_h_phi3phi3_tautautautau_CMS8(mH3);
    }

    THoEX_pp_h_phi3phi3_bbmumu_CMS8 = 0.0;

    if(mH3 >= 25.0 && mH3 <= 62.5)
    {
        THoEX_pp_h_phi3phi3_bbmumu_CMS8 = (2.0 * pph8 * GTHDM_BR_h_AA * Br_phi3tobb * Br_phi3tomumu) / ip_low_pp_h_phi3phi3_bbmumu_CMS8(mH3);
    }

    THoEX_pp_h_phi3phi3_mumutautau_CMS8 = 0.0;

    if(mH3 >= 15.0 && mH3 <= 62.5)
    {
        THoEX_pp_h_phi3phi3_mumutautau_CMS8 = (2.0 * pph8 * GTHDM_BR_h_AA * Br_phi3totautau * Br_phi3tomumu) / ip_low_pp_h_phi3phi3_mumutautau_CMS8(mH3);
    }

    THoEX_pp_bbphi3_bbtautau_CMS13 = 0.0;

    if(mH3 >= 25.0 && mH3 <= 70.0)
    {
        THoEX_pp_bbphi3_bbtautau_CMS13 = (SigmabbF_phi3_13 * Br_phi3totautau) / ip_low_pp_bbphi3_bbtautau_CMS13(mH3);
    }

    THoEX_pp_bbphi3_bbtautau_CMS8 = 0.0;

    if(mH3 >= 25.0 && mH3 <= 80.0)
    {
        THoEX_pp_bbphi3_bbtautau_CMS8 = (SigmabbF_phi3_8 * Br_phi3totautau) / ip_low_pp_bbphi3_bbtautau_CMS8(mH3);
    }

    THoEX_pp_ttphi3_ttmumu_ATLAS13 = 0.0;

    if(mH3 >= 15.1 && mH3 <= 73.3)
    {
        THoEX_pp_ttphi3_ttmumu_ATLAS13 = (SigmattF_phi3_13 * Br_phi3tomumu) / ip_low_pp_ttphi3_ttmumu_ATLAS13(mH3);
    }

    THoEX_pp_bbphi3_bbmumu_CMS8 = 0.0;

    if(mH3 >= 25.0 && mH3 <= 60.0)
    {
        THoEX_pp_bbphi3_bbmumu_CMS8 = (SigmabbF_phi3_8 * Br_phi3tomumu) / ip_low_pp_bbphi3_bbmumu_CMS8(mH3);
    }

    THoEX_pp_h_phi3phi3_gagagaga_ATLAS13 = 0.0;

    if(mH3 >= 10.0 && mH3 <= 62.0)
    {
        THoEX_pp_h_phi3phi3_gagagaga_ATLAS13 = (GTHDM_BR_h_AA * Br_phi3togaga * Br_phi3togaga) / ip_low_pp_h_phi3phi3_gagagaga_ATLAS13(mH3);
    }

    THoEX_pp_h_phi3phi3_bbtautau_ATLAS13 = 0.0;

    if(mH3 >= 12.0 && mH3 <= 60.0)
    {
        THoEX_pp_h_phi3phi3_bbtautau_ATLAS13 = (2.0 * pph13 * GTHDM_BR_h_AA * Br_phi3tobb * Br_phi3totautau) / ip_low_pp_h_phi3phi3_bbtautau_ATLAS13(mH3);
    }

    /*********************************/
    /* Observables with phi_2, i.e H */
    /*********************************/

    THoEX_pp_h_phi2Z_mumull_CMS13 = 0.0;

    if(mH2 >= 4.2045 && mH2 <= 34.843)
    {
        THoEX_pp_h_phi2Z_mumull_CMS13 = (GTHDM_BR_h_HZ * Br_phi2tomumu) / ip_low_pp_h_phi23Z_mumull_CMS13(mH2);
    }

    THoEX_pp_h_phi2phi2_mumumumu_CMS13 = 0.0;

    if(mH2 >= 4.2 && mH2 <= 59.946)
    {
        THoEX_pp_h_phi2phi2_mumumumu_CMS13 = (GTHDM_BR_h_HH * Br_phi2tomumu * Br_phi2tomumu) / ip_low_pp_h_phi23phi23_mumumumu_CMS13(mH2);
    }

    THoEX_gg_h_phi2phi2_mumumumu_ATLAS13 = 0.0;

    if((mH2 >= 1.5 && mH2 <= 2.0)||(mH2 >= 4.5 && mH2 <= 8.5)||(mH2 >= 11.5 && mH2 <= 14.999)||(mH2 >= 15.0 && mH2 <= 60.0))
    {
        //The experimental values are in fb whereas the theoretical value is calculated in pb.
        THoEX_gg_h_phi2phi2_mumumumu_ATLAS13 = (1.0e3) * ((SigmaggF13 * rh_gg) * GTHDM_BR_h_HH * Br_phi2tomumu * Br_phi2tomumu) / ip_low_gg_h_phi23phi23_mumumumu_ATLAS13(mH2);
    }

    THoEX_gg_h_phi2Z_mumull_ATLAS13 = 0.0;

    if(mH2 >= 15.0 && mH2 <= 30.0)
    {
        //The experimental values are in fb whereas the theoretical value is calculated in pb.
        THoEX_gg_h_phi2Z_mumull_ATLAS13 = (1.0e3) * ((SigmaggF13 * rh_gg) * GTHDM_BR_h_HZ * Br_phi2tomumu * BrSM_Ztoll) / ip_low_gg_h_phi23Z_mumull_ATLAS13(mH2);
    }

    THoEX_Vh_h_phi2phi2_bbbb_ATLAS13 = 0.0;

    if(mH2 >= 20.0 && mH2 <= 60.0)
    {
        THoEX_Vh_h_phi2phi2_bbbb_ATLAS13 = ((SigmaVh13 * rh_VV) * GTHDM_BR_h_HH * Br_phi2tobb * Br_phi2tobb) / ip_low_Vh_h_phi23phi23_bbbb_ATLAS13(mH2);
    }

    THoEX_Zh_h_phi2phi2_bbbb_ATLAS13 = 0.0;

    if(mH2 >= 15.0 && mH2 <= 30.0)
    {
        THoEX_Zh_h_phi2phi2_bbbb_ATLAS13 = ((SigmaZh13 * rh_VV) * GTHDM_BR_h_HH * Br_phi2tobb * Br_phi2tobb) / ip_low_Zh_h_phi23phi23_bbbb_ATLAS13(mH2);
    }

    THoEX_pp_h_phi2phi2_bbmumu_ATLAS13_old = 0.0;

    if(mH2 >= 20.0 && mH2 <= 60.0)
    {
        THoEX_pp_h_phi2phi2_bbmumu_ATLAS13_old = (2.0 * ggF_VBF_Vh13 * GTHDM_BR_h_AA * Br_phi2tobb * Br_phi2tomumu) / ip_low_pp_h_phi23phi23_bbmumu_ATLAS13_old(mH2);
    }

    THoEX_pp_h_phi2phi2_gagagg_ATLAS13 = 0.0;

    if((mH2 >= 20.0 && mH2 <= 24.999)||(mH2 >= 25.0 && mH2 <= 34.999)||(mH2 >= 35.0 && mH2 <= 44.999)||(mH2 >= 45.0 && mH2 <= 54.999)||(mH2 >= 55.0 && mH2 <= 60.0))
    {
        THoEX_pp_h_phi2phi2_gagagg_ATLAS13 = (2.0 * pph13 * GTHDM_BR_h_HH * Br_phi2togaga * Br_phi2togg) / ip_low_pp_h_phi23phi23_gagagg_ATLAS13(mH2);
    }

    THoEX_pp_phi2_gaga_CMS8 = 0.0;

    if(mH2 >= 80.0 && mH2 <= 109.9)
    {
        THoEX_pp_phi2_gaga_CMS8 = (SigmaSumphi2_8 * Br_phi2togaga) / ip_low_pp_phi2_gaga_CMS8(mH2);
    }

    THoEX_pp_phi2_gaga_CMS13 = 0.0;

    if(mH2 >= 70.0 && mH2 <= 110.0)
    {
        THoEX_pp_phi2_gaga_CMS13 = (SigmaSumphi2_13 * Br_phi2togaga) / ip_low_pp_phi2_gaga_CMS13(mH2);
    }

    THoEX_pp_phi2_gaga_ATLAS13_low = 0.0;

    if(mH2 >= 66.0 && mH2 <= 110.0)
    {
        THoEX_pp_phi2_gaga_ATLAS13_low = (SigmaSumphi2_13 * Br_phi2togaga) / ip_low_pp_phi2_gaga_ATLAS13(mH2);
    }

    THoEX_phi2Z_gagaZ_LEP209 = 0.0;

    if(mH2 >= 20.0 && mH2 <= 116.0)
    {
        THoEX_phi2Z_gagaZ_LEP209 = (Br_phi2togaga) / ip_low_phi2Z_gagaZ_LEP209(mH2);
    }

    THoEX_phi2Z_bbZ_LEP209 = 0.0;

    if(mH2 >= 12.0 && mH2 <= 120.0)
    {
        THoEX_phi2Z_bbZ_LEP209 = (rphi2_VV * Br_phi2tobb) / ip_low_phi2Z_bbZ_LEP209(mH2);
    }

    THoEX_phi2Z_tautauZ_LEP209 = 0.0;

    if(mH2 >= 4.0 && mH2 <= 120.0)
    {
        THoEX_phi2Z_tautauZ_LEP209 = (rphi2_VV * Br_phi2totautau) / ip_low_phi2Z_tautauZ_LEP209(mH2);
    }

    /**********************************/
    /* Observables with phi_i + phi_j */
    /**********************************/

    // To generalise for CPV analyses, also the Z-h-H combination has to be considered with the LEP data below

    THoEX_phi2phi3_bbbb_LEP209 = 0.0;
    THoEX_phi2phi3_tautautautau_LEP209 = 0.0;

    if(mH2 <= mH3)
    {
        if(mH2 >= 10.0 && mH3 >= 15.0 && mH3 <= 145.0 && mH2+mH3 <= 205.0)
            THoEX_phi2phi3_bbbb_LEP209 = (R11*R11 * Br_phi2tobb * Br_phi3tobb) / ip_low_phiEphi3_bbbb_LEP209(mH3,mH2,5.,5.);

        if(mH2 >=  5.0 && mH3 >=  5.0 && mH3 <= 150.0 && mH2+mH3 <= 205.0)
            THoEX_phi2phi3_tautautautau_LEP209 = (R11*R11 * Br_phi2totautau * Br_phi3totautau) / ip_low_phiEphi3_tautautautau_LEP209(mH3,mH2,5.,5.);
    }
    else
    {
        if(mH3 >= 10.0 && mH2 >= 15.0 && mH2 <= 145.0 && mH2+mH3 <= 205.0)
            THoEX_phi2phi3_bbbb_LEP209 = (R11*R11 * Br_phi2tobb * Br_phi3tobb) / ip_low_phiEphi3_bbbb_LEP209(mH2,mH3,5.,5.);

        if(mH3 >=  5.0 && mH2 >=  5.0 && mH2 <= 150.0 && mH2+mH3 <= 205.0)
            THoEX_phi2phi3_tautautautau_LEP209 = (R11*R11 * Br_phi2totautau * Br_phi3totautau) / ip_low_phiEphi3_tautautautau_LEP209(mH2,mH3,5.,5.);
    }

    THoEX_phi1phi3_bbbb_LEP209 = 0.0;

    if(mH3 >= 10.0 && mHl+mH3 <= 205.0)
    {
        THoEX_phi1phi3_bbbb_LEP209 = (R12*R12 * GTHDM_BR_h_bb * Br_phi3tobb) / ip_low_phiEphi3_bbbb_LEP209(mHl,mH3,5.,5.);
    }

    THoEX_phi1phi3_tautautautau_LEP209 = 0.0;

    if(mH3 >=  5.0 && mHl+mH3 <= 205.0)
    {
        THoEX_phi1phi3_tautautautau_LEP209 = (R12*R12 * GTHDM_BR_h_tautau * Br_phi3totautau) / ip_low_phiEphi3_tautautautau_LEP209(mHl,mH3,5.,5.);
    }

    /***********************/
    /* Observables with Hp */
    /***********************/

    THoEX_t_Hpb_csb_CMS8 = 0.0;

    if(mHp >= 90.0 && mHp <= 160.0)
    {
        THoEX_t_Hpb_csb_CMS8 = (Br_ttoHpb * Br_Hptocs) / ip_low_t_Hpb_csb_CMS8(mHp);
    }

    THoEX_t_Hpb_taunub_CMS8 = 0.0;

    if(mHp >= 80.0 && mHp <= 160.0)
    {
        THoEX_t_Hpb_taunub_CMS8 = (Br_ttoHpb * Br_Hptotaunu) / ip_low_t_Hpb_taunub_CMS8(mHp);
    }

    THoEX_t_Hpb_cbb_CMS8 = 0.0;

    if(mHp >= 90.0 && mHp <= 150.0)
    {
        THoEX_t_Hpb_cbb_CMS8 = (Br_ttoHpb * Br_Hptocb) / ip_low_t_Hpb_cbb_CMS8(mHp);
    }

    THoEX_t_Hpb_WAb_Wmumub_CMS13 = 0.0;

    if(mH3 >= 15.0 && mH3 <= 75.0 && mHp >= (mH3+85.0) && mHp <= 160.)
    {
        THoEX_t_Hpb_WAb_Wmumub_CMS13 = (Br_ttoHpb * Br_Hptophi3W * Br_phi3tomumu) / ip_low_t_Hpb_WAb_Wmumub_CMS13(mH3);
    }

    THoEX_t_Hpb_csb_CMS13 = 0.0;

    if(mHp >= 80.0 && mHp <= 160.0)
    {
        THoEX_t_Hpb_csb_CMS13 = (Br_ttoHpb * Br_Hptocs) / ip_low_t_Hpb_csb_CMS13(mHp);
    }

    THoEX_t_Hpb_taunub_ATLAS8 = 0.0;

    if(mHp >= 80.0 && mHp <= 160.0)
    {
        THoEX_t_Hpb_taunub_ATLAS8 = (Br_ttoHpb * Br_Hptotaunu) / ip_low_t_Hpb_taunub_ATLAS8(mHp);
    }

    THoEX_t_Hpb_cbb_ATLAS13 = 0.0;

    if(mHp >= 90.0 && mHp <= 150.0)
    {
        THoEX_t_Hpb_cbb_ATLAS13 = (Br_ttoHpb * Br_Hptocb) / ip_low_t_Hpb_cbb_ATLAS13(mHp);
    }

    THoEX_t_Hpb_WAb_Wmumub_ATLAS13 = 0.0;

    if(mH3 >= 15.1 && mH3 <= 70.2 && mHp >= 120.0 && mHp <= 160.0)
    {
        THoEX_t_Hpb_WAb_Wmumub_ATLAS13 = (Br_ttoHpb * Br_Hptophi3W * Br_phi3tomumu) / ip_low_t_Hpb_WAb_Wmumub_ATLAS13(mH3);
    }

    THoEX_t_Hpb_csb_ATLAS13 = 0.0;

    if(mHp >= 60.0 && mHp <= 168.0)
    {
        THoEX_t_Hpb_csb_ATLAS13 = (Br_ttoHpb * Br_Hptocs) / ip_low_t_Hpb_csb_ATLAS13(mHp);
    }

    // GeV^-2 -> pb
    double GeV2pb = 0.389379e9;

    THoEX_HpHm_taunutaunu_LEP209 = 0.0;

    if(mHp >= 43.0 && mHp <= 95.0)
    {
        THoEX_HpHm_taunutaunu_LEP209 = (GeV2pb * Sigma_HpHm_LEP(mHp2,209.) * Br_Hptotaunu * Br_Hptotaunu) / ip_low_HpHm_taunutaunu_LEP209(mHp);
    }

    THoEX_HpHm_qqqq_LEP209 = 0.0;

    if(mHp >= 43.0 && mHp <= 95.0)
    {
        THoEX_HpHm_qqqq_LEP209 = (GeV2pb * Sigma_HpHm_LEP(mHp2,209.) * (Br_Hptocs + Br_Hptocb) * (Br_Hptocs + Br_Hptocb)) / ip_low_HpHm_qqqq_LEP209(mHp);
    }

    THoEX_HpHm_qqtaunu_OPAL209 = 0.0;

    if(mHp >= 50.0 && mHp <= 93.0)
    {
        THoEX_HpHm_qqtaunu_OPAL209 = (2.0 * Br_Hptotaunu * (Br_Hptocs + Br_Hptocb)) / ip_low_HpHm_qqtaunu_OPAL209(mHp);
    }

    THoEX_HpHm_qqtaunu_OPAL172 = 0.0;

    if(mHp >= 40.0 && mHp < 50.0)
    {
        THoEX_HpHm_qqtaunu_OPAL172 = (2.0 * Br_Hptotaunu * (Br_Hptocs + Br_Hptocb)) / ip_low_HpHm_qqtaunu_OPAL172(mHp);
    }

    THoEX_HpHm_AWstAWst_bbbb_OPAL209 = 0.0;
    THoEX_HpHm_taunuAWst_taunubb_OPAL209 = 0.0;

    if(mHp >= 40.0 && mHp <= 93.0 && mH3 >= 12.0 && mH3 <= (mHp - 3.))
    {
        THoEX_HpHm_AWstAWst_bbbb_OPAL209 = (Br_Hptophi3Wst * Br_Hptophi3Wst * Br_phi3tobb * Br_phi3tobb) / ip_low_HpHm_AWstAWst_bbbb_OPAL209(mHp,mH3,1.,1.);

        THoEX_HpHm_taunuAWst_taunubb_OPAL209 = (4.0 * Br_Hptotaunu * Br_Hptophi3Wst * Br_phi3tobb) / ip_low_HpHm_taunuAWst_taunubb_OPAL209(mHp,mH3,1.,1.);
    }

    THoEX_pp_HpHm_taunutaunu_ATLAS13 = 0.0;

    if(mHp >= 80.0 && mHp < 440.0)
    {
        //One can try with LO cross-section too.
        THoEX_pp_HpHm_taunutaunu_ATLAS13 = (ip_csr_HpHm_pp_13_NLO(mHp) * Br_Hptotaunu * Br_Hptotaunu) / ip_susy_pp_HpHm_taunutaunu_ATLAS13(mHp);
    }

    THoEX_pp_HpHm_taunutaunu_CMS13 = 0.0;

    if(mHp >= 90.0 && mHp < 500.0)
    {
        //One can try with LO cross-section too.
        THoEX_pp_HpHm_taunutaunu_CMS13 = (ip_csr_HpHm_pp_13_NLO(mHp) * Br_Hptotaunu * Br_Hptotaunu) / ip_susy_pp_HpHm_taunutaunu_CMS13(mHp);
    }

    THoEX_pp_HpHm_munumunu_ATLAS13 = 0.0;

    if(mHp >= 90.0 && mHp < 350.0)
    {
        //One can try with LO cross-section too.
        THoEX_pp_HpHm_munumunu_ATLAS13 = (ip_csr_HpHm_pp_13_NLO(mHp) * Br_Hptomunu * Br_Hptomunu) / ip_susy_pp_HpHm_munumunu_ATLAS13(mHp);
    }

    THoEX_pp_HpHm_munumunu_CMS13 = 0.0;

    if(mHp >= 100.0 && mHp < 800.0)
    {
        //One can try with LO cross-section too.
        THoEX_pp_HpHm_munumunu_CMS13 = (ip_csr_HpHm_pp_13_NLO(mHp) * Br_Hptomunu * Br_Hptomunu) / ip_susy_pp_HpHm_munumunu_CMS13(mHp);
    }

    THoEX_HpHm_munumunu_LEP208 = 0.0;

    if(mHp >= 45.0 && mHp <= 102.0)
    {
        THoEX_HpHm_munumunu_LEP208 = (GeV2pb * Sigma_HpHm_LEP(mHp2,208.) * Br_Hptomunu * Br_Hptomunu) / ip_susy_HpHm_munumunu_LEP208(mHp);
    }
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


//////// Here we update and define the values of the constants and parameters!!!!
double GeneralTHDMcache::updateCache()
{      
    
    
    //The Standard Model h branching ratios
    //Could be made into mass-dependent formulae in the future
    //SEE IF WE CAN UPDATE THESE NUMBERS. THEY MAY BE OLD!!!!!!!!!
    BrSM_htobb = 5.77e-1;
    BrSM_htotautau = 6.32e-2;
    BrSM_htogaga = 2.28e-3;
    BrSM_htoWW = 2.15e-1;
    BrSM_htoZZ = 2.64e-2;
    BrSM_htogg = 8.57e-2;
    BrSM_htoZga = 1.54e-3;
    BrSM_htocc = 2.91e-2;

    //The Standard Model Z to ll branching ratio
    BrSM_Ztoll = 6.73e-2;


    mHl   = myGTHDM->getMHl();
    m1    = mHl;
    mH1sq = myGTHDM->getmH1sq();
    mH2sq = myGTHDM->getmH2sq();
    mH3sq = myGTHDM->getmH3sq();
    mH2   = sqrt(mH2sq);
    mH3   = sqrt(mH3sq);
    mHp2  = myGTHDM->getmHp2();
    mHp   = sqrt(mHp2);
    vev   = myGTHDM->v();

    m2 = mH2;
    m3 = mH3;

    m1_2 = mH1sq;
    m2_2 = mH2sq;
    m3_2 = mH3sq;

    if (m1_2 < 0. || m2_2 < 0. || m3_2 < 0.)
        return std::numeric_limits<double>::quiet_NaN();


    cosa1=myGTHDM->getcosalpha1();
    sina1=myGTHDM->getsinalpha1();
    tana1=myGTHDM->gettanalpha1();
    
    cosa2=myGTHDM->getcosalpha2();
    sina2=myGTHDM->getsinalpha2();
    cosa3=myGTHDM->getcosalpha3();
    sina3=myGTHDM->getsinalpha3();
    
    lambda1 = myGTHDM->getlambda1();   
    lambda2 = myGTHDM->getlambda2();    
    lambda3 = myGTHDM->getlambda3();    
    lambda4 = myGTHDM->getlambda4();
    Relambda5=myGTHDM->getRelambda5();
    Imlambda5=myGTHDM->getImlambda5();
    Relambda6=myGTHDM->getRelambda6();
    Imlambda6=myGTHDM->getImlambda6();
    Relambda7=myGTHDM->getRelambda7();
    
    Rij_GTHDM.assign(0, 0, cosa1*cosa2);
    Rij_GTHDM.assign(0, 1, sina1*cosa2);
    Rij_GTHDM.assign(0, 2, -sina2);
    Rij_GTHDM.assign(1, 0, cosa1*sina2*sina3 - sina1*cosa3);
    Rij_GTHDM.assign(1, 1, sina1*sina2*sina3 + cosa1*cosa3);
    Rij_GTHDM.assign(1, 2, cosa2*sina3);
    Rij_GTHDM.assign(2, 0, cosa2*sina2*cosa3 + sina1*sina3);
    Rij_GTHDM.assign(2, 1, sina1*sina2*cosa3 - cosa1*sina3);
    Rij_GTHDM.assign(2, 2, cosa2*cosa3);
    
    //What about the case "Heavy case" in which the SM is Higgs is not the light one, for the g-2
    //in that case the \pi/2 phase is included so that the angles are still small perturbations but
    //that is not being done here
    
    // Variables without 'GTHDM' label to simplify expressions with lambdaijk and KaellenFunction
    R11 = Rij_GTHDM(0,0);
    R12 = Rij_GTHDM(0,1);
    R13 = Rij_GTHDM(0,2);
    R21 = Rij_GTHDM(1,0);
    R22 = Rij_GTHDM(1,1);
    R23 = Rij_GTHDM(1,2);
    R31 = Rij_GTHDM(2,0);
    R32 = Rij_GTHDM(2,1);
    R33 = Rij_GTHDM(2,2);
    
    
    
    //std::cout<<"\033[1;33m R11 = \033[0m "<<R11<<std::endl;
    //std::cout<<"\033[1;33m R12 = \033[0m "<<R12<<std::endl;
    //std::cout<<"\033[1;33m R13 = \033[0m "<<R13<<std::endl;
    //std::cout<<"\033[1;33m R21 = \033[0m "<<R21<<std::endl;
    //std::cout<<"\033[1;33m R22 = \033[0m "<<R22<<std::endl;
    //std::cout<<"\033[1;33m R23 = \033[0m "<<R23<<std::endl;
    //std::cout<<"\033[1;33m R31 = \033[0m "<<R31<<std::endl;
    //std::cout<<"\033[1;33m R32 = \033[0m "<<R32<<std::endl;
    //std::cout<<"\033[1;33m R33 = \033[0m "<<R33<<std::endl;

    
    
    //LOOK AT THIS, NOT SURE IF IT MAKES SENSE NOW
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
    /*
    m11sq     = M11_2 - M33_2 - M12_2*tanb + 0.5*Relambda5*vev*vev
                + (M33_2-0.5*Relambda5*vev*vev)*(cosb*cosb-sinb*sinb)
                + 0.5*vev*vev*((Relambda6-Relambda7)*sinb*cosb+Relambda7*tanb);

    m22sq     = M11_2 - M33_2 + M12_2/tanb + 0.5*Relambda5*vev*vev
                - (M33_2-0.5*Relambda5*vev*vev)*(cosb*cosb-sinb*sinb)
                + 0.25*vev*vev*(Relambda6+Relambda7+(Relambda6-Relambda7)*(cosb*cosb-sinb*sinb))/tanb;

    Rem12sq   = 0.25*vev*vev*(Relambda6+Relambda7+(Relambda6-Relambda7)*(cosb*cosb-sinb*sinb))
                - (2.0*M33_2-Relambda5*vev*vev)*sinb*cosb;

    Imm12sq   = M13_2;
    */
    
    //Let's the define the parameters of the potential in terms of the linear independent parameters
    //Here we'll assume CP-conservation and therefore we'll take the imaginary parts to zero explicitly
    //For the CP-violating case we should decide which parameters are the best, probably the easiest is to remove one phase of lambda5 and lambda6 and keep the one for lambda7
    //doing we'll simplify the rotation matrix although we'll keep one parameter that seems to be less physical (ImLambda7)
    
    /*
    mH1sq=mHl*mHl;
    mH2sq=myGTHDM->getmH2sq();
    mH3sq=myGTHDM->getmH3sq();
    mHp2=myGTHDM->getmHp2();
    */
    
   
    
    //std::cout<<"\033[1;33m vev = \033[0m "<<vev<<std::endl;
    
    // This was clearly wrong, keep it to try to find where exactly it was coming from
    //lambda1   = (-2.0*(M11_2-M22_2+M33_2) + Relambda5*vev*vev
    //             - (2.0*M22_2-2.0*M33_2+Relambda5*vev*vev)/(cosb*cosb)
    //             + (4.0*M12_2-2.0*Relambda6*vev*vev)*tanb)/(vev*vev);
    
    //lambda2   = (-2.0*(M11_2-M22_2+M33_2) + Relambda5*vev*vev
    //             - (2.0*M22_2-2.0*M33_2+Relambda5*vev*vev)/(sinb*sinb)
    //             - (4.0*M12_2+2.0*Relambda7*vev*vev)/tanb)/(vev*vev);
    
    //lambda3   = -(2.0*(M11_2-M22_2-M33_2-mHp2) + Relambda5*vev*vev
    //              + (2.0*M12_2+Relambda6*vev*vev)/tanb
    //              - (2.0*M12_2-Relambda7*vev*vev)*tanb)/(vev*vev);
    
    //lambda4   = Relambda5 - (2.0*mHp2+4.0*M33_2)/(vev*vev);
    
    //Imlambda6 = (2.0*M13_2-(2.0*M23_2+0.5*Imlambda5*vev*vev)*tanb)/(vev*vev);

    //Imlambda7 = 2.0*M13_2/(vev*vev) + (-0.5*Imlambda5+(2.0*M23_2)/(vev*vev))/tanb;

    //Higgs potential parameters
    /*
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
    */
   
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
        
        Nd_GTHDM.assign(0,0, sigmad_ATHDM*Md_GTHDM(0,0));
        Nd_GTHDM.assign(1,1, sigmad_ATHDM*Md_GTHDM(1,1));
        Nd_GTHDM.assign(2,2, sigmad_ATHDM*Md_GTHDM(2,2));
        
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
    //All this matrices where suppose to be defined in the basis where we have two vevs, 
    //such a basis doesn't make sense in the Align (or general Align) THDM, we remove them
    /*
    Yu1_GTHDM = (cosb*Mu_GTHDM - sinb*Nu_GTHDM)*sqrt(2.)/vev;
    Yu2_GTHDM = (cosb*Nu_GTHDM + sinb*Mu_GTHDM)*sqrt(2.)/vev;
    Yd1_GTHDM = (cosb*Md_GTHDM - sinb*Nd_GTHDM)*sqrt(2.)/vev;
    Yd2_GTHDM = (cosb*Nd_GTHDM + sinb*Md_GTHDM)*sqrt(2.)/vev;
    Yl1_GTHDM = (cosb*Ml_GTHDM - sinb*Nl_GTHDM)*sqrt(2.)/vev;
    Yl2_GTHDM = (cosb*Nl_GTHDM + sinb*Ml_GTHDM)*sqrt(2.)/vev;
    */ 
    
    /*up, down and leptonic couplings */
    su = myGTHDM->getNu_11();
    sd = myGTHDM->getNd_11();
    sl = myGTHDM->getNl_11();

    //if (yu1R*yu1R >4. ||  yd1R*yd1R >4. ||  yl1R*yl1R >4.) 
         //       return std::numeric_limits<double>::quiet_NaN();


    Q_GTHDM=myGTHDM->getQ_GTHDM();
    Ale=myGTHDM->getAle();
    Als=myGTHDM->getAlsMz();
    MZ=myGTHDM->getMz();
    MW=MWGTHDM(myGTHDM->Mw_tree());
    cW2=cW2GTHDM(myGTHDM->c02());
    sW2=1.0-cW2;
    GF=1.0/(sqrt(2.0)*vev*vev);

    Mt   = myGTHDM->getQuarks(QCD::TOP).getMass();
    Mb   = myGTHDM->getQuarks(QCD::BOTTOM).getMass();
    Mc   = myGTHDM->getQuarks(QCD::CHARM).getMass();
    Ms   = myGTHDM->getQuarks(QCD::STRANGE).getMass();
    Mu   = myGTHDM->getQuarks(QCD::UP).getMass();
    Md   = myGTHDM->getQuarks(QCD::DOWN).getMass();
    Mtau = myGTHDM->getLeptons(StandardModel::TAU).getMass();
    Mmu  = myGTHDM->getLeptons(StandardModel::MU).getMass();
    Me   = myGTHDM->getLeptons(StandardModel::ELECTRON).getMass();


    runGeneralTHDMparameters();
    computeSignalStrengths();
    computephi2quantities();
    computephi3quantities();
    computeHpquantities();
    computeWZquantities();
    computeHeavyHiggs();
    computeLowMass();

    return mH1sq;
}
