/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THDMCACHE_H
#define	THDMCACHE_H

#include <cmath>
#include "PVfunctions.h"
#include "THDM.h"
#include "THDMquantities.h"
#include "Runner.h"

#include <stdexcept>
#include "gslpp.h"

/**
 * @class THDMcache
 * @ingroup THDM
 * @brief A class for the caching of some %THDM objects.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details At the moment the Passarino-Veltman functions for %STU and the calculation of the Higgs observables (@f$h@f$ signal strengths and searches for @f$H@f$ and @f$A@f$) are cached.
 * The tables are also read here.
 */
class THDMcache {
    
public:

    /**
     * @brief THDMcache constructor.
     * @details Reads all the tables values and stores them in the memory.
     */
    THDMcache(const StandardModel& SM_i);

    /**
     * @brief THDMcache destructor.
     */
    ~THDMcache();
    
    /**
     * @brief Cache size.
     * @details Determines the size of the cache. If it is set to 5, the cache will remember the last five function calls and store their results.
     */
    static const int CacheSize = 5;

    /**
     * @brief Check whether for the latest set of parameters a value is in the cache.
     * @details Takes a complex value.
     */
    int CacheCheck(const gslpp::complex cache[][CacheSize],
                   const int NumPar, const double params[]) const;

    /**
     * @brief Check whether for the latest set of parameters a value is in the cache.
     * @details Takes a real value.
     */
    int CacheCheckReal(const double cache[][CacheSize],
                   const int NumPar, const double params[]) const;

    /**
     * @brief Adds a new result and its parameters into the cache.
     * @details The new values are added on top. The oldest set on the stack is deleted. Takes a complex value.
     */
    void CacheShift(gslpp::complex cache[][CacheSize], const int NumPar,
                    const double params[], const gslpp::complex newResult) const; 

    /**
     * @brief Adds a new result and its parameters into the cache.
     * @details The new values are added on top. The oldest set on the stack is deleted. Takes a real value.
     */
    void CacheShiftReal(double cache[][CacheSize], const int NumPar,
                    const double params[], const double newResult) const; 

    /**
     * @return an object of PVfunctions class
     */
    const PVfunctions getPV() const {
        return PV;
    }

    /**
     * @brief This function reads values from a table and returns them as an array.
     * @return the tabled values
     */
    gslpp::matrix<double> readTable(std::string filename, int rowN, int colN);

    /**
     * @brief Linearly interpolates a table with one parameter dimension.
     * @return the interpolated value
     */
    double interpolate (gslpp::matrix<double> arrayTab, double x);

    /**
     * @brief Linearly interpolates a table with two parameter dimensions.
     * @return the interpolated value
     */
    double interpolate2D (gslpp::matrix<double> arrayTab, double x, double y);

    /**
     * @brief Fills all required arrays with the values read from the tables.
     */
    void read();

    /**
     * @brief SM Higgs branching ratio tables (obtained with HDECAY 6.10), depending on the Higgs mass.
     */
    gslpp::matrix<double> br_tt, br_bb, br_tautau, br_cc, br_mumu, br_ZZ, br_WW;

    /**
     * @brief Total SM decay width (obtained with HDECAY 6.10), depending on the Higgs mass.
     */
    gslpp::matrix<double> GammaHtot_SM;

    /**
     * @brief SM Higgs production cross section tables at 8 TeV from the LHC Higgs Cross Section Working Group, depending on the Higgs mass.
     */
    gslpp::matrix<double> log_cs_ggH_8, log_cs_VBF_8, log_cs_WH_8, log_cs_ZH_8;

    /**
     * @brief SM Higgs production cross section tables at 13 TeV from the LHC Higgs Cross Section Working Group, depending on the Higgs mass.
     */
    gslpp::matrix<double> log_cs_ggH_13, log_cs_VBF_13, log_cs_WH_13, log_cs_ZH_13;

    /**
     * @brief SM Higgs production cross section table at 8 TeV obtained with MadGraph 5, depending on the Higgs mass.
     */
    gslpp::matrix<double> log_cs_ttH_8;

    /**
     * @brief SM Higgs production cross section table at 13 TeV obtained with MadGraph 5, depending on the Higgs mass.
     */
    gslpp::matrix<double> log_cs_ttH_13;

    /**
     * @brief SM Higgs production cross section table at 8 TeV obtained with SusHi 1.5, depending on the Higgs mass.
     */
    gslpp::matrix<double> log_cs_bbH_8;

    /**
     * @brief SM Higgs production cross section table at 13 TeV obtained with SusHi 1.5, depending on the Higgs mass.
     */
    gslpp::matrix<double> log_cs_bbH_13;

    /**
     * @brief CP-odd Higgs production cross section tables at 8 TeV obtained with HIGLU 4.34, depending on the Higgs mass.
     */
    gslpp::matrix<double> log_cs_ggA_8, log_cs_ttA_8, log_cs_bbA_8;

    /**
     * @brief CP-odd Higgs production cross section tables at 13 TeV obtained with HIGLU 4.34, depending on the Higgs mass.
     */
    gslpp::matrix<double> log_cs_ggA_13, log_cs_ttA_13, log_cs_bbA_13;

    /**
     * @brief Charged Higgs production cross section table at 8 TeV from LHCHXSWGMSSMCharged, depending on the charged Higgs mass and logtb.
     */
    gslpp::matrix<double> log_cs_ggHp_8;

    /**
     * @brief Charged Higgs production cross section table at 13 TeV from LHCHXSWGMSSMCharged, depending on the charged Higgs mass and logtb.
     */
    gslpp::matrix<double> log_cs_ggHp_13;

    /**
     * @brief Production cross section ratio tables at 8 TeV obtained with HIGLU 4.34, depending on the Higgs mass.
     */
    gslpp::matrix<double> csrH_top_8, csrH_bottom_8, csrA_top_8, csrA_bottom_8;

    /**
     * @brief Production cross section ratio tables at 13 TeV obtained with HIGLU 4.34, depending on the Higgs mass.
     */
    gslpp::matrix<double> csrH_top_13, csrH_bottom_13, csrA_top_13, csrA_bottom_13;

    /**
     * @brief ATLAS observed @f$95\%@f$ upper cross section limits at 8 TeV, depending on the Higgs mass.
     */
    gslpp::matrix<double> ATLAS8_pp_phi_gaga, ATLAS8_pp_phi_Zga_llga, ATLAS8_gg_phi_tautau, ATLAS8_bb_phi_tautau, ATLAS8_gg_A_hZ_tautauZ, ATLAS8_gg_A_hZ_bbZ, ATLAS8_gg_phi_tt, ATLAS8_gg_H_WW, ATLAS8_VBF_H_WW, ATLAS8_gg_H_ZZ, ATLAS8_VBF_H_ZZ, ATLAS8_gg_H_hh;

    /**
     * @brief ATLAS expected @f$95\%@f$ upper cross section limits at 8 TeV, depending on the Higgs mass.
     */
    gslpp::matrix<double> ATLAS8_pp_phi_gaga_e, ATLAS8_pp_phi_Zga_llga_e, ATLAS8_gg_phi_tautau_e, ATLAS8_bb_phi_tautau_e, ATLAS8_gg_A_hZ_tautauZ_e, ATLAS8_gg_A_hZ_bbZ_e, ATLAS8_gg_phi_tt_e, ATLAS8_gg_H_WW_e, ATLAS8_VBF_H_WW_e, ATLAS8_gg_H_ZZ_e, ATLAS8_VBF_H_ZZ_e, ATLAS8_gg_H_hh_e;

    /**
     * @brief CMS observed @f$95\%@f$ upper signal strength limits at 8 TeV, depending on the Higgs mass.
     */
    gslpp::matrix<double> CMS8_mu_pp_H_VV;

    /**
     * @brief CMS expected @f$95\%@f$ upper signal strength limits at 8 TeV, depending on the Higgs mass.
     */
    gslpp::matrix<double> CMS8_mu_pp_H_VV_e;

    /**
     * @brief CMS observed @f$95\%@f$ upper cross section limits at 8 TeV, depending on the Higgs mass.
     */
    gslpp::matrix<double> CMS8_gg_A_hZ_bbll, CMS8_pp_H_hh, CMS8_pp_H_hh_gagabb, CMS8_pp_H_hh_bbbb, CMS8_bb_phi_bb, CMS8_gg_phi_tautau, CMS8_bb_phi_tautau, CMS8_gg_phi_gaga, CMS8_pp_A_Zga_llga, CMS8_pp_phi_Zga, CMS8_gg_H_hh_bbtautau, CMS8_gg_A_hZ_tautaull;

    /**
     * @brief CMS expected @f$95\%@f$ upper cross section limits at 8 TeV, depending on the Higgs mass.
     */
    gslpp::matrix<double> CMS8_gg_A_hZ_bbll_e, CMS8_pp_H_hh_e, CMS8_pp_H_hh_gagabb_e, CMS8_pp_H_hh_bbbb_e, CMS8_bb_phi_bb_e, CMS8_gg_phi_tautau_e, CMS8_bb_phi_tautau_e, CMS8_gg_phi_gaga_e, CMS8_pp_A_Zga_llga_e, CMS8_gg_H_hh_bbtautau_e, CMS8_gg_A_hZ_tautaull_e;

//    gslpp::matrix<double> CMS_ggF_phi_gaga_ep2, CMS_ggF_phi_gaga_em2;

    /**
     * @brief ATLAS observed @f$95\%@f$ upper cross section limits at 13 TeV, depending on the Higgs mass.
     */
    gslpp::matrix<double> ATLAS13_bb_phi_tt, ATLAS13_tt_phi_tt, ATLAS13_gg_phi_tautau, ATLAS13_bb_phi_tautau,\
                          ATLAS13_pp_phi_gaga, ATLAS13_pp_phi_Zga, ATLAS13_gg_H_ZZ_llnunu, ATLAS13_gg_H_ZZ_llll,\
                          ATLAS13_VBF_H_ZZ_llll, ATLAS13_gg_H_ZZ_llqq, ATLAS13_VBF_H_ZZ_llqq, ATLAS13_gg_H_ZZ_nunuqq,\
                          ATLAS13_gg_H_WW_enumumu, ATLAS13_VBF_H_WW_enumumu, ATLAS13_gg_H_WW_lnuqq, ATLAS13_pp_H_hh_bbbb,\
                          ATLAS13_pp_H_hh_gagabb, ATLAS13_pp_H_hh_gagaWW, ATLAS13_gg_A_Zh_Zbb, ATLAS13_bb_A_Zh_Zbb;

    /**
     * @brief ATLAS expected @f$95\%@f$ upper cross section limits at 13 TeV, depending on the Higgs mass.
     */
    gslpp::matrix<double> ATLAS13_bb_phi_tt_e, ATLAS13_tt_phi_tt_e, ATLAS13_gg_phi_tautau_e, ATLAS13_bb_phi_tautau_e,\
                          ATLAS13_pp_phi_gaga_e, ATLAS13_pp_phi_Zga_e, ATLAS13_gg_H_ZZ_llnunu_e, ATLAS13_gg_H_ZZ_llll_e,\
                          ATLAS13_VBF_H_ZZ_llll_e, ATLAS13_gg_H_ZZ_llqq_e, ATLAS13_VBF_H_ZZ_llqq_e, ATLAS13_gg_H_ZZ_nunuqq_e,\
                          ATLAS13_gg_H_WW_enumumu_e, ATLAS13_VBF_H_WW_enumumu_e, ATLAS13_gg_H_WW_lnuqq_e, ATLAS13_pp_H_hh_bbbb_e,\
                          ATLAS13_pp_H_hh_gagabb_e, ATLAS13_pp_H_hh_gagaWW_e, ATLAS13_gg_A_Zh_Zbb_e, ATLAS13_bb_A_Zh_Zbb_e;

    /**
     * @brief CMS observed @f$95\%@f$ upper cross section limits at 13 TeV, depending on the Higgs mass.
     */
    gslpp::matrix<double> CMS13_pp_phi_bb, CMS13_gg_phi_tautau, CMS13_bb_phi_tautau, CMS13_gg_phi_gaga, CMS13_pp_phi_Zga_llga,\
                          CMS13_pp_phi_Zga_qqga, CMS13_pp_H_ZZ_llll, CMS13_VBFVH_H_ZZ_llll, CMS13_pp_H_ZZ_llqq, CMS13_ggFVBF_H_WW_lnulnu,\
                          CMS13_pp_H_hh_bbbb, CMS13_pp_H_hh_gagabb, CMS13_pp_H_hh_bbtautau, CMS13_pp_H_hh_bbtautau1, CMS13_pp_H_hh_bblnulnu, CMS13_pp_H_hh_bbVV;

    /**
     * @brief CMS expected @f$95\%@f$ upper cross section limits at 13 TeV, depending on the Higgs mass.
     */
    gslpp::matrix<double> CMS13_pp_phi_bb_e, CMS13_gg_phi_tautau_e, CMS13_bb_phi_tautau_e, CMS13_gg_phi_gaga_e, CMS13_pp_phi_Zga_llga_e,\
                          CMS13_pp_phi_Zga_qqga_e, CMS13_pp_H_ZZ_llll_e, CMS13_VBFVH_H_ZZ_llll_e, CMS13_pp_H_ZZ_llqq_e, CMS13_ggFVBF_H_WW_lnulnu_e,\
                          CMS13_pp_H_hh_bbbb_e, CMS13_pp_H_hh_gagabb_e, CMS13_pp_H_hh_bbtautau_e, CMS13_pp_H_hh_bbtautau1_e, CMS13_pp_H_hh_bblnulnu_e, CMS13_pp_H_hh_bbVV_e;

    /**
     * @brief data templates
     */
    gslpp::matrix<double> temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10;
    gslpp::matrix<double> temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19, temp20;
    gslpp::matrix<double> temp21, temp22, temp23, temp24, temp25, temp26, temp27, temp28, temp29, temp30;
    gslpp::matrix<double> temp31, temp32, temp33, temp34, temp35, temp36, temp37, temp38, temp39, temp40;
    gslpp::matrix<double> temp1e, temp2e, temp3e, temp4e, temp5e, temp6e, temp7e, temp8e, temp9e, temp10e;
    gslpp::matrix<double> temp11e, temp12e, temp13e, temp14e, temp15e, temp16e, temp17e, temp18e, temp19e, temp20e;
    gslpp::matrix<double> temp21e, temp22e, temp23e, temp24e, temp25e, temp26e, temp27e, temp28e, temp29e, temp30e;
    gslpp::matrix<double> temp31e, temp32e, temp33e, temp34e, temp35e, temp36e, temp37e, temp38e, temp39e, temp40e;

    /**
     * @brief ATLAS observed @f$95\%@f$ upper cross section limits at 8 TeV, depending on the charged Higgs mass.
     */
    gslpp::matrix<double> ATLAS8_pp_Hpm_taunu, ATLAS8_pp_Hp_tb;

    /**
     * @brief ATLAS expected @f$95\%@f$ upper cross section limits at 8 TeV, depending on the charged Higgs mass.
     */
    gslpp::matrix<double> ATLAS8_pp_Hpm_taunu_e, ATLAS8_pp_Hp_tb_e;

    /**
     * @brief CMS observed @f$95\%@f$ upper cross section limits at 8 TeV, depending on the charged Higgs mass.
     */
    gslpp::matrix<double> CMS8_pp_Hp_taunu, CMS8_pp_Hp_tb;

    /**
     * @brief CMS expected @f$95\%@f$ upper cross section limits at 8 TeV, depending on the charged Higgs mass.
     */
    gslpp::matrix<double> CMS8_pp_Hp_taunu_e, CMS8_pp_Hp_tb_e;

    /**
     * @brief ATLAS observed @f$95\%@f$ upper cross section limits at 13 TeV, depending on the charged Higgs mass.
     */
    gslpp::matrix<double> ATLAS13_pp_Hpm_taunu, ATLAS13_pp_Hp_tb1, ATLAS13_pp_Hp_tb2;

    /**
     * @brief ATLAS expected @f$95\%@f$ upper cross section limits at 13 TeV, depending on the charged Higgs mass.
     */
    gslpp::matrix<double> ATLAS13_pp_Hpm_taunu_e, ATLAS13_pp_Hp_tb1_e, ATLAS13_pp_Hp_tb2_e;

    /**
     * @brief CMS observed @f$95\%@f$ upper cross section limits at 13 TeV, depending on the charged Higgs mass.
     */
    gslpp::matrix<double> CMS13_pp_Hpm_taunu;

    /**
     * @brief CMS expected @f$95\%@f$ upper cross section limits at 13 TeV, depending on the charged Higgs mass.
     */
    gslpp::matrix<double> CMS13_pp_Hpm_taunu_e;

    /**
     * @brief @f$b\to s \gamma@f$ table, depending on logtb and the logarithm of the charged Higgs mass.
     */
    gslpp::matrix<double> arraybsgamma;

    /**
     * @brief Interpolating function for the SM branching ratio to two top quarks.
     * @return @f$BR^{\text{SM}}(H\to t\bar t)@f$
     */
    double ip_Br_HPtott(double mass);

    /**
     * @brief Interpolating function for the SM branching ratio to two bottom quarks.
     * @return @f$BR^{\text{SM}}(H\to b\bar b)@f$
     */
    double ip_Br_HPtobb(double mass);

    /**
     * @brief Interpolating function for the SM branching ratio to two tau leptons.
     * @return @f$BR^{\text{SM}}(H\to \tau\tau)@f$
     */
    double ip_Br_HPtotautau(double mass);

    /**
     * @brief Interpolating function for the SM branching ratio to two charm quarks.
     * @return @f$BR^{\text{SM}}(H\to c\bar c)@f$
     */
    double ip_Br_HPtocc(double mass);

    /**
     * @brief Interpolating function for the SM branching ratio to two muons.
     * @return @f$BR^{\text{SM}}(H\to \mu \mu)@f$
     */
    double ip_Br_HPtomumu(double mass);

    /**
     * @brief Interpolating function for the SM branching ratio to two @f$Z@f$ bosons.
     * @return @f$BR^{\text{SM}}(H\to ZZ)@f$
     */
    double ip_Br_HPtoZZ(double mass);

    /**
     * @brief Interpolating function for the SM branching ratio to two @f$W@f$ bosons.
     * @return @f$BR^{\text{SM}}(H\to WW)@f$
     */
    double ip_Br_HPtoWW(double mass);

    /**
     * @brief Interpolating function for the total SM Higgs decay width.
     * @return @f$\Gamma^{\text{tot}}_H@f$
     */
    double ip_GammaHPtotSM(double mass);

    /**
     * @brief Interpolating function for the H production cross section via gluon-gluon fusion at 8 TeV.
     * @return @f$\sigma(gg\to H)@f$
     */
    double ip_cs_ggtoH_8(double mass);

    /**
     * @brief Interpolating function for the H production cross section via gluon-gluon fusion at 13 TeV.
     * @return @f$\sigma(gg\to H)@f$
     */
    double ip_cs_ggtoH_13(double mass);

    /**
     * @brief Interpolating function for the H production cross section via vector boson fusion at 8 TeV.
     * @return @f$\sigma(VV\to H)@f$
     */
    double ip_cs_VBFtoH_8(double mass);

    /**
     * @brief Interpolating function for the H production cross section via vector boson fusion at 13 TeV.
     * @return @f$\sigma(VV\to H)@f$
     */
    double ip_cs_VBFtoH_13(double mass);

    /**
     * @brief Interpolating function for the W associated H production cross section at 8 TeV.
     * @return @f$\sigma(W\to WH)@f$
     */
    double ip_cs_WtoWH_8(double mass);

    /**
     * @brief Interpolating function for the W associated H production cross section at 13 TeV.
     * @return @f$\sigma(W\to WH)@f$
     */
    double ip_cs_WtoWH_13(double mass);

    /**
     * @brief Interpolating function for the Z associated H production cross section at 8 TeV.
     * @return @f$\sigma(Z\to ZH)@f$
     */
    double ip_cs_ZtoZH_8(double mass);

    /**
     * @brief Interpolating function for the Z associated H production cross section at 13 TeV.
     * @return @f$\sigma(Z\to ZH)@f$
     */
    double ip_cs_ZtoZH_13(double mass);

    /**
     * @brief Interpolating function for the top associated H production cross section at 8 TeV.
     * @return @f$\sigma(pp\to t \bar t H)@f$
     */
    double ip_cs_pptottH_8(double mass);

    /**
     * @brief Interpolating function for the top associated H production cross section at 13 TeV.
     * @return @f$\sigma(pp\to t \bar t H)@f$
     */
    double ip_cs_pptottH_13(double mass);

    /**
     * @brief Interpolating function for the bottom associated H production cross section at 8 TeV.
     * @return @f$\sigma(pp\to b \bar b H)@f$
     */
    double ip_cs_pptobbH_8(double mass);

    /**
     * @brief Interpolating function for the bottom associated H production cross section at 13 TeV.
     * @return @f$\sigma(pp\to b \bar b H)@f$
     */
    double ip_cs_pptobbH_13(double mass);

    /**
     * @brief Interpolating function for the A production cross section via gluon-gluon fusion at 8 TeV.
     * @return @f$\sigma(gg\to A)@f$
     */
    double ip_cs_ggtoA_8(double mass);

    /**
     * @brief Interpolating function for the A production cross section via gluon-gluon fusion at 13 TeV.
     * @return @f$\sigma(gg\to A)@f$
     */
    double ip_cs_ggtoA_13(double mass);

    /**
     * @brief Interpolating function for the top associated A production cross section at 8 TeV.
     * @return @f$\sigma(pp\to t \bar t A)@f$
     */
    double ip_cs_pptottA_8(double mass);

    /**
     * @brief Interpolating function for the top associated A production cross section at 13 TeV.
     * @return @f$\sigma(pp\to t \bar t A)@f$
     */
    double ip_cs_pptottA_13(double mass);

    /**
     * @brief Interpolating function for the bottom associated A production cross section at 8 TeV.
     * @return @f$\sigma(pp\to b \bar b A)@f$
     */
    double ip_cs_pptobbA_8(double mass);

    /**
     * @brief Interpolating function for the bottom associated A production cross section at 13 TeV.
     * @return @f$\sigma(pp\to b \bar b A)@f$
     */
    double ip_cs_pptobbA_13(double mass);

    /**
     * @brief Interpolating function for the H+ production cross section from two gluons at 8 TeV.
     * @return @f$\sigma(gg\to H^+)@f$
     */
    double ip_cs_ggtoHp_8(double mHp, double logtb);

    /**
     * @brief Interpolating function for the H+ production cross section from two gluons at 13 TeV.
     * @return @f$\sigma(gg\to H^+)@f$
     */
    double ip_cs_ggtoHp_13(double mHp, double logtb);

    /**
     * @brief Interpolating function for the gluon-gluon fusion H cross section ratio of the top-loop and the total contribution at 8 TeV.
     * @return @f$\sigma_t(gg\to H)/\sigma(gg\to H)@f$
     */
    double ip_csr_ggH_t_8(double mass);

    /**
     * @brief Interpolating function for the gluon-gluon fusion H cross section ratio of the top-loop and the total contribution at 13 TeV.
     * @return @f$\sigma_t(gg\to H)/\sigma(gg\to H)@f$
     */
    double ip_csr_ggH_t_13(double mass);

    /**
     * @brief Interpolating function for the gluon-gluon fusion H cross section ratio of the bottom-loop and the total contribution at 8 TeV.
     * @return @f$\sigma_b(gg\to H)/\sigma(gg\to H)@f$
     */
    double ip_csr_ggH_b_8(double mass);

    /**
     * @brief Interpolating function for the gluon-gluon fusion H cross section ratio of the bottom-loop and the total contribution at 13 TeV.
     * @return @f$\sigma_b(gg\to H)/\sigma(gg\to H)@f$
     */
    double ip_csr_ggH_b_13(double mass);

    /**
     * @brief Interpolating function for the gluon-gluon fusion A cross section ratio of the top-loop and the total contribution at 8 TeV.
     * @return @f$\sigma_t(gg\to A)/\sigma(gg\to A)@f$
     */
    double ip_csr_ggA_t_8(double mass);

    /**
     * @brief Interpolating function for the gluon-gluon fusion A cross section ratio of the top-loop and the total contribution at 13 TeV.
     * @return @f$\sigma_t(gg\to A)/\sigma(gg\to A)@f$
     */
    double ip_csr_ggA_t_13(double mass);

    /**
     * @brief Interpolating function for the gluon-gluon fusion A cross section ratio of the bottom-loop and the total contribution at 8 TeV.
     * @return @f$\sigma_b(gg\to A)/\sigma(gg\to A)@f$
     */
    double ip_csr_ggA_b_8(double mass);

    /**
     * @brief Interpolating function for the gluon-gluon fusion A cross section ratio of the bottom-loop and the total contribution at 13 TeV.
     * @return @f$\sigma_b(gg\to A)/\sigma(gg\to A)@f$
     */
    double ip_csr_ggA_b_13(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a scalar resonance decaying to two photons.
     * @return @f$[\sigma_{pp\to \phi}\cdot BR(\phi\to \gamma \gamma)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1407.6583, Figure 4 @cite Aad:2014ioa.
     */
    double ip_ex_pp_phi_gaga_ATLAS8(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a scalar resonance decaying to two photons.
     * @return @f$[\sigma_{pp\to \phi}\cdot BR(\phi\to \gamma \gamma)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1407.6583, Figure 4 @cite Aad:2014ioa.
     */
    double ip_ex_pp_phi_gaga_ATLAS8_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a pseudoscalar resonance decaying to a Z boson and a photon which further decay into two leptons and a photon.
     * @return @f$[\sigma_{pp\to \phi}\cdot BR(\phi\to Z\gamma \to \ell \ell \gamma)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1407.8150, Figure 3c @cite Aad:2014fha.
     */
    double ip_ex_pp_phi_Zga_llga_ATLAS8(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a pseudoscalar resonance decaying to a Z boson and a photon which further decay into two leptons and a photon.
     * @return @f$[\sigma_{pp\to \phi}\cdot BR(\phi\to Z\gamma \to \ell \ell \gamma)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1407.8150, Figure 3c @cite Aad:2014fha.
     */
    double ip_ex_pp_phi_Zga_llga_ATLAS8_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two tau leptons.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to \tau \tau)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1409.6064, Figure 11a @cite Aad:2014vgg.
     */
    double ip_ex_gg_phi_tautau_ATLAS8(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two tau leptons.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to \tau \tau)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1409.6064, Figure 11a @cite Aad:2014vgg.
     */
    double ip_ex_gg_phi_tautau_ATLAS8_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a bottom quark produced scalar resonance decaying to two tau leptons.
     * @return @f$[\sigma_{bb\to \phi}\cdot BR(\phi\to \tau \tau)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1409.6064, Figure 11b @cite Aad:2014vgg.
     */
    double ip_ex_bb_phi_tautau_ATLAS8(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a bottom quark produced scalar resonance decaying to two tau leptons.
     * @return @f$[\sigma_{bb\to \phi}\cdot BR(\phi\to \tau \tau)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1409.6064, Figure 11b @cite Aad:2014vgg.
     */
    double ip_ex_bb_phi_tautau_ATLAS8_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a gluon-gluon produced pseudoscalar resonance decaying to @f$hZ@f$ of which the Higgs further decays to a @f$\tau@f$ lepton pair.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to hZ\to \tau \tau Z)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1502.04478, Figure 3a @cite Aad:2015wra.
     */
    double ip_ex_gg_A_hZ_tautauZ_ATLAS8(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a gluon-gluon produced pseudoscalar resonance decaying to @f$hZ@f$ of which the Higgs further decays to a @f$\tau@f$ lepton pair.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to hZ\to \tau \tau Z)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1502.04478, Figure 3a @cite Aad:2015wra.
     */
    double ip_ex_gg_A_hZ_tautauZ_ATLAS8_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a gluon-gluon produced pseudoscalar resonance decaying to @f$hZ@f$ of which the Higgs further decays to a bottom quark pair.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to hZ\to b\bar b Z)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1502.04478, Figure 3b @cite Aad:2015wra.
     */
    double ip_ex_gg_A_hZ_bbZ_ATLAS8(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a gluon-gluon produced pseudoscalar resonance decaying to @f$hZ@f$ of which the Higgs further decays to a bottom quark pair.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to hZ\to b\bar b Z)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1502.04478, Figure 3b @cite Aad:2015wra.
     */
    double ip_ex_gg_A_hZ_bbZ_ATLAS8_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to a top quark pair.
     * @return @f$[\sigma_{pp\to \phi}\cdot BR(\phi\to t\bar t)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1505.07018, Figure 11d @cite Aad:2015fna.
     */
    double ip_ex_gg_phi_tt_ATLAS8(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to a top quark pair.
     * @return @f$[\sigma_{pp\to \phi}\cdot BR(\phi\to t\bar t)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1505.07018, Figure 11d @cite Aad:2015fna.
     */
    double ip_ex_gg_phi_tt_ATLAS8_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two @f$W@f$ bosons.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to WW)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1509.00389, Figure 13, left @cite Aad:2015agg.
     */
    double ip_ex_gg_H_WW_ATLAS8(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two @f$W@f$ bosons.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to WW)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1509.00389, Figure 13, left @cite Aad:2015agg.
     */
    double ip_ex_gg_H_WW_ATLAS8_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a vector boson fusion produced scalar resonance decaying to two @f$W@f$ bosons.
     * @return @f$[\sigma_{VV\to \phi}\cdot BR(\phi\to WW)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1509.00389, Figure 13, right @cite Aad:2015agg.
     */
    double ip_ex_VBF_H_WW_ATLAS8(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a vector boson fusion produced scalar resonance decaying to two @f$W@f$ bosons.
     * @return @f$[\sigma_{VV\to \phi}\cdot BR(\phi\to WW)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1509.00389, Figure 13, right @cite Aad:2015agg.
     */
    double ip_ex_VBF_H_WW_ATLAS8_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two @f$Z@f$ bosons.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to ZZ)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1507.05930, Figure 12(a) @cite Aad:2015kna.
     */
    double ip_ex_gg_H_ZZ_ATLAS8(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two @f$Z@f$ bosons.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to ZZ)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1507.05930, Figure 12(a) @cite Aad:2015kna.
     */
    double ip_ex_gg_H_ZZ_ATLAS8_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a vector boson fusion produced scalar resonance decaying to two @f$Z@f$ bosons.
     * @return @f$[\sigma_{VV\to \phi}\cdot BR(\phi\to ZZ)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1507.05930, Figure 12(b) @cite Aad:2015kna.
     */
    double ip_ex_VBF_H_ZZ_ATLAS8(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a vector boson fusion produced scalar resonance decaying to two @f$Z@f$ bosons.
     * @return @f$[\sigma_{VV\to \phi}\cdot BR(\phi\to ZZ)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1507.05930, Figure 12(b) @cite Aad:2015kna.
     */
    double ip_ex_VBF_H_ZZ_ATLAS8_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two @f$h@f$ bosons.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to hh)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1509.04670, Figure 6 @cite Aad:2015xja.
     */
    double ip_ex_gg_H_hh_ATLAS8(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two @f$h@f$ bosons.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to hh)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1509.04670, Figure 6 @cite Aad:2015xja.
     */
    double ip_ex_gg_H_hh_ATLAS8_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a scalar resonance decaying to two massive vector bosons.
     * @return @f$[\mu_H(H\to VV)]_{\text{CMS,95\%}}@f$
     * @details Taken from arXiv:1504.00936, Figure 7, bottom right @cite Khachatryan:2015cwa.
     */
    double ip_ex_mu_pp_H_VV_CMS8(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a scalar resonance decaying to two massive vector bosons.
     * @return @f$[\mu_H(H\to VV)]_{\text{CMS,95\%}}@f$
     * @details Taken from arXiv:1504.00936, Figure 7, bottom right @cite Khachatryan:2015cwa.
     */
    double ip_ex_mu_pp_H_VV_CMS8_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a gluon-gluon produced pseudoscalar resonance decaying to @f$hZ@f$ which further decay to a bottom quark pair and a light lepton pair.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to hZ\to b\bar b \ell \ell)]_{\text{CMS,95\%}}@f$
     * @details Taken from arXiv:1504.04710, Figure 3 @cite Khachatryan:2015lba.
     */
    double ip_ex_gg_A_hZ_bbll_CMS8(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a gluon-gluon produced pseudoscalar resonance decaying to @f$hZ@f$ which further decay to a bottom quark pair and a light lepton pair.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to hZ\to b\bar b \ell \ell)]_{\text{CMS,95\%}}@f$
     * @details Taken from arXiv:1504.04710, Figure 3 @cite Khachatryan:2015lba.
     */
    double ip_ex_gg_A_hZ_bbll_CMS8_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a scalar resonance decaying to two @f$h@f$ bosons.
     * @return @f$[\sigma_{pp\to H}\cdot BR(H\to hh)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-PAS-HIG-15-013, Figure 5a @cite CMS:2016zxv.
     */
    double ip_ex_pp_H_hh_CMS8(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a scalar resonance decaying to two @f$h@f$ bosons.
     * @return @f$[\sigma_{pp\to H}\cdot BR(H\to hh)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-PAS-HIG-15-013, Figure 5a @cite CMS:2016zxv.
     */
    double ip_ex_pp_H_hh_CMS8_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a scalar resonance decaying to two @f$h@f$ bosons which further decay to a photon pair and a bottom quark pair.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to hh\to \gamma \gamma b\bar b)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-PAS-HIG-13-032, Figure 8 @cite CMS:2014ipa.
     */
    double ip_ex_pp_phi_hh_gagabb_CMS8(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a scalar resonance decaying to two @f$h@f$ bosons which further decay to a photon pair and a bottom quark pair.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to hh\to \gamma \gamma b\bar b)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-PAS-HIG-13-032, Figure 8 @cite CMS:2014ipa.
     */
    double ip_ex_pp_phi_hh_gagabb_CMS8_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a scalar resonance decaying to two @f$h@f$ bosons which further decay to two bottom quark pairs.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to hh\to b\bar b b\bar b)]_{\text{CMS,95\%}}@f$
     * @details Taken from arXiv:1503.04114, Figure 5, left @cite Khachatryan:2015yea.
     */
    double ip_ex_pp_phi_hh_bbbb_CMS8(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a scalar resonance decaying to two @f$h@f$ bosons which further decay to two bottom quark pairs.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to hh\to b\bar b b\bar b)]_{\text{CMS,95\%}}@f$
     * @details Taken from arXiv:1503.04114, Figure 5, left @cite Khachatryan:2015yea.
     */
    double ip_ex_pp_phi_hh_bbbb_CMS8_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a bottom quark produced scalar resonance decaying to two bottom quarks.
     * @return @f$[\sigma_{bb\to \phi}\cdot BR(\phi\to b\bar b)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-14-017, Figure 6 @cite Khachatryan:2015tra.
     */
    double ip_ex_bb_phi_bb_CMS8(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a bottom quark produced scalar resonance decaying to two bottom quarks.
     * @return @f$[\sigma_{bb\to \phi}\cdot BR(\phi\to b\bar b)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-14-017, Figure 6 @cite Khachatryan:2015tra.
     */
    double ip_ex_bb_phi_bb_CMS8_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a gluon-gluon produced scalar resonance decaying to two tau leptons.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to \tau \tau)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-PAS-HIG-14-029, Figure 10-a @cite CMS:2015mca.
     */
    double ip_ex_gg_phi_tautau_CMS8(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a gluon-gluon produced scalar resonance decaying to two tau leptons.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to \tau \tau)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-PAS-HIG-14-029, Figure 10-a @cite CMS:2015mca.
     */
    double ip_ex_gg_phi_tautau_CMS8_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a bottom quark produced scalar resonance decaying to two tau leptons.
     * @return @f$[\sigma_{bb\to \phi}\cdot BR(\phi\to \tau \tau)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-PAS-HIG-14-029, Figure 10-b @cite CMS:2015mca.
     */
    double ip_ex_bb_phi_tautau_CMS8(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a bottom quark produced scalar resonance decaying to two tau leptons.
     * @return @f$[\sigma_{bb\to \phi}\cdot BR(\phi\to \tau \tau)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-PAS-HIG-14-029, Figure 10-b @cite CMS:2015mca.
     */
    double ip_ex_bb_phi_tautau_CMS8_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a gluon-gluon produced scalar resonance decaying to two photons.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to \gamma \gamma)]_{\text{CMS,95\%}}@f$
     * @details Taken from arXiv:1506.02301, Figure 7, left @cite Khachatryan:2015qba.
     */
    double ip_ex_gg_phi_gaga_CMS8(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a gluon-gluon produced scalar resonance decaying to two photons.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to \gamma \gamma)]_{\text{CMS,95\%}}@f$
     * @details Taken from arXiv:1506.02301, Figure 7, left @cite Khachatryan:2015qba.
     */
    double ip_ex_gg_phi_gaga_CMS8_e(double mass);

//        double ip_ex_ggF_phi_gaga_CMS_ep2(double mass);
//
//        double ip_ex_ggF_phi_gaga_CMS_em2(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a pseudoscalar resonance decaying to a Z boson and a photon which further decay into two leptons and a photon.
     * @return @f$[\sigma_{pp\to A}\cdot BR(A\to Z\gamma \to \ell \ell \gamma)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-PAS-HIG-16-014, Figure 2 @cite CMS:2016all.
     */
    double ip_ex_pp_A_Zga_llga_CMS8(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a pseudoscalar resonance decaying to a Z boson and a photon which further decay into two leptons and a photon.
     * @return @f$[\sigma_{pp\to A}\cdot BR(A\to Z\gamma \to \ell \ell \gamma)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-PAS-HIG-16-014, Figure 2 @cite CMS:2016all.
     */
    double ip_ex_pp_A_Zga_llga_CMS8_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a gluon-gluon produced scalar resonance decaying to two @f$h@f$ bosons which further decay to a bottom quark pair and a @f$\tau@f$ lepton pair.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to hh\to b\bar b \tau \tau)]_{\text{CMS,95\%}}@f$
     * @details Taken from arXiv:1510.01181, Figure 8, bottom right @cite Khachatryan:2015tha.
     */
    double ip_ex_gg_H_hh_bbtautau_CMS8(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a gluon-gluon produced scalar resonance decaying to two @f$h@f$ bosons which further decay to a bottom quark pair and a @f$\tau@f$ lepton pair.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to hh\to b\bar b \tau \tau)]_{\text{CMS,95\%}}@f$
     * @details Taken from arXiv:1510.01181, Figure 8, bottom right @cite Khachatryan:2015tha.
     */
    double ip_ex_gg_H_hh_bbtautau_CMS8_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a gluon-gluon produced pseudoscalar resonance decaying to @f$hZ@f$ which further decay to a @f$\tau@f$ lepton pair and a light lepton pair.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to hZ\to \tau \tau \ell \ell)]_{\text{CMS,95\%}}@f$
     * @details Taken from arXiv:1510.01181, Figure 10, left @cite Khachatryan:2015tha.
     */
    double ip_ex_gg_A_hZ_tautaull_CMS8(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a gluon-gluon produced pseudoscalar resonance decaying to @f$hZ@f$ which further decay to a @f$\tau@f$ lepton pair and a light lepton pair.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to hZ\to \tau \tau \ell \ell)]_{\text{CMS,95\%}}@f$
     * @details Taken from arXiv:1510.01181, Figure 10, left @cite Khachatryan:2015tha.
     */
    double ip_ex_gg_A_hZ_tautaull_CMS8_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a bb associated scalar resonance decaying to t quarks.
     * @return @f$[\sigma_{bb\to \phi}\cdot BR(\phi\to t\bar t)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-104, Figure 21 @cite TheATLAScollaboration:2016loc.
     */
    double ip_ex_bb_phi_tt_ATLAS13(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a bb associated scalar resonance decaying to t quarks.
     * @return @f$[\sigma_{bb\to \phi}\cdot BR(\phi\to t\bar t)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-104, Figure 21 @cite TheATLAScollaboration:2016loc.
     */
    double ip_ex_bb_phi_tt_ATLAS13_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a tt associated scalar resonance decaying to t quarks.
     * @return @f$[\sigma_{tt\to \phi}\cdot BR(\phi\to t\bar t)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-104, Figure 22 @cite TheATLAScollaboration:2016loc.
     */
    double ip_ex_tt_phi_tt_ATLAS13(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a tt associated scalar resonance decaying to t quarks.
     * @return @f$[\sigma_{tt\to \phi}\cdot BR(\phi\to t\bar t)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-104, Figure 22 @cite TheATLAScollaboration:2016loc.
     */
    double ip_ex_tt_phi_tt_ATLAS13_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two tau leptons.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to \tau \tau)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-085, Table 5 @cite ATLAS:2016fpj.
     */
    double ip_ex_gg_phi_tautau_ATLAS13(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two tau leptons.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to \tau \tau)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-085, Table 5 @cite ATLAS:2016fpj.
     */
    double ip_ex_gg_phi_tautau_ATLAS13_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a bb associated scalar resonance decaying to two tau leptons.
     * @return @f$[\sigma_{bb\to \phi}\cdot BR(\phi\to \tau \tau)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-085, Table 5 @cite ATLAS:2016fpj.
     */
    double ip_ex_bb_phi_tautau_ATLAS13(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a bb associated scalar resonance decaying to two tau leptons.
     * @return @f$[\sigma_{bb\to \phi}\cdot BR(\phi\to \tau \tau)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-085, Table 5 @cite ATLAS:2016fpj.
     */
    double ip_ex_bb_phi_tautau_ATLAS13_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a scalar resonance decaying to two photons.
     * @return @f$[\sigma_{pp\to \phi}\cdot BR(\phi\to \gamma \gamma)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-059, Figure 7a @cite ATLAS:2016eeo.
     */
    double ip_ex_pp_phi_gaga_ATLAS13(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a scalar resonance decaying to two photons.
     * @return @f$[\sigma_{pp\to \phi}\cdot BR(\phi\to \gamma \gamma)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-059, Figure 7a @cite ATLAS:2016eeo.
     */
    double ip_ex_pp_phi_gaga_ATLAS13_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a scalar resonance decaying to a @f$Z@f$ boson and a photon, of which the @f$Z@f$ further decays to two leptons.
     * @return @f$[\sigma_{pp\to \phi}\cdot BR(\phi\to Z \gamma)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-044, Figure 6 @cite ATLAS:2016lri.
     */
    double ip_ex_pp_phi_Zga_llga_ATLAS13(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a scalar resonance decaying to a @f$Z@f$ boson and a photon, of which the @f$Z@f$ further decays to two leptons.
     * @return @f$[\sigma_{pp\to \phi}\cdot BR(\phi\to Z \gamma)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-044, Figure 6 @cite ATLAS:2016lri.
     */
    double ip_ex_pp_phi_Zga_llga_ATLAS13_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two @f$Z@f$ bosons which further decay to two leptons and two neutrinos.
     * @return @f$[\sigma_{pp\to H}\cdot BR(H\to ZZ)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-056, Figure 7a @cite ATLAS:2016bza.
     */
    double ip_ex_gg_H_ZZ_llnunu_ATLAS13(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two @f$Z@f$ bosons which further decay to two leptons and two neutrinos.
     * @return @f$[\sigma_{pp\to H}\cdot BR(H\to ZZ)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-056, Figure 7a @cite ATLAS:2016bza.
     */
    double ip_ex_gg_H_ZZ_llnunu_ATLAS13_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two @f$Z@f$ bosons which further decay to four leptons.
     * @return @f$[\sigma_{gg\to H}\cdot BR(H\to ZZ \to \ell \ell \ell \ell)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-079, Figure 11a @cite ATLAS:2016oum.
     */
    double ip_ex_gg_H_ZZ_llll_ATLAS13(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two @f$Z@f$ bosons which further decay to four leptons.
     * @return @f$[\sigma_{gg\to H}\cdot BR(H\to ZZ \to \ell \ell \ell \ell)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-079, Figure 11a @cite ATLAS:2016oum.
     */
    double ip_ex_gg_H_ZZ_llll_ATLAS13_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a vector boson fusion produced scalar resonance decaying to two @f$Z@f$ bosons which further decay to four leptons.
     * @return @f$[\sigma_{VV\to H}\cdot BR(H\to ZZ \to \ell \ell \ell \ell)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-079, Figure 11b @cite ATLAS:2016oum.
     */
    double ip_ex_VBF_H_ZZ_llll_ATLAS13(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a vector boson fusion produced scalar resonance decaying to two @f$Z@f$ bosons which further decay to four leptons.
     * @return @f$[\sigma_{VV\to H}\cdot BR(H\to ZZ \to \ell \ell \ell \ell)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-079, Figure 11b @cite ATLAS:2016oum.
     */
    double ip_ex_VBF_H_ZZ_llll_ATLAS13_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two @f$Z@f$ bosons which further decay to two leptons and two quarks.
     * @return @f$[\sigma_{gg\to H}\cdot BR(H\to ZZ)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-082, Figure 10a @cite ATLAS:2016npe.
     */
    double ip_ex_gg_H_ZZ_llqq_ATLAS13(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two @f$Z@f$ bosons which further decay to two leptons and two quarks.
     * @return @f$[\sigma_{gg\to H}\cdot BR(H\to ZZ)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-082, Figure 10a @cite ATLAS:2016npe.
     */
    double ip_ex_gg_H_ZZ_llqq_ATLAS13_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a vector boson fusion produced scalar resonance decaying to two @f$Z@f$ bosons which further decay to two leptons and two quarks.
     * @return @f$[\sigma_{VV\to H}\cdot BR(H\to ZZ)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-082, Figure 10b @cite ATLAS:2016npe.
     */
    double ip_ex_VBF_H_ZZ_llqq_ATLAS13(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a vector boson fusion produced scalar resonance decaying to two @f$Z@f$ bosons which further decay to two leptons and two quarks.
     * @return @f$[\sigma_{VV\to H}\cdot BR(H\to ZZ)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-082, Figure 10b @cite ATLAS:2016npe.
     */
    double ip_ex_VBF_H_ZZ_llqq_ATLAS13_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two @f$Z@f$ bosons which further decay to two neutrinos and two quarks.
     * @return @f$[\sigma_{gg\to H}\cdot BR(H\to ZZ)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-082, Figure 12c @cite ATLAS:2016npe.
     */
    double ip_ex_gg_H_ZZ_nunuqq_ATLAS13(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two @f$Z@f$ bosons which further decay to two neutrinos and two quarks.
     * @return @f$[\sigma_{gg\to H}\cdot BR(H\to ZZ)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-082, Figure 12c @cite ATLAS:2016npe.
     */
    double ip_ex_gg_H_ZZ_nunuqq_ATLAS13_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two @f$W@f$ bosons which further decay to two lepton-neutrino pairs.
     * @return @f$[\sigma_{gg\to H}\cdot BR(H\to WW)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-074, Figure 4a @cite ATLAS:2016kjy.
     */
    double ip_ex_gg_H_WW_enumunu_ATLAS13(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two @f$W@f$ bosons which further decay to two lepton-neutrino pairs.
     * @return @f$[\sigma_{gg\to H}\cdot BR(H\to WW)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-074, Figure 4a @cite ATLAS:2016kjy.
     */
    double ip_ex_gg_H_WW_enumunu_ATLAS13_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a vector boson fusion produced scalar resonance decaying to two @f$W@f$ bosons which further decay to two lepton-neutrino pairs.
     * @return @f$[\sigma_{VV\to H}\cdot BR(H\to WW)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-074, Figure 4b @cite ATLAS:2016kjy.
     */
    double ip_ex_VBF_H_WW_enumunu_ATLAS13(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a vector boson fusion produced scalar resonance decaying to two @f$W@f$ bosons which further decay to two lepton-neutrino pairs.
     * @return @f$[\sigma_{VV\to H}\cdot BR(H\to WW)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-074, Figure 4b @cite ATLAS:2016kjy.
     */
    double ip_ex_VBF_H_WW_enumunu_ATLAS13_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two @f$W@f$ bosons which further decay to a lepton-neutrino pair and a pair of quarks.
     * @return @f$[\sigma_{gg\to H}\cdot BR(H\to WW)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-062, Figure 5 @cite ATLAS:2016cwq.
     */
    double ip_ex_gg_H_WW_lnuqq_ATLAS13(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two @f$W@f$ bosons which further decay to a lepton-neutrino pair and a pair of quarks.
     * @return @f$[\sigma_{gg\to H}\cdot BR(H\to WW)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-062, Figure 5 @cite ATLAS:2016cwq.
     */
    double ip_ex_gg_H_WW_lnuqq_ATLAS13_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a spin-2 resonance decaying to two @f$h@f$ bosons which further decay to four b quarks.
     * @return @f$[\sigma_{pp\to H}\cdot BR(H\to hh\to b\bar b b\bar b)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-049, Figure 11 @cite ATLAS:2016ixk.
     */
    double ip_ex_pp_H_hh_bbbb_ATLAS13(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a spin-2 resonance decaying to two @f$h@f$ bosons which further decay to four b quarks.
     * @return @f$[\sigma_{pp\to H}\cdot BR(H\to hh\to b\bar b b\bar b)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-049, Figure 11 @cite ATLAS:2016ixk.
     */
    double ip_ex_pp_H_hh_bbbb_ATLAS13_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a scalar resonance decaying to two @f$h@f$ bosons which further decay to two photons and a b quark pair.
     * @return @f$[\sigma_{pp\to H}\cdot BR(H\to hh)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-004, Figure 6a @cite TheATLAScollaboration:2016ibb.
     */
    double ip_ex_pp_H_hh_gagabb_ATLAS13(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a scalar resonance decaying to two @f$h@f$ bosons which further decay to two photons and a b quark pair.
     * @return @f$[\sigma_{pp\to H}\cdot BR(H\to hh)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-004, Figure 6a @cite TheATLAScollaboration:2016ibb.
     */
    double ip_ex_pp_H_hh_gagabb_ATLAS13_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two @f$h@f$ bosons which further decay to two photons and two @f$W@f$ bosons.
     * @return @f$[\sigma_{gg\to H}\cdot BR(H\to hh)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-071, Figure 3a @cite ATLAS:2016qmt.
     */
    double ip_ex_pp_H_hh_gagaWW_ATLAS13(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two @f$h@f$ bosons which further decay to two photons and two @f$W@f$ bosons.
     * @return @f$[\sigma_{gg\to H}\cdot BR(H\to hh)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-071, Figure 3a @cite ATLAS:2016qmt.
     */
    double ip_ex_pp_H_hh_gagaWW_ATLAS13_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a gluon-gluon produced pseudoscalar resonance decaying to a @f$Z@f$ and an @f$h@f$ boson, of which the latter further decays to a b quark pair.
     * @return @f$[\sigma_{gg\to A}\cdot BR(A\to Zh\to Zb\bar b)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-015, Figure 8a @cite TheATLAScollaboration:2016loc.
     */
    double ip_ex_gg_A_Zh_Zbb_ATLAS13(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a gluon-gluon produced pseudoscalar resonance decaying to a @f$Z@f$ and an @f$h@f$ boson, of which the latter further decays to a b quark pair.
     * @return @f$[\sigma_{gg\to A}\cdot BR(A\to Zh\to Zb\bar b)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-015, Figure 8a @cite TheATLAScollaboration:2016loc.
     */
    double ip_ex_gg_A_Zh_Zbb_ATLAS13_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a bb associated pseudoscalar resonance decaying to a @f$Z@f$ and an @f$h@f$ boson, of which the latter further decays to a b quark pair.
     * @return @f$[\sigma_{bb\to A}\cdot BR(A\to Zh\to Zb\bar b)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-015, Figure 8b @cite TheATLAScollaboration:2016loc.
     */
    double ip_ex_bb_A_Zh_Zbb_ATLAS13(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a bb associated pseudoscalar resonance decaying to a @f$Z@f$ and an @f$h@f$ boson, of which the latter further decays to a b quark pair.
     * @return @f$[\sigma_{bb\to A}\cdot BR(A\to Zh\to Zb\bar b)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-015, Figure 8b @cite TheATLAScollaboration:2016loc.
     */
    double ip_ex_bb_A_Zh_Zbb_ATLAS13_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a scalar resonance decaying to a b quark pair.
     * @return @f$[\sigma_{pp\to \phi}\cdot BR(\phi \to b\bar b)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-PAS-16-025, Figure 5 @cite CMS:2016ncz.
     */
    double ip_ex_pp_phi_bb_CMS13(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a scalar resonance decaying to a b quark pair.
     * @return @f$[\sigma_{pp\to \phi}\cdot BR(\phi \to b\bar b)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-PAS-16-025, Figure 5 @cite CMS:2016ncz.
     */
    double ip_ex_pp_phi_bb_CMS13_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a gluon-gluon produced scalar resonance decaying to two tau leptons.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to \tau \tau)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-PAS-16-037, Figure 9-a @cite CMS:2016rjp.
     */
    double ip_ex_gg_phi_tautau_CMS13(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a gluon-gluon produced scalar resonance decaying to two tau leptons.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to \tau \tau)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-PAS-16-037, Figure 9-a @cite CMS:2016rjp.
     */
    double ip_ex_gg_phi_tautau_CMS13_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a b associated scalar resonance decaying to two tau leptons.
     * @return @f$[\sigma_{bb\to \phi}\cdot BR(\phi\to \tau \tau)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-PAS-16-037, Figure 9-b @cite CMS:2016rjp.
     */
    double ip_ex_bb_phi_tautau_CMS13(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a b associated scalar resonance decaying to two tau leptons.
     * @return @f$[\sigma_{bb\to \phi}\cdot BR(\phi\to \tau \tau)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-PAS-16-037, Figure 9-b @cite CMS:2016rjp.
     */
    double ip_ex_bb_phi_tautau_CMS13_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a gluon-gluon produced scalar resonance decaying to two photons.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to \gamma \gamma)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-EXO-PAS-16-027, Figure 7, top @cite CMS:2016crm.
     */
    double ip_ex_gg_phi_gaga_CMS13(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a gluon-gluon produced scalar resonance decaying to two photons.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to \gamma \gamma)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-EXO-PAS-16-027, Figure 7, top @cite CMS:2016crm.
     */
    double ip_ex_gg_phi_gaga_CMS13_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a scalar resonance decaying to a @f$Z@f$ boson and a photon, of which the @f$Z@f$ further decays to two leptons.
     * @return @f$[\sigma_{pp\to \phi}\cdot BR(\phi\to Z \gamma)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-EXO-PAS-16-034, Figure 4 @cite CMS:2016pax.
     */
    double ip_ex_pp_phi_Zga_llga_CMS13(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a scalar resonance decaying to a @f$Z@f$ boson and a photon, of which the @f$Z@f$ further decays to two leptons.
     * @return @f$[\sigma_{pp\to \phi}\cdot BR(\phi\to Z \gamma)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-EXO-PAS-16-034, Figure 4 @cite CMS:2016pax.
     */
    double ip_ex_pp_phi_Zga_llga_CMS13_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a scalar resonance decaying to a @f$Z@f$ boson and a photon, of which the @f$Z@f$ further decays to two quarks.
     * @return @f$[\sigma_{pp\to \phi}\cdot BR(\phi\to Z \gamma)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-EXO-PAS-16-035, Figure 5-a @cite CMS:2016cbb.
     */
    double ip_ex_pp_phi_Zga_qqga_CMS13(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a scalar resonance decaying to a @f$Z@f$ boson and a photon, of which the @f$Z@f$ further decays to two quarks.
     * @return @f$[\sigma_{pp\to \phi}\cdot BR(\phi\to Z \gamma)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-EXO-PAS-16-035, Figure 5-a @cite CMS:2016cbb.
     */
    double ip_ex_pp_phi_Zga_qqga_CMS13_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a gluon-gluon produced scalar resonance decaying to two @f$Z@f$ bosons which further decay to four leptons.
     * @return @f$[\sigma_{gg\to H}\cdot BR(H\to ZZ\to \ell \ell \ell \ell)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-PAS-16-033, Figure 16-a @cite CMS:2016ilx.
     */
    double ip_ex_pp_H_ZZ_llll_CMS13(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a gluon-gluon produced scalar resonance decaying to two @f$Z@f$ bosons which further decay to four leptons.
     * @return @f$[\sigma_{gg\to H}\cdot BR(H\to ZZ\to \ell \ell \ell \ell)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-PAS-16-033, Figure 16-a @cite CMS:2016ilx.
     */
    double ip_ex_pp_H_ZZ_llll_CMS13_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a vector boson fusion produced scalar resonance decaying to two @f$Z@f$ bosons which further decay to four leptons.
     * @return @f$[\sigma_{VV\to H}\cdot BR(H\to ZZ\to \ell \ell \ell \ell)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-PAS-16-033, Figure 16-b @cite CMS:2016ilx.
     */
    double ip_ex_VBF_VH_H_ZZ_llll_CMS13(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a vector boson fusion produced scalar resonance decaying to two @f$Z@f$ bosons which further decay to four leptons.
     * @return @f$[\sigma_{VV\to H}\cdot BR(H\to ZZ\to \ell \ell \ell \ell)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-PAS-16-033, Figure 16-b @cite CMS:2016ilx.
     */
    double ip_ex_VBF_VH_H_ZZ_llll_CMS13_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a scalar resonance decaying to two @f$Z@f$ bosons which further decay to two leptons and two quarks.
     * @return @f$[\sigma_{pp\to H}\cdot BR(H\to ZZ)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-PAS-16-034, Figure 4-a @cite CMS:2017sbi.
     */
    double ip_ex_pp_H_ZZ_llqq_CMS13(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a scalar resonance decaying to two @f$Z@f$ bosons which further decay to two leptons and two quarks.
     * @return @f$[\sigma_{pp\to H}\cdot BR(H\to ZZ)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-PAS-16-034, Figure 4-a @cite CMS:2017sbi.
     */
    double ip_ex_pp_H_ZZ_llqq_CMS13_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a gluon-gluon or vector boson fusion produced scalar resonance decaying to two @f$W@f$ bosons which further decay to two lepton-neutrino pairs.
     * @return @f$[\sigma_{gg/VV\to H}\cdot BR(H\to WW\to \ell \nu \ell \nu)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-PAS-16-023, Figure 5-a @cite CMS:2016jpd.
     */
    double ip_ex_ggVV_H_WW_lnulnu_CMS13(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a gluon-gluon or vector boson fusion produced scalar resonance decaying to two @f$W@f$ bosons which further decay to two lepton-neutrino pairs.
     * @return @f$[\sigma_{gg/VV\to H}\cdot BR(H\to WW\to \ell \nu \ell \nu)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-PAS-16-023, Figure 5-a @cite CMS:2016jpd.
     */
    double ip_ex_ggVV_H_WW_lnulnu_CMS13_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a scalar resonance decaying to two @f$h@f$ bosons which further decay to four b quarks.
     * @return @f$[\sigma_{pp\to H}\cdot BR(H\to hh\to b\bar b b\bar b)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-PAS-16-002, Figure 7 @cite CMS:2016tlj.
     */
    double ip_ex_pp_H_hh_bbbb_CMS13(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a scalar resonance decaying to two @f$h@f$ bosons which further decay to four b quarks.
     * @return @f$[\sigma_{pp\to H}\cdot BR(H\to hh\to b\bar b b\bar b)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-PAS-16-002, Figure 7 @cite CMS:2016tlj.
     */
    double ip_ex_pp_H_hh_bbbb_CMS13_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a scalar resonance decaying to two @f$h@f$ bosons which further decay to two photons and a b quark pair.
     * @return @f$[\sigma_{pp\to H}\cdot BR(H\to hh\to \gamma \gamma b\bar b)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-PAS-16-032, Figure 14-a @cite CMS:2016vpz.
     */
    double ip_ex_pp_H_hh_gagabb_CMS13(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a scalar resonance decaying to two @f$h@f$ bosons which further decay to two photons and a b quark pair.
     * @return @f$[\sigma_{pp\to H}\cdot BR(H\to hh\to \gamma \gamma b\bar b)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-PAS-16-032, Figure 14-a @cite CMS:2016vpz.
     */
    double ip_ex_pp_H_hh_gagabb_CMS13_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a scalar resonance decaying to two @f$h@f$ bosons which further decay to a b quark pair and two tau leptons.
     * @return @f$[\sigma_{pp\to H}\cdot BR(H\to hh\to b\bar b \tau \tau)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-PAS-16-029, Figure 3 @cite CMS:2016knm.
     */
    double ip_ex_pp_H_hh_bbtautau_CMS13(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a scalar resonance decaying to two @f$h@f$ bosons which further decay to a b quark pair and two tau leptons.
     * @return @f$[\sigma_{pp\to H}\cdot BR(H\to hh\to b\bar b \tau \tau)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-PAS-16-029, Figure 3 @cite CMS:2016knm.
     */
    double ip_ex_pp_H_hh_bbtautau_CMS13_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a scalar resonance decaying to two @f$h@f$ bosons which further decay to a b quark pair and two tau leptons.
     * @return @f$[\sigma_{pp\to H}\cdot BR(H\to hh\to b\bar b \tau \tau)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-PAS-17-002, Figure 5-a @cite CMS:2017orf.
     */
    double ip_ex_pp_H_hh_bbtautau1_CMS13(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a scalar resonance decaying to two @f$h@f$ bosons which further decay to a b quark pair and two tau leptons.
     * @return @f$[\sigma_{pp\to H}\cdot BR(H\to hh\to b\bar b \tau \tau)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-PAS-17-002, Figure 5-a @cite CMS:2017orf.
     */
    double ip_ex_pp_H_hh_bbtautau1_CMS13_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a scalar resonance decaying to two @f$h@f$ bosons which further decay to a b quark pair and two lepton-neutrino pairs.
     * @return @f$[\sigma_{pp\to H}\cdot BR(H\to hh\to b\bar b \ell \nu \ell \nu)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-PAS-16-011, Figure 5-a @cite CMS:2016rec.
     */
    double ip_ex_pp_H_hh_bblnulnu_CMS13(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a scalar resonance decaying to two @f$h@f$ bosons which further decay to a b quark pair and two lepton-neutrino pairs.
     * @return @f$[\sigma_{pp\to H}\cdot BR(H\to hh\to b\bar b \ell \nu \ell \nu)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-PAS-16-011, Figure 5-a @cite CMS:2016rec.
     */
    double ip_ex_pp_H_hh_bblnulnu_CMS13_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a scalar resonance decaying to two @f$h@f$ bosons which further decay to a b quark pair, two leptons and two neutrinos.
     * @return @f$[\sigma_{pp\to H}\cdot BR(H\to hh\to b\bar b \ell \nu \ell \nu)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-PAS-17-006, Figure 6-a @cite CMS:2017ums.
     */
    double ip_ex_pp_H_hh_bbVV_CMS13(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a scalar resonance decaying to two @f$h@f$ bosons which further decay to a b quark pair, two leptons and two neutrinos.
     * @return @f$[\sigma_{pp\to H}\cdot BR(H\to hh\to b\bar b \ell \nu \ell \nu)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-PAS-17-006, Figure 6-a @cite CMS:2017ums.
     */
    double ip_ex_pp_H_hh_bbVV_CMS13_e(double mass);

    /*Empty interpolation templates*/
    double ie1(double mass);
    double ie1e(double mass);
    double ie2(double mass);
    double ie2e(double mass);
    double ie3(double mass);
    double ie3e(double mass);
    double ie4(double mass);
    double ie4e(double mass);
    double ie5(double mass);
    double ie5e(double mass);
    double ie6(double mass);
    double ie6e(double mass);
    double ie7(double mass);
    double ie7e(double mass);
    double ie8(double mass);
    double ie8e(double mass);
    double ie9(double mass);
    double ie9e(double mass);
    double ie10(double mass);
    double ie10e(double mass);
    double ie11(double mass);
    double ie11e(double mass);
    double ie12(double mass);
    double ie12e(double mass);
    double ie13(double mass);
    double ie13e(double mass);
    double ie14(double mass);
    double ie14e(double mass);
    double ie15(double mass);
    double ie15e(double mass);
    double ie16(double mass);
    double ie16e(double mass);
    double ie17(double mass);
    double ie17e(double mass);
    double ie18(double mass);
    double ie18e(double mass);
    double ie19(double mass);
    double ie19e(double mass);
    double ie20(double mass);
    double ie20e(double mass);
    double ie21(double mass);
    double ie21e(double mass);
    double ie22(double mass);
    double ie22e(double mass);
    double ie23(double mass);
    double ie23e(double mass);
    double ie24(double mass);
    double ie24e(double mass);
    double ie25(double mass);
    double ie25e(double mass);
    double ie26(double mass);
    double ie26e(double mass);
    double ie27(double mass);
    double ie27e(double mass);
    double ie28(double mass);
    double ie28e(double mass);
    double ie29(double mass);
    double ie29e(double mass);
    double ie30(double mass);
    double ie30e(double mass);
    double ie31(double mass);
    double ie31e(double mass);
    double ie32(double mass);
    double ie32e(double mass);
    double ie33(double mass);
    double ie33e(double mass);
    double ie34(double mass);
    double ie34e(double mass);
    double ie35(double mass);
    double ie35e(double mass);
    double ie36(double mass);
    double ie36e(double mass);
    double ie37(double mass);
    double ie37e(double mass);
    double ie38(double mass);
    double ie38e(double mass);
    double ie39(double mass);
    double ie39e(double mass);
    double ie40(double mass);
    double ie40e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a singly charged scalar resonance decaying to a @f$\tau@f$ lepton and a neutrino.
     * @return @f$[\sigma_{pp\to H^\pm}\cdot BR(H^\pm\to \tau \nu)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1412.6663, Figure 7-b @cite Aad:2014kga.
     */
    double ip_ex_pp_Hpm_taunu_ATLAS8(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a singly charged scalar resonance decaying to a @f$\tau@f$ lepton and a neutrino.
     * @return @f$[\sigma_{pp\to H^\pm}\cdot BR(H^\pm\to \tau \nu)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1412.6663, Figure 7-b @cite Aad:2014kga.
     */
    double ip_ex_pp_Hpm_taunu_ATLAS8_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a singly charged scalar resonance decaying to a @f$t@f$ quark and a @f$b@f$ quark.
     * @return @f$[\sigma_{pp\to H^+}\cdot BR(H^+\to tb)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1512.03704, Figure 6 @cite Aad:2015typ.
     */
    double ip_ex_pp_Hp_tb_ATLAS8(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a singly charged scalar resonance decaying to a @f$t@f$ quark and a @f$b@f$ quark.
     * @return @f$[\sigma_{pp\to H^+}\cdot BR(H^+\to tb)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1512.03704, Figure 6 @cite Aad:2015typ.
     */
    double ip_ex_pp_Hp_tb_ATLAS8_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a singly charged scalar resonance decaying to a @f$\tau@f$ lepton and a neutrino.
     * @return @f$[\sigma_{pp\to H^+}]_{\text{CMS,95\%}}@f$
     * @details Taken from arXiv:1508.07774, Table 10 bottom @cite Khachatryan:2015qxa.
     */
    double ip_ex_pp_Hp_taunu_CMS8(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a singly charged scalar resonance decaying to a @f$\tau@f$ lepton and a neutrino.
     * @return @f$[\sigma_{pp\to H^+}]_{\text{CMS,95\%}}@f$
     * @details Taken from arXiv:1508.07774, Table 10 bottom @cite Khachatryan:2015qxa.
     */
    double ip_ex_pp_Hp_taunu_CMS8_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a singly charged scalar resonance decaying to a @f$t@f$ quark and a @f$b@f$ quark.
     * @return @f$[\sigma_{pp\to H^+}]_{\text{CMS,95\%}}@f$
     * @details Taken from arXiv:1508.07774, Table 11 @cite Khachatryan:2015qxa.
     */
    double ip_ex_pp_Hp_tb_CMS8(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a singly charged scalar resonance decaying to a @f$t@f$ quark and a @f$b@f$ quark.
     * @return @f$[\sigma_{pp\to H^+}]_{\text{CMS,95\%}}@f$
     * @details Taken from arXiv:1508.07774, Table 11 @cite Khachatryan:2015qxa.
     */
    double ip_ex_pp_Hp_tb_CMS8_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a singly charged scalar resonance decaying to a @f$\tau@f$ lepton and a neutrino.
     * @return @f$[\sigma_{pp\to H^\pm}\cdot BR(H^\pm\to \tau \nu)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-088, Figure 6 @cite ATLAS:2016grc.
     */
    double ip_ex_pp_Hpm_taunu_ATLAS13(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a singly charged scalar resonance decaying to a @f$\tau@f$ lepton and a neutrino.
     * @return @f$[\sigma_{pp\to H^\pm}\cdot BR(H^\pm\to \tau \nu)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-088, Figure 6 @cite ATLAS:2016grc.
     */
    double ip_ex_pp_Hpm_taunu_ATLAS13_e(double mass);

    /**
     * @brief Interpolating function for the observed CMS upper limit on a singly charged scalar resonance decaying to a @f$\tau@f$ lepton and a neutrino.
     * @return @f$[\sigma_{pp\to H^\pm}\cdot BR(H^\pm\to \tau \nu)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-PAS-HIG-16-031, Figure 6 right @cite CMS:2016szv.
     */
    double ip_ex_pp_Hpm_taunu_CMS13(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a singly charged scalar resonance decaying to a @f$\tau@f$ lepton and a neutrino.
     * @return @f$[\sigma_{pp\to H^\pm}\cdot BR(H^\pm\to \tau \nu)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-PAS-HIG-16-031, Figure 6 right @cite CMS:2016szv.
     */
    double ip_ex_pp_Hpm_taunu_CMS13_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a singly charged scalar resonance decaying to a @f$t@f$ quark and a @f$b@f$ quark.
     * @return @f$[\sigma_{pp\to H^+}\cdot BR(H^+\to tb)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-089, Figure 12 @cite ATLAS:2016qiq.
     */
    double ip_ex_pp_Hp_tb_ATLAS13_1(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a singly charged scalar resonance decaying to a @f$t@f$ quark and a @f$b@f$ quark.
     * @return @f$[\sigma_{pp\to H^+}\cdot BR(H^+\to tb)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-089, Figure 12 @cite ATLAS:2016qiq.
     */
    double ip_ex_pp_Hp_tb_ATLAS13_1_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a singly charged scalar resonance decaying to a @f$t@f$ quark and a @f$b@f$ quark.
     * @return @f$[\sigma_{pp\to H^+}\cdot BR(H^+\to tb)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-104, Figure 23 @cite ATLAS:2016btu.
     */
    double ip_ex_pp_Hp_tb_ATLAS13_2(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a singly charged scalar resonance decaying to a @f$t@f$ quark and a @f$b@f$ quark.
     * @return @f$[\sigma_{pp\to H^+}\cdot BR(H^+\to tb)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from ATLAS-CONF-2016-104, Figure 23 @cite ATLAS:2016btu.
     */
    double ip_ex_pp_Hp_tb_ATLAS13_2_e(double mass);

    /**
     * @brief Interpolating function for the NNLO value for the branching ratio of @f$b\to s \gamma@f$ decays in the THDM.
     * @return @f$BR(B\to X_s \gamma)@f$
     * @details Values derived with the help of the authors of @cite Misiak:2015xwa.
     */
    double ip_ex_bsgamma(double logtb, double logmHp);

    /*One-loop functions*/

    gslpp::complex B0_MZ2_0_MW2_mHh2(const double MZ2, const double MW2, const double mHh2) const;
    gslpp::complex B0_MZ2_0_MW2_mHl2(const double MZ2, const double MW2, const double mHl2) const;
    gslpp::complex B0_MZ2_0_MZ2_mHh2(const double MZ2, const double mHh2) const;
    gslpp::complex B0_MZ2_0_MZ2_mHl2(const double MZ2, const double mHl2) const;
    gslpp::complex B0_MZ2_MW2_MW2_mHh2(const double MZ2, const double MW2, const double mHh2) const;
    gslpp::complex B0_MZ2_MW2_MW2_mHl2(const double MZ2, const double MW2, const double mHl2) const;
    gslpp::complex B0_MZ2_MZ2_MZ2_mHh2(const double MZ2, const double mHh2) const;
    gslpp::complex B0_MZ2_MZ2_MZ2_mHl2(const double MZ2, const double mHl2) const;

    gslpp::complex B0_MZ2_0_0_mHl2(const double MZ2, const double mHl2) const;
    gslpp::complex B0_MZ2_0_0_mHh2(const double MZ2, const double mHh2) const;
    gslpp::complex B0_MZ2_0_mHp2_mHl2(const double MZ2, const double mHp2, const double mHl2) const;
    gslpp::complex B0_MZ2_0_mHp2_mHh2(const double MZ2, const double mHp2, const double mHh2) const;
    gslpp::complex B0_MZ2_0_mA2_mHl2(const double MZ2, const double mA2, const double mHl2) const;
    gslpp::complex B0_MZ2_0_mA2_mHh2(const double MZ2, const double mA2, const double mHh2) const;
    gslpp::complex B0_MZ2_mHl2_0_0(const double MZ2, const double mHl2) const;
    gslpp::complex B0_MZ2_mHl2_0_mHp2(const double MZ2, const double mHl2, const double mHp2) const;
    gslpp::complex B0_MZ2_mHl2_0_mA2(const double MZ2, const double mHl2, const double mA2) const;
    gslpp::complex B0_MZ2_mHl2_mHl2_mHl2(const double MZ2, const double mHl2) const;
    gslpp::complex B0_MZ2_mHl2_mHh2_mHl2(const double MZ2, const double mHl2, const double mHh2) const;
    gslpp::complex B0_MZ2_mHl2_mHh2_mHh2(const double MZ2, const double mHl2, const double mHh2) const;
    gslpp::complex B0_MZ2_mHl2_mHp2_mHp2(const double MZ2, const double mHl2, const double mHp2) const;
    gslpp::complex B0_MZ2_mHl2_mA2_mA2(const double MZ2, const double mHl2, const double mA2) const;
    gslpp::complex B0_MZ2_mHh2_0_0(const double MZ2, const double mHh2) const;
    gslpp::complex B0_MZ2_mHh2_0_mHp2(const double MZ2, const double mHh2, const double mHp2) const;
    gslpp::complex B0_MZ2_mHh2_0_mA2(const double MZ2, const double mHh2, const double mA2) const;
    gslpp::complex B0_MZ2_mHh2_mHl2_mHl2(const double MZ2, const double mHh2, const double mHl2) const;
    gslpp::complex B0_MZ2_mHh2_mHh2_mHl2(const double MZ2, const double mHh2, const double mHl2) const;
    gslpp::complex B0_MZ2_mHh2_mHh2_mHh2(const double MZ2, const double mHh2) const;
    gslpp::complex B0_MZ2_mHh2_mHp2_mHp2(const double MZ2, const double mHh2, const double mHp2) const;
    gslpp::complex B0_MZ2_mHh2_mA2_mA2(const double MZ2, const double mHh2, const double mA2) const;
    gslpp::complex B0_MZ2_mHp2_0_mHl2(const double MZ2, const double mHp2, const double mHl2) const;
    gslpp::complex B0_MZ2_mHp2_0_mHh2(const double MZ2, const double mHp2, const double mHh2) const;
    gslpp::complex B0_MZ2_mHp2_mHp2_mHl2(const double MZ2, const double mHp2, const double mHl2) const;
    gslpp::complex B0_MZ2_mHp2_mHp2_mHh2(const double MZ2, const double mHp2, const double mHh2) const;
    gslpp::complex B0_MZ2_mA2_0_mHl2(const double MZ2, const double mA2, const double mHl2) const;
    gslpp::complex B0_MZ2_mA2_0_mHh2(const double MZ2, const double mA2, const double mHh2) const;
    gslpp::complex B0_MZ2_mA2_mA2_mHl2(const double MZ2, const double mA2, const double mHl2) const;
    gslpp::complex B0_MZ2_mA2_mA2_mHh2(const double MZ2, const double mA2, const double mHh2) const;

    gslpp::complex B0p_MZ2_0_0_mHl2(const double MZ2, const double mHl2) const;
    gslpp::complex B0p_MZ2_0_0_mHh2(const double MZ2, const double mHh2) const;
    gslpp::complex B0p_MZ2_0_mHp2_mHl2(const double MZ2, const double mHp2, const double mHl2) const;
    gslpp::complex B0p_MZ2_0_mHp2_mHh2(const double MZ2, const double mHp2, const double mHh2) const;
    gslpp::complex B0p_MZ2_0_mHp2_mA2(const double MZ2, const double mHp2, const double mA2) const;
    gslpp::complex B0p_MZ2_0_mA2_mHl2(const double MZ2, const double mA2, const double mHl2) const;
    gslpp::complex B0p_MZ2_0_mA2_mHh2(const double MZ2, const double mA2, const double mHh2) const;
    gslpp::complex B0p_MZ2_mHl2_0_0(const double MZ2, const double mHl2) const;
    gslpp::complex B0p_MZ2_mHl2_0_mHp2(const double MZ2, const double mHl2, const double mHp2) const;
    gslpp::complex B0p_MZ2_mHl2_0_mA2(const double MZ2, const double mHl2, const double mA2) const;
    gslpp::complex B0p_MZ2_mHl2_mHl2_mHl2(const double MZ2, const double mHl2) const;
    gslpp::complex B0p_MZ2_mHl2_mHh2_mHl2(const double MZ2, const double mHl2, const double mHh2) const;
    gslpp::complex B0p_MZ2_mHl2_mHh2_mHh2(const double MZ2, const double mHl2, const double mHh2) const;
    gslpp::complex B0p_MZ2_mHl2_mHp2_mHp2(const double MZ2, const double mHl2, const double mHp2) const;
    gslpp::complex B0p_MZ2_mHl2_mA2_mA2(const double MZ2, const double mHl2, const double mA2) const;
    gslpp::complex B0p_MZ2_mHh2_0_0(const double MZ2, const double mHh2) const;
    gslpp::complex B0p_MZ2_mHh2_0_mHp2(const double MZ2, const double mHh2, const double mHp2) const;
    gslpp::complex B0p_MZ2_mHh2_0_mA2(const double MZ2, const double mHh2, const double mA2) const;
    gslpp::complex B0p_MZ2_mHh2_mHl2_mHl2(const double MZ2, const double mHh2, const double mHl2) const;
    gslpp::complex B0p_MZ2_mHh2_mHh2_mHl2(const double MZ2, const double mHh2, const double mHl2) const;
    gslpp::complex B0p_MZ2_mHh2_mHh2_mHh2(const double MZ2, const double mHh2) const;
    gslpp::complex B0p_MZ2_mHh2_mHp2_mHp2(const double MZ2, const double mHh2, const double mHp2) const;
    gslpp::complex B0p_MZ2_mHh2_mA2_mA2(const double MZ2, const double mHh2, const double mA2) const;
    gslpp::complex B0p_MZ2_mHp2_0_mHl2(const double MZ2, const double mHp2, const double mHl2) const;
    gslpp::complex B0p_MZ2_mHp2_0_mHh2(const double MZ2, const double mHp2, const double mHh2) const;
    gslpp::complex B0p_MZ2_mHp2_0_mA2(const double MZ2, const double mHp2, const double mA2) const;
    gslpp::complex B0p_MZ2_mHp2_mHp2_mHl2(const double MZ2, const double mHp2, const double mHl2) const;
    gslpp::complex B0p_MZ2_mHp2_mHp2_mHh2(const double MZ2, const double mHp2, const double mHh2) const;
    gslpp::complex B0p_MZ2_mA2_0_mHl2(const double MZ2, const double mA2, const double mHl2) const;
    gslpp::complex B0p_MZ2_mA2_0_mHh2(const double MZ2, const double mA2, const double mHh2) const;
    gslpp::complex B0p_MZ2_mA2_0_mHp2(const double MZ2, const double mA2, const double mHp2) const;
    gslpp::complex B0p_MZ2_mA2_mA2_mHl2(const double MZ2, const double mA2, const double mHl2) const;
    gslpp::complex B0p_MZ2_mA2_mA2_mHh2(const double MZ2, const double mA2, const double mHh2) const;

    gslpp::complex B00_MZ2_0_mA2_mHp2(const double MZ2, const double mA2, const double mHp2) const;
    gslpp::complex B00_MZ2_0_mHh2_mA2(const double MZ2, const double mHh2, const double mA2) const;
    gslpp::complex B00_MZ2_0_mHh2_mHp2(const double MZ2, const double mHh2, const double mHp2) const;
    gslpp::complex B00_MZ2_0_mHl2_mA2(const double MZ2, const double mHl2, const double mA2) const;
    gslpp::complex B00_MZ2_0_mHl2_mHp2(const double MZ2, const double mHl2, const double mHp2) const;
    gslpp::complex B00_MZ2_0_mHp2_mHp2(const double MZ2, const double mHp2) const;
    gslpp::complex B00_MZ2_0_MW2_mHh2(const double MZ2, const double MW2, const double mHh2) const;
    gslpp::complex B00_MZ2_0_MW2_mHl2(const double MZ2, const double MW2, const double mHl2) const;
    gslpp::complex B00_MZ2_0_MZ2_mHh2(const double MZ2, const double mHh2) const;
    gslpp::complex B00_MZ2_0_MZ2_mHl2(const double MZ2, const double mHl2) const;
    gslpp::complex B00_MZ2_MW2_mA2_mHp2(const double MZ2, const double MW2, const double mA2, const double mHp2) const;
    gslpp::complex B00_MZ2_MW2_mHh2_mHp2(const double MZ2, const double MW2, const double mHh2, const double mHp2) const;
    gslpp::complex B00_MZ2_MW2_mHl2_mHp2(const double MZ2, const double MW2, const double mHl2, const double mHp2) const;
    gslpp::complex B00_MZ2_MW2_mHp2_mHp2(const double MZ2, const double MW2, const double mHp2) const;
    gslpp::complex B00_MZ2_MW2_MW2_mHh2(const double MZ2, const double MW2, const double mHh2) const;
    gslpp::complex B00_MZ2_MW2_MW2_mHl2(const double MZ2, const double MW2, const double mHl2) const;
    gslpp::complex B00_MZ2_MZ2_mHh2_mA2(const double MZ2, const double mHh2, const double mA2) const;
    gslpp::complex B00_MZ2_MZ2_mHl2_mA2(const double MZ2, const double mHl2, const double mA2) const;
    gslpp::complex B00_MZ2_MZ2_mHp2_mHp2(const double MZ2, const double mHp2) const;
    gslpp::complex B00_MZ2_MZ2_MZ2_mHh2(const double MZ2, const double mHh2) const;
    gslpp::complex B00_MZ2_MZ2_MZ2_mHl2(const double MZ2, const double mHl2) const;

    double ghHpHm(const double mHp2, const double tanb, const double m12_2, const double bma, const double mHl2, const double vev) const;
    double g_HH_HpHm(const double mHp2, const double mHh2, const double tanb, const double m12_2, const double bma, const double vev) const;

    gslpp::complex I_h_U(const double mHl2, const double Mu, const double Mc, const double Mt) const;
    gslpp::complex I_HH_U(const double mHh2, const double Mc, const double Mt) const;
    gslpp::complex I_A_U(const double mA2, const double Mc, const double Mt) const;
    gslpp::complex I_h_D(const double mHl2, const double Md, const double Ms, const double Mb) const;
    gslpp::complex I_HH_D(const double mHh2, const double Ms, const double Mb) const;
    gslpp::complex I_A_D(const double mA2, const double Ms, const double Mb) const;
    gslpp::complex I_h_L(const double mHl2, const double Me, const double Mmu, const double Mtau) const;
    gslpp::complex I_HH_L(const double mHh2, const double Mmu, const double Mtau) const;
    gslpp::complex I_A_L(const double mA2, const double Mmu, const double Mtau) const;
    gslpp::complex I_H_W(const double mH, const double MW) const;
    gslpp::complex I_H_Hp(const double mHp2, const double mH) const;

    gslpp::complex A_h_U(const double mHl2, const double cW2, const double Mu, const double Mc, const double Mt, const double MZ) const;
    gslpp::complex A_HH_U(const double mHh2, const double cW2, const double Mc, const double Mt, const double MZ) const;
    gslpp::complex A_A_U(const double mA2, const double cW2, const double Mc, const double Mt, const double MZ) const;
    gslpp::complex A_h_D(const double mHl2, const double cW2, const double Md, const double Ms, const double Mb, const double MZ) const;
    gslpp::complex A_HH_D(const double mHh2, const double cW2, const double Ms, const double Mb, const double MZ) const;
    gslpp::complex A_A_D(const double mA2, const double cW2, const double Ms, const double Mb, const double MZ) const;
    gslpp::complex A_h_L(const double mHl2, const double cW2, const double Me, const double Mmu, const double Mtau, const double MZ) const;
    gslpp::complex A_HH_L(const double mHh2, const double cW2, const double Mmu, const double Mtau, const double MZ) const;
    gslpp::complex A_A_L(const double mA2, const double cW2, const double Mmu, const double Mtau, const double MZ) const;
    gslpp::complex A_H_W(const double mH, const double cW2, const double MW, const double MZ) const;
    gslpp::complex A_H_Hp(const double mHp2, const double mH, const double cW2, const double MZ) const;

    /**
     * @brief Kaellen function
     * @return @f$\kappa(a,b,c)=\frac{1}{2a}\sqrt{a^2+b^a+c^2-2ab-2ac-2bc}@f$
     */
    double KaellenFunction(const double a2, const double b2, const double c2) const;

    double cW2THDM(const double c02) const;
    double MWTHDM(const double MW) const;

    void computeSignalStrengthQuantities();

    void computeHHquantities();

    void computeAquantities();

    void computeHpquantities();

    ////////////////////////////////////////////////////////////////////////////

    void runTHDMparameters();

    void computeWFRcombinations();

    ////////////////////////////////////////////////////////////////////////////
    double updateCache();

    /**
     * @brief SM branching ratio of @f$h\to b \bar b@f$.
     * @return @f$BR{\text SM}(h\to b \bar b)@f$
     */
    double BrSM_htobb;

    /**
     * @brief SM branching ratio of @f$h\to \gamma \gamma@f$.
     * @return @f$BR{\text SM}(h\to \gamma \gamma)@f$
     */
    double BrSM_htogaga;

    /**
     * @brief SM branching ratio of @f$h\to \tau \tau@f$.
     * @return @f$BR{\text SM}(h\to \tau \tau)@f$
     */
    double BrSM_htotautau;

    /**
     * @brief Squared relative coupling of @f$h@f$ to two down type quarks.
     * @return @f$r^{(h)}_{Q_dQ_d}@f$
     * @details Depends on the type of @f$Z_2@f$ symmetry.
     */
    double rh_QdQd;

    /**
     * @brief Squared relative coupling of @f$h@f$ to two massive vector bosons.
     * @return @f$r^{(h)}_{WW}=r^{(h)}_{ZZ}@f$
     */
    double rh_VV;

    /**
     * @brief Squared relative coupling of @f$h@f$ to two charged leptons.
     * @return @f$r^{(h)}_{\ell \ell}@f$
     * @details Depends on the type of @f$Z_2@f$ symmetry.
     */
    double rh_ll;

    /**
     * @brief Squared relative coupling of @f$h@f$ to two photons.
     * @return @f$r^{(h)}_{\gamma \gamma}@f$
     * @details Depends on the type of @f$Z_2@f$ symmetry.
     */
    double rh_gaga;

    /**
     * @brief Squared relative coupling of @f$h@f$ to two gluons.
     * @return @f$r^{(h)}_{gg}@f$
     * @details Depends on the type of @f$Z_2@f$ symmetry.
     */
    double rh_gg;

    /**
     * @brief Ratio of THDM and SM cross sections for ggF and tth production of h at 8 TeV.
     * @return @f$\sigma^{\text THDM}_{\text ggF+tth}/\sigma^{\text SM}_{\text ggF+tth}@f$
     */
    double ggF_tth8;

    /**
     * @brief Ratio of THDM and SM cross sections for ggF and tth production of h at 13 TeV.
     * @return @f$\sigma^{\text THDM}_{\text ggF+tth}/\sigma^{\text SM}_{\text ggF+tth}@f$
     */
    double ggF_tth13;

    /**
     * @brief Ratio of THDM and SM cross sections for the production of h at 13 TeV.
     * @return @f$\sigma^{\text THDM}_{\text ggF+VBF+Vh+tth}/\sigma^{\text SM}_{\text ggF+VBF+Vh+tth}@f$
     */
    double pph13;

    /**
     * @brief Ratio of THDM and SM cross sections for VBF and Vh production of h.
     * @return @f$\sigma^{\text THDM}_{\text VBF+Vh}/\sigma^{\text SM}_{\text VBF+Vh}@f$
     */
    double VBF_Vh;

    /**
     * @brief Sum of the modified branching ratios.
     * @return @f$\sum _i r^{(h)}_{i} BR^{\text SM}(h\to i)@f$
     */
    double sumModBRs;

    /**
     * @brief Total h decay rate in the THDM.
     * @return @f$\Gamma_h@f$
     */
    double Gamma_h;

    /**
     * @brief Squared relative coupling of @f$h@f$ to a @f$Z@f$ boson and a photon.
     * @return @f$r^{(h)}_{Z\gamma}@f$
     * @details Depends on the type of @f$Z_2@f$ symmetry.
     */
    double rh_Zga;

    /**
     * @brief Squared relative coupling of @f$h@f$ to two up type quarks.
     * @return @f$r^{(h)}_{Q_uQ_u}@f$
     */
    double rh_QuQu;

    /**
     * @brief @f$h@f$ branching ratio to two @f$b@f$ quarks in the %THDM.
     * @return @f$BR^{\text{THDM}}(h\to b \bar b)@f$
     */
    double THDM_BR_h_bb;

    /**
     * @brief @f$h@f$ branching ratio to two photons in the %THDM.
     * @return @f$BR^{\text{THDM}}(h\to \gamma \gamma)@f$
     */
    double THDM_BR_h_gaga;

    /**
     * @brief @f$h@f$ branching ratio to two @f$\tau@f$ leptons in the %THDM.
     * @return @f$BR^{\text{THDM}}(h\to \tau\tau )@f$
     */
    double THDM_BR_h_tautau;

    /**
     * @brief @f$h@f$ branching ratio to two @f$W@f$ bosons in the %THDM.
     * @return @f$BR^{\text{THDM}}(h\to WW)@f$
     */
    double THDM_BR_h_WW;

    /**
     * @brief @f$h@f$ branching ratio to two @f$Z@f$ bosons in the %THDM.
     * @return @f$BR^{\text{THDM}}(h\to ZZ)@f$
     */
    double THDM_BR_h_ZZ;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to H\to \tau\tau@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)@f$
     */
    double ggF_H_tautau_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$b\bar b\to H\to \tau\tau@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)@f$
     */
    double bbF_H_tautau_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to H\to \gamma\gamma@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to \gamma\gamma)@f$
     */
    double pp_H_gaga_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to H\to \gamma\gamma@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \gamma\gamma)@f$
     */
    double ggF_H_gaga_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to H\to Z\gamma \to \ell \ell \gamma@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to Z\gamma \to \ell \ell \gamma)@f$
     */
    double pp_H_Zga_llga_TH8;

    /**
     * @brief Signal strength for the process @f$pp\to H\to VV@f$ with $VV=WW,ZZ$ at the LHC with 8 TeV.
     * @return @f$\mu_H^{\text{THDM}}(H\to VV)=[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to VV)] / [\sigma^{\text{SM}}_{pp\to H}\cdot BR^{\text{SM}}(H\to VV)]@f$
     */
    double mu_pp_H_VV_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to H\to ZZ@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to ZZ)@f$
     */
    double ggF_H_ZZ_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$VV\to H\to ZZ@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to ZZ)@f$
     */
    double VBF_H_ZZ_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to H\to WW@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to WW)@f$
     */
    double ggF_H_WW_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$VV\to H\to WW@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to WW)@f$
     */
    double VBF_H_WW_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to H\to hh@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to hh)@f$
     */
    double ggF_H_hh_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to H\to hh@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh)@f$
     */
    double pp_H_hh_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to H\to hh\to b\bar b \tau\tau@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b \tau\tau)@f$
     */
    double ggF_H_hh_bbtautau_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to H\to hh\to b\bar b b\bar b@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b b\bar b)@f$
     */
    double pp_H_hh_bbbb_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to H\to hh\to \gamma\gamma b\bar b@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to \gamma\gamma b\bar b)@f$
     */
    double pp_H_hh_gagabb_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to H\to t\bar t@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to t\bar t)@f$
     */
    double ggF_H_tt_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$b\bar b\to H\to b\bar b@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to b\bar b)@f$
     */
    double bbF_H_bb_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to H\to \tau\tau@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)@f$
     */
    double ggF_H_tautau_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$b\bar b\to H\to \tau\tau@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)@f$
     */
    double bbF_H_tautau_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to H\to \gamma\gamma@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to \gamma\gamma)@f$
     */
    double pp_H_gaga_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to H\to \gamma\gamma@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \gamma\gamma)@f$
     */
    double ggF_H_gaga_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to H\to Z\gamma@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to Z\gamma)@f$
     */
    double pp_H_Zga_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to H\to ZZ@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to ZZ)@f$
     */
    double pp_H_ZZ_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to H\to ZZ@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to ZZ)@f$
     */
    double ggF_H_ZZ_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$VV\to H\to ZZ@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to ZZ)@f$
     */
    double VBF_H_ZZ_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to H\to ZZ\to 4\ell@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to ZZ\to 4\ell)@f$
     */
    double ggF_H_ZZ_llll_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$VV\to H\to ZZ\to 4\ell@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to ZZ\to 4\ell)@f$
     */
    double VBF_H_ZZ_llll_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to H\to ZZ\to 4\ell@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to ZZ\to 4\ell)@f$
     */
    double pp_H_ZZ_llll_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$VV+VH\to H\to ZZ\to 4\ell@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{VV+VH\to H}\cdot BR^{\text{THDM}}(H\to ZZ\to 4\ell)@f$
     */
    double VBF_VH_H_ZZ_llll_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to H\to WW@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to WW)@f$
     */
    double ggF_H_WW_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$VV\to H\to WW@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to WW)@f$
     */
    double VBF_H_WW_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$(gg+VV)\to H\to WW\to \ell \nu \ell \nu@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{(gg+VV)\to H}\cdot BR^{\text{THDM}}(H\to WW\to \ell \nu \ell \nu)@f$
     */
    double ggF_VBF_H_WW_lnulnu_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to H\to hh@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to hh)@f$
     */
    double ggF_H_hh_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to H\to hh@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh)@f$
     */
    double pp_H_hh_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to H\to hh\to b\bar b b\bar b@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b b\bar b)@f$
     */
    double pp_H_hh_bbbb_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to H\to hh\to \gamma\gamma b\bar b@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to \gamma\gamma b\bar b)@f$
     */
    double pp_H_hh_gagabb_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to H\to hh\to b\bar b \tau\tau@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b \tau\tau)@f$
     */
    double pp_H_hh_bbtautau_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to H\to hh\to b\bar b WW\to b\bar b \ell \nu \ell \nu@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b WW\to b\bar b \ell \nu \ell \nu)@f$
     */
    double pp_H_hh_bblnulnu_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to H\to hh\to b\bar b VV(\ell\ell \nu\nu)@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}[H\to hh\to b\bar b VV(\ell\ell \nu\nu)]@f$
     */
    double pp_H_hh_bbVV_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$t\bar t\to H\to t\bar t@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{t\bar t\to H}\cdot BR^{\text{THDM}}(H\to t\bar t)@f$
     */
    double ttF_H_tt_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$b\bar b\to H\to t\bar t@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to t\bar t)@f$
     */
    double bbF_H_tt_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to H\to b\bar b@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to b\bar b)@f$
     */
    double pp_H_bb_TH13;

    /**
     * @brief Total decay width of the heavy CP-even Higgs @f$H@f$.
     * @return @f$\Gamma_H@f$
     */
    double GammaHtot;

//    /**
//     * @brief Squared relative coupling of @f$H@f$ to two photons.
//     * @return @f$r^{(H)}_{\gamma \gamma}@f$
//     * @details Depends on the type of @f$Z_2@f$ symmetry.
//     */
//    double rHH_gaga;

    /**
     * @brief Squared relative coupling of @f$H@f$ to two gluons.
     * @return @f$r^{(H)}_{gg}@f$
     * @details Depends on the type of @f$Z_2@f$ symmetry.
     */
    double rHH_gg;

    /**
     * @brief @f$H@f$ branching ratio to two light Higgs bosons.
     * @return @f$BR(H\to hh)@f$
     */
    double Br_Htohh;

    /**
     * @brief @f$H@f$ branching ratio to two CP-odd Higgs bosons.
     * @return @f$BR(H\to AA)@f$
     */
    double Br_HtoAA;

    /**
     * @brief @f$H@f$ branching ratio to two charged Higgs bosons.
     * @return @f$BR(H\to H^\pm H^\mp)@f$
     */
    double Br_HtoHpHm;

    /**
     * @brief @f$H@f$ branching ratio to two light Higgs bosons.
     * @return @f$BR(H\to hh)@f$
     */
    double Br_HtoAZ;

    /**
     * @brief @f$H@f$ branching ratio to two light Higgs bosons.
     * @return @f$BR(H\to hh)@f$
     */
    double Br_HtoHpW;

    double THoEX_ggF_H_tautau_ATLAS8;
    double R_ggF_H_tautau_ATLAS8;
    double THoEX_ggF_H_tautau_CMS8;
    double R_ggF_H_tautau_CMS8;
    double THoEX_bbF_H_tautau_ATLAS8;
    double R_bbF_H_tautau_ATLAS8;
    double THoEX_bbF_H_tautau_CMS8;
    double R_bbF_H_tautau_CMS8;
    double THoEX_pp_H_gaga_ATLAS8;
    double R_pp_H_gaga_ATLAS8;
    double THoEX_ggF_H_gaga_CMS8;
    double R_ggF_H_gaga_CMS8;
//    double LIMIT_ggF_H_gaga_CMS8;
//    double LIMEST_ggF_H_gaga_CMS8;
//    double DEVIATION_ggF_H_gaga_CMS8;
//    double BANDSIZE_ggF_H_gaga_CMS8;
    double THoEX_pp_H_Zga_llga_CMS8;
    double R_pp_H_Zga_llga_CMS8;
    double THoEX_pp_H_Zga_llga_ATLAS8;
    double R_pp_H_Zga_llga_ATLAS8;
    double THoEX_mu_pp_H_VV_CMS8;
    double R_mu_pp_H_VV_CMS8;
    double THoEX_ggF_H_WW_ATLAS8;
    double R_ggF_H_WW_ATLAS8;
    double THoEX_VBF_H_WW_ATLAS8;
    double R_VBF_H_WW_ATLAS8;
    double THoEX_ggF_H_ZZ_ATLAS8;
    double R_ggF_H_ZZ_ATLAS8;
    double THoEX_VBF_H_ZZ_ATLAS8;
    double R_VBF_H_ZZ_ATLAS8;
    double THoEX_ggF_H_hh_ATLAS8;
    double R_ggF_H_hh_ATLAS8;
    double THoEX_pp_H_hh_CMS8;
    double R_pp_H_hh_CMS8;
    double THoEX_ggF_H_hh_bbtautau_CMS8;
    double R_ggF_H_hh_bbtautau_CMS8;
    double THoEX_pp_H_hh_bbbb_CMS8;
    double R_pp_H_hh_bbbb_CMS8;
    double THoEX_pp_H_hh_gagabb_CMS8;
    double R_pp_H_hh_gagabb_CMS8;
    double THoEX_ggF_H_tt_ATLAS8;
    double R_ggF_H_tt_ATLAS8;
    double THoEX_bbF_H_bb_CMS8;
    double R_bbF_H_bb_CMS8;

    double THoEX_ttF_H_tt_ATLAS13;
    double R_ttF_H_tt_ATLAS13;
    double THoEX_bbF_H_tt_ATLAS13;
    double R_bbF_H_tt_ATLAS13;
    double THoEX_ggF_H_tautau_ATLAS13;
    double R_ggF_H_tautau_ATLAS13;
    double THoEX_bbF_H_tautau_ATLAS13;
    double R_bbF_H_tautau_ATLAS13;
    double THoEX_ggF_H_tautau_CMS13;
    double R_ggF_H_tautau_CMS13;
    double THoEX_bbF_H_tautau_CMS13;
    double R_bbF_H_tautau_CMS13;
    double THoEX_pp_H_gaga_ATLAS13;
    double R_pp_H_gaga_ATLAS13;
    double THoEX_ggF_H_gaga_CMS13;
    double R_ggF_H_gaga_CMS13;
    double THoEX_pp_H_Zga_ATLAS13;
    double R_pp_H_Zga_ATLAS13;
    double THoEX_pp_H_Zga_llga_CMS13;
    double R_pp_H_Zga_llga_CMS13;
    double THoEX_pp_H_Zga_qqga_CMS13;
    double R_pp_H_Zga_qqga_CMS13;
    double THoEX_ggF_H_ZZ_llnunu_ATLAS13;
    double R_ggF_H_ZZ_llnunu_ATLAS13;
    double THoEX_ggF_H_ZZ_llll_ATLAS13;
    double R_ggF_H_ZZ_llll_ATLAS13;
    double THoEX_VBF_H_ZZ_llll_ATLAS13;
    double R_VBF_H_ZZ_llll_ATLAS13;
    double THoEX_pp_H_ZZ_llll_CMS13;
    double R_pp_H_ZZ_llll_CMS13;
    double THoEX_VBF_VH_H_ZZ_llll_CMS13;
    double R_VBF_VH_H_ZZ_llll_CMS13;
    double THoEX_ggF_H_ZZ_llqq_ATLAS13;
    double R_ggF_H_ZZ_llqq_ATLAS13;
    double THoEX_VBF_H_ZZ_llqq_ATLAS13;
    double R_VBF_H_ZZ_llqq_ATLAS13;
    double THoEX_ggF_H_ZZ_nunuqq_ATLAS13;
    double R_ggF_H_ZZ_nunuqq_ATLAS13;
    double THoEX_pp_H_ZZ_llqq_CMS13;
    double R_pp_H_ZZ_llqq_CMS13;
    double THoEX_ggF_H_WW_lnuqq_ATLAS13;
    double R_ggF_H_WW_lnuqq_ATLAS13;
    double THoEX_ggF_H_WW_enumunu_ATLAS13;
    double R_ggF_H_WW_enumunu_ATLAS13;
    double THoEX_VBF_H_WW_enumunu_ATLAS13;
    double R_VBF_H_WW_enumunu_ATLAS13;
    double THoEX_ggF_VBF_H_WW_lnulnu_CMS13;
    double R_ggF_VBF_H_WW_lnulnu_CMS13;
    double THoEX_pp_H_hh_bbgaga_ATLAS13;
    double R_pp_H_hh_bbgaga_ATLAS13;
    double THoEX_pp_H_hh_bbgaga_CMS13;
    double R_pp_H_hh_bbgaga_CMS13;
    double THoEX_pp_H_hh_bbbb_ATLAS13;
    double R_pp_H_hh_bbbb_ATLAS13;
    double THoEX_pp_H_hh_bbbb_CMS13;
    double R_pp_H_hh_bbbb_CMS13;
    double THoEX_ggF_H_hh_gagaWW_ATLAS13;
    double R_ggF_H_hh_gagaWW_ATLAS13;
    double THoEX_pp_H_hh_bbtautau_CMS13;
    double R_pp_H_hh_bbtautau_CMS13;
    double THoEX_pp_H_hh_bbtautau1_CMS13;
    double R_pp_H_hh_bbtautau1_CMS13;
    double THoEX_pp_H_hh_bblnulnu_CMS13;
    double R_pp_H_hh_bblnulnu_CMS13;
    double THoEX_pp_H_hh_bbVV_CMS13;
    double R_pp_H_hh_bbVV_CMS13;
    double THoEX_pp_H_bb_CMS13;
    double R_pp_H_bb_CMS13;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to A\to \tau\tau@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to \tau\tau)@f$
     */
    double ggF_A_tautau_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$b\bar b\to A\to \tau\tau@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{b\bar b\to A}\cdot BR^{\text{THDM}}(A\to \tau\tau)@f$
     */
    double bbF_A_tautau_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to A\to \gamma\gamma@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to A}\cdot BR^{\text{THDM}}(A\to \gamma\gamma)@f$
     */
    double pp_A_gaga_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to A\to Z\gamma \to \ell \ell \gamma@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to A}\cdot BR^{\text{THDM}}(A\to Z\gamma \to \ell \ell \gamma)@f$
     */
    double pp_A_Zga_llga_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to A\to \gamma\gamma@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to \gamma\gamma)@f$
     */
    double ggF_A_gaga_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to A\to hZ@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ)@f$
     */
    double ggF_A_hZ_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to A\to hZ \to b\bar b \ell \ell@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ \to b\bar b \ell \ell)@f$
     */
    double ggF_A_hZ_bbll_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to A\to hZ \to b\bar b Z@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ \to b\bar b Z)@f$
     */
    double ggF_A_hZ_bbZ_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to A\to hZ \to \tau \tau \ell \ell@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ \to \tau \tau \ell \ell)@f$
     */
    double ggF_A_hZ_tautaull_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to A\to hZ \to \tau \tau Z@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ \to \tau \tau Z)@f$
     */
    double ggF_A_hZ_tautauZ_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to A\to t\bar t@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to t\bar t)@f$
     */
    double ggF_A_tt_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$b\bar b\to A\to b\bar b@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{b\bar b\to A}\cdot BR^{\text{THDM}}(A\to b\bar b)@f$
     */
    double bbF_A_bb_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to A\to \tau \tauZ@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to \tau \tau)@f$
     */
    double ggF_A_tautau_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$b\bar b\to A\to \tau \tauZ@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{b\bar b\to A}\cdot BR^{\text{THDM}}(A\to \tau \tau)@f$
     */
    double bbF_A_tautau_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to A\to \gamma \gamma@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to A}\cdot BR^{\text{THDM}}(A\to \gamma \gamma)@f$
     */
    double pp_A_gaga_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to A\to \gamma \gamma@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to \gamma \gamma)@f$
     */
    double ggF_A_gaga_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to A\to Z\gamma@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to A}\cdot BR^{\text{THDM}}(A\to Z\gamma)@f$
     */
    double pp_A_Zga_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to A\to hZ\to b\bar b Z@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ\to b\bar b Z)@f$
     */
    double ggF_A_hZ_bbZ_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$b\bar b\to A\to hZ\to b\bar b Z@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{b\bar b\to A}\cdot BR^{\text{THDM}}(A\to hZ\to b\bar b Z)@f$
     */
    double bbF_A_hZ_bbZ_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$t\bar t\to A\to t\bar t@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{t\bar t\to A}\cdot BR^{\text{THDM}}(A\to t\bar t)@f$
     */
    double ttF_A_tt_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$b\bar b\to A\to t\bar t@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{b\bar b\to A}\cdot BR^{\text{THDM}}(A\to t\bar t)@f$
     */
    double bbF_A_tt_TH13;
    
    /**
     * @brief Cross section times branching ratio for the process @f$pp\to A\to b\bar b@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to A}\cdot BR^{\text{THDM}}(A\to b\bar b)@f$
     */
    double pp_A_bb_TH13;

    /**
     * @brief Total decay width of the CP-odd Higgs @f$A@f$.
     * @return @f$\Gamma_A@f$
     */
    double GammaAtot;

//    /**
//     * @brief Squared relative coupling of @f$A@f$ to two photons.
//     * @return @f$r^{(A)}_{\gamma \gamma}@f$
//     * @details Depends on the type of @f$Z_2@f$ symmetry.
//     */
//    double rA_gaga;

    /**
     * @brief Squared relative coupling of @f$A@f$ to two gluons.
     * @return @f$r^{(A)}_{gg}@f$
     * @details Depends on the type of @f$Z_2@f$ symmetry.
     */
    double rA_gg;

    /**
     * @brief @f$A@f$ branching ratio to an @f$H@f$ and a @f$Z@f$ boson.
     * @return @f$BR(A\to HZ)@f$
     */
    double Br_AtoHZ;

    /**
     * @brief @f$A@f$ branching ratio to an @f$h@f$ and a @f$Z@f$ boson.
     * @return @f$BR(A\to hZ)@f$
     */
    double Br_AtohZ;

    /**
     * @brief @f$A@f$ branching ratio to a charged Higgs and a @f$W^\mp@f$ boson.
     * @return @f$BR(A\to H^\pm W^\mp)@f$
     */
    double Br_AtoHpW;

    double THoEX_ggF_A_tautau_ATLAS8;
    double R_ggF_A_tautau_ATLAS8;
    double THoEX_ggF_A_tautau_CMS8;
    double R_ggF_A_tautau_CMS8;
    double THoEX_bbF_A_tautau_ATLAS8;
    double R_bbF_A_tautau_ATLAS8;
    double THoEX_bbF_A_tautau_CMS8;
    double R_bbF_A_tautau_CMS8;
    double THoEX_pp_A_gaga_ATLAS8;
    double R_pp_A_gaga_ATLAS8;
    double THoEX_ggF_A_gaga_CMS8;
    double R_ggF_A_gaga_CMS8;
    double THoEX_pp_A_Zga_llga_CMS8;
    double R_pp_A_Zga_llga_CMS8;
    double THoEX_pp_A_Zga_llga_ATLAS8;
    double R_pp_A_Zga_llga_ATLAS8;
    double THoEX_ggF_A_hZ_bbll_CMS8;
    double R_ggF_A_hZ_bbll_CMS8;
    double THoEX_ggF_A_hZ_bbZ_ATLAS8;
    double R_ggF_A_hZ_bbZ_ATLAS8;
    double THoEX_ggF_A_hZ_tautaull_CMS8;
    double R_ggF_A_hZ_tautaull_CMS8;
    double THoEX_ggF_A_hZ_tautauZ_ATLAS8;
    double R_ggF_A_hZ_tautauZ_ATLAS8;
    double THoEX_ggF_A_tt_ATLAS8;
    double R_ggF_A_tt_ATLAS8;
    double THoEX_bbF_A_bb_CMS8;
    double R_bbF_A_bb_CMS8;

    double THoEX_ttF_A_tt_ATLAS13;
    double R_ttF_A_tt_ATLAS13;
    double THoEX_bbF_A_tt_ATLAS13;
    double R_bbF_A_tt_ATLAS13;
    double THoEX_ggF_A_tautau_ATLAS13;
    double R_ggF_A_tautau_ATLAS13;
    double THoEX_bbF_A_tautau_ATLAS13;
    double R_bbF_A_tautau_ATLAS13;
    double THoEX_ggF_A_tautau_CMS13;
    double R_ggF_A_tautau_CMS13;
    double THoEX_bbF_A_tautau_CMS13;
    double R_bbF_A_tautau_CMS13;
    double THoEX_pp_A_gaga_ATLAS13;
    double R_pp_A_gaga_ATLAS13;
    double THoEX_ggF_A_gaga_CMS13;
    double R_ggF_A_gaga_CMS13;
    double THoEX_pp_A_Zga_llga_ATLAS13;
    double R_pp_A_Zga_llga_ATLAS13;
    double THoEX_pp_A_Zga_llga_CMS13;
    double R_pp_A_Zga_llga_CMS13;
    double THoEX_pp_A_Zga_qqga_CMS13;
    double R_pp_A_Zga_qqga_CMS13;
    double THoEX_ggF_A_hZ_bbZ_ATLAS13;
    double R_ggF_A_hZ_bbZ_ATLAS13;
    double THoEX_bbF_A_hZ_bbZ_ATLAS13;
    double R_bbF_A_hZ_bbZ_ATLAS13;
    double THoEX_pp_A_bb_CMS13;
    double R_pp_A_bb_CMS13;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to H^\pm\to \tau^\pm \nu@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to H^\pm}\cdot BR^{\text{THDM}}(H^\pm\to \tau^\pm \nu)@f$
     */
    double pp_Hpm_taunu_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to H^+\to t\bar b@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to H^+}\cdot BR^{\text{THDM}}(H^+\to t\bar b)@f$
     */
    double pp_Hp_tb_TH8;

    /**
     * @brief Cross section for the process @f$pp\to H^+@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to H^+}@f$
     */
    double pp_Hp_TH8;

    /**
     * @brief Cross section for the process @f$pp\to H^+\to \tau^+ \nu@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to H^+\to \tau^+ \nu}@f$
     */
    double pp_Hp_taunu_TH8;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to H^\pm\to \tau^\pm \nu@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to H^\pm}\cdot BR^{\text{THDM}}(H^\pm\to \tau^\pm \nu)@f$
     */
    double pp_Hpm_taunu_TH13;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to H^+\to t\bar b@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to H^+}\cdot BR^{\text{THDM}}(H^+\to t\bar b)@f$
     */
    double pp_Hp_tb_TH13;

    double THoEX_pp_Hpm_taunu_ATLAS8;
    double R_pp_Hpm_taunu_ATLAS8;
    double THoEX_pp_Hp_taunu_CMS8;
    double R_pp_Hp_taunu_CMS8;
    double THoEX_pp_Hp_tb_CMS8;
    double R_pp_Hp_tb_CMS8;
    double THoEX_pp_Hp_tb_ATLAS8;
    double R_pp_Hp_tb_ATLAS8;
    double THoEX_pp_Hpm_taunu_ATLAS13;
    double R_pp_Hpm_taunu_ATLAS13;
    double THoEX_pp_Hpm_taunu_CMS13;
    double R_pp_Hpm_taunu_CMS13;
    double THoEX_pp_Hp_tb_ATLAS13_1;
    double R_pp_Hp_tb_ATLAS13_1;
    double THoEX_pp_Hp_tb_ATLAS13_2;
    double R_pp_Hp_tb_ATLAS13_2;
    double THoEX_pp_Hp_tb_ATLAS13;
    double R_pp_Hp_tb_ATLAS13;

    /**
     * @brief Total decay width of the charged Higgs @f$H^+@f$.
     * @return @f$\Gamma_{H^+}@f$
     */
    double GammaHptot;

    double Q_cutoff;
    double g1_at_Q;
    double g2_at_Q;
    double g3_at_Q;
    double Ytop_at_Q;
    double Ybottom1_at_Q;
    double Ybottom2_at_Q;
    double Ytau1_at_Q;
    double Ytau2_at_Q;
    double m11_2_at_Q;
    double m22_2_at_Q;
    double m12_2_at_Q;
    double lambda1_at_Q;
    double lambda2_at_Q;
    double lambda3_at_Q;
    double lambda4_at_Q;
    double lambda5_at_Q;
    double WFRcomb1;
    double WFRcomb2;
    double WFRcomb3;
    double WFRcomb4;
    double Rpeps;

private:

    const THDM * myTHDM;
    const PVfunctions PV;
    m11_2 * mym11_2;
    m22_2 * mym22_2;
    lambda1 * mylambda1;
    lambda2 * mylambda2;
    lambda3 * mylambda3;
    lambda4 * mylambda4;
    lambda5 * mylambda5;
    Runner * myRunner;

    double Q_THDM;
    double bma;
    double logtb;
    double tanb;
    double m12_2;
    double BDtaunu_SM;
    double BDtaunu_A;
    double BDtaunu_B;
    double BDstartaunu_SM;
    double BDstartaunu_A;
    double BDstartaunu_B;
    double mHl;
    double mHl2;
    double mHh2;
    double mA2;
    double mHp2;
    double MW;
    double cW2;
    double vev;
    double Ale;
    double Als;
    double Mt;
    double Mb;
    double Mtau;
    double Mc;
    double Ms;
    double Mmu;
    double Mu;
    double Md;
    double Me;
    double MZ;
    std::string modelflag;
    bool WFRflag;

    ////////////////////////////////////////////////////////////////////////////
    //Caches

    /*One-loop functions*/

    mutable gslpp::complex B0_MZ2_0_MW2_mHh2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_0_MW2_mHl2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_0_MZ2_mHh2_cache[3][CacheSize];
    mutable gslpp::complex B0_MZ2_0_MZ2_mHl2_cache[3][CacheSize];
    mutable gslpp::complex B0_MZ2_MW2_MW2_mHh2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_MW2_MW2_mHl2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_MZ2_MZ2_mHh2_cache[3][CacheSize];
    mutable gslpp::complex B0_MZ2_MZ2_MZ2_mHl2_cache[3][CacheSize];

    mutable gslpp::complex B0_MZ2_0_0_mHl2_cache[3][CacheSize];
    mutable gslpp::complex B0_MZ2_0_0_mHh2_cache[3][CacheSize];
    mutable gslpp::complex B0_MZ2_0_mHp2_mHl2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_0_mHp2_mHh2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_0_mA2_mHl2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_0_mA2_mHh2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_mHl2_0_0_cache[3][CacheSize];
    mutable gslpp::complex B0_MZ2_mHl2_0_mHp2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_mHl2_0_mA2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_mHl2_mHl2_mHl2_cache[3][CacheSize];
    mutable gslpp::complex B0_MZ2_mHl2_mHh2_mHl2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_mHl2_mHh2_mHh2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_mHl2_mHp2_mHp2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_mHl2_mA2_mA2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_mHh2_0_0_cache[3][CacheSize];
    mutable gslpp::complex B0_MZ2_mHh2_0_mHp2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_mHh2_0_mA2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_mHh2_mHl2_mHl2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_mHh2_mHh2_mHl2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_mHh2_mHh2_mHh2_cache[3][CacheSize];
    mutable gslpp::complex B0_MZ2_mHh2_mHp2_mHp2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_mHh2_mA2_mA2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_mHp2_0_mHl2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_mHp2_0_mHh2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_mHp2_mHp2_mHl2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_mHp2_mHp2_mHh2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_mA2_0_mHl2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_mA2_0_mHh2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_mA2_mA2_mHl2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_mA2_mA2_mHh2_cache[4][CacheSize];

    mutable gslpp::complex B0p_MZ2_0_0_mHl2_cache[3][CacheSize];
    mutable gslpp::complex B0p_MZ2_0_0_mHh2_cache[3][CacheSize];
    mutable gslpp::complex B0p_MZ2_0_mHp2_mHl2_cache[4][CacheSize];
    mutable gslpp::complex B0p_MZ2_0_mHp2_mHh2_cache[4][CacheSize];
    mutable gslpp::complex B0p_MZ2_0_mHp2_mA2_cache[4][CacheSize];
    mutable gslpp::complex B0p_MZ2_0_mA2_mHl2_cache[4][CacheSize];
    mutable gslpp::complex B0p_MZ2_0_mA2_mHh2_cache[4][CacheSize];
    mutable gslpp::complex B0p_MZ2_mHl2_0_0_cache[3][CacheSize];
    mutable gslpp::complex B0p_MZ2_mHl2_0_mHp2_cache[4][CacheSize];
    mutable gslpp::complex B0p_MZ2_mHl2_0_mA2_cache[4][CacheSize];
    mutable gslpp::complex B0p_MZ2_mHl2_mHl2_mHl2_cache[3][CacheSize];
    mutable gslpp::complex B0p_MZ2_mHl2_mHh2_mHl2_cache[4][CacheSize];
    mutable gslpp::complex B0p_MZ2_mHl2_mHh2_mHh2_cache[4][CacheSize];
    mutable gslpp::complex B0p_MZ2_mHl2_mHp2_mHp2_cache[4][CacheSize];
    mutable gslpp::complex B0p_MZ2_mHl2_mA2_mA2_cache[4][CacheSize];
    mutable gslpp::complex B0p_MZ2_mHh2_0_0_cache[3][CacheSize];
    mutable gslpp::complex B0p_MZ2_mHh2_0_mHp2_cache[4][CacheSize];
    mutable gslpp::complex B0p_MZ2_mHh2_0_mA2_cache[4][CacheSize];
    mutable gslpp::complex B0p_MZ2_mHh2_mHl2_mHl2_cache[4][CacheSize];
    mutable gslpp::complex B0p_MZ2_mHh2_mHh2_mHl2_cache[4][CacheSize];
    mutable gslpp::complex B0p_MZ2_mHh2_mHh2_mHh2_cache[3][CacheSize];
    mutable gslpp::complex B0p_MZ2_mHh2_mHp2_mHp2_cache[4][CacheSize];
    mutable gslpp::complex B0p_MZ2_mHh2_mA2_mA2_cache[4][CacheSize];
    mutable gslpp::complex B0p_MZ2_mHp2_0_mHl2_cache[4][CacheSize];
    mutable gslpp::complex B0p_MZ2_mHp2_0_mHh2_cache[4][CacheSize];
    mutable gslpp::complex B0p_MZ2_mHp2_0_mA2_cache[4][CacheSize];
    mutable gslpp::complex B0p_MZ2_mHp2_mHp2_mHl2_cache[4][CacheSize];
    mutable gslpp::complex B0p_MZ2_mHp2_mHp2_mHh2_cache[4][CacheSize];
    mutable gslpp::complex B0p_MZ2_mA2_0_mHl2_cache[4][CacheSize];
    mutable gslpp::complex B0p_MZ2_mA2_0_mHh2_cache[4][CacheSize];
    mutable gslpp::complex B0p_MZ2_mA2_0_mHp2_cache[4][CacheSize];
    mutable gslpp::complex B0p_MZ2_mA2_mA2_mHl2_cache[4][CacheSize];
    mutable gslpp::complex B0p_MZ2_mA2_mA2_mHh2_cache[4][CacheSize];

    mutable gslpp::complex B00_MZ2_0_mA2_mHp2_cache[4][CacheSize];
    mutable gslpp::complex B00_MZ2_0_mHh2_mA2_cache[4][CacheSize];
    mutable gslpp::complex B00_MZ2_0_mHh2_mHp2_cache[4][CacheSize];
    mutable gslpp::complex B00_MZ2_0_mHl2_mA2_cache[4][CacheSize];
    mutable gslpp::complex B00_MZ2_0_mHl2_mHp2_cache[4][CacheSize];
    mutable gslpp::complex B00_MZ2_0_mHp2_mHp2_cache[3][CacheSize];
    mutable gslpp::complex B00_MZ2_0_MW2_mHh2_cache[4][CacheSize];
    mutable gslpp::complex B00_MZ2_0_MW2_mHl2_cache[4][CacheSize];
    mutable gslpp::complex B00_MZ2_0_MZ2_mHh2_cache[3][CacheSize];
    mutable gslpp::complex B00_MZ2_0_MZ2_mHl2_cache[3][CacheSize];
    mutable gslpp::complex B00_MZ2_MW2_mA2_mHp2_cache[5][CacheSize];
    mutable gslpp::complex B00_MZ2_MW2_mHh2_mHp2_cache[5][CacheSize];
    mutable gslpp::complex B00_MZ2_MW2_mHl2_mHp2_cache[5][CacheSize];
    mutable gslpp::complex B00_MZ2_MW2_mHp2_mHp2_cache[4][CacheSize];
    mutable gslpp::complex B00_MZ2_MW2_MW2_mHh2_cache[4][CacheSize];
    mutable gslpp::complex B00_MZ2_MW2_MW2_mHl2_cache[4][CacheSize];
    mutable gslpp::complex B00_MZ2_MZ2_mHh2_mA2_cache[4][CacheSize];
    mutable gslpp::complex B00_MZ2_MZ2_mHl2_mA2_cache[4][CacheSize];
    mutable gslpp::complex B00_MZ2_MZ2_mHp2_mHp2_cache[3][CacheSize];
    mutable gslpp::complex B00_MZ2_MZ2_MZ2_mHh2_cache[3][CacheSize];
    mutable gslpp::complex B00_MZ2_MZ2_MZ2_mHl2_cache[3][CacheSize];

    mutable double ip_Br_HPtott_cache[2][CacheSize];
    mutable double ip_Br_HPtobb_cache[2][CacheSize];
    mutable double ip_Br_HPtotautau_cache[2][CacheSize];
    mutable double ip_Br_HPtocc_cache[2][CacheSize];
    mutable double ip_Br_HPtomumu_cache[2][CacheSize];
    mutable double ip_Br_HPtoZZ_cache[2][CacheSize];
    mutable double ip_Br_HPtoWW_cache[2][CacheSize];
    mutable double ip_GammaHPtotSM_cache[2][CacheSize];
    mutable double ip_cs_ggtoH_8_cache[2][CacheSize];
    mutable double ip_cs_ggtoH_13_cache[2][CacheSize];
    mutable double ip_cs_VBFtoH_8_cache[2][CacheSize];
    mutable double ip_cs_VBFtoH_13_cache[2][CacheSize];
    mutable double ip_cs_WtoWH_8_cache[2][CacheSize];
    mutable double ip_cs_WtoWH_13_cache[2][CacheSize];
    mutable double ip_cs_ZtoZH_8_cache[2][CacheSize];
    mutable double ip_cs_ZtoZH_13_cache[2][CacheSize];
    mutable double ip_cs_pptottH_8_cache[2][CacheSize];
    mutable double ip_cs_pptottH_13_cache[2][CacheSize];
    mutable double ip_cs_pptobbH_8_cache[2][CacheSize];
    mutable double ip_cs_pptobbH_13_cache[2][CacheSize];
    mutable double ip_cs_ggtoA_8_cache[2][CacheSize];
    mutable double ip_cs_ggtoA_13_cache[2][CacheSize];
    mutable double ip_cs_pptottA_8_cache[2][CacheSize];
    mutable double ip_cs_pptottA_13_cache[2][CacheSize];
    mutable double ip_cs_pptobbA_8_cache[2][CacheSize];
    mutable double ip_cs_pptobbA_13_cache[2][CacheSize];
    mutable double ip_cs_ggtoHp_8_cache[3][CacheSize];
    mutable double ip_cs_ggtoHp_13_cache[3][CacheSize];
    mutable double ip_csr_ggH_t_8_cache[2][CacheSize];
    mutable double ip_csr_ggH_t_13_cache[2][CacheSize];
    mutable double ip_csr_ggH_b_8_cache[2][CacheSize];
    mutable double ip_csr_ggH_b_13_cache[2][CacheSize];
    mutable double ip_csr_ggA_t_8_cache[2][CacheSize];
    mutable double ip_csr_ggA_t_13_cache[2][CacheSize];
    mutable double ip_csr_ggA_b_8_cache[2][CacheSize];
    mutable double ip_csr_ggA_b_13_cache[2][CacheSize];
    mutable double ip_ex_pp_phi_gaga_ATLAS8_cache[2][CacheSize];
    mutable double ip_ex_pp_phi_gaga_ATLAS8_cache_e[2][CacheSize];
    mutable double ip_ex_pp_phi_Zga_llga_ATLAS8_cache[2][CacheSize];
    mutable double ip_ex_pp_phi_Zga_llga_ATLAS8_cache_e[2][CacheSize];
    mutable double ip_ex_gg_phi_tautau_ATLAS8_cache[2][CacheSize];
    mutable double ip_ex_gg_phi_tautau_ATLAS8_cache_e[2][CacheSize];
    mutable double ip_ex_bb_phi_tautau_ATLAS8_cache[2][CacheSize];
    mutable double ip_ex_bb_phi_tautau_ATLAS8_cache_e[2][CacheSize];
    mutable double ip_ex_gg_A_hZ_tautauZ_ATLAS8_cache[2][CacheSize];
    mutable double ip_ex_gg_A_hZ_tautauZ_ATLAS8_cache_e[2][CacheSize];
    mutable double ip_ex_gg_A_hZ_bbZ_ATLAS8_cache[2][CacheSize];
    mutable double ip_ex_gg_A_hZ_bbZ_ATLAS8_cache_e[2][CacheSize];
    mutable double ip_ex_gg_phi_tt_ATLAS8_cache[2][CacheSize];
    mutable double ip_ex_gg_phi_tt_ATLAS8_cache_e[2][CacheSize];
    mutable double ip_ex_gg_H_WW_ATLAS8_cache[2][CacheSize];
    mutable double ip_ex_gg_H_WW_ATLAS8_cache_e[2][CacheSize];
    mutable double ip_ex_VBF_H_WW_ATLAS8_cache[2][CacheSize];
    mutable double ip_ex_VBF_H_WW_ATLAS8_cache_e[2][CacheSize];
    mutable double ip_ex_gg_H_ZZ_ATLAS8_cache[2][CacheSize];
    mutable double ip_ex_gg_H_ZZ_ATLAS8_cache_e[2][CacheSize];
    mutable double ip_ex_VBF_H_ZZ_ATLAS8_cache[2][CacheSize];
    mutable double ip_ex_VBF_H_ZZ_ATLAS8_cache_e[2][CacheSize];
    mutable double ip_ex_gg_H_hh_ATLAS8_cache[2][CacheSize];
    mutable double ip_ex_gg_H_hh_ATLAS8_cache_e[2][CacheSize];
    mutable double ip_ex_mu_pp_H_VV_CMS8_cache[2][CacheSize];
    mutable double ip_ex_mu_pp_H_VV_CMS8_cache_e[2][CacheSize];
    mutable double ip_ex_gg_A_hZ_bbll_CMS8_cache[2][CacheSize];
    mutable double ip_ex_gg_A_hZ_bbll_CMS8_cache_e[2][CacheSize];
    mutable double ip_ex_pp_H_hh_CMS8_cache[2][CacheSize];
    mutable double ip_ex_pp_H_hh_CMS8_cache_e[2][CacheSize];
    mutable double ip_ex_pp_phi_hh_gagabb_CMS8_cache[2][CacheSize];
    mutable double ip_ex_pp_phi_hh_gagabb_CMS8_cache_e[2][CacheSize];
    mutable double ip_ex_pp_phi_hh_bbbb_CMS8_cache[2][CacheSize];
    mutable double ip_ex_pp_phi_hh_bbbb_CMS8_cache_e[2][CacheSize];
    mutable double ip_ex_bb_phi_bb_CMS8_cache[2][CacheSize];
    mutable double ip_ex_bb_phi_bb_CMS8_cache_e[2][CacheSize];
    mutable double ip_ex_gg_phi_tautau_CMS8_cache[2][CacheSize];
    mutable double ip_ex_gg_phi_tautau_CMS8_cache_e[2][CacheSize];
    mutable double ip_ex_bb_phi_tautau_CMS8_cache[2][CacheSize];
    mutable double ip_ex_bb_phi_tautau_CMS8_cache_e[2][CacheSize];
    mutable double ip_ex_gg_phi_gaga_CMS8_cache[2][CacheSize];
    mutable double ip_ex_gg_phi_gaga_CMS8_cache_e[2][CacheSize];
//            mutable double ip_ex_ggF_phi_gaga_CMS_cache_ep2[2][CacheSize];
//            mutable double ip_ex_ggF_phi_gaga_CMS_cache_em2[2][CacheSize];
    mutable double ip_ex_pp_A_Zga_llga_CMS8_cache[2][CacheSize];
    mutable double ip_ex_pp_A_Zga_llga_CMS8_cache_e[2][CacheSize];
    mutable double ip_ex_gg_H_hh_bbtautau_CMS8_cache[2][CacheSize];
    mutable double ip_ex_gg_H_hh_bbtautau_CMS8_cache_e[2][CacheSize];
    mutable double ip_ex_gg_A_hZ_tautaull_CMS8_cache[2][CacheSize];
    mutable double ip_ex_gg_A_hZ_tautaull_CMS8_cache_e[2][CacheSize];

    mutable double ip_ex_bb_phi_tt_ATLAS13_cache[2][CacheSize];
    mutable double ip_ex_bb_phi_tt_ATLAS13_cache_e[2][CacheSize];
    mutable double ip_ex_tt_phi_tt_ATLAS13_cache[2][CacheSize];
    mutable double ip_ex_tt_phi_tt_ATLAS13_cache_e[2][CacheSize];
    mutable double ip_ex_gg_phi_tautau_ATLAS13_cache[2][CacheSize];
    mutable double ip_ex_gg_phi_tautau_ATLAS13_cache_e[2][CacheSize];
    mutable double ip_ex_bb_phi_tautau_ATLAS13_cache[2][CacheSize];
    mutable double ip_ex_bb_phi_tautau_ATLAS13_cache_e[2][CacheSize];
    mutable double ip_ex_pp_phi_gaga_ATLAS13_cache[2][CacheSize];
    mutable double ip_ex_pp_phi_gaga_ATLAS13_cache_e[2][CacheSize];
    mutable double ip_ex_pp_phi_Zga_llga_ATLAS13_cache[2][CacheSize];
    mutable double ip_ex_pp_phi_Zga_llga_ATLAS13_cache_e[2][CacheSize];
    mutable double ip_ex_gg_H_ZZ_llnunu_ATLAS13_cache[2][CacheSize];
    mutable double ip_ex_gg_H_ZZ_llnunu_ATLAS13_cache_e[2][CacheSize];
    mutable double ip_ex_gg_H_ZZ_llll_ATLAS13_cache[2][CacheSize];
    mutable double ip_ex_gg_H_ZZ_llll_ATLAS13_cache_e[2][CacheSize];
    mutable double ip_ex_VBF_H_ZZ_llll_ATLAS13_cache[2][CacheSize];
    mutable double ip_ex_VBF_H_ZZ_llll_ATLAS13_cache_e[2][CacheSize];
    mutable double ip_ex_gg_H_ZZ_llqq_ATLAS13_cache[2][CacheSize];
    mutable double ip_ex_gg_H_ZZ_llqq_ATLAS13_cache_e[2][CacheSize];
    mutable double ip_ex_VBF_H_ZZ_llqq_ATLAS13_cache[2][CacheSize];
    mutable double ip_ex_VBF_H_ZZ_llqq_ATLAS13_cache_e[2][CacheSize];
    mutable double ip_ex_gg_H_ZZ_nunuqq_ATLAS13_cache[2][CacheSize];
    mutable double ip_ex_gg_H_ZZ_nunuqq_ATLAS13_cache_e[2][CacheSize];
    mutable double ip_ex_pp_H_ZZ_llqq_CMS13_cache[2][CacheSize];
    mutable double ip_ex_pp_H_ZZ_llqq_CMS13_cache_e[2][CacheSize];
    mutable double ip_ex_gg_H_WW_enumunu_ATLAS13_cache[2][CacheSize];
    mutable double ip_ex_gg_H_WW_enumunu_ATLAS13_cache_e[2][CacheSize];
    mutable double ip_ex_VBF_H_WW_enumunu_ATLAS13_cache[2][CacheSize];
    mutable double ip_ex_VBF_H_WW_enumunu_ATLAS13_cache_e[2][CacheSize];
    mutable double ip_ex_gg_H_WW_lnuqq_ATLAS13_cache[2][CacheSize];
    mutable double ip_ex_gg_H_WW_lnuqq_ATLAS13_cache_e[2][CacheSize];
    mutable double ip_ex_pp_H_hh_bbbb_ATLAS13_cache[2][CacheSize];
    mutable double ip_ex_pp_H_hh_bbbb_ATLAS13_cache_e[2][CacheSize];
    mutable double ip_ex_pp_H_hh_gagabb_ATLAS13_cache[2][CacheSize];
    mutable double ip_ex_pp_H_hh_gagabb_ATLAS13_cache_e[2][CacheSize];
    mutable double ip_ex_pp_H_hh_gagaWW_ATLAS13_cache[2][CacheSize];
    mutable double ip_ex_pp_H_hh_gagaWW_ATLAS13_cache_e[2][CacheSize];
    mutable double ip_ex_gg_A_Zh_Zbb_ATLAS13_cache[2][CacheSize];
    mutable double ip_ex_gg_A_Zh_Zbb_ATLAS13_cache_e[2][CacheSize];
    mutable double ip_ex_bb_A_Zh_Zbb_ATLAS13_cache[2][CacheSize];
    mutable double ip_ex_bb_A_Zh_Zbb_ATLAS13_cache_e[2][CacheSize];
    mutable double ip_ex_pp_phi_bb_CMS13_cache[2][CacheSize];
    mutable double ip_ex_pp_phi_bb_CMS13_cache_e[2][CacheSize];
    mutable double ip_ex_gg_phi_tautau_CMS13_cache[2][CacheSize];
    mutable double ip_ex_gg_phi_tautau_CMS13_cache_e[2][CacheSize];
    mutable double ip_ex_bb_phi_tautau_CMS13_cache[2][CacheSize];
    mutable double ip_ex_bb_phi_tautau_CMS13_cache_e[2][CacheSize];
    mutable double ip_ex_gg_phi_gaga_CMS13_cache[2][CacheSize];
    mutable double ip_ex_gg_phi_gaga_CMS13_cache_e[2][CacheSize];
    mutable double ip_ex_pp_phi_Zga_llga_CMS13_cache[2][CacheSize];
    mutable double ip_ex_pp_phi_Zga_llga_CMS13_cache_e[2][CacheSize];
    mutable double ip_ex_pp_phi_Zga_qqga_CMS13_cache[2][CacheSize];
    mutable double ip_ex_pp_phi_Zga_qqga_CMS13_cache_e[2][CacheSize];
    mutable double ip_ex_pp_H_ZZ_llll_CMS13_cache[2][CacheSize];
    mutable double ip_ex_pp_H_ZZ_llll_CMS13_cache_e[2][CacheSize];
    mutable double ip_ex_VBF_VH_H_ZZ_llll_CMS13_cache[2][CacheSize];
    mutable double ip_ex_VBF_VH_H_ZZ_llll_CMS13_cache_e[2][CacheSize];
    mutable double ip_ex_ggVV_H_WW_lnulnu_CMS13_cache[2][CacheSize];
    mutable double ip_ex_ggVV_H_WW_lnulnu_CMS13_cache_e[2][CacheSize];
    mutable double ip_ex_pp_H_hh_bbbb_CMS13_cache[2][CacheSize];
    mutable double ip_ex_pp_H_hh_bbbb_CMS13_cache_e[2][CacheSize];
    mutable double ip_ex_pp_H_hh_gagabb_CMS13_cache[2][CacheSize];
    mutable double ip_ex_pp_H_hh_gagabb_CMS13_cache_e[2][CacheSize];
    mutable double ip_ex_pp_H_hh_bbtautau_CMS13_cache[2][CacheSize];
    mutable double ip_ex_pp_H_hh_bbtautau_CMS13_cache_e[2][CacheSize];
    mutable double ip_ex_pp_H_hh_bbtautau1_CMS13_cache[2][CacheSize];
    mutable double ip_ex_pp_H_hh_bbtautau1_CMS13_cache_e[2][CacheSize];
    mutable double ip_ex_pp_H_hh_bblnulnu_CMS13_cache[2][CacheSize];
    mutable double ip_ex_pp_H_hh_bblnulnu_CMS13_cache_e[2][CacheSize];
    mutable double ip_ex_pp_H_hh_bbVV_CMS13_cache[2][CacheSize];
    mutable double ip_ex_pp_H_hh_bbVV_CMS13_cache_e[2][CacheSize];

    mutable double ie1c[2][CacheSize];
    mutable double ie1ec[2][CacheSize];
    mutable double ie2c[2][CacheSize];
    mutable double ie2ec[2][CacheSize];
    mutable double ie3c[2][CacheSize];
    mutable double ie3ec[2][CacheSize];
    mutable double ie4c[2][CacheSize];
    mutable double ie4ec[2][CacheSize];
    mutable double ie5c[2][CacheSize];
    mutable double ie5ec[2][CacheSize];
    mutable double ie6c[2][CacheSize];
    mutable double ie6ec[2][CacheSize];
    mutable double ie7c[2][CacheSize];
    mutable double ie7ec[2][CacheSize];
    mutable double ie8c[2][CacheSize];
    mutable double ie8ec[2][CacheSize];
    mutable double ie9c[2][CacheSize];
    mutable double ie9ec[2][CacheSize];
    mutable double ie10c[2][CacheSize];
    mutable double ie10ec[2][CacheSize];
    mutable double ie11c[2][CacheSize];
    mutable double ie11ec[2][CacheSize];
    mutable double ie12c[2][CacheSize];
    mutable double ie12ec[2][CacheSize];
    mutable double ie13c[2][CacheSize];
    mutable double ie13ec[2][CacheSize];
    mutable double ie14c[2][CacheSize];
    mutable double ie14ec[2][CacheSize];
    mutable double ie15c[2][CacheSize];
    mutable double ie15ec[2][CacheSize];
    mutable double ie16c[2][CacheSize];
    mutable double ie16ec[2][CacheSize];
    mutable double ie17c[2][CacheSize];
    mutable double ie17ec[2][CacheSize];
    mutable double ie18c[2][CacheSize];
    mutable double ie18ec[2][CacheSize];
    mutable double ie19c[2][CacheSize];
    mutable double ie19ec[2][CacheSize];
    mutable double ie20c[2][CacheSize];
    mutable double ie20ec[2][CacheSize];
    mutable double ie21c[2][CacheSize];
    mutable double ie21ec[2][CacheSize];
    mutable double ie22c[2][CacheSize];
    mutable double ie22ec[2][CacheSize];
    mutable double ie23c[2][CacheSize];
    mutable double ie23ec[2][CacheSize];
    mutable double ie24c[2][CacheSize];
    mutable double ie24ec[2][CacheSize];
    mutable double ie25c[2][CacheSize];
    mutable double ie25ec[2][CacheSize];
    mutable double ie26c[2][CacheSize];
    mutable double ie26ec[2][CacheSize];
    mutable double ie27c[2][CacheSize];
    mutable double ie27ec[2][CacheSize];
    mutable double ie28c[2][CacheSize];
    mutable double ie28ec[2][CacheSize];
    mutable double ie29c[2][CacheSize];
    mutable double ie29ec[2][CacheSize];
    mutable double ie30c[2][CacheSize];
    mutable double ie30ec[2][CacheSize];
    mutable double ie31c[2][CacheSize];
    mutable double ie31ec[2][CacheSize];
    mutable double ie32c[2][CacheSize];
    mutable double ie32ec[2][CacheSize];
    mutable double ie33c[2][CacheSize];
    mutable double ie33ec[2][CacheSize];
    mutable double ie34c[2][CacheSize];
    mutable double ie34ec[2][CacheSize];
    mutable double ie35c[2][CacheSize];
    mutable double ie35ec[2][CacheSize];
    mutable double ie36c[2][CacheSize];
    mutable double ie36ec[2][CacheSize];
    mutable double ie37c[2][CacheSize];
    mutable double ie37ec[2][CacheSize];
    mutable double ie38c[2][CacheSize];
    mutable double ie38ec[2][CacheSize];
    mutable double ie39c[2][CacheSize];
    mutable double ie39ec[2][CacheSize];
    mutable double ie40c[2][CacheSize];
    mutable double ie40ec[2][CacheSize];

    mutable double ip_ex_pp_Hpm_taunu_ATLAS8_cache[2][CacheSize];
    mutable double ip_ex_pp_Hpm_taunu_ATLAS8_cache_e[2][CacheSize];
    mutable double ip_ex_pp_Hp_taunu_CMS8_cache[2][CacheSize];
    mutable double ip_ex_pp_Hp_taunu_CMS8_cache_e[2][CacheSize];
    mutable double ip_ex_pp_Hp_tb_ATLAS8_cache[2][CacheSize];
    mutable double ip_ex_pp_Hp_tb_ATLAS8_cache_e[2][CacheSize];
    mutable double ip_ex_pp_Hp_tb_CMS8_cache[2][CacheSize];
    mutable double ip_ex_pp_Hp_tb_CMS8_cache_e[2][CacheSize];
    mutable double ip_ex_pp_Hpm_taunu_ATLAS13_cache[2][CacheSize];
    mutable double ip_ex_pp_Hpm_taunu_ATLAS13_cache_e[2][CacheSize];
    mutable double ip_ex_pp_Hpm_taunu_CMS13_cache[2][CacheSize];
    mutable double ip_ex_pp_Hpm_taunu_CMS13_cache_e[2][CacheSize];
    mutable double ip_ex_pp_Hp_tb_ATLAS13_1_cache[2][CacheSize];
    mutable double ip_ex_pp_Hp_tb_ATLAS13_1_cache_e[2][CacheSize];
    mutable double ip_ex_pp_Hp_tb_ATLAS13_2_cache[2][CacheSize];
    mutable double ip_ex_pp_Hp_tb_ATLAS13_2_cache_e[2][CacheSize];

    mutable double ip_ex_bsgamma_cache[3][CacheSize];

    mutable double ghHpHm_cache[7][CacheSize];
    mutable double g_HH_HpHm_cache[7][CacheSize];

    mutable gslpp::complex I_h_U_cache[5][CacheSize];
    mutable gslpp::complex I_HH_U_cache[4][CacheSize];
    mutable gslpp::complex I_A_U_cache[4][CacheSize];
    mutable gslpp::complex I_h_D_cache[5][CacheSize];
    mutable gslpp::complex I_HH_D_cache[4][CacheSize];
    mutable gslpp::complex I_A_D_cache[4][CacheSize];
    mutable gslpp::complex I_h_L_cache[5][CacheSize];
    mutable gslpp::complex I_HH_L_cache[4][CacheSize];
    mutable gslpp::complex I_A_L_cache[4][CacheSize];
    mutable gslpp::complex I_H_W_cache[3][CacheSize];
    mutable gslpp::complex I_H_Hp_cache[3][CacheSize];

    mutable gslpp::complex A_h_U_cache[7][CacheSize];
    mutable gslpp::complex A_HH_U_cache[6][CacheSize];
    mutable gslpp::complex A_A_U_cache[6][CacheSize];
    mutable gslpp::complex A_h_D_cache[7][CacheSize];
    mutable gslpp::complex A_HH_D_cache[6][CacheSize];
    mutable gslpp::complex A_A_D_cache[6][CacheSize];
    mutable gslpp::complex A_h_L_cache[7][CacheSize];
    mutable gslpp::complex A_HH_L_cache[6][CacheSize];
    mutable gslpp::complex A_A_L_cache[6][CacheSize];
    mutable gslpp::complex A_H_W_cache[5][CacheSize];
    mutable gslpp::complex A_H_Hp_cache[5][CacheSize];

    mutable double KaellenFunction_cache[4][CacheSize];

    /**
     * @brief f function for the gamma gamma coupling to h, H and A
     * @return @f$f(x)@f$
     * @details The definition can be found in (2.19) of @cite Gunion:1989we.
     */
    gslpp::complex f_func(const double x) const;

    /**
     * @brief g function for the Int1 function
     * @return @f$g(x)@f$
     * @details The definition can be found in (2.24) of @cite Gunion:1989we.
     */
    gslpp::complex g_func(const double x) const;

    /**
     * @brief @f$I_1@f$ function for Z gamma coupling to h, H and A
     * @return @f$I_1(\tau,\lambda)@f$
     * @details The definition can be found in (2.24) of @cite Gunion:1989we.
     */
    gslpp::complex Int1(const double tau, const double lambda) const;

    /**
     * @brief @f$I_2@f$ function for Z gamma coupling to h, H and A
     * @return @f$I_2(\tau,\lambda)@f$
     * @details The definition can be found in (2.24) of @cite Gunion:1989we.
     */
    gslpp::complex Int2(const double tau, const double lambda) const;

    /**
     * @brief Heaviside @f$\Theta@f$ function
     * @return @f$\Theta(x)@f$
     * @details Gives 1 for @f$x\geq 0@f$ and 0 for @f$x<0@f$.
     */
    int HSTheta (const double x) const;
};

#endif	/* THDMCACHE_H */
