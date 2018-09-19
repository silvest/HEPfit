/* 
 * Copyright (C) 2017 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THDMWCACHE_H
#define	THDMWCACHE_H

#include "THDMW.h"
#include "RunnerTHDMW.h"

/**
 * @class THDMWcache
 * @ingroup THDMW
 * @brief A class for the caching of some %THDMW objects.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class THDMWcache {
    
public:

    /**
     * @brief THDMWcache constructor.
     * @details Reads all the tables values and stores them in the memory.
     */
    THDMWcache(const StandardModel& SM_i);

    /**
     * @brief THDMWcache destructor.
     */
    ~THDMWcache();
    
    void updateCache();
    double setOtherParameters();

    double Q_cutoff;
//    double g1_at_Q;
//    double g2_at_Q;
//    double g3_at_Q;
//    double Ytop_at_Q;
//    double Ybottom1_at_Q;
//    double Ybottom2_at_Q;
//    double Ytau1_at_Q;
//    double Ytau2_at_Q;
    double lambda1_at_Q;
    double lambda2_at_Q;
    double lambda3_at_Q;
    double lambda4_at_Q;
    double mu1_at_Q;
    double mu3_at_Q;
    double mu4_at_Q;
    double nu1_at_Q;
    double omega1_at_Q;
    double kappa1_at_Q;
    double nu2_at_Q;
    double omega2_at_Q;
    double kappa2_at_Q;
    double nu4_at_Q;
    double omega4_at_Q;
    double nu3_at_Q;
    double nu5_at_Q;
    double mu2_at_Q;
    double mu5_at_Q;
    double mu6_at_Q;
    double m12sq;
    double m11sq;
    double m22sq;
    double mhsq;
    double mHsq;
    double mAsq;
    double mSRsq;
    double mSIsq;
    double mHpsq;
    double mSpsq;


    
     /**
     * @brief Cross section times branching ratio for the process @f$pp\to Sr\to t\bar t@f$ at the LHC with 13 TeV.
     * @return @f$\sigma^{\text{THDMW}}_{pp\to Sr}\cdot BR^{\text{THDMW}}(Sr\to t\bar t)@f$
     */
    double pp_Sr_tt_TH13;
    
    
    
    
    
    double THoEX_pp_Sr_tt;
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    double RpepsTHDMW;
    gslpp::vector<gslpp::complex> unitarityeigenvalues;
    gslpp::vector<gslpp::complex> NLOunitarityeigenvalues;

    double rh_QuQu, rh_VV, rh_gg, rh_QdQd, rh_ll, rh_gaga, rh_Zga;
    double sumModBRs, Gamma_h, THDM_BR_h_bb, THDM_BR_h_gaga, THDM_BR_h_tautau, THDM_BR_h_WW, THDM_BR_h_ZZ;
    
protected:

private:

    const THDMW * myTHDMW;
    RunnerTHDMW * myRunnerTHDMW;

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

    
    void computeHHlimits();
    
    
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
    
//    mutable double ip_cs_ppto2Sto4t_13_cache[5][CacheSize];
    mutable double ip_th_pp_Sr_tt_cache[5][CacheSize];
    mutable double ip_ex_pp_phi_hh_bbbb_CMS8_cache[2][CacheSize];
    mutable double ip_ex_pp_phi_hh_bbbb_CMS8_cache_e[2][CacheSize];
    mutable double ip_ex_bb_phi_bb_CMS8_cache[2][CacheSize];
    mutable double ip_ex_bb_phi_bb_CMS8_cache_e[2][CacheSize];
    mutable double ip_ex_gg_phi_tt_ATLAS8_cache[2][CacheSize];
    mutable double ip_ex_gg_phi_tt_ATLAS8_cache_e[2][CacheSize];
    mutable double ip_ex_bb_phi_tt_ATLAS13_cache[2][CacheSize];
    mutable double ip_ex_bb_phi_tt_ATLAS13_cache_e[2][CacheSize];
    mutable double ip_ex_tt_phi_tt_ATLAS13_cache[2][CacheSize];
    mutable double ip_ex_tt_phi_tt_ATLAS13_cache_e[2][CacheSize];
    mutable double ip_ex_pp_H_hh_bbbb_ATLAS13_cache[2][CacheSize];
    mutable double ip_ex_pp_H_hh_bbbb_ATLAS13_cache_e[2][CacheSize];
    mutable double ip_ex_pp_phi_bb_CMS13_cache[2][CacheSize];
    mutable double ip_ex_pp_phi_bb_CMS13_cache_e[2][CacheSize];
    mutable double ip_ex_pp_H_hh_bbbb_CMS13_cache[2][CacheSize];
    mutable double ip_ex_pp_H_hh_bbbb_CMS13_cache_e[2][CacheSize];
    mutable double ip_ex_pp_Hpm_tb_ATLAS8_cache[2][CacheSize];
    mutable double ip_ex_pp_Hpm_tb_ATLAS8_cache_e[2][CacheSize];
    mutable double ip_ex_pp_Hp_tb_CMS8_cache[2][CacheSize];
    mutable double ip_ex_pp_Hp_tb_CMS8_cache_e[2][CacheSize];
    mutable double ip_ex_pp_Hp_tb_ATLAS13_1_cache[2][CacheSize];
    mutable double ip_ex_pp_Hp_tb_ATLAS13_1_cache_e[2][CacheSize];
    mutable double ip_ex_pp_Hp_tb_ATLAS13_2_cache[2][CacheSize];
    mutable double ip_ex_pp_Hp_tb_ATLAS13_2_cache_e[2][CacheSize];
    mutable double ip_ex_ggF_H_hh_bbbb_CMS13_cache[2][CacheSize];
    mutable double ip_ex_ggF_H_hh_bbbb_CMS13_cache_e[2][CacheSize];
    mutable double ip_ex_pp_Gkk_tt_ATLAS13_cache[2][CacheSize];

    
    
    
    
        /**
     * @brief Fills all required arrays with the values read from the tables.
     */
    void read();
    
        /**
     * @brief This function reads values from a table and returns them as an array.
     * @return the tabled values
     */
    gslpp::matrix<double> readTable(std::string filename, int rowN, int colN);
    
    



    /**
     * @brief ATLAS observed @f$95\%@f$ upper cross section limits at 8 TeV, depending on the Higgs mass.
     */
    gslpp::matrix<double>  ATLAS8_gg_phi_tt;

    /**
     * @brief ATLAS expected @f$95\%@f$ upper cross section limits at 8 TeV, depending on the Higgs mass.
     */
    gslpp::matrix<double>   ATLAS8_gg_phi_tt_e;



    /**
     * @brief CMS observed @f$95\%@f$ upper cross section limits at 8 TeV, depending on the Higgs mass.
     */
    gslpp::matrix<double>  CMS8_pp_H_hh_bbbb, CMS8_bb_phi_bb;

    /**
     * @brief CMS expected @f$95\%@f$ upper cross section limits at 8 TeV, depending on the Higgs mass.
     */
    gslpp::matrix<double>  CMS8_pp_H_hh_bbbb_e, CMS8_bb_phi_bb_e;

//    gslpp::matrix<double> CMS_ggF_phi_gaga_ep2, CMS_ggF_phi_gaga_em2;

    /**
     * @brief ATLAS observed @f$95\%@f$ upper cross section limits at 13 TeV, depending on the Higgs mass.
     */
    gslpp::matrix<double> ATLAS13_bb_phi_tt, ATLAS13_tt_phi_tt, ATLAS13_pp_H_hh_bbbb;

    /**
     * @brief ATLAS expected @f$95\%@f$ upper cross section limits at 13 TeV, depending on the Higgs mass.
     */
    gslpp::matrix<double> ATLAS13_bb_phi_tt_e, ATLAS13_tt_phi_tt_e, ATLAS13_pp_H_hh_bbbb_e;

    /**
     * @brief CMS observed @f$95\%@f$ upper cross section limits at 13 TeV, depending on the Higgs mass.
     */
    gslpp::matrix<double> CMS13_pp_phi_bb, CMS13_pp_H_hh_bbbb, CMS13_ggF_H_hh_bbbb;

    /**
     * @brief CMS expected @f$95\%@f$ upper cross section limits at 13 TeV, depending on the Higgs mass.
     */
    gslpp::matrix<double> CMS13_pp_phi_bb_e, CMS13_pp_H_hh_bbbb_e, CMS13_ggF_H_hh_bbbb_e;


    /**
     * @brief ATLAS observed @f$95\%@f$ upper cross section limits at 8 TeV, depending on the charged Higgs mass.
     */
    gslpp::matrix<double> ATLAS8_pp_Hpm_tb;

    /**
     * @brief ATLAS expected @f$95\%@f$ upper cross section limits at 8 TeV, depending on the charged Higgs mass.
     */
    gslpp::matrix<double> ATLAS8_pp_Hpm_tb_e;

    /**
     * @brief CMS observed @f$95\%@f$ upper cross section limits at 8 TeV, depending on the charged Higgs mass.
     */
    gslpp::matrix<double> CMS8_pp_Hp_tb;

    /**
     * @brief CMS expected @f$95\%@f$ upper cross section limits at 8 TeV, depending on the charged Higgs mass.
     */
    gslpp::matrix<double> CMS8_pp_Hp_tb_e;

    /**
     * @brief ATLAS observed @f$95\%@f$ upper cross section limits at 13 TeV, depending on the charged Higgs mass.
     */
    gslpp::matrix<double> ATLAS13_pp_Hp_tb1, ATLAS13_pp_Hp_tb2;

    /**
     * @brief ATLAS expected @f$95\%@f$ upper cross section limits at 13 TeV, depending on the charged Higgs mass.
     */
    gslpp::matrix<double> ATLAS13_pp_Hp_tb1_e, ATLAS13_pp_Hp_tb2_e;


        /**
     * @brief ATLAS expected @f$95\%@f$ upper cross section limits at 13 TeV, depending on the Kaluza-Klein Graviton mass. Process pp -> Gkk -> t tbar.
     */
    gslpp::matrix<double> ATLAS_13_pp_Gkk_tt;
    
    
    
        /**
     * @brief Table for xsection times branching ratio for p p -> Sr -> t tbar generated with Madgraph
     */
    gslpp::matrix<double> MadGraph_pp_Sr_tt;
    
    /**
     * @brief @f$b\to s \gamma@f$ table, depending on logtb and the logarithm of the charged Higgs mass.
     */
    gslpp::matrix<double> arraybsgamma;
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    
    
    
    
    

    
    
    
    
    
    
    
    
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
     * @brief Interpolating function for the observed CMS upper limit on a scalar resonance decaying to a b quark pair.
     * @return @f$[\sigma_{pp\to \phi}\cdot BR(\phi \to b\bar b)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-PAS-HIG-16-025, Figure 5 @cite CMS:2016ncz.
     */
    double ip_ex_pp_phi_bb_CMS13(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a scalar resonance decaying to a b quark pair.
     * @return @f$[\sigma_{pp\to \phi}\cdot BR(\phi \to b\bar b)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-PAS-HIG-16-025, Figure 5 @cite CMS:2016ncz.
     */
    double ip_ex_pp_phi_bb_CMS13_e(double mass);

       /**
     * @brief Interpolating function for the observed CMS upper limit on a scalar resonance decaying to two @f$h@f$ bosons which further decay to four b quarks.
     * @return @f$[\sigma_{pp\to H}\cdot BR(H\to hh\to b\bar b b\bar b)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-PAS-HIG-16-002, Figure 7 @cite CMS:2016tlj.
     */
    double ip_ex_pp_H_hh_bbbb_CMS13(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a scalar resonance decaying to two @f$h@f$ bosons which further decay to four b quarks.
     * @return @f$[\sigma_{pp\to H}\cdot BR(H\to hh\to b\bar b b\bar b)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-PAS-HIG-16-002, Figure 7 @cite CMS:2016tlj.
     */
    double ip_ex_pp_H_hh_bbbb_CMS13_e(double mass);

    /**
     * @brief Interpolating function for the observed ATLAS upper limit on a singly charged scalar resonance decaying to a @f$t@f$ quark and a @f$b@f$ quark.
     * @return @f$[2\sigma_{pp\to H^+}\cdot BR(H^+\to tb)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1512.03704, Figure 6 @cite Aad:2015typ.
     */
    double ip_ex_pp_Hpm_tb_ATLAS8(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on a singly charged scalar resonance decaying to a @f$t@f$ quark and a @f$b@f$ quark.
     * @return @f$[2\sigma_{pp\to H^+}\cdot BR(H^+\to tb)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1512.03704, Figure 6 @cite Aad:2015typ.
     */
    double ip_ex_pp_Hpm_tb_ATLAS8_e(double mass);

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
     * @brief Interpolating function for the observed CMS upper limit on a scalar resonance decaying to two @f$h@f$ bosons which further decay to four b quarks.
     * @return @f$[\sigma_{gg\to H}\cdot BR(H\to hh\to b\bar b b\bar b)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-PAS-B2G-16-026, Figure 9 left @cite CMS:2017gxe.
     */
    double ip_ex_ggF_H_hh_bbbb_CMS13(double mass);

    /**
     * @brief Interpolating function for the expected CMS upper limit on a scalar resonance decaying to two @f$h@f$ bosons which further decay to four b quarks.
     * @return @f$[\sigma_{gg\to H}\cdot BR(H\to hh\to b\bar b b\bar b)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-PAS-B2G-16-026, Figure 9 left @cite CMS:2017gxe.
     */
    double ip_ex_ggF_H_hh_bbbb_CMS13_e(double mass);

    /**
     * @brief Interpolating function for the expected ATLAS upper limit on pp -> Gkk (Kaluza-Klein graviton) -> t tbar
     * @return xsection times branching ratio
     * @details ATLAS arXiv:1804.01939, Fig. 16, 36.1 fb^-1, 13 TeV
     */
    double ip_ex_pp_Gkk_tt_ATLAS13(double mass);    
    
    /**
     * @brief Interpolating function for the theoretical value of p p -> Sr -> t tbar
     * @return xsection times branching ratio of pp -> Sr -> t tbar
     * @details Generated with MadGraph
     */
    double ip_th_pp_Sr_tt(double etaD, double etaU, double Lambda4, double mass);
    
    
    
    
    
    
    
    
    

    
    
    
    
    
    
    
    gslpp::complex f_func(const double x) const;
    gslpp::complex g_func(const double x) const;

    gslpp::complex Int1(const double tau, const double lambda) const;
    gslpp::complex Int2(const double tau, const double lambda) const;

    void runTHDMWparameters();

    void computeUnitarity();
    gslpp::vector<gslpp::complex> betaeigenvalues;

    void computeSignalStrengthQuantities();

    std::string THDMWmodel;
    double Q_THDMW;
    double MZ;
    double vev;
    double tanb;
    double sinb;
    double cosb;
    double bma;
    double sina;
    double cosa;
    double lambda1;
    double lambda2;
    double lambda3;
    double lambda4;
    double lambda5;
    double mSsq;
    double mu1;
    double mu2;
    double mu3;
    double mu4;
    double mu5;
    double mu6;
    double nu1;
    double nu2;
    double nu3;
    double nu4;
    double nu5;
    double omega1;
    double omega2;
    double omega3;
    double omega4;
    double kappa1;
    double kappa2;
    double kappa3;
    double etaU;
    double etaD;
    
    
        /**
     * @brief Linearly interpolates a table with one parameter dimension.
     * @return the interpolated value
     */
    double interpolate (gslpp::matrix<double> arrayTab, double x);
        
    /**
     * @brief Linearly interpolates a table with four parameter dimensions.
     * @return the interpolated value
     */
    double interpolate4D (gslpp::matrix<double> arrayTab, double x, double y, double z, double v);

    /**
     * @brief Interpolating function for the H production cross section via gluon-gluon fusion at 8 TeV.
     * @return @f$\sigma(gg\to H)@f$
     */
    double ip_cs_ppto2Sto4t_13(double etaD, double etaU, double THDMW_nu4, double mSR);
    
    
};

#endif	/* THDMWCACHE_H */
