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
//#include "THDM.h"

#include <stdexcept>
#include "gslpp.h"

/**
 * @class THDMcache
 * @ingroup THDM
 * @brief A class for the caching of some THDM objects.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details At the moment only the Passarino-Veltman functions for STU are cached.
 * The tables are also read here.
 */
class THDMcache {
    
public:

    /**
     * @brief THDMcache constructor.
     * @details Reads all the tables values and stores them in the memory.
     */
    THDMcache();
    
    /**
     * @brief Cache size.
     * @details Determines the size of the cache. If it is set to 5, the cache will remember the last five function calls and store their results.
     */
    static const int CacheSize = 5;

    /**
     * @brief Check whether for the latest set of parameters a value is in the cache.
     */
    int CacheCheck(const gslpp::complex cache[][CacheSize],
                   const int NumPar, const double params[]) const;

    /**
     * @brief Adds a new result and its parameters into the cache.
     * @details The new values are added on top. The oldest set on the stack is deleted.
     */
    void CacheShift(gslpp::complex cache[][CacheSize], const int NumPar,
                    const double params[], const gslpp::complex newResult) const; 

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
     * @brief LHC production cross section percentage contributions (derived from the cross sections provided by the LHC Higgs Cross Section Working Group), depending on the Higgs mass.
     */
    gslpp::matrix<double> pc_ggF, pc_VBF, pc_WH, pc_ZH, pc_ttH;

    /**
     * @brief Total SM decay width (obtained with HDECAY 6.10), depending on the Higgs mass.
     */
    gslpp::matrix<double> GammaHtotSM;

    /**
     * @brief HIGLU v4.00 cross sections, depending on the Higgs mass.
     */
    gslpp::matrix<double> cs_ggH, cs_ggH_tt, cs_ggH_bb, cs_ggA, cs_ggA_tt, cs_ggA_bb;

    /**
     * @brief HiggsBounds 4.1.0 @f$b\bar b \to H@f$ cross sections, depending on the Higgs mass.
     */
    gslpp::matrix<double> cs_bbFtoHP;

    /**
     * @brief ATLAS @f$95\%@f$ upper cross section limits, depending on the Higgs mass.
     */
    gslpp::matrix<double> ATLAS_ggF_phi_gaga, ATLAS_ggF_phi_tautau, ATLAS_bbF_phi_tautau, ATLAS_ggF_A_hZ_tautauZ, ATLAS_ggF_A_hZ_bbZ, ATLAS_pp_phi_tt, ATLAS_ggF_H_WW, ATLAS_VBF_H_WW, ATLAS_ggF_H_hh;

    /**
     * @brief CMS @f$95\%@f$ upper signal strength limits, depending on the Higgs mass.
     */
    gslpp::matrix<double> CMS_pp_H_ZZ;

    /**
     * @brief CMS @f$95\%@f$ upper cross section limits, depending on the Higgs mass.
     */
    gslpp::matrix<double> CMS_ggF_A_hZ_bbll, CMS_pp_H_hh_gagabb, CMS_pp_H_hh_bbbb, CMS_bbF_phi_bb, CMS_ggF_phi_tautau, CMS_bbF_phi_tautau, CMS_ggF_phi_gaga, CMS_ggF_H_hh_bbtautau, CMS_ggF_A_hZ_tautaull;

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
     * @brief Interpolating function for the SM percentage contribution of gluon-gluon fusion to the total Higgs production cross section.
     * @return @f$pc(gg\to H)@f$
     */
    double ip_pc_ggFtoHP(double mass);

    /**
     * @brief Interpolating function for the SM percentage contribution of vector boson fusion to the total Higgs production cross section.
     * @return @f$pc(VV\to H)@f$
     */
    double ip_pc_VBFtoHP(double mass);

    /**
     * @brief Interpolating function for the SM percentage contribution of @f$W@f$ Higgsstrahlung to the total Higgs production cross section.
     * @return @f$pc(W\to WH)@f$
     */
    double ip_pc_WHP_HP(double mass);

    /**
     * @brief Interpolating function for the SM percentage contribution of @f$Z@f$ Higgsstrahlung to the total Higgs production cross section.
     * @return @f$pc(Z\to ZH)@f$
     */
    double ip_pc_ZHP_HP(double mass);

    /**
     * @brief Interpolating function for the SM percentage contribution of @f$t\bar t@f$ associated production to the total Higgs production cross section.
     * @return @f$pc(t\bar t\to H)@f$
     */
    double ip_pc_ttFtoHP(double mass);

    /**
     * @brief Interpolating function for the total SM Higgs decay width.
     * @return @f$\Gamma^{\text{tot}}_H@f$
     */
    double ip_GammaHPtotSM(double mass);

    /**
     * @brief Interpolating function for the SM Higgs production cross section via gluon-gluon fusion.
     * @return @f$\sigma(gg\to H)@f$
     */
    double ip_cs_ggFtoHP(double mass);

    /**
     * @brief Interpolating function for the SM Higgs production cross section via gluon-gluon fusion (top-loop only).
     * @return @f$\sigma_t(gg\to H)@f$
     */
    double ip_cs_ggHP_tt(double mass);

    /**
     * @brief Interpolating function for the SM Higgs production cross section via gluon-gluon fusion (bottom-loop only).
     * @return @f$\sigma_b(gg\to H)@f$
     */
    double ip_cs_ggHP_bb(double mass);

    /**
     * @brief Interpolating function for the production cross section of a pseudoscalar via gluon-gluon fusion.
     * @return @f$\sigma(gg\to A)@f$
     */
    double ip_cs_ggA(double mass);

    /**
     * @brief Interpolating function for the production cross section of a pseudoscalar via gluon-gluon fusion (top-loop only).
     * @return @f$\sigma_t(gg\to A)@f$
     */
    double ip_cs_ggA_tt(double mass);

    /**
     * @brief Interpolating function for the production cross section of a pseudoscalar via gluon-gluon fusion (bottom-loop only).
     * @return @f$\sigma_b(gg\to A)@f$
     */
    double ip_cs_ggA_bb(double mass);

    /**
     * @brief Interpolating function for the bottom quark associated production cross section of a Higgs.
     * @return @f$\sigma(b\bar b\to H)@f$
     */
    double ip_cs_bbFtoHP(double mass);

    /**
     * @brief Interpolating function for the ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two photons.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to \gamma \gamma)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1407.6583, Figure 4 @cite Aad:2014ioa.
     */
    double ip_ex_ggF_phi_gaga_ATLAS(double mass);

    /**
     * @brief Interpolating function for the ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two tau leptons.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to \tau \tau)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1409.6064, Figure 11a @cite Aad:2014vgg.
     */
    double ip_ex_ggF_phi_tautau_ATLAS(double mass);

    /**
     * @brief Interpolating function for the ATLAS upper limit on a bottom quark produced scalar resonance decaying to two tau leptons.
     * @return @f$[\sigma_{bb\to \phi}\cdot BR(\phi\to \tau \tau)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1409.6064, Figure 11b @cite Aad:2014vgg.
     */
    double ip_ex_bbF_phi_tautau_ATLAS(double mass);

    /**
     * @brief Interpolating function for the ATLAS upper limit on a gluon-gluon produced pseudoscalar resonance decaying to @f$hZ@f$ of which the Higgs further decays to a @f$\tau@f$ lepton pair.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to hZ\to \tau \tau Z)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1502.04478, Figure 3a @cite Aad:2015wra.
     */
    double ip_ex_ggF_A_hZ_tautauZ_ATLAS(double mass);

    /**
     * @brief Interpolating function for the ATLAS upper limit on a gluon-gluon produced pseudoscalar resonance decaying to @f$hZ@f$ of which the Higgs further decays to a bottom quark pair.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to hZ\to b\bar b Z)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1502.04478, Figure 3b @cite Aad:2015wra.
     */
    double ip_ex_ggF_A_hZ_bbZ_ATLAS(double mass);

    /**
     * @brief Interpolating function for the ATLAS upper limit on scalar resonance decaying to a top quark pair.
     * @return @f$[\sigma_{pp\to \phi}\cdot BR(\phi\to t\bar t)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1505.07018, Figure 11d @cite Aad:2015fna.
     */
    double ip_ex_pp_phi_tt_ATLAS(double mass);

    /**
     * @brief Interpolating function for the ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two @f$W@f$ bosons.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to WW)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1509.00389, Figure 13, left @cite Aad:2015agg.
     */
    double ip_ex_ggF_H_WW_ATLAS(double mass);

    /**
     * @brief Interpolating function for the ATLAS upper limit on a vector boson fusion produced scalar resonance decaying to two @f$W@f$ bosons.
     * @return @f$[\sigma_{VV\to \phi}\cdot BR(\phi\to WW)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1509.00389, Figure 13, right @cite Aad:2015agg.
     */
    double ip_ex_VBF_H_WW_ATLAS(double mass);

    /**
     * @brief Interpolating function for the ATLAS upper limit on a gluon-gluon produced scalar resonance decaying to two @f$h@f$ bosons.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to hh)]_{\text{ATLAS,95\%}}@f$
     * @details Taken from arXiv:1509.04670, Figure 6 @cite Aad:2015xja.
     */
    double ip_ex_ggF_H_hh_ATLAS(double mass);

    /**
     * @brief Interpolating function for the CMS upper limit on a scalar resonance decaying to two @f$Z@f$ bosons.
     * @return @f$[\mu_H(H\to ZZ)]_{\text{CMS,95\%}}@f$
     * @details Taken from arXiv:1504.00936, Figure 7, bottom right @cite Khachatryan:2015cwa.
     */
    double ip_ex_pp_H_ZZ_CMS(double mass);

    /**
     * @brief Interpolating function for the CMS upper limit on a gluon-gluon produced pseudoscalar resonance decaying to @f$hZ@f$ which further decay to a bottom quark pair and a light lepton pair.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to hZ\to b\bar b \ell \ell)]_{\text{CMS,95\%}}@f$
     * @details Taken from arXiv:1504.04710, Figure 3 @cite Khachatryan:2015lba.
     */
    double ip_ex_ggF_A_hZ_bbll_CMS(double mass);

    /**
     * @brief Interpolating function for the CMS upper limit on a scalar resonance decaying to two @f$h@f$ bosons which further decay to a photon pair and a bottom quark pair.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to hh\to \gamma \gamma b\bar b)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-PAS-HIG-13-032, Figure 8 @cite CMS:2014ipa.
     */
    double ip_ex_pp_phi_hh_gagabb_CMS(double mass);

    /**
     * @brief Interpolating function for the CMS upper limit on a scalar resonance decaying to two @f$h@f$ bosons which further decay to two bottom quark pairs.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to hh\to b\bar b b\bar b)]_{\text{CMS,95\%}}@f$
     * @details Taken from arXiv:1503.04114, Figure 5, left @cite Khachatryan:2015yea.
     */
    double ip_ex_pp_phi_hh_bbbb_CMS(double mass);

    /**
     * @brief Interpolating function for the CMS upper limit on a bottom quark produced scalar resonance decaying to two bottom quarks.
     * @return @f$[\sigma_{bb\to \phi}\cdot BR(\phi\to b\bar b)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-HIG-14-017, Figure 6 @cite Khachatryan:2015tra.
     */
    double ip_ex_bbF_phi_bb_CMS(double mass);

    /**
     * @brief Interpolating function for the CMS upper limit on a gluon-gluon produced scalar resonance decaying to two tau leptons.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to \tau \tau)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-PAS-HIG-14-029, Figure 10-a @cite CMS:2015mca.
     */
    double ip_ex_ggF_phi_tautau_CMS(double mass);

    /**
     * @brief Interpolating function for the CMS upper limit on a bottom quark produced scalar resonance decaying to two tau leptons.
     * @return @f$[\sigma_{bb\to \phi}\cdot BR(\phi\to \tau \tau)]_{\text{CMS,95\%}}@f$
     * @details Taken from CMS-PAS-HIG-14-029, Figure 10-b @cite CMS:2015mca.
     */
    double ip_ex_bbF_phi_tautau_CMS(double mass);

    /**
     * @brief Interpolating function for the CMS upper limit on a gluon-gluon produced scalar resonance decaying to two photons.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to \gamma \gamma)]_{\text{CMS,95\%}}@f$
     * @details Taken from arXiv:1506.02301, Figure 7, left @cite Khachatryan:2015qba.
     */
    double ip_ex_ggF_phi_gaga_CMS(double mass);

    /**
     * @brief Interpolating function for the CMS upper limit on a gluon-gluon produced scalar resonance decaying to two @f$h@f$ bosons which further decay to a bottom quark pair and a @f$\tau@f$ lepton pair.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to hh\to b\bar b \tau \tau)]_{\text{CMS,95\%}}@f$
     * @details Taken from arXiv:1510.01181, Figure 8, bottom right @cite Khachatryan:2015tha.
     */
    double ip_ex_ggF_H_hh_bbtautau_CMS(double mass);

    /**
     * @brief Interpolating function for the CMS upper limit on a gluon-gluon produced pseudoscalar resonance decaying to @f$hZ@f$ which further decay to a @f$\tau@f$ lepton pair and a light lepton pair.
     * @return @f$[\sigma_{gg\to \phi}\cdot BR(\phi\to hZ\to \tau \tau \ell \ell)]_{\text{CMS,95\%}}@f$
     * @details Taken from arXiv:1510.01181, Figure 10, left @cite Khachatryan:2015tha.
     */
    double ip_ex_ggF_A_hZ_tautaull_CMS(double mass);

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

 ///////////////////////////////////////////////////////////////////////////////

private:

    const PVfunctions PV;
    //const THDM myTHDM;

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

};

#endif	/* THDMCACHE_H */
