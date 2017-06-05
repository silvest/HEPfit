/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSMTWOFERMIONSLEP2_H
#define	EWSMTWOFERMIONSLEP2_H

#include <string>
#include <gslpp.h>
#include <Polylogarithms.h>
#include <PVfunctions.h>
#include "EWSMOneLoopEW.h"

/**
 * @class EWSMTwoFermionsLEP2
 * @ingroup StandardModel
 * @brief A class for the form factors @f$G_1@f$, @f$G_2@f$ and @f$G_3@f$ in the processes @f$e^+e^-\to f\bar{f}@f$ at LEP-II. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class EWSMTwoFermionsLEP2 {
public:

    /**
     * @brief Constructor. 
     * @param[in] cache_i reference to an EWSMcommon object
     * @param[in] bKeepNonUnitary_i true if keeping non-unitary terms
     */
    EWSMTwoFermionsLEP2(const EWSMcache& cache_i,
            const bool bKeepNonUnitary_i = false);

    ////////////////////////////////////////////////////////////////////////  

    void setBDebug(bool bDebug)
    {
        this->bDebug = bDebug;
    }

    ////////////////////////////////////////////////////////////////////////

    double G_1(const double s, const double t, const double Mw,
            const double GammaZ, const double I3f, const double Qf,
            const double mf, const double mfp, const bool bWeak,
            const bool bWWbox, const bool bZZbox) const;
    double G_2(const double s, const double t, const double Mw,
            const double GammaZ, const double I3f, const double Qf,
            const double mf, const double mfp, const bool bWeak,
            const bool bWWbox, const bool bZZbox) const;
    double G_3(const double s, const double t, const double Mw,
            const double GammaZ, const double I3f, const double Qf,
            const double mf, const double mfp, const bool bWeak,
            const bool bWWbox, const bool bZZbox) const;

    double G_1_noBox(const double s, const double Mw, const double GammaZ,
            const double I3f, const double Qf, const double mf,
            const double mfp, const bool bWeak) const;
    double G_2_noBox(const double s, const double Mw, const double GammaZ,
            const double I3f, const double Qf, const double mf,
            const double mfp, const bool bWeak) const;
    double G_3_noBox(const double s, const double Mw, const double GammaZ,
            const double I3f, const double Qf, const double mf,
            const double mfp, const bool bWeak) const;

    double G_1_box(const double s, const double t, const double Mw,
            const double GammaZ, const double I3f, const double Qf,
            const double mf, const double mfp, const bool bWWbox = true,
            const bool bZZbox = true) const;
    double G_2_box(const double s, const double t, const double Mw,
            const double GammaZ, const double I3f, const double Qf,
            const double mf, const double mfp, const bool bWWbox = true,
            const bool bZZbox = true) const;
    double G_3_box(const double s, const double t, const double Mw,
            const double GammaZ, const double I3f, const double Qf,
            const double mf, const double mfp, const bool bWWbox = true,
            const bool bZZbox = true) const;

    ////////////////////////////////////////////////////////////////////////  

    gslpp::complex V_pol(const double s) const;

    gslpp::complex chi_Z(const double s, const double Mw, const double GammaZ) const;

    gslpp::complex G_e(const double s, const double t, const double Mw,
            const double I3f, const double Qf, const double mf,
            const double mfp, const bool bWeak, const bool bWWbox,
            const bool bZZbox) const;
    gslpp::complex G_f(const double s, const double t, const double Mw,
            const double I3f, const double Qf, const double mf,
            const double mfp, const bool bWeak, const bool bWWbox,
            const bool bZZbox) const;
    gslpp::complex G_ef(const double s, const double t, const double Mw,
            const double I3f, const double Qf, const double mf,
            const double mfp, const bool bWeak, const bool bWWbox,
            const bool bZZbox) const;

    gslpp::complex rho_ef(const double s, const double t, const double Mw,
            const double I3f, const double Qf, const double mf,
            const double mfp, const bool bWeak, const bool bWWbox,
            const bool bZZbox) const;
    gslpp::complex kappa_e(const double s, const double t, const double Mw,
            const double I3f, const double Qf, const double mf,
            const double mfp, const bool bWeak, const bool bWWbox,
            const bool bZZbox) const;
    gslpp::complex kappa_f(const double s, const double t, const double Mw,
            const double I3f, const double Qf, const double mf,
            const double mfp, const bool bWeak, const bool bWWbox,
            const bool bZZbox) const;
    gslpp::complex kappa_ef(const double s, const double t, const double Mw,
            const double I3f, const double Qf, const double mf,
            const double mfp, const bool bWeak, const bool bWWbox,
            const bool bZZbox) const;

    gslpp::complex Delta_rho_ef_TOP(const double s, const double t, const double u,
            const double Mw, const bool bWWbox) const;
    gslpp::complex Delta_kappa_e_TOP(const double s, const double t, const double u,
            const double Mw, const bool bWWbox) const;
    gslpp::complex Delta_kappa_f_TOP(const double s, const double t, const double u,
            const double Mw, const bool bWWbox) const;
    gslpp::complex Delta_kappa_ef_TOP(const double s, const double t, const double u,
            const double Mw, const bool bWWbox) const;

    gslpp::complex Delta_rho_ef_WW_hat(const double s, const double t, const double u,
            const double Mw, const double I3f) const;
    gslpp::complex Delta_kappa_e_WW_hat(const double s, const double t, const double u,
            const double Mw, const double I3f) const;
    gslpp::complex Delta_kappa_f_WW_hat(const double s, const double t, const double u,
            const double Mw, const double I3f) const;
    gslpp::complex Delta_kappa_ef_WW_hat(const double s, const double t, const double u,
            const double Mw, const double I3f) const;

    gslpp::complex Delta_rho_ef_WW_TOP_hat(const double s, const double t, const double u,
            const double Mw) const;
    gslpp::complex Delta_kappa_e_WW_TOP_hat(const double s, const double t, const double u,
            const double Mw) const;
    gslpp::complex Delta_kappa_f_WW_TOP_hat(const double s, const double t, const double u,
            const double Mw) const;
    gslpp::complex Delta_kappa_ef_WW_TOP_hat(const double s, const double t, const double u,
            const double Mw) const;

    gslpp::complex Delta_rho_ef_ZZ(const double mu, const double s, const double t,
            const double u, const double Mw, const double I3f,
            const double Qf) const;
    gslpp::complex Delta_kappa_e_ZZ(const double mu, const double s, const double t,
            const double u, const double Mw, const double I3f,
            const double Qf) const;
    gslpp::complex Delta_kappa_f_ZZ(const double mu, const double s, const double t,
            const double u, const double Mw, const double I3f,
            const double Qf) const;
    gslpp::complex Delta_kappa_ef_ZZ(const double mu, const double s, const double t,
            const double u, const double Mw, const double I3f,
            const double Qf) const;

    ////////////////////////////////////////////////////////////////////////  

    // Weak corrections
    //gslpp::complex I2e(const double s, const double Mw, const bool bWeak) const;
    //gslpp::complex I2f(const double s, const double Mw, const bool bWeak) const;
    gslpp::complex DeltaRhobar(const double mu, const double Mw) const;
    gslpp::complex DeltaRhobarZ(const double mu, const double Mw) const;
    gslpp::complex D_Z(const double mu, const double s, const double Mw) const;
    gslpp::complex Pibar_Zgamma(const double mu, const double s, const double Mw) const;
    gslpp::complex Pibar_gg_bos(const double mu, const double s, const double Mw) const;
    gslpp::complex F_za_0(const double s, const double Mw) const;
    gslpp::complex F_Wa_0(const double s, const double Mw) const;
    gslpp::complex F_Wa_t(const double s, const double Mw) const;
    gslpp::complex F_Wn_0(const double s, const double Mw) const;
    gslpp::complex F_Wn_t(const double s, const double Mw) const;
    gslpp::complex F_W_0(const double s, const double Mw) const;
    gslpp::complex F_W_t(const double s, const double Mw) const;

    // WW box
    gslpp::complex B_WW_d_0(const double mu, const double s, const double t,
            const double u, const double Mw) const;
    gslpp::complex B_WW_d(const double mu, const double s, const double t,
            const double u, const double Mw) const;
    gslpp::complex Delta_B_WW_d(const double mu, const double s, const double t,
            const double u, const double Mw) const;
    gslpp::complex B_WW_c_0(const double mu, const double s, const double t,
            const double u, const double Mw) const;

    // ZZ box
    gslpp::complex B_ZZ_0(const double mu, const double s, const double t,
            const double u) const;

    // weak corrections without non-unitary terms
    gslpp::complex D_Z_hat(const double s, const double Mw) const;
    gslpp::complex Pibar_Zgamma_hat(const double s, const double Mw) const;
    gslpp::complex Pibar_gg_bos_hat(const double s, const double Mw) const;
    gslpp::complex F_Wn_0_hat(const double s, const double Mw) const;
    gslpp::complex F_Wn_t_hat(const double s, const double Mw) const;
    gslpp::complex F_W_0_hat(const double s, const double Mw) const;
    gslpp::complex F_W_t_hat(const double s, const double Mw) const;
    gslpp::complex B_WW_d_0_hat(const double s, const double t, const double u,
            const double Mw) const;
    gslpp::complex B_WW_d_0_hat_TEST(const double s, const double t, const double u,
            const double Mw) const;
    gslpp::complex Delta_B_WW_d_hat(const double s, const double t, const double u,
            const double Mw) const;
    gslpp::complex B_WW_c_0_hat(const double s, const double t, const double u,
            const double Mw) const;

    
    
        /**
     * @param[in] l lepton in the final state
     * @param[in] mf the mass of the final-state fermion
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] cosTheta cosine of the scattering angle theta
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width
     * @param[in] bWeak flag to control weak corrections (not including box diagrams)
     * @return the differential cross section d sigma(e^+ e^- -> l lbar)/d cosTheta in GeV^{-2}
     */
    double dsigma_l(const QCD::lepton l, const double mf, 
                    const double s, const double cosTheta,
                    const double Mw, const double GammaZ, const bool bWeak) const;

    /**
     * @param[in] q quark in the final state
     * @param[in] mf the mass of the final-state fermion
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] cosTheta cosine of the scattering angle theta
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width
     * @param[in] bWeak flag to control weak corrections (not including box diagrams)
     * @return the differential cross section d sigma(e^+ e^- -> q qbar)/d cosTheta in GeV^{-2}
     */
    double dsigma_q(const QCD::quark q, const double mf, 
                    const double s, const double cosTheta,
                    const double Mw, const double GammaZ, const bool bWeak) const;

    /**
     * @param[in] l lepton in the final state
     * @param[in] mf the mass of the final-state fermion
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] cosTheta cosine of the scattering angle theta
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width
     * @return the box contribution to the differential cross section d sigma(e^+ e^- -> l lbar)/d cosTheta in GeV^{-2}
     */
    double dsigma_l_box(const QCD::lepton l, const double mf, 
                        const double s, const double cosTheta,
                        const double Mw, const double GammaZ) const;
      
    /**
     * @param[in] q quark in the final state
     * @param[in] mf the mass of the final-state fermion
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] cosTheta cosine of the scattering angle theta
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width
     * @return the box contribution to the differential cross section d sigma(e^+ e^- -> q qbar)/d cosTheta in GeV^{-2}
     */
    double dsigma_q_box(const QCD::quark q, const double mf, 
                        const double s, const double cosTheta,
                        const double Mw, const double GammaZ) const;
    
    /**
     * @param[in] l lepton in the final state
     * @param[in] mf the mass of the final-state fermion
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width
     * @param[in] bWeak flag to control weak corrections (not including box diagrams)
     * @return the total cross section for e^+ e^- -> l lbar in GeV^{-2}
     */
    double sigma_l(const QCD::lepton l, const double mf, const double s, 
                   const double Mw, const double GammaZ, const bool bWeak) const;
    
    /**
     * @param[in] q quark in the final state
     * @param[in] mf the mass of the final-state fermion
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width
     * @param[in] bWeak flag to control weak corrections (not including box diagrams)
     * @return the total cross section for e^+ e^- -> q qbar in GeV^{-2}
     */
    double sigma_q(const QCD::quark q, const double mf, const double s, 
                   const double Mw, const double GammaZ, const bool bWeak) const;

    /**
     * @param[in] l lepton in the final state
     * @param[in] mf the mass of the final-state fermion
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width
     * @param[in] bWeak flag to control weak corrections (not including box diagrams)
     * @return the forward-backward asymmetry for e^+ e^- -> l lbar
     */
    double AFB_l(const QCD::lepton l, const double mf, const double s, 
                 const double Mw, const double GammaZ, const bool bWeak) const;
    
    /**
     * @param[in] q quark in the final state
     * @param[in] mf the mass of the final-state fermion
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width
     * @param[in] bWeak flag to control weak corrections (not including box diagrams)
     * @return the forward-backward asymmetry for e^+ e^- -> q qbar
     */
    double AFB_q(const QCD::quark q, const double mf, const double s, 
                 const double Mw, const double GammaZ, const bool bWeak) const;

    /**
     * @param[in] s the invariant mass squared of the initial-state e^+ e^- pair
     * @return the final-state QCD corrections to cross sections 
     */
    double QCD_FSR_forSigma(const double s) const;
    
    /**
     * @param[in] q quark in the final state
     * @param[in] mf the mass of the final-state fermion
     * @param[in] s the invariant mass squared of the initial-state e^+ e^- pair
     * @return the final-state QCD corrections to forward-backward asymmetries
     */
    double QCD_FSR_forAFB(const QCD::quark q, const double mf, 
                          const double s) const;

    /**
     * @param[in] s the invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Qf the electromagnetic charge of the final-state fermion
     * @return the final-state QED corrections to cross sections 
     */
    double QED_FSR_forSigma(const double s, const double Qf) const;
    
    /**
     * @param[in] x s'=(1-x)s
     * @param[in] s the invariant mass squared of the initial-state e^+ e^- pair
     * @return the additive radiator of initial-state radiations in cross sections
     */
    double H_ISR(const double x, const double s) const;

    /**
     * @param[in] x s'=(1-x)s
     * @param[in] s the invariant mass squared of the initial-state e^+ e^- pair
     * @return the additive radiator of initial-state radiations in forward-backward asysmmetries
     */
    double H_ISR_FB(const double x, const double s) const;

    /**
     * @param[in] l lepton in the final state
     * @param[in] mf the mass of the final-state fermion
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width
     * @param[in] bWeak flag to control weak corrections (not including box diagrams)
     * @return the form factor beta_f^2*G_3(s) for e^+ e^- -> l lbar
     */
    double G_3prime_l(const QCD::lepton l, const double mf, 
                      const double s, const double Mw, const double GammaZ, 
                      const bool bWeak) const;
    
    /**
     * @param[in] q quark in the final state
     * @param[in] mf the mass of the final-state fermion
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width
     * @param[in] bWeak flag to control weak corrections (not including box diagrams)
     * @return the form factor beta_f^2*G_3(s) for e^+ e^- -> q qbar
     */
    double G_3prime_q(const QCD::quark q, const double mf, 
                      const double s, const double Mw, const double GammaZ, 
                      const bool bWeak) const;   

    
    

    ////////////////////////////////////////////////////////////////////////  
private:
    
    double alpha_at_s(const double s) const;

    double dsigma(const double s, const double cosTheta, 
                  const double Mw, const double GammaZ, 
                  const double I3f, const double Qf, const double mf,
                  const double mfp, const double Ncf, const bool bWeak) const;

    double dsigma_box(const double s, const double cosTheta, 
                      const double Mw, const double GammaZ, 
                      const double I3f, const double Qf, const double mf,
                      const double mfp, const double Ncf) const;
    
    double sigma(const double s, const double Mw, const double GammaZ, 
                 const double I3f, const double Qf, const double mf,
                 const double mfp, const double Ncf, const bool bWeak) const;
    
    double AFB(const double s, const double Mw, const double GammaZ, 
               const double I3f, const double Qf, const double mf,
               const double mfp, const bool bWeak) const;    
    
    
    bool bDebug; // for debug
    bool bKeepNonUnitary; // true if keeping non-unitary terms

    const EWSMcache& cache; ///< A reference to an object of type EWSMcache.
    const EWSMOneLoopEW myOneLoopEW; ///< An object of type EWSMOneLoopEW.

};

#endif	/* EWSMTWOFERMIONSLEP2_H */

