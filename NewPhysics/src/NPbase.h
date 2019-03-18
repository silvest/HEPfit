/*
 * Copyright (C) 2013 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPBASE_H
#define	NPBASE_H

#include "StandardModel.h"

/**
 * @class NPbase
 * @brief The auxiliary base model class for other model classes.
 * @ingroup NewPhysics
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 *  
 * @details This is an auxiliary Model class containing the basic structure 
 * to compute new physics (NP) contributions to the electroweak precision
 * observables. The NP contributions are described by the following quantities:
 *
 * @li @f$\Delta G@f$&nbsp;&nbsp; (with DeltaGF()),
 * @li @f$S@f$, @f$T@f$ and @f$U@f$&nbsp;&nbsp;
 * (with obliqueS(), obliqueT() and obliqueU()),
 *
 * Using these quantities, the mass and width of the @f$W@f$ boson and the NP
 * contributions to the effective vector and axial-vector couplings of the
 * @f$Z@f$-boson to leptons and quarks are computed:
 *
 * @li @f$M_W@f$, @f$c_W^2@f$ and @f$s_W^2@f$&nbsp;&nbsp; (with Mw(), cW2() and sW2()),
 * @li @f$\Gamma_W@f$&nbsp;&nbsp; (with GammaW()),
 * @li @f$\delta g_V^f@f$&nbsp;&nbsp;(with deltaGVl() and deltaGVq()),
 * @li @f$\delta g_A^f@f$&nbsp;&nbsp;(with deltaGAl() and deltaGAq()). 
 *
 * In the model classes that are inherited from the current class (see the
 * inheritance diagram above), some of these methods are reimplemented to
 * account for the details of more specific scenarios.
 *
 * 
 * @anchor NPbaseInitialization
 * <h3>Initialization</h3>
 *
 * This class is intended to be used with an inherited model class.
 *
 *
 * @anchor NPbaseParameters
 * <h3>%Model parameters</h3>
 *
 * There is no model parameter in the current class.
 *
 * 
 * @anchor NPbaseFlags
 * <h3>%Model flags</h3>
 *
 * There is no model flag in the current class.
 *
 *
 * @anchor NPbaseFunctions
 * <h3>Important member functions</h3>
 *
 * The functions are explained above. 
 *
 */
class NPbase : public StandardModel {
public:

    /**
     * @brief The default constructor.
     */
    NPbase();

    /**
     * @brief The update method for %NPbase.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    /**
     * @brief The postupdate method for %NPbase.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PostUpdate();
    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief A method to return a StandardModel object from %NPbase.
     * @return an StandardModel object
     */
    virtual StandardModel getTrueSM() const
    {
        return trueSM;
    }

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief New physics contribution to the Fermi constant.
     * @details The new physics contribution @f$\Delta G@f$ is defined as
     * @f[
     * G_\mu = G_{\mu,\mathrm{SM}}(1+\Delta G)\,,
     * @f]
     * where @f$G_\mu@f$ is the experimental value measured through muon decays, 
     * and @f$G_{\mu,\mathrm{SM}}@f$ is the Fermi constant in the SM.
     * @return @f$\Delta G@f$
     */
    virtual double DeltaGF() const
    {
        return 0.;
    }

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The oblique parameter \f$S\f$.
     * @return the value of @f$S@f$
     */
    virtual double obliqueS() const
    {
        return 0.;
    }

    /**
     * @brief The oblique parameter \f$T\f$.
     * @return the value of @f$T@f$
     */
    virtual double obliqueT() const
    {
        return 0.;
    }

    /**
     * @brief The oblique parameter \f$U\f$.
     * @return the value of @f$U@f$
     */
    virtual double obliqueU() const
    {
        return 0.;
    }
    
    /**
     * @brief The oblique parameter \f$W\f$.
     * @return the value of @f$W@f$
     */
    virtual double obliqueW() const
    {
        return 0.;
    }
    
    /**
     * @brief The oblique parameter \f$Y\f$.
     * @return the value of @f$Y@f$
     */
    virtual double obliqueY() const
    {
        return 0.;
    }

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The mass of the @f$W@f$ boson, @f$M_W@f$.
     * @details
     * The @f$W@f$-boson mass receives the new physics
     * contribution via the oblique parameters @f$S@f$, @f$T@f$ and @f$U@f$ and
     * the shift in the Fermi constant, @f$\Delta G@f$:
     * @f[
     * M_W = M_{W,\mathrm{SM}}
     * \left[
     * 1 - \frac{\alpha(M_Z^2)}{4(c_W^2-s_W^2)}
     * \left( S - 2c_W^2\,T - \frac{c_W^2-s_W^2}{2s_W^2}\,U \right)
     * - \frac{s_W^2}{2(c_W^2-s_W^2)}\,\Delta G
     * \right].
     * @f]
     *
     * See @cite Ciuchini:2013pca and references therein.
     * @return @f$M_W@f$ in GeV
     */
    virtual double Mw() const;
    
    /**
     * @brief The new physics contribution to the decay width of the @f$W@f$ boson into a given fermion pair, @f$\delta \Gamma_Z^{f}@f$.
     * @param[in] fi a lepton or quark
     * @param[in] fj a lepton or quark
     * @return @f$\delta \Gamma_W^{ff}@f$ in GeV
     */
    virtual double deltaGamma_Wff(const Particle fi, const Particle fj) const
    {
        return 0.0;
    };
    
    /**
     * @brief A partial decay width of the @f$W@f$ boson decay into a SM fermion pair.
     * @details
     * The partial @f$W@f$-boson widths receives the new physics
     * contribution via the oblique parameters @f$S@f$, @f$T@f$ and @f$U@f$ and
     * the shift in the Fermi constant, @f$\Delta G@f$:
     * @f[
     * \Gamma_W^{ij} = \Gamma_{W,\mathrm{SM}}
     * \left[ 1
     * - \frac{3\alpha(M_Z^2)}{4(c_W^2-s_W^2)}
     *  \left( S - 2c_W^2\,T - \frac{c_W^2-s_W^2}{2s_W^2}\,U \right)
     * - \frac{1+c_W^2}{2(c_W^2-s_W^2)}\, \Delta G
     * \right].
     * @f]
     * @param[in] fi a lepton or quark
     * @param[in] fj a lepton or quark
     * @return @f$\Gamma^W_{ij}@f$
     *
     * @attention Fermion masses are neglected.
     */
    virtual double GammaW(const Particle fi, const Particle fj) const;
    
    /**
     * @brief The new physics contribution to the total decay width of the @f$W@f$ boson, @f$\delta \Gamma_W@f$.
     * @return @f$\delta \Gamma_W@f$ in GeV
     */
    virtual double deltaGamma_W() const
    {
        return 0.0;
    };

    /**
     * @brief The total width of the @f$W@f$ boson, @f$\Gamma_W@f$.
     * @details
     * The @f$W@f$-boson width receives the new physics
     * contribution via the oblique parameters @f$S@f$, @f$T@f$ and @f$U@f$ and
     * the shift in the Fermi constant, @f$\Delta G@f$:
     * @f[
     * \Gamma_W = \Gamma_{W,\mathrm{SM}}
     * \left[ 1
     * - \frac{3\alpha(M_Z^2)}{4(c_W^2-s_W^2)}
     *  \left( S - 2c_W^2\,T - \frac{c_W^2-s_W^2}{2s_W^2}\,U \right)
     * - \frac{1+c_W^2}{2(c_W^2-s_W^2)}\, \Delta G
     * \right].
     * @f]
     *
     * See @cite Ciuchini:2013pca and references therein.
     * @return @f$\Gamma_W@f$ in GeV
     */
    virtual double GammaW() const;

    /**
     * @brief New physics contribution to the neutral-current vector coupling @f$g_V^f@f$.
     * @details
     * The neutral-current vector coupling @f$g_V^f@f$ receives the new physics
     * contribution via the oblique parameters @f$S@f$ and @f$T@f$ and the shift
     * in the Fermi constant, @f$\Delta G@f$:
     * @f[
     * \delta g_V^f =
     * \frac{g_{V,\mathrm{SM}}^f}{2}
     * \left[ \alpha(M_Z^2)\, T - \Delta G \right]
     * +
     * \frac{\big( g_{V,\mathrm{SM}}^f - g_{A,\mathrm{SM}}^f \big)
     * \left[
     * \alpha(M_Z^2)\left( S - 4\,c_W^2s_W^2\, T \right)
     * + 4\,c_W^2s_W^2\, \Delta G
     * \right]}{4s_W^2\,(c_W^2-s_W^2)}\,.
     * @f]
     *
     * See @cite Ciuchini:2013pca and references therein.
     * @param[in] f a lepton or quark
     * @return @f$\delta g_V^f@f$
     */
    virtual double deltaGV_f(const Particle f) const;

    /**
     * @brief The total (SM+NP) contribution to the neutral-current vector coupling @f$g_V^f@f$.
     * @param[in] f a lepton or quark
     * @return @f$g_V^f@f$, including SM plus NP contributions
     */
    virtual gslpp::complex gV_f(const Particle f) const;

    /**
     * @brief New physics contribution to the neutral-current axial-vector coupling @f$g_A^f@f$.
     * @details
     * The neutral-current axial-vector coupling @f$g_A^f@f$ receives the new
     * physics contribution via the oblique parameter @f$T@f$ and the shift in
     * the Fermi constant, @f$\Delta G@f$:
     * @f[
     * \delta g_A^f
     * = \frac{g_{A,\mathrm{SM}}^f}{2} \left[ \alpha(M_Z^2)\, T - \Delta G \right].
     * @f]
     *
     * See @cite Ciuchini:2013pca and references therein.
     * @param[in] f a lepton or quark
     * @return @f$\delta g_A^f@f$
     */
    virtual double deltaGA_f(const Particle f) const;

    /**
     * @brief The total (SM+NP) contribution to the neutral-current axial-vector coupling @f$g_A^f@f$.
     * @param[in] f a lepton or quark
     * @return @f$g_A^f@f$, including SM plus NP contributions
     */
    virtual gslpp::complex gA_f(const Particle f) const;

    /**
     * @brief The effective neutral-current coupling @f$\rho_Z^f@f$ including SM plus NP contributions.
     * @param[in] f a lepton or quark
     * @return @f$\rho_Z^f@f$, including SM plus NP contributions
     */
    virtual gslpp::complex rhoZ_f(const Particle f) const;

    /**
     * @brief The effective neutral-current coupling @f$\kappa_Z^f@f$ including SM plus NP contributions.
     * @param[in] f a lepton or quark
     * @return @f$\kappa_Z^f@f$, including SM plus NP contributions
     */
    virtual gslpp::complex kappaZ_f(const Particle f) const;
    
    /**
     * @brief The new physics contribution to the decay width of the @f$Z@f$ boson into a given fermion pair, @f$\delta \Gamma_Z^{f}@f$.
     * @param[in] f a lepton or quark
     * @return @f$\delta \Gamma_Z^{f}@f$ in GeV
     */
    virtual double deltaGamma_Zf(const Particle f) const;
    
    /**
     * @brief The decay width of the @f$Z@f$ boson into a given fermion pair, @f$\Gamma_Z^{f}@f$.
     * @param[in] f a lepton or quark
     * @return @f$\Gamma_Z^{f}@f$ in GeV, including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    virtual double Gamma_Zf(const Particle f) const;

    /**
     * @brief The new physics contribution to the total decay width of the @f$Z@f$ boson, @f$\delta \Gamma_Z@f$.
     * @return @f$\delta \Gamma_Z@f$ in GeV
     */
    virtual double deltaGamma_Z() const;

    /**
     * @brief The total decay width of the @f$Z@f$ boson, @f$\Gamma_Z@f$.
     * @return @f$\Gamma_Z@f$ in GeV, including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    virtual double Gamma_Z() const;
    
    /**
     * @brief The Branching ratio of the @f$Z@f$ boson into a given fermion pair, @f$BR_Z^{f}@f$.
     * @param[in] f a lepton or quark
     * @return @f$BR_Z^{f}@f$ including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    virtual double BR_Zf(const Particle f) const;

    /**
     * @brief The new physics contribution to the cross section for the process @f$e^+ e^-\to Z\to \mathrm{hadrons}@f$
     * at the @f$Z@f$ pole, @f$\delta \sigma_h^0@f$.
     * @return @f$\delta \sigma_h^0@f$ in GeV@f$^{-2}@f$
     */
    virtual double deltaSigmaHadron() const;

    /**
     * @brief The cross section for the process @f$e^+ e^-\to Z\to \mathrm{hadrons}@f$
     * at the @f$Z@f$ pole, @f$\sigma_h^0@f$.
     * @return @f$\sigma_h^0@f$ in GeV@f$^{-2}@f$, including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    virtual double sigma0_had() const;

    /**
     * @brief The new physics contribution to the effective electron/leptonic weak angle @f$\delta \sin^2\theta_{\rm eff}^{\rm lept}@f$
     * at the @f$Z@f$ pole.
     * @return @f$\delta \sin^2\theta_{\rm eff}^{\rm lept}@f$
     */
    virtual double deltaSin2thetaEff_e() const;
    
    /**
     * @brief The new physics contribution to the effective muonic weak angle @f$\delta \sin^2\theta_{\rm eff}^{\mu\mu}@f$
     * at the @f$Z@f$ pole.
     * @return @f$\delta \sin^2\theta_{\rm eff}^{\mu\mu}@f$
     */
    virtual double deltaSin2thetaEff_mu() const;

    /**
     * @brief @copybrief sin2thetaEff::computeThValue()
     * @param[in] f a lepton or quark
     * @return @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$, including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    virtual double sin2thetaEff(const Particle f) const;

    /**
     * @brief The new physics contribution to the left-right asymmetry in @f$e^+e^-\to Z\to f \bar{f}@f$ at the
     * @f$Z@f$-pole, @f$\delta \mathcal{A}_f@f$. 
     * @param[in] f a lepton or quark
     * @return @f$\delta \mathcal{A}_f@f$
     */
    virtual double deltaA_f(const Particle f) const;

    /**
     * @brief The left-right asymmetry in @f$e^+e^-\to Z\to f \bar{f}@f$ at the
     * @f$Z@f$-pole, @f$\mathcal{A}_f@f$.
     * @param[in] f a lepton or quark
     * @return @f$\mathcal{A}_f@f$, including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    virtual double A_f(const Particle f) const;

    /**
     * @brief The new physics contribution to the forward-backward asymmetry in @f$e^+e^-\to Z\to f \bar{f}@f$ at the
     * @f$Z@f$-pole, @f$\delta A^f_{FB}@f$. 
     * @param[in] f a lepton or quark
     * @return @f$\delta A^f_{FB}@f$
     */
    virtual double deltaAFB(const Particle f) const;

    /**
     * @brief The forward-backward asymmetry in @f$e^+e^-\to Z\to f \bar{f}@f$ at the
     * @f$Z@f$-pole, @f$A^f_{FB}@f$.
     * @param[in] f a lepton or quark
     * @return @f$A^f_{FB}@f$, including SM plus NP contributions
     */
    virtual double AFB(const Particle f) const;

    /**
     * @brief The new physics contribution to the ratio @f$R_\ell^0=\Gamma_{\mathrm{had}}/\Gamma_\ell@f$
     * or @f$R_q^0=\Gamma_q/\Gamma_{\mathrm{had}}@f$, for leptons or quarks, respectively.
     * @param f a lepton or quark
     * @return @f$\delta R_f^0@f$
     */
    virtual double deltaR0_f(const Particle f) const;

    /**
     * @brief The ratio @f$R_\ell^0=\Gamma_{\mathrm{had}}/\Gamma_\ell@f$
     * or @f$R_q^0=\Gamma_q/\Gamma_{\mathrm{had}}@f$, for leptons or quarks, respectively. 
     * @param[in] f a lepton or quark
     * @return @f$R_f^0@f$, including SM plus NP contributions
     */
    virtual double R0_f(const Particle f) const;
    
    /**
     * @brief The new physics contribution to the ratio of invisible and leptonic (electron) decay widths of the @f$Z@f$ boson, @f$\delta R_{inv}@f$.
     * @return @f$\delta R_{inv}@f$
     */
    virtual double deltaR_inv() const;

    /**
     * @brief The ratio of the invisible and leptonic (electron) decay widths of the @f$Z@f$ boson, @f$R_{inv}@f$.
     * @return @f$R_{inv}@f$, including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    virtual double R_inv() const;
    
    /**
     * @brief The new physics contribution to the number of neutrinos dervied from the @f$Z@f$ pole measurements.
     * @return @f$\delta N_{\nu}@f$
     */
    virtual double deltaN_nu() const;

    /**
     * @brief The number of neutrinos dervied from the @f$Z@f$ pole measurements, @f$N_{\nu}@f$.
     * @return @f$N_{\nu}@f$, including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    virtual double N_nu() const;

    ////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief New physics contribution to the charged current coupling @f$W_\mu \bar{f_L}\gamma^mu f_L@f$.
     * @param[in] pbar a lepton or quark
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{Wff}^{L}@f$
     */
    // no generation mixing
    virtual gslpp::complex deltaGL_Wff(const Particle pbar, const Particle p) const
    {
        return 0.0;
    };
    /**
     * @brief New physics contribution to the charged current coupling @f$W_\mu \bar{f_R}\gamma^mu f_R@f$.
     * @param[in] pbar a lepton or quark
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{Wff}^{R}@f$
     */
    // no generation mixing
    virtual gslpp::complex deltaGR_Wff(const Particle pbar, const Particle p) const
    {
        return 0.0;
    };

    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H G_{\mu\nu}^AG^{A \mu\nu}@f$.
     * @return @f$\delta g_{HGG}@f$
     */
    virtual double deltaG_hgg() const
    {
        return 0.0;
    };
    /**
     * @brief The full new physics contribution to the coupling of the effective interaction @f$H G_{\mu\nu}^AG^{A \mu\nu}@f$,
     * including new local terms and modifications on the SM-loops. Normalized to the SM value.
     * @return @f$\delta g_{HGG}/g_{HGG}^SM}@f$
     */
    virtual double deltaG_hggRatio() const
    {
        return 0.0;
    };
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H W_{\mu\nu}^\dagger W^{\mu\nu}@f$.
     * @return @f$\delta g_{HWW}^{(1)}@f$
     */
    virtual double deltaG1_hWW() const
    {
        return 0.0;
    };
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H W_{\nu}^\dagger \partial^\mu W^{\mu\nu}@f$.
     * @return @f$\delta g_{HWW}^{(2)}@f$
     */
    virtual double deltaG2_hWW() const
    {
        return 0.0;
    };
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H W_{\mu}^\dagger W^{\mu}@f$.
     * @return @f$\delta g_{HWW}^{(3)}@f$
     */
    virtual double deltaG3_hWW() const
    {
        return 0.0;
    };
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\mu\nu} Z^{\mu\nu}@f$.
     * @return @f$\delta g_{HZZ}^{(1)}@f$
     */
    virtual double deltaG1_hZZ() const
    {
        return 0.0;
    };
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\nu} \partial^\mu Z^{\mu\nu}@f$.
     * @return @f$\delta g_{HZZ}^{(2)}@f$
     */
    virtual double deltaG2_hZZ() const
    {
        return 0.0;
    };
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\mu} Z^{\mu}@f$.
     * @return @f$\delta g_{HZZ}^{(3)}@f$
     */
    virtual double deltaG3_hZZ() const
    {
        return 0.0;
    };
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\mu\nu} F^{\mu\nu}@f$.
     * @return @f$\delta g_{HZA}^{(1)}@f$
     */
    virtual double deltaG1_hZA() const
    {
        return 0.0;
    };
    /**
     * @brief The full new physics contribution to the coupling of the effective interaction @f$H Z_{\mu\nu} F^{A \mu\nu}@f$,
     * including new local terms and modifications on the SM-loops. Normalized to the SM value.
     * @return @f$\delta g_{HZA}^{(1)}/g_{HZA}^{(1),SM}@f$
     */
    virtual double deltaG1_hZARatio() const
    {
        return 0.0;
    };
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\nu} \partial^\mu F^{\mu\nu}@f$.
     * @return @f$\delta g_{HZA}^{(2)}@f$
     */
    virtual double deltaG2_hZA() const
    {
        return 0.0;
    };
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H F_{\mu\nu} F^{\mu\nu}@f$.
     * @return @f$\delta g_{HAA}@f$
     */
    virtual double deltaG_hAA() const
    {
        return 0.0;
    };
    /**
     * @brief The full new physics contribution to the coupling of the effective interaction @f$H F_{\mu\nu} F^{\mu\nu}@f$,
     * including new local terms and modifications on the SM-loops. Normalized to the SM value.
     * @return @f$\delta g_{HAA}/g_{HAA}^SM}@f$
     */
    virtual double deltaG_hAARatio() const
    {
        return 0.0;
    };
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H f\bar{f}@f$.
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{Hff}@f$
     */
    // no generation mixing
    virtual gslpp::complex deltaG_hff(const Particle p) const
    {
        return 0.0;
    };
    /**
     * @brief The new physics contribution to the Higgs self-coupling @f$ H H H@f$. Normalized to the SM value.
     * @return @f$\delta g_{HHH}/g_{HHH}^SM}@f$
     */
    virtual double deltaG_hhhRatio() const
    {
        return 0.0;
    };

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The ratio @f$\mu_{ggH}@f$ between the gluon-gluon fusion Higgs
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH}@f$
     */
    virtual double muggH(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{ggHH}@f$ between the gluon-gluon fusion di-Higgs
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggHH}@f$
     */
    virtual double muggHH(const double sqrt_s) const
    {
        return 1.0;
    }

    /**
     * @brief The ratio @f$\mu_{VBF}@f$ between the vector-boson fusion Higgs
     * production cross-section in the current model and in the Standard Model. 
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF}@f$
     */
    virtual double muVBF(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{VBF+\gamma}@f$ between the vector-boson fusion Higgs
     * production cross-section in association with a hard photon in the current model
     * and in the Standard Model. 
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF+\gamma}@f$
     */
    virtual double muVBFgamma(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{eeWBF}@f$ between the 
     * @f$ e^{+}e^{-}\to \nu\bar{\nu} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeWBF}@f$
     */
    virtual double mueeWBF(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{eeWBF}@f$ between the 
     * @f$ e^{+}e^{-}\to \nu\bar{\nu} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively
     * @return @f$\mu_{eeWBF}@f$
     */
    virtual double mueeWBFPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$ between the 
     * @f$ e^+e^- \to H\nu\bar{\nu} @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$
     */
    virtual double mueeHvv(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$ between the 
     * @f$ e^+e^- \to H\nu\bar{\nu} @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$
     */
    virtual double mueeHvvPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{eeZBF}@f$ between the 
     * @f$ e^{+}e^{-}\to e^{+}e^{-} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeZBF}@f$
     */
    virtual double mueeZBF(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{eeZBF}@f$ between the 
     * @f$ e^{+}e^{-}\to e^{+}e^{-} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeZBF}@f$
     */
    virtual double mueeZBFPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{epWBF}@f$ between the 
     * @f$ e^{-} p\to \nu j H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{epWBF}@f$
     */
    virtual double muepWBF(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{epZBF}@f$ between the 
     * @f$ e^{-} p\to e^{-} j H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{epZBF}@f$
     */
    virtual double muepZBF(const double sqrt_s) const
    {
        return 1.0;
    }

    /**
     * @brief The ratio @f$\mu_{WH}@f$ between the W-Higgs associated production
     * cross-section in the current model and in the Standard Model. 
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH}@f$
     */
    virtual double muWH(const double sqrt_s) const
    {
        return 1.0;
    }

    /**
     * @brief The ratio @f$\mu_{ZH}@f$ between the Z-Higgs associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH}@f$
     */
    virtual double muZH(const double sqrt_s) const
    {
        return 1.0;
    }

    /**
     * @brief The ratio @f$\mu_{eeZH}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeZH}@f$
     */
    virtual double mueeZH(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{eeZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH, Z \to e^+ e^-, \mu^+ \mu^- @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$
     */
    virtual double mueeZllH(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{eeZH, Z \to q \bar{q}}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH, Z \to q \bar{q} @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeZH, Z \to q \bar{q}}@f$
     */
    virtual double mueeZqqH(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{eeZH}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeZH}@f$
     */
    virtual double mueeZHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{eeZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH, Z \to e^+ e^-, \mu^+ \mu^- @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$
     */
    virtual double mueeZllHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{eeZH, Z \to q \bar{q}}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH, Z \to q \bar{q} @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeZH, Z \to q \bar{q}}@f$
     */
    virtual double mueeZqqHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
    {
        return 1.0;
    }

    /**
     * @brief The ratio @f$\mu_{VH}@f$ between the WH+ZH associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH}@f$
     */
    virtual double muVH(const double sqrt_s) const
    {
        return 1.0;
    }

    /**
     * @brief The ratio @f$\mu_{VBF+VH}@f$ between the sum of VBF and WH+ZH associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF+VH}@f$
     */
    virtual double muVBFpVH(const double sqrt_s) const
    {
        return 1.0;
    }

    /**
     * @brief The ratio @f$\mu_{ttH}@f$ between the t-tbar-Higgs associated 
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH}@f$
     */
    virtual double muttH(const double sqrt_s) const
    {
        return 1.0;
    }

    /**
     * @brief The ratio @f$\mu_{ggH+ttH}@f$ between the sum of gluon-gluon fusion
     * and t-tbar-Higgs associated 
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH+ttH}@f$
     */
    virtual double muggHpttH(const double sqrt_s) const
    {
        return 1.0;
    }
        
    /**
     * @brief The ratio @f$\mu_{eettH}@f$ between the 
     * @f$ e^{+}e^{-}\to t\bar{t} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eettH}@f$
     */
    virtual double mueettH(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{eettH}@f$ between the 
     * @f$ e^{+}e^{-}\to t\bar{t} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively
     * @return @f$\mu_{eettH}@f$
     */
    virtual double mueettHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{\mu\mu H}@f$ between the @f$\sigma(\mu \mu \to H)}@f$
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{\mu\mu H}@f$
     */
    virtual double mummH(const double sqrt_s) const
    {
        return 1.0;
    }

    /**
     * @brief The ratio of the Br@f$(H\to gg)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to gg)@f$/Br@f$(H\to gg)_{\mathrm{SM}}@f$
     */
    virtual double BrHggRatio() const
    {
        return 1.0;
    }

    /**
     * @brief The ratio of the Br@f$(H\to WW)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to WW)@f$/Br@f$(H\to WW)_{\mathrm{SM}}@f$
     */
    virtual double BrHWWRatio() const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio of the Br@f$(H\to W l\nu)@f$ (@f$l=e,\mu @f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Wl\nu)@f$/Br@f$(H\to Wl\nu)_{\mathrm{SM}}@f$
     */
    virtual double BrHWlvRatio() const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio of the Br@f$(H\to WW^*\to l\nu l\nu)@f$ (@f$l=e,\mu @f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to WW^*\to l\nu l\nu)@f$/Br@f$(H\to WW^*\to l\nu l\nu)_{\mathrm{SM}}@f$
     */
    virtual double BrHWW2l2vRatio() const
    {
        return 1.0;
    }

    /**
     * @brief The ratio of the Br@f$(H\to ZZ)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ)@f$/Br@f$(H\to ZZ)_{\mathrm{SM}}@f$
     */
    virtual double BrHZZRatio() const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio of the Br@f$(H\to Zll)@f$ (@f$l=e,\mu @f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Zll)@f$/Br@f$(H\to Zll)_{\mathrm{SM}}@f$
     */
    virtual double BrHZllRatio() const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio of the Br@f$(H\to ZZ* \to 4l)@f$ (@f$l=e,\mu @f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ* \to 4l)@f$/Br@f$(H\to ZZ* \to 4l)_{\mathrm{SM}}@f$
     */
    virtual double BrHZZ4lRatio() const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio of the Br@f$(H\to ZZ* \to 4e)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ* \to 4e)@f$/Br@f$(H\to ZZ* \to 4e)_{\mathrm{SM}}@f$
     */
    virtual double BrHZZ4eRatio() const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio of the Br@f$(H\to ZZ* \to 2e 2\mu)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ* \to 2e 2\mu)@f$/Br@f$(H\to ZZ* \to 2e 2\mu)_{\mathrm{SM}}@f$
     */
    virtual double BrHZZ2e2muRatio() const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio of the Br@f$(H\to ZZ* \to 4\mu)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ* \to 4\mu)@f$/Br@f$(H\to ZZ* \to 4\mu)_{\mathrm{SM}}@f$
     */
    virtual double BrHZZ4muRatio() const
    {
        return 1.0;
    }

    /**
     * @brief The ratio of the Br@f$(H\to Z\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z\gamma)@f$/Br@f$(H\to Z\gamma)_{\mathrm{SM}}@f$
     */
    virtual double BrHZgaRatio() const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio of the Br@f$(H\to Z\gamma\to ll\gamma)@f$ (@f$l=e,\mu @f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z\gamma\to ll\gamma)@f$/Br@f$(H\to Z\gamma\to ll\gamma)_{\mathrm{SM}}@f$
     */
    virtual double BrHZgallRatio() const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio of the Br@f$(H\to Z\gamma\to ee\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z\gamma\to ee\gamma)@f$/Br@f$(H\to Z\gamma\to ee\gamma)_{\mathrm{SM}}@f$
     */
    virtual double BrHZgaeeRatio() const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio of the Br@f$(H\to Z\gamma\to \mu\mu\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z\gamma\to \mu\mu\gamma)@f$/Br@f$(H\to Z\gamma\to \mu\mu\gamma)_{\mathrm{SM}}@f$
     */
    virtual double BrHZgamumuRatio() const
    {
        return 1.0;
    }

    /**
     * @brief The ratio of the Br@f$(H\to \gamma\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to \gamma\gamma)@f$/Br@f$(H\to \gamma\gamma)_{\mathrm{SM}}@f$
     */
    virtual double BrHgagaRatio() const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio of the Br@f$(H\to \mu^+\mu^-)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to \mu^+\mu^-)@f$/Br@f$(H\to \mu^+\mu^-)_{\mathrm{SM}}@f$
     */
    virtual double BrHmumuRatio() const
    {
        return 1.0;
    }

    /**
     * @brief The ratio of the Br@f$(H\to \tau^+\tau^-)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to \tau^+\tau^-)@f$/Br@f$(H\to \tau^+\tau^-)_{\mathrm{SM}}@f$
     */
    virtual double BrHtautauRatio() const
    {
        return 1.0;
    }

    /**
     * @brief The ratio of the Br@f$(H\to c\bar{c})@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to c\bar{c})@f$/Br@f$(H\to c\bar{c})_{\mathrm{SM}}@f$
     */
    virtual double BrHccRatio() const
    {
        return 1.0;
    }

    /**
     * @brief The ratio of the Br@f$(H\to b\bar{b})@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to b\bar{b})@f$/Br@f$(H\to b\bar{b})_{\mathrm{SM}}@f$
     */
    virtual double BrHbbRatio() const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\sigma(ttH)/\sigma(ttZ)@f$ 
     * in the @f$H,Z\to b\bar{b}@f$ channel in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\sigma(ttH)/\sigma(ttZ)@f$ normalized to the SM
     */
    virtual double muttHZbbboost(const double sqrt_s) const
    {
        return 1.0;
    }
        
    virtual double muggHgaga(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{ggH,\gamma\gamma}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model. Includes interference effects
     * with the background, following arXiv:1704.08259
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,\gamma\gamma}@f$
     */
    virtual double muggHgagaInt(const double sqrt_s) const
    {
        return 1.0;        
    };
    
    virtual double muVBFHgaga(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muZHgaga(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muWHgaga(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVHgaga(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muttHgaga(const double sqrt_s) const
    {
        return 1.0;
    }   
    virtual double muggHZga(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVBFHZga(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muZHZga(const double sqrt_s) const
    {
        return 1.0;
    }    
    virtual double muWHZga(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVHZga(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muttHZga(const double sqrt_s) const
    {
        return 1.0;
    }  
    virtual double muggHZZ(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVBFHZZ(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muZHZZ(const double sqrt_s) const
    {
        return 1.0;
    }    
    virtual double muWHZZ(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVHZZ(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muttHZZ(const double sqrt_s) const
    {
        return 1.0;
    }
    
    virtual double muggHZZ4l(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVBFHZZ4l(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muZHZZ4l(const double sqrt_s) const
    {
        return 1.0;
    }    
    virtual double muWHZZ4l(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVHZZ4l(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muttHZZ4l(const double sqrt_s) const
    {
        return 1.0;
    }
    
    
    virtual double muggHWW(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVBFHWW(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muZHWW(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muWHWW(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVHWW(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muttHWW(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muggHWW2l2v(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVBFHWW2l2v(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muZHWW2l2v(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muWHWW2l2v(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVHWW2l2v(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muttHWW2l2v(const double sqrt_s) const
    {
        return 1.0;
    }   
    virtual double muggHmumu(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVBFHmumu(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muZHmumu(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muWHmumu(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVHmumu(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muttHmumu(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muggHtautau(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVBFHtautau(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muZHtautau(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muWHtautau(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVHtautau(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muttHtautau(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muggHbb(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVBFHbb(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muZHbb(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muWHbb(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVHbb(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muttHbb(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muppHmumu(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muppHZga(const double sqrt_s) const
    {
        return 1.0;
    }

    /**
     * @brief The ratio of the @f$\Gamma(H)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H)@f$/@f$\Gamma(H)_{\mathrm{SM}}@f$
     */
    virtual double computeGammaTotalRatio() const
    {
        return 1.0;
    }
    
    /**
     * @brief The branching ratio of the of the Higgs into exotic particles.
     * @return Br@f$(H\to exotic)@f$
     */
    virtual double Br_H_exo() const
    {
        return 0.0;
    };
    
    /**
     * @brief The branching ratio of the of the Higgs into invisible particles.
     * @return Br@f$(H\to invisible)@f$
     */
    virtual double Br_H_inv() const
    {
        return 0.0;
    };
    
    /**
     * @brief The branching ratio of the of the Higgs into invisible particles 
     * (only invisible new particles).
     * @return Br@f$(H\to invisible,NP)@f$
     */
    virtual double Br_H_inv_NP() const
    {
        return 0.0;
    };
    
    /**
     * @brief The ratio of the Br@f$(H\to visible)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to visible)@f$/Br@f$(H\to visible)_{\mathrm{SM}}@f$
     */
    virtual double BrHvisRatio() const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio of the Br@f$(H\to invisible)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to invisible)@f$/Br@f$(H\to ZZ \to invisible)_{\mathrm{SM}}@f$
     */
    virtual double BrHtoinvRatio() const
    {
        return 1.0;
    }

    /**
     * @brief 
     * @return 
     */
    virtual double UpperLimitZgammaA13(const double sqrt_s) const
    {
        return 1.0;
    }

    /**
     * @brief 
     * @return 
     */
    virtual double UpperLimitZgammaC13(const double sqrt_s) const
    {
        return 1.0;
    }

    /**
     * @brief 
     * @return 
     */
    virtual double UpperLimitZgammaA(const double sqrt_s) const
    {
        return 1.0;
    }

    /**
     * @brief 
     * @return 
     */
    virtual double UpperLimitZgammaC(const double sqrt_s) const
    {
        return 1.0;
    }

    /**
     * @brief 
     * @return 
     */
    virtual double cgplusct() const
    {
        return 1.0;
    }

    /**
     * @brief 
     * @return 
     */
    virtual double cgaplusct() const
    {
        return 1.0;
    }

    /**
     * @brief 
     * @return 
     */
    virtual double cgminuscga() const
    {
        return 0.0;
    }

    /**
     * @brief 
     * @return 
     */
    virtual double cVpluscb() const
    {
        return 2.0;
    }

    /**
     * @brief 
     * @return 
     */
    virtual double cVplusctau() const
    {
        return 2.0;
    }

    /**
     * @brief 
     * @return 
     */
    virtual double cbminuscc() const
    {
        return 0.0;
    }

    /**
     * @brief 
     * @return 
     */
    virtual double cbminusctau() const
    {
        return 0.0;
    }

    /**
     * @brief 
     * @return 
     */
    virtual double ccminusctau() const
    {
        return 0.0;
    }

    ////////////////////////////////////////////////////////////////////////
      
    /**
     * @brief The new physics contribution to the anomalous triple gauge coupling @f$g_{1,Z}@f$.
     * @return @f$\delta g_{1,Z}@f$
     */
    virtual double deltag1ZNP() const
    {
        return 0.0;
    }
      
    /**
     * @brief The new physics contribution to the anomalous triple gauge coupling @f$\kappa_{\gamma}@f$.
     * @return @f$\delta \kappa_{\gamma}@f$
     */
    virtual double deltaKgammaNP() const
    {
        return 0.0;
    }
      
    /**
     * @brief The new physics contribution to the anomalous triple gauge coupling @f$\lambda_{Z}@f$.
     * @return @f$\lambda_{Z}@f$
     */
    virtual double lambdaZNP() const
    {
        return 0.0;
    }
    
    
    ////////////////////////////////////////////////////////////////////////
      
    /**
     * @brief The new physics contribution to the effective anomalous triple 
     * gauge coupling @f$g_{1,Z}^{Eff}@f$ from arXiv: 1708.09079 [hep-ph].
     * @return @f$\delta g_{1,Z}@f$
     */
    virtual double deltag1ZNPEff() const
    {
        return 0.0;
    }
      
    /**
     * @brief The new physics contribution to the effective anomalous triple 
     * gauge coupling @f$\kappa_{\gamma}^{Eff}@f$ from arXiv: 1708.09079 [hep-ph].
     * @return @f$\delta \kappa_{\gamma}@f$
     */
    virtual double deltaKgammaNPEff() const
    {
        return 0.0;
    }
    
    ////////////////////////////////////////////////////////////////////////
    /**
     * @brief The differential distribution for @f$e^+ e^- \to W^+ W^- \to jj \ell \nu@f$, 
     * with @f$\ell= e, \mu@f$, as a function of the @f$W@f$ polar angle.
     * @return @f$d\sigma/d\cos{\theta}@f$
     */
    virtual double dxseeWWdcos(const double sqrt_s, const double cos) const
    {
        return 0.0;
    }
    
    /**
     * @brief The integral of differential distribution for @f$e^+ e^- \to W^+ W^- \to jj \ell \nu@f$, 
     * with @f$\ell= e, \mu@f$ in a given bin of the @f$W@f$ polar angle.
     * @return @f$\int_{\cos{\theta_1}}^{\cos{\theta_2}} d\sigma/d\cos{\theta}@f$
     */
    virtual double dxseeWWdcosBin(const double sqrt_s, const double cos1, const double cos2) const
    {
        return 0.0;
    }
    
    /**
     * @brief Total @f$e^+ e^- \to W^+ W^- \to jj \ell \nu@f$ cross section in pb, 
     * with @f$\ell= e, \mu@f$.
     * @return @f$\sigma(e^+ e^- \to W^+ W^- \to jj \ell \nu) @f$
     */
    virtual double xseeWW(const double sqrt_s) const
    {
        return 0.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{eeWW}@f$ between the 
     * @f$ e^{+}e^{-}\to W^{+}W^{-} @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeWW}@f$
     */
    virtual double mueeWW(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{eeWW}@f$ between the 
     * @f$ e^{+}e^{-}\to W^{+}W^{-} @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeWW}@f$
     */
    virtual double mueeWWPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
    {
        return 1.0;
    }
    
    ////////////////////////////////////////////////////////////////////////
    
    //----- High Energy diboson observables at hadron colliders

    /**
     * @brief The direction constrained by @f$ p p \to Z H@f$ in the boosted regime, @f$g_p^Z@f$.
     * From arXiv:1807.01796 and the contribution to FCC CDR Vol 1. Implemented only in NPSMEFTd6 class.
     * @return @f$g_p^Z@f$
     */
    virtual double ppZHprobe(const double sqrt_s) const
    {
        return 0.0;
    }
    
    /**
     * @brief The number of events in  @f$ p p \to WZ@f$
     * in a given @f$p_{TV}@f$ bin, normalized to the SM prediction.
     * From arXiv: 1712.01310 [hep-ph] and private communication.
     * Implemented only in NPSMEFTd6 class.
     * @return @f$N_{ev}^{p_{TV}}/N_{ev,SM}^{p_{TV}}@f$
     */
    virtual double mupTVppWZ(const double sqrt_s, const double pTV1, const double pTV2) const
    {
        return 1.0;
    }
 
    ////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief The effective coupling @f$\kappa_{\mu,eff}=\sqrt{\Gamma_{H\mu\mu}/\Gamma_{H\mu\mu}^{SM}}@f$.
     * @return @f$\kappa_{\mu,eff}@f$
     */
    virtual double kappamueff() const
    {
        return 1.0;
    }
    
    /**
     * @brief The effective coupling @f$\kappa_{\tau,eff}=\sqrt{\Gamma_{H\tau\tau}/\Gamma_{H\tau\tau}^{SM}}@f$.
     * @return @f$\kappa_{\tau,eff}@f$
     */
    virtual double kappataueff() const
    {
        return 1.0;
    }
    
    /**
     * @brief The effective coupling @f$\kappa_{c,eff}=\sqrt{\Gamma_{Hcc}/\Gamma_{Hcc}^{SM}}@f$.
     * @return @f$\kappa_{c,eff}@f$
     */
    virtual double kappaceff() const
    {
        return 1.0;
    }
    
    /**
     * @brief The effective coupling @f$\kappa_{b,eff}=\sqrt{\Gamma_{Hbb}/\Gamma_{Hbb}^{SM}}@f$.
     * @return @f$\kappa_{b,eff}@f$
     */
    virtual double kappabeff() const
    {
        return 1.0;
    }
    
    /**
     * @brief The effective coupling @f$\kappa_{G,eff}=\sqrt{\Gamma_{HGG}/\Gamma_{HGG}^{SM}}@f$.
     * @return @f$\kappa_{G,eff}@f$
     */
    virtual double kappaGeff() const
    {
        return 1.0;
    }
    
    /**
     * @brief The effective coupling @f$\kappa_{Z,eff}=\sqrt{\Gamma_{HZZ}/\Gamma_{HZZ}^{SM}}@f$.
     * @return @f$\kappa_{Z,eff}@f$
     */
    virtual double kappaZeff() const
    {
        return 1.0;
    }
    
    /**
     * @brief The effective coupling @f$\kappa_{W,eff}=\sqrt{\Gamma_{HWW}/\Gamma_{HWW}^{SM}}@f$.
     * @return @f$\kappa_{W,eff}@f$
     */
    virtual double kappaWeff() const
    {
        return 1.0;
    }
    
    /**
     * @brief The effective coupling @f$\kappa_{A,eff}=\sqrt{\Gamma_{HAA}/\Gamma_{HAA}^{SM}}@f$.
     * @return @f$\kappa_{A,eff}@f$
     */
    virtual double kappaAeff() const
    {
        return 1.0;
    }
    
    /**
     * @brief The effective coupling @f$\kappa_{ZA,eff}=\sqrt{\Gamma_{HZA}/\Gamma_{HZA}^{SM}}@f$.
     * @return @f$\kappa_{ZA,eff}@f$
     */
    virtual double kappaZAeff() const
    {
        return 1.0;
    }
    
    /////////////Basic interactions of the so-called Higgs basis////////////////
    
    /**
     * @brief The Higgs-basis coupling @f$\delta y_t@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * @return @f$\delta y_t@f$
     */
    virtual double deltayt_HB() const
    {
        return 0.0;
    }
    
    /**
     * @brief The Higgs-basis coupling @f$\delta y_b@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * @return @f$\delta y_b@f$
     */
    virtual double deltayb_HB() const
    {
        return 0.0;
    }
    
    /**
     * @brief The Higgs-basis coupling @f$\delta y_\tau@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * @return @f$\delta y_\tau@f$
     */
    virtual double deltaytau_HB() const
    {
        return 0.0;
    }
    
    /**
     * @brief The Higgs-basis coupling @f$\delta y_c@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * @return @f$\delta y_c@f$
     */
    virtual double deltayc_HB() const
    {
        return 0.0;
    }
    
    
    /**
     * @brief The Higgs-basis coupling @f$\delta y_\mu@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * @return @f$\delta y_\mu@f$
     */
    virtual double deltaymu_HB() const
    {
        return 0.0;
    }
    
    /**
     * @brief The Higgs-basis coupling @f$\delta c_z@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * @return @f$\delta c_z@f$
     */
    virtual double deltacZ_HB() const
    {
        return 0.0;
    }
    
    /**
     * @brief The Higgs-basis coupling @f$c_{z\Box}@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * @return @f$c_{z\Box}@f$
     */
    virtual double cZBox_HB() const
    {
        return 0.0;
    }
    
    /**
     * @brief The Higgs-basis coupling @f$c_{zz}@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * @return @f$c_{zz}@f$
     */
    virtual double cZZ_HB() const
    {
        return 0.0;
    }
    
    /**
     * @brief The Higgs-basis coupling @f$c_{z\gamma}@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * @return @f$c_{z\gamma}@f$
     */
    virtual double cZga_HB() const
    {
        return 0.0;
    }
    
    /**
     * @brief The Higgs-basis coupling @f$c_{\gamma\gamma}@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * @return @f$c_{\gamma\gamma}@f$
     */
    virtual double cgaga_HB() const
    {
        return 0.0;
    }
    
    /**
     * @brief The Higgs-basis coupling @f$c_{gg}@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * @return @f$c_{gg}@f$
     */
    virtual double cgg_HB() const
    {
        return 0.0;
    }
    
    /**
     * @brief The Higgs-basis coupling @f$\lambda_{z}@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * @return @f$\lambda_{z}@f$
     */
    virtual double lambz_HB() const
    {
        return 0.0;
    }
    
    /////////////Auxiliary observables////////////////
    
    /**
     * @brief Auxiliary observable AuxObs_NP1
     * @return AuxObs_NP1
     */
    virtual double AuxObs_NP1() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP2
     * @return AuxObs_NP2
     */
    virtual double AuxObs_NP2() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP3
     * @return AuxObs_NP3
     */
    virtual double AuxObs_NP3() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP4
     * @return AuxObs_NP4
     */
    virtual double AuxObs_NP4() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP5
     * @return AuxObs_NP5
     */
    virtual double AuxObs_NP5() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP6
     * @return AuxObs_NP6
     */
    virtual double AuxObs_NP6() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP7
     * @return AuxObs_NP7
     */
    virtual double AuxObs_NP7() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP8
     * @return AuxObs_NP8
     */
    virtual double AuxObs_NP8() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP9
     * @return AuxObs_NP9
     */
    virtual double AuxObs_NP9() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP10
     * @return AuxObs_NP10
     */
    virtual double AuxObs_NP10() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP11
     * @return AuxObs_NP11
     */
    virtual double AuxObs_NP11() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP12
     * @return AuxObs_NP12
     */
    virtual double AuxObs_NP12() const
    {
        return 0.0;
    }
      
    ////////////////////////////////////////////////////////////////////////
protected:
    StandardModel trueSM;
};

#endif	/* NPBASE_H */

