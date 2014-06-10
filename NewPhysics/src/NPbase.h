/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPBASE_H
#define	NPBASE_H

#include <StandardModel.h>

/**
 * @class NPbase
 * @brief The auxiliary base model class for other model classes.
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
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
class NPbase : virtual public StandardModel {
public:

    /**
     * @brief Th default constructor.
     */
    NPbase();

    virtual bool PostUpdate();

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
    virtual double DeltaGF() const {
        return 0.; 
    }

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The oblique parameter \f$S\f$.
     * @return the value of @f$S@f$
     */
    virtual double obliqueS() const {
        return 0.; 
    }

    /**
     * @brief The oblique parameter \f$T\f$.
     * @return the value of @f$T@f$
     */
    virtual double obliqueT() const {
        return 0.; 
    }

    /**
     * @brief The oblique parameter \f$U\f$.
     * @return the value of @f$U@f$
     */
    virtual double obliqueU() const {
        return 0.; 
    }

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief @copybrief StandardModel::Mw()
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
     * @brief @copybrief StandardModel::GammaW()
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
     * @brief New physics contribution to @f$g_V^l@f$.
     * @details
     * The neutral-current vector coupling @f$g_V^l@f$ receives the new physics
     * contribution via the oblique parameters @f$S@f$ and @f$T@f$ and the shift
     * in the Fermi constant, @f$\Delta G@f$:
     * @f[
     * \delta g_V^l =
     * \frac{g_{V,\mathrm{SM}}^l}{2}
     * \left[ \alpha(M_Z^2)\, T - \Delta G \right]
     * +
     * \frac{\big( g_{V,\mathrm{SM}}^l - g_{A,\mathrm{SM}}^l \big)
     * \left[
     * \alpha(M_Z^2)\left( S - 4\,c_W^2s_W^2\, T \right)
     * + 4\,c_W^2s_W^2\, \Delta G
     * \right]}{4s_W^2\,(c_W^2-s_W^2)}\,.
     * @f]
     *
     * See @cite Ciuchini:2013pca and references therein.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$\delta g_V^l@f$
     */
    virtual double deltaGV_f(const Particle p) const;
    
    virtual complex gV_f(const Particle p) const;

    /**
     * @brief New physics contribution to @f$g_A^l@f$.
     * @details
     * The neutral-current axial-vector coupling @f$g_A^l@f$ receives the new
     * physics contribution via the oblique parameter @f$T@f$ and the shift in
     * the Fermi constant, @f$\Delta G@f$:
     * @f[
     * \delta g_A^l
     * = \frac{g_{A,\mathrm{SM}}^l}{2} \left[ \alpha(M_Z^2)\, T - \Delta G \right].
     * @f]
     *
     * See @cite Ciuchini:2013pca and references therein.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$\delta g_A^l@f$
     */   
    virtual double deltaGA_f(const Particle p) const;
    
    virtual complex gA_f(const Particle p) const;
    
       ////////////////////////////////////////////////////////////////////////

    virtual double deltaGamma_Z() const;
    
    /**
     * @brief The total decay width of the @f$Z@f$ boson, @f$\Gamma_Z@f$.
     * @param[in] GammaZ_SM the SM prediction for @f$\Gamma_Z@f$ in GeV
     * @return @f$\Gamma_Z@f$ in GeV, including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    
    virtual double Gamma_Z() const;

    /**
     * @brief The cross section for the process @f$e^+ e^-\to Z\to \mathrm{hadrons}@f$
     * at the @f$Z@f$ pole, @f$\sigma_h^0@f$.
     * @param[in] sigmaHadron_SM the SM prediction for @f$\sigma_h^0@f$ in GeV@f$^{-2}@f$
     * @return @f$\sigma_h^0@f$ in GeV@f$^{-2}@f$, including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    virtual double deltaSigmaHadron() const;

    virtual double sigma0_had() const;
    
    /**
     * @brief @copybrief sin2thetaEff::computeThValue()
     * @param[in] sin2thetaEff_SM the SM prediction for @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$
     * @return @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$, including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    virtual double deltaSin2thetaEff_e() const;
    
    virtual double sin2thetaEff(const Particle p) const;

    /**
     * @brief @copybrief PtauPol::computeThValue()
     * @param[in] PtauPol_SM the SM prediction for @f$P_\tau^{\mathrm{pol}}@f$
     * @return @f$P_\tau^{\mathrm{pol}}@f$, including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    virtual double deltaA_f(const Particle p) const;
    
    virtual double A_f(const Particle p) const;

    virtual double deltaAFB(const Particle p) const;
    virtual double AFB(const Particle p) const;


    /**
     * @brief @copybrief Rlepton::computeThValue()
     * @param[in] Rlepton_SM the SM prediction for @f$R_\ell^0@f$
     * @return @f$R_\ell^0@f$, including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    virtual double deltaR0_f(const Particle p) const;
    virtual double R0_f(const Particle p) const;

    ////////////////////////////////////////////////////////////////////////
protected:
        StandardModel trueSM;
};

#endif	/* NPBASE_H */

