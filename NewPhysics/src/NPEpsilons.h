/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPEPSILONS_H
#define	NPEPSILONS_H

#include "NPbase.h"

/**
 * @class NPEpsilons
 * @brief A model class for new physics in the form of contributions to the
 * \f$\varepsilon_{1,2,3,b}\f$ parameters.
 * @ingroup NewPhysics
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This is a Model class containing parameters and functions to work
 * with the epsilon paramterization 
 * \cite Altarelli:1990zd, \cite Altarelli:1991fk,\cite Altarelli:1993sz
 * for new physics contributions to the electroweak precision observables.
 * It is noted that the \f$\varepsilon_{1,2,3,b}\f$ parameters include both the
 * SM and new physics contributions.
 *
 *
 * @anchor NPEpsilonsInitialization
 * <h3>Initialization</h3>
 *
 * The constructor NPEpsilons() initializes the model flags explained below to
 * their default values. After creating an instance of the current class,
 * it is required to call the initialization method InitializeModel(). 
 * This pointer is then used in computing the \f$W\f$-boson mass and the
 * fermionic neutral-current couplings with the epsilon parameters.
 * In the Monte Carlo run, the constructor as well as the initialization
 * method are called in InputParser::ReadParameters().
 *
 *
 * @anchor NPEpsilonsParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %NPEpsilons are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%epsilon_1 </td>
 *   <td class="mod_symb">\f$\varepsilon_1\f$</td>
 *   <td class="mod_desc">The parameter \f$\varepsilon_1\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%epsilon_2 </td>
 *   <td class="mod_symb">\f$\varepsilon_2\f$</td>
 *   <td class="mod_desc">The parameter \f$\varepsilon_2\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%epsilon_3 </td>
 *   <td class="mod_symb">\f$\varepsilon_3\f$</td>
 *   <td class="mod_desc">The parameter \f$\varepsilon_3\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%epsilon_b </td>
 *   <td class="mod_symb">\f$\varepsilon_b\f$</td>
 *   <td class="mod_desc">The parameter \f$\varepsilon_b\f$.</td>
 * </tr>
 * </table>
 * 
 * @anchor NPEpsilonsFlags
 * <h3>%Model Flags</h3>
 *
 * The flags of %NPEpsilons are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>Value</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%epsilon1SM</td>
 *   <td class="mod_desc">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if only the SM contribution
 *   is considered for \f$\varepsilon_1\f$. The default value is FALSE.</td>
 * <tr>
 *   <td class="mod_name">%epsilon2SM</td>
 *   <td class="mod_desc">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if only the SM contribution
 *   is considered for \f$\varepsilon_2\f$. The default value is FALSE.</td>
 * <tr>
 *   <td class="mod_name">%epsilon3SM</td>
 *   <td class="mod_desc">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if only the SM contribution
 *   is considered for \f$\varepsilon_3\f$. The default value is FALSE.</td>
 * <tr>
 *   <td class="mod_name">%epsilonbSM</td>
 *   <td class="mod_desc">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if only the SM contribution
 *   is considered for \f$\varepsilon_b\f$. The default value is FALSE.</td>
 * </table>
 *
 *
 * @anchor NPEpsilonsFunctions
 * <h3>Important member functions</h3>
 *
 * Compared to the base class NPbase, the functions for the following quantities
 * are reimplemented in the current class:
 *
 * @li @f$M_W@f$&nbsp;&nbsp; (with Mw()),
 * @li @f$\Gamma_W@f$&nbsp;&nbsp; (with GammaW()),
 *
 * where the latter is not available, since it cannot be simply parameterized
 * by the epsilon parameters.
 * In addition, the functions for the epsilon parameters are also provided:
 * 
 * @li @f$\varepsilon_1@f$, @f$\varepsilon_2@f$, @f$\varepsilon_3@f$ and
 * @f$\varepsilon_b@f$&nbsp;&nbsp;
 * (with epsilon1(), epsilon2(), epsilon3() and epsilonb()).
 *
 */
class NPEpsilons : public NPbase {
public:

    /**
     * @brief The number of the model parameters in %NPEpsilons.
     */
    static const int NEPSILONvars = 4;

    /**
     * @brief A string array containing the labels of the model parameters in %NPEpsilons.
     */
    static const std::string EPSILONvars[NEPSILONvars];

    /**
     * @brief The default constructor.
     */
    NPEpsilons();

    /**
     * @brief The post-update method for %NPEpsilons.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PostUpdate();

    /**
     * @brief @copybrief Model::CheckParameters()
     * @copydetails Model::CheckParameters()
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    /**
     * @brief @copybrief Model::setFlag()
     * @copydetails Model::setFlag()
     */
    virtual bool setFlag(const std::string name, const bool value);

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The parameter \f$\varepsilon_1\f$.
     * @return the SM value (FlagEpsilon1SM=true) or the SM plus new physics
     * value (FlagEpsilon1SM=false) of \f$\varepsilon_1\f$
     */
    virtual double epsilon1() const;

    /**
     * @brief The parameter \f$\varepsilon_2\f$.
     * @return the SM value (FlagEpsilon2SM=true) or the SM plus new physics
     * value (FlagEpsilon2SM=false) of \f$\varepsilon_2\f$
     */
    virtual double epsilon2() const;

    /**
     * @brief The parameter \f$\varepsilon_3\f$.
     * @return the SM value (FlagEpsilon3SM=true) or the SM plus new physics
     * value (FlagEpsilon3SM=false) of \f$\varepsilon_3\f$
     */
    virtual double epsilon3() const;

    /**
     * @brief The parameter \f$\varepsilon_b\f$.
     * @return the SM value (FlagEpsilonbSM=true) or the SM plus new physics
     * value (FlagEpsilonbSM=false) of \f$\varepsilon_b\f$
     */
    virtual double epsilonb() const;


    ////////////////////////////////////////////////////////////////////////     

    /**
     * @brief The mass of the @f$W@f$ boson, @f$M_W@f$.
     * @details This function calls EWNPEpsilons::Mw() via
     * EWNPEpsilons::Mw_NPEpsilons().
     * @return @f$M_W@f$ in GeV
     */
    virtual double Mw() const;

    /** 
     * @brief The total width of the @f$W@f$ boson, @f$\Gamma_W@f$.
     * 
     * @warning This function is not available.
     */
    virtual double GammaW() const;

    /** 
     * @brief @copybrief StandardModel::Gamma_Z()
     * @return @f$\Gamma_Z@f$ in GeV
     */
    virtual double Gamma_Z() const;

    /** 
     * @brief @copybrief StandardModel::sigma0_had()
     * @return @f$\sigma_h^0@f$ in GeV@f$^{-2}@f$
     */
    virtual double sigma0_had() const;

    /** 
     * @brief @copybrief StandardModel::sin2thetaEff()
     * @param[in] f a lepton or quark
     * @return @f$\sin^2\theta_{\rm eff}^{\,\ell}@f$
     */
    virtual double sin2thetaEff(const Particle p) const;

    /** 
     * @brief @copybrief StandardModel::A_f()
     * @param[in] f a lepton or quark
     * @return @f$\mathcal{A}_\ell@f$
     */
    virtual double A_f(const Particle p) const;

    /** 
     * @brief The forward-backward asymmetry in @f$e^+e^-\to Z\to f \bar{f}@f$ at the
     * @f$Z@f$-pole, @f$A^f_{FB}@f$.
     * @param[in] f a lepton or quark
     * @return @f$A^f_{FB}@f$
     */
    virtual double AFB(const Particle p) const;

    /** 
     * @brief The ratio @f$R_\ell^0=\Gamma_{\mathrm{had}}/\Gamma_\ell@f$
     * or @f$R_q^0=\Gamma_q/\Gamma_{\mathrm{had}}@f$, for leptons or quarks, respectively. 
     * @param[in] f a lepton or quark
     * @return @f$R_f^0@f$
     */
    virtual double R0_f(const Particle p) const;

    /**
     * @brief The \f$W\f$ boson mass @f$M_W@f$. 
     * @return @f$M_W@f$ in GeV
     */
    double Mw_NPEpsilons() const;


    /**
     * @brief The @f$W@f$-boson mass @f$M_W@f$.
     * @details The radiative corrections to @f$M_W@f$ is parameterized in terms
     * of the quantity @f$\Delta r@f$:
     * @f[
     * M_W^2 = \frac{M_Z^2}{2} \left( 1+\sqrt{1-\frac{4\pi\alpha}{\sqrt{2}G_\mu M_Z^2(1-\Delta r)}}\ \right)\,,
     * @f]
     * where @f$\Delta r@f$ contains both SM and NP contributions, and is resummed. 
     * In @cite Altarelli:1990zd and @cite Altarelli:1991fk, @f$\Delta r@f$
     * is given in terms of @f$\Delta r_W@f$:
     * @f[
     * \Delta r = 1 - [1 - \Delta\alpha(M_Z^2)][1-\Delta r_W(\varepsilon_1,\varepsilon_2,\varepsilon_3)]\,,
     * @f]
     * where @f$\Delta r_W@f$ is a function of @f$\varepsilon_1@f$,
     * @f$\varepsilon_2@f$ and @f$\varepsilon_3@f$. See Delta_rW().
     * @param[in] eps1 the @f$\varepsilon_1@f$ parameter
     * @param[in] eps2 the @f$\varepsilon_2@f$ parameter
     * @param[in] eps3 the @f$\varepsilon_3@f$ parameter
     * @return @f$M_W@f$ in GeV
     */
    double Mw_eps(const double eps1, const double eps2, const double eps3) const;

    /**
     * @brief The effective neutral-current coupling @f$\rho_Z^f@f$.
     * @details The @f$Zf\bar{f}@f$ effective coupling @f$\rho_Z^f@f$ for
     * @f$f\neq b,t@f$ is flavour universal in the original version of the epsilon parameterization:
     * @f[
     *   \rho_Z^f = \rho_Z^e(\varepsilon_1).
     * @f]
     * See rhoZ_e() for @f$\rho_Z^e(\varepsilon_1)@f$.
     * When StandardModel::FlagWithoutNonUniversalVC is true, the flavour
     * non-universal vertex corrections are also taken into account:
     * @f[
     *   \rho_Z^f = \rho_Z^e(\varepsilon_1) + \Delta\rho_Z^f,
     * @f]
     * where @f$\Delta\rho_Z^f@f$ denotes the non-universal corrections
     * given by StandardModel::rhoZ_f_SM_FlavorDep().
     * @param[in] f a lepton or quark
     * @param[in] eps1 the @f$\varepsilon_1@f$ parameter
     * @param[in] epsb the @f$\varepsilon_b@f$ parameter
     * @return @f$\rho_Z^f@f$
     */
    gslpp::complex rhoZ_f_eps(const Particle f, const double eps1, const double epsb = 0.) const;

    /**
     * @brief @copybrief NPbase::rhoZ_f()
     * @details NP contribution is included via the \f$\varepsilon_i\f$ parameter
     * @copydetails NPbase::rhoZ_f()
     */
    virtual gslpp::complex rhoZ_f(const Particle f) const;

    /**
     * @brief The effective neutral-current coupling @f$\kappa_Z^f@f$. 
     * @details The @f$Zf\bar{f}@f$ effective coupling @f$\kappa_Z^f@f$ for
     * @f$f\neq b,t@f$ is flavour universal in the original version of the
     * epsilon parameterization:
     * @f[
     *   \kappa_Z^f = \kappa_Z^e(\varepsilon_1,\varepsilon_3).
     * @f]
     * See kappaZ_e() for @f$\kappa_Z^e(\varepsilon_1,\varepsilon_3)@f$.
     * When StandardModel::FlagWithoutNonUniversalVC is true, the flavour
     * non-universal vertex corrections are also taken into account:
     * @f[
     *   \kappa_Z^f = \kappa_Z^e(\varepsilon_1,\varepsilon_3) + \Delta\kappa_Z^f,
     * @f]
     * where @f$\Delta\kappa_Z^f@f$ denotes the non-universal corrections
     * given by StandardModel::kappaZ_f_SM_FlavorDep().
     * @param[in] f a lepton or quark
     * @param[in] eps1 the @f$\varepsilon_1@f$ parameter
     * @param[in] eps3 the @f$\varepsilon_3@f$ parameter
     * @param[in] epsb the @f$\varepsilon_b@f$ parameter
     * @return @f$\kappa_Z^f@f$
     */
    gslpp::complex kappaZ_f_eps(const Particle f, const double eps1, const double eps3, const double epsb = 0.) const;

    /**
     * @brief @copybrief NPbase::kappaZ_f()
     * @details NP contribution is included via the \f$\varepsilon_i\f$ parameter
     * @param[in] f a lepton or quark
     * @copydetails NPbase::kappaZ_f()
     */
    virtual gslpp::complex kappaZ_f(const Particle f) const;

    /**
     * @brief The effective neutral-current vector coupling @f$g_V^f@f$.
     * @details The coupling @f$g_V^e@f$ is given in terms of the epsilon
     * parameters @f$\varepsilon_1@f$ and @f$\varepsilon_3@f$:
     * @f[
     * g_V^e(\varepsilon_1,\varepsilon_3)
     * = \left\{ 1 - 4|Q_e|\,[1+\Delta\kappa'(\varepsilon_1,\varepsilon_3)] s_0^2 \right\}
     *   g_A^e(\varepsilon_1)\,.
     * @f]
     * See @cite Altarelli:1990zd and @cite Altarelli:1991fk.
     * @param[in] f a lepton or quark
     * @param[in] eps1 the @f$\varepsilon_1@f$ parameter
     * @param[in] eps3 the @f$\varepsilon_3@f$ parameter
     * @param[in] epsb the @f$\varepsilon_b@f$ parameter
     * @return @f$g_V^f@f$
     */
    gslpp::complex gV_f_eps(const Particle f, const double eps1, const double eps3, const double epsb = 0.) const;

    /**
     * @brief @copybrief NPbase::gV_f()
     * @details NP contribution is included via the \f$\varepsilon_i\f$ parameter
     * @copydetails NPbase::gV_f()
     */
    virtual gslpp::complex gV_f(const Particle f) const;

    /**
     * @brief The effective neutral-current axial-vector coupling @f$g_A^f@f$. 
     * @details The coupling @f$g_A^e@f$ is given in terms of the
     * @f$\varepsilon_1@f$ parameter:
     * @f[
     * g_A^e(\varepsilon_1)
     * = - \frac{1}{2}\left( 1 + \frac{\varepsilon_1}{2} \right).
     * @f]
     * See @cite Altarelli:1990zd and @cite Altarelli:1991fk.
     * @param[in] f a lepton or quark
     * @param[in] eps1 the @f$\varepsilon_1@f$ parameter
     * @param[in] epsb the @f$\varepsilon_b@f$ parameter
     * @return @f$g_A^f@f$
     */
    gslpp::complex gA_f_eps(const Particle f, const double eps1, const double epsb = 0.) const;

    /**
     * @brief @copybrief NPbase::gA_f()
     * @details NP contribution is included via the \f$\varepsilon_i\f$ parameter
     * @param[in] f a lepton or quark
     * @copydetails NPbase::gA_f()
     */
    virtual gslpp::complex gA_f(const Particle f) const;


    ////////////////////////////////////////////////////////////////////////   
protected:

    double myEpsilon_1; ///< The parameter \f$\varepsilon_1\f$.
    double myEpsilon_2; ///< The parameter \f$\varepsilon_2\f$.
    double myEpsilon_3; ///< The parameter \f$\varepsilon_3\f$.
    double myEpsilon_b; ///< The parameter \f$\varepsilon_b\f$.

    /**
     * @brief @copybrief Model::setParameter()
     * @copydetails Model::setParameter()
     */
    virtual void setParameter(const std::string name, const double& value);


    ////////////////////////////////////////////////////////////////////////         
private:

    bool FlagEpsilon1SM; ///< A boolean flag that is true if only the SM contribution is considered for \f$\varepsilon_1\f$.
    bool FlagEpsilon2SM; ///< A boolean flag that is true if only the SM contribution is considered for \f$\varepsilon_2\f$.
    bool FlagEpsilon3SM; ///< A boolean flag that is true if only the SM contribution is considered for \f$\varepsilon_3\f$.
    bool FlagEpsilonbSM; ///< A boolean flag that is true if only the SM contribution is considered for \f$\varepsilon_b\f$.

    /**
     * @brief The auxiliary function @f$\Delta r_W@f$.
     * @details The function @f$\Delta r_W@f$ is given in terms of the epsilon
     * parameters:
     * @f[
     * \Delta r_W(\varepsilon_1,\varepsilon_2,\varepsilon_3)
     *   = \frac{c_0^2-s_0^2}{s_0^2}
     *     \left[ \varepsilon_2 - c_0^2\,\varepsilon_1
     *       + 2s_0^2\Delta\kappa'(\varepsilon_1,\varepsilon_3) \right],
     * @f]
     * where @f$\Delta\kappa'(\varepsilon_1,\varepsilon_3)@f$ is defined as
     * the function Delta_kappaPrime().
     * See @cite Altarelli:1990zd and @cite Altarelli:1991fk.
     * @param[in] eps1 the @f$\varepsilon_1@f$ parameter
     * @param[in] eps2 the @f$\varepsilon_2@f$ parameter
     * @param[in] eps3 the @f$\varepsilon_3@f$ parameter
     * @return @f$\Delta r_W@f$
     */
    double Delta_rW(const double eps1, const double eps2, const double eps3) const;

    /**
     * @brief The auxiliary function @f$\Delta\kappa'@f$.
     * @details The function @f$\Delta\kappa'@f$ is given in terms of the epsilon
     * parameters:
     * @f[
     * \Delta\kappa'(\varepsilon_1,\varepsilon_3)
     * = \frac{\epsilon_3-c_0^2\,\epsilon_1}{c_0^2-s_0^2}\,.
     * @f]
     * See @cite Altarelli:1990zd and @cite Altarelli:1991fk.
     * @param[in] eps1 the @f$\varepsilon_1@f$ parameter
     * @param[in] eps3 the @f$\varepsilon_3@f$ parameter
     * @return @f$\Delta\kappa'@f$
     */
    double Delta_kappaPrime(const double eps1, const double eps3) const;


};

#endif	/* NPEPSILONS_H */

