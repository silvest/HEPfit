/*
 * Copyright (C) 2013-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWNPEPSILONS_H
#define	EWNPEPSILONS_H

#include <EWSM.h>

/**
 * @addtogroup NewPhysics
 * @brief A module for model-independent analysis of new physics.
 * @{
 */

/**
 * @class EWNPEpsilons
 * @brief A class for the \f$W\f$-boson mass and the fermionic neutral-current
 * couplings with new physics in the form of contributions to the
 * \f$\varepsilon_{1,2,3,b}\f$ parameters.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class contains the necessary functions to compute new physics 
 * contributions to electroweak precision observables with the
 * \f$\varepsilon\f$ parameterization \cite Altarelli:1990zd,
 * \cite Altarelli:1991fk, \cite Altarelli:1993sz :
 *
 * @li @f$M_W@f$&nbsp;&nbsp; (with Mw_NPEpsilons()),
 * @li @f$\rho_Z^f@f$&nbsp;&nbsp; (with rhoZ_l() and rhoZ_q()),
 * @li @f$\kappa_Z^f@f$&nbsp;&nbsp; (with kappaZ_l() and kappaZ_q()),
 * @li @f$g_V^f@f$&nbsp;&nbsp; (with gVl() and gVq()),
 * @li @f$g_A^f@f$&nbsp;&nbsp; (with gAl() and gAq()). 
 *
 * These functions call the other member functions in the current class which
 * are parameterized in terms of the \f$\varepsilon_i\f$ parameters.
 * The definitions of the effective fermionic neutral-current couplings 
 * @f$\rho_Z^f@f$, @f$\kappa_Z^f@f$, @f$g_V^f@f$ and @f$g_A^f@f$ can be found 
 * in the description of EWSM class.
 *
 * When the model flag @ref StandardModelFlags "WithoutNonUniversalVC" defined 
 * in StandardModel is set to true, the flavour non-universal vertex corrections 
 * which are not considered in the original version of the epsilon parameterization
 * are also taken into account. See @cite Ciuchini:2013pca for detail.
 */
class EWNPEpsilons : public EWSM {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    EWNPEpsilons(const StandardModel& SM_i);


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The \f$W\f$ boson mass @f$M_W@f$. 
     * @return @f$M_W@f$ in GeV
     */
    double Mw_NPEpsilons() const;


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief @copybrief EWSM::rhoZ_l()
     * @details NP contribution is included via the \f$\varepsilon_i\f$ parameter
     * @copydetails EWSM::rhoZ_l()
     */
    virtual complex rhoZ_l(const StandardModel::lepton l) const;

    /**
     * @brief @copybrief EWSM::rhoZ_q()
     * @details NP contribution is included via the \f$\varepsilon_i\f$ parameter
     * @copydetails EWSM::rhoZ_q()
     */
    virtual complex rhoZ_q(const StandardModel::quark q) const;

    /**
     * @brief @copybrief EWSM::kappaZ_l()
     * @details NP contribution is included via the \f$\varepsilon_i\f$ parameter
     * @copydetails EWSM::kappaZ_l()
     */
    virtual complex kappaZ_l(const StandardModel::lepton l) const;

    /**
     * @brief @copybrief EWSM::kappaZ_q()
     * @details NP contribution is included via the \f$\varepsilon_i\f$ parameter
     * @copydetails EWSM::kappaZ_q()
     */
    virtual complex kappaZ_q(const StandardModel::quark q) const;

    /**
     * @brief @copybrief EWSM::gVl()
     * @details NP contribution is included via the \f$\varepsilon_i\f$ parameter
     * @copydetails EWSM::gVl()
     */
    virtual complex gVl(const StandardModel::lepton l) const;

    /**
     * @brief @copybrief EWSM::gVq()
     * @details NP contribution is included via the \f$\varepsilon_i\f$ parameter
     * @copydetails EWSM::gVq()
     */
    virtual complex gVq(const StandardModel::quark q) const;

    /**
     * @brief @copybrief EWSM::gAl()
     * @details NP contribution is included via the \f$\varepsilon_i\f$ parameter
     * @copydetails EWSM::gAl()
     */
    virtual complex gAl(const StandardModel::lepton l) const;

    /**
     * @brief @copybrief EWSM::gAq()
     * @details NP contribution is included via the \f$\varepsilon_i\f$ parameter
     * @copydetails EWSM::gAq()
     */
    virtual complex gAq(const StandardModel::quark q) const;


    ////////////////////////////////////////////////////////////////////////

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
    double Mw(const double eps1, const double eps2, const double eps3) const;

    /**
     * @brief The effective leptonic neutral-current coupling @f$\rho_Z^l@f$.
     * @details The @f$Zl\bar{l}@f$ effective coupling @f$\rho_Z^l@f$ is flavour
     * universal in the original version of the epsilon parameterization:
     * @f[
     *   \rho_Z^l = \rho_Z^e(\varepsilon_1).
     * @f]
     * See rhoZ_e() for @f$\rho_Z^e(\varepsilon_1)@f$.
     * When StandardModel::FlagWithoutNonUniversalVC is true, the flavour
     * non-universal vertex corrections are also taken into account:
     * @f[
     *   \rho_Z^l = \rho_Z^e(\varepsilon_1) + \Delta\rho_Z^l,
     * @f]
     * where @f$\Delta\rho_Z^l@f$ denotes the non-universal corrections
     * given by EWSM::rhoZ_l_SM_FlavorDep().
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @param[in] eps1 the @f$\varepsilon_1@f$ parameter
     * @return @f$\rho_Z^l@f$
     */
    complex rhoZ_l(const StandardModel::lepton l, const double eps1) const;

    /**
     * @brief The effective quark neutral-current coupling @f$\rho_Z^q@f$ for @f$q\neq b,t@f$.
     * @details The @f$Zq\bar{q}@f$ effective coupling @f$\rho_Z^q@f$ for
     * @f$q\neq b,t@f$ is flavour universal in the original version of the
     * epsilon parameterization:
     * @f[
     *   \rho_Z^q = \rho_Z^e(\varepsilon_1).
     * @f]
     * See rhoZ_e() for @f$\rho_Z^e(\varepsilon_1)@f$.
     * When StandardModel::FlagWithoutNonUniversalVC is true, the flavour
     * non-universal vertex corrections are also taken into account:
     * @f[
     *   \rho_Z^q = \rho_Z^e(\varepsilon_1) + \Delta\rho_Z^q,
     * @f]
     * where @f$\Delta\rho_Z^q@f$ denotes the non-universal corrections
     * given by EWSM::rhoZ_q_SM_FlavorDep().
     * @param[in] q name of a quark (see QCD::quark); @f$q\neq b,t@f$
     * @param[in] eps1 the @f$\varepsilon_1@f$ parameter
     * @return @f$\rho_Z^q@f$
     */
    complex rhoZ_q(const StandardModel::quark q, const double eps1) const;

    /**
     * @brief The effective leptonic neutral-current coupling @f$\kappa_Z^l@f$.
     * @details The @f$Zl\bar{l}@f$ effective coupling @f$\kappa_Z^l@f$ is flavour
     * universal in the original version of the epsilon parameterization:
     * @f[
     *   \kappa_Z^l = \kappa_Z^e(\varepsilon_1,\varepsilon_3).
     * @f]
     * See kappaZ_e() for @f$\kappa_Z^e(\varepsilon_1,\varepsilon_3)@f$.
     * When StandardModel::FlagWithoutNonUniversalVC is true, the flavour
     * non-universal vertex corrections are also taken into account:
     * @f[
     *   \kappa_Z^l = \kappa_Z^e(\varepsilon_1,\varepsilon_3) + \Delta\kappa_Z^l,
     * @f]
     * where @f$\Delta\kappa_Z^l@f$ denotes the non-universal corrections
     * given by EWSM::kappaZ_l_SM_FlavorDep().
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @param[in] eps1 the @f$\varepsilon_1@f$ parameter
     * @param[in] eps3 the @f$\varepsilon_3@f$ parameter
     * @return @f$\kappa_Z^l@f$
     */
    complex kappaZ_l(const StandardModel::lepton l,
                     const double eps1, const double eps3) const;

    /**
     * @brief The effective quark neutral-current coupling @f$\kappa_Z^q@f$ for @f$q\neq b,t@f$.
     * @details The @f$Zq\bar{q}@f$ effective coupling @f$\kappa_Z^q@f$ for
     * @f$q\neq b,t@f$ is flavour universal in the original version of the
     * epsilon parameterization:
     * @f[
     *   \kappa_Z^q = \kappa_Z^e(\varepsilon_1,\varepsilon_3).
     * @f]
     * See kappaZ_e() for @f$\kappa_Z^e(\varepsilon_1,\varepsilon_3)@f$.
     * When StandardModel::FlagWithoutNonUniversalVC is true, the flavour
     * non-universal vertex corrections are also taken into account:
     * @f[
     *   \kappa_Z^q = \kappa_Z^e(\varepsilon_1,\varepsilon_3) + \Delta\kappa_Z^q,
     * @f]
     * where @f$\Delta\kappa_Z^q@f$ denotes the non-universal corrections
     * given by EWSM::kappaZ_q_SM_FlavorDep().
     * @param[in] q name of a quark (see QCD::quark); @f$q\neq b,t@f$
     * @param[in] eps1 the @f$\varepsilon_1@f$ parameter
     * @param[in] eps3 the @f$\varepsilon_3@f$ parameter
     * @return @f$\kappa_Z^q@f$
     */
    complex kappaZ_q(const StandardModel::quark q,
                     const double eps1, const double eps3) const;

    /**
     * @brief The effective leptonic neutral-current vector coupling @f$g_V^l@f$.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @param[in] eps1 the @f$\varepsilon_1@f$ parameter
     * @param[in] eps3 the @f$\varepsilon_3@f$ parameter
     * @return @f$g_V^l@f$
     */
    complex gVl(const StandardModel::lepton l,
                const double eps1, const double eps3) const;

    /**
     * @brief The effective quark neutral-current vector coupling @f$g_V^q@f$ for @f$q\neq b,t@f$.
     * @param[in] q name of a quark (see QCD::quark); @f$q\neq b,t@f$
     * @param[in] eps1 the @f$\varepsilon_1@f$ parameter
     * @param[in] eps3 the @f$\varepsilon_3@f$ parameter
     * @return @f$g_V^q@f$
     */
    complex gVq(const StandardModel::quark q,
                const double eps1, const double eps3) const;

    /**
     * @brief The effective leptonic neutral-current axial-vector coupling @f$g_A^l@f$.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @param[in] eps1 the @f$\varepsilon_1@f$ parameter
     * @return @f$g_A^l@f$
     */
    complex gAl(const StandardModel::lepton l, const double eps1) const;

    /**
     * @brief The effective quark neutral-current axial-vector coupling @f$g_A^q@f$ for @f$q\neq b,t@f$.
     * @param[in] q name of a quark (see QCD::quark); @f$q\neq b,t@f$
     * @param[in] eps1 the @f$\varepsilon_1@f$ parameter
     * @return @f$g_A^q@f$
     */
    complex gAq(const StandardModel::quark q, const double eps1) const;

    /**
     * @brief The effective neutral-current coupling @f$\rho_Z^b@f$.
     * @details In the original version of the epsilon parameterization
     * @cite Altarelli:1993sz, the @f$Zb\bar{b}@f$ effective coupling
     * @f$\rho_Z^b@f$ is given in terms of @f$\varepsilon_b@f$:
     * @f[
     *   \rho_Z^b = \rho_Z^{e}(\varepsilon_1)(1 + \varepsilon_b)^2\,,
     * @f]
     * See rhoZ_e() for @f$\rho_Z^e(\varepsilon_1)@f$.
     * When StandardModel::FlagWithoutNonUniversalVC is true, the additional
     * flavour non-universal vertex corrections are taken into account:
     * @f[
     *   \rho_Z^b = \left[\rho_Z^{e}(\varepsilon_1) + \Delta\rho_Z^b \right]
     *   (1 + \varepsilon_b)^2\,,
     * @f]
     * where @f$\Delta\rho_Z^b@f$ denotes the non-universal corrections which are
     * not considered in @f$\varepsilon_b@f$.
     * @param[in] eps1 the @f$\varepsilon_1@f$ parameter
     * @param[in] epsb the @f$\varepsilon_b@f$ parameter
     * @return @f$\rho_Z^b@f$
     */
    complex rhoZ_b(const double eps1, const double epsb) const;

    /**
     * @brief The effective neutral-current coupling @f$\kappa_Z^b@f$.
     * @details In the original version of the epsilon parameterization
     * @cite Altarelli:1993sz, the @f$Zb\bar{b}@f$ effective coupling
     * @f$\kappa_Z^b@f$ is given in terms of @f$\varepsilon_b@f$:
     * @f[
     *   \kappa_Z^b = \frac{\kappa_Z^{e}(\varepsilon_1,\varepsilon_3)}{1+\varepsilon_b}\,,
     * @f]
     * See kappaZ_e() for @f$\kappa_Z^e(\varepsilon_1,\varepsilon_3)@f$.
     * When StandardModel::FlagWithoutNonUniversalVC is true, the additional
     * flavour non-universal vertex corrections are taken into account:
     * @f[
     *   \kappa_Z^b = \frac{\kappa_Z^{e}(\varepsilon_1,\varepsilon_3) + \Delta\kappa_Z^b}{1+\varepsilon_b}\,,
     * @f]
     * where @f$\Delta\kappa_Z^b@f$ denotes the non-universal corrections which are
     * not considered in @f$\varepsilon_b@f$.
     * @param[in] eps1 the @f$\varepsilon_1@f$ parameter
     * @param[in] eps3 the @f$\varepsilon_3@f$ parameter
     * @param[in] epsb the @f$\varepsilon_b@f$ parameter
     * @return @f$\rho_Z^b@f$
     */
    complex kappaZ_b(const double eps1, const double eps3, const double epsb) const;

    /**
     * @brief The effective neutral-current vector coupling @f$g_V^b@f$.
     * @param[in] eps1 the @f$\varepsilon_1@f$ parameter
     * @param[in] eps3 the @f$\varepsilon_3@f$ parameter
     * @param[in] epsb the @f$\varepsilon_b@f$ parameter
     * @return @f$g_V^b@f$
     */
    complex gVb(const double eps1, const double eps3, const double epsb) const;

    /**
     * @brief The effective neutral-current axial-vector coupling @f$g_A^b@f$.
     * @param[in] eps1 the @f$\varepsilon_1@f$ parameter
     * @param[in] epsb the @f$\varepsilon_b@f$ parameter
     * @return @f$g_A^b@f$
     */
    complex gAb(const double eps1, const double epsb) const;


    ////////////////////////////////////////////////////////////////////////
private:

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

    /**
     * @brief The effective neutral-current coupling @f$\rho_Z^e@f$.
     * @param[in] eps1 the @f$\varepsilon_1@f$ parameter
     * @return @f$\rho_Z^e@f$
     */
    complex rhoZ_e(const double eps1) const;

    /**
     * @brief The effective coupling @f$\kappa_Z^e@f$.
     * @param[in] eps1 the @f$\varepsilon_1@f$ parameter
     * @param[in] eps3 the @f$\varepsilon_3@f$ parameter
     * @return @f$\kappa_Z^e@f$
     */
    complex kappaZ_e(const double eps1, const double eps3) const;

    /**
     * @brief The effective neutral-current vector coupling @f$g_V^e@f$.
     * @details The coupling @f$g_V^e@f$ is given in terms of the epsilon
     * parameters @f$\varepsilon_1@f$ and @f$\varepsilon_3@f$:
     * @f[
     * g_V^e(\varepsilon_1,\varepsilon_3)
     * = \left\{ 1 - 4|Q_e|\,[1+\Delta\kappa'(\varepsilon_1,\varepsilon_3)] s_0^2 \right\}
     *   g_A^e(\varepsilon_1)\,.
     * @f]
     * See @cite Altarelli:1990zd and @cite Altarelli:1991fk.
     * @param[in] eps1 the @f$\varepsilon_1@f$ parameter
     * @param[in] eps3 the @f$\varepsilon_3@f$ parameter
     * @return @f$g_V^e@f$
     */
    complex gVe(const double eps1, const double eps3) const;

    /**
     * @brief The effective neutral-current axial-vector coupling @f$g_A^e@f$.
     * @details The coupling @f$g_A^e@f$ is given in terms of the
     * @f$\varepsilon_1@f$ parameter:
     * @f[
     * g_A^e(\varepsilon_1)
     * = - \frac{1}{2}\left( 1 + \frac{\varepsilon_1}{2} \right).
     * @f]
     * See @cite Altarelli:1990zd and @cite Altarelli:1991fk.
     * @param[in] eps1 the @f$\varepsilon_1@f$ parameter
     * @return @f$g_A^e@f$
     */
    complex gAe(const double eps1) const;


};

/**
 * @}
 */

#endif	/* EWNPEPSILONS_H */

