/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSMTWOLOOPEW_H
#define	EWSMTWOLOOPEW_H

#include "EWSMcache.h"
#include "EWSMOneLoopEW.h"


/**
 * @class EWSMTwoLoopEW
 * @ingroup StandardModel
 * @brief A class for @f$O(\alpha^2)@f$ two-loop corrections to the %EW
 * precision observables.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class handles two-loop %EW contributions of
 * @f$O(\alpha^2)@f$ to the following quantities, which are relevant to the %EW
 * precision observables: 
 *
 * @li @f$\Delta\alpha_{\mathrm{lept}}(M_Z^2)@f$&nbsp;&nbsp; (with DeltaAlpha_l()),
 * @li @f$\Delta\alpha_{\mathrm{top}}(M_Z^2)@f$&nbsp;&nbsp; (with DeltaAlpha_t()),
 * @li @f$\Delta\rho@f$&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; (with DeltaRho()),
 * @li @f$\Delta r_{\mathrm{rem}}@f$&nbsp;&nbsp; (with DeltaR_rem()),
 * @li @f$\delta\rho_{\mathrm{rem}}^{f}@f$&nbsp;&nbsp; (with deltaRho_rem_l() and deltaRho_rem_q()),
 * @li @f$\delta\kappa_{\mathrm{rem}}^{f}@f$&nbsp;&nbsp; (with deltaKappa_rem_l() and deltaKappa_rem_q()),
 *
 * and the @f$O(\alpha^2)@f$ corrections to @f$\Delta\rho@f$ and to @f$Zb\bar{b}@f$:
 *
 * @li @f$\rho^{(2)}@f$&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; (with rho_2()),
 * @li @f$\tau^{(2)}@f$&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; (with tau_2()).
 *
 * See also the description of EWSM class for their definitions.
 * The @f$O(\alpha^2)@f$ two-loop %EW contributions to the vacuum
 * polarization amplitudes of the gauge bosons were calculated in
 * @cite Barbieri:1992nz, @cite Barbieri:1992dq, @cite Fleischer:1993ub,
 * @cite Fleischer:1994cb, @cite Degrassi:1996mg, @cite Degrassi:1996ps and
 * @cite Degrassi:1999jd with large-@f$m_t@f$ expansion.
 * In the current class, the @f$O(\alpha^2(m_t^4/m_Z^4 + m_t^2/M_Z^2))@f$
 * corrections to @f$\Delta\rho@f$, @f$\Delta r_{\mathrm{rem}}@f$,
 * @f$\delta\rho_{\mathrm{rem}}^{f}@f$ and @f$\delta\kappa_{\mathrm{rem}}^{f}@f$
 * and the @f$O(\alpha^2 m_t^4/M_Z^4)@f$ corrections to @f$Zb\bar{b}@f$,
 * denoted by @f$\rho^{(2)}@f$ and @f$\tau^{(2)}@f$, are computed with the
 * auxiliary functions defined as private members.
 * In @cite Degrassi:1996mg, the former corrections were calculated in
 * the MSbar scheme in order to undertake resummations correctly.
 * In subsequent papers @cite Degrassi:1996ps and @cite Degrassi:1999jd,
 * the resultant two-loop contributions were rewritten in terms of parameters
 * in the on-shell scheme by taking into account additional contributions, 
 * which correspond to the member functions with the word "Add".
 */
class EWSMTwoLoopEW {
public:

    /**
     * @brief Constructor.
     * @param[in] cache_i a reference to an object of type EWSMcache
     */
    EWSMTwoLoopEW(const EWSMcache& cache_i);


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief Leptonic contribution of @f$O(\alpha^2)@f$
     * to the electromagnetic coupling @f$\alpha@f$,
     * denoted as @f$\Delta\alpha_{\mathrm{lept}}^{\alpha^2}(s)@f$.
     * @details The expressions used here can be found in @cite Steinhauser:1998rq.
     * @param[in] s invariant mass squared
     * @return @f$\Delta\alpha_{\mathrm{lept}}^{\alpha^2}(s)@f$
     *
     * @attention This function is valid in the limit of @f$s\gg m_l^2@f$.
     */
    double DeltaAlpha_l(const double s) const;

    /**
     * @brief Top-quark contribution of @f$O(\alpha^2)@f$
     * to the electromagnetic coupling @f$\alpha@f$,
     * denoted as @f$\Delta\alpha_{\mathrm{top}}^{\alpha^2}(s)@f$.
     * @details This contribution is not implemented, since it is tiny and negligible.
     * @param[in] s invariant mass squared
     * @return @f$\Delta\alpha_{\mathrm{top}}^{\alpha^2}(s)=0@f$
     */
    double DeltaAlpha_t(const double s) const;

    /**
     * @brief Leading two-loop contribution of @f$O(\alpha^2)@f$
     * to @f$\Delta\rho@f$, denoted as @f$\Delta\rho^{\alpha^2}@f$. 
     * @details This function handles the leading irreducible two-loop %EW
     * contribution of @f$O(\alpha^2(m_t^4/m_Z^4 + m_t^2/M_Z^2))@f$ to @f$\Delta\rho@f$
     * in the on-shell scheme. 
     * The expression can be found in
     * @cite Degrassi:1996mg and @cite Degrassi:1996ps :
     * @f[
     * \Delta\rho^{\alpha^2} =
     *  3 (X_t^{\alpha})^2
     *  \left( \Delta\hat{\rho}^{(2)}
     *  + 4\, {\it zt}\, c_W^2 \Delta\bar{\rho}_{\mathrm{add}}^{(2)} \right)
     *  - \left(\frac{\alpha}{4\pi}\right)^2 \frac{c_W^2}{s_W^2}
     *    \left[ \mathrm{Re}\Pi^{\mathrm{fer}}_{Z\gamma}(M_Z^2) \right]^2,
     * @f]
     * where the definitions of the symbols can be read from the codes below,
     * and the last term originates from the @f$Z@f$-@f$\gamma@f$ mixing
     * (see, e.g., Chapter 6 of @cite Bardin:1999ak).
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\Delta\rho^{\alpha^2}@f$
     */
    double DeltaRho(const double Mw_i) const;

    /**
     * @brief Remainder contribution of @f$O(\alpha^2)@f$ to @f$\Delta r@f$,
     * denoted as @f$\Delta r_{\mathrm{rem}}^{\alpha^2}@f$. 
     * @details This function handles the remainder two-loop %EW
     * contribution of @f$O(\alpha^2(m_t^4/m_Z^4 + m_t^2/M_Z^2))@f$ to @f$\Delta r@f$
     * in the on-shell scheme.
     * The expression can be found in
     * @cite Degrassi:1996mg, @cite Degrassi:1996ps and @cite Degrassi:1999jd :
     * @f[
     * \Delta r_{\rm rem}^{\alpha^2}
     * =
     * 3\left(\frac{\alpha}{4\pi s_W^2}\right)^2 \frac{m_t^2}{M_W^2}
     * \left[ \Delta \hat{r}_W^{(2)}
     * + s_W^2 \bigg(\frac{\delta e}{e}\bigg)^{(2)}
     * + \frac{1}{4}\, \bar{f}_{\rm add}^{(2)} \right], 
     * @f]
     * where the definitions of the symbols can be read from the codes below. 
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\Delta r_{\mathrm{rem}}^{\alpha^2}@f$
     */
    double DeltaR_rem(const double Mw_i) const;

    /**
     * @brief Remainder contribution of @f$O(\alpha^2)@f$ to the effective
     * couplings @f$\rho_Z^f@f$,
     * denoted as @f$\delta\rho_{\mathrm{rem}}^{f,\, \alpha^2}@f$.
     * @details This function handles the @f$O(\alpha^2)@f$ remainder contribution
     * to @f$\rho_{Z}^{f}@f$ in the on-shell scheme, which was calculated
     * in @cite Degrassi:1999jd :
     * @f[
     * \delta\rho_{\rm rem}^{f,\, \alpha^2} = 3 (X_t^{\alpha})^2
     * \left[ 16\, {\it zt}\,c_W^2\, \Delta\hat{\eta}^{(2)}
     * +  4\, {\it zt}\,c_W^2\, \Delta\bar{\eta}_{\rm add}^{(2)} \right],
     * @f]
     * where the definitions of the symbols can be read from the codes below. 
     * @param[in] f a lepton or quark
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\delta\rho_{\mathrm{rem}}^{f,\, \alpha^2}@f$
     */
    gslpp::complex deltaRho_rem_f(const Particle f, const double Mw_i) const;

    /**
     * @brief Remainder contribution of @f$O(\alpha^2)@f$ to the effective
     * couplings @f$\kappa_Z^f@f$,
     * denoted as @f$\delta\kappa_{\mathrm{rem}}^{f,\, \alpha^2}@f$.
     * @details This function handles the @f$O(\alpha^2)@f$ remainder contribution
     * to @f$\kappa_{Z}^{f}@f$ in the on-shell scheme, which was calculated
     * in @cite Degrassi:1999jd :
     * @f[
     * \delta\kappa_{\rm rem}^{f,\, \alpha^2} = 3 (X_t^{\alpha})^2
     * \left[ 16\, {\it zt}\,c_W^2\, \Delta\hat{k}^{(2)}
     * + 4\, {\it zt}\,c_W^2\, \Delta\bar{k}_{\rm add}^{(2)} \right],
     * @f]
     * where the definitions of the symbols can be read from the codes below.
     * @param[in] f a lepton or quark
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\delta\kappa_{\mathrm{rem}}^{f,\, \alpha^2}@f$
     */
    gslpp::complex deltaKappa_rem_f(const Particle f, const double Mw_i) const;

    ////////////////////////////////////////////////////////////////////////        
    // O(GF^2 Mt^2) contributions

    /**
     * @brief The function @f$\rho^{(2)}@f$. 
     * @details This function parameterize the @f$O(\alpha^2 m_t^4/M_Z^4)@f$
     * contribution to @f$\Delta\rho@f$: 
     * @f[
     * \Delta\rho^{\alpha^2} = 3 (X_t^{\alpha})^2\rho^{(2)},
     * @f]
     * where the expression of @f$\rho^{(2)}@f$ can be found in
     * @cite Fleischer:1993ub and @cite Fleischer:1994cb
     * (see also @cite Barbieri:1992nz and @cite Barbieri:1992dq). 
     * @return @f$\rho^{(2)}@f$
     */
    double rho_2() const;

    /**
     * @brief The function @f$\tau^{(2)}@f$.
     * @details This function parmeterize the @f$O(\alpha^2 m_t^4/M_Z^4)@f$
     * contribution to the @f$Zb\bar{b}@f$ vertex (see EWSM::taub()),
     * where the expression of @f$\tau^{(2)}@f$ can be found in
     * @cite Fleischer:1993ub and @cite Fleischer:1994cb
     * (see also @cite Barbieri:1992nz and @cite Barbieri:1992dq).
     * @return @f$\tau^{(2)}@f$
     */
    double tau_2() const;


    ////////////////////////////////////////////////////////////////////////        

private:
    const EWSMcache& cache; ///< A reference to an object of type EWSMcache.
    const EWSMOneLoopEW myOneLoopEW; ///< An object of type EWSMOneLoopEW.

    /**
     * @brief The auxiliary function @f$g(a)@f$.
     * @details See @cite Fleischer:1993ub and @cite Fleischer:1994cb.
     * @param[in] a the ratio @f$a=(m_h/m_t)^2@f$
     * @return @f$g(a)@f$
     */
    double g(const double a) const;

    /**
     * @brief The auxiliary function @f$f(a,0)@f$.
     * @details See @cite Fleischer:1993ub and @cite Fleischer:1994cb.
     * @param[in] a the ratio @f$a=(m_h/m_t)^2@f$
     * @return @f$f(a,0)@f$
     */
    double f0(const double a) const;

    /**
     * @brief The auxiliary function @f$f(a,1)@f$.
     * @details See @cite Fleischer:1993ub and @cite Fleischer:1994cb.
     * @param[in] a the ratio @f$a=(m_h/m_t)^2@f$
     * @return @f$f(a,1)@f$
     */
    double f1(const double a) const;


    ////////////////////////////////////////////////////////////////////////        

    /**
     * @brief The auxiliary function @f$\Delta\hat{\rho}^{(2)}@f$.
     * @details See @cite Degrassi:1996mg. 
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\Delta\hat{\rho}^{(2)}@f$
     */
    double DeltaRho2(const double Mw_i) const;

    /**
     * @brief The auxiliary function @f$\Delta\bar{\rho}_{\mathrm{add}}^{(2)}@f$.
     * @details See @cite Degrassi:1996ps.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\Delta\bar{\rho}_{\mathrm{add}}^{(2)}@f$
     */
    double DeltaRho2Add(const double Mw_i) const;

    /**
     * @brief The auxiliary function @f$\Delta \hat{r}_W^{(2)}@f$.
     * @details See @cite Degrassi:1996mg.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\Delta \hat{r}_W^{(2)}@f$
     */
    double DeltaRw2(const double Mw_i) const;

    /**
     * @brief The auxiliary function @f$(\delta e/e)^{(2)}@f$.
     * @details See @cite Degrassi:1999jd. 
     * @return @f$\displaystyle\bigg(\frac{\delta e}{e}\bigg)^{(2)}@f$
     */
    double deltaEoverE2() const;

    /**
     * @brief The auxiliary function @f$\bar{f}_{\rm add}^{(2)}@f$.
     * @details See @cite Degrassi:1996ps.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\bar{f}_{\rm add}^{(2)}@f$
     */
    double f2Add(const double Mw_i) const;

    /**
     * @brief The auxiliary function @f$\Delta\hat{\eta}^{(2)}@f$.
     * @details See @cite Degrassi:1999jd.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\Delta\hat{\eta}^{(2)}@f$
     */
    double DeltaEta2(const double Mw_i) const;

    /**
     * @brief The auxiliary function @f$\Delta\bar{\eta}_{\rm add}^{(2)}@f$.
     * @details This functions is used in DeltaEta2Add_l() and DeltaEta2Add_q().
     * See @cite Degrassi:1999jd.
     * @param[in] I3f the isospin of a final-state fermion
     * @param[in] Qf the electric charge of a final-state fermion
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\Delta\bar{\eta}_{\rm add}^{(2)}@f$
     */
    gslpp::complex DeltaEta2Add_tmp(const double I3f, const double Qf, const double Mw_i) const;

    /**
     * @brief The auxiliary function @f$\Delta\bar{\eta}_{\rm add}^{(2)}@f$ for @f$Z\to f\bar{f}@f$.
     * @details See @cite Degrassi:1999jd.
     * @param[in] f a lepton or quark
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\Delta\bar{\eta}_{\rm add}^{(2)}@f$
     */
    gslpp::complex DeltaEta2Add_f(const Particle f, const double Mw_i) const;

    /**
     * @brief The auxiliary function @f$\Delta\hat{\kappa}^{(2)}@f$.
     * @details See @cite Degrassi:1999jd.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\Delta\hat{\kappa}^{(2)}@f$
     */
    double DeltaKappa2(const double Mw_i) const;

    /**
     * @brief The auxiliary function @f$\Delta\bar{\kappa}_{\rm add}^{(2)}@f$.
     * @details This functions is used in DeltaKappa2Add_l() and DeltaKappa2Add_q().
     * See @cite Degrassi:1999jd.
     * @param[in] I3f the isospin of a final-state fermion
     * @param[in] Qf the electric charge of a final-state fermion
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\Delta\bar{\kappa}_{\rm add}^{(2)}@f$
     */
    gslpp::complex DeltaKappa2Add_tmp(const double I3f, const double Qf, const double Mw_i) const;

    /**
     * @brief The auxiliary function @f$\Delta\bar{\kappa}_{\rm add}^{(2)}@f$ for @f$Z\to f\bar{f}@f$.
     * @details See @cite Degrassi:1999jd.
     * @param[in] f a lepton or quark
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\Delta\bar{\kappa}_{\rm add}^{(2)}@f$
     */
    gslpp::complex DeltaKappa2Add_f(const Particle f, const double Mw_i) const;

    /**
     * @brief The auxiliary function @f$V_{\rm add}@f$.
     * @details This functions is used in DeltaEta2Add_tmp().
     * See @cite Degrassi:1999jd.
     * @param[in] I3f the isospin of a final-state fermion
     * @param[in] Qf the electric charge of a final-state fermion
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$V_{\rm add}@f$
     */
    gslpp::complex Vadd(const double I3f, const double Qf, const double Mw_i) const;

    /**
     * @brief The auxiliary function @f$\Delta\bar{\eta}^{(1)}_f@f$.
     * @details This functions is used in DeltaEta2Add_tmp().
     * See @cite Degrassi:1990ec and @cite Degrassi:1999jd.
     * @param[in] I3f the isospin of a final-state fermion
     * @param[in] Qf the electric charge of a final-state fermion
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\Delta\bar{\eta}^{(1)}_f@f$
     */
    gslpp::complex DeltaEtaf1(const double I3f, const double Qf, const double Mw_i) const;

    /**
     * @brief The auxiliary function @f${\cal V}_{fi}(q^2)@f$. 
     * @details This function is used in DeltaEtaf1().
     * See @cite Degrassi:1990ec.
     * @param[in] I3f the isospin of a final-state fermion
     * @param[in] Qf the electric charge of a final-state fermion
     * @param[in] q2 invariant mass squared
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f${\cal V}_{fi}(q^2)@f$
     */
    gslpp::complex Vfi(const double I3f, const double Qf, const double q2, const double Mw_i) const;

    /**
     * @brief The auxiliary function @f$\Lambda(x)@f$.
     * @details This functions is used in DeltaRho2(), DeltaRw2() and DeltaEta2().
     * See @cite Degrassi:1996mg and @cite Degrassi:1999jd.
     * @param[in] x a real variable
     * @return @f$\Lambda(x)@f$
     *
     * @attention This function is valid for @f$x\geq 0@f$.
     */
    double Lambda(const double x) const;

    /**
     * @brief The auxiliary function @f$\phi(x)@f$.
     * @details This functions is used in DeltaRho2(), DeltaRw2(), deltaEoverE2(), 
     * DeltaEta2() and DeltaKappa2().
     * See @cite Degrassi:1996mg and @cite Degrassi:1999jd.
     * @param[in] x a real variable
     * @return @f$\phi(x)@f$
     *
     * @attention This function is valid for @f$x\geq 0@f$.
     */
    double phi(const double x) const;

    /**
     * @brief The auxiliary function @f$f(x)@f$.
     * @details This function is used in Vadd() and Vfi().
     * See @cite Degrassi:1990ec.
     * @param[in] x a real variable 
     * @return @f$f(x)@f$
     *
     * @attention This function is valid for @f$x>0@f$.
     */
    gslpp::complex FV(const double x) const;

    /**
     * @brief The auxiliary function @f$g(x)@f$.
     * @details This function is used in Vadd() and Vfi().
     * See @cite Degrassi:1990ec.
     * @param[in] x a real variable
     * @return @f$g(x)@f$
     *
     * @attention This function is valid for @f$0<x<4@f$.
     */
    gslpp::complex GV(const double x) const;

};

#endif	/* EWSMTWOLOOPEW_H */

