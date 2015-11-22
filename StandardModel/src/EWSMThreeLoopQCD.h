/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSMTHREELOOPQCD_H
#define	EWSMTHREELOOPQCD_H

#include "EWSMcache.h"

/**
 * @class EWSMThreeLoopQCD
 * @ingroup StandardModel
 * @brief A class for @f$O(\alpha\alpha_s^2)@f$ three-loop corrections to the
 * %EW precision observables.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class handles three-loop %QCD contributions of
 * @f$O(\alpha\alpha_s^2)@f$ to the following quantities, which are relevant to
 * the %EW precision observables:
 *
 * @li @f$\Delta\alpha_{\mathrm{lept}}(M_Z^2)@f$&nbsp;&nbsp; (with DeltaAlpha_l()),
 * @li @f$\Delta\alpha_{\mathrm{top}}(M_Z^2)@f$&nbsp;&nbsp; (with DeltaAlpha_t()),
 * @li @f$\Delta\rho@f$&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; (with DeltaRho()),
 * @li @f$\Delta r_{\mathrm{rem}}@f$&nbsp;&nbsp; (with DeltaR_rem()),
 * @li @f$\delta\rho_{\mathrm{rem}}^{f}@f$&nbsp;&nbsp; (with deltaRho_rem_l() and deltaRho_rem_q()),
 * @li @f$\delta\kappa_{\mathrm{rem}}^{f}@f$&nbsp;&nbsp; (with deltaKappa_rem_l() and deltaKappa_rem_q()).
 *
 * See also the description of EWSM class for their definitions.
 */
class EWSMThreeLoopQCD {
public:

    /**
     * @brief Constructor. 
     * @param[in] cache_i a reference to an object of type EWSMcache
     */
    EWSMThreeLoopQCD(const EWSMcache& cache_i);


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief Leptonic contribution of @f$O(\alpha\alpha_s^2)@f$
     * to the electromagnetic coupling @f$\alpha@f$,
     * denoted as @f$\Delta\alpha_{\mathrm{lept}}^{\alpha\alpha_s^2}(s)@f$.
     * @details This contribution vanishes at @f$O(\alpha\alpha_s^2)@f$. 
     * @param[in] s invariant mass squared
     * @return @f$\Delta\alpha_{\mathrm{lept}}^{\alpha\alpha_s^2}(s)=0@f$
     */
    double DeltaAlpha_l(const double s) const;

    /**
     * @brief Top-quark contribution of @f$O(\alpha\alpha_s^2)@f$
     * to the electromagnetic coupling @f$\alpha@f$,
     * denoted as @f$\Delta\alpha_{\mathrm{top}}^{\alpha\alpha_s^2}(s)@f$.
     * @details A simple numerical formula presented in @cite Kuhn:1998ze has
     * been employed. See also @cite Chetyrkin:1995ii, @cite Chetyrkin:1996cf
     * and @cite Chetyrkin:1997mb.
     * @param[in] s invariant mass squared
     * @return @f$\Delta\alpha_{\mathrm{top}}^{\alpha\alpha_s^2}(s)@f$
     */
    double DeltaAlpha_t(const double s) const;

    /**
     * @brief Leading three-loop %QCD contribution of @f$O(\alpha\alpha_s^2)@f$
     * to @f$\Delta\rho@f$, denoted as @f$\Delta\rho^{\alpha\alpha_s^2}@f$.
     * @details The formula used here is given by
     * @f[
     * \Delta\rho^{\alpha\alpha_s^2} 
     * = 3\,X_t^\alpha\biggl(\frac{\alpha_s(m_t^2)}{\pi}\biggr)^2 \delta^{\mathrm{QCD}}_3,
     * @f]
     * where @f$X_t^\alpha = \alpha\, m_t^2/(16\pi s_W^2 M_W^2)@f$, and
     * @f$\delta^{\mathrm{QCD}}_3@f$ is computed via deltaQCD_3().
     * See @cite Avdeev:1994db, @cite Chetyrkin:1995ix, @cite Chetyrkin:1995js
     * and Chapter 8 of @cite Bardin:1999ak.
     * This quantity contributes to @f$\Delta r@f$ and the @f$Zf\bar{f}@f$
     * effective couplings @f$\rho_Z^f@f$ and @f$\kappa_Z^f@f$.
     * See also the description of EWSM class.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\Delta\rho^{\alpha\alpha_s^2}@f$
     */
    double DeltaRho(const double Mw_i) const;

    /**
     * @brief Remainder contribution of @f$O(\alpha\alpha_s^2)@f$ to @f$\Delta r@f$,
     * denoted as @f$\Delta r_{\mathrm{rem}}^{\alpha\alpha_s^2}@f$.
     * @details The three-loop remainder contribution of @f$O(\alpha\alpha_s^2)@f$
     * is obtained from the @f$O(\alpha)@f$ remainder contribution @cite Halzen:1990je :
     * @f[
     * \Delta r_{\mathrm{rem},ud}^{\alpha}
     * \biggl[1+\frac{\alpha_s(M_Z^2)}{\pi}
     * + 1.4097\, \biggl(\frac{\alpha_s(M_Z^2)}{\pi}\biggr)^2\biggr]
     * = \Delta r_{\mathrm{rem},ud}^{\alpha} 
     *   + \Delta r_{\mathrm{rem}}^{\alpha\alpha_s}
     *   + \Delta r_{\mathrm{rem}}^{\alpha\alpha_s^2},
     * @f]
     * where @f$\Delta r_{\mathrm{rem},ud}^{\alpha}@f$ is the one-loop light-quark
     * contribution to @f$\Delta r_{\mathrm{rem}}^{\alpha}@f$ and given by
     * @f[
     * \Delta r_{\mathrm{rem},ud}^{\alpha}
     * = - \frac{\alpha}{\pi} \frac{c_W^2 - s_W^2}{4s_W^4}\,\ln c_W^2,
     * @f]
     * and the %QCD corrections are associated with the @f$R@f$ ratio
     * (see, e.g., @cite Baikov:2008jh):
     * @f[
     * R = \frac{11}{3}\biggl[ 1 + \frac{\alpha_s(M_Z^2)}{\pi}
     * + 1.4097\, \biggl(\frac{\alpha_s(M_Z^2)}{\pi}\biggr)^2 + \cdots \biggr].
     * @f]
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\Delta r_{\mathrm{rem}}^{\alpha\alpha_s^2}@f$ 
     */
    double DeltaR_rem(const double Mw_i) const;

    /**
     * @brief Remainder contribution of @f$O(\alpha\alpha_s^2)@f$ 
     * to the effective couplings @f$\rho_Z^f@f$,
     * denoted as @f$\delta\rho_{\mathrm{rem}}^{f,\, \alpha\alpha_s^2}@f$.
     * @details This contribution is not implemented, since it is tiny and negligible.
     * @param[in] f a lepton or quark
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\delta\rho_{\mathrm{rem}}^{f,\, \alpha\alpha_s^2}=0@f$
     */
    gslpp::complex deltaRho_rem_f(const Particle f, const double Mw_i) const;

    /**
     * @brief Remainder contribution of @f$O(\alpha\alpha_s^2)@f$ 
     * to the effective couplings @f$\kappa_Z^f@f$,
     * denoted as @f$\delta\kappa_{\mathrm{rem}}^{f,\, \alpha\alpha_s^2}@f$.
     * @details The formula used here is given by
     * @f[
     * \delta\kappa_{\mathrm{rem}}^{f,\alpha\alpha_s^2}
     * = - 3\,X_t^\alpha \frac{c_W^2}{s_W^2}
     * \biggl(\frac{\alpha_s(m_t^2)}{\pi}\biggr)^2
     * \bigl( \delta^{\mathrm{QCD}}_3
     * + \mathrm{Re}\,[\delta^{\mathrm{QCD}}_{\kappa,\,3}]\bigr),
     * @f]
     * where @f$\delta^{\mathrm{QCD}}_3@f$ and @f$\delta^{\mathrm{QCD}}_3@f$
     * are computed via deltaQCD_3() and deltaQCD_kappa3(), respectively.
     * See @cite Avdeev:1994db, @cite Chetyrkin:1995ix and @cite Chetyrkin:1995js.
     * @param[in] f a lepton or quark
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\delta\kappa_{\mathrm{rem}}^{f,\, \alpha\alpha_s^2}@f$
     */
    gslpp::complex deltaKappa_rem_f(const Particle f, const double Mw_i) const;


    ////////////////////////////////////////////////////////////////////////        

private:
    const EWSMcache& cache; ///< A reference to an object of type EWSMcache.


    ////////////////////////////////////////////////////////////////////////        

    /**
     * @brief The function @f$\delta^{\mathrm{QCD}}_3@f$. 
     * @details This function is associated with the leading three-loop %QCD
     * contribution of @f$O(\alpha\alpha_s^2(m_t^2/M_Z^2+1+M_Z^2/m_t^2))@f$
     * to @f$\Delta\rho@f$, as explained in the description of DeltaRho(). 
     * See @cite Avdeev:1994db, @cite Chetyrkin:1995ix, @cite Chetyrkin:1995js
     * and Chapter 8 of @cite Bardin:1999ak.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\delta^{\mathrm{QCD}}_3@f$
     */
    double deltaQCD_3(const double Mw_i) const;

    /**
     * @brief The function @f$\delta^{\mathrm{QCD}}_{\kappa, 3}@f$.
     * @details The sum
     * @f$\delta^{\mathrm{QCD}}_3 + \delta^{\mathrm{QCD}}_{\kappa, 3}@f$
     * corresponds to the
     * @f$O(\alpha\alpha_s^2(m_t^2/M_Z^2+1+M_Z^2/m_t^2))@f$ contribution
     * to @f$\delta\kappa_{\mathrm{rem}}^{f}@f$.
     * See the arXiv version of @cite Chetyrkin:1995js 
     * (and also @cite Avdeev:1994db, @cite Chetyrkin:1995ix
     * and Chapter 8 of @cite Bardin:1999ak).
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\delta^{\mathrm{QCD}}_{\kappa,3}@f$
     */
    gslpp::complex deltaQCD_kappa3(const double Mw_i) const;


};

#endif	/* EWSMTHREELOOPQCD_H */

