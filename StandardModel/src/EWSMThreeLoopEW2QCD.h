/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSMTHREELOOPEW2QCD_H
#define	EWSMTHREELOOPEW2QCD_H

#include "EWSMcache.h"

/**
 * @class EWSMThreeLoopEW2QCD
 * @ingroup StandardModel
 * @brief A class for @f$O(\alpha^2\alpha_s)@f$ three-loop corrections to the
 * %EW precision observables.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class handles three-loop mixed %EW-%QCD contributions of
 * @f$O(\alpha^2\alpha_s)@f$ to the following quantities, which are relevant to 
 * the %EW precision observables:
 *
 * @li @f$\Delta\alpha_{\mathrm{lept}}(M_Z^2)@f$&nbsp;&nbsp; (with DeltaAlpha_l()),
 * @li @f$\Delta\alpha_{\mathrm{top}}(M_Z^2)@f$&nbsp;&nbsp; (with DeltaAlpha_t()),
 * @li @f$\Delta\rho@f$&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; (with DeltaRho()),
 * @li @f$\Delta r_{\mathrm{rem}}@f$&nbsp;&nbsp; (with DeltaR_rem()),
 * @li @f$\delta\rho_{\mathrm{rem}}^{f}@f$&nbsp;&nbsp; (with deltaRho_rem_l() and deltaRho_rem_q()),
 * @li @f$\delta\kappa_{\mathrm{rem}}^{f}@f$&nbsp;&nbsp; (with deltaKappa_rem_l() and deltaKappa_rem_q()),
 *
 * where only @f$\Delta\rho@f$ is non-zero in the current class.
 * See also the description of EWSM class for their definitions. 
 */
class EWSMThreeLoopEW2QCD {
public:

    /**
     * @brief Constructor.
     * @param[in] cache_i a reference to an object of type EWSMcache
     */
    EWSMThreeLoopEW2QCD(const EWSMcache& cache_i);


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief Leptonic contribution of @f$O(\alpha^2\alpha_s)@f$
     * to the electromagnetic coupling @f$\alpha@f$,
     * denoted as @f$\Delta\alpha_{\mathrm{lept}}^{\alpha^2\alpha_s}(s)@f$.
     * @details This contribution vanishes at @f$O(\alpha^2\alpha_s)@f$.
     * @param[in] s invariant mass squared
     * @return @f$\Delta\alpha_{\mathrm{lept}}^{\alpha^2\alpha_s}(s)=0@f$
     */
    double DeltaAlpha_l(const double s) const;

    /**
     * @brief Top-quark contribution of @f$O(\alpha^2\alpha_s)@f$
     * to the electromagnetic coupling @f$\alpha@f$,
     * denoted as @f$\Delta\alpha_{\mathrm{top}}^{\alpha^2\alpha_s}(s)@f$.
     * @details This contribution is not implemented, since it is tiny and negligible.
     * @param[in] s invariant mass squared
     * @return @f$\Delta\alpha_{\mathrm{top}}^{\alpha^2\alpha_s}(s)=0@f$
     */
    double DeltaAlpha_t(const double s) const;

    /**
     * @brief Leading three-loop contribution of @f$O(\alpha^2\alpha_s)@f$
     * to @f$\Delta\rho@f$, denoted as @f$\Delta\rho^{\alpha^2\alpha_s}@f$. 
     * @details This function represents the leading three-loop mixed %EW-%QCD
     * contribution of @f$O(\alpha^2\alpha_s m_t^4/M_Z^4)@f$ to @f$\Delta\rho@f$.
     * Expressions are available for @f$m_h=0@f$ in @cite vanderBij:2000cg
     * and for @f$m_h\approx m_t@f$ and @f$m_h\gg m_t@f$ in @cite Faisst:2003px.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\Delta\rho^{\alpha^2\alpha_s}@f$
     */
    double DeltaRho(const double Mw_i) const;

    /**
     * @brief Remainder contribution of @f$O(\alpha^2\alpha_s)@f$ to @f$\Delta r@f$,
     * denoted as @f$\Delta r_{\mathrm{rem}}^{\alpha^2\alpha_s}@f$.
     * @details This contribution is not implemented, since it is tiny and negligible.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\Delta r_{\mathrm{rem}}^{\alpha^2\alpha_s}=0@f$
     */
    double DeltaR_rem(const double Mw_i) const;

    /**
     * @brief Remainder contribution of @f$O(\alpha^2\alpha_s)@f$
     * to the effective couplings @f$\rho_Z^f@f$,
     * denoted as @f$\delta\rho_{\mathrm{rem}}^{f,\, \alpha^2\alpha_s}@f$.
     * @details This contribution is not implemented, since it is tiny and negligible.
     * @param[in] f a lepton or quark
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\delta\rho_{\mathrm{rem}}^{f,\, \alpha^2\alpha_s}=0@f$
     */
    gslpp::complex deltaRho_rem_f(const Particle f, const double Mw_i) const;

    /**
     * @brief Remainder contribution of @f$O(\alpha^2\alpha_s)@f$
     * to the effective couplings @f$\kappa_Z^f@f$,
     * denoted as @f$\delta\kappa_{\mathrm{rem}}^{f,\, \alpha^2\alpha_s}@f$.
     * @details This contribution is not implemented, since it is tiny and negligible.
     * @param[in] f a lepton or quark
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\delta\kappa_{\mathrm{rem}}^{f,\, \alpha^2\alpha_s}=0@f$
     */
    gslpp::complex deltaKappa_rem_f(const Particle f, const double Mw_i) const;


    ////////////////////////////////////////////////////////////////////////        

private:
    const EWSMcache& cache; ///< A reference to an object of type EWSMcache.


};

#endif	/* EWSMTHREELOOPEW2QCD_H */

