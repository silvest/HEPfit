/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSMTWOLOOPQCD_H
#define	EWSMTWOLOOPQCD_H

#include "EWSMcache.h"

/**
 * @class EWSMTwoLoopQCD
 * @ingroup StandardModel
 * @brief A class for @f$O(\alpha\alpha_s)@f$ two-loop corrections to the %EW
 * precision observables.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class handles two-loop %QCD contributions of
 * @f$O(\alpha\alpha_s)@f$ to the following quantities, which are relevant to
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
 * The above quantities are computed with the auxiliary functions: 
 *
 * @li @f$\delta_2^{\mathrm{QCD}}@f$&nbsp;&nbsp; (with deltaQCD_2()),
 * @li @f$\Delta r^{ud}@f$&nbsp;&nbsp; (with DeltaR_ud()),
 * @li @f$\Delta r^{tb}@f$&nbsp;&nbsp; (with DeltaR_tb()),
 * @li @f$\Delta \rho^{ud}@f$&nbsp;&nbsp; (with DeltaRho_ud()),
 * @li @f$\Delta \rho^{tb}@f$&nbsp;&nbsp; (with DeltaRho_tb()),
 * @li @f$\Delta \kappa^{ud}@f$&nbsp;&nbsp; (with DeltaKappa_ud()),
 * @li @f$\Delta \kappa^{tb}@f$&nbsp;&nbsp; (with DeltaKappa_tb()),
 *
 * and
 *
 * @li @f$F_1(x)@f$&nbsp;&nbsp; (with F1()),
 * @li @f$V_1(r)@f$&nbsp;&nbsp; (with V1()),
 * @li @f$V'_1(r)@f$&nbsp;&nbsp; (with V1prime()),
 * @li @f$A_1(r)@f$&nbsp;&nbsp; (with A1()),
 * @li @f$A'_1(r)@f$&nbsp;&nbsp; (with A1prime()).
 *
 * The @f$O(\alpha\alpha_s)@f$ two-loop %QCD contributions to the vacuum
 * polarization amplitudes of the gauge bosons were calculated in
 * @cite Djouadi:1987gn, @cite Djouadi:1987di, @cite Kniehl:1989yc, 
 * @cite Halzen:1990je, @cite Kniehl:1991gu, @cite Kniehl:1992dx and @cite Djouadi:1993ss.
 */
class EWSMTwoLoopQCD {
public:

    /**
     * @brief Constructor.
     * @param[in] cache_i a reference to an object of type EWSMcache
     */
    EWSMTwoLoopQCD(const EWSMcache& cache_i);


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief Leptonic contribution of @f$O(\alpha\alpha_s)@f$
     * to the electromagnetic coupling @f$\alpha@f$,
     * denoted as @f$\Delta\alpha_{\mathrm{lept}}^{\alpha\alpha_s}(s)@f$.
     * @details This contribution vanishes at @f$O(\alpha\alpha_s)@f$.
     * @param[in] s invariant mass squared
     * @return @f$\Delta\alpha_{\mathrm{lept}}^{\alpha\alpha_s}(s)=0@f$
     */
    double DeltaAlpha_l(const double s) const;

    /**
     * @brief Top-quark contribution of @f$O(\alpha\alpha_s)@f$
     * to the electromagnetic coupling @f$\alpha@f$,
     * denoted as @f$\Delta\alpha_{\mathrm{top}}^{\alpha\alpha_s}(s)@f$.
     * @details A simple numerical formula presented in @cite Kuhn:1998ze has
     * been employed.
     * @param[in] s invariant mass squared
     * @return @f$\Delta\alpha_{\mathrm{top}}^{\alpha\alpha_s}(s)@f$
     */
    double DeltaAlpha_t(const double s) const;

    /**
     * @brief Leading two-loop %QCD contribution of @f$O(\alpha\alpha_s)@f$
     * to @f$\Delta\rho@f$, denoted as @f$\Delta\rho^{\alpha\alpha_s}@f$.
     * @details The formula used here is given by
     * @f[
     * \Delta\rho^{\alpha\alpha_s}
     * = 3\,X_t^\alpha\frac{\alpha_s(m_t^2)}{\pi} \delta^{\mathrm{QCD}}_2,
     * @f]
     * where @f$X_t^\alpha = \alpha\, m_t^2/(16\pi s_W^2 M_W^2)@f$, and
     * @f$\delta^{\mathrm{QCD}}_2@f$ is computed via deltaQCD_2().
     * See, e.g., Chapter 8 of @cite Bardin:1999ak.
     * This quantity contributes to @f$\Delta r@f$ and the @f$Zf\bar{f}@f$
     * effective couplings @f$\rho_Z^f@f$ and @f$\kappa_Z^f@f$.
     * See also the description of EWSM class.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\Delta\rho^{\alpha\alpha_s}@f$
     */
    double DeltaRho(const double Mw_i) const;

    /**
     * @brief Remainder contribution of @f$O(\alpha\alpha_s)@f$ to @f$\Delta r@f$,
     * denoted as @f$\Delta r_{\mathrm{rem}}^{\alpha\alpha_s^2}@f$. 
     * @details The @f$O(\alpha\alpha_s)@f$ remainder contribution to
     * @f$\Delta r@f$ is decomposed as
     * @f[
     * \Delta r_{\mathrm{rem}}^{\alpha\alpha_s}
     *  = 2 \Delta r^{ud} + \Delta r^{tb}
     *    + \frac{c_W^2}{s_W^2} \Delta\rho^{\alpha\alpha_s},
     * @f]
     * where @f$\Delta r^{ud}@f$ and @f$\Delta r^{tb}@f$ are associated with 
     * corrections to the self-energies of the gauge bosons with loops of
     * a light-quark doublet and with those of the @f$t@f$-@f$b@f$ doublet, respectively,
     * and @f$\Delta\rho^{\alpha\alpha_s}@f$ is the leading contribution.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\Delta r_{\mathrm{rem}}^{\alpha\alpha_s^2}@f$ 
     */
    double DeltaR_rem(const double Mw_i) const;

    /**
     * @brief Remainder contribution of @f$O(\alpha\alpha_s)@f$
     * to the effective couplings @f$\rho_Z^f@f$,
     * denoted as @f$\delta\rho_{\mathrm{rem}}^{f,\, \alpha\alpha_s}@f$.
     * @details The @f$O(\alpha\alpha_s)@f$ remainder contribution to
     * @f$\rho_{Z}^{f}@f$ is decomposed as
     * @f[
     * \delta\rho_{\mathrm{rem}}^{f,\,\alpha\alpha_s}
     *  = 2 \Delta\rho^{ud} + \Delta\rho^{tb} - \Delta\rho^{\alpha\alpha_s},
     * @f]
     * where @f$\Delta\rho^{ud}@f$ and @f$\Delta\rho^{tb}@f$ are associated with
     * corrections to the self-energies of the gauge bosons with loops of
     * a light-quark doublet and with those of the @f$t@f$-@f$b@f$ doublet, respectively,
     * and @f$\Delta\rho^{\alpha\alpha_s}@f$ is the leading contribution
     * of @f$O(\alpha\alpha_s)@f$ to @f$\rho_{Z}^{f}@f$.
     * @param[in] f a lepton or quark
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\delta\rho_{\mathrm{rem}}^{f,\, \alpha\alpha_s}@f$
     */
    gslpp::complex deltaRho_rem_f(const Particle f, const double Mw_i) const;

    /**
     * @brief Remainder contribution of @f$O(\alpha\alpha_s)@f$
     * to the effective couplings @f$\kappa_Z^f@f$,
     * denoted as @f$\delta\kappa_{\mathrm{rem}}^{f,\, \alpha\alpha_s}@f$.
     * @details The @f$O(\alpha\alpha_s)@f$ remainder contribution to
     * @f$\kappa_{Z}^{f}@f$ is decomposed as
     * @f[
     * \delta\kappa_{\mathrm{rem}}^{f,\,\alpha\alpha_s}
     *  = 2 \Delta\kappa^{ud} + \Delta\kappa^{tb}
     * - \frac{c_W^2}{s_W^2}\Delta\rho^{\alpha\alpha_s},
     * @f]
     * where @f$\Delta\kappa^{ud}@f$ and @f$\Delta\kappa^{tb}@f$ are associated with
     * corrections to the self-energies of the gauge bosons with loops of
     * a light-quark doublet and with those of the @f$t@f$-@f$b@f$ doublet, respectively,
     * and @f$(c_W^2/s_W^2)\Delta\rho^{\alpha\alpha_s}@f$ is the leading contribution
     * of @f$O(\alpha\alpha_s)@f$ to @f$\kappa_{Z}^{f}@f$.
     * @param[in] f a lepton or quark
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\delta\kappa_{\mathrm{rem}}^{f,\, \alpha\alpha_s}@f$
     */
    gslpp::complex deltaKappa_rem_f(const Particle f, const double Mw_i) const;


    ////////////////////////////////////////////////////////////////////////        

    /**
     * @brief The function @f$\delta^{\mathrm{QCD}}_2@f$.
     * @details This function is associated with the leading two-loop %QCD
     * contribution of @f$O(\alpha\alpha_s m_t^2/M_Z^2)@f$ to @f$\Delta\rho@f$,
     * as explained in the description of DeltaRho().
     * See @cite Kniehl:1989yc, @cite Halzen:1990je 
     * and Chapter 8 of @cite Bardin:1999ak.
     * @return @f$\delta^{\mathrm{QCD}}_2@f$
     */
    double deltaQCD_2() const;

    /**
     * @brief The function @f$F_1(x)@f$.
     * @details The expression for @f$F_1(x)@f$ can be found in @cite Kniehl:1989yc
     * and @cite Halzen:1990je. See also Chapter 8 of @cite Bardin:1999ak.
     * @param[in] x the ratio @f$x=s/m_t^2@f$
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$F_1(x)@f$
     *
     * @attention This function is valid for @f$0 \leq x < 1@f$.
     */
    double F1(const double x, const double Mw_i) const;

    /**
     * @brief The function @f$V_1(r)@f$. 
     * @details The expression for @f$V_1(r)@f$ can be found in @cite Kniehl:1989yc
     * and @cite Halzen:1990je. See also Chapter 8 of @cite Bardin:1999ak.
     * @param[in] r the ratio @f$r=s/(4m_t^2)@f$
     * @return @f$V_1(r)@f$
     *
     * @attention This function is valid for @f$0 \leq r < 1@f$. 
     */
    double V1(const double r) const;

    /**
     * @brief The function @f$A_1(r)@f$.
     * @details The expression for @f$A_1(r)@f$ can be found in @cite Kniehl:1989yc
     * and @cite Halzen:1990je. See also Chapter 8 of @cite Bardin:1999ak.
     * @param[in] r the ratio @f$r=s/(4m_t^2)@f$
     * @return @f$A_1(r)@f$
     *
     * @attention This function is valid for @f$0 \leq r < 1@f$.
     */
    double A1(const double r) const;

    /**
     * @brief The derivative of the function @f$V_1(r)@f$.
     * @details The expression for @f$V'_1(r)@f$ has been derived from @f$V_1(r)@f$
     * in @cite Kniehl:1989yc and @cite Halzen:1990je.
     * @param[in] r the ratio @f$r=s/(4m_t^2)@f$
     * @return @f$V'_1(r)@f$
     *
     * @attention This function is valid for @f$0 \leq r < 1@f$.
     */
    double V1prime(const double r) const;

    /**
     * @brief The derivative of the function @f$A_1(r)@f$.
     * @details The expression for @f$A'_1(r)@f$ has been derived from @f$A_1(r)@f$
     * in @cite Kniehl:1989yc and @cite Halzen:1990je.
     * @param[in] r the ratio @f$r=s/(4m_t^2)@f$
     * @return @f$A'_1(r)@f$
     * 
     * @attention This function is valid for @f$0 \leq r < 1@f$. 
     */
    double A1prime(const double r) const;

    /**
     * @brief Light-quark contribution to @f$\Delta r@f$, not including
     * @f$\Delta\alpha^{l+5q}(M_Z^2)@f$, denoted as @f$\Delta r^{ud}@f$. 
     * @details The quantity @f$\Delta r^{ud}@f$ is associated with
     * @f$O(\alpha\alpha_s)@f$ corrections to the self-energies of the 
     * gauge bosons with loops of the light-quark doublets.
     * The expression of @f$\Delta r^{ud}@f$ is given by
     * @f[
     * \Delta r^{ud} 
     * = - \frac{\alpha\alpha_s(M_Z^2)}{\pi^2}\frac{c_W^2 - s_W^2}{4s_W^4}\,\ln c_W^2.
     * @f]
     * See @cite Kniehl:1989yc, @cite Halzen:1990je, @cite Kniehl:1991gu and
     * Chapter 8 of @cite Bardin:1999ak. 
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\Delta r^{ud}@f$
     */
    double DeltaR_ud(const double Mw_i) const;

    /**
     * @brief Heavy-quark contribution to @f$\Delta r@f$,
     * denoted as @f$\Delta r^{tb}@f$.
     * @details The quantity @f$\Delta r^{tb}@f$ is associated with
     * @f$O(\alpha\alpha_s)@f$ corrections to the self-energies of the 
     * gauge bosons with loops of the @f$t@f$-@f$b@f$ doublet.
     * The expression of @f$\Delta r^{tb}@f$ is given by
     * @f[
     * \Delta r^{tb}
     * = \frac{\alpha\alpha_s(M_t^2)}{\pi^2}
     * \Bigg\{ Q_t^2 V_1'(0)
     *   + \frac{c_W^2}{s_W^4}\frac{w_t}{4}\left[ \zeta(2) + \frac{1}{2} \right]
     *   - \frac{z_t}{s_W^4}
     *     \Big[ v_t^2 V_1(r^Z_{4t}) + a_t^2\left[A_1(r^Z_{4t}) - A_1(0)\right] \Big]
     *   + \frac{c_W^2 - s_W^2}{s_W^4} w_t \Big[ F_1(x^W_{t}) -  F_1(0) \Big]
     *   - \frac{v_ba_b}{2 s_W^4} \ln z_t
     * \Bigg\},
     * @f]
     * where the definitions of the symbols can be read from the codes below. 
     * See @cite Kniehl:1989yc, @cite Halzen:1990je, @cite Kniehl:1991gu and
     * Chapter 8 of @cite Bardin:1999ak.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\Delta r^{tb}@f$
     */
    double DeltaR_tb(const double Mw_i) const;

    /**
     * @brief Light-quark contribution to @f$\rho_Z^f@f$,
     * denoted as @f$\Delta\rho^{ud}@f$. 
     * @details The quantity @f$\Delta\rho^{ud}@f$ is associated with
     * @f$O(\alpha\alpha_s)@f$ corrections to the self-energies of the 
     * gauge bosons with loops of the light-quark doublets.
     * The expression of @f$\Delta\rho^{ud}@f$ is given by
     * @f[
     * \Delta\rho^{ud}
     * = \frac{\alpha\alpha_s(M_Z^2)}{\pi^2} \frac{1}{4s_W^2 c_W^2}
     *   \left( v_u^2 + v_d^2 + a_u^2 + a_d^2 \right),
     * @f]
     * where the definitions of the symbols can be read from the codes below. 
     * See @cite Kniehl:1989yc, @cite Halzen:1990je, @cite Kniehl:1991gu and
     * Chapter 8 of @cite Bardin:1999ak.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\Delta\rho^{ud}@f$
     */
    double DeltaRho_ud(const double Mw_i) const;

    /**
     * @brief Heavy-quark contribution to @f$\rho_Z^f@f$,
     * denoted as @f$\Delta\rho^{tb}@f$.
     * @details The quantity @f$\Delta\rho^{tb}@f$ is associated with
     * @f$O(\alpha\alpha_s)@f$ corrections to the self-energies of the 
     * gauge bosons with loops of the @f$t@f$-@f$b@f$ doublet.
     * The expression of @f$\Delta\rho^{tb}@f$ is given by
     * @f[
     * \Delta\rho^{tb}
     * = \frac{\alpha\alpha_s(M_t^2)}{\pi^2} \frac{1}{4s_W^2 c_W^2}
     *   \Big\{ - \left[ v_t^2 V_1'(r^Z_{4t}) +  a_t^2 A_1'(r^Z_{4t}) \right]
     *          + 4 z_t \left[  v_t^2 V_1(r^Z_{4t}) + a_t^2 A_1(r^Z_{4t}) \right]
     *          + v_b^2 + a_b^2 - 4 z_t\, F_1(0) \Big\},
     * @f]
     * where the definitions of the symbols can be read from the codes below. 
     * See @cite Kniehl:1989yc, @cite Halzen:1990je, @cite Kniehl:1991gu and
     * Chapter 8 of @cite Bardin:1999ak.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\Delta\rho^{tb}@f$
     */
    double DeltaRho_tb(const double Mw_i) const;

    /**
     * @brief Light-quark contribution to @f$\kappa_Z^f@f$,
     * denoted as @f$\Delta\kappa^{ud}@f$.
     * @details The quantity @f$\Delta\kappa^{ud}@f$ is associated with
     * @f$O(\alpha\alpha_s)@f$ corrections to the self-energies of the 
     * gauge bosons with loops of the light-quark doublets.
     * The expression of @f$\Delta\kappa^{ud}@f$ is given by
     * @f[
     * \Delta\kappa^{ud}
     * = \frac{\alpha\alpha_s(M_Z^2)}{\pi^2} \frac{c_W^2}{4 s_W^4}\ln c_W^2
     *   + i \frac{\alpha\alpha_s(M_Z^2)}{4\pi s_W^2}
     *     \left( 1 - \frac{20}{9}s_W^2 \right),
     * @f]
     * where the definitions of the symbols can be read from the codes below. 
     * See @cite Kniehl:1989yc, @cite Halzen:1990je, @cite Kniehl:1991gu and
     * Chapter 8 of @cite Bardin:1999ak.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\Delta\kappa^{ud}@f$
     */
    gslpp::complex DeltaKappa_ud(const double Mw_i) const;

    /**
     * @brief Heavy-quark contribution to @f$\kappa_Z^f@f$,
     * denoted as @f$\Delta\kappa^{tb}@f$. 
     * @details The quantity @f$\Delta\kappa^{tb}@f$ is associated with
     * @f$O(\alpha\alpha_s)@f$ corrections to the self-energies of the 
     * gauge bosons with loops of the @f$t@f$-@f$b@f$ doublet.
     * The expression of @f$\Delta\kappa^{tb}@f$ is given by
     * @f[
     * \Delta\kappa^{tb}
     * = \frac{\alpha\alpha_s(M_t^2)}{\pi^2} \frac{1}{4s_W^4}
     *   \Bigg\{
     *     4 c_W^2 w_t \Big[ v_t^2 V_1(r^Z_{4t}) + a_t^2A_1(r^Z_{4t}) - F_1(x^W_{t}) \Big]
     *    + 4s_W^2 \left( |Q_t| - 4s_W^2 Q_t^2 \right) z_t V_1(r^Z_{4t})
     *    + \Big[ v_b^2 + a_b ^2 + s_W^2 \left( |Q_b| - 4s_W^2 Q_b^2 \right) \Big]\ln z_t
     *   \Bigg\}
     *   + i \frac{\alpha\alpha_s(M_Z^2)}{4\pi s_W^2}
     *     \left( \frac{1}{3} - \frac{4}{9}s_W^2 \right),
     * @f]
     * where the definitions of the symbols can be read from the codes below. 
     * See @cite Kniehl:1989yc, @cite Halzen:1990je, @cite Kniehl:1991gu and
     * Chapter 8 of @cite Bardin:1999ak.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\Delta\kappa^{tb}@f$
     */
    gslpp::complex DeltaKappa_tb(const double Mw_i) const;


    ////////////////////////////////////////////////////////////////////////        

private:
    const EWSMcache& cache; ///< A reference to an object of type EWSMcache.


};

#endif	/* EWSMTWOLOOPQCD_H */

