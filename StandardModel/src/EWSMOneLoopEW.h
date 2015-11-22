/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSMONELOOPEW_H
#define	EWSMONELOOPEW_H

#include "EWSMcache.h"

/**
 * @class EWSMOneLoopEW
 * @ingroup StandardModel
 * @brief A class for @f$O(\alpha)@f$ one-loop corrections to the %EW
 * precision observables.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class handles one-loop %EW contributions of
 * @f$O(\alpha)@f$ to the following quantities, which are relevant to the %EW
 * precision observables: 
 *
 * @li @f$\Delta\alpha_{\mathrm{lept}}(M_Z^2)@f$&nbsp;&nbsp; (with DeltaAlpha_l()),
 * @li @f$\Delta\alpha_{\mathrm{top}}(M_Z^2)@f$&nbsp;&nbsp; (with DeltaAlpha_t()),
 * @li @f$\Delta\rho@f$&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; (with DeltaRho()),
 * @li @f$\Delta r_{\mathrm{rem}}@f$&nbsp;&nbsp; (with DeltaR_rem()),
 * @li @f$\delta\rho_{\mathrm{rem}}^{f}@f$&nbsp;&nbsp; (with deltaRho_rem_f()),
 * @li @f$\delta\kappa_{\mathrm{rem}}^{f}@f$&nbsp;&nbsp; (with deltaKappa_rem_f()),
 * @li @f$\rho_{ij}^W@f$&nbsp;&nbsp; (with rho_GammaW()). 
 *
 * See also the description of EWSM class for their definitions.
 *
 * The above quantities are computed with the help of the auxiliary member functions
 * defined in the current class.
 * It is noted that gauge-dependent quantities in the current class are given
 * in the Unitary gauge. 
 * For example, the functions
 * SigmabarWW_bos(), SigmabarWW_bos(), etc. represent the self-energies of the
 * gauge bosons.
 * @anchor SeflEnergies
 * The following table shows the comparisons of the definitions of the 
 * self-energies in the current class to those in @cite Bardin:1999ak and
 * @cite Hollik:1988ii :
 * <table border="1" align="center">
 * <tr>
 *   <td align="center">Current class</td>
 *   <td align="center">Bardin & Passarino @cite Bardin:1999ak</td>
 *   <td align="center">Hollik @cite Hollik:1988ii</td>
 * </tr>
 * <tr>
 *   <td align="center">@f$\displaystyle\Sigma_{WW}(s)=\frac{\alpha}{4\pi s_W^2 } \overline{\Sigma}_{WW}(s)@f$</td>
 *   <td align="center">@f$\displaystyle\frac{\alpha}{4\pi s_W^2 } \Sigma_{WW}(s)@f$</td>
 *   <td align="center">@f$\displaystyle\Sigma^W(s)@f$</td>
 * </tr>
 * <tr>
 *   <td align="center">@f$\displaystyle\Sigma_{ZZ}(s)=\frac{\alpha}{4\pi s_W^2  c_W^2 } \overline{\Sigma}_{ZZ}(s)@f$</td>
 *   <td align="center">@f$\displaystyle\frac{\alpha}{4\pi s_W^2  c_W^2 } \Sigma_{ZZ}(s)@f$</td>
 *   <td align="center">@f$\displaystyle\Sigma^Z(s)@f$</td>
 * </tr>
 * <tr>
 *   <td align="center">@f$\displaystyle\Pi_{\gamma\gamma}(s)=\frac{\alpha}{4\pi}\,\overline{\Pi}_{\gamma\gamma}(s)@f$</td>
 *   <td align="center">@f$\displaystyle\frac{\alpha}{4\pi}\,\left[ -\Pi_{\gamma\gamma}(s) \right]@f$</td>
 *   <td align="center">@f$\displaystyle\Pi^\gamma(s)=\frac{\partial\Sigma^\gamma(s)}{\partial s}@f$</td>
 * </tr>
 * <tr>
 *   <td align="center">@f$\displaystyle\Pi_{Z\gamma}(s) = \frac{\alpha}{4\pi s_W c_W}\overline{\Pi}_{Z\gamma}(s)@f$</td>
 *   <td align="center">@f$\displaystyle\frac{\alpha}{4\pi s_W c_W}\left[ - \Pi_{ZA}(s)\right]@f$</td>
 *   <td align="center"></td>
 * </tr>
 * <tr>
 *   <td align="center">@f$\displaystyle\overline{\Sigma}^{\prime}_{WW}(s)
 *       = \frac{\partial\,\overline{\Sigma}_{WW}(s)}{\partial s}@f$</td>
 *   <td align="center">@f$\displaystyle- \frac{\alpha}{4\pi s_W^2 } {\cal W}(s)@f$</td>
 *   <td align="center"></td>
 * </tr>
 * <tr>
 *   <td align="center">@f$\displaystyle\overline{\Sigma}^{\prime}_{ZZ}(s)
 *       =\frac{\partial\,\overline{\Sigma}_{ZZ}(s)}{\partial s}@f$</td>
 *   <td align="center">@f$\displaystyle- \frac{\alpha}{4\pi s_W^2  c_W^2 } c_W^2  {\cal Z}(s)@f$</td>
 *   <td align="center"></td>
 * </tr>
 * </table>
 * These self-energies are decomposed into the sum of bosonic and fermionic
 * contributions:
 * \f{align*}{
 * \overline{\Sigma}_{WW}(s)
 * = \overline{\Sigma}^{\rm bos}_{WW}(s) + \overline{\Sigma}^{\rm fer}_{WW}(s)\,,
 * \\
 * \overline{\Sigma}_{ZZ}(s)
 * = \overline{\Sigma}^{\rm bos}_{ZZ}(s) + \overline{\Sigma}^{\rm fer}_{ZZ}(s)\,,
 * \\
 * \overline{\Pi}_{\gamma\gamma}(s)
 * = \overline{\Pi}_{\gamma\gamma}^{\rm bos}(s) + \overline{\Pi}_{\gamma\gamma}^{\rm fer}(s)\,,
 * \\
 * \overline{\Pi}_{Z\gamma}(s)
 * = \overline{\Pi}_{Z\gamma}^{\rm bos}(s) + \overline{\Pi}_{Z\gamma}^{\rm fer}(s)\,.
 * \f}
 *
 * Most of the formulae used in the current class have been copied or derived 
 * from from those in Bardin and Passarino's book @cite Bardin:1999ak as well as
 * from @cite Sirlin:1980nh, @cite Marciano:1980pb, @cite Bardin:1981sv,
 * @cite Akhundov:1985fc, @cite Bardin:1986fi and @cite Bardin:1989di.
 */
class EWSMOneLoopEW {
public:

    /**
     * @brief Constructor.
     * @param[in] cache_i a reference to an object of type EWSMcache
     */
    EWSMOneLoopEW(const EWSMcache& cache_i);


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief Leptonic contribution of @f$O(\alpha)@f$
     * to the electromagnetic coupling @f$\alpha(s)@f$,
     * denoted as @f$\Delta\alpha_{\mathrm{lept}}^{\alpha}@f$.
     * @details This function uses the function PiGammaGamma_fer_l() for the
     * fermionic contribution to the photon self-energy.
     * @param[in] s invariant mass squared
     * @return @f$\Delta\alpha_{\mathrm{lept}}^{\alpha}(s)@f$
     */
    double DeltaAlpha_l(const double s) const;

    /**
     * @brief Top-quark contribution of @f$O(\alpha)@f$
     * to the electromagnetic coupling @f$\alpha(s)@f$,
     * denoted as @f$\Delta\alpha_{\mathrm{top}}^{\alpha}@f$.
     * @details A simple numerical formula presented in @cite Kuhn:1998ze has
     * been employed.
     * @param[in] s invariant mass squared
     * @return @f$\Delta\alpha_{\mathrm{top}}^{\alpha}(s)@f$
     */
    double DeltaAlpha_t(const double s) const;

    /**
     * @brief Leading one-loop contribution of @f$O(\alpha)@f$
     * to @f$\Delta\rho@f$, denoted as @f$\Delta\rho^{\alpha}@f$.
     * @details The leading one-loop contribution is written in terms of the 
     * quantity @f$\Delta\overline{\rho}@f$:
     * @f[
     * \Delta\rho^{\alpha} = - \frac{\alpha}{4\pi s_W^2}\Delta\overline{\rho}\big|_{\mu=M_Z},
     * @f]
     * where @f$\Delta\overline{\rho}@f$ is renormalized at the scale @f$M_Z@f$, 
     * and computed with the function DeltaRhobar().
     * @param[in] Mw_i the @f$W@f$-boson mass @f$M_W@f$
     * @return @f$\Delta\rho^{\alpha}@f$
     */
    double DeltaRho(const double Mw_i) const;

    /**
     * @brief Remainder contribution of @f$O(\alpha)@f$ to @f$\Delta r@f$,
     * denoted as @f$\Delta r_{\rm rem}^{\alpha}@f$.
     * @details The @f$O(\alpha)@f$ remainder contribution
     * @f$\Delta r_{\mathrm{rem}}^{\alpha}@f$ is given by
     * @f[
     * \Delta r_{\rm rem}^\alpha
     * =
     * \frac{\alpha}{4\pi s_W^2}
     * \bigg[
     * - \frac{2}{3} s_W^2
     * + s_W^2\overline{\Pi}_{\gamma\gamma}^{t}(0)\big|_{\mu=M_Z}
     * + s_W^2{\rm Re}\Big[\overline{\Pi}_{\gamma\gamma}^{\ell+5q}(M_Z^2)\big|_{\mu=M_Z}\Big]
     * + \Delta\overline{\rho}_W\big|_{\mu=M_Z}
     * + \left( 4 - \frac{25}{4}c_W^2 + \frac{3}{4} c_W^4 + \frac{9 c_W^2}{4 s_W^2} \right)\ln c_W^2
     * + \frac{11}{2} - \frac{5}{8} c_W^2 (1 + c_W^2) \bigg].
     * @f]
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\Delta r_{\mathrm{rem}}^{\alpha}@f$
     */
    double DeltaR_rem(const double Mw_i) const;

    /**
     * @brief @f$\Delta \bar{r}_{\rm rem}^{\alpha}@f$.
     * @details The quantity @f$\Delta \bar{r}_{\rm rem}^{\alpha}@f$, which
     * contributes to @f$\rho_Z^f@f$ and @f$\kappa_Z^f@f$, is given by
     * @f[
     * \Delta \bar{r}_{\rm rem}^{\alpha}
     * = \bigg[
     * \Delta r_{\rm rem}^\alpha
     * - \frac{\alpha}{4\pi} \overline{\Pi}_{\gamma\gamma}^{t}(0)\big|_{\mu=M_Z}
     * - \frac{\alpha}{4\pi}
     * {\rm Re}\left[ \overline{\Pi}_{\gamma\gamma}^{\ell+5q}(M_Z^2)\big|_{\mu=M_Z} \right] \bigg],
     * @f]
     * where the definition of @f$\Delta r_{\rm rem}^\alpha@f$ is given in
     * DeltaR_rem().
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\Delta \bar{r}_{\rm rem}^{\alpha}@f$
     */
    double DeltaRbar_rem(const double Mw_i) const;

    /**
     * @brief Remainder contribution of @f$O(\alpha)@f$ to the effective
     * couplings @f$\rho_Z^f@f$,
     * denoted as @f$\delta\rho_{\rm rem}^{f,\alpha}@f$.
     * @details This function represents @f$\delta\rho_{\rm rem}^{f,\alpha}@f$,
     * given in terms of the flavour-dependent quantity @f$u_f@f$:
     * @f[
     * \delta\rho_{\rm rem}^{f,\alpha}
     * = \frac{\alpha}{4\pi s_W^2}
     * \left\{
     *   - \frac{1}{ c_W^2}
     *     {\rm Re}\big[\overline{\Sigma}^{\prime}_{ZZ}(M_Z^2)\big|_{\mu=M_Z} \big]
     *   - \Delta\overline{\rho}_W\big|_{\mu=M_Z}
     *   + 2u_f
     *   - \left[ \frac{1}{6 c_W^2} - \frac{1}{3} + \frac{3}{4} c_W^2 (1+ c_W^2)
     *            + \frac{9 c_W^2}{4 s_W^2} \right]\ln c_W^2
     *   - \frac{11}{2} + \frac{5}{8} c_W^2(1+ c_W^2)
     * \right\},
     * @f]
     * where @f$u_f@f$ is defined as 
     * @f[
     * u_f = \frac{3v_f^2+a_f^2}{4 c_W^2}{\cal F}_Z(M_Z^2) + {\cal F}_W(M_Z^2)
     * @f]
     * with the so-called unified form factors @f${\cal F}_Z(M_Z^2)@f$ and
     * @f${\cal F}_W(M_Z^2)@f$. 
     * @param[in] u_f the quantity @f$u_f@f$
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\delta\rho_{\rm rem}^{f,\alpha}@f$
     *
     * @sa deltaRho_rem_l() and deltaRho_rem_q() as well as
     * FZ(), FW_l() and FW_q()
     */
    gslpp::complex deltaRho_rem_tmp(const gslpp::complex u_f, const double Mw_i) const;

    /**
     * @brief Remainder contribution of @f$O(\alpha)@f$ to the effective
     * couplings @f$\rho_Z^f@f$,
     * denoted as @f$\delta\rho_{\mathrm{rem}}^{f,\, \alpha}@f$.
     * @details This function handles the remainder contribution
     * @f$\delta\rho_{\mathrm{rem}}^{f,\, \alpha}@f$ for @f$Z\to f\bar{f}@f$.
     * @param[in] f a lepton or quark
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\delta\rho_{\mathrm{rem}}^{f,\, \alpha}@f$
     *
     * @sa deltaRho_rem_tmp()
     */
    gslpp::complex deltaRho_rem_f(const Particle f, const double Mw_i) const;

    /**
     * @brief Remainder contribution of @f$O(\alpha)@f$ to the effective
     * couplings @f$\kappa_Z^f@f$,
     * denoted as @f$\delta\kappa_{\rm rem}^{f,\,\alpha}@f$.
     * @details This function represents @f$\delta\kappa_{\rm rem}^{f,\,\alpha}@f$,
     * given in terms of the flavour-dependent quantities @f$\delta_f@f$ and @f$u_f@f$:
     * @f[
     * \delta\kappa_{\rm rem}^{f,\alpha}
     * = \frac{\alpha}{4\pi s_W^2}
     * \left[ \overline{\Pi}_{Z\gamma}(M_Z^2)\big|_{\mu=M_Z}
     *        + \frac{\delta_f^2}{4 c_W^2}{\cal F}_Z(M_Z^2) - u_f
     *        + \left( \frac{1}{12 c_W^2} + \frac{4}{3} \right)\ln c_W^2
     * \right]
     * @f]
     * where @f$\delta_f@f$ and @f$u_f@f$ are defined as
     * @f{align}{
     * \delta_f &= v_f - a_f = -2Q_f s_W^2\,,
     * \\
     * u_f &= \frac{3v_f^2+a_f^2}{4 c_W^2}{\cal F}_Z(M_Z^2) + {\cal F}_W(M_Z^2)
     * @f}
     * with the so-called unified form factors @f${\cal F}_Z(M_Z^2)@f$ and
     * @f${\cal F}_W(M_Z^2)@f$.
     * @param[in] deltaf the quantity @f$\delta_f@f$
     * @param[in] uf the quantity @f$u_f@f$
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\delta\kappa_{\rm rem}^{f,\,\alpha}@f$
     *
     * @sa deltaKappa_rem_l() and deltaKappa_rem_q() as well as
     * FZ(), FW_l() and FW_q()

     */
    gslpp::complex deltaKappa_rem_tmp(const double deltaf, const gslpp::complex uf,
            const double Mw_i) const;

    /**
     * @brief Remainder contribution of @f$O(\alpha)@f$ to the effective
     * couplings @f$\kappa_Z^f@f$,
     * denoted as @f$\delta\kappa_{\mathrm{rem}}^{f,\, \alpha}@f$.
     * @details This function handles the remainder contribution
     * @f$\delta\kappa_{\mathrm{rem}}^{f,\, \alpha}@f$ for @f$Z\to f\bar{f}@f$.
     * @param[in] f a lepton or quark
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\delta\kappa_{\mathrm{rem}}^{f,\, \alpha}@f$
     *
     * @sa deltaKappa_rem_tmp()
     */
    gslpp::complex deltaKappa_rem_f(const Particle f, const double Mw_i) const;

    /**
     * @brief %EW radiative corrections to the widths of @f$W \to f_i \bar{f}_j@f$, 
     * denoted as @f$\rho^W_{ij}@f$. 
     * @details The factor @f$\rho^W_{ij}@f$ is decomposed as
     * @cite Bardin:1986fi, @cite Bardin:1999ak
     * @f[
     * \rho^W_{ij} = 1 + \delta f^W_{ij} + \delta f^{\rm QED}_{ij},
     * @f]
     * where @f$\delta f^W_{ij}@f$ and @f$\delta f^{\rm QED}_{ij}@f$ are given by
     * @f{align}{
     * \delta f^W_{ij}
     * &=
     * \frac{\alpha}{4\pi s_W^2}
     * \Bigg[
     * - \Delta\overline{\rho}_W \big|_{\mu=M_W}
     * - {\rm Re}\big[\overline{\Sigma}^{\prime}_{WW}(M_W^2)\big]\big|_{\mu=M_W}
     * + \frac{5}{8} c_W^2 (1 + c_W^2)
     * - \frac{11}{2}
     * - \frac{9}{4}\frac{c_W^2}{s_W^2}\ln c_W^2
     * \\
     * &\qquad\qquad\qquad
     * + \left( -1 + \frac{1}{2\,c_W^2} + \frac{2 s_W^4}{c_W^2} Q_i Q_j \right)
     *   \left( V_1(M_W^2,M_Z^2) + \frac{3}{2} \right)
     * + 2\,c_W^2
     * \left( V_2(M_W^2,M_W^2,M_Z^2) + \frac{3}{2} \right)
     * \Bigg]\,,
     * \\
     * \delta f^{\rm QED}_{ij}
     * &=
     * \frac{\alpha}{\pi}
     * \left[
     *   \frac{85}{18} - \frac{\pi^2}{3} + \frac{3}{4}\,Q_iQ_j
     * \right].
     * @f}
     * The functions @f$V_1@f$ and @f$V_2@f$ are given in @cite Bardin:1981sv :
     * @f{align}{
     * V_1(M_W^2,M_Z^2)
     * &=
     * {\rm Re}\left[{\cal F}_{Za}^0(M_W^2)\right] - \frac{3}{2}\,,
     * \\
     * V_2(M_W^2,M_W^2,M_Z^2)
     * &=
     * - 2(2 + c_W^2)M_Z^2\,{\rm Re}\left[C_0(M_W^2;\,M_W,0,M_Z)\right]
     * - \left(\frac{1}{12 c_W^4} + \frac{5}{3 c_W^2} + 1 \right)
     *   {\rm Re}\left[B_0^F(M_W^2;\,M_Z,M_W)\Big|_{\mu=M_W}\right]
     * \\
     * &\quad
     * + \left( \frac{1}{12 c_W^4} + \frac{1}{c_W^2} + 1 \right)\log c_W^2
     * + \frac{1}{12 c_W^4} + \frac{13}{12 c_W^2} + \frac{59}{18}\,,
     * @f}
     * where the imaginary parts have been neglected. 
     * @param[in] Qi the electric charge of f_i
     * @param[in] Qj the electric charge of f_j
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\rho^W_{ij}@f$
     *
     * @attention The masses of virtual fermions are neglected.
     */
    double rho_GammaW_tmp(const double Qi, const double Qj,
            const double Mw_i) const;

    /**
     * @brief %EW radiative corrections to the width of @f$W \to f_i \bar{f}_j@f$,
     * denoted as @f$\rho^W_{ij}@f$.
     * @param[in] fi a lepton or quark
     * @param[in] fj a lepton or quark
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\rho^W_{ij}@f$
     *
     * @sa rho_GammaW_tmp()
     */
    double rho_GammaW(const Particle fi, const Particle fj, const double Mw_i) const;


    ////////////////////////////////////////////////////////////////////////    

    /**
     * @brief The bosonic contribution to the self-energy of the @f$W@f$ boson
     * in the Unitary gauge, @f$\overline{\Sigma}^{\mathrm{bos}}_{WW}(s)@f$.
     * @details This function represents the @f$O(\alpha)@f$ bosonic contribution
     * to @f$\overline{\Sigma}_{WW}(s)@f$,
     * whose definition is given in @ref SeflEnergies "the detailed description"
     * of the current class.
     * @param[in] mu renormalization scale @f$\mu@f$
     * @param[in] s momentum squared @f$s@f$
     * @param[in] Mw_i the @f$W@f$-boson mass @f$M_W@f$
     * @return @f$\overline{\Sigma}^{\mathrm{bos}}_{WW}(s)@f$
     */
    gslpp::complex SigmabarWW_bos(const double mu, const double s, const double Mw_i) const;

    /**
     * @brief The fermionic contribution to the self-energy of the @f$W@f$ boson
     * in the Unitary gauge, @f$\overline{\Sigma}^{\mathrm{fer}}_{WW}(s)@f$.
     * @details This function represents the @f$O(\alpha)@f$ fermionic contribution
     * to @f$\overline{\Sigma}_{WW}(s)@f$,
     * whose definition is given in @ref SeflEnergies "the detailed description"
     * of the current class.
     * @param[in] mu renormalization scale @f$\mu@f$
     * @param[in] s momentum squared @f$s@f$
     * @param[in] Mw_i the @f$W@f$-boson mass @f$M_W@f$
     * @return @f$\overline{\Sigma}^{\mathrm{fer}}_{WW}(s)@f$
     */
    gslpp::complex SigmabarWW_fer(const double mu, const double s, const double Mw_i) const;

    /**
     * @brief The bosonic contribution to the self-energy of the @f$Z@f$ boson
     * in the Unitary gauge, @f$\overline{\Sigma}^{\mathrm{bos}}_{ZZ}(s)@f$.
     * @details This function represents the @f$O(\alpha)@f$ bosonic contribution
     * to @f$\overline{\Sigma}_{ZZ}(s)@f$,
     * whose definition is given in @ref SeflEnergies "the detailed description"
     * of the current class.
     * @param[in] mu renormalization scale @f$\mu@f$
     * @param[in] s momentum squared @f$s@f$
     * @param[in] Mw_i the @f$W@f$-boson mass @f$M_W@f$
     * @return @f$\overline{\Sigma}^{\mathrm{bos}}_{ZZ}(s)@f$
     */
    gslpp::complex SigmabarZZ_bos(const double mu, const double s, const double Mw_i) const;

    /**
     * @brief The fermionic contribution to the self-energy of the @f$Z@f$ boson
     * in the Unitary gauge, @f$\overline{\Sigma}^{\mathrm{fer}}_{ZZ}(s)@f$.
     * @details This function represents the @f$O(\alpha)@f$ fermionic contribution
     * to @f$\overline{\Sigma}_{ZZ}(s)@f$,
     * whose definition is given in @ref SeflEnergies "the detailed description"
     * of the current class.
     * @param[in] mu renormalization scale @f$\mu@f$
     * @param[in] s momentum squared @f$s@f$
     * @param[in] Mw_i the @f$W@f$-boson mass @f$M_W@f$
     * @return @f$\overline{\Sigma}^{\mathrm{fer}}_{ZZ}(s)@f$
     */
    gslpp::complex SigmabarZZ_fer(const double mu, const double s, const double Mw_i) const;

    /**
     * @brief The bosonic contribution to the self-energy of the photon
     * in the Unitary gauge, @f$\overline{\Pi}^{\mathrm{bos}}_{\gamma\gamma}(s)@f$.
     * @details This function represents the @f$O(\alpha)@f$ bosonic contribution
     * to @f$\overline{\Pi}_{\gamma\gamma}(s)@f$,
     * whose definition is given in @ref SeflEnergies "the detailed description"
     * of the current class.
     * @param[in] mu renormalization scale @f$\mu@f$
     * @param[in] s momentum squared @f$s@f$
     * @param[in] Mw_i the @f$W@f$-boson mass @f$M_W@f$
     * @return @f$\overline{\Pi}^{\mathrm{bos}}_{\gamma\gamma}(s)@f$
     */
    gslpp::complex PibarGammaGamma_bos(const double mu, const double s, const double Mw_i) const;

    /**
     * @brief The fermionic contribution to the self-energy of the photon
     * in the Unitary gauge, associated with loops of a lepton or quark. 
     * @f$\overline{\Pi}^{\mathrm{fer},f}_{\gamma\gamma}(s)@f$.
     * @details This function represents the @f$O(\alpha)@f$ fermionic contribution
     * to @f$\overline{\Pi}_{\gamma\gamma}(s)@f$,
     * whose definition is given in @ref SeflEnergies "the detailed description"
     * of the current class.
     * @param[in] mu renormalization scale @f$\mu@f$
     * @param[in] s momentum squared @f$s@f$
     * @param[in] f a lepton or quark
     * @return @f$\overline{\Pi}^{\mathrm{fer},f}_{\gamma\gamma}(s)@f$
     */
    gslpp::complex PibarGammaGamma_fer(const double mu, const double s, const Particle f) const;

    /**
     * @brief The fermionic contribution to the self-energy of the photon
     * in the Unitary gauge,
     * @f$\overline{\Pi}^{\mathrm{fer}}_{\gamma\gamma}(s)@f$.
     * @details This function represents the @f$O(\alpha)@f$ fermionic contribution
     * to @f$\overline{\Pi}_{\gamma\gamma}(s)@f$,
     * whose definition is given in @ref SeflEnergies "the detailed description"
     * of the current class.
     * @param[in] mu renormalization scale @f$\mu@f$
     * @param[in] s momentum squared @f$s@f$
     * @return @f$\overline{\Pi}^{\mathrm{fer}}_{\gamma\gamma}(s)@f$
     */
    gslpp::complex PibarGammaGamma_fer(const double mu, const double s) const;

    /**
     * @brief The bosonic contribution to the self-energy of the @f$Z@f$-@f$\gamma@f$
     * mixing in the Unitary gauge, @f$\overline{\Pi}^{\mathrm{bos}}_{Z\gamma}(s)@f$.
     * @details This function represents the @f$O(\alpha)@f$ bosonic contribution
     * to @f$\overline{\Pi}_{Z\gamma}(s)@f$,
     * whose definition is given in @ref SeflEnergies "the detailed description"
     * of the current class.
     * @param[in] mu renormalization scale @f$\mu@f$
     * @param[in] s momentum squared @f$s@f$
     * @param[in] Mw_i the @f$W@f$-boson mass @f$M_W@f$
     * @return @f$\overline{\Pi}^{\mathrm{bos}}_{Z\gamma}(s)@f$
     */
    gslpp::complex PibarZgamma_bos(const double mu, const double s, const double Mw_i) const;

    /**
     * @brief The fermionic contribution to the self-energy of the @f$Z@f$-@f$\gamma@f$
     * mixing in the Unitary gauge, @f$\overline{\Pi}^{\mathrm{fer}}_{Z\gamma}(s)@f$.
     * @details This function represents the @f$O(\alpha)@f$ fermionic contribution
     * to @f$\overline{\Pi}_{Z\gamma}(s)@f$,
     * whose definition is given in @ref SeflEnergies "the detailed description"
     * of the current class.
     * @param[in] mu renormalization scale @f$\mu@f$
     * @param[in] s momentum squared @f$s@f$
     * @param[in] Mw_i the @f$W@f$-boson mass @f$M_W@f$
     * @return @f$\overline{\Pi}^{\mathrm{fer}}_{Z\gamma}(s)@f$
     */
    gslpp::complex PibarZgamma_fer(const double mu, const double s, const double Mw_i) const;


    ////////////////////////////////////////////////////////////////////////   

    /**
     * @brief The derivative of the bosonic contribution to the self-energy of
     * the @f$W@f$ boson for @f$s=M_W^2@f$ in the Unitary gauge, 
     * @f$\overline{\Sigma}^{\prime,\mathrm{bos}}_{WW}(M_W^2)@f$.
     * @details See also the definition of the self-energy given in
     * @ref SeflEnergies "the detailed description" of the current class.
     * @param[in] mu renormalization scale @f$\mu@f$
     * @param[in] Mw_i the @f$W@f$-boson mass @f$M_W@f$
     * @return @f$\overline{\Sigma}^{\prime,\mathrm{bos}}_{WW}(M_W^2)@f$
     */
    gslpp::complex SigmabarPrime_WW_bos_Mw2(const double mu, const double Mw_i) const;

    /**
     * @brief The derivative of the fermionic contribution to the self-energy of
     * the @f$W@f$ boson for @f$s=M_W^2@f$ in the Unitary gauge,
     * @f$\overline{\Sigma}^{\prime,\mathrm{fer}}_{WW}(M_W^2)@f$.
     * @details See also the definition of the self-energy given in
     * @ref SeflEnergies "the detailed description" of the current class.
     * @param[in] mu renormalization scale @f$\mu@f$
     * @param[in] Mw_i the @f$W@f$-boson mass @f$M_W@f$
     * @return @f$\overline{\Sigma}^{\prime,\mathrm{fer}}_{WW}(M_W^2)@f$
     */
    gslpp::complex SigmabarPrime_WW_fer_Mw2(const double mu, const double Mw_i) const;

    /**
     * @brief The derivative of the bosonic contribution to the self-energy of
     * the @f$Z@f$ boson for @f$s=M_Z^2@f$ in the Unitary gauge,
     * @f$\overline{\Sigma}^{\prime,\mathrm{bos}}_{ZZ}(M_Z^2)@f$.
     * @details See also the definition of the self-energy given in
     * @ref SeflEnergies "the detailed description" of the current class.
     * @param[in] mu renormalization scale @f$\mu@f$
     * @param[in] Mw_i the @f$W@f$-boson mass @f$M_W@f$
     * @return @f$\overline{\Sigma}^{\prime,\mathrm{bos}}_{ZZ}(M_Z^2)@f$
     */
    gslpp::complex SigmabarPrime_ZZ_bos_Mz2(const double mu, const double Mw_i) const;

    /**
     * @brief The derivative of the fermionic contribution to the self-energy of
     * the @f$Z@f$ boson for @f$s=M_Z^2@f$ in the Unitary gauge,
     * @f$\overline{\Sigma}^{\prime,\mathrm{fer}}_{ZZ}(M_Z^2)@f$.
     * @details See also the definition of the self-energy given in
     * @ref SeflEnergies "the detailed description" of the current class.
     * @param[in] mu renormalization scale @f$\mu@f$
     * @param[in] Mw_i the @f$W@f$-boson mass @f$M_W@f$
     * @return @f$\overline{\Sigma}^{\prime,\mathrm{fer}}_{ZZ}(M_Z^2)@f$
     */
    gslpp::complex SigmabarPrime_ZZ_fer_Mz2(const double mu, const double Mw_i) const;

    ////////////////////////////////////////////////////////////////////////       

    /**
     * @brief @f$\Delta\overline{\rho}@f$.
     * @details The quantity @f$\Delta\overline{\rho}@f$, which is associated with
     * @f$\Delta\rho@f$ as explained in the description of DeltaRho(), is
     * defined as
     * @f[
     * \Delta\overline{\rho}\big|_{\mu}
     * =
     * \frac{1}{M_W^2}\left[ {\rm Re}\,\overline{\Sigma}_{WW}(M_W^2)\big|_{\mu}
     *  - {\rm Re}\,\overline{\Sigma}_{ZZ}(M_Z^2)\big|_{\mu} \right],
     * @f]
     * where @f$\mu@f$ denotes the renormalization scale.
     * @param[in] mu renormalization scale @f$\mu@f$
     * @param[in] Mw_i the @f$W@f$-boson mass @f$M_W@f$
     * @return @f$\Delta\overline{\rho}\big|_{\mu}@f$
     */
    double DeltaRhobar(const double mu, const double Mw_i) const;

    /**
     * @brief @f$\Delta\overline{\rho}_W@f$.
     * @details The quantity @f$\Delta\overline{\rho}_W@f$ is defined as
     * @f[
     * \Delta\overline{\rho}_W\big|_{\mu}
     * =
     * \frac{1}{M_W^2}\left[ \overline{\Sigma}_{WW}(0)\big|_{\mu}
     *  - {\rm Re}\,\overline{\Sigma}_{WW}(M_W^2)\big|_{\mu} \right],
     * @f]
     * where @f$\mu@f$ denotes the renormalization scale.
     * @param[in] mu renormalization scale @f$\mu@f$
     * @param[in] Mw_i the @f$W@f$-boson mass @f$M_W@f$
     * @return @f$\Delta\overline{\rho}_W\big|_{\mu}@f$
     */
    double DeltaRhobarW(const double mu, const double Mw_i) const;


    ////////////////////////////////////////////////////////////////////////   

    /**
     * @brief A test function.
     * @details @f$\Delta\overline{\rho}^{\mathrm{bos}}@f$ is given without
     * the use of the self-energies of the gauge bosons.
     * @param[in] Mw_i the @f$W@f$-boson mass @f$M_W@f$
     * @return @f$\Delta\overline{\rho}^{\mathrm{bos}}@f$
     *
     * @attention The renormalization scale is fixed to be @f$\mu=M_W@f$.
     */
    double TEST_DeltaRhobar_bos(const double Mw_i) const;

    /**
     * @brief A test function. 
     * @details @f$\Delta\overline{\rho}_W^{\mathrm{bos}}@f$ is given without
     * the use of the self-energies of the gauge bosons.
     * @param[in] Mw_i the @f$W@f$-boson mass @f$M_W@f$
     * @return @f$\Delta\overline{\rho}_W^{\mathrm{bos}}@f$
     *
     * @attention The renormalization scale is fixed to be @f$\mu=M_W@f$.
     */
    double TEST_DeltaRhobarW_bos(const double Mw_i) const;


    ////////////////////////////////////////////////////////////////////////    

    /**
     * @brief The form factor @f$\mathcal{F}_{Za}^0@f$.
     * @details The form factor @f$\mathcal{F}_{Za}^0@f$,
     * associated with abelian-type diagrams of @f$Z\to f\bar{f}@f$
     * with a virtual @f$Z@f$ boson, is given in the chiral limit by
     * @f[
     *  {\cal F}_{Za}^0(s)
     * = 2(R_Z + 1)^2 s\, C_0(s;0,(\widetilde{M}_Z^2)^{1/2},0)
     *   - (2R_Z+3) \ln\left(-\frac{\widetilde{M}_Z^2}{s}\right)
     *   - 2R_Z - \frac{7}{2}\,,
     * @f]
     * where @f$\widetilde{M}_Z^2 \equiv M_Z^2 - i\, M_Z\Gamma_Z\approx M_Z^2-i\epsilon@f$,
     * and the definitions of the other symbols can be read from the codes below.
     * 
     * See @cite Bardin:1999ak.
     * @param[in] s momentum squared @f$s@f$
     * @param[in] Mw_i the @f$W@f$-boson mass @f$M_W@f$
     * @return @f$\mathcal{F}_{Za}^0@f$
     */
    gslpp::complex FZa_0(const double s, const double Mw_i) const;

    /**
     * @brief The form factor @f$\mathcal{F}_{Wa}^0@f$.
     * @details The form factor @f$\mathcal{F}_{Wa}^0@f$,
     * associated with abelian-type diagrams of @f$Z\to f\bar{f}@f$
     * with a virtual @f$W@f$ boson, is given in the chiral limit by
     * @f[
     * {\cal F}_{Wa}^0(s)
     * = 2(R_W + 1)^2 s\, C_0(s;0,(\widetilde{M}_W^2)^{1/2},0)
     *   - (2R_W+3) \ln\left(-\frac{\widetilde{M}_W^2}{s}\right)
     *   - 2R_W - \frac{7}{2}\,,
     * @f]
     * where @f$\widetilde{M}_W^2 \equiv M_W^2 - i\, M_W\Gamma_W\approx M_W^2-i\epsilon@f$,
     * and the definitions of the other symbols can be read from the codes below.
     *
     * See @cite Bardin:1999ak.
     * @param[in] s momentum squared @f$s@f$
     * @param[in] Mw_i the @f$W@f$-boson mass @f$M_W@f$
     * @return @f$\mathcal{F}_{Wa}^0@f$
     */
    gslpp::complex FWa_0(const double s, const double Mw_i) const;

    /**
     * @brief The form factor @f$\overline{\mathcal{F}}_{Wa}^0@f$.
     * @details The form factor @f$\overline{\mathcal{F}}_{Wa}^0@f$,
     * associated with abelian-type diagrams of @f$Z\to f\bar{f}@f$
     * with a virtual @f$W@f$ boson, is given in the chiral limit by
     * @f[
     * \overline{{\cal F}}_{Wa}^0(s) = 0\,. 
     * @f]
     *
     * See @cite Bardin:1999ak.
     * @param[in] s momentum squared @f$s@f$
     * @return @f$\overline{\mathcal{F}}_{Wa}^0@f$
     */
    gslpp::complex FbarWa_0(const double s) const;

    /**
     * @brief The form factor @f$\mathcal{F}_{Wn}^0@f$.
     * @details The form factor @f$\mathcal{F}_{Wn}^0@f$,
     * associated with nonabelian-type diagrams of @f$Z\to f\bar{f}@f$
     * with virtual @f$W@f$ bosons, is given in the chiral limit by
     * @f[
     * {\cal F}_{Wn}^0(s)
     * = - 2 (R_W + 2) M_W^2 C_0(s;M_W,0,M_W)
     *   - \left( 2R_W + \frac{7}{3} - \frac{3}{2R_W} - \frac{1}{12R_W^2} \right)
     *     B_0(s; M_W, M_W)\big|_{\mu=M_W}
     *   + 2R_W + \frac{9}{2} - \frac{11}{18R_W} + \frac{1}{18R_W^2}\,,
     * @f]
     * where the definitions of the symbols can be read from the codes below.
     *
     * See @cite Bardin:1999ak.
     * @param[in] s momentum squared @f$s@f$
     * @param[in] Mw_i the @f$W@f$-boson mass @f$M_W@f$
     * @return @f$\mathcal{F}_{Wn}^0@f$
     */
    gslpp::complex FWn_0(const double s, const double Mw_i) const;

    /**
     * @brief The form factor @f$\mathcal{F}_{Wa}^t@f$.
     * @details The form factor @f$\mathcal{F}_{Wa}^t@f$,
     * associated with abelian-type diagrams of @f$Z\to f\bar{f}@f$
     * with a virtual @f$W@f$ boson and the heavy top quark, is given by
     * @f{align}{
     * {\cal F}_{Wa}^t(s)
     * &=
     * 2(R_W +1)^2 s
     * \left[ C_0(s;M_t,M_W,M_t) - C_0(s;0,(\widetilde{M}_W^2)^{1/2},0) \right]
     * + (2R_W + 3) \left[ - B_0(s;M_t,M_t)\big|_{\mu=M_W}
     *   + \ln\left(-\frac{\widetilde{M}_W^2}{s}\right) + 2 \right]
     * \\
     * &\quad
     * - w_t
     * \Bigg\{
     *   \left(  3R_W +2-w_t - w_t^2 R_W  \right) M_W^2 C_0(s;M_t,M_W,M_t)
     *   + \left( R_W  + \frac{1}{2} + w_t R_W \right)
     *     \left[ 1 - B_0(s;M_t,M_t)\big|_{\mu=M_W} \right]
     * \\
     * &\qquad\qquad
     *   - \left( 2 R_W + \frac{1}{2} - \frac{2}{w_t - 1}
     *   + \frac{3}{2}\frac{1}{(w_t - 1)^2} + w_t R_W \right) \ln w_t
     *   + \frac{3}{2}\frac{1}{w_t - 1} + \frac{3}{4}
     * \Bigg\}\,,
     * @f}
     * where @f$\widetilde{M}_W^2 \equiv M_W^2 - i\, M_W\Gamma_W\approx M_W^2-i\epsilon@f$,
     * and the definitions of the other symbols can be read from the codes below.
     *
     * See @cite Bardin:1999ak.
     * @param[in] s momentum squared @f$s@f$
     * @param[in] Mw_i the @f$W@f$-boson mass @f$M_W@f$
     * @return @f$\mathcal{F}_{Wa}^t@f$
     */
    gslpp::complex FWa_t(const double s, const double Mw_i) const;

    /**
     * @brief The form factor @f$\overline{\mathcal{F}}_{Wa}^t@f$.
     * @details The form factor @f$\overline{\mathcal{F}}_{Wa}^t@f$,
     * associated with abelian-type diagrams of @f$Z\to f\bar{f}@f$
     * with a virtual @f$W@f$ boson and the heavy top quark, is given by
     * @f[
     * \overline{{\cal F}}_{Wa}^t(s)
     * =
     * - w_t
     * \Bigg\{ \left[ R_W + 2 - w_t(2 - w_t)R_W  \right]
     *   M_W^2C_0(s;M_t,M_W,M_t)
     *   - \left( \frac{1}{2} - R_W + w_t R_W \right)
     *     \left[ - B_0(s;M_t,M_t)\big|_{\mu=M_W} + 1 \right] + w_t R_W \ln w_t
     * \Bigg\}\,,
     * @f]
     * where the definitions of the symbols can be read from the codes below.
     *
     * See @cite Bardin:1999ak.
     * @param[in] s momentum squared @f$s@f$
     * @param[in] Mw_i the @f$W@f$-boson mass @f$M_W@f$
     * @return @f$\overline{\mathcal{F}}_{Wa}^t@f$
     */
    gslpp::complex FbarWa_t(const double s, const double Mw_i) const;

    /**
     * @brief The form factor @f$\mathcal{F}_{Wn}^t@f$.
     * @details The form factor @f$\mathcal{F}_{Wn}^t@f$,
     * associated with nonabelian-type diagrams of @f$Z\to f\bar{f}@f$
     * with virtual @f$W@f$ bosons and the heavy top quark, is given by
     * @f{align}{
     * {\cal F}_{Wn}^t(s)
     * &=
     * -2 (R_W +2) M_W^2
     * \left[ C_0(s;M_W,M_t,M_W) - C_0(s;M_W,0,M_W) \right]
     * \\
     * &\quad
     * + w_t
     * \Bigg\{
     *   \Bigg[ 3 R_W + \frac{5}{2} - \frac{2}{R_W }
     *     - w_t \left( 2 - \frac{1}{2R_W} \right)
     *  + w_t ^2 \left( \frac{1}{2} - R_W \right)
     *   \Bigg] M_W^2 C_0(s;M_W,M_t,M_W)
     * \\
     * &\qquad\qquad
     * +\left[ R_W + 1 - \frac{1}{4 R_W}
     *   - w_t \left( \frac{1}{2} - R_W \right) \right]
     *   \left[ B_0(s;M_W,M_W)\big|_{\mu=M_W} - 1 \right]
     * \\
     * &\qquad\qquad
     *   + \left[ 2 R_W + \frac{1}{2}
     *     - \frac{2}{w_t - 1} + \frac{3}{2}\frac{1}{(w_t - 1)^2}
     *     - w_t \left( \frac{1}{2} - R_W \right) \right] \ln w_t
     *   - \frac{3}{2}\frac{1}{w_t -1} + \frac{1}{4}
     *   - \frac{1}{2 R_W}
     * \Bigg\}\,,
     * @f}
     * where the definitions of the symbols can be read from the codes below.
     *
     * See @cite Bardin:1999ak.
     * @param[in] s momentum squared @f$s@f$
     * @param[in] Mw_i the @f$W@f$-boson mass @f$M_W@f$
     * @return @f$\mathcal{F}_{Wn}^0@f$
     */
    gslpp::complex FWn_t(const double s, const double Mw_i) const;

    /**
     * @brief The unified form factor @f$\mathcal{F}_Z@f$.
     * @details The so-called unified form factor @f$\mathcal{F}_Z@f$, associated
     * with radiative corrections to the @f$Z\to f\bar{f}@f$ vertex with a virtual
     * @f$Z@f$ boson, is given by
     * @f[
     * {\cal F}_Z(s) = {\cal F}_{Za}^0(s)\,,
     * @f]
     * where @f${\cal F}_{Za}^0(s)@f$ corresponds to the function FZa_0().
     * 
     * See @cite Bardin:1999ak.
     * @param[in] s momentum squared @f$s@f$
     * @param[in] Mw_i the @f$W@f$-boson mass @f$M_W@f$
     * @return @f$\mathcal{F}_Z@f$
     */
    gslpp::complex FZ(const double s, const double Mw_i) const;

    /**
     * @brief The unified form factor @f$\mathcal{F}_W@f$ for @f$Z\to f\bar{f}@f$.
     * @details The so-called unified form factor @f$\mathcal{F}_W@f$, associated
     * with radiative corrections to the @f$Z\to f\bar{f}@f$ vertex with a virtual
     * @f$W@f$ boson as well as with virtual @f$W@f$ bosons, is given by
     * @f[
     * {\cal F}_W(s) = c_W^2 {\cal F}_{Wn}^0(s)
     * - \frac{1}{2}\sigma_{l'}^a {\cal F}_{Wa}^0(s)
     * - \frac{1}{2}\overline{{\cal F}}_{Wa}^0(s)\,,
     * @f]
     * where the suprescripts "0" denote the chiral limit,
     * @f$\sigma_{f'}^a = |v_{f'} + a_{f'}|
     * = 1 - 2|Q_{f'}|s_W^2 = 2c_W^2 - 1 + 2|Q_{f}| s_W^2@f$
     * with @f$f'@f$ being the partner of @f$f@f$ in the @f$SU(2)_L@f$ doublet,
     * and @f${\cal F}_{Wn}^0(s)@f$, @f${\cal F}_{Wa}^0(s)@f$
     * and @f$\overline{{\cal F}}_{Wa}^0(s)@f$ correspond to the functions 
     * FWn_0(), FWa_0() and FbarWa_0(), respectively.
     *
     * See @cite Bardin:1999ak.
     * @param[in] s momentum squared @f$s@f$
     * @param[in] f a lepton or quark
     * @param[in] Mw_i the @f$W@f$-boson mass @f$M_W@f$
     * @return @f$\mathcal{F}_W@f$
     */
    gslpp::complex FW(const double s, const Particle f, const double Mw_i) const;

    ////////////////////////////////////////////////////////////////////////        

    /**
     * @brief A test function for @f$\mathcal{F}_{Wn}@f$ with a finite fermion mass.
     * @param[in] s momentum squared @f$s@f$
     * @param[in] mf the mass of the fermion in the loop
     * @param[in] Mw_i the @f$W@f$-boson mass @f$M_W@f$
     * @return @f$\mathcal{F}_{Wn}@f$
     */
    gslpp::complex TEST_FWn(const double s, const double mf, const double Mw_i) const;


    ////////////////////////////////////////////////////////////////////////    

private:
    const EWSMcache& cache; ///< A reference to an object of type EWSMcache.

};

#endif	/* EWSMONELOOPEW_H */

