/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef PVFUNCTIONS_H
#define	PVFUNCTIONS_H

// set in case where LoopTools library is employed.
//#define USE_LOOPTOOLS

#include "gslpp.h"
#include "Polylogarithms.h"
#ifdef USE_LOOPTOOLS
#include "LoopToolsWrapper.h"
#endif

/**
 * @class PVfunctions
 * @ingroup LoopFunctions 
 * @brief A class for Passarino-Veltman functions.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class handles the so-called Passarino-Veltman (PV) functions,
 * which appear in one-loop amplitudes.
 * The definitions of the two-point and four-point functions used in the current 
 * class are identical to those in
 * <a href="http://www.feynarts.de/looptools/" target=blank>LoopTools</a> library
 * @cite Hahn:1998yk.
 * The one-point and three-point functions are identical to those in 
 * <a href="http://www.feynarts.de/looptools/" target=blank>LoopTools</a> library
 * when bExtraMinusSign is set to false. On the other hand, when bExtraMinusSign
 * is set to true, an extra minus sign is added to them in order to match their
 * definitions to those in @cite Bardin:1999ak.
 * If the preprocessor macro USE_LOOPTOOLS is defined in PVfunctions.h or
 * Makefile, the functions in LoopTools library, called via LoopToolsWrapper
 * class, are employed instead of those defined in the current class.
 *
 * See, e.g., @cite tHooft:1978xw, @cite Passarino:1978jh, @cite VeltBook,
 * @cite Denner:1991kt and @cite Bardin:1999ak
 */
class PVfunctions {
public:

    /**
     * @brief Constructor.
     * @details The boolean argument bExtraMinusSign controls whether an
     * extra overall minus sign is added to the one-point and three-point
     * functions or not. See also the detailed description of the current class.
     * @param[in] bExtraMinusSign a flag to control whether an extra overall
     * minus sign is added to the one-point and three-point functions or not
     */
    PVfunctions(const bool bExtraMinusSign);

    /**
     * @brief @f$A_0(m^2)@f$.
     * @details The scalar one-point function @f$A_0(m^2)@f$ is defined as
     * @f[
     * A_0(m^2)
     * = \frac{(2\pi\mu)^{4-d}}{i\pi^2}\int d^dk\, \frac{1}{k^2-m^2+i\varepsilon},
     * @f]
     * where the UV divergence is regularized with the dimensional regularization.
     * When bExtraMinusSign=true is passed to the constructor, an extra overall
     * minus sign is added to the above definition. 
     * @param[in] mu2 the renormalization scale squared, @f$\mu^2@f$
     * @param[in] m2 mass squared, @f$m^2@f$
     * @return the finite part of @f$A_0(m^2)@f$
     * in the sense of the @f$\overline{\mathrm{MS}}@f$ scheme
     */
    double A0(const double mu2, const double m2) const;
    
    /**
     * @brief @f$B_0(p^2; m_0^2, m_1^2)@f$.
     * @details The scalar two-point function @f$B_0(p^2; m_0^2, m_1^2)@f$
     * is defined as
     * @f[
     * B_0(p^2;m_0^2,m_1^2)
     * = \frac{(2\pi\mu)^{4-d}}{i\pi^2}\int d^dk\,
     * \frac{1}{(k^2-m_0^2+i\varepsilon)\left[(k+p)^2-m_1^2+i\varepsilon\right]}\,,
     * @f]
     * where the UV divergence is regularized with the dimensional regularization.
     * @param[in] mu2 the renormalization scale squared, @f$\mu^2@f$
     * @param[in] p2 momentum squared, @f$p^2@f$
     * @param[in] m02, m12 mass squared, @f$m_0^2@f$ and @f$m_1^2@f$
     * @return the finite part of @f$B_0(p^2; m_0^2, m_1^2)@f$
     * in the sense of the @f$\overline{\mathrm{MS}}@f$ scheme
     */
    gslpp::complex B0(const double mu2, const double p2,
               const double m02, const double m12) const;
    
    /**
     * @brief @f$B_1(p^2; m_0^2, m_1^2)@f$.
     * @details The vector two-point PV function @f$B_1(p^2; m_0^2, m_1^2)@f$
     * is defined as
     * @f[
     * p_\mu B_1(p^2;m_0^2,m_1^2)
     * = \frac{(2\pi\mu)^{4-d}}{i\pi^2}\int d^dk\,
     * \frac{k_\mu}{(k^2-m_0^2+i\varepsilon)\left[(k+p)^2-m_1^2+i\varepsilon\right]},
     * @f]
     * where the UV divergence is regularized with the dimensional regularization.
     * @param[in] mu2 the renormalization scale squared, @f$\mu^2@f$
     * @param[in] p2 momentum squared, @f$p^2@f$
     * @param[in] m02, m12 mass squared, @f$m_0^2@f$ and @f$m_1^2@f$
     * @return the finite part of @f$B_1(p^2; m_0^2, m_1^2)@f$
     * in the sense of the @f$\overline{\mathrm{MS}}@f$ scheme
     */
    gslpp::complex B1(const double mu2, const double p2,
               const double m02, const double m12) const;
    
    /**
     * @brief @f$B_{11}(p^2; m_0^2, m_1^2)@f$.
     * @details The tensor two-point PV function @f$B_{11}(p^2; m_0^2, m_1^2)@f$
     * is defined as
     * @f[
     * g_{\mu\nu} B_{00}(p^2;m_0^2,m_1^2) + p_\mu p_\nu B_{11}(p^2;m_0^2,m_1^2)
     * =
     * \frac{(2\pi\mu)^{4-d}}{i\pi^2}\int d^dk\,
     * \frac{k_\mu k_\nu}{(k^2-m_0^2+i\varepsilon)
     * \left[(k+p)^2-m_1^2+i\varepsilon\right]}, 
     * @f]
     * where the UV divergence is regularized with the dimensional regularization.
     * @param[in] mu2 the renormalization scale squared, @f$\mu^2@f$
     * @param[in] p2 momentum squared, @f$p^2@f$
     * @param[in] m02, m12 mass squared, @f$m_0^2@f$ and @f$m_1^2@f$
     * @return the finite part of @f$B_{11}(p^2; m_0^2, m_1^2)@f$
     * in the sense of the @f$\overline{\mathrm{MS}}@f$ scheme
     */
    gslpp::complex B11(const double mu2, const double p2,
                const double m02, const double m12) const;
    
    /**
     * @brief @f$B_{00}(p^2; m_0^2, m_1^2)@f$.
     * @details The tensor two-point PV function @f$B_{00}(p^2; m_0^2, m_1^2)@f$
     * is defined as
     * @f[
     * g_{\mu\nu} B_{00}(p^2;m_0^2,m_1^2) + p_\mu p_\nu B_{11}(p^2;m_0^2,m_1^2)
     * =
     * \frac{(2\pi\mu)^{4-d}}{i\pi^2}\int d^dk\,
     * \frac{k_\mu k_\nu}{(k^2-m_0^2+i\varepsilon)
     * \left[(k+p)^2-m_1^2+i\varepsilon\right]}, 
     * @f]
     * where the UV divergence is regularized with the dimensional regularization.
     * @param[in] mu2 the renormalization scale squared, @f$\mu^2@f$
     * @param[in] p2 momentum squared, @f$p^2@f$
     * @param[in] m02, m12 mass squared, @f$m_0^2@f$ and @f$m_1^2@f$
     * @return the finite part of @f$B_{00}(p^2; m_0^2, m_1^2)@f$
     * in the sense of the @f$\overline{\mathrm{MS}}@f$ scheme
     */
    gslpp::complex B00(const double mu2, const double p2,
                const double m02, const double m12) const;
    
    /**
     * @brief @f$B_{f}(p^2; m_0^2, m_1^2)@f$.
     * @details The function @f$B_{f}(p^2; m_0^2, m_1^2)@f$ is defined as a sum
     * of the two PV functions:
     * @f[
     * B_f(p^2;m_0^2,m_1^2)
     * = 2 \left[ B_{11}(p^2;m_0^2,m_1^2) + B_{1}(p^2;m_0^2,m_1^2) \right],
     * @f]
     * where the UV divergence is regularized with the dimensional regularization.
     * @param[in] mu2 the renormalization scale squared, @f$\mu^2@f$
     * @param[in] p2 momentum squared, @f$p^2@f$
     * @param[in] m02, m12 mass squared, @f$m_0^2@f$ and @f$m_1^2@f$
     * @return the finite part of @f$B_{f}(p^2; m_0^2, m_1^2)@f$
     * in the sense of the @f$\overline{\mathrm{MS}}@f$ scheme
     */
    gslpp::complex Bf(const double mu2, const double p2,
               const double m02, const double m12) const;
    
    /**
     * @brief @f$B_{0p}(p^2; m_0^2, m_1^2)@f$.
     * @details The function @f$B_{0p}(p^2; m_0^2, m_1^2)@f$ is defined as
     * @f[
     * B_{0p}(p^2;m_0^2,m_1^2) = \frac{\partial}{\partial p^2} B_0(p^2;m_0^2,m_1^2)\,,
     * @f]
     * which is UV finite, while @f$B_{0p}(m^2; 0, m^2)@f$ is IR divergent.
     * The IR divergence is regularized with the dimensional regularization. 
     * @param[in] muIR2 the renormalization scale squared for the IR divergence, @f$\mu_{\mathrm{IR}}^2@f$
     * @param[in] p2 momentum squared, @f$p^2@f$
     * @param[in] m02, m12 mass squared, @f$m_0^2@f$ and @f$m_1^2@f$
     * @return the finite part of @f$B_{0p}(p^2; m_0^2, m_1^2)@f$
     * in the sense of the @f$\overline{\mathrm{MS}}@f$ scheme
     */
    gslpp::complex B0p(const double muIR2, const double p2,
                const double m02, const double m12) const;
    
    /**
     * @brief @f$B_{1p}(p^2; m_0^2, m_1^2)@f$.
     * @details The function @f$B_{1p}(p^2; m_0^2, m_1^2)@f$ is defined as
     * @f[
     * B_{1p}(p^2;m_0^2,m_1^2) = \frac{\partial}{\partial p^2} B_1(p^2;m_0^2,m_1^2)\,,
     * @f]
     * where the UV divergence is regularized with the dimensional regularization.
     * @param[in] mu2 the renormalization scale squared, @f$\mu^2@f$
     * @param[in] p2 momentum squared, @f$p^2@f$
     * @param[in] m02, m12 mass squared, @f$m_0^2@f$ and @f$m_1^2@f$
     * @return the finite part of @f$B_{1p}(p^2; m_0^2, m_1^2)@f$
     * in the sense of the @f$\overline{\mathrm{MS}}@f$ scheme
     */
    gslpp::complex B1p(const double mu2, const double p2,
                const double m02, const double m12) const;
    
    /**
     * @brief @f$B_{11p}(p^2; m_0^2, m_1^2)@f$.
     * @details The function @f$B_{11p}(p^2; m_0^2, m_1^2)@f$ is defined as
     * @f[
     * B_{11p}(p^2;m_0^2,m_1^2) = \frac{\partial}{\partial p^2} B_{11}(p^2;m_0^2,m_1^2)\,,
     * @f]
     * where the UV divergence is regularized with the dimensional regularization.
     * @param[in] mu2 the renormalization scale squared, @f$\mu^2@f$
     * @param[in] p2 momentum squared, @f$p^2@f$
     * @param[in] m02, m12 mass squared, @f$m_0^2@f$ and @f$m_1^2@f$
     * @return the finite part of @f$B_{11p}(p^2; m_0^2, m_1^2)@f$
     * in the sense of the @f$\overline{\mathrm{MS}}@f$ scheme
     */
    gslpp::complex B11p(const double mu2, const double p2,
                 const double m02, const double m12) const;

    /**
     * @brief @f$B_{00p}(p^2; m_0^2, m_1^2)@f$.
     * @details The function @f$B_{00p}(p^2; m_0^2, m_1^2)@f$ is defined as
     * @f[
     * B_{00p}(p^2;m_0^2,m_1^2) = \frac{\partial}{\partial p^2} B_{00}(p^2;m_0^2,m_1^2)\,,
     * @f]
     * where the UV divergence is regularized with the dimensional regularization.
     * @param[in] mu2 the renormalization scale squared, @f$\mu^2@f$
     * @param[in] p2 momentum squared, @f$p^2@f$
     * @param[in] m02, m12 mass squared, @f$m_0^2@f$ and @f$m_1^2@f$
     * @return the finite part of @f$B_{00p}(p^2; m_0^2, m_1^2)@f$
     * in the sense of the @f$\overline{\mathrm{MS}}@f$ scheme
     */
    gslpp::complex B00p(const double mu2, const double p2,
                 const double m02, const double m12) const;
    
    /**
     * @brief @f$B_{fp}(p^2; m_0^2, m_1^2)@f$.
     * @details The function @f$B_{fp}(p^2; m_0^2, m_1^2)@f$ is defined as
     * @f[
     * B_{fp}(p^2;m_0^2,m_1^2) = \frac{\partial}{\partial p^2} B_{f}(p^2;m_0^2,m_1^2)\,,
     * @f]
     * where the UV divergence is regularized with the dimensional regularization.
     * @param[in] mu2 the renormalization scale squared, @f$\mu^2@f$
     * @param[in] p2 momentum squared, @f$p^2@f$
     * @param[in] m02, m12 mass squared, @f$m_0^2@f$ and @f$m_1^2@f$
     * @return the finite part of @f$B_{fp}(p^2; m_0^2, m_1^2)@f$
     * in the sense of the @f$\overline{\mathrm{MS}}@f$ scheme
     */
    gslpp::complex Bfp(const double mu2, const double p2,
                const double m02, const double m12) const;

    /**
     * @brief @f$C_{0}(0,0,p^2; m_0^2, m_1^2, m_2^2)@f$.
     * @details The scalar three-point function 
     * @f$C_{0}(p_1^2,p_2^2,(p_1+p_2)^2; m_0^2, m_1^2, m_2^2)@f$ is defined as
     * @f[
     * C_0(p_1^2,p_2^2,(p_1+p_2)^2; m_0^2,m_1^2,m_2^2)
     * = \frac{1}{i\pi^2}\int d^4k\,
     * \frac{1}{(k^2-m_0^2+i\varepsilon)
     * \left[(k+p_1)^2-m_1^2+i\varepsilon\right]
     * \left[(k+p_1+p_2)^2-m_2^2+i\varepsilon\right]}\,,
     * @f]
     * The current functions handles only the special case of @f$p_1^2=p_2^2=0@f$.
     * When bExtraMinusSign=true is passed to the constructor, an extra overall
     * minus sign is added to the above definition. 
     * @param[in] p2 momentum squared, @f$p^2@f$
     * @param[in] m02, m12, m22 mass squared, @f$m_0^2@f$, @f$m_1^2@f$ and @f$m_2^2@f$
     * @return @f$C_{0}(0,0,p^2; m_0^2, m_1^2, m_2^2)@f$
     */
    gslpp::complex C0(const double p2, 
               const double m02, const double m12, const double m22) const;

    /**
     * @brief @f$C_{11}(m_1^2, m_2^2, m_3^2)@f$.
     * @details The function
     * @f$C_{11}(m_1^2, m_2^2, m_3^2)@f$ is defined as
     * @f[
     * C_{11}(m_1^2,m_2^2,m_3^2)
     * = \frac{m_1^4 m_2^2 (2 m_1^2-m_2^2) \log \left(\frac{m_1^2}{m_2^2}\right)
     * +m_1^4 m_3^2 (m_3^2-2 m_1^2) \log \left(\frac{m_1^2}{m_3^2}\right)
     * -m_1^2 (m_1^2-m_2^2) (m_1^2-m_3^2) (m_2^2-m_3^2)
     * +m_2^2 m_3^2 (m_2^2-2 m_1^2) (m_3^2-2 m_1^2) \log \left(\frac{m_2^2}{m_3^2}\right)}
     * {2 (m_1^2-m_2^2)^2 (m_1^2-m_3^2)^2 (m_2^2-m_3^2)}.
     * @f]
     * The definition is taken from Equation (B8) in @cite Arganda:2005ji.
     * @param[in] m12, m22, m32 mass squared, @f$m_1^2@f$, @f$m_2^2@f$ and @f$m_3^2@f$
     * @return @f$C_{11}(m_1^2, m_2^2, m_3^2)@f$
     */
    double C11(const double m12, const double m22, const double m32) const;

    /**
     * @brief @f$C_{12}(m_1^2, m_2^2, m_3^2)@f$.
     * @details The function
     * @f$C_{12}(m_1^2, m_2^2, m_3^2)@f$ is defined as
     * @f[
     * C_{12}(m_1^2,m_2^2,m_3^2)
     * = \frac{m_1^4 \left(m_2^4 \log \left(\frac{m_1^2}{m_2^2}\right)
     * +m_3^2 (m_3^2-2 m_2^2) \log \left(\frac{m_1^2}{m_3^2}\right)\right)
     * +m_2^4 m_3^2 (2 m_1^2-m_3^2) \log \left(\frac{m_2^2}{m_3^2}\right)
     * +m_3^2 (m_1^2-m_2^2) (m_1^2-m_3^2) (m_2^2-m_3^2)}
     * {2 (m_1^2-m_2^2) (m_1^2-m_3^2)^2 (m_2^2-m_3^2)^2}.
     * @f]
     * The definition is taken from Equation (B9) in @cite Arganda:2005ji.
     * @param[in] m12, m22, m32 mass squared, @f$m_1^2@f$, @f$m_2^2@f$ and @f$m_3^2@f$
     * @return @f$C_{12}(m_1^2, m_2^2, m_3^2)@f$
     */
    double C12(const double m12, const double m22, const double m32) const;

    /**
     * @brief @f$D_{0}(0,0,0,0,s,t; m_0^2, m_1^2, m_2^2, m_3^2)@f$.
     * @details The scalar four-point function 
     * @f$D_{0}(p_1^2,p_2^2,p_3^2,p_4^2,(p_1+p_2)^2,(p_2+p_3)^2; m_0^2, m_1^2, m_2^2, m_3^2)@f$
     * is defined as
     * @f{eqnarray*}{
     * &&D_0(p_1^2,p_2^2,p_3^2,p_4^2,(p_1+p_2)^2,(p_2+p_3)^2; m_0^2,m_1^2,m_2^2,m_3^2)
     * \\
     * &&\quad
     * =
     * \frac{1}{i\pi^2}\int d^4k\,
     * \frac{1}{(k^2-m_0^2+i\varepsilon)
     * \left[(k+p_1)^2-m_1^2+i\varepsilon\right]
     * \left[(k+p_1+p_2)^2-m_2^2+i\varepsilon\right]
     * \left[(k+p_1+p_2+p_3)^2-m_2^2+i\varepsilon\right]}\,,
     * @f}
     * where @f$p_1+p_2+p_3+p_4=0@f$. 
     * The current functions handles only the special case of 
     * @f$p_1^2=p_2^2=p_3^2=p_4^2=0@f$.
     * @param[in] s,t momentum squared, @f$s@f$ and @f$t@f$
     * @param[in] m02, m12, m22, m32 mass squared, @f$m_0^2@f$, @f$m_1^2@f$, @f$m_2^2@f$ and @f$m_3^2@f$
     * @return @f$D_{0}(0,0,0,0,s,t; m_0^2, m_1^2, m_2^2, m_3^2)@f$
     *
     * @warning Only the case of @f$s=t=0@f$ has been implemented. Other cases
     * can be computed with the help of LoopTools library, by setting the
     * preprocessor macro USE_LOOPTOOLS.
     */
    gslpp::complex D0(const double s, const double t, const double m02, const double m12,
               const double m22, const double m32) const;

    /**
     * @brief @f$D_{00}(0,0,0,0,s,t; m_0^2, m_1^2, m_2^2, m_3^2)@f$.
     * @details The tensor four-point function
     * @f$D_{0}(p_1^2,p_2^2,p_3^2,p_4^2,(p_1+p_2)^2,(p_2+p_3)^2; m_0^2, m_1^2, m_2^2, m_3^2)@f$
     * is defined as
     * @f{eqnarray*}{
     * &&g_{\mu\nu}
     * D_{00}(p_1^2,p_2^2,p_3^2,p_4^2,(p_1+p_2)^2,(p_2+p_3)^2; m_0^2,m_1^2,m_2^2,m_3^2)
     * + \sum_{i,j=1}^{3}q_{i\mu}q_{j\mu}
     * D_{ij}(p_1^2,p_2^2,p_3^2,p_4^2,(p_1+p_2)^2,(p_2+p_3)^2; m_0^2,m_1^2,m_2^2,m_3^2)
     * \\
     * &&\quad
     * =
     * \frac{1}{i\pi^2}\int d^4k\,
     * \frac{k_\mu k_\nu}{(k^2-m_0^2+i\varepsilon)
     * \left[(k+p_1)^2-m_1^2+i\varepsilon\right]
     * \left[(k+p_1+p_2)^2-m_2^2+i\varepsilon\right]
     * \left[(k+p_1+p_2+p_3)^2-m_2^2+i\varepsilon\right]}\,,
     * @f}
     * where @f$q_N=\sum_{i=1}^N p_i@f$ and @f$p_1+p_2+p_3+p_4=0@f$.
     * The current functions handles only the special case of
     * @f$p_1^2=p_2^2=p_3^2=p_4^2=0@f$.
     * @param[in] s,t momentum squared, @f$s@f$ and @f$t@f$
     * @param[in] m02, m12, m22, m32 mass squared, @f$m_0^2@f$, @f$m_1^2@f$, @f$m_2^2@f$ and @f$m_3^2@f$
     * @return @f$D_{00}(0,0,0,0,s,t; m_0^2, m_1^2, m_2^2, m_3^2)@f$
     *
     * @warning Only the case of @f$s=t=0@f$ has been implemented. Other cases
     * can be computed with the help of LoopTools library, by setting the
     * preprocessor macro USE_LOOPTOOLS.
     */
    gslpp::complex D00(const double s, const double t, const double m02, const double m12,
                const double m22, const double m32) const;

private:
    double ExtraMinusSign;///< An overall factor for the one-point and three-point functions, initialized in PVfunctions().
    Polylogarithms myPolylog;///< An object of type Polylogarithms.
#ifdef USE_LOOPTOOLS
    LoopToolsWrapper myLT;///< An object of type LoopToolsWrapper.
#endif
};

#endif	/* PVFUNCTIONS_H */

