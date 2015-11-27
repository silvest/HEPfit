/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LOOPTOOLSWRAPPER_H
#define	LOOPTOOLSWRAPPER_H

#include "gslpp.h"
   
/**
 * @class LoopToolsWrapper
 * @ingroup LoopFunctions
 * @brief A wrapper class for LoopTools library.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is responsible for the interface to 
 * <a href="http://www.feynarts.de/looptools/" target=blank>LoopTools</a> library
 * @cite Hahn:1998yk.
 *
 * See @cite Hahn:1998yk and
 * <a href="http://www.feynarts.de/looptools/" target=blank>the webpage of LoopTools</a>.
 *
 * See also the documentation of PVfunctions class for the definitions of the
 * loop functions. 
 *
 * @attention The preprocessor macro USE_LOOPTOOLS has to be set in order to
 * use the current wrapper class. See also PVfunctions.h.
 */
class LoopToolsWrapper {
public:

    /**
     * @brief The default constructor.
     */
    LoopToolsWrapper();

    /**
     * @brief The default destructor.
     */
    virtual ~LoopToolsWrapper();
    
#ifdef USE_LOOPTOOLS

    /**
     * @brief @f$A_0(m^2)@f$.
     * @param[in] mu2 the renormalization scale squared, @f$\mu^2@f$
     * @param[in] m2 mass squared, @f$m^2@f$
     * @return the finite part of @f$A_0(m^2)@f$
     * in the sense of the @f$\overline{\mathrm{MS}}@f$ scheme
     */
    double PV_A0(const double mu2, const double m2) const;
    
    /**
     * @brief @f$B_0(p^2; m_0^2, m_1^2)@f$.
     * @param[in] mu2 the renormalization scale squared, @f$\mu^2@f$
     * @param[in] p2 momentum squared, @f$p^2@f$
     * @param[in] m02, m12 mass squared, @f$m_0^2@f$ and @f$m_1^2@f$
     * @return the finite part of @f$B_0(p^2; m_0^2, m_1^2)@f$
     * in the sense of the @f$\overline{\mathrm{MS}}@f$ scheme
     */
    gslpp::complex PV_B0(const double mu2, const double p2,
                  const double m02, const double m12) const;
    
    /**
     * @brief @f$B_1(p^2; m_0^2, m_1^2)@f$.
     * @param[in] mu2 the renormalization scale squared, @f$\mu^2@f$
     * @param[in] p2 momentum squared, @f$p^2@f$
     * @param[in] m02, m12 mass squared, @f$m_0^2@f$ and @f$m_1^2@f$
     * @return the finite part of @f$B_1(p^2; m_0^2, m_1^2)@f$
     * in the sense of the @f$\overline{\mathrm{MS}}@f$ scheme
     */
    gslpp::complex PV_B1(const double mu2, const double p2,
                  const double m02, const double m12) const;

    /**
     * @brief @f$B_{11}(p^2; m_0^2, m_1^2)@f$.
     * @param[in] mu2 the renormalization scale squared, @f$\mu^2@f$
     * @param[in] p2 momentum squared, @f$p^2@f$
     * @param[in] m02, m12 mass squared, @f$m_0^2@f$ and @f$m_1^2@f$
     * @return the finite part of @f$B_{11}(p^2; m_0^2, m_1^2)@f$
     * in the sense of the @f$\overline{\mathrm{MS}}@f$ scheme
     */
    gslpp::complex PV_B11(const double mu2, const double p2,
                   const double m02, const double m12) const;

    /**
     * @brief @f$B_{00}(p^2; m_0^2, m_1^2)@f$.
     * @param[in] mu2 the renormalization scale squared, @f$\mu^2@f$
     * @param[in] p2 momentum squared, @f$p^2@f$
     * @param[in] m02, m12 mass squared, @f$m_0^2@f$ and @f$m_1^2@f$
     * @return the finite part of @f$B_{00}(p^2; m_0^2, m_1^2)@f$
     * in the sense of the @f$\overline{\mathrm{MS}}@f$ scheme
     */
    gslpp::complex PV_B00(const double mu2, const double p2,
                   const double m02, const double m12) const;

    /**
     * @brief @f$B_{0p}(p^2; m_0^2, m_1^2)@f$.
     * @param[in] muIR2 the renormalization scale squared for the IR divergence, @f$\mu_{\mathrm{IR}}^2@f$
     * @param[in] p2 momentum squared, @f$p^2@f$
     * @param[in] m02, m12 mass squared, @f$m_0^2@f$ and @f$m_1^2@f$
     * @return the finite part of @f$B_{0p}(p^2; m_0^2, m_1^2)@f$
     * in the sense of the @f$\overline{\mathrm{MS}}@f$ scheme
     */
    gslpp::complex PV_B0p(const double muIR2, const double p2,
                   const double m02, const double m12) const;
    
    /**
     * @brief @f$B_{1p}(p^2; m_0^2, m_1^2)@f$.
     * @param[in] mu2 the renormalization scale squared, @f$\mu^2@f$
     * @param[in] p2 momentum squared, @f$p^2@f$
     * @param[in] m02, m12 mass squared, @f$m_0^2@f$ and @f$m_1^2@f$
     * @return the finite part of @f$B_{1p}(p^2; m_0^2, m_1^2)@f$
     * in the sense of the @f$\overline{\mathrm{MS}}@f$ scheme
     */
    gslpp::complex PV_B1p(const double mu2, const double p2,
                   const double m02, const double m12) const;
    
    /**
     * @brief @f$B_{11p}(p^2; m_0^2, m_1^2)@f$.
     * @param[in] mu2 the renormalization scale squared, @f$\mu^2@f$
     * @param[in] p2 momentum squared, @f$p^2@f$
     * @param[in] m02, m12 mass squared, @f$m_0^2@f$ and @f$m_1^2@f$
     * @return the finite part of @f$B_{11p}(p^2; m_0^2, m_1^2)@f$
     * in the sense of the @f$\overline{\mathrm{MS}}@f$ scheme
     */
    gslpp::complex PV_B11p(const double mu2, const double p2,
                    const double m02, const double m12) const;

    /**
     * @brief @f$B_{00p}(p^2; m_0^2, m_1^2)@f$.
     * @param[in] mu2 the renormalization scale squared, @f$\mu^2@f$
     * @param[in] p2 momentum squared, @f$p^2@f$
     * @param[in] m02, m12 mass squared, @f$m_0^2@f$ and @f$m_1^2@f$
     * @return the finite part of @f$B_{00p}(p^2; m_0^2, m_1^2)@f$
     * in the sense of the @f$\overline{\mathrm{MS}}@f$ scheme
     */
    gslpp::complex PV_B00p(const double mu2, const double p2,
                    const double m02, const double m12) const;

    /**
     * @brief @f$C_{0}(0,0,p^2; m_0^2, m_1^2, m_2^2)@f$.
     * @param[in] p2 momentum squared, @f$p^2@f$
     * @param[in] m02, m12, m22 mass squared, @f$m_0^2@f$, @f$m_1^2@f$ and @f$m_2^2@f$
     * @return @f$C_{0}(0,0,p^2; m_0^2, m_1^2, m_2^2)@f$
     */
    gslpp::complex PV_C0(const double p2, 
                  const double m02, const double m12, const double m22) const;
    
    /**
     * @brief @f$D_{0}(0,0,0,0,s,t; m_0^2, m_1^2, m_2^2, m_3^2)@f$.
     * @param[in] s,t momentum squared, @f$s@f$ and @f$t@f$
     * @param[in] m02, m12, m22, m32 mass squared, @f$m_0^2@f$, @f$m_1^2@f$, @f$m_2^2@f$ and @f$m_3^2@f$
     * @return @f$D_{0}(0,0,0,0,s,t; m_0^2, m_1^2, m_2^2, m_3^2)@f$
     */
    gslpp::complex PV_D0(const double s, const double t, const double m02, const double m12,
                  const double m22, const double m32) const;
    
    /**
     * @brief @f$D_{00}(0,0,0,0,s,t; m_0^2, m_1^2, m_2^2, m_3^2)@f$.
     * @param[in] s,t momentum squared, @f$s@f$ and @f$t@f$
     * @param[in] m02, m12, m22, m32 mass squared, @f$m_0^2@f$, @f$m_1^2@f$, @f$m_2^2@f$ and @f$m_3^2@f$
     * @return @f$D_{00}(0,0,0,0,s,t; m_0^2, m_1^2, m_2^2, m_3^2)@f$
     *
     * @warning This function does not work.
     */
    gslpp::complex PV_D00(const double s, const double t, const double m02, const double m12,
                   const double m22, const double m32) const;

#endif

};

#endif	/* LOOPTOOLSWRAPPER_H */

