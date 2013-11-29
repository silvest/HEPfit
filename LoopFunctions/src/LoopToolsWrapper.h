/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LOOPTOOLSWRAPPER_H
#define	LOOPTOOLSWRAPPER_H

/**
 * @class LoopToolsWrapper
 * @ingroup LoopFunctions 
 * @brief A wrapper class for LoopTools library.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
#include <gslpp.h>
using namespace gslpp;
   
/**
 * @class LoopToolsWrapper
 * @brief C++ interface class for LoopTools library. 
 */
class LoopToolsWrapper {
public:

    LoopToolsWrapper();

    virtual ~LoopToolsWrapper();
    
    /**
     * @brief The scalar one-point Passarino-Veltman function.
     * @param[in] mu2 The renormalization scale squared.
     * @param[in] m2 Mass squared.
     * @return The finite part of the scalar one-point PV function.
     */
    double PV_A0(const double mu2, const double m2) const;
    
    /**
     * @brief The scalar two-point Passarino-Veltman function.
     * @param[in] mu2 The renormalization scale squared.
     * @param[in] p2 Momentum squared.
     * @param[in] m02, m12 Mass squared.
     * @return The finite part of the scalar two-point PV function.
     */
    complex PV_B0(const double mu2, const double p2,
                  const double m02, const double m12) const;
    
    /**
     * @brief The vector two-point Passarino-Veltman function,
     * @param[in] mu2 The renormalization scale squared.
     * @param[in] p2 Momentum squared.
     * @param[in] m02, m12 Mass squared.
     * @return The finite part of the vector two-point PV function.
     */
    complex PV_B1(const double mu2, const double p2,
                  const double m02, const double m12) const;

    /**
     * @brief A tensor two-point Passarino-Veltman function.
     * @param[in] mu2 The renormalization scale squared.
     * @param[in] p2 Momentum squared.
     * @param[in] m02, m12 Mass squared.
     * @return The finite part of a tensor two-point PV function B_{21}.
     */
    complex PV_B21(const double mu2, const double p2,
                   const double m02, const double m12) const;

    /**
     * @brief A tensor two-point Passarino-Veltman function.
     * @param[in] mu2 The renormalization scale squared.
     * @param[in] p2 Momentum squared.
     * @param[in] m02, m12 Mass squared.
     * @return The finite part of a tensor two-point PV function B_{22}.
     */
    complex PV_B22(const double mu2, const double p2,
                   const double m02, const double m12) const;

    /**
     * @brief The derivative of B_0.
     * @param[in] muIR2 The IR renormalization scale squared.
     * @param[in] p2 Momentum squared.
     * @param[in] m02, m12 Mass squared.
     * @return The finite part of B_{0p}.
     */
    complex PV_B0p(const double muIR2, const double p2,
                   const double m02, const double m12) const;
    
    /**
     * @brief The derivative of B_1.
     * @param[in] mu2 The renormalization scale squared.
     * @param[in] p2 Momentum squared.
     * @param[in] m02, m12 Mass squared.
     * @return The finite part of B_{1p}.
     */
    complex PV_B1p(const double mu2, const double p2,
                   const double m02, const double m12) const;
    
    /**
     * @brief The derivative of B_{21}.
     * @param[in] mu2 The renormalization scale squared.
     * @param[in] p2 Momentum squared.
     * @param[in] m02, m12 Mass squared.
     * @return The finite part of B_{21p}.
     */
    complex PV_B21p(const double mu2, const double p2,
                    const double m02, const double m12) const;

    /**
     * @brief The derivative of B_{22}.
     * @param[in] mu2 The renormalization scale squared.
     * @param[in] p2 Momentum squared.
     * @param[in] m02, m12 Mass squared.
     * @return The finite part of B_{22p}.
     */
    complex PV_B22p(const double mu2, const double p2,
                    const double m02, const double m12) const;

    /**
     * @brief The scalar three-point Passarino-Veltman function C_0(0,0,p2;m02,m12,m22).
     * @param[in] p2 Momentum squared.
     * @param[in] m02, m12, m22 Mass squared.
     * @return The scalar three-point PV function C_0(0,0,p2;m02,m12,m22).
     */
    complex PV_C0(const double p2, 
                  const double m02, const double m12, const double m22) const;
    
    /**
     * @brief The scalar four-point Passarino-Veltman function D_0(0,0,0,0,s,t;m02,m12,m22,m32).
     * @param[in] s, t Momentum squared.
     * @param[in] m02, m12, m22, m32 Mass squared.
     * @return The scalar four-point PV function D_0(0,0,0,0,s,t;m02,m12,m22,m32).
     */
    complex PV_D0(const double s, const double t, const double m02, const double m12,
                  const double m22, const double m32) const;
    
    /**
     * @brief A tensor four-point Passarino-Veltman function D_22(0,0,0,0,s,t;m02,m12,m22,m32).
     * @param[in] s, t Momentum squared.
     * @param[in] m02, m12, m22, m32 Mass squared.
     * @return A tensor four-point PV function D_0(0,0,0,0,s,t;m02,m12,m22,m32).
     */
    complex PV_D22(const double s, const double t, const double m02, const double m12,
                   const double m22, const double m32) const;

private:

};

#endif	/* LOOPTOOLSWRAPPER_H */

