/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef PVFUNCTIONS_H
#define	PVFUNCTIONS_H

#include <gslpp_complex.h>
#include "Polylogarithms.h"
#include "LoopToolsWrapper.h"

using namespace gslpp;

// set in case where LoopTools library is employed. 
//#define USE_LOOPTOOLS


/**
 * @class PVfunctions
 * @ingroup LoopFunctions 
 * @brief A class for Passarino-Veltman functions.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details The definitions of the two-point and four-point functions are identical 
 * to those in LoopTools library. The functions A_0 and C_0 are identical to those
 * in LoopTools library when the argument passes to the constructor is "false",
 * while an extra minus sign is added to them when the argument is "true".
 * If the preprocessor macro "USE_LOOPTOOLS" is defined in PVdunctions.h or 
 * Makefile, the functions in LoopTools library, called via LoopToolsWrapper class,
 * are employed instead of those defined in the current class.
 */
class PVfunctions {
public:

    /**
     * @brief PVfunctions constructor. 
     * @param[in] bExtraMinusSign true if adding an overall extra minus sign
     * to A_0 and C_0, compared to the ones in LoopTools library.
     */
    PVfunctions(const bool bExtraMinusSign);

    /**
     * @brief The scalar one-point Passarino-Veltman function. 
     * @param[in] mu2 The renormalization scale squared. 
     * @param[in] m2 Mass squared.
     * @return The finite part of the scalar one-point PV function.
     */
    double A0(const double mu2, const double m2) const;
    
    /**
     * @brief The scalar two-point Passarino-Veltman function.
     * @param[in] mu2 The renormalization scale squared.
     * @param[in] p2 Momentum squared.
     * @param[in] m02, m12 Mass squared.
     * @return The finite part of the scalar two-point PV function.
     */
    complex B0(const double mu2, const double p2,
               const double m02, const double m12) const;
    
    /**
     * @brief The vector two-point Passarino-Veltman function,
     * @param[in] mu2 The renormalization scale squared.
     * @param[in] p2 Momentum squared.
     * @param[in] m02, m12 Mass squared.
     * @return The finite part of the vector two-point PV function.
     */
    complex B1(const double mu2, const double p2,
               const double m02, const double m12) const;
    
    /**
     * @brief A tensor two-point Passarino-Veltman function.
     * @param[in] mu2 The renormalization scale squared.
     * @param[in] p2 Momentum squared.
     * @param[in] m02, m12 Mass squared.
     * @return The finite part of a tensor two-point PV function B_{11}.
     */
    complex B11(const double mu2, const double p2,
                const double m02, const double m12) const;
    
    /**
     * @brief A tensor two-point Passarino-Veltman function.
     * @param[in] mu2 The renormalization scale squared.
     * @param[in] p2 Momentum squared.
     * @param[in] m02, m12 Mass squared.
     * @return The finite part of a tensor two-point PV function B_{00}.
     */
    complex B00(const double mu2, const double p2,
                const double m02, const double m12) const;
    
    /**
     * @brief A sum of two-point Passarino-Veltman functions
     * @param[in] mu2 The renormalization scale squared.
     * @param[in] p2 Momentum squared.
     * @param[in] m02, m12 Mass squared.
     * @return The finite part of a sum of two-point PV function B_f.
     */
    complex Bf(const double mu2, const double p2,
               const double m02, const double m12) const;
    
    /**
     * @brief The derivative of B_0.
     * @param[in] muIR2 The IR renormalization scale squared. 
     * @param[in] p2 Momentum squared.
     * @param[in] m02, m12 Mass squared.
     * @return The finite part of B_{0p}.
     */
    complex B0p(const double muIR2, const double p2,
                const double m02, const double m12) const;
    
    /**
     * @brief The derivative of B_1.
     * @param[in] mu2 The renormalization scale squared.
     * @param[in] p2 Momentum squared.
     * @param[in] m02, m12 Mass squared.
     * @return The finite part of B_{1p}. 
     */
    complex B1p(const double mu2, const double p2,
                const double m02, const double m12) const;
    
    /**
     * @brief The derivative of B_{11}.
     * @param[in] mu2 The renormalization scale squared.
     * @param[in] p2 Momentum squared.
     * @param[in] m02, m12 Mass squared.
     * @return The finite part of B_{11p}.
     */
    complex B11p(const double mu2, const double p2,
                 const double m02, const double m12) const;

    /**
     * @brief The derivative of B_{00}.
     * @param[in] mu2 The renormalization scale squared.
     * @param[in] p2 Momentum squared.
     * @param[in] m02, m12 Mass squared.
     * @return The finite part of B_{00p}.
     */
    complex B00p(const double mu2, const double p2,
                 const double m02, const double m12) const;
    
    /**
     * @brief The derivative of B_{f}.
     * @param[in] mu2 The renormalization scale squared.
     * @param[in] p2 Momentum squared.
     * @param[in] m02, m12 Mass squared.
     * @return The finite part of B_{fp}. 
     */
    complex Bfp(const double mu2, const double p2,
                const double m02, const double m12) const;
    
    /**
     * @brief The scalar three-point Passarino-Veltman function C_0(0,0,p2;m02,m12,m22).
     * @param[in] p2 Momentum squared.
     * @param[in] m02, m12, m22 Mass squared.
     * @return The scalar three-point PV function C_0(0,0,p2;m02,m12,m22).
     */
    complex C0(const double p2, 
               const double m02, const double m12, const double m22) const;
    
    /**
     * @brief The scalar four-point Passarino-Veltman function D_0(0,0,0,0,s,t;m02,m12,m22,m32).
     * @param[in] s, t Momentum squared.
     * @param[in] m02, m12, m22, m32 Mass squared.
     * @return The scalar four-point PV function D_0(0,0,0,0,s,t;m02,m12,m22,m32).
     */
    complex D0(const double s, const double t, const double m02, const double m12,
               const double m22, const double m32) const;

    /**
     * @brief A tensor four-point Passarino-Veltman function D_00(0,0,0,0,s,t;m02,m12,m22,m32).
     * @param[in] s, t Momentum squared.
     * @param[in] m02, m12, m22, m32 Mass squared.
     * @return A tensor four-point PV function D_00(0,0,0,0,s,t;m02,m12,m22,m32).
     */
    complex D00(const double s, const double t, const double m02, const double m12,
                const double m22, const double m32) const;

private:
    double ExtraMinusSign;
    Polylogarithms myPolylog;
    LoopToolsWrapper myLT;

};

#endif	/* PVFUNCTIONS_H */

