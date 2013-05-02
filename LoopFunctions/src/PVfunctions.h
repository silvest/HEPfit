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
#include "LoopTools.h"

using namespace gslpp;

// set in case where LoopTools library is employed. 
//#define USE_LOOPTOOLS


/**
 * @class PVfunctions
 * @ingroup LoopFunctions 
 * @brief A class for Passarino-Veltman function. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class PVfunctions {
public:

    PVfunctions() 
    {};

    /**
     * @brief the scalar one-point Passarino-Veltman function
     * @param[in] mu renormalization scale
     * @param[in] m mass m
     * @return the finite part of the scalar one-point PV function at the scale mu
     */
    double A0(const double mu, const double m) const;
    
    /**
     * @brief the scalar two-point Passarino-Veltman function
     * @param[in] mu renormalization scale
     * @param[in] p2 p squared
     * @param[in] m0 mass m_0
     * @param[in] m1 mass m_1
     * @return the finite part of the scalar two-point PV function at the scale mu
     */
    complex B0(const double mu, const double p2, 
               const double m0, const double m1) const; 
    
    /**
     * @brief the vector two-point Passarino-Veltman function
     * @param[in] mu renormalization scale
     * @param[in] p2 p squared
     * @param[in] m0 mass m_0
     * @param[in] m1 mass m_1
     * @return the finite part of the vector two-point PV function at the scale mu
     */
    complex B1(const double mu, const double p2, 
               const double m0, const double m1) const; 
    
    /**
     * @brief a tensor two-point Passarino-Veltman function
     * @param[in] mu renormalization scale
     * @param[in] p2 p squared
     * @param[in] m0 mass m_0
     * @param[in] m1 mass m_1
     * @return the finite part of a tensor two-point PV function B_{21} at the scale mu
     */
    complex B21(const double mu, const double p2, 
                const double m0, const double m1) const;     
    
    /**
     * @brief a tensor two-point Passarino-Veltman function
     * @param[in] mu renormalization scale
     * @param[in] p2 p squared
     * @param[in] m0 mass m_0
     * @param[in] m1 mass m_1
     * @return the finite part of a tensor two-point PV function B_{22} at the scale mu
     */
    complex B22(const double mu, const double p2, 
                const double m0, const double m1) const;       
    
    /**
     * @brief a sum of two-point Passarino-Veltman functions
     * @param[in] mu renormalization scale
     * @param[in] p2 p squared
     * @param[in] m0 mass m_0
     * @param[in] m1 mass m_1
     * @return the finite part of a sum of two-point PV function at the scale mu
     */
    complex Bf(const double mu, const double p2, 
               const double m0, const double m1) const;       
    
    /**
     * @brief the derivative of B_0
     * @param[in] muIR renormalization scale for an IR divergence
     * @param[in] p2 p squared
     * @param[in] m0 mass m_0
     * @param[in] m1 mass m_1
     * @return the finite part of B_{0p}
     */
    complex B0p(const double muIR, const double p2, 
                const double m0, const double m1) const; 
    
    /**
     * @brief the derivative of B_1
     * @param[in] mu renormalization scale
     * @param[in] p2 p squared
     * @param[in] m0 mass m_0
     * @param[in] m1 mass m_1
     * @return the finite part of B_{1p}
     */
    complex B1p(const double mu, const double p2, 
                const double m0, const double m1) const; 
    
    /**
     * @brief the derivative of B_{21}
     * @param[in] mu renormalization scale
     * @param[in] p2 p squared
     * @param[in] m0 mass m_0
     * @param[in] m1 mass m_1
     * @return the finite part of B_{21p}
     */
    complex B21p(const double mu, const double p2, 
                 const double m0, const double m1) const;     
    
    /**
     * @brief the derivative of B_{22}
     * @param[in] mu renormalization scale
     * @param[in] p2 p squared
     * @param[in] m0 mass m_0
     * @param[in] m1 mass m_1
     * @return the finite part of B_{22p}
     */
    //complex B22p(const double mu, const double p2, 
    //             const double m0, const double m1);       
    
    /**
     * @brief the derivative of B_{f}
     * @param[in] mu renormalization scale
     * @param[in] p2 p squared
     * @param[in] m0 mass m_0
     * @param[in] m1 mass m_1
     * @return the finite part of B_{fp}
     */
    complex Bfp(const double mu, const double p2, 
                const double m0, const double m1) const;       
    
    /**
     * @brief the scalar three-point Passarino-Veltman function C_0(0,0,p2;m0,m1,m2)
     * @param[in] p2 p squared
     * @param[in] m0 mass m_0
     * @param[in] m1 mass m_1
     * @param[in] m2 mass m_2
     * @return the scalar three-point PV function C_0(0,0,p2;m0,m1,m2)
     */
    complex C0(const double p2, 
               const double m0, const double m1, const double m2) const; 
    
    /**
     * @brief the scalar four-point Passarino-Veltman function D_0(0,0,0,0,s,t;m0,m1,m2,m3)
     * @param[in] s 
     * @param[in] t 
     * @param[in] m0 mass m_0
     * @param[in] m1 mass m_1
     * @param[in] m2 mass m_2
     * @param[in] m3 mass m_3
     * @return the scalar four-point PV function D_0(0,0,0,0,s,t;m0,m1,m2,m3)
     */
    complex D0(const double s, const double t, const double m0, const double m1, 
               const double m2, const double m3) const;

private:
    Polylogarithms myPolylog;
#ifdef USE_LOOPTOOLS
    LoopTools myLT;
#endif
    
};

#endif	/* PVFUNCTIONS_H */

