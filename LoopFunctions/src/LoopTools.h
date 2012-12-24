/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LOOPTOOLS_H
#define	LOOPTOOLS_H

#include <gslpp.h>
using namespace gslpp;
   
/**
 * @class LoopTools
 * @brief C++ interface class for LoopTools library
 */
class LoopTools {
public:

    LoopTools();

    virtual ~LoopTools();
    
    /**
     * @brief the scalar one-point Passarino-Veltman function
     * @param[in] mu renormalization scale
     * @param[in] m mass m
     * @return the finite part of the scalar one-point PV function at the scale mu
     */
    double PV_A0(const double mu, const double m) const;
    
    /**
     * @brief the scalar two-point Passarino-Veltman function
     * @param[in] mu renormalization scale
     * @param[in] p2 p squared
     * @param[in] m0 mass m_0
     * @param[in] m1 mass m_1
     * @return the finite part of the scalar two-point PV function at the scale mu
     */
    complex PV_B0(const double mu, const double p2, 
                  const double m0, const double m1) const; 
    
    /**
     * @brief the scalar three-point Passarino-Veltman function C_0(0,0,p2;m0,m1,m2)
     * @param[in] p2 p squared
     * @param[in] m0 mass m_0
     * @param[in] m1 mass m_1
     * @param[in] m2 mass m_2
     * @return the scalar three-point PV function C_0(0,0,p2;m0,m1,m2)
     */
    complex PV_C0(const double p2, 
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
    complex PV_D0(const double s, const double t, const double m0, const double m1, 
                  const double m2, const double m3) const;    
    
private:

};

#endif	/* LOOPTOOLS_H */

