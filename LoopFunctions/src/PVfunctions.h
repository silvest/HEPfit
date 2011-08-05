/* 
 * File:   PVfunctions.h
 * Author: mishima
 *
 * Created on June 10, 2011, 2:47 PM
 */

#ifndef PVFUNCTIONS_H
#define	PVFUNCTIONS_H

#include <gslpp_complex.h>

using namespace gslpp;

class PVfunctions {
public:

    PVfunctions();

    PVfunctions(const PVfunctions& orig);

    virtual ~PVfunctions();

    /**
     * @brief the scalar one-point Passarino-Veltman function
     * @param[in] mu renormalization scale
     * @param[in] m mass m
     * @return the finite part of the scalar one-point PV function at the scale mu
     */
    double A0(const double mu, const double m);
    
    /**
     * @brief the scalar two-point Passarino-Veltman function
     * @param[in] mu renormalization scale
     * @param[in] p2 p squared
     * @param[in] m0 mass m_0
     * @param[in] m1 mass m_1
     * @return the finite part of the scalar two-point PV function at the scale mu
     */
    complex B0(const double mu, const double p2, 
               const double m0, const double m1); 
    
    /**
     * @brief the scalar three-point Passarino-Veltman function
     * @param[in] p2 p squared
     * @param[in] m0 mass m_0
     * @param[in] m1 mass m_1
     * @param[in] m2 mass m_2
     * @return the scalar three-point PV function 
     */
    complex C0(const double p2, 
               const double m0, const double m1, const double m2); 
    
private:

};

#endif	/* PVFUNCTIONS_H */

