/* 
 * File:   PVfunctions.h
 * Author: mishima
 */

#ifndef PVFUNCTIONS_H
#define	PVFUNCTIONS_H

#include <gslpp_complex.h>

using namespace gslpp;

class PVfunctions {
public:

    PVfunctions();

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
     * @param[in] mu renormalization scale
     * @param[in] p2 p squared
     * @param[in] m0 mass m_0
     * @param[in] m1 mass m_1
     * @return the finite part of B_{0p}
     */
    complex B0p(const double mu, const double p2, 
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
     * @brief the scalar three-point Passarino-Veltman function
     * @param[in] p2 p squared
     * @param[in] m0 mass m_0
     * @param[in] m1 mass m_1
     * @param[in] m2 mass m_2
     * @return the scalar three-point PV function 
     */
    complex C0(const double p2, 
               const double m0, const double m1, const double m2) const; 
    
    /**
     * @brief the scalar three-point Passarino-Veltman function C_11
     * @param[in] mu the renrmalization scale
     * @param[in] p12 p1 squared
     * @param[in] p22 p2 squared
     * @param[in] s (p1+p2) squared
     * @param[in] m1 mass m_1
     * @param[in] m2 mass m_2
     * @param[in] m3 mass m_3
     * @return the scalar three-point PV function 
     */
    complex C11(const double mu,const double p12, const double p22, const double s,
               const double m1, const double m2, const double m3) const; 
    
    /**
     * @brief the scalar three-point Passarino-Veltman function C_12
     * @param[in] s momentum squared
     * @param[in] m1 mass m_1
     * @param[in] m2 mass m_2
     * @param[in] m3 mass m_3
     * @return the scalar three-point PV function 
     */
    complex C12(const double mu,const double p12, const double p22, const double s,
               const double m1, const double m2, const double m3) const;
    
    
    /**
     * @brief function F(m0,m1) used for THDM. Remember that this function is
     * defined for THDM while for SUSY we have a multiplicative factor 2.
     * @param[in] m0 mass m_0
     * @param[in] m1 mass m_1
     * @return the function F for THDM 
     */
    double F(const double m0, const double m1) const;
    
    
    /**
     * @brief the scalar four-point Passarino-Veltman function
     * @param[in] p12 p1 squared
     * @param[in] p22 p2 squared
     * @param[in] p32 p3 squared
     * @param[in] p42 p4 squared
     * @param[in] m1 mass m_1
     * @param[in] m2 mass m_2
     * @param[in] m3 mass m_3
     * @param[in] m4 mass m_4
     * @return the function D0
     */
    complex D0(const double p12,const double p22,const double p32,const double p42,
               const double s, const double t,const double m1, const double m2, 
               const double m3, const double m4) const;
    
    
    /**
     * @brief the scalar four-point Passarino-Veltman function D11
     * @param[in] p12 p1 squared
     * @param[in] p22 p2 squared
     * @param[in] p32 p3 squared
     * @param[in] p42 p4 squared
     * @param[in] m1 mass m_1
     * @param[in] m2 mass m_2
     * @param[in] m3 mass m_3
     * @param[in] m4 mass m_4
     * @return the function D11
     */
    complex D11(const double p12, const double p22,
                         const double p32,const double p42, const double s,
                         const double m1, const double m2, const double m3, 
                         const double m4, const double theta) const;
    
    /**
     * @brief the scalar four-point Passarino-Veltman function D11
     * @param[in] p12 p1 squared
     * @param[in] p22 p2 squared
     * @param[in] p32 p3 squared
     * @param[in] p42 p4 squared
     * @param[in] m1 mass m_1
     * @param[in] m2 mass m_2
     * @param[in] m3 mass m_3
     * @param[in] m4 mass m_4
     * @return the function D11
     */
    complex D12(const double p12, const double p22,
                         const double p32,const double p42, const double s,
                         const double m1, const double m2, const double m3, 
                         const double m4, const double theta) const;
    
    /**
     * @brief the scalar four-point Passarino-Veltman function D11
     * @param[in] p12 p1 squared
     * @param[in] p22 p2 squared
     * @param[in] p32 p3 squared
     * @param[in] p42 p4 squared
     * @param[in] m1 mass m_1
     * @param[in] m2 mass m_2
     * @param[in] m3 mass m_3
     * @param[in] m4 mass m_4
     * @return the function D11
     */
    complex D13(const double p12, const double p22,
                         const double p32,const double p42, const double s,
                         const double m1, const double m2, const double m3, 
                         const double m4, const double theta) const;
    
    /**
     * @brief the scalar four-point Passarino-Veltman function D11
     * @param[in] p12 p1 squared
     * @param[in] p22 p2 squared
     * @param[in] p32 p3 squared
     * @param[in] p42 p4 squared
     * @param[in] m1 mass m_1
     * @param[in] m2 mass m_2
     * @param[in] m3 mass m_3
     * @param[in] m4 mass m_4
     * @return the function D11
     */
    complex D27(const double p12, const double p22,
                         const double p32,const double p42, const double s,
                         const double m1, const double m2, const double m3, 
                         const double m4, const double theta) const;
    
    /**
     * @brief the scalar four-point Passarino-Veltman function D11
     * @param[in] p12 p1 squared
     * @param[in] p22 p2 squared
     * @param[in] p32 p3 squared
     * @param[in] p42 p4 squared
     * @param[in] m1 mass m_1
     * @param[in] m2 mass m_2
     * @param[in] m3 mass m_3
     * @param[in] m4 mass m_4
     * @return the function D11
     */
    complex D21(const double mu,const double p12, const double p22,
                         const double p32,const double p42, const double s,
                         const double m1, const double m2, const double m3, 
                         const double m4, const double theta) const;
    
    /**
     * @brief the scalar four-point Passarino-Veltman function D11
     * @param[in] p12 p1 squared
     * @param[in] p22 p2 squared
     * @param[in] p32 p3 squared
     * @param[in] p42 p4 squared
     * @param[in] m1 mass m_1
     * @param[in] m2 mass m_2
     * @param[in] m3 mass m_3
     * @param[in] m4 mass m_4
     * @return the function D11
     */
    complex D24(const double mu,const double p12, const double p22,
                         const double p32,const double p42, const double s,
                         const double m1, const double m2, const double m3, 
                         const double m4, const double theta) const;
    
    /**
     * @brief the scalar four-point Passarino-Veltman function D11
     * @param[in] p12 p1 squared
     * @param[in] p22 p2 squared
     * @param[in] p32 p3 squared
     * @param[in] p42 p4 squared
     * @param[in] m1 mass m_1
     * @param[in] m2 mass m_2
     * @param[in] m3 mass m_3
     * @param[in] m4 mass m_4
     * @return the function D11
     */
    complex D25(const double mu,const double p12, const double p22,
                         const double p32,const double p42, const double s,
                         const double m1, const double m2, const double m3, 
                         const double m4, const double theta) const;
    
    
    /**
     * @brief the scalar four-point Passarino-Veltman function D11
     * @param[in] p12 p1 squared
     * @param[in] p22 p2 squared
     * @param[in] p32 p3 squared
     * @param[in] p42 p4 squared
     * @param[in] m1 mass m_1
     * @param[in] m2 mass m_2
     * @param[in] m3 mass m_3
     * @param[in] m4 mass m_4
     * @return the function D11
     */
    complex D22(const double mu,const double p12, const double p22,
                         const double p32,const double p42, const double s,
                         const double m1, const double m2, const double m3, 
                         const double m4, const double theta) const;
    
    
    /**
     * @brief the scalar four-point Passarino-Veltman function D11
     * @param[in] p12 p1 squared
     * @param[in] p22 p2 squared
     * @param[in] p32 p3 squared
     * @param[in] p42 p4 squared
     * @param[in] m1 mass m_1
     * @param[in] m2 mass m_2
     * @param[in] m3 mass m_3
     * @param[in] m4 mass m_4
     * @return the function D11
     */
    complex D26(const double mu,const double p12, const double p22,
                         const double p32,const double p42, const double s,
                         const double m1, const double m2, const double m3, 
                         const double m4, const double theta) const;
    /**
     * @brief the scalar four-point Passarino-Veltman function D11
     * @param[in] p12 p1 squared
     * @param[in] p22 p2 squared
     * @param[in] p32 p3 squared
     * @param[in] p42 p4 squared
     * @param[in] m1 mass m_1
     * @param[in] m2 mass m_2
     * @param[in] m3 mass m_3
     * @param[in] m4 mass m_4
     * @return the function D11
     */
    complex D23(const double mu,const double p12, const double p22,
                         const double p32,const double p42, const double s,
                         const double m1, const double m2, const double m3, 
                         const double m4, const double theta) const;
    
    
    
    
private:

};

#endif	/* PVFUNCTIONS_H */

