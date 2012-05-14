/* 
 * File:   EWSMTwoLoopEW.h
 * Author: mishima
 */

#ifndef EWSMTWOLOOPEW_H
#define	EWSMTWOLOOPEW_H

#include "EWSMcache.h"
using namespace gslpp;


class EWSMTwoLoopEW {

public:

    /**
     * @brief TwoLoopEW constructor
     * @param[in] cache_i reference to an EWSMcommon object
     */
    EWSMTwoLoopEW(const EWSMcache& cache_i);

    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief leptonic contribution to alpha
     * @return Delta alpha_{lept}^{alpha^2}
     */
    double DeltaAlpha_l() const;

    /**
     * @brief top-quark contribution to alpha
     * @return Delta alpha_{top}^{alpha^2}
     */
    double DeltaAlpha_t() const;    
    
    /**
     * @brief leading contribution to Delta r
     * @param[in] Mw_i the W-boson mass
     * @return Delta rho^{alpha^2}
     */
    double DeltaRho(const double Mw_i) const;

    /**
     * @brief remainder contribution to Delta r
     * @param[in] Mw_i the W-boson mass
     * @return Delta r_{rem}^{alpha^2}
     */
    double DeltaR_rem(const double Mw_i) const;

    /**
     * @brief remainder contribution to rho_Z^f
     * @param[in] f StandardModel::quark or StandardModel::lepton  
     * @param[in] Mw_i the W-boson mass
     * @return delta rho_{rem}^{f, alpha^2}
     */
    template<typename T> complex deltaRho_rem_f(const T f, const double Mw_i) const;

    /**
     * @brief remainder contribution to kappa_Z^f
     * @param[in] f StandardModel::quark or StandardModel::lepton 
     * @param[in] Mw_i the W-boson mass
     * @return delta kappa_{rem}^{f, alpha^2}
     */
    template<typename T> complex deltaKappa_rem_f(const T f, const double Mw_i) const;

    
    ////////////////////////////////////////////////////////////////////////        
    
    /**
     * @param[in] Mw_i the W-boson mass
     * @return O(alpha^2 Mt^4/M_Z^4) contribution to Delta rho
     */
    complex rho_2(const double Mw_i) const;
    
    
    /**
     * @param[in] Mw_i the W-boson mass
     * @return O(alpha^2 Mt^4/M_Z^4) contribution to the Z-b-bbar vertex
     */
    complex tau_2(const double Mw_i) const;
    
    
    ////////////////////////////////////////////////////////////////////////        
    
private:
    const EWSMcache& cache;
     

};

#endif	/* EWSMTWOLOOPEW_H */

