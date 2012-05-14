/* 
 * File:   EWSMThreeLoopEW.h
 * Author: mishima
 */

#ifndef EWSMTHREELOOPEW_H
#define	EWSMTHREELOOPEW_H

#include "EWSMcache.h"
using namespace gslpp;


class EWSMThreeLoopEW {

public:

    /**
     * @brief EWSMThreeLoopEW constructor
     * @param[in] cache_i reference to an EWSMcommon object
     */
    EWSMThreeLoopEW(const EWSMcache& cache_i);

    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief leptonic contribution to alpha
     * @return Delta alpha_{lept}^{alpha^3}
     */
    double DeltaAlpha_l() const;

    /**
     * @brief top-quark contribution to alpha
     * @return Delta alpha_{top}^{alpha^3}
     */
    double DeltaAlpha_t() const;
    
    /**
     * @brief leading contribution to Delta r
     * @param[in] Mw_i the W-boson mass
     * @return Delta rho^{alpha^3}
     */
    double DeltaRho(const double Mw_i) const;

    /**
     * @brief remainder contribution to Delta r
     * @param[in] Mw_i the W-boson mass
     * @return Delta r_{rem}^{alpha^3}
     */
    double DeltaR_rem(const double Mw_i) const;

    /**
     * @brief remainder contribution to rho_Z^f
     * @param[in] f StandardModel::quark or StandardModel::lepton 
     * @param[in] Mw_i the W-boson mass
     * @return delta rho_{rem}^{f, alpha^3}
     */
    template<typename T> complex deltaRho_rem_f(const T f, const double Mw_i) const;

    /**
     * @brief remainder contribution to kappa_Z^f
     * @param[in] f StandardModel::quark or StandardModel::lepton
     * @param[in] Mw_i the W-boson mass
     * @return delta kappa_{rem}^{f, alpha^3}
     */
    template<typename T> complex deltaKappa_rem_f(const T f, const double Mw_i) const;

    
    ////////////////////////////////////////////////////////////////////////        
    
private:
    const EWSMcache& cache;
    

    ////////////////////////////////////////////////////////////////////////    
     
    
};

#endif	/* EWSMTHREELOOPEW_H */

