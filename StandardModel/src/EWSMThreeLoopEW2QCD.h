/* 
 * File:   EWSMThreeLoopEW2QCD.h
 * Author: mishima
 */

#ifndef EWSMTHREELOOPEW2QCD_H
#define	EWSMTHREELOOPEW2QCD_H

#include "EWSMcache.h"
using namespace gslpp;


class EWSMThreeLoopEW2QCD {

public:
 
    /**
     * @brief EWSMThreeLoopEW2QCD constructor
     * @param[in] cache_i reference to an EWSMcommon object
     */
    EWSMThreeLoopEW2QCD(const EWSMcache& cache_i);

    
    ////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief leptonic contribution to alpha
     * @return Delta alpha_{lept}^{alpha^2 alpha_s}
     */
    double DeltaAlpha_l() const;

    /**
     * @brief top-quark contribution to alpha
     * @return Delta alpha_{top}^{alpha^2 alpha_s}
     */
    double DeltaAlpha_t() const;
    
    /**
     * @brief leading contribution to Delta r
     * @param[in] Mw_i the W-boson mass
     * @return Delta rho^{alpha^2 alpha_s}
     */
    double DeltaRho(const double Mw_i) const;

    /**
     * @brief remainder contribution to Delta r
     * @param[in] Mw_i the W-boson mass
     * @return Delta r_{rem}^{alpha^2 alpha_s}
     */
    double DeltaR_rem(const double Mw_i) const;

    /**
     * @brief remainder contribution to rho_Z^f
     * @param[in] f StandardModel::quark or StandardModel::lepton 
     * @param[in] Mw_i the W-boson mass
     * @return delta rho_{rem}^{f, alpha^2 alpha_s}
     */
    template<typename T> complex deltaRho_rem_f(const T f, const double Mw_i) const;

    /**
     * @brief remainder contribution to kappa_Z^f
     * @param[in] f StandardModel::quark or StandardModel::lepton 
     * @param[in] Mw_i the W-boson mass
     * @return delta kappa_{rem}^{f, alpha^2 alpha_s}
     */
    template<typename T> complex deltaKappa_rem_f(const T f, const double Mw_i) const;
    
    
    ////////////////////////////////////////////////////////////////////////        
    
private:
    const EWSMcache& cache;
    

    ////////////////////////////////////////////////////////////////////////    
     
    
};

#endif	/* EWSMTHREELOOPEW2QCD_H */

