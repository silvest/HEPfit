/* 
 * File:   EWSMThreeLoopQCD.h
 * Author: mishima
 */

#ifndef EWSMTHREELOOPQCD_H
#define	EWSMTHREELOOPQCD_H

#include "EWSMcache.h"
using namespace gslpp;


class EWSMThreeLoopQCD {

public:

    /**
     * @brief EWSMThreeLoopQCD constructor
     * @param[in] cache_i reference to an EWSMcommon object
     */
    EWSMThreeLoopQCD(const EWSMcache& cache_i);


    ////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief leptonic contribution to alpha
     * @return Delta alpha_{lept}^{alpha alpha_s^2}
     */
    double DeltaAlpha_l() const;

    /**
     * @brief top-quark contribution to alpha
     * @return Delta alpha_{top}^{alpha alpha_s^2}
     */
    double DeltaAlpha_t() const;
    
    /**
     * @brief leading contribution to Delta r
     * @param[in] Mw_i the W-boson mass
     * @return Delta rho^{alpha alpha_s^2}
     */
    double DeltaRho(const double Mw_i) const;

    /**
     * @brief remainder contribution to Delta r
     * @param[in] Mw_i the W-boson mass
     * @return Delta r_{rem}^{alpha alpha_s^2}
     */
    double DeltaR_rem(const double Mw_i) const;

    /**
     * @brief remainder contribution to rho_Z^f
     * @param[in] f StandardModel::quark or StandardModel::lepton 
     * @param[in] Mw_i the W-boson mass
     * @return delta rho_{rem}^{f, alpha alpha_s^2}
     */
    template<typename T> complex deltaRho_rem_f(const T f, const double Mw_i) const;

    /**
     * @brief remainder contribution to kappa_Z^f
     * @param[in] f StandardModel::quark or StandardModel::lepton 
     * @param[in] Mw_i the W-boson mass
     * @return delta kappa_{rem}^{f, alpha alpha_s^2}
     */
    template<typename T> complex deltaKappa_rem_f(const T f, const double Mw_i) const;
    
    
    ////////////////////////////////////////////////////////////////////////        
    
private:
    const EWSMcache& cache;

    
    ////////////////////////////////////////////////////////////////////////        
    
    /**
     * @param[in] Mw_i the W-boson mass
     * @return delta^QCD_3 for the leading contribution to Delta rho
     */
    double deltaQCD_3(const double Mw_i) const;    
    
    /**
     * @param[in] Mw_i the W-boson mass
     * @return delta^QCD_{kappa,3} for the remainder contribution
     */
    complex deltaQCD_kappa3(const double Mw_i) const;
    
    
};

#endif	/* EWSMTHREELOOPQCD_H */

