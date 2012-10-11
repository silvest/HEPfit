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
     * @param[in] s invariant mass squared 
     * @return Delta alpha_{lept}^{alpha^3}
     */
    double DeltaAlpha_l(const double s) const;

    /**
     * @brief top-quark contribution to alpha
     * @param[in] s invariant mass squared 
     * @return Delta alpha_{top}^{alpha^3}
     */
    double DeltaAlpha_t(const double s) const;
    
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
     * @brief remainder contribution to rho_Z^l
     * @param[in] l name of lepton
     * @param[in] Mw_i the W-boson mass
     * @return delta rho_{rem}^{l, alpha^3}
     */
    complex deltaRho_rem_l(const StandardModel::lepton l, const double Mw_i) const;

    /**
     * @brief remainder contribution to rho_Z^q
     * @param[in] q name of quark
     * @param[in] Mw_i the W-boson mass
     * @return delta rho_{rem}^{q, alpha^3}
     */
    complex deltaRho_rem_q(const StandardModel::quark q, const double Mw_i) const;

    /**
     * @brief remainder contribution to kappa_Z^l
     * @param[in] l name of lepton
     * @param[in] Mw_i the W-boson mass
     * @return delta kappa_{rem}^{l, alpha^3}
     */
    complex deltaKappa_rem_l(const StandardModel::lepton l, const double Mw_i) const;
                                                  
    /**
     * @brief remainder contribution to kappa_Z^q
     * @param[in] q name of quark
     * @param[in] Mw_i the W-boson mass
     * @return delta kappa_{rem}^{q, alpha^3}
     */
    complex deltaKappa_rem_q(const StandardModel::quark q, const double Mw_i) const;

    
    ////////////////////////////////////////////////////////////////////////        
    
private:
    const EWSMcache& cache;
    

    ////////////////////////////////////////////////////////////////////////    
     
    
};

#endif	/* EWSMTHREELOOPEW_H */

