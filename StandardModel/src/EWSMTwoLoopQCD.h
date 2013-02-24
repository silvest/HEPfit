/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSMTWOLOOPQCD_H
#define	EWSMTWOLOOPQCD_H

#include "EWSMcache.h"
using namespace gslpp;

/**
 * @class EWSMTwoLoopQCD
 * @ingroup StandardModel
 * @brief A class for @f$O(\alpha_s^2)@f$ two-loop radiative corrections to the EW precision observables.  
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class EWSMTwoLoopQCD {

public:

    /**
     * @brief EWSMTwoLoopQCD constructor
     * @param[in] cache_i reference to an EWSMcommon object
     */
    EWSMTwoLoopQCD(const EWSMcache& cache_i);


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief leptonic contribution to alpha
     * @param[in] s invariant mass squared 
     * @return Delta alpha_{lept}^{alpha alpha_s}
     */
    double DeltaAlpha_l(const double s) const;

    /**
     * @brief top-quark contribution to alpha
     * @param[in] s invariant mass squared 
     * @return Delta alpha_{top}^{alpha alpha_s}
     */
    double DeltaAlpha_t(const double s) const;
    
    /**
     * @brief leading contribution to Delta r
     * @param[in] Mw_i the W-boson mass
     * @return Delta rho^{alpha alpha_s}
     */
    double DeltaRho(const double Mw_i) const;

    /**
     * @brief remainder contribution to Delta r
     * @param[in] Mw_i the W-boson mass
     * @return Delta r_{rem}^{alpha alpha_s}
     */
    double DeltaR_rem(const double Mw_i) const;

    /**
     * @brief remainder contribution to rho_Z^l
     * @param[in] l name of lepton
     * @param[in] Mw_i the W-boson mass
     * @return delta rho_{rem}^{l, alpha alpha_s}
     */
    complex deltaRho_rem_l(const StandardModel::lepton l, const double Mw_i) const;

    /**
     * @brief remainder contribution to rho_Z^q
     * @param[in] q name of quark
     * @param[in] Mw_i the W-boson mass
     * @return delta rho_{rem}^{q, alpha alpha_s}
     */
    complex deltaRho_rem_q(const StandardModel::quark q, const double Mw_i) const;

    /**
     * @brief remainder contribution to kappa_Z^l
     * @param[in] l name of lepton
     * @param[in] Mw_i the W-boson mass
     * @return delta kappa_{rem}^{l, alpha alpha_s}
     */
    complex deltaKappa_rem_l(const StandardModel::lepton l, const double Mw_i) const;
                                                  
    /**
     * @brief remainder contribution to kappa_Z^q
     * @param[in] q name of quark
     * @param[in] Mw_i the W-boson mass
     * @return delta kappa_{rem}^{q, alpha alpha_s}
     */
    complex deltaKappa_rem_q(const StandardModel::quark q, const double Mw_i) const;  
    

    ////////////////////////////////////////////////////////////////////////        
    
    /**
     * @return delta^QCD_2 for the leading contribution to Delta rho
     */
    double deltaQCD_2() const;

    /**
     * @param[in] x x=s/M_t^2
     * @param[in] Mw_i the W-boson mass
     * @return F_1(x)
     * @attention valid for 0<=x<1
     */
    double F1(const double x, const double Mw_i) const;
    
    /**
     * @param[in] r r=s/(4M_t^2)
     * @return V_1(r)
     * @attention valid for 0<=r<1
     */
    double V1(const double r) const;
        
    /**
     * @param[in] r r=s/(4M_t^2)
     * @return FA_1(r)
     * @attention valid for 0<=r<1
     */
    double A1(const double r) const;

    /**
     * @param[in] r r=s/(4M_t^2)
     * @return the derivative of V_1(r)
     * @attention valid for r << 1
     */
    double V1prime(const double r) const;
    
    /**
     * @param[in] r r=s/(4M_t^2)
     * @return the derivative of A_1(r)
     * @attention valid for r << 1
     */
    double A1prime(const double r) const;
    
    /**
     * @param[in] Mw_i the W-boson mass
     * @return Delta r^{ud}, contribution to Delta r from the light-quark doublets
     */
    double DeltaR_ud(const double Mw_i) const;

    /**
     * @param[in] Mw_i the W-boson mass
     * @return Delta r^{tb}, contribution to Delta r from the top-bottom doublet
     */
    double DeltaR_tb(const double Mw_i) const;
        
    /**
     * @param[in] Mw_i the W-boson mass
     * @return Delta rho^{ud}, contribution to rho_Z^f from the light-quark doublets
     */
    double DeltaRho_ud(const double Mw_i) const;

    /**
     * @param[in] Mw_i the W-boson mass
     * @return Delta rho^{tb}, contribution to rho_Z^f from the top-bottom doublet
     */
    double DeltaRho_tb(const double Mw_i) const;
    
    /**
     * @param[in] Mw_i the W-boson mass
     * @return Delta kappa^{ud}, contribution to kappa_Z^f from the light-quark doublets
     */
    complex DeltaKappa_ud(const double Mw_i) const;
    
    /**
     * @param[in] Mw_i the W-boson mass
     * @return Delta kappa^{tb}, contribution to kappa_Z^f from the top-bottom doublet
     */
    complex DeltaKappa_tb(const double Mw_i) const;
     
    
    ////////////////////////////////////////////////////////////////////////        
    
private:
    const EWSMcache& cache;
    
    
};

#endif	/* EWSMTWOLOOPQCD_H */

