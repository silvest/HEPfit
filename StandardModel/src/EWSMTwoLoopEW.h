/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSMTWOLOOPEW_H
#define	EWSMTWOLOOPEW_H

#include "EWSMcache.h"
#include "EWSMOneLoopEW.h"
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
     * @param[in] s invariant mass squared 
     * @return Delta alpha_{lept}^{alpha^2}
     */
    double DeltaAlpha_l(const double s) const;

    /**
     * @brief top-quark contribution to alpha
     * @param[in] s invariant mass squared 
     * @return Delta alpha_{top}^{alpha^2}
     */
    double DeltaAlpha_t(const double s) const;    
    
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
     * @brief remainder contribution to rho_Z^l
     * @param[in] l name of lepton
     * @param[in] Mw_i the W-boson mass
     * @return delta rho_{rem}^{l, alpha^2}
     */
    complex deltaRho_rem_l(const StandardModel::lepton l, const double Mw_i) const;

    /**
     * @brief remainder contribution to rho_Z^q
     * @param[in] q name of quark
     * @param[in] Mw_i the W-boson mass
     * @return delta rho_{rem}^{q, alpha^2}
     */
    complex deltaRho_rem_q(const StandardModel::quark q, const double Mw_i) const;

    /**
     * @brief remainder contribution to kappa_Z^l
     * @param[in] l name of lepton
     * @param[in] Mw_i the W-boson mass
     * @return delta kappa_{rem}^{l, alpha^2}
     */
    complex deltaKappa_rem_l(const StandardModel::lepton l, const double Mw_i) const;
                                                  
    /**
     * @brief remainder contribution to kappa_Z^q
     * @param[in] q name of quark
     * @param[in] Mw_i the W-boson mass
     * @return delta kappa_{rem}^{q, alpha^2}
     */
    complex deltaKappa_rem_q(const StandardModel::quark q, const double Mw_i) const;

    
    ////////////////////////////////////////////////////////////////////////        
    // O(GF^2 Mt^2) contributions
    
    /**
     * @return O(alpha^2 Mt^4/M_Z^4) contribution to Delta rho
     */
    double rho_2() const;    
    
    /**
     * @return O(alpha^2 Mt^4/M_Z^4) contribution to the Z-b-bbar vertex
     */
    double tau_2() const;
    
        
    ////////////////////////////////////////////////////////////////////////        
    
private:
    const EWSMcache& cache;
    const EWSMOneLoopEW myOneLoopEW;
    
    /**
     * @param[in] a a=(m_h/M_t)^2
     * @return g(a)
     */
    double g(const double a) const;

    /**
     * @param[in] a a=(m_h/M_t)^2
     * @return f(a,0)
     */
    double f0(const double a) const;

    /**
     * @param[in] a a=(m_h/M_t)^2
     * @return f(a,1)
     */
    double f1(const double a) const;    

    
    ////////////////////////////////////////////////////////////////////////        
    // O(alpha^2 Mt^4/M_Z^4 + alpha^2 Mt^2/M_Z^2) contributions

    double DeltaRho2(const double Mw_i) const;
    double DeltaRho2Add(const double Mw_i) const;
    double DeltaRw2(const double Mw_i) const;
    double deltaEoverE2() const;
    double f2Add(const double Mw_i) const;
    double DeltaEta2(const double Mw_i) const;
    complex DeltaEta2Add_tmp(const double I3f, const double Qf, const double Mw_i) const;
    complex DeltaEta2Add_l(const StandardModel::lepton l, const double Mw_i) const;
    complex DeltaEta2Add_q(const StandardModel::quark q, const double Mw_i) const;
    double DeltaKappa2(const double Mw_i) const;
    complex DeltaKappa2Add_tmp(const double I3f, const double Qf, const double Mw_i) const;
    complex DeltaKappa2Add_l(const StandardModel::lepton l, const double Mw_i) const;
    complex DeltaKappa2Add_q(const StandardModel::quark q, const double Mw_i) const;

    complex Vadd(const double I3f, const double Qf, const double Mw_i) const;
    complex DeltaEtaf1(const double I3f, const double Qf, const double Mw_i) const;
    complex Vfi(const double I3f, const double Qf, const double q2, const double Mw_i) const;

    double Lambda(const double x) const;
    double phi(const double x) const;
    complex FV(const double x) const;
    complex GV(const double x) const;    
    
};

#endif	/* EWSMTWOLOOPEW_H */

