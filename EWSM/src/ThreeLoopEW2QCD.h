/* 
 * File:   ThreeLoopEW2QCD.h
 * Author: mishima
 */

#ifndef THREELOOPEW2QCD_H
#define	THREELOOPEW2QCD_H

#include <StandardModel.h>
#include "EWSMcommon.h"

using namespace gslpp;


class ThreeLoopEW2QCD {

public:
 
    /**
     * @brief ThreeLoopEW2QCD constructor
     * @param[in] EWSMC_i reference to an EWSMcommon object
     */
    ThreeLoopEW2QCD(const EWSMcommon& EWSMC_i);

    /**
     * @brief ThreeLoopEW2QCD copy constructor
     * @param[in] orig reference to a ThreeLoopEW2QCD object
     */
    //ThreeLoopEW2QCD(const ThreeLoopEW2QCD& orig);

    /**
     * @brief ThreeLoopEW2QCD destructor
     */
    virtual ~ThreeLoopEW2QCD();

    
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
     * @return Delta rho^{alpha^2 alpha_s}
     */
    double DeltaRho() const;

    /**
     * @brief remainder contribution to Delta r
     * @return Delta r_{rem}^{alpha^2 alpha_s}
     */
    double DeltaR_rem() const;

    /**
     * @brief remainder contribution to rho_Z^l
     * @param[in] l name of a lepton 
     * @return delta rho_{rem}^{l, alpha^2 alpha_s}
     */
    complex deltaRho_rem_l(const StandardModel::lepton l) const;

    /**
     * @brief remainder contribution to rho_Z^q
     * @param[in] q name of a quark 
     * @return delta rho_{rem}^{q, alpha^2 alpha_s}
     */
    complex deltaRho_rem_q(const StandardModel::quark q) const;

    /**
     * @brief remainder contribution to kappa_Z^l
     * @param[in] l name of a lepton 
     * @return delta kappa_{rem}^{l, alpha^2 alpha_s}
     */
    complex deltaKappa_rem_l(const StandardModel::lepton l) const;

    /**
     * @brief remainder contribution to kappa_Z^q
     * @param[in] q name of a quark 
     * @return delta kappa_{rem}^{q, alpha^2 alpha_s}
     */
    complex deltaKappa_rem_q(const StandardModel::quark q) const;
    
    
    ////////////////////////////////////////////////////////////////////////        
    
private:
    const EWSMcommon& EWSMC;
    

    ////////////////////////////////////////////////////////////////////////    
     
    
};

#endif	/* THREELOOPEW2QCD_H */

