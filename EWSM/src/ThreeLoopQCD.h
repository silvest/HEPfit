/* 
 * File:   ThreeLoopQCD.h
 * Author: mishima
 */

#ifndef THREELOOPQCD_H
#define	THREELOOPQCD_H

#include <StandardModel.h>
#include "EWSMcommon.h"

using namespace gslpp;


class ThreeLoopQCD {

public:

    /**
     * @brief ThreeLoopQCD constructor
     * @param[in] EWSMC_i reference to an EWSMcommon object
     */
    ThreeLoopQCD(const EWSMcommon& EWSMC_i);


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
     * @return Delta rho^{alpha alpha_s^2}
     */
    double DeltaRho() const;

    /**
     * @brief remainder contribution to Delta r
     * @return Delta r_{rem}^{alpha alpha_s^2}
     */
    double DeltaR_rem() const;

    /**
     * @brief remainder contribution to rho_Z^l
     * @param[in] l name of a lepton 
     * @return delta rho_{rem}^{l, alpha alpha_s^2}
     */
    complex deltaRho_rem_l(const StandardModel::lepton l) const;

    /**
     * @brief remainder contribution to rho_Z^q
     * @param[in] q name of a quark 
     * @return delta rho_{rem}^{q, alpha alpha_s^2}
     */
    complex deltaRho_rem_q(const StandardModel::quark q) const;

    /**
     * @brief remainder contribution to kappa_Z^l
     * @param[in] l name of a lepton 
     * @return delta kappa_{rem}^{l, alpha alpha_s^2}
     */
    complex deltaKappa_rem_l(const StandardModel::lepton l) const;

    /**
     * @brief remainder contribution to kappa_Z^q
     * @param[in] q name of a quark 
     * @return delta kappa_{rem}^{q, alpha alpha_s^2}
     */
    complex deltaKappa_rem_q(const StandardModel::quark q) const;
    
    
    ////////////////////////////////////////////////////////////////////////        
    
private:
    const EWSMcommon& EWSMC;

    
    ////////////////////////////////////////////////////////////////////////        
    
    /**
     * @return delta^QCD_3 for the leading contribution to Delta rho
     */
    double deltaQCD_3() const;    
    
    /**
     * @return delta^QCD_{kappa,3} for the remainder contribution
     */
    complex deltaQCD_kappa3() const;
    
    
};

#endif	/* THREELOOPQCD_H */

