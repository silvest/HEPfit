/* 
 * File:   TwoLoopQCD.h
 * Author: mishima
 */

#ifndef TWOLOOPQCD_H
#define	TWOLOOPQCD_H

#include <StandardModel.h>
#include "EWSMcommon.h"

using namespace gslpp;


class TwoLoopQCD : public EWSMcommon {

public:

    /**
     * @brief TwoLoopQCD constructor
     * @param[in] EWSMC_i reference to an EWSMcommon object
     */
    TwoLoopQCD(const EWSMcommon& EWSMC_i);

    /**
     * @brief TwoLoopQCD copy constructor
     * @param[in] orig reference to a TwoLoopQCD object
     */
    //TwoLoopQCD(const TwoLoopQCD& orig);

    /**
     * @brief TwoLoopQCD destructor
     */
    virtual ~TwoLoopQCD();


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief leptonic contribution to alpha
     * @return Delta alpha_{lept}^{alpha alpha_s}
     */
    double DeltaAlpha_l() const;

    /**
     * @brief top-quark contribution to alpha
     * @return Delta alpha_{top}^{alpha alpha_s}
     */
    double DeltaAlpha_t() const;
    
    /**
     * @brief leading contribution to Delta r
     * @return Delta rho^{G alpha_s}
     */
    double DeltaRho() const;

    /**
     * @brief remainder contribution to Delta r
     * @return Delta r_{rem}^{alpha alpha_s}
     */
    double DeltaR_rem() const;

    /**
     * @brief remainder contribution to rho_Z^l
     * @param[in] l name of a lepton 
     * @return delta rho_{rem}^{l, G alpha_s}
     */
    complex deltaRho_rem_l(const StandardModel::lepton l) const;

    /**
     * @brief remainder contribution to rho_Z^q
     * @param[in] q name of a quark 
     * @return delta rho_{rem}^{q, G alpha_s}
     */
    complex deltaRho_rem_q(const StandardModel::quark q) const;

    /**
     * @brief remainder contribution to kappa_Z^l
     * @param[in] l name of a lepton 
     * @return delta kappa_{rem}^{l, G alpha_s}
     */
    complex deltaKappa_rem_l(const StandardModel::lepton l) const;

    /**
     * @brief remainder contribution to kappa_Z^q
     * @param[in] q name of a quark 
     * @return delta kappa_{rem}^{q, G alpha_s}
     */
    complex deltaKappa_rem_q(const StandardModel::quark q) const;
    
    
    ////////////////////////////////////////////////////////////////////////        
    
private:


    ////////////////////////////////////////////////////////////////////////    
     

};

#endif	/* TWOLOOPQCD_H */

