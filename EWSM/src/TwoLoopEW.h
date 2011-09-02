/* 
 * File:   TwoLoopEW.h
 * Author: mishima
 */

#ifndef TWOLOOPEW_H
#define	TWOLOOPEW_H

#include <StandardModel.h>
#include "EWSMcommon.h"

using namespace gslpp;


class TwoLoopEW : public EWSMcommon {

public:

    /**
     * @brief TwoLoopEW constructor
     * @param[in] EWSMC_i reference to an EWSMcommon object
     */
    TwoLoopEW(const EWSMcommon& EWSMC_i);

    /**
     * @brief TwoLoopEW copy constructor
     * @param[in] orig reference to a TwoLoopEW object
     */
    //TwoLoopEW(const TwoLoopEW& orig);

    /**
     * @brief TwoLoopEW destructor
     */
    virtual ~TwoLoopEW();

    
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
     * @return Delta rho^{G^2}
     */
    double DeltaRho() const;

    /**
     * @brief remainder contribution to Delta r
     * @return Delta r_{rem}^{alpha^2}
     */
    double DeltaR_rem() const;

    /**
     * @brief remainder contribution to rho_Z^l
     * @param[in] l name of a lepton 
     * @return delta rho_{rem}^{l, G^2}
     */
    complex deltaRho_rem_l(const StandardModel::lepton l) const;

    /**
     * @brief remainder contribution to rho_Z^q
     * @param[in] q name of a quark 
     * @return delta rho_{rem}^{q, G^2}
     */
    complex deltaRho_rem_q(const StandardModel::quark q) const;

    /**
     * @brief remainder contribution to kappa_Z^l
     * @param[in] l name of a lepton 
     * @return delta kappa_{rem}^{l, G^2}
     */
    complex deltaKappa_rem_l(const StandardModel::lepton l) const;

    /**
     * @brief remainder contribution to kappa_Z^q
     * @param[in] q name of a quark 
     * @return delta kappa_{rem}^{q, G^2}
     */
    complex deltaKappa_rem_q(const StandardModel::quark q) const;
    
    
    ////////////////////////////////////////////////////////////////////////        
    
private:


    ////////////////////////////////////////////////////////////////////////    
     

};

#endif	/* TWOLOOPEW_H */

