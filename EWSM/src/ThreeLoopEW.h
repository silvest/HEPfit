/* 
 * File:   ThreeLoopEW.h
 * Author: mishima
 */

#ifndef THREELOOPEW_H
#define	THREELOOPEW_H

#include <StandardModel.h>
#include "EWSMcommon.h"

using namespace gslpp;


class ThreeLoopEW : public EWSMcommon {

public:

    /**
     * @brief ThreeLoopEW constructor
     * @param[in] EWSMC_i reference to an EWSMcommon object
     */
    ThreeLoopEW(const EWSMcommon& EWSMC_i);

    /**
     * @brief ThreeLoopEW copy constructor
     * @param[in] orig reference to a ThreeLoopEW object
     */
    //ThreeLoopEW(const ThreeLoopEW& orig);

    /**
     * @brief ThreeLoopEW destructor
     */
    virtual ~ThreeLoopEW();

    
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
     * @return Delta rho^{G^3}
     */
    double DeltaRho() const;

    /**
     * @brief remainder contribution to Delta r
     * @return Delta r_{rem}^{alpha^3}
     */
    double DeltaR_rem() const;

    /**
     * @brief remainder contribution to rho_Z^l
     * @param[in] l name of a lepton 
     * @return delta rho_{rem}^{l, G^3}
     */
    complex deltaRho_rem_l(const StandardModel::lepton l) const;

    /**
     * @brief remainder contribution to rho_Z^q
     * @param[in] q name of a quark 
     * @return delta rho_{rem}^{q, G^3}
     */
    complex deltaRho_rem_q(const StandardModel::quark q) const;

    /**
     * @brief remainder contribution to kappa_Z^l
     * @param[in] l name of a lepton 
     * @return delta kappa_{rem}^{l, G^3}
     */
    complex deltaKappa_rem_l(const StandardModel::lepton l) const;

    /**
     * @brief remainder contribution to kappa_Z^q
     * @param[in] q name of a quark 
     * @return delta kappa_{rem}^{q, G^3}
     */
    complex deltaKappa_rem_q(const StandardModel::quark q) const;

    
    ////////////////////////////////////////////////////////////////////////        
    
private:


    ////////////////////////////////////////////////////////////////////////    
     
    
};

#endif	/* THREELOOPEW_H */

