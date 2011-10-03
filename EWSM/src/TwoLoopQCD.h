/* 
 * File:   TwoLoopQCD.h
 * Author: mishima
 */

#ifndef TWOLOOPQCD_H
#define	TWOLOOPQCD_H

#include <StandardModel.h>
#include "EWSMcommon.h"

using namespace gslpp;


class TwoLoopQCD {

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
     * @return Delta rho^{alpha alpha_s}
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
     * @return delta rho_{rem}^{l, alpha alpha_s}
     */
    complex deltaRho_rem_l(const StandardModel::lepton l) const;

    /**
     * @brief remainder contribution to rho_Z^q
     * @param[in] q name of a quark 
     * @return delta rho_{rem}^{q, alpha alpha_s}
     */
    complex deltaRho_rem_q(const StandardModel::quark q) const;

    /**
     * @brief remainder contribution to kappa_Z^l
     * @param[in] l name of a lepton 
     * @return delta kappa_{rem}^{l, alpha alpha_s}
     */
    complex deltaKappa_rem_l(const StandardModel::lepton l) const;

    /**
     * @brief remainder contribution to kappa_Z^q
     * @param[in] q name of a quark 
     * @return delta kappa_{rem}^{q, alpha alpha_s}
     */
    complex deltaKappa_rem_q(const StandardModel::quark q) const;
    

    ////////////////////////////////////////////////////////////////////////        
    
    /**
     * @return delta^QCD_2 for the leading contribution to Delta rho
     */
    double deltaQCD_2() const;

    /**
     * @param[in] x x=s/M_t^2
     * @return F_1(x)
     * @attention valid for 0<=x<1
     */
    double F1(const double x) const;
    
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
     * @return Delta r^{ud}, contribution to Delta r from the light-quark doublets
     */
    double DeltaR_ud() const;

    /**
     * @return Delta r^{tb}, contribution to Delta r from the top-bottom doublet
     */
    double DeltaR_tb() const;
        
    /**
     * @return Delta rho^{ud}, contribution to rho_Z^f from the light-quark doublets
     */
    double DeltaRho_ud() const;

    /**
     * @return Delta rho^{tb}, contribution to rho_Z^f from the top-bottom doublet
     */
    complex DeltaRho_tb() const;
    
    /**
     * @return Delta kappa^{ud}, contribution to kappa_Z^f from the light-quark doublets
     */
    double DeltaKappa_ud() const;
    
    /**
     * @return Delta kappa^{tb}, contribution to kappa_Z^f from the top-bottom doublet
     */
    complex DeltaKappa_tb() const;
     
    
    ////////////////////////////////////////////////////////////////////////        
    
private:
    const EWSMcommon& EWSMC;
    
    
};

#endif	/* TWOLOOPQCD_H */

