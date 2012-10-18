/* 
 * File:   EWSMApproximateFormulae.h
 * Author: mishima
 */

#ifndef EWSMAPPROXIMATEFORMULAE_H
#define	EWSMAPPROXIMATEFORMULAE_H

#include "StandardModel.h"


class EWSMApproximateFormulae {
public:
      
    /**
     * @brief  EWSMApproximateFormulae constructor
     * @param[in] SM_i reference to a StandardModel object
     */
    EWSMApproximateFormulae(const StandardModel& SM_i);    

    
    ////////////////////////////////////////////////////////////////////////

    /** 
     * @param[in] DeltaAlpha_i the sum of the leptonic and hadronic corrections to alpha at Mz
     * @return the W-boson mass from an approximate two-loop formula
     */
    double Mw(const double DeltaAlphaL5q_i) const;
    
    /**
     * @param[in] l name of lepton
     * @param[in] DeltaAlphaL5q_i the sum of the leptonic and hadronic corrections to alpha at Mz
     * @return the effective weak mixing angle for Z->l lbar from an approximate two-loop formula
     */
    double sin2thetaEff_l(const StandardModel::lepton l, const double DeltaAlphaL5q_i) const;

    /**
     * @param[in] q name of quark
     * @param[in] DeltaAlphaL5q_i the sum of the leptonic and hadronic corrections to alpha at Mz
     * @return the effective weak mixing angle for Z->q qbar from an approximate two-loop formula
     */
    double sin2thetaEff_q(const StandardModel::quark q, const double DeltaAlphaL5q_i) const;    
    
    /** 
     * @param[in] DeltaAlphaL5q_i the sum of the leptonic and hadronic corrections to alpha at Mz
     * @return EW two-loop contribution to Delta r from an approximate formula
     */
    double DeltaR_TwoLoopEW(const double DeltaAlphaL5q_i) const;

    /**
     * @brief EW two-loop contribution to Delta kappa_Z^l = kappaZ^l - 1  
     * @param[in] DeltaAlphaL5q_i the sum of the leptonic and hadronic corrections to alpha at Mz
     * @return Delta kappa^{l, alpha^2} from an approximate formula
     */
    double DeltaKappa_l_TwoLoopEW(const double DeltaAlphaL5q_i) const;    

    /**
     * @brief EW two-loop contribution to Delta kappa_Z^b = kappaZ^b - 1 
     * @param[in] DeltaAlphaL5q_i the sum of the leptonic and hadronic corrections to alpha at Mz
     * @return Delta kappa^{b, alpha^2} from an approximate formula
     */
    double DeltaKappa_b_TwoLoopEW(const double DeltaAlphaL5q_i) const;    
    
    /**
     * @brief R_b^0 with the complete fermionic EW two-loop corrections
     * @param[in] DeltaAlphaL5q_i the sum of the leptonic and hadronic corrections to alpha at Mz
     * @return R_b^0 from an approximate two-loop formula
     */
    double R0_bottom(const double DeltaAlphaL5q_i) const;
    
    
    ////////////////////////////////////////////////////////////////////////
    
private:
    bool bDebug; // true for debugging    
    
    const StandardModel& SM;
    
};

#endif	/* EWSMAPPROXIMATEFORMULAE_H */

