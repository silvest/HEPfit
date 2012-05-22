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
     * @param[in] bDebug_i boolean value for debugging (true for debugging)
     */
    EWSMApproximateFormulae(const StandardModel& SM_i, const bool bDebug_i=false);    

    
    ////////////////////////////////////////////////////////////////////////

    /** 
     * @param[in] DeltaAlpha_i the total radiative corrections to alpha at Mz
     * @return the W-boson mass from an approximate two-loop formula
     */
    double Mw(const double DeltaAlpha_i) const;
    
    /**
     * @param[in] l name of lepton
     * @param[in] DeltaAlpha_i the total radiative corrections to alpha at Mz
     * @return the effective weak mixing angle for Z->l lbar from an approximate two-loop formula
     */
    double sin2thetaEff_l(const StandardModel::lepton l, const double DeltaAlpha_i) const;

    /**
     * @param[in] q name of quark
     * @param[in] DeltaAlpha_i the total radiative corrections to alpha at Mz
     * @return the effective weak mixing angle for Z->q qbar from an approximate two-loop formula
     */
    double sin2thetaEff_q(const StandardModel::quark q, const double DeltaAlpha_i) const;    
    
    /** 
     * @param[in] DeltaAlpha_i the total radiative corrections to alpha at Mz
     * @return EW two-loop contribution to Delta r from an approximate formula
     */
    double DeltaR_TwoLoopEW(const double DeltaAlpha_i) const;

    /**
     * @brief EW two-loop contribution to Delta kappa_Z^l = kappaZ^l - 1  
     * @param[in] DeltaAlpha_i the total radiative corrections to alpha at Mz
     * @return Delta kappa^{l, alpha^2} from an approximate formula
     */
    double DeltaKappa_l_TwoLoopEW(const double DeltaAlpha_i) const;    

    /**
     * @brief EW two-loop contribution to Delta kappa_Z^b = kappaZ^b - 1 
     * @param[in] DeltaAlpha_i the total radiative corrections to alpha at Mz
     * @return Delta kappa^{b, alpha^2} from an approximate formula
     */
    double DeltaKappa_b_TwoLoopEW(const double DeltaAlpha_i) const;    
    
    
    ////////////////////////////////////////////////////////////////////////
    
private:
    bool bDebug; // true for debugging    
    
    const StandardModel& SM;
    
};

#endif	/* EWSMAPPROXIMATEFORMULAE_H */

