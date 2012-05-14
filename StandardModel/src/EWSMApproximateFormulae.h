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
     * @param[in] DeltaAlpha_i the total radiative corrections to alpha at Mz
     * @return the W-boson mass from an approximate two-loop formula
     */
    double Mw(const double DeltaAlpha_i) const;
    
    /**
     * @param[in] f StandardModel::quark or StandardModel::lepton 
     * @param[in] DeltaAlpha_i the total radiative corrections to alpha at Mz
     * @return the effective weak mixing angle for Z->f\bar{f} from an approximate two-loop formula
     */
    template<typename T> 
    double sin2thetaEff(const T f, const double DeltaAlpha_i) const;
    
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
    const StandardModel& SM;
    
};

#endif	/* EWSMAPPROXIMATEFORMULAE_H */

