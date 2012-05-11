/* 
 * File:   ApproximateFormulae.h
 * Author: mishima
 */

#ifndef APPROXIMATEFORMULAE_H
#define	APPROXIMATEFORMULAE_H

#include <StandardModel.h>


class ApproximateFormulae {
public:
      
    /**
     * @brief  ApproximateFormulae constructor
     * @param[in] SM_i reference to a StandardModel object
     */
    ApproximateFormulae(const StandardModel& SM_i);    

    
    ////////////////////////////////////////////////////////////////////////

    /** 
     * @param[in] DeltaAlpha_i the total radiative corrections to alpha at Mz
     * @return the W-boson mass from an approximate two-loop formula
     */
    double Mw(const double DeltaAlpha_i) const;
    
    /**
     * @param[in] l name of a lepton 
     * @param[in] DeltaAlpha_i the total radiative corrections to alpha at Mz
     * @return the effective weak mixing angle for Z->l\bar{l} from an approximate two-loop formula
     */
    double sin2thetaEff(const StandardModel::lepton l, 
                        const double DeltaAlpha_i) const;

    /**
     * @param[in] q name of a quark 
     * @param[in] DeltaAlpha_i the total radiative corrections to alpha at Mz
     * @return the effective weak mixing angle for Z->q\bar{q} from an approximate two-loop formula
     */
    double sin2thetaEff(const StandardModel::quark q, 
                        const double DeltaAlpha_i) const;    
    
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

#endif	/* APPROXIMATEFORMULAE_H */

