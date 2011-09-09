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
     * @param[in] DeltaAlpha_i radiative corrections to alpha
     */
    ApproximateFormulae(const StandardModel& SM_i, const double DeltaAlpha_i);    
    
    /**
     * @brief ApproximateFormulae copy constructor
     * @param[in] orig reference to an ApproximateFormulae object
     */    
    //ApproximateFormulae(const ApproximateFormulae& orig);
    
    /**
     * @brief ApproximateFormulae destructor
     */
    virtual ~ApproximateFormulae();

    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @return the W-boson mass from an approximate two-loop formula
     */
    double Mw() const;
    
    /**
     * @param[in] l name of a lepton
     * @return the effective weak mixing angle for Z->l\bar{l} from an approximate two-loop formula
     */
    double sin2thetaEff(const StandardModel::lepton l) const;

    /**
     * @param[in] q name of a quark
     * @return the effective weak mixing angle for Z->q\bar{q} from an approximate two-loop formula
     */
    double sin2thetaEff(const StandardModel::quark q) const;    
    
    /**
     * @return EW two-loop contribution to Delta r from an approximate formula
     */
    double DeltaR_TwoLoopEW() const;

    /**
     * @brief EW two-loop contribution to Delta kappa_Z^l = kappaZ^l - 1 
     * @return Delta kappa^{l, alpha^2} from an approximate formula
     */
    double DeltaKappa_l_TwoLoopEW() const;    

    /**
     * @brief EW two-loop contribution to Delta kappa_Z^b = kappaZ^b - 1 
     * @return Delta kappa^{b, alpha^2} from an approximate formula
     */
    double DeltaKappa_b_TwoLoopEW() const;    
    
    
    ////////////////////////////////////////////////////////////////////////
    
private:
    const StandardModel& SM;
    
    double myDeltaAlpha;
    
};

#endif	/* APPROXIMATEFORMULAE_H */

