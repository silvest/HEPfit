/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef SIN2THETAEFF_H
#define	SIN2THETAEFF_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class sin2thetaEff 
 * @ingroup EW 
 * @brief An observable class for the leptonic effective weak mixing angle
 * @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$ at the @f$Z@f$ pole.
 * To be used for the electron or in case of asumming lepton universality.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the leptonic effective weak mixing
 * angle @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$ at the @f$Z@f$ pole.
 *
 * @sa EW_NPZff::sin2thetaEff() and the detailed description of EW class
 * for the inclusion of new physics contribution
 *
 */
class sin2thetaEff : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    sin2thetaEff(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The leptonic effective weak mixing angle 
     * @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$ at the the @f$Z@f$ pole.
     * @return @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$
     */
    double computeThValue();

    
private:


};


/**
 * @class sin2thetaEffel 
 * @ingroup EW 
 * @brief An observable class for the electron effective weak mixing angle
 * @f$\sin^2\theta_{\rm eff}^{ee}@f$ at the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the electron effective weak mixing
 * angle @f$\sin^2\theta_{ee}^{\rm lept}@f$ at the @f$Z@f$ pole.
 *
 */
class sin2thetaEffel : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    sin2thetaEffel(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The leptonic effective weak mixing angle 
     * @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$ at the the @f$Z@f$ pole.
     * @return @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$
     */
    double computeThValue();

    
private:


};



/**
 * @class sin2thetaEffmu 
 * @ingroup EW 
 * @brief An observable class for the muon effective weak mixing angle
 * @f$\sin^2\theta_{\rm eff}^{\mu}@f$ at the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the muon effective weak mixing
 * angle @f$\sin^2\theta_{\rm eff}^{\mu}@f$ at the @f$Z@f$ pole.
 *
 */
class sin2thetaEffmu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    sin2thetaEffmu(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The leptonic effective weak mixing angle 
     * @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$ at the the @f$Z@f$ pole.
     * @return @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$
     */
    double computeThValue();

    
private:


};

#endif	/* SIN2THETAEFF_H */

