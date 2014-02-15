/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef SIN2THETAEFF_H
#define	SIN2THETAEFF_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"

/**
 * @class sin2thetaEff 
 * @ingroup EW 
 * @brief An observable class for the leptonic effective weak mixing angle
 * @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$ at the @f$Z@f$ pole.
 * @author SusyFit Collaboration
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
     * @param[in] EW_i a reference to an object of type EW
     */
    sin2thetaEff(const EW& EW_i) 
    : ThObservable(EW_i), myEW(EW_i)
    {
    };

    /**
     * @brief The leptonic effective weak mixing angle 
     * @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$ at the the @f$Z@f$ pole.
     * @return @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$
     */
    double computeThValue();

    
private:
    const EW& myEW;///< A reference to an object of type EW.

};

#endif	/* SIN2THETAEFF_H */

