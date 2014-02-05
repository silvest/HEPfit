/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
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
 * @brief A class for the effective weak mixing angle @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$ 
 * for a leptonic channel at the @f$Z@f$ pole.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the effective weak mixing angle 
 * @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$ for a leptonic channel at the @f$Z@f$ pole.
 */
class sin2thetaEff : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] EW_i A reference to an object of EW class, which is the base 
     * class of the electroweak precision observables.
     */
    sin2thetaEff(const EW& EW_i) 
    : ThObservable(EW_i), myEW(EW_i)
    {
    };

    /**
     * @return The effective weak mixing angle @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$ 
     * for a leptonic channel at the @f$Z@f$ pole.
     */
    double computeThValue();

    
private:

    /**
     * A reference to an object of EW class, which is the base class of the 
     * electroweak precision observables.
     */
    const EW& myEW;///< A reference to an object of the EW class.
};

#endif	/* SIN2THETAEFF_H */

