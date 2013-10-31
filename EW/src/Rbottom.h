/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef RBOTTOM_H
#define	RBOTTOM_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"

/**
 * @class Rbottom
 * @ingroup EW 
 * @brief A class for @f$R_b=\Gamma(Z\to b\bar{b})/\Gamma(Z\to {\rm hadrons})@f$ 
 * at the @f$Z@f$ pole.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the ratio of the @f$Z\to b\bar{b}@f$ 
 * width to the @f$Z@f$-boson hadronic width at the @f$Z@f$ pole. 
 */
class Rbottom : public ThObservable {
public:

    /**
     * A constructor.
     * @param[in] EW_i A reference to an object of EW class, which is the base 
     * class of the electroweak precision observables.
     */
    Rbottom(const EW& EW_i) 
    : ThObservable(EW_i), myEW(EW_i)
    {
    };

    /**
     * @return The ratio of the @f$Z\to b\bar{b}@f$ width to the @f$Z@f$-boson 
     * hadronic width at the @f$Z@f$ pole. 
     */
    double computeThValue();

    
private:

    /**
     * A reference to an object of EW class, which is the base class of the 
     * electroweak precision observables.
     */
    const EW& myEW;
};

#endif	/* RBOTTOM_H */

