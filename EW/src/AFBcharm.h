/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef AFBCHARM_H
#define	AFBCHARM_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"

/**
 * @class AFBcharm
 * @ingroup EW 
 * @brief A class for the forward-backward asymmetry of @f$Z\to c\bar{c}@f$ 
 * at the @f$Z@f$ pole.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the forward-backward asymmetry of 
 * the @f$Z\to c\bar{c}@f$ channel  at the @f$Z@f$ pole. 
 */
class AFBcharm : public ThObservable {
public:

    /**
     * A constructor.
     * @param[in] EW_i A reference to an object of EW class, which is the base class of 
     * the electroweak precision observables.
     */
    AFBcharm(const EW& EW_i) 
    : ThObservable(EW_i), myEW(EW_i), myEWTYPE(EW_i.getEWTYPE()) 
    {
    };
    
    /**
     * @return The forward-backward asymmetry of @f$Z\to c\bar{c}@f$ at the 
     * @f$Z@f$ pole.
     */
    double getThValue();

    
private:

    /**
     * A reference to an object of EW class, which is the base class of the electroweak 
     * precision observables.
     */
    const EW& myEW;

    /**
     * An enumerator controlling the formulae used in the computation.
     */
    const EW::EWTYPE myEWTYPE;
};

#endif	/* AFBCHARM_H */

