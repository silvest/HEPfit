/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ACHARM_H
#define	ACHARM_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"

/**
 * @class Acharm
 * @ingroup EW 
 * @brief A class for the left-right asymmetry of @f$Z\to c\bar{c}@f$ at the 
 * @f$Z@f$ pole.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the left-right asymmetry of the 
 * @f$Z\to c\bar{c}@f$ channel at the @f$Z@f$ pole. 
 */
class Acharm : public ThObservable {
public:

    /**
     * A constructor.
     * @param[in] EW_i A reference to an object of EW class, which is the base class of 
     * the electroweak precision observables.
     */
    Acharm(const EW& EW_i) 
    : ThObservable(EW_i), myEW(EW_i)
    {
    };

    /**
     * @return the left-right asymmetry of the c-cbar channel
     */
    double getThValue();

    
private:

    /**
     * A reference to an object of EW class, which is the base class of the electroweak 
     * precision observables.
     */
    const EW& myEW;
};

#endif	/* ACHARM_H */

