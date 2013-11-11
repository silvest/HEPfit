/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MW_H
#define	MW_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"

/**
 * @class Mw 
 * @ingroup EW 
 * @brief A class for the W-boson mass
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the W-boson mass. 
 */
class Mw : public ThObservable {
public:

    /**
     * A constructor.
     * @param[in] EW_i A reference to an object of EW class, which is the base class of 
     * the electroweak precision observables.
     */
    Mw(const EW& EW_i) 
    : ThObservable(EW_i), myEW(EW_i)
    {
    };

    /**
     * @return The W-boson mass
     */
    double computeThValue();

    
private:

    /**
     * A reference to an object of EW class, which is the base class of the electroweak 
     * precision observables.
     */
    const EW& myEW;
};

#endif	/* MW_H */

