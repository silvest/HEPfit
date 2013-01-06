/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef OBLIQUEU_H
#define	OBLIQUEU_H

#include <ThObservable.h>
#include "EW.h"

/**
 * @class obliqueU
 * @ingroup EW 
 * @brief A class for the oblique parameter U
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the oblique parameter U.  
 */
class obliqueU : public ThObservable {
public:

    /**
     * A constructor.
     * @param[in] EW_i A reference to an object of EW class, which is the base 
     * class of the electroweak precision observables.
     */
    obliqueU(const EW& EW_i) 
    : ThObservable(EW_i), myEW(EW_i) 
    {
    };

    /**
     * @return The oblique parameter U
     */
    double getThValue();

    
private:

    /**
     * A reference to an object of EW class, which is the base class of the 
     * electroweak precision observables.
     */
    const EW& myEW;
};

#endif	/* OBLIQUEU_H */

