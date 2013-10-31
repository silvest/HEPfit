/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef OBLIQUES_H
#define	OBLIQUES_H

#include <ThObservable.h>
#include "EW.h"

/**
 * @class obliqueS
 * @ingroup EW 
 * @brief A class for the oblique parameter S
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the oblique parameter S.  
 */
class obliqueS : public ThObservable {
public:

    /**
     * A constructor.
     * @param[in] EW_i A reference to an object of EW class, which is the base 
     * class of the electroweak precision observables.
     */
    obliqueS(const EW& EW_i) 
    : ThObservable(EW_i)//, myEW(EW_i) 
    {
    };

    /**
     * @return The oblique parameter S. 
     */
    double computeThValue();


private:

    /**
     * A reference to an object of EW class, which is the base class of the 
     * electroweak precision observables.
     */
    //const EW& myEW;
};

#endif	/* OBLIQUES_H */

