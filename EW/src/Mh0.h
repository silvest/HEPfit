/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MH0_H
#define	MH0_H

#include <ThObservable.h>
#include "EW.h"

/**
 * @class Mh0
 * @ingroup EW 
 * @brief A class for the mass of the SM-like Higgs.  
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class Mh0 : public ThObservable {
public: 
    
    /**
     * A constructor.
     * @param[in] EW_i A reference to an object of EW class, which is the base 
     * class of the electroweak precision observables.
     */
    Mh0(const EW& EW_i) 
    : ThObservable(EW_i)
    { 
    };
    
    /**
     * @return The mass of the SM-like Higgs. 
     */
    double getThValue()
    {
        return SM.getMHl();
    }

    
private:

};

#endif	/* MH0_H */

