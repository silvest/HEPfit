/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ACHARM_H
#define	ACHARM_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"


class Acharm : public ThObservable {
public:

    /**
     * @brief Acharm constructor
     * @param[in] EW_i an object of EW class
     */
    Acharm(const EW& EW_i) : ThObservable(EW_i), myEW(EW_i), 
            myEWTYPE(EW_i.getEWTYPE()) {
    };

    /**
     * @return the left-right asymmetry of the c-cbar channel
     */
    double getThValue();

    
private:
    const EW& myEW;
    const EW::EWTYPE myEWTYPE;
};

#endif	/* ACHARM_H */

