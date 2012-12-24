/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef AFBCHARM_H
#define	AFBCHARM_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"


class AFBcharm : public ThObservable {
public:

    /**
     * @brief AFBcharm constructor
     * @param[in] EW_i an object of EW class
     */
    AFBcharm(const EW& EW_i) : ThObservable(EW_i), myEW(EW_i), 
            myEWTYPE(EW_i.getEWTYPE()) {
    };
    /**
     * @return the forward-backward asymmetry of the c-cbar channel
     */
    double getThValue();

    
private:
    const EW& myEW;
    const EW::EWTYPE myEWTYPE;
};

#endif	/* AFBCHARM_H */

