/* 
 * Copyright (C) 2012 SUSYfit Collaboration
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
     * @param[in] type EWDEFAULT(default), EWCHMN, EWBURGESS, EWABC or EWABC2
     */
    AFBcharm(const EW& EW_i, const EW::EWTYPE type=EW::EWDEFAULT) : ThObservable(EW_i), 
            myEW(EW_i), myEWTYPE(type) {
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

