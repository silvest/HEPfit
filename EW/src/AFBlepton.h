/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef AFBLEPTON_H
#define	AFBLEPTON_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"


class AFBlepton : public ThObservable {
public:
 
    /**
     * @brief AFBlepton constructor
     * @param[in] EW_i an object of EW class
     * @param[in] type EWDEFAULT(default), EWCHMN, EWBURGESS, EWABC or EWABC2
     */
    AFBlepton(const EW& EW_i, const EW::EWTYPE type=EW::EWDEFAULT) : ThObservable(EW_i), 
            myEW(EW_i), myEWTYPE(type) {
    };

    /**
     * @return the forward-backward asymmetry of a leptonic channel
     */
    double getThValue();

    
private:
    const EW& myEW;
    const EW::EWTYPE myEWTYPE;
};

#endif	/* AFBLEPTON_H */

