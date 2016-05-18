/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ZFMW_H
#define	ZFMW_H

#include <ThObservable.h>
#include <StandardModel.h>
#include "ZFitterWrapper.h"


class ZFMw : public ThObservable {
public:

    /**
     * @brief ZFMw constructor
     * @param[in] ZF_i an object of ZFitter class
     */
    ZFMw(const StandardModel& SM_i) : ThObservable(SM_i), myZF(SM_i) {};

    /**
     * @return the W-boson mass
     */
    double computeThValue();
    
private:
    const ZFitterWrapper& myZF;
};

#endif	/* ZFMW_H */

