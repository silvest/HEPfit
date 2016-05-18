/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ZFMH0_H
#define	ZFMH0_H

#include <ThObservable.h>
#include <StandardModel.h>
#include "ZFitterWrapper.h"


class ZFMh0 : public ThObservable {
public:

    /**
     * @brief ZFMh0 constructor
     * @param[in] ZF_i an object of ZFitter class
     */
    ZFMh0(const StandardModel& SM_i) : ThObservable(SM_i), myZF(SM_i) {};

    /**
     * @return the Higgs mass
     */
    double computeThValue();
    
private:
    const ZFitterWrapper& myZF;
};

#endif	/* ZFMH0_H */

