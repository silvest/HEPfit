/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ZFABOTTOM_H
#define	ZFABOTTOM_H

#include <ThObservable.h>
#include "ZFitterWrapper.h"


class ZFAbottom : public ThObservable {
public:

    /**
     * @brief ZFAbottom constructor
     * @param[in] ZF_i an object of ZFitter class
     */
    ZFAbottom(const StandardModel& SM_i) : ThObservable(SM_i), myZF(SM_i) {};

    /**
     * @return the left-right asymmetry of the b-bbar channel
     */
    double computeThValue();

    
private:
    const ZFitterWrapper& myZF;
};

#endif	/* ZFABOTTOM_H */

