/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ZFAFBCHARM_H
#define	ZFAFBCHARM_H

#include <ThObservable.h>
#include "ZFitterWrapper.h"


class ZFAFBcharm : public ThObservable {
public:

    /**
     * @brief ZFAFBcharm constructor
     * @param[in] ZF_i an object of ZFitter class
     */
    ZFAFBcharm(const ZFitterWrapper& ZF_i) : ThObservable(ZF_i), myZF(ZF_i) {};

    /**
     * @return the forward-backward asymmetry of the b-bar channel
     */
    double computeThValue();

    
private:
    const ZFitterWrapper& myZF;
};

#endif	/* ZFAFBCHARM_H */

