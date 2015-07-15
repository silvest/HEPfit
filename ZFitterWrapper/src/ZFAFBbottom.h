/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ZFAFBBOTTOM_H
#define	ZFAFBBOTTOM_H

#include <ThObservable.h>
#include "ZFitterWrapper.h"


class ZFAFBbottom : public ThObservable {
public:

    /**
     * @brief ZFAFBbottom constructor
     * @param[in] ZF_i an object of ZFitter class
     */
    ZFAFBbottom(const ZFitterWrapper& ZF_i) : ThObservable(ZF_i), myZF(ZF_i) {};

    /**
     * @return the forward-backward asymmetry of the b-bar channel
     */
    double computeThValue();

    
private:
    const ZFitterWrapper& myZF;
};

#endif	/* ZFAFBBOTTOM_H */

