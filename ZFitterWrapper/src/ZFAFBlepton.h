/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ZFAFBLEPTON_H
#define	ZFAFBLEPTON_H

#include <ThObservable.h>
#include "ZFitterWrapper.h"


class ZFAFBlepton : public ThObservable {
public:

    /**
     * @brief ZFAFBlepton constructor
     * @param[in] ZF_i an object of ZFitter class
     */
    ZFAFBlepton(const ZFitterWrapper& ZF_i) : ThObservable(ZF_i), myZF(ZF_i) {};

    /**
     * @return the forward-backward asymmetry of a leptonic channel
     */
    double computeThValue();

    
private:
    const ZFitterWrapper& myZF;
};

#endif	/* ZFAFBLEPTON_H */

