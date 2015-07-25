/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ZFALEPTON_H
#define	ZFALEPTON_H

#include <ThObservable.h>
#include "ZFitterWrapper.h"


class ZFAlepton : public ThObservable {
public:

    /**
     * @brief ZFAlepton constructor
     * @param[in] ZF_i an object of ZFitter class
     */
    ZFAlepton(const ZFitterWrapper& ZF_i) : ThObservable(ZF_i), myZF(ZF_i) {};

    /**
     * @return the left-right asymmetry of a leptonic channel
     */
    double computeThValue();

    
private:
    const ZFitterWrapper& myZF;
};

#endif	/* ZFALEPTON_H */

