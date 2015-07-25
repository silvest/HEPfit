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
    ZFAbottom(const ZFitterWrapper& ZF_i) : ThObservable(ZF_i), myZF(ZF_i) {};

    /**
     * @return the left-right asymmetry of the b-bbar channel
     */
    double computeThValue();

    
private:
    const ZFitterWrapper& myZF;
};

#endif	/* ZFABOTTOM_H */

