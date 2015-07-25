/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ZFRCHARM_H
#define	ZFRCHARM_H

#include <ThObservable.h>
#include "ZFitterWrapper.h"


class ZFRcharm : public ThObservable {
public:

    /**
     * @brief ZFRcharm constructor
     * @param[in] ZF_i an object of ZFitter class
     */
    ZFRcharm(const ZFitterWrapper& ZF_i) : ThObservable(ZF_i), myZF(ZF_i) {};

    /**
     * @return the ratio of the c-cbar width to the hadronic width
     */
    double computeThValue();

    
private:
    const ZFitterWrapper& myZF;
};

#endif	/* ZFRCHARM_H */

