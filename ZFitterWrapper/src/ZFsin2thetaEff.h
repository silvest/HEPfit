/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ZFSIN2THETAEFF_H
#define	ZFSIN2THETAEFF_H

#include <ThObservable.h>
#include "ZFitterWrapper.h"


class ZFsin2thetaEff : public ThObservable {
public:

    /**
     * @brief ZFsin2thetaEff constructor
     * @param[in] ZF_i an object of ZFitter class
     */
    ZFsin2thetaEff(const ZFitterWrapper& ZF_i) : ThObservable(ZF_i), myZF(ZF_i) {};

    /**
     * @return the effective weak mixing angle for a leptonic channel
     */
    double computeThValue();

    
private:
    const ZFitterWrapper& myZF;
};


#endif	/* ZFSIN2THETAEFF_H */

