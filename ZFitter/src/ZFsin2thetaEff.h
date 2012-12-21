/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ZFSIN2THETAEFF_H
#define	ZFSIN2THETAEFF_H

#include <ThObservable.h>
#include "ZFitter.h"


class ZFsin2thetaEff : public ThObservable {
public:

    /**
     * @brief ZFsin2thetaEff constructor
     * @param[in] ZF_i an object of ZFitter class
     */
    ZFsin2thetaEff(const ZFitter& ZF_i) : ThObservable(ZF_i), myZF(ZF_i) {};

    /**
     * @return the effective weak mixing angle for a leptonic channel
     */
    double getThValue();

    
private:
    const ZFitter& myZF;
};


#endif	/* ZFSIN2THETAEFF_H */

