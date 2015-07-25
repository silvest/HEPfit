/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ZFSIGMAQUARKSLEP2_H
#define	ZFSIGMAQUARKSLEP2_H

#include <ThObservable.h>
#include "ZFitterWrapper.h"


class ZFsigmaQuarksLEP2 : public ThObservable {
public:

    /**
     * @brief ZFsigmaQuarksLEP2 constructor
     * @param[in] ZF_i an object of ZFitter class
     * @param[in] sqrt_s_i \sqrt{s} of the e+ e- pair in the initial state
     */
    ZFsigmaQuarksLEP2(const ZFitterWrapper& ZF_i, const double sqrt_s_i) : ThObservable(ZF_i), 
            myZF(ZF_i), sqrt_s(sqrt_s_i) {};

    /**
     * @return the hadronic cross section in pb at sqrt_s
     */
    double computeThValue();

    
private:
    const ZFitterWrapper& myZF;
    const double sqrt_s;
};

#endif	/* ZFSIGMAQUARKSLEP2_H */

