/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ZFAFBMULEP2_H
#define	ZFAFBMULEP2_H

#include <ThObservable.h>
#include <StandardModel.h>
#include "ZFitterWrapper.h"


class ZFAFBmuLEP2 : public ThObservable {
public:

    /**
     * @brief ZFAFBmuLEP2 constructor
     * @param[in] ZF_i an object of ZFitter class
     * @param[in] sqrt_s_i \sqrt{s} of the e+ e- pair in the initial state
     */
    ZFAFBmuLEP2(const ZFitterWrapper& ZF_i, const double sqrt_s_i) : ThObservable(ZF_i), 
            myZF(ZF_i), sqrt_s(sqrt_s_i) {};

    /**
     * @return the forward-backward asymmetry for e^+ e^- -> mu^+ mu^- 
     */
    double computeThValue();
    
private:
    const ZFitterWrapper& myZF;
    const double sqrt_s;
};

#endif	/* ZFAFBMULEP2_H */

