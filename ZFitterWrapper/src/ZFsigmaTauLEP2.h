/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ZFSIGMATAULEP2_H
#define	ZFSIGMATAULEP2_H

#include <ThObservable.h>
#include <StandardModel.h>
#include "ZFitterWrapper.h"


class ZFsigmaTauLEP2 : public ThObservable {
public:

    /**
     * @brief ZFsigmaTauLEP2 constructor
     * @param[in] ZF_i an object of ZFitter class
     * @param[in] sqrt_s_i \sqrt{s} of the e+ e- pair in the initial state
     */
    ZFsigmaTauLEP2(const ZFitterWrapper& ZF_i, const double sqrt_s_i) : ThObservable(ZF_i), 
            myZF(ZF_i), sqrt_s(sqrt_s_i) {};

    /**
     * @return the cross section for e^+ e^- -> tau^+ tau^- in pb at sqrt_s
     */
    double computeThValue();
    
private:
    const ZFitterWrapper& myZF;
    const double sqrt_s;
};

#endif	/* ZFSIGMATAULEP2_H */

