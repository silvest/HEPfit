/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ZFGAMMAZ_H
#define	ZFGAMMAZ_H

#include <ThObservable.h>
#include "ZFitterWrapper.h"


class ZFGammaZ : public ThObservable {
public:

    /**
     * @brief ZFGammaZ constructor
     * @param[in] ZF_i an object of ZFitter class
     */
    ZFGammaZ(const StandardModel& SM_i) : ThObservable(SM_i), myZF(SM_i) {};

    /**
     * @return the total width of the Z boson
     */
    double computeThValue();

    
private:
    const ZFitterWrapper& myZF;
};

#endif	/* ZFGAMMAZ_H */

