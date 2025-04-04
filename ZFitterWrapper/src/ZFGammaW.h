/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ZFGAMMAW_H
#define	ZFGAMMAW_H

#include <ThObservable.h>
#include "ZFitterWrapper.h"


class ZFGammaW : public ThObservable {
public:

    /**
     * @brief ZFGammaW constructor
     * @param[in] ZF_i an object of ZFitter class
     */
    ZFGammaW(const StandardModel& SM_i) : ThObservable(SM_i), myZF(SM_i) {};

    /**
     * @return the total width of the W boson
     */
    double computeThValue();

    
private:
    const ZFitterWrapper& myZF;
};


#endif	/* ZFGAMMAW_H */

