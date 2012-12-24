/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ZFMH0_H
#define	ZFMH0_H

#include <ThObservable.h>
#include <StandardModel.h>
#include "ZFitter.h"


class ZFMh0 : public ThObservable {
public:

    /**
     * @brief ZFMh0 constructor
     * @param[in] ZF_i an object of ZFitter class
     */
    ZFMh0(const ZFitter& ZF_i) : ThObservable(ZF_i), myZF(ZF_i) {};

    /**
     * @return the Higgs mass
     */
    double getThValue();
    
private:
    const ZFitter& myZF;
};

#endif	/* ZFMH0_H */

