/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ZFGAMMAZ_H
#define	ZFGAMMAZ_H

#include <ThObservable.h>
#include "ZFitter.h"


class ZFGammaZ : public ThObservable {
public:

    /**
     * @brief ZFGammaZ constructor
     * @param[in] ZF_i an object of ZFitter class
     */
    ZFGammaZ(const ZFitter& ZF_i) : ThObservable(ZF_i), myZF(ZF_i) {};

    /**
     * @return the total width of the Z boson
     */
    double getThValue();

    
private:
    const ZFitter& myZF;
};

#endif	/* ZFGAMMAZ_H */

