/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ZFABOTTOM_H
#define	ZFABOTTOM_H

#include <ThObservable.h>
#include "ZFitter.h"


class ZFAbottom : public ThObservable {
public:

    /**
     * @brief ZFAbottom constructor
     * @param[in] ZF_i an object of ZFitter class
     */
    ZFAbottom(const ZFitter& ZF_i) : ThObservable(ZF_i), myZF(ZF_i) {};

    /**
     * @return the left-right asymmetry of the b-bbar channel
     */
    double getThValue();

    
private:
    const ZFitter& myZF;
};

#endif	/* ZFABOTTOM_H */

