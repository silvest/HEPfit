/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ZFRLEPTON_H
#define	ZFRLEPTON_H

#include <ThObservable.h>
#include "ZFitterWrapper.h"


class ZFRlepton : public ThObservable {
public:

    /**
     * @brief ZFRlepton constructor
     * @param[in] ZF_i an object of ZFitter class
     */
    ZFRlepton(const StandardModel& SM_i) : ThObservable(SM_i), myZF(SM_i) {};

    /**
     * @return the ratio of the hadronic width to the leptonic width
     */
    double computeThValue();

    
private:
    const ZFitterWrapper& myZF;
};

#endif	/* ZFRLEPTON_H */

