/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ZFRBOTTOM_H
#define	ZFRBOTTOM_H

#include <ThObservable.h>
#include "ZFitterWrapper.h"


class ZFRbottom : public ThObservable {
public:

    /**
     * @brief ZFRbottom constructor
     * @param[in] ZF_i an object of ZFitter class
     */
    ZFRbottom(const ZFitterWrapper& ZF_i) : ThObservable(ZF_i), myZF(ZF_i) {};

    /**
     * @return the ratio of the b-bbar width to the hadronic width
     */
    double computeThValue();

    
private:
    const ZFitterWrapper& myZF;
};

#endif	/* ZFRBOTTOM_H */

