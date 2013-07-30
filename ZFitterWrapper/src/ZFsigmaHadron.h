/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ZFSIGMAHADRON_H
#define	ZFSIGMAHADRON_H

#include <ThObservable.h>
#include "ZFitterWrapper.h"


class ZFsigmaHadron : public ThObservable {
public:

    /**
     * @brief ZFsigmaHadron constructor
     * @param[in] ZF_i an object of ZFitter class
     */
    ZFsigmaHadron(const ZFitterWrapper& ZF_i) : ThObservable(ZF_i), myZF(ZF_i) {};

    /**
     * @return the hadronic cross section in nb
     */
    double getThValue();

    
private:
    const ZFitterWrapper& myZF;
};

#endif	/* ZFSIGMAHADRON_H */

