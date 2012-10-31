/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ZFPTAUPOL_H
#define	ZFPTAUPOL_H

#include <ThObservable.h>
#include "ZFitter.h"


class ZFPtauPol : public ThObservable {
public:

    /**
     * @brief ZFPtauPol constructor
     * @param[in] ZF_i an object of ZFitter class
     */
    ZFPtauPol(const ZFitter& ZF_i) : ThObservable(ZF_i), myZF(ZF_i) {};

    /**
     * @return the longitudinal polarization of the tau-taubar channel
     */
    double getThValue();

    
private:
    const ZFitter& myZF;
};

#endif	/* ZFPTAUPOL_H */

