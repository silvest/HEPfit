/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ZFSIGMAMULEP2_H
#define	ZFSIGMAMULEP2_H

#include <ThObservable.h>
#include <StandardModel.h>
#include "ZFitter.h"


class ZFsigmaMuLEP2 : public ThObservable {
public:

    /**
     * @brief ZFsigmaMuLEP2 constructor
     * @param[in] ZF_i an object of ZFitter class
     * @param[in] sqrt_s_i \sqrt{s} of the e+ e- pair in the initial state
     */
    ZFsigmaMuLEP2(const ZFitter& ZF_i, const double sqrt_s_i) : ThObservable(ZF_i), 
            myZF(ZF_i), sqrt_s(sqrt_s_i) {};

    /**
     * @return the cross section for e^+ e^- -> mu^+ mu^- in pb at sqrt_s
     */
    double getThValue();
    
private:
    const ZFitter& myZF;
    const double sqrt_s;
};

#endif	/* ZFSIGMAMULEP2_H */

