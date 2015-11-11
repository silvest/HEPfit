/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ZFDSIGMATAULEP2_H
#define	ZFDSIGMATAULEP2_H

#include <ThObservable.h>
#include <StandardModel.h>
#include "ZFitterWrapper.h"


class ZFDsigmaTauLEP2 : public ThObservable {
public:

    /**
     * @brief ZFDsigmaTauLEP2 constructor
     * @param[in] ZF_i an object of ZFitter class
     * @param[in] sqrt_s_i \sqrt{s} of the e+ e- pair in the initial state
     * @param[in] cos_theta_i cos(theta)
     */
    ZFDsigmaTauLEP2(const ZFitterWrapper& ZF_i, const double sqrt_s_i,  
                    const double cos_theta_i) : ThObservable(ZF_i), myZF(ZF_i), 
                    sqrt_s(sqrt_s_i), cos_theta(cos_theta_i) {};

    /**
     * @return the differential cross section for e^+ e^- -> tau^+ tau^- in pb at sqrt_s and cos_theta
     */
    double computeThValue();
    
private:
    const ZFitterWrapper& myZF;
    const double sqrt_s, cos_theta;
};


#endif	/* ZFDSIGMATAULEP2_H */

