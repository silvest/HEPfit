/* 
 * File:   ZFAFBmuLEP2.h
 * Author: mishima
 */

#ifndef ZFAFBMULEP2_H
#define	ZFAFBMULEP2_H

#include <ThObservable.h>
#include <StandardModel.h>
#include "ZFitter.h"


class ZFAFBmuLEP2 : public ThObservable {
public:

    /**
     * @brief ZFAFBmuLEP2 constructor
     * @param[in] ZF_i an object of ZFitter class
     * @param[in] sqrt_s_i \sqrt{s} of the e+ e- pair in the initial state
     */
    ZFAFBmuLEP2(const ZFitter& ZF_i, const double sqrt_s_i) : ThObservable(ZF_i), 
            myZF(ZF_i), sqrt_s(sqrt_s_i) {};

    /**
     * @return the forward-backward asymmetry for e^+ e^- -> mu^+ mu^- 
     */
    double getThValue();
    
private:
    const ZFitter& myZF;
    const double sqrt_s;
};

#endif	/* ZFAFBMULEP2_H */

