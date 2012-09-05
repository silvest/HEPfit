/* 
 * File:   ZFAFBcharm.h
 * Author: mishima
 */

#ifndef ZFAFBCHARM_H
#define	ZFAFBCHARM_H

#include <ThObservable.h>
#include "ZFitter.h"


class ZFAFBcharm : public ThObservable {
public:

    /**
     * @brief ZFAFBcharm constructor
     * @param[in] ZF_i an object of ZFitter class
     */
    ZFAFBcharm(const ZFitter& ZF_i) : ThObservable(ZF_i), myZF(ZF_i) {};

    /**
     * @return the forward-backward asymmetry of the b-bar channel
     */
    double getThValue();

    
private:
    const ZFitter& myZF;
};

#endif	/* ZFAFBCHARM_H */

