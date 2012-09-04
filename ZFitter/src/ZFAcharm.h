/* 
 * File:   ZFAcharm.h
 * Author: mishima
 */

#ifndef ZFACHARM_H
#define	ZFACHARM_H

#include <ThObservable.h>
#include "ZFitter.h"


class ZFAcharm : public ThObservable {
public:

    /**
     * @brief ZFAcharm constructor
     * @param[in] ZF_i an object of ZFitter class
     */
    ZFAcharm(const ZFitter& ZF_i) : ThObservable(ZF_i), myZF(ZF_i) {};

    /**
     * @return the left-right asymmetry of the c-cbar channel
     */
    double getThValue();

    
private:
    const ZFitter& myZF;
};

#endif	/* ZFACHARM_H */

