/* 
 * File:   ZFMw.h
 * Author: mishima
 */

#ifndef ZFMW_H
#define	ZFMW_H

#include <ThObservable.h>
#include <StandardModel.h>
#include "ZFitter.h"


class ZFMw : public ThObservable {
public:

    /**
     * @brief ZFMw constructor
     * @param[in] ZF_i an object of ZFitter class
     */
    ZFMw(const ZFitter& ZF_i) : ThObservable(ZF_i), myZF(ZF_i) {};

    /**
     * @return the W-boson mass
     */
    double getThValue();
    
private:
    const ZFitter& myZF;
};

#endif	/* ZFMW_H */

