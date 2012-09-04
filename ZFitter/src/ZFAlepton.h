/* 
 * File:   ZFAlepton.h
 * Author: mishima
 */

#ifndef ZFALEPTON_H
#define	ZFALEPTON_H

#include <ThObservable.h>
#include "ZFitter.h"


class ZFAlepton : public ThObservable {
public:

    /**
     * @brief ZFAlepton constructor
     * @param[in] ZF_i an object of ZFitter class
     */
    ZFAlepton(const ZFitter& ZF_i) : ThObservable(ZF_i), myZF(ZF_i) {};

    /**
     * @return the left-right asymmetry of a leptonic channel
     */
    double getThValue();

    
private:
    const ZFitter& myZF;
};

#endif	/* ZFALEPTON_H */

