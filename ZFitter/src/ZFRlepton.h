/* 
 * File:   ZFRlepton.h
 * Author: mishima
 */

#ifndef ZFRLEPTON_H
#define	ZFRLEPTON_H

#include <ThObservable.h>
#include "ZFitter.h"


class ZFRlepton : public ThObservable {
public:

    /**
     * @brief ZFRlepton constructor
     * @param[in] ZF_i an object of ZFitter class
     */
    ZFRlepton(const ZFitter& ZF_i) : ThObservable(ZF_i), myZF(ZF_i) {};

    /**
     * @return the ratio of the hadronic width to the leptonic width
     */
    double getThValue();

    
private:
    const ZFitter& myZF;
};

#endif	/* ZFRLEPTON_H */

