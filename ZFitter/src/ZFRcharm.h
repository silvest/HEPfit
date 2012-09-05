/* 
 * File:   ZFRcharm.h
 * Author: mishima
 */

#ifndef ZFRCHARM_H
#define	ZFRCHARM_H

#include <ThObservable.h>
#include "ZFitter.h"


class ZFRcharm : public ThObservable {
public:

    /**
     * @brief ZFRcharm constructor
     * @param[in] ZF_i an object of ZFitter class
     */
    ZFRcharm(const ZFitter& ZF_i) : ThObservable(ZF_i), myZF(ZF_i) {};

    /**
     * @return the ratio of the c-cbar width to the hadronic width
     */
    double getThValue();

    
private:
    const ZFitter& myZF;
};

#endif	/* ZFRCHARM_H */

