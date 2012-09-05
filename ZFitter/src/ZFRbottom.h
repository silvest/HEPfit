/* 
 * File:   ZFRbottom.h
 * Author: mishima
 */

#ifndef ZFRBOTTOM_H
#define	ZFRBOTTOM_H

#include <ThObservable.h>
#include "ZFitter.h"


class ZFRbottom : public ThObservable {
public:

    /**
     * @brief ZFRbottom constructor
     * @param[in] ZF_i an object of ZFitter class
     */
    ZFRbottom(const ZFitter& ZF_i) : ThObservable(ZF_i), myZF(ZF_i) {};

    /**
     * @return the ratio of the b-bbar width to the hadronic width
     */
    double getThValue();

    
private:
    const ZFitter& myZF;
};

#endif	/* ZFRBOTTOM_H */

