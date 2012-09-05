/* 
 * File:   ZFsigmaTauLEP2.h
 * Author: mishima
 */

#ifndef ZFSIGMATAULEP2_H
#define	ZFSIGMATAULEP2_H

#include <ThObservable.h>
#include <StandardModel.h>
#include "ZFitter.h"


class ZFsigmaTauLEP2 : public ThObservable {
public:

    /**
     * @brief ZFsigmaTauLEP2 constructor
     * @param[in] ZF_i an object of ZFitter class
     * @param[in] sqrt_s_i \sqrt{s} of the e+ e- pair in the initial state
     */
    ZFsigmaTauLEP2(const ZFitter& ZF_i, const double sqrt_s_i) : ThObservable(ZF_i), 
            myZF(ZF_i), sqrt_s(sqrt_s_i) {};

    /**
     * @return the cross section for e^+ e^- -> tau^+ tau^- in pb at sqrt_s
     */
    double getThValue();
    
private:
    const ZFitter& myZF;
    const double sqrt_s;
};

#endif	/* ZFSIGMATAULEP2_H */

