/* 
 * File:   ZFRbottomLEP2.h
 * Author: mishima
 */

#ifndef ZFRBOTTOMLEP2_H
#define	ZFRBOTTOMLEP2_H

#include <ThObservable.h>
#include <StandardModel.h>
#include "ZFitter.h"


class ZFRbottomLEP2 : public ThObservable {
public:

    /**
     * @brief ZFRbottomLEP2 constructor
     * @param[in] ZF_i an object of ZFitter class
     * @param[in] sqrt_s_i \sqrt{s} of the e+ e- pair in the initial state
     */
    ZFRbottomLEP2(const ZFitter& ZF_i, const double sqrt_s_i) : ThObservable(ZF_i), 
            myZF(ZF_i), sqrt_s(sqrt_s_i) {};

    /**
     * @return the ratio of the b-bbar cross section to the hadronic cross section at sqrt_s
     */
    double getThValue();
    
private:
    const ZFitter& myZF;
    const double sqrt_s;
};

#endif	/* ZFRBOTTOMLEP2_H */

