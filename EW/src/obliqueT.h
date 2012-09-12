/* 
 * File:   obliqueT.h
 * Author: mishima
 */

#ifndef OBLIQUET_H
#define	OBLIQUET_H

#include <ThObservable.h>
#include "EW.h"


class obliqueT : public ThObservable {
public:

    /**
     * @brief obliqueT constructor
     * @param[in] EW_i an object of EW class
     */
    obliqueT(const EW& EW_i) : ThObservable(EW_i), myEW(EW_i) {};

    /**
     * @return the oblique parameter T
     */
    double getThValue();


private:
    const EW& myEW;
};

#endif	/* OBLIQUET_H */

