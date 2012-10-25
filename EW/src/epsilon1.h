/* 
 * File:   epsilon1.h
 * Author: mishima
 */

#ifndef EPSILON1_H
#define	EPSILON1_H

#include <ThObservable.h>
#include "EW.h"


class epsilon1 : public ThObservable {
public:

    /**
     * @brief epsilon1 constructor
     * @param[in] EW_i an object of EW class
     */
    epsilon1(const EW& EW_i) : ThObservable(EW_i), myEW(EW_i) {};

    /**
     * @return the oblique parameter S
     */
    double getThValue();


private:
    const EW& myEW;
};

#endif	/* EPSILON1_H */

