/* 
 * File:   AFBcharm.h
 * Author: mishima
 */

#ifndef AFBCHARM_H
#define	AFBCHARM_H

#include <ThObservable.h>
#include "EW.h"


class AFBcharm : public ThObservable {
public:

    /**
     * @brief AFBcharm constructor
     * @param[in] EW_i an object of EW class
     */
    AFBcharm(const EW& EW_i) : ThObservable(EW_i), myEW(EW_i) {};

    /**
     * @return the forward-backward asymmetry of the c-cbar channel
     */
    double getThValue();

    
private:
    const EW& myEW;
};

#endif	/* AFBCHARM_H */

