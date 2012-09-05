/* 
 * File:   Rcharm.h
 * Author: mishima
 */

#ifndef RCHARM_H
#define	RCHARM_H

#include <ThObservable.h>
#include "EW.h"


class Rcharm : public ThObservable {
public:

    /**
     * @brief Rcharm constructor
     * @param[in] EW_i an object of EW class
     */
    Rcharm(const EW& EW_i) : ThObservable(EW_i), myEW(EW_i) {};

    /**
     * @return the ratio of the c-cbar width to the hadronic width
     */
    double getThValue();

    
private:
    const EW& myEW;
};

#endif	/* RCHARM_H */

