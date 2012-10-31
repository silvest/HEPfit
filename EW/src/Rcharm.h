/* 
 * File:   Rcharm.h
 * Author: mishima
 */

#ifndef RCHARM_H
#define	RCHARM_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"


class Rcharm : public ThObservable {
public:

    /**
     * @brief Rcharm constructor
     * @param[in] EW_i an object of EW class
     * @param[in] type EWDEFAULT(default), EWCHMN, EWBURGESS or EWABC
     */
    Rcharm(const EW& EW_i, const EW::EWTYPE type=EW::EWDEFAULT) : ThObservable(EW_i), 
            myEW(EW_i), myEWTYPE(type) {
    };

    /**
     * @return the ratio of the c-cbar width to the hadronic width
     */
    double getThValue();

    
private:
    const EW& myEW;
    const EW::EWTYPE myEWTYPE;
};

#endif	/* RCHARM_H */

