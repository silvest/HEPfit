/* 
 * File:   Rbottom.h
 * Author: mishima
 */

#ifndef RBOTTOM_H
#define	RBOTTOM_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"


class Rbottom : public ThObservable {
public:

    /**
     * @brief Rbottom constructor
     * @param[in] EW_i an object of EW class
     * @param[in] type EWDEFAULT(default), EWCHMN, EWBURGESS or EWABC
     */
    Rbottom(const EW& EW_i, const EW::EWTYPE type=EW::EWDEFAULT) : ThObservable(EW_i), 
            myEW(EW_i), myEWTYPE(type) {
    };

    /**
     * @return the ratio of the b-bbar width to the hadronic width
     */
    double getThValue();

    
private:
    const EW& myEW;
    const EW::EWTYPE myEWTYPE;
};

#endif	/* RBOTTOM_H */

