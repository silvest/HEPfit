/* 
 * File:   AFBbottom.h
 * Author: mishima
 */

#ifndef AFBBOTTOM_H
#define	AFBBOTTOM_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"


class AFBbottom : public ThObservable {
public:

    /**
     * @brief AFBbottom constructor
     * @param[in] EW_i an object of EW class
     * @param[in] type EWDEFAULT(default), EWCHMN, EWBURGESS or EWABC
     */
    AFBbottom(const EW& EW_i, const EW::EWTYPE type=EW::EWDEFAULT) : ThObservable(EW_i), 
            myEW(EW_i), myEWTYPE(type) {
    };

    /**
     * @return the forward-backward asymmetry of the b-bar channel
     */
    double getThValue();

    
private:
    const EW& myEW;
    const EW::EWTYPE myEWTYPE;
};

#endif	/* AFBBOTTOM_H */
