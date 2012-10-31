/* 
 * File:   Mw.h
 * Author: mishima
 */

#ifndef MW_H
#define	MW_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"


class Mw : public ThObservable {
public:

    /**
     * @brief Mw constructor
     * @param[in] EW_i an object of EW class
     * @param[in] type EWDEFAULT(default), EWCHMN, EWBURGESS or EWABC
     */
    Mw(const EW& EW_i, const EW::EWTYPE type=EW::EWDEFAULT) : ThObservable(EW_i), 
            myEW(EW_i), myEWTYPE(type) {
    };

    /**
     * @return the W-boson mass
     */
    double getThValue();

    
private:
    const EW& myEW;
    const EW::EWTYPE myEWTYPE;
};

#endif	/* MW_H */

