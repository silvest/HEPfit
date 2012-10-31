/* 
 * File:   Rlepton.h
 * Author: mishima
 */

#ifndef RLEPTON_H
#define	RLEPTON_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"


class Rlepton : public ThObservable {
public:

    /**
     * @brief Rlepton constructor
     * @param[in] EW_i an object of EW class
     * @param[in] type EWDEFAULT(default), EWCHMN, EWBURGESS or EWABC
     */
    Rlepton(const EW& EW_i, const EW::EWTYPE type=EW::EWDEFAULT) : ThObservable(EW_i), 
            myEW(EW_i), myEWTYPE(type) {
    };

    /**
     * @return the ratio of the hadronic width to the leptonic width
     */
    double getThValue();

    
private:
    const EW& myEW; 
    const EW::EWTYPE myEWTYPE;
};

#endif	/* RLEPTON_H */

