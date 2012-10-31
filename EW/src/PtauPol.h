/* 
 * File:   PtauPol.h
 * Author: mishima
 */

#ifndef PTAUPOL_H
#define	PTAUPOL_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"


class PtauPol : public ThObservable {
public:

    /**
     * @brief PtauPol constructor
     * @param[in] EW_i an object of EW class
     * @param[in] type EWDEFAULT(default), EWCHMN, EWBURGESS or EWABC
     */
    PtauPol(const EW& EW_i, const EW::EWTYPE type=EW::EWDEFAULT) : ThObservable(EW_i), 
            myEW(EW_i), myEWTYPE(type) {
    };

    /**
     * @return the longitudinal polarization of the tau-taubar channel
     */
    double getThValue();

    
private:
    const EW& myEW;
    const EW::EWTYPE myEWTYPE;
};

#endif	/* PTAUPOL_H */

