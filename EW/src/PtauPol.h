/* 
 * File:   PtauPol.h
 * Author: mishima
 */

#ifndef PTAUPOL_H
#define	PTAUPOL_H

#include <ThObservable.h>
#include "EW.h"


class PtauPol : public ThObservable {
public:

    /**
     * @brief PtauPol constructor
     * @param[in] EW_i an object of EW class
     */
    PtauPol(const EW& EW_i) : ThObservable(EW_i), myEW(EW_i) {};

    /**
     * @return the longitudinal polarization of the tau-taubar channel
     */
    double getThValue();

    
private:
    const EW& myEW;
};

#endif	/* PTAUPOL_H */

