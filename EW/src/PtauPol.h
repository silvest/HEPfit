/* 
 * File:   PtauPol.h
 * Author: mishima
 */

#ifndef PTAUPOL_H
#define	PTAUPOL_H

#include <ThObservable.h>
#include "EW.h"
#include "EW_CHMN.h"


class PtauPol : public ThObservable {
public:

    /**
     * @brief PtauPol constructor
     * @param[in] EW_i an object of EW class
     * @param[in] bCHMN_i true if using EW_CHMN class 
     */
    PtauPol(const EW& EW_i, const bool bCHMN_i=false) : ThObservable(EW_i), 
            myEW(EW_i), myEW_CHMN(EW_i.getSM()), bCHMN(bCHMN_i) {};

    /**
     * @return the longitudinal polarization of the tau-taubar channel
     */
    double getThValue();

    
private:
    const EW& myEW;
    const EW_CHMN myEW_CHMN;
    const bool bCHMN;
};

#endif	/* PTAUPOL_H */

