/* 
 * File:   PtauPol.h
 * Author: mishima
 *
 * Created on June 9, 2011, 3:42 PM
 */

#ifndef PTAUPOL_H
#define	PTAUPOL_H

#include <ThObservable.h>
#include "EW.h"

class PtauPol : public ThObservable {
public:

    /**
     * @brief PtauPol constructor
     * @param[in] myEW an object of EW class
     */
    PtauPol(const EW& myEW);

    /**
     * @return the longitudinal polarization of the tau-taubar channel
     */
    double getThValue();

private:

};

#endif	/* PTAUPOL_H */

