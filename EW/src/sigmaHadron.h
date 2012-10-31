/* 
 * File:   sigmaHadron.h
 * Author: mishima
 */

#ifndef SIGMAHADRON_H
#define	SIGMAHADRON_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"


class sigmaHadron : public ThObservable {
public:

    /**
     * @brief sigmaHadron constructor
     * @param[in] EW_i an object of EW class
     * @param[in] type EWDEFAULT(default), EWCHMN, EWBURGESS or EWABC
     */
    sigmaHadron(const EW& EW_i, const EW::EWTYPE type=EW::EWDEFAULT) : ThObservable(EW_i), 
            myEW(EW_i), myEWTYPE(type) {
    };

    /**
     * @return the hadronic cross section in nb
     */
    double getThValue();


private:
    const EW& myEW;
    const EW::EWTYPE myEWTYPE;
};

#endif	/* SIGMAHADRON_H */

