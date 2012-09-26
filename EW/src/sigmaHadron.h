/* 
 * File:   sigmaHadron.h
 * Author: mishima
 */

#ifndef SIGMAHADRON_H
#define	SIGMAHADRON_H

#include <ThObservable.h>
#include "EW.h"


class sigmaHadron : public ThObservable {
public:

    /**
     * @brief sigmaHadron constructor
     * @param[in] EW_i an object of EW class
     */
    sigmaHadron(const EW& EW_i) : ThObservable(EW_i), myEW(EW_i) {};

    /**
     * @return the hadronic cross section in nb
     */
    double getThValue();


private:
    const EW& myEW;
};

#endif	/* SIGMAHADRON_H */

