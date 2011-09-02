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
    sigmaHadron(const EW& EW_i);

    /**
     * @return the hadronic cross section 
     */
    double getThValue();

private:

};

#endif	/* SIGMAHADRON_H */

