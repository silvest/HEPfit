/* 
 * File:   sigmaHadron.h
 * Author: mishima
 *
 * Created on June 9, 2011, 3:44 PM
 */

#ifndef SIGMAHADRON_H
#define	SIGMAHADRON_H

#include <ThObservable.h>
#include "EW.h"

class sigmaHadron : public ThObservable {
public:

    /**
     * @brief sigmaHadron constructor
     * @param[in] myEW an object of EW class
     */
    sigmaHadron(const EW& myEW);

    /**
     * @return the hadronic cross section 
     */
    double getThValue();

private:

};

#endif	/* SIGMAHADRON_H */

