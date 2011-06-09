/* 
 * File:   Alepton.h
 * Author: mishima
 *
 * Created on June 9, 2011, 3:43 PM
 */

#ifndef ALEPTON_H
#define	ALEPTON_H

#include <ThObservable.h>
#include "EW.h"

class Alepton : public ThObservable {
public:

    /**
     * @brief Alepton constructor
     * @param[in] myEW an object of EW class
     */
    Alepton(const EW& myEW);

    /**
     * @return the left-right asymmetry of a leptonic channel
     */
    double getThValue();

private:

};

#endif	/* ALEPTON_H */

