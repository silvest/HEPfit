/* 
 * File:   Alepton.h
 * Author: mishima
 */

#ifndef ALEPTON_H
#define	ALEPTON_H

#include <ThObservable.h>
#include "EW.h"

class Alepton : public ThObservable {
public:

    /**
     * @brief Alepton constructor
     * @param[in] EW_i an object of EW class
     */
    Alepton(const EW& EW_i);

    /**
     * @return the left-right asymmetry of a leptonic channel
     */
    double getThValue();

private:

};

#endif	/* ALEPTON_H */

