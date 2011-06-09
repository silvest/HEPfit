/* 
 * File:   obliqueU.h
 * Author: mishima
 *
 * Created on June 9, 2011, 3:46 PM
 */

#ifndef OBLIQUEU_H
#define	OBLIQUEU_H

#include <ThObservable.h>
#include "EW.h"

class obliqueU : public ThObservable {
public:

    /**
     * @brief obliqueU constructor
     * @param[in] myEW an object of EW class
     */
    obliqueU(const EW& myEW);

    /**
     * @return the oblique parameter U
     */
    double getThValue();

private:

};

#endif	/* OBLIQUEU_H */

