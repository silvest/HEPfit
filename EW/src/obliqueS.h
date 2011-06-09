/* 
 * File:   obliqueS.h
 * Author: mishima
 *
 * Created on June 9, 2011, 3:46 PM
 */

#ifndef OBLIQUES_H
#define	OBLIQUES_H

#include <ThObservable.h>
#include "EW.h"

class obliqueS : public ThObservable {
public:

    /**
     * @brief obliqueS constructor
     * @param[in] myEW an object of EW class
     */
    obliqueS(const EW& myEW);

    /**
     * @return the oblique parameter S
     */
    double getThValue();

private:

};

#endif	/* OBLIQUES_H */

