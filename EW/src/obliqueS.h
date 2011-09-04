/* 
 * File:   obliqueS.h
 * Author: mishima
 */

#ifndef OBLIQUES_H
#define	OBLIQUES_H

#include <ThObservable.h>
#include "EW.h"
#include "obliqueEpsilon3.h"


class obliqueS : public ThObservable {
public:

    /**
     * @brief obliqueS constructor
     * @param[in] EW_i an object of EW class
     */
    obliqueS(const EW& EW_i);

    /**
     * @return the oblique parameter S
     */
    double getThValue();


private:
    obliqueEpsilon3 epsilon_3;
    double S;
    
};

#endif	/* OBLIQUES_H */

