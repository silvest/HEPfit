/* 
 * File:   obliqueEpsilon2.h
 * Author: mishima
 */

#ifndef OBLIQUEEPSILON2_H
#define	OBLIQUEEPSILON2_H

#include <ThObservable.h>
#include "EW.h"

class obliqueEpsilon2 : public ThObservable {
public:

    /**
     * @brief obliqueEpsilon2 constructor
     * @param[in] EW_i an object of EW class
     */
    obliqueEpsilon2(const EW& EW_i);

    /**
     * @return the oblique parameter epsilon_2
     */
    double getThValue();

private:

};

#endif	/* OBLIQUEEPSILON2_H */

