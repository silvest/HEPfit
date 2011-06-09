/* 
 * File:   obliqueEpsilon2.h
 * Author: mishima
 *
 * Created on June 9, 2011, 3:45 PM
 */

#ifndef OBLIQUEEPSILON2_H
#define	OBLIQUEEPSILON2_H

#include <ThObservable.h>
#include "EW.h"

class obliqueEpsilon2 : public ThObservable {
public:

    /**
     * @brief obliqueEpsilon2 constructor
     * @param[in] myEW an object of EW class
     */
    obliqueEpsilon2(const EW& myEW);

    /**
     * @return the oblique parameter epsilon_2
     */
    double getThValue();

private:

};

#endif	/* OBLIQUEEPSILON2_H */

