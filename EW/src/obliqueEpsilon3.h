/* 
 * File:   obliqueEpsilon3.h
 * Author: mishima
 *
 * Created on June 9, 2011, 3:45 PM
 */

#ifndef OBLIQUEEPSILON3_H
#define	OBLIQUEEPSILON3_H

#include <ThObservable.h>
#include "EW.h"

class obliqueEpsilon3 : public ThObservable {
public:

    /**
     * @brief obliqueEpsilon3 constructor
     * @param[in] myEW an object of EW class
     */
    obliqueEpsilon3(const EW& myEW);

    /**
     * @return the oblique parameter epsilon_3
     */
    double getThValue();

private:

};

#endif	/* OBLIQUEEPSILON3_H */

