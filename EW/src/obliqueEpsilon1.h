/* 
 * File:   obliqueEpsilon1.h
 * Author: mishima
 *
 * Created on June 9, 2011, 3:45 PM
 */

#ifndef OBLIQUEEPSILON1_H
#define	OBLIQUEEPSILON1_H

#include <ThObservable.h>
#include "EW.h"

class obliqueEpsilon1 : public ThObservable {
public:

    /**
     * @brief obliqueEpsilon1 constructor
     * @param[in] myEW an object of EW class
     */
    obliqueEpsilon1(const EW& myEW);

    /**
     * @return the oblique parameter epsilon_1
     */
    double getThValue();

private:

};

#endif	/* OBLIQUEEPSILON1_H */

