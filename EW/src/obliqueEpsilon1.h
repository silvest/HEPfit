/* 
 * File:   obliqueEpsilon1.h
 * Author: mishima
 */

#ifndef OBLIQUEEPSILON1_H
#define	OBLIQUEEPSILON1_H

#include <ThObservable.h>
#include "EW.h"

class obliqueEpsilon1 : public ThObservable {
public:

    /**
     * @brief obliqueEpsilon1 constructor
     * @param[in] EW_i an object of EW class
     */
    obliqueEpsilon1(const EW& EW_i);

    /**
     * @return the oblique parameter epsilon_1
     */
    double getThValue();

private:

};

#endif	/* OBLIQUEEPSILON1_H */

