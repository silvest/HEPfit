/* 
 * File:   obliqueEpsilon3.h
 * Author: mishima
 */

#ifndef OBLIQUEEPSILON3_H
#define	OBLIQUEEPSILON3_H

#include <ThObservable.h>
#include "EW.h"


class obliqueEpsilon3 : public ThObservable {
public:

    /**
     * @brief obliqueEpsilon3 constructor
     * @param[in] EW_i an object of EW class
     */
    obliqueEpsilon3(const EW& EW_i);

    /**
     * @return the oblique parameter epsilon_3
     */
    double getThValue();

    
private:
    double epsilon_3;
    
};

#endif	/* OBLIQUEEPSILON3_H */

