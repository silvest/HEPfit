/* 
 * File:   obliqueY.h
 * Author: mishima
 */

#ifndef OBLIQUEY_H
#define	OBLIQUEY_H

#include <ThObservable.h>
#include "EW.h"

class obliqueY : public ThObservable {
public:

    /**
     * @brief obliqueY constructor
     * @param[in] EW_i an object of EW class
     */
    obliqueY(const EW& EW_i);

    /**
     * @return the oblique parameter Y
     */
    double getThValue();

private:

};

#endif	/* OBLIQUEY_H */

