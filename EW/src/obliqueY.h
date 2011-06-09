/* 
 * File:   obliqueY.h
 * Author: mishima
 *
 * Created on June 9, 2011, 3:47 PM
 */

#ifndef OBLIQUEY_H
#define	OBLIQUEY_H

#include <ThObservable.h>
#include "EW.h"

class obliqueY : public ThObservable {
public:

    /**
     * @brief obliqueY constructor
     * @param[in] myEW an object of EW class
     */
    obliqueY(const EW& myEW);

    /**
     * @return the oblique parameter Y
     */
    double getThValue();

private:

};

#endif	/* OBLIQUEY_H */

