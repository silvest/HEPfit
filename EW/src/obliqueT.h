/* 
 * File:   obliqueT.h
 * Author: mishima
 *
 * Created on June 9, 2011, 3:46 PM
 */

#ifndef OBLIQUET_H
#define	OBLIQUET_H

#include <ThObservable.h>
#include "EW.h"

class obliqueT : public ThObservable {
public:

    /**
     * @brief obliqueT constructor
     * @param[in] myEW an object of EW class
     */
    obliqueT(const EW& myEW);

    /**
     * @return the oblique parameter T
     */
    double getThValue();

private:

};

#endif	/* OBLIQUET_H */

