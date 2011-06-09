/* 
 * File:   obliqueW.h
 * Author: mishima
 *
 * Created on June 9, 2011, 3:47 PM
 */

#ifndef OBLIQUEW_H
#define	OBLIQUEW_H

#include <ThObservable.h>
#include "EW.h"

class obliqueW : public ThObservable {
public:

    /**
     * @brief obliqueW constructor
     * @param[in] myEW an object of EW class
     */
    obliqueW(const EW& myEW);

    /**
     * @return the oblique parameter W
     */
    double getThValue();

private:

};

#endif	/* OBLIQUEW_H */

