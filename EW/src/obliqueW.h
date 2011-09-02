/* 
 * File:   obliqueW.h
 * Author: mishima
 */

#ifndef OBLIQUEW_H
#define	OBLIQUEW_H

#include <ThObservable.h>
#include "EW.h"

class obliqueW : public ThObservable {
public:

    /**
     * @brief obliqueW constructor
     * @param[in] EW_i an object of EW class
     */
    obliqueW(const EW& EW_i);

    /**
     * @return the oblique parameter W
     */
    double getThValue();

private:

};

#endif	/* OBLIQUEW_H */

