/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef OBLIQUES_H
#define	OBLIQUES_H

#include <ThObservable.h>
#include "EW.h"


class obliqueS : public ThObservable {
public:

    /**
     * @brief obliqueS constructor
     * @param[in] EW_i an object of EW class
     */
    obliqueS(const EW& EW_i) : ThObservable(EW_i), myEW(EW_i) {};

    /**
     * @return the oblique parameter S
     */
    double getThValue();


private:
    const EW& myEW;
};

#endif	/* OBLIQUES_H */

