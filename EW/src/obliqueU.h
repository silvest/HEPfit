/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef OBLIQUEU_H
#define	OBLIQUEU_H

#include <ThObservable.h>
#include "EW.h"


class obliqueU : public ThObservable {
public:

    /**
     * @brief obliqueU constructor
     * @param[in] EW_i an object of EW class
     */
    obliqueU(const EW& EW_i) : ThObservable(EW_i), myEW(EW_i) {};

    /**
     * @return the oblique parameter U
     */
    double getThValue();

    
private:
    const EW& myEW; 
};

#endif	/* OBLIQUEU_H */

