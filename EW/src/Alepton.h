/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ALEPTON_H
#define	ALEPTON_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"


class Alepton : public ThObservable {
public:

    /**
     * @brief Alepton constructor
     * @param[in] EW_i an object of EW class
     */
    Alepton(const EW& EW_i) : ThObservable(EW_i), myEW(EW_i), 
            myEWTYPE(EW_i.getEWTYPE()) {
    };

    /**
     * @return the left-right asymmetry of a leptonic channel
     */
    double getThValue();

    
private:
    const EW& myEW;
    const EW::EWTYPE myEWTYPE;
};

#endif	/* ALEPTON_H */

