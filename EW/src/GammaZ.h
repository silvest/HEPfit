/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GAMMAZ_H
#define	GAMMAZ_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"


class GammaZ : public ThObservable {
public:

    /**
     * @brief GammaZ constructor
     * @param[in] EW_i an object of EW class
     */
    GammaZ(const EW& EW_i) : ThObservable(EW_i), myEW(EW_i), 
            myEWTYPE(EW_i.getEWTYPE()) {
    };

    /**
     * @return the total width of the Z boson
     */
    double getThValue();

    
private:
    const EW& myEW;
    const EW::EWTYPE myEWTYPE;
};

#endif	/* GAMMAZ_H */

