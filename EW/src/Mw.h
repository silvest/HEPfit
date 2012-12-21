/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MW_H
#define	MW_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"


class Mw : public ThObservable {
public:

    /**
     * @brief Mw constructor
     * @param[in] EW_i an object of EW class
     */
    Mw(const EW& EW_i) : ThObservable(EW_i), myEW(EW_i), 
            myEWTYPE(EW_i.getEWTYPE()) {
    };

    /**
     * @return the W-boson mass
     */
    double getThValue();

    
private:
    const EW& myEW;
    const EW::EWTYPE myEWTYPE;
};

#endif	/* MW_H */

