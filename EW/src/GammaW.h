/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GAMMAW_H
#define	GAMMAW_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"


class GammaW : public ThObservable {
public:
    
    /**
     * @brief GammaW constructor
     * @param[in] EW_i an object of EW class
     */
    GammaW(const EW& EW_i) : ThObservable(EW_i), myEW(EW_i), 
            myEWTYPE(EW_i.getEWTYPE()) {
    };

    /**
     * @return the total width of the W boson 
     */
    double getThValue();

    
private:
    const EW& myEW;
    const EW::EWTYPE myEWTYPE;
};

#endif	/* GAMMAW_H */

