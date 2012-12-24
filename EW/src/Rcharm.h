/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef RCHARM_H
#define	RCHARM_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"


class Rcharm : public ThObservable {
public:

    /**
     * @brief Rcharm constructor
     * @param[in] EW_i an object of EW class
     */
    Rcharm(const EW& EW_i) : ThObservable(EW_i), myEW(EW_i), 
            myEWTYPE(EW_i.getEWTYPE()) {
    };

    /**
     * @return the ratio of the c-cbar width to the hadronic width
     */
    double getThValue();

    
private:
    const EW& myEW;
    const EW::EWTYPE myEWTYPE;
};

#endif	/* RCHARM_H */

