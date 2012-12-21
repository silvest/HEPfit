/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef RBOTTOM_H
#define	RBOTTOM_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"


class Rbottom : public ThObservable {
public:

    /**
     * @brief Rbottom constructor
     * @param[in] EW_i an object of EW class
     */
    Rbottom(const EW& EW_i) : ThObservable(EW_i), myEW(EW_i), 
            myEWTYPE(EW_i.getEWTYPE()) {
    };

    /**
     * @return the ratio of the b-bbar width to the hadronic width
     */
    double getThValue();

    
private:
    const EW& myEW;
    const EW::EWTYPE myEWTYPE;
};

#endif	/* RBOTTOM_H */

