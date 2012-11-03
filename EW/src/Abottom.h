/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ABOTTOM_H
#define	ABOTTOM_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"


class Abottom : public ThObservable {
public:

    /**
     * @brief Abottom constructor
     * @param[in] EW_i an object of EW class
     */
    Abottom(const EW& EW_i) : ThObservable(EW_i), myEW(EW_i), 
            myEWTYPE(EW_i.getEWTYPE()) {
    };

    /**
     * @return the left-right asymmetry of the b-bbar channel
     */
    double getThValue();

    
private:
    const EW& myEW;
    const EW::EWTYPE myEWTYPE;
};

#endif	/* ABOTTOM_H */

