/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef SIN2THETAEFF_H
#define	SIN2THETAEFF_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"


class sin2thetaEff : public ThObservable {
public:

    /**
     * @brief sin2thetaEff constructor
     * @param[in] EW_i an object of EW class
     * @param[in] type EWDEFAULT(default), EWCHMN, EWBURGESS, EWABC or EWABC2
     */
    sin2thetaEff(const EW& EW_i, const EW::EWTYPE type=EW::EWDEFAULT) : ThObservable(EW_i), 
            myEW(EW_i), myEWTYPE(type) {
    };

    /**
     * @return the effective weak mixing angle for a leptonic channel
     */
    double getThValue();

    
private:
    const EW& myEW;
    const EW::EWTYPE myEWTYPE;
};

#endif	/* SIN2THETAEFF_H */

