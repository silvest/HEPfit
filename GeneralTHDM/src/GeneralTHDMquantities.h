/* 
 * Copyright (C) 2016 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDMQUANTITIES_H
#define	GENERALTHDMQUANTITIES_H

//#include <stdexcept>
#include "ThObservable.h"
#include "GeneralTHDM.h"

/**
 * @class 
 * @ingroup GeneralTHDM 
 * @brief .
 */
class Re_sigma_u: public ThObservable {
public:

    /**
     * @brief Re_sigma_u constructor.
     */
    Re_sigma_u(const StandardModel& SM_i);

    /**
     * @return Re_sigma_u
     */
    double computeThValue();

    const GeneralTHDM * myGeneralTHDM;
};

#endif	/* GENERALTHDMQUANTITIES_H */
