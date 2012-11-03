/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EPSILON1_H
#define	EPSILON1_H

#include <ThObservable.h>
#include "EW.h"


class epsilon1 : public ThObservable {
public:

    /**
     * @brief epsilon1 constructor
     * @param[in] EW_i an object of EW class
     */
    epsilon1(const EW& EW_i) : ThObservable(EW_i), myEW(EW_i) {};

    /**
     * @return epsilon_1
     */
    double getThValue();


private:
    const EW& myEW;
};

#endif	/* EPSILON1_H */

