/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EPSILON2_H
#define	EPSILON2_H

#include <ThObservable.h>
#include "EW.h"


class epsilon2 : public ThObservable {
public:

    /**
     * @brief epsilon2 constructor
     * @param[in] EW_i an object of EW class
     */
    epsilon2(const EW& EW_i) : ThObservable(EW_i), myEW(EW_i) {};

    /**
     * @return epsilon_2
     */
    double getThValue();


private:
    const EW& myEW;
};

#endif	/* EPSILON2_H */

