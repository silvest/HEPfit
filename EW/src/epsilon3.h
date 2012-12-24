/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EPSILON3_H
#define	EPSILON3_H

#include <ThObservable.h>
#include "EW.h"


class epsilon3 : public ThObservable {
public:

    /**
     * @brief epsilon3 constructor
     * @param[in] EW_i an object of EW class
     */
    epsilon3(const EW& EW_i) : ThObservable(EW_i), myEW(EW_i) {};

    /**
     * @return epsilon_3
     */
    double getThValue();


private:
    const EW& myEW;
};

#endif	/* EPSILON3_H */

