/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EPSILONB_H
#define	EPSILONB_H

#include <ThObservable.h>
#include "EW.h"


class epsilonb : public ThObservable {
public:

    /**
     * @brief epsilonb constructor
     * @param[in] EW_i an object of EW class
     */
    epsilonb(const EW& EW_i) : ThObservable(EW_i), myEW(EW_i) {};

    /**
     * @return epsilon_b
     */
    double getThValue();


private:
    const EW& myEW;
};

#endif	/* EPSILONB_H */

