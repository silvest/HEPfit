/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BSGAMMATHDM_H
#define	BSGAMMATHDM_H

#include <stdexcept>
#include <ThObservable.h>
#include "THDM.h"
#include "THDMcache.h"

/**
 * @class bsgammaTHDM
 * @brief .
 */
class bsgammaTHDM : public ThObservable {
public:
    bsgammaTHDM(const StandardModel& SM_i);
    virtual ~bsgammaTHDM();
    double computeThValue();

protected:
    THDMcache * mycache;

private:
    const THDM * myTHDM;
    const StandardModel& mySM;
};

#endif	/* BSGAMMATHDM_H */
