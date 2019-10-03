/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ALPHA_2A_H
#define	ALPHA_2A_H

#include "ThObservable.h"
#include "AmpDB2.h"

class Alpha_2a : public ThObservable, AmpDB2 {
public:
    Alpha_2a(const StandardModel& SM_i) : ThObservable(SM_i), AmpDB2(SM_i) {};

    double computeThValue();
};

#endif	/* ALPHA_2A_H */

