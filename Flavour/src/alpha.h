/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ALPHA_H
#define	ALPHA_H

#include "ThObservable.h"
#include "AmpDB2.h"

class Alpha : public ThObservable, AmpDB2 {
public:
    Alpha(const StandardModel& SM_i);

    double computeThValue();
};

#endif	/* ALPHA_H */
