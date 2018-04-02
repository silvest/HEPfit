/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GAMMA_H
#define	GAMMA_H

#include "ThObservable.h"

class CKMGamma : public ThObservable {
public:
    CKMGamma(const StandardModel& SM_i);

    double computeThValue();
};

#endif	/* GAMMA_H */
