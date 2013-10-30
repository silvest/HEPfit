/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ALPHA_2A_H
#define	ALPHA_2A_H

#include <ThObservable.h>
#include <ThObsType.h>

class Alpha_2a : public ThObservable {
public:
    Alpha_2a(const ThObsType& ObsType) : ThObservable(ObsType) {};

    double computeThValue();
};

#endif	/* ALPHA_2A_H */

