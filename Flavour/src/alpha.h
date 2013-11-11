/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ALPHA_H
#define	ALPHA_H

#include <ThObservable.h>
#include <ThObsType.h>

class Alpha : public ThObservable {
public:
    Alpha(const ThObsType& ObsType) : ThObservable(ObsType) {};

    double computeThValue();
};

#endif	/* ALPHA_H */
