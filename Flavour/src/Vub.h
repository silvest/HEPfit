/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef VUB_H
#define	VUB_H

#include <ThObservable.h>
#include <ThObsType.h>

class Vub : public ThObservable {
public:
    Vub(const ThObsType& ObsType) : ThObservable(ObsType) {};

    double computeThValue();
};

#endif	/* VUB_H */
