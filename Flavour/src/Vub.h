/* 
 * Copyright (C) 2012 SUSYfit Collaboration
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

    double getThValue();
};

#endif	/* VUB_H */
