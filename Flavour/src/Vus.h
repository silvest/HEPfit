/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef VUS_H
#define	VUS_H

#include <ThObservable.h>
#include <ThObsType.h>

class Vus : public ThObservable {
public:
    Vus(const ThObsType& ObsType) : ThObservable(ObsType) {};

    double getThValue();
};

#endif	/* VUS_H */
