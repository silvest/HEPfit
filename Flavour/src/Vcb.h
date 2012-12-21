/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef VCB_H
#define	VCB_H

#include <ThObservable.h>
#include <ThObsType.h>

class Vcb : public ThObservable {
public:
    Vcb(const ThObsType& ObsType) : ThObservable(ObsType) {};

    double getThValue();
};

#endif	/* VCB_H */
