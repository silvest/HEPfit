/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef VUD_H
#define	VUD_H

#include <ThObservable.h>
#include <ThObsType.h>

class Vud : public ThObservable {
public:
    Vud(const ThObsType& ObsType) : ThObservable(ObsType) {};

    double getThValue();
};

#endif	/* VUD_H */
