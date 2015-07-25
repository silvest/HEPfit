/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ALPHA_H
#define	ALPHA_H

#include <ThObservable.h>



class Alpha : public ThObservable {
public:
    Alpha(const StandardModel& SM_i) : ThObservable(SM_i) {};

    double computeThValue();
};

#endif	/* ALPHA_H */
